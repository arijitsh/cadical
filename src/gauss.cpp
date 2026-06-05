#include "internal.hpp"

#include "gauss/packed_matrix.hpp"

namespace CaDiCaL {

// Build the Gauss-Jordan matrix from the stored external XOR constraints,
// reduce it to row echelon form, deal with the resulting unit / empty rows at
// the root level, and keep the multi-variable rows for use during search.
//
// Each external XOR encodes 'XOR(literals) == true' (a negative literal flips
// the right-hand-side parity).  We translate every external literal to its
// internal variable via the external 'e2i' mapping; the sign of the literal is
// preserved by 'internalize' (which 'External::add_xor_clause' already called),
// so a negative literal simply flips the row's parity.
void Internal::init_gauss () {

  if (gauss)
    return;
  if (!opts.gauss)
    return;
  if (external->xors.empty ())
    return;
  // The engine does not (yet) emit proof steps for its derived clauses, so it
  // only runs when no proof / LRAT is being produced.  In proof mode the
  // 'xorblast' fallback keeps the solver correct.
  if (proof || lrat) {
    LOG ("not activating Gauss-Jordan engine: proof/LRAT enabled");
    return;
  }
  assert (!level);

  // Collect the involved internal variables and assign them matrix columns.
  // 'var_to_col' is indexed by internal variable index.
  std::vector<int> var_to_col (max_var + 1, -1);
  std::vector<int> col_to_var;
  std::vector<std::vector<int>> rows;   // internal var indices per input row
  std::vector<signed char> rows_rhs;

  for (const auto &xlits : external->xors) {
    std::vector<int> vars;
    int rhs = 1; // 'XOR(literals) == true'
    bool ok = true;
    for (const int elit : xlits) {
      const int eidx = (elit < 0) ? -elit : elit;
      if (eidx >= (int) external->e2i.size ()) {
        ok = false;
        break;
      }
      // 'e2i' maps an external variable to an internal *literal* which may be
      // negative if the variable was substituted by an equivalent one during
      // inprocessing.  Compute the internal literal of this XOR literal and
      // fold its sign into the right-hand side.
      int ilit = external->e2i[eidx];
      if (!ilit) {
        ok = false;
        break;
      }
      if (elit < 0)
        ilit = -ilit;
      if (ilit < 0)
        rhs ^= 1;
      const int ivar = (ilit < 0) ? -ilit : ilit;
      if (ivar > max_var) {
        ok = false;
        break;
      }
      vars.push_back (ivar);
      if (var_to_col[ivar] < 0) {
        var_to_col[ivar] = (int) col_to_var.size ();
        col_to_var.push_back (ivar);
      }
    }
    if (!ok)
      continue;
    rows.push_back (std::move (vars));
    rows_rhs.push_back ((signed char) rhs);
  }

  const uint32_t num_rows = (uint32_t) rows.size ();
  const uint32_t num_cols = (uint32_t) col_to_var.size ();
  if (!num_rows || !num_cols)
    return;

  // Fill the bit-matrix.
  Gauss::PackedMatrix mat;
  mat.resize (num_rows, num_cols);
  for (uint32_t r = 0; r < num_rows; r++) {
    Gauss::PackedRow pr = mat[r];
    pr.setZero ();
    pr.rhs () = rows_rhs[r];
    for (const int v : rows[r])
      pr.flipBit ((uint32_t) var_to_col[v]); // XOR-toggle: duplicates cancel
  }

  // Gaussian elimination to reduced row echelon form.
  uint32_t pivot = 0;
  for (uint32_t c = 0; c < num_cols && pivot < num_rows; c++) {
    int sel = -1;
    for (uint32_t r = pivot; r < num_rows; r++)
      if (mat[r][c]) {
        sel = (int) r;
        break;
      }
    if (sel < 0)
      continue;
    if ((uint32_t) sel != pivot)
      mat[pivot].swapBoth (mat[sel]);
    for (uint32_t r = 0; r < num_rows; r++)
      if (r != pivot && mat[r][c])
        mat[r].xor_in (mat[pivot]);
    pivot++;
  }

  // Extract rows: handle units / empty rows now, keep multi-variable rows.
  gauss = new GaussMatrix ();
  for (uint32_t r = 0; r < num_rows; r++) {
    Gauss::PackedRow pr = mat[r];
    const uint32_t pc = pr.popcnt ();
    const int rhs = (int) (pr.rhs () & 1);
    if (pc == 0) {
      if (rhs) {
        LOG ("Gauss-Jordan derived empty XOR row -> UNSAT");
        learn_empty_clause ();
        return;
      }
      continue; // 0 == 0, redundant
    }
    // Gather the variables of this row.
    std::vector<int> vars;
    for (uint32_t c = 0; c < num_cols; c++)
      if (pr[c])
        vars.push_back (col_to_var[c]);
    assert (vars.size () == pc);
    if (pc == 1) {
      const int v = vars[0];
      const int forced = rhs ? v : -v; // value making XOR(v)==rhs hold
      const signed char cur = val (v);
      if (!cur)
        assign_unit (forced);
      else if ((cur > 0) != (forced > 0)) {
        LOG ("Gauss-Jordan unit row conflicts with assignment -> UNSAT");
        learn_empty_clause ();
        return;
      }
      continue;
    }
    gauss->row_vars.push_back (std::move (vars));
    gauss->row_rhs.push_back ((signed char) rhs);
  }

  // Build the variable -> rows occurrence index for incremental evaluation.
  gauss->var_rows.assign (max_var + 1, {});
  for (size_t r = 0; r < gauss->row_vars.size (); r++)
    for (const int v : gauss->row_vars[r])
      gauss->var_rows[v].push_back ((int) r);
  gauss->qhead = 0;

  // After preprocessing the watch lists may be disconnected.  Reconnect them
  // so that (a) the unit literals just assigned by 'assign_unit' are actually
  // propagated against the CNF clauses once search starts (connect_watches
  // rewinds 'propagated' across falsified watches at level 0) and (b) the
  // reason clauses created later by 'gauss_round' interact correctly with the
  // existing clauses.
  clear_watches ();
  connect_watches ();

  VERBOSE (2,
           "Gauss-Jordan: %zu input rows over %u variables, %zu active rows "
           "after reduction",
           (size_t) num_rows, num_cols, gauss->row_vars.size ());
}

void Internal::reset_gauss () {
  if (!gauss)
    return;
  delete gauss;
  gauss = nullptr;
}

// Create an implied clause for row 'r' under the current assignment.  The
// clause buffer 'gauss->reason' is filled with 'forced' (the propagated
// literal, or 0 for a pure conflict clause) followed by the currently-false
// literal of every other variable in the row.  The two watched literals are
// placed at positions 0 and 1: for a propagation, position 0 is the propagated
// literal and position 1 is the highest-level false literal; for a conflict the
// two highest-level literals are used.  The clause is exactly one clause of the
// row's CNF expansion and is therefore implied by the original XOR system.
Clause *Internal::gauss_build_clause (const std::vector<int> &vars, int forced,
                                      int forced_var) {
  auto &c = gauss->reason;
  c.clear ();
  if (forced)
    c.push_back (forced);
  for (const int v : vars) {
    if (v == forced_var)
      continue;
    const signed char s = val (v);
    assert (s); // every other variable is assigned
    c.push_back ((s > 0) ? -v : v); // the currently-false literal
  }
  assert (c.size () >= 2);

  // Order watches.  Position 0 must be the propagated literal (if any), which
  // is either unassigned (propagation) or will be the asserting literal.
  // Position 1 must be a highest-level *false* literal so the two-watched
  // invariant holds once 'forced' becomes true.
  size_t first_false = forced ? 1 : 0;
  size_t best = first_false;
  int best_level = -1;
  for (size_t i = first_false; i < c.size (); i++) {
    const int lvl = var (c[i]).level;
    if (lvl > best_level) {
      best_level = lvl;
      best = i;
    }
  }
  std::swap (c[first_false], c[best]);
  if (!forced) {
    // Pure conflict clause: also pull the next highest-level literal to pos 1.
    size_t best2 = 1;
    int best2_level = -1;
    for (size_t i = 1; i < c.size (); i++) {
      const int lvl = var (c[i]).level;
      if (lvl > best2_level) {
        best2_level = lvl;
        best2 = i;
      }
    }
    std::swap (c[1], c[best2]);
  }

  // Build the clause through the standard staging vector.  The reason/conflict
  // clauses are redundant (and thus collectable by 'reduce' once they are no
  // longer reasons) so that a search with many Gauss-Jordan propagations does
  // not accumulate unboundedly many permanent clauses.
  clause = c;
  const int glue = (int) clause.size ();
  Clause *res = new_clause (true, glue);
  clause.clear ();
  watch_clause (res);
  return res;
}

// Verify a candidate complete assignment against every active XOR row.  The
// incremental 'gauss_round' may miss propagations that only become unit after
// backtracking, so before accepting a model we rescan all rows.  If a row is
// falsified, install it as a conflict and return false so search continues;
// returns true when every row is satisfied.  Called only at candidate models,
// so the full scan is cheap.
bool Internal::gauss_check_model () {
  if (!gauss)
    return true;
  for (size_t r = 0; r < gauss->row_vars.size (); r++) {
    const std::vector<int> &vars = gauss->row_vars[r];
    int parity = 0;
    bool full = true;
    for (const int v : vars) {
      const signed char s = val (v);
      if (!s) {
        full = false;
        break;
      }
      if (s > 0)
        parity ^= 1;
    }
    if (full && parity != gauss->row_rhs[r]) {
      conflict = gauss_build_clause (vars, 0, 0);
      gauss->conflicts++;
      return false;
    }
  }
  return true;
}

// Evaluate a single active row against the current assignment.  Propagates the
// unique unassigned literal (if any) or sets 'conflict' if the row is fully
// assigned with the wrong parity.  Returns true if it propagated a literal.
bool Internal::gauss_eval_row (size_t r) {
  const std::vector<int> &vars = gauss->row_vars[r];
  int count_unassigned = 0;
  int unassigned_var = 0;
  int assigned_parity = 0;
  for (const int v : vars) {
    const signed char s = val (v);
    if (!s) {
      if (++count_unassigned >= 2)
        return false; // not unit yet, nothing to do
      unassigned_var = v;
    } else if (s > 0)
      assigned_parity ^= 1;
  }
  if (count_unassigned == 1) {
    const int want = gauss->row_rhs[r] ^ assigned_parity; // value of last var
    const int forced = want ? unassigned_var : -unassigned_var;
    Clause *reason = gauss_build_clause (vars, forced, unassigned_var);
    search_assign_driving (forced, reason);
    gauss->propagations++;
    return true;
  }
  // fully assigned
  if (assigned_parity != gauss->row_rhs[r]) {
    conflict = gauss_build_clause (vars, 0, 0);
    gauss->conflicts++;
  }
  return false;
}

// Incrementally feed newly assigned trail literals to the matrix: for every
// assignment since the last call, examine only the rows that contain that
// variable.  Returns true if at least one literal was propagated (so the
// caller should re-run BCP); sets 'conflict' on a falsified row.
bool Internal::gauss_round () {

  if (!gauss)
    return false;

  // The trail may have shrunk due to backtracking since the last call.
  if (gauss->qhead > trail.size ())
    gauss->qhead = trail.size ();

  bool progress = false;
  while (!conflict && gauss->qhead < trail.size ()) {
    const int lit = trail[gauss->qhead++];
    const int v = (lit < 0) ? -lit : lit;
    if (v >= (int) gauss->var_rows.size ())
      continue;
    for (const int r : gauss->var_rows[v]) {
      if (gauss_eval_row ((size_t) r))
        progress = true;
      if (conflict)
        break;
    }
  }
  return progress;
}

} // namespace CaDiCaL
