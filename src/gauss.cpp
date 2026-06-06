#include "internal.hpp"

namespace CaDiCaL {

// ===========================================================================
// Incremental "watched-column" Gauss-Jordan XOR engine (port of CryptoMiniSat
// EGaussian onto CaDiCaL).  See 'gauss.hpp' for the high-level idea.  Reason
// and conflict clauses are real CaDiCaL clauses (one clause of the responsible
// row's CNF expansion), so conflict analysis is reused unchanged; the matrix
// rows always stay valid GF(2) combinations of the original XORs, so conflicts
// are sound, and the 'gauss_check_model' gate guarantees correct models.
// ===========================================================================

void Internal::reset_gauss () {
  if (!gauss)
    return;
  delete gauss;
  gauss = nullptr;
}

// Called from 'backtrack' on every real backtrack: the matrix and watches are
// valid GF(2) state and persist, but the column-assignment bitsets and the
// per-row 'satisfied' flags must be invalidated (CryptoMiniSat's 'canceling').
void Internal::gauss_notify_backtrack () {
  if (!gauss)
    return;
  gauss->cancelled_since_val_update = true;
  std::fill (gauss->satisfied_xors.begin (), gauss->satisfied_xors.end (),
             (char) 0);
}

// Build the matrix from the external XOR constraints, reduce it to RREF (which
// fixes each row's basic/responsible variable), deal with unit/empty rows at
// the root level, and set up two watches per remaining row.
void Internal::init_gauss () {

  if (!opts.gauss)
    return;
  if (external->xors.empty ())
    return;
  if (proof || lrat) {
    LOG ("not activating Gauss-Jordan engine: proof/LRAT enabled");
    return;
  }
  // Reuse an existing matrix only if no XOR clauses were added since it was
  // built; otherwise rebuild (model counters add random XOR hashes between
  // solves).
  if (gauss) {
    if (gauss->num_input_xors == external->xors.size ())
      return;
    reset_gauss ();
  }
  assert (!level);

  std::vector<int> var_to_col (max_var + 1, -1);
  std::vector<int> col_to_var;
  std::vector<std::vector<int>> rows;
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

  gauss = new GaussMatrix ();
  gauss->num_input_xors = external->xors.size ();
  gauss->num_cols = num_cols;
  gauss->col_to_var = col_to_var;
  gauss->var_to_col = var_to_col;
  gauss->var_has_resp_row.assign (max_var + 1, 0);
  gauss->gwatch.assign (max_var + 1, {});

  Gauss::PackedMatrix &mat = gauss->mat;
  mat.resize (num_rows, num_cols);
  for (uint32_t r = 0; r < num_rows; r++) {
    Gauss::PackedRow pr = mat[r];
    pr.setZero ();
    pr.rhs () = rows_rhs[r];
    for (const int v : rows[r])
      pr.flipBit ((uint32_t) var_to_col[v]); // XOR-toggle: duplicates cancel
  }

  // Gauss-Jordan elimination (reduced).  Establish each row's basic variable.
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
    gauss->var_has_resp_row[col_to_var[c]] = 1; // basic variable of this row
    if ((uint32_t) sel != pivot)
      mat[pivot].swapBoth (mat[sel]);
    for (uint32_t r = 0; r < num_rows; r++)
      if (r != pivot && mat[r][c])
        mat[r].xor_in (mat[pivot]);
    pivot++;
  }

  gauss->num_rows = num_rows;
  gauss->row_to_var_non_resp.assign (num_rows, GaussMatrix::no_var);
  gauss->satisfied_xors.assign (num_rows, 0);
  gauss->aux.resize (4, num_cols); // cols_vals, cols_unset, tmp_col, tmp_col2

  // init_adjust_matrix: handle unit/empty rows, set up watches for the rest.
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
      gauss->satisfied_xors[r] = 1;
      continue;
    }
    if (pc == 1) {
      const int col = pr.first_set_bit ();
      const int v = col_to_var[col];
      const int forced = rhs ? v : -v;
      const signed char cur = val (v);
      if (!cur)
        assign_unit (forced);
      else if ((cur > 0) != (forced > 0)) {
        LOG ("Gauss-Jordan unit row conflicts with assignment -> UNSAT");
        learn_empty_clause ();
        return;
      }
      pr.setZero ();
      gauss->var_has_resp_row[v] = 0;
      gauss->satisfied_xors[r] = 1;
      continue;
    }
    // pc >= 2: find the basic variable and one non-basic variable to watch.
    uint32_t basic_var = GaussMatrix::no_var, non_resp = GaussMatrix::no_var;
    pr.for_each_set_bit ([&] (uint32_t col) {
      const int v = col_to_var[col];
      if (gauss->var_has_resp_row[v])
        basic_var = v;
      else
        non_resp = v;
    });
    assert (basic_var != GaussMatrix::no_var);
    assert (non_resp != GaussMatrix::no_var);
    gauss->gwatch[basic_var].push_back (r);
    gauss->gwatch[non_resp].push_back (r);
    gauss->row_to_var_non_resp[r] = non_resp;

    // Snapshot the original active row for the model-check gate.
    std::vector<int> vs;
    pr.for_each_set_bit (
        [&] (uint32_t col) { vs.push_back (col_to_var[col]); });
    gauss->row_vars.push_back (std::move (vs));
    gauss->row_rhs.push_back ((signed char) rhs);
  }

  gauss->cancelled_since_val_update = true;
  gauss->last_val_update = trail.size ();
  gauss->qhead = 0;
  gauss->last_trail = trail.size ();

  // Watch lists may be disconnected after preprocessing; reconnect so the
  // unit literals assigned above are propagated against the CNF and the reason
  // clauses created later interact correctly with the existing clauses.
  clear_watches ();
  connect_watches ();

  VERBOSE (2, "Gauss-Jordan: %u rows over %u variables, %zu active",
           num_rows, num_cols, gauss->row_vars.size ());
}

// --------------------------------------------------------------------------
// Maintain the column assignment bitsets aux[0]=cols_vals, aux[1]=cols_unset.
// 'just_assigned_var >= 0': reflect that single variable immediately.
// 'just_assigned_var < 0' : lazily catch up to the current trail (a full
// rebuild after backtracking, otherwise incrementally).
// --------------------------------------------------------------------------
void Internal::gauss_update_cols (int just_assigned_var) {
  Gauss::PackedRow cols_vals = gauss->aux[0];
  Gauss::PackedRow cols_unset = gauss->aux[1];

  if (just_assigned_var >= 0) {
    const int c = gauss->var_to_col[just_assigned_var];
    if (c < 0)
      return;
    cols_unset.clearBit ((uint32_t) c);
    if (val (just_assigned_var) > 0)
      cols_vals.setBit ((uint32_t) c);
    else
      cols_vals.clearBit ((uint32_t) c);
    return;
  }

  if (gauss->cancelled_since_val_update) {
    cols_vals.setZero ();
    cols_unset.setOne ();
    for (uint32_t c = 0; c < gauss->num_cols; c++) {
      const signed char s = val (gauss->col_to_var[c]);
      if (!s)
        continue;
      cols_unset.clearBit (c);
      if (s > 0)
        cols_vals.setBit (c);
    }
    gauss->last_val_update = trail.size ();
    gauss->cancelled_since_val_update = false;
    return;
  }

  for (size_t i = gauss->last_val_update; i < trail.size (); i++) {
    const int lit = trail[i];
    const int v = (lit < 0) ? -lit : lit;
    // Variables added after this matrix was built (e.g. by an incremental
    // model counter between solves) are not part of it; index by the engine's
    // own table size, not the possibly-larger current 'max_var'.
    if ((size_t) v >= gauss->var_to_col.size ())
      continue;
    const int c = gauss->var_to_col[v];
    if (c < 0)
      continue;
    cols_unset.clearBit ((uint32_t) c);
    if (val (v) > 0)
      cols_vals.setBit ((uint32_t) c);
    else
      cols_vals.clearBit ((uint32_t) c);
  }
  gauss->last_val_update = trail.size ();
}

// Build a reason/conflict clause from the current bit pattern of 'row'.
Clause *Internal::gauss_row_clause (uint32_t row, int forced, int forced_var) {
  std::vector<int> &vs = gauss->rowbuf;
  vs.clear ();
  gauss->mat[row].for_each_set_bit (
      [&] (uint32_t col) { vs.push_back (gauss->col_to_var[col]); });
  if (vs.size () < 2)
    return nullptr; // degenerate row we cannot explain (see gauss.hpp / round)
  return gauss_build_clause (vs, forced, forced_var);
}

// Evaluate a single row against the current column assignment (propGauss).
GaussRet Internal::gauss_prop_row (uint32_t row, uint32_t &new_resp_var,
                                   int &ret_lit_prop) {
  Gauss::PackedRow pr = gauss->mat[row];
  Gauss::PackedRow cols_vals = gauss->aux[0];
  Gauss::PackedRow cols_unset = gauss->aux[1];
  Gauss::PackedRow tmp_col = gauss->aux[2];
  Gauss::PackedRow tmp_col2 = gauss->aux[3];

  const uint32_t pop = tmp_col.set_and_popcnt (pr, cols_unset);
  if (pop >= 2) {
    // Find an unassigned non-basic variable to watch.
    new_resp_var = GaussMatrix::no_var;
    tmp_col.for_each_set_bit ([&] (uint32_t col) {
      if (new_resp_var != GaussMatrix::no_var)
        return;
      const int v = gauss->col_to_var[col];
      if (!gauss->var_has_resp_row[v])
        new_resp_var = v;
    });
    assert (new_resp_var != GaussMatrix::no_var);
    return GaussRet::new_watch;
  }

  const uint32_t pop_t =
      tmp_col2.set_and_popcnt (pr, cols_vals) + (uint32_t) (pr.rhs () & 1);
  if (pop == 1) {
    const int col = tmp_col.first_set_bit ();
    const int v = gauss->col_to_var[col];
    ret_lit_prop = (pop_t & 1) ? v : -v;
    return GaussRet::prop;
  }
  assert (pop == 0);
  return (pop_t & 1) ? GaussRet::confl : GaussRet::satisfied;
}

// Process row 'row' (watched on variable 'p') after 'p' was assigned.
// Returns 0 = keep the watch, 1 = drop it (re-watched elsewhere),
// 2 = conflict (keep the watch).  May schedule a column elimination.
int Internal::gauss_find_truths (uint32_t row, uint32_t p,
                                 uint32_t &new_resp_var, uint32_t &new_resp_row,
                                 bool &do_eliminate) {
  if (gauss->satisfied_xors[row])
    return 0;

  bool was_resp = false;
  if (gauss->var_has_resp_row[p]) {
    was_resp = true;
    gauss->var_has_resp_row[gauss->row_to_var_non_resp[row]] = 1;
    gauss->var_has_resp_row[p] = 0;
  }

  uint32_t nrv = GaussMatrix::no_var;
  int ret_lit = 0;
  const GaussRet ret = gauss_prop_row (row, nrv, ret_lit);

  switch (ret) {
  case GaussRet::confl: {
    Clause *c = gauss_row_clause (row, 0, 0);
    if (was_resp) {
      gauss->var_has_resp_row[gauss->row_to_var_non_resp[row]] = 0;
      gauss->var_has_resp_row[p] = 1;
    }
    if (!c)
      return 0; // degenerate: leave to the model-check gate
    conflict = c;
    gauss->conflicts++;
    return 2;
  }
  case GaussRet::prop: {
    const int fv = (ret_lit < 0) ? -ret_lit : ret_lit;
    Clause *reason = gauss_row_clause (row, ret_lit, fv);
    if (reason) {
      search_assign_driving (ret_lit, reason);
      gauss_update_cols (fv);
      gauss->propagations++;
    }
    gauss->satisfied_xors[row] = 1;
    if (was_resp) {
      gauss->var_has_resp_row[gauss->row_to_var_non_resp[row]] = 0;
      gauss->var_has_resp_row[p] = 1;
    }
    return 0;
  }
  case GaussRet::new_watch: {
    if (was_resp) {
      gauss->gwatch[nrv].clear (); // its column becomes basic -> rebuilt below
      gauss->gwatch[nrv].push_back (row);
      gauss->var_has_resp_row[gauss->row_to_var_non_resp[row]] = 0;
      gauss->var_has_resp_row[nrv] = 1;
      new_resp_var = nrv;
      new_resp_row = row;
      do_eliminate = true;
      return 1;
    }
    gauss->gwatch[nrv].push_back (row);
    gauss->row_to_var_non_resp[row] = nrv;
    return 1;
  }
  case GaussRet::satisfied:
    gauss->satisfied_xors[row] = 1;
    if (was_resp) {
      gauss->var_has_resp_row[gauss->row_to_var_non_resp[row]] = 0;
      gauss->var_has_resp_row[p] = 1;
    }
    return 0;
  }
  return 0;
}

static inline void gauss_remove_watch (std::vector<uint32_t> &ws,
                                       uint32_t row) {
  for (size_t k = 0; k < ws.size (); k++)
    if (ws[k] == row) {
      ws[k] = ws.back ();
      ws.pop_back ();
      return;
    }
}

// Make 'new_resp_var' the basic variable of 'new_resp_row' by XORing that row
// into every other row that has a 1 in the new basic column, fixing up watches
// (and propagating/conflicting) for rows whose non-basic watch is eliminated.
// Returns false if a conflict was installed.
bool Internal::gauss_eliminate_col (uint32_t p, uint32_t new_resp_var,
                                    uint32_t new_resp_row) {
  gauss->eliminations++;
  const int new_resp_col = gauss->var_to_col[new_resp_var];
  bool conflicted = false;

  for (uint32_t row = 0; row < gauss->num_rows; row++) {
    if (row == new_resp_row)
      continue;
    if (!gauss->mat[row][new_resp_col])
      continue;

    const uint32_t orig_non_resp = gauss->row_to_var_non_resp[row];
    const int orig_non_resp_col = gauss->var_to_col[orig_non_resp];

    gauss->mat[row].xor_in (gauss->mat[new_resp_row]);

    if (gauss->mat[row][orig_non_resp_col])
      continue; // still watchable on the same non-basic variable

    if (orig_non_resp != new_resp_var)
      gauss_remove_watch (gauss->gwatch[orig_non_resp], row);

    if (conflicted) {
      gauss->gwatch[p].push_back (row);
      gauss->row_to_var_non_resp[row] = p;
      continue;
    }

    uint32_t nnrv = GaussMatrix::no_var;
    int ret_lit = 0;
    const GaussRet ret = gauss_prop_row (row, nnrv, ret_lit);
    switch (ret) {
    case GaussRet::confl: {
      gauss->gwatch[p].push_back (row);
      gauss->row_to_var_non_resp[row] = p;
      Clause *c = gauss_row_clause (row, 0, 0);
      if (c) {
        conflict = c;
        gauss->conflicts++;
        conflicted = true;
      }
      break;
    }
    case GaussRet::prop: {
      gauss->gwatch[p].push_back (row);
      gauss->row_to_var_non_resp[row] = p;
      const int fv = (ret_lit < 0) ? -ret_lit : ret_lit;
      Clause *reason = gauss_row_clause (row, ret_lit, fv);
      if (reason) {
        search_assign_driving (ret_lit, reason);
        gauss_update_cols (fv);
        gauss->propagations++;
      }
      gauss->satisfied_xors[row] = 1;
      break;
    }
    case GaussRet::new_watch:
      gauss->gwatch[nnrv].push_back (row);
      gauss->row_to_var_non_resp[row] = nnrv;
      break;
    case GaussRet::satisfied:
      gauss->gwatch[p].push_back (row);
      gauss->row_to_var_non_resp[row] = p;
      gauss->satisfied_xors[row] = 1;
      break;
    }
  }
  return !conflicted;
}

// Driver: feed newly assigned trail literals to the matrix.  For each, walk its
// gauss watch list (find_truths, with watch-list compaction) and then, if a
// basic variable moved, run the column elimination.
bool Internal::gauss_round () {
  if (!gauss || !gauss->num_rows)
    return false;

  const size_t tsz = trail.size ();
  // 'gauss_notify_backtrack' already invalidated the column bitsets on any
  // backtrack; here we only clamp the processed-trail pointer.
  if (gauss->qhead > tsz)
    gauss->qhead = tsz;
  if (gauss->qhead >= tsz && !gauss->cancelled_since_val_update)
    return false;

  gauss->rounds++;
  gauss_update_cols (-1);

  const int64_t props_before = gauss->propagations;

  while (!conflict && gauss->qhead < trail.size ()) {
    const int lit = trail[gauss->qhead++];
    const uint32_t p = (lit < 0) ? (uint32_t) -lit : (uint32_t) lit;
    if ((size_t) p >= gauss->var_to_col.size () || gauss->var_to_col[p] < 0)
      continue;

    std::vector<uint32_t> &ws = gauss->gwatch[p];
    uint32_t new_resp_var = GaussMatrix::no_var;
    uint32_t new_resp_row = GaussMatrix::no_var;
    bool do_eliminate = false;
    bool row_confl = false;

    size_t i = 0, j = 0;
    for (; i < ws.size (); i++) {
      const uint32_t row = ws[i];
      const int r =
          gauss_find_truths (row, p, new_resp_var, new_resp_row, do_eliminate);
      if (r == 0)
        ws[j++] = ws[i];
      else if (r == 2) {
        ws[j++] = ws[i];
        i++;
        row_confl = true;
        break;
      }
      // r == 1: watch dropped (re-watched elsewhere)
    }
    for (; i < ws.size (); i++)
      ws[j++] = ws[i];
    ws.resize (j);

    if (row_confl)
      break;
    if (do_eliminate && !gauss_eliminate_col (p, new_resp_var, new_resp_row))
      break; // conflict installed
  }

  return conflict || gauss->propagations != props_before;
}

// --------------------------------------------------------------------------
// Reason/conflict clause builder shared by the engine.  Fills the clause with
// 'forced' (the propagated literal, or 0 for a conflict) followed by the
// currently-false literal of every other variable; orders the two watched
// literals so the two-watched invariant holds.  Returns a redundant clause
// (collectable once it is no longer a reason).
// --------------------------------------------------------------------------
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
    assert (s);
    c.push_back ((s > 0) ? -v : v);
  }
  if (c.size () < 2)
    return nullptr;

  const size_t first_false = forced ? 1 : 0;
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

  clause = c;
  const int glue = (int) clause.size ();
  Clause *res = new_clause (true, glue);
  clause.clear ();
  watch_clause (res);
  return res;
}

// Backstop: verify a candidate complete assignment satisfies every original
// active XOR row.  If one is falsified, install it as a conflict so search
// continues; this guarantees correct models regardless of the incremental
// engine's watch bookkeeping.
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
      Clause *c = gauss_build_clause (vars, 0, 0);
      if (c) {
        conflict = c;
        gauss->conflicts++;
        return false;
      }
    }
  }
  return true;
}

} // namespace CaDiCaL
