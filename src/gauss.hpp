#ifndef CADICAL_GAUSS_HPP
#define CADICAL_GAUSS_HPP

#include <cstdint>
#include <limits>
#include <vector>

#include "gauss/packed_matrix.hpp"

namespace CaDiCaL {

// State of the incremental Gauss-Jordan XOR engine (a port of CryptoMiniSat's
// EGaussian "watched-column" scheme onto CaDiCaL).
//
// The active XOR rows are kept in 'mat' in reduced row echelon form.  Each row
// owns a unique "responsible" (basic) variable -- the only row with a 1 in that
// column -- recorded by 'var_has_resp_row'; it also watches one non-basic
// variable, 'row_to_var_non_resp[row]'.  'gwatch[var]' lists the rows watching
// 'var'.  When a watched variable is assigned, only that row is examined: it
// either finds a fresh non-basic variable to watch, propagates its last
// unassigned variable, or conflicts.  When a *basic* variable is assigned the
// responsibility moves to the freshly watched variable and 'eliminate_col'
// XORs that row into the others sharing the new basic column (the only costly
// step, amortised over the search path).  The matrix elimination is a valid
// GF(2) row operation, so it is never undone on backtracking; only the
// assignment bitsets 'cols_vals'/'cols_unset' and the per-row 'satisfied_xors'
// flags are refreshed (see 'canceling').
struct GaussMatrix {
  Gauss::PackedMatrix mat; // active rows, num_rows x num_cols, kept in RREF
  uint32_t num_rows = 0;
  uint32_t num_cols = 0;

  std::vector<int> col_to_var;   // column -> internal variable index
  std::vector<int> var_to_col;   // internal var -> column (-1 if none)

  std::vector<char> var_has_resp_row;      // per internal var: is it basic?
  std::vector<uint32_t> row_to_var_non_resp; // per row: its watched non-basic var
  std::vector<char> satisfied_xors;        // per row

  // Immutable snapshot of the original active rows (variables + parity), used
  // by the 'gauss_check_model' correctness gate.
  std::vector<std::vector<int>> row_vars;
  std::vector<signed char> row_rhs;

  // gwatch[v] = rows that watch internal variable v.
  std::vector<std::vector<uint32_t>> gwatch;

  // aux rows: 0=cols_vals (bit set => var TRUE), 1=cols_unset (bit set =>
  // var UNASSIGNED), 2=tmp_col, 3=tmp_col2.
  Gauss::PackedMatrix aux;

  bool cancelled_since_val_update = true;
  size_t last_val_update = 0;

  std::vector<int> reason; // clause buffer reused by gauss_build_clause
  std::vector<int> rowbuf; // reused by gauss_row_clause to gather a row's vars

  // Driver bookkeeping: trail position consumed and last seen trail size.
  size_t qhead = 0;
  size_t last_trail = 0;

  // Number of input XOR constraints ('external->xors') this matrix was built
  // from.  Used to detect incrementally-added XOR clauses so the matrix is
  // rebuilt on the next solve instead of silently ignoring the new rows.
  size_t num_input_xors = 0;

  // Statistics.
  int64_t rounds = 0;
  int64_t propagations = 0;
  int64_t conflicts = 0;
  int64_t eliminations = 0;

  static constexpr uint32_t no_var = std::numeric_limits<uint32_t>::max ();
};

// Result of evaluating a single row under the current assignment.
enum class GaussRet { confl, prop, new_watch, satisfied };

} // namespace CaDiCaL

#endif // CADICAL_GAUSS_HPP
