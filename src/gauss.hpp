#ifndef CADICAL_GAUSS_HPP
#define CADICAL_GAUSS_HPP

#include <cstdint>
#include <vector>

#include "gauss/packed_matrix.hpp"

namespace CaDiCaL {

// State of the Gauss-Jordan XOR engine.
//
// The input XOR constraints are reduced once to row echelon form at the start
// of solving; unit / empty rows are dealt with at the root level and the
// remaining multi-variable rows are kept (their bit pattern in 'base', their
// variables in 'row_vars', their target parity in 'row_rhs').
//
// During search 'Internal::gauss_round' performs *full* Gaussian elimination
// of the rows over the currently unassigned columns (a fresh copy in 'work').
// Reducing over the unassigned columns exposes every variable that the linear
// system forces under the current partial assignment, and every linear
// combination of rows that has become contradictory -- i.e. it has the full
// propagation strength of Gaussian reasoning, not just per-row XOR-unit
// propagation.  The computation is stateless (it reads 'val' and rebuilds
// 'work' each time), so no backtrack bookkeeping is needed; it is only redone
// when an assignment to a matrix variable has changed since the last round.
//
// Propagations and conflicts are explained with a real CaDiCaL clause (one
// clause of the responsible row's CNF expansion, selected by the assignment),
// so conflict analysis / backtracking are reused unchanged.
struct GaussMatrix {
  // Active rows after the initial reduction: bit matrix and metadata.
  Gauss::PackedMatrix base; // num_rows x num_cols, initial reduced rows
  uint32_t num_rows = 0;
  uint32_t num_cols = 0;

  std::vector<int> col_to_var;        // column -> internal variable index
  std::vector<int> var_to_col;        // internal var -> column (-1 if none)
  std::vector<std::vector<int>> row_vars; // variables of each active row
  std::vector<signed char> row_rhs;       // target parity of each active row

  // Scratch reused across rounds.
  Gauss::PackedMatrix work; // working copy that gets eliminated
  Gauss::PackedMatrix aux;  // two rows: aux[0]=cols_vals, aux[1]=cols_unset
  std::vector<int> reason;  // clause buffer

  // Change detection: trail position consumed and last seen trail size.
  size_t qhead = 0;
  size_t last_trail = 0;
  bool need_round = true;

  // Statistics.
  int64_t rounds = 0;
  int64_t propagations = 0;
  int64_t conflicts = 0;
};

} // namespace CaDiCaL

#endif // CADICAL_GAUSS_HPP
