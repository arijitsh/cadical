#ifndef CADICAL_GAUSS_HPP
#define CADICAL_GAUSS_HPP

#include <cstdint>
#include <vector>

namespace CaDiCaL {

// Search-time state of the Gauss-Jordan XOR engine.
//
// The input XOR constraints are reduced once to row echelon form at the start
// of solving (in 'Internal::init_gauss').  After that reduction every row is a
// linear (GF(2)) combination of the original XORs and is therefore itself a
// valid XOR constraint.  We keep, for every row with at least two variables,
// the list of its internal variable indices and its target parity ('rhs').
//
// During search 'Internal::gauss_round' simply re-evaluates each stored row
// against the current assignment ('val'): a row with exactly one unassigned
// variable propagates it; a fully assigned row with the wrong parity is a
// conflict.  Propagations and conflicts are explained with a real CaDiCaL
// clause (one clause of the row's CNF expansion, selected by the assignment),
// so the existing conflict-analysis / backtracking machinery is reused
// unchanged.  Because the rows are static, no backtrack hook is required.
//
// This is the correct-but-unoptimised engine; the watched-column incremental
// elimination (Han-Jiang / CryptoMiniSat EGaussian) is a later performance
// layer on top of the same validated matrix core.
struct GaussMatrix {
  // For each active row (>= 2 variables): its internal variable indices.
  std::vector<std::vector<int>> row_vars;
  // Target parity (0/1) for each active row, aligned with 'row_vars'.
  std::vector<signed char> row_rhs;

  // Occurrence index: 'var_rows[v]' lists the active rows containing internal
  // variable 'v' (sized max_var+1).  Used to examine, on each new assignment,
  // only the rows that variable participates in instead of all rows.
  std::vector<std::vector<int>> var_rows;

  // Position on the solver trail up to which assignments have already been fed
  // to the engine.  Clamped down on backtrack (the trail shrinks); evaluation
  // is stateless (reads 'val'), so re-assigned variables are simply
  // re-examined.
  size_t qhead = 0;

  // Scratch buffer reused while building reason / conflict clauses.
  std::vector<int> reason;

  // Statistics.
  int64_t propagations = 0;
  int64_t conflicts = 0;
};

} // namespace CaDiCaL

#endif // CADICAL_GAUSS_HPP
