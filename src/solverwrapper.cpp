#include "gaussian/solverwrapper.h"
#include "internal.hpp"

#include <algorithm>
#include <cstdint>

namespace CMSat {

Solver::Solver(CaDiCaL::Internal *i) : internal(i) {}

bool Solver::init_gauss() { return true; }

// Simple Gaussian elimination over GF(2) that detects inconsistencies
// between the current assignment in 'internal' and the stored XOR system.
Solver::gauss_ret Solver::run_gauss() {
  if (xors.empty())
    return gauss_ret::g_nothing;

  // determine matrix width
  unsigned max_var = 0;
  for (const auto &xr : xors)
    for (unsigned v : xr.vars)
      max_var = std::max(max_var, v);

  const unsigned nvars = max_var + 1;
  const unsigned chunks = (nvars + 63) / 64;

  struct Row {
    std::vector<uint64_t> bits;
    bool rhs;
  };

  std::vector<Row> rows;
  rows.reserve(xors.size());

  // Substitute assigned variables and build matrix of unassigned ones.
  for (const auto &xr : xors) {
    bool parity = xr.rhs;
    Row row;
    row.bits.assign(chunks, 0);

    for (unsigned v : xr.vars) {
      const int val = internal->val((int)v + 1);
      if (!val) {
        row.bits[v >> 6] ^= (uint64_t)1 << (v & 63);
      } else if (val > 0) {
        parity = !parity;
      }
    }

    row.rhs = parity;

    bool empty = true;
    for (auto b : row.bits)
      if (b) {
        empty = false;
        break;
      }

    if (empty) {
      if (row.rhs)
        return gauss_ret::g_false; // 0 = 1 -> conflict
      continue;                     // redundant 0 = 0 row
    }

    rows.push_back(std::move(row));
  }

  if (rows.empty())
    return gauss_ret::g_nothing;

  // Perform elimination to detect 0 = 1 rows.
  size_t rank = 0;
  for (unsigned col = 0; col < nvars && rank < rows.size(); ++col) {
    const unsigned chunk = col >> 6;
    const uint64_t mask = (uint64_t)1 << (col & 63);

    size_t pivot = rank;
    while (pivot < rows.size() && !(rows[pivot].bits[chunk] & mask))
      ++pivot;
    if (pivot == rows.size())
      continue;
    std::swap(rows[rank], rows[pivot]);

    for (size_t i = 0; i < rows.size(); ++i) {
      if (i == rank)
        continue;
      if (rows[i].bits[chunk] & mask) {
        for (unsigned k = 0; k < chunks; ++k)
          rows[i].bits[k] ^= rows[rank].bits[k];
        rows[i].rhs ^= rows[rank].rhs;
      }
    }
    ++rank;
  }

  for (const auto &row : rows) {
    bool empty = true;
    for (auto b : row.bits)
      if (b) {
        empty = false;
        break;
      }
    if (empty && row.rhs)
      return gauss_ret::g_false;
  }

  return gauss_ret::g_nothing;
}

void Solver::cancel() {}

void Solver::run_top_level_gauss() {}

bool Solver::add_xor_clause(const std::vector<unsigned> &vars, bool rhs) {
  xors.push_back({vars, rhs});
  return true;
}

} // namespace CMSat
