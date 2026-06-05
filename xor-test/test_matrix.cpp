// Standalone unit test for the ported GF(2) matrix core (PackedRow /
// PackedMatrix).  It builds many random linear systems over GF(2), runs
// Gauss-Jordan elimination using the ported arithmetic, and cross-checks the
// SAT/UNSAT (consistency) verdict against a completely independent reference
// eliminator.  No CaDiCaL internals are involved.
//
// Build:  g++ -O2 -std=c++17 -I../src test_matrix.cpp -o test_matrix

#include "gauss/packed_matrix.hpp"

#include <cstdint>
#include <iostream>
#include <random>
#include <vector>

using CaDiCaL::Gauss::PackedMatrix;
using CaDiCaL::Gauss::PackedRow;

// ---- Reference GF(2) eliminator (independent, simple, vector<uint64_t>) ----
struct RefRow {
  std::vector<uint64_t> bits;
  int rhs;
};

static bool ref_consistent (std::vector<RefRow> rows, int ncols) {
  int words = (ncols + 63) / 64;
  for (auto &r : rows)
    r.bits.resize (words, 0);
  auto getbit = [&] (const RefRow &r, int c) {
    return (r.bits[c / 64] >> (c % 64)) & 1ULL;
  };
  auto xorinto = [&] (RefRow &dst, const RefRow &src) {
    for (int i = 0; i < words; i++)
      dst.bits[i] ^= src.bits[i];
    dst.rhs ^= src.rhs;
  };
  int pivot = 0;
  for (int c = 0; c < ncols && pivot < (int) rows.size (); c++) {
    int sel = -1;
    for (int r = pivot; r < (int) rows.size (); r++)
      if (getbit (rows[r], c)) {
        sel = r;
        break;
      }
    if (sel < 0)
      continue;
    std::swap (rows[pivot], rows[sel]);
    for (int r = 0; r < (int) rows.size (); r++)
      if (r != pivot && getbit (rows[r], c))
        xorinto (rows[r], rows[pivot]);
    pivot++;
  }
  for (auto &r : rows) {
    bool zero = true;
    for (int i = 0; i < words; i++)
      if (r.bits[i]) {
        zero = false;
        break;
      }
    if (zero && r.rhs)
      return false; // 0 == 1  -> inconsistent
  }
  return true;
}

// ---- Eliminator using the ported PackedMatrix ----
static bool packed_consistent (const std::vector<RefRow> &rows, int ncols) {
  const uint32_t nrows = rows.size ();
  if (nrows == 0)
    return true;
  PackedMatrix m;
  m.resize (nrows, ncols);
  for (uint32_t r = 0; r < nrows; r++) {
    PackedRow pr = m[r];
    pr.setZero ();
    pr.rhs () = rows[r].rhs;
    for (int c = 0; c < ncols; c++)
      if ((rows[r].bits[c / 64] >> (c % 64)) & 1ULL)
        pr.setBit (c);
  }
  uint32_t pivot = 0;
  for (int c = 0; c < ncols && pivot < nrows; c++) {
    int sel = -1;
    for (uint32_t r = pivot; r < nrows; r++)
      if (m[r][c]) {
        sel = r;
        break;
      }
    if (sel < 0)
      continue;
    if ((uint32_t) sel != pivot)
      m[pivot].swapBoth (m[sel]);
    for (uint32_t r = 0; r < nrows; r++)
      if (r != pivot && m[r][c])
        m[r].xor_in (m[pivot]);
    pivot++;
  }
  for (uint32_t r = 0; r < nrows; r++)
    if (m[r].isZero () && m[r].rhs ())
      return false;
  return true;
}

int main () {
  std::mt19937 rng (12345);
  int mismatches = 0, total = 0, sat = 0, unsat = 0;
  for (int trial = 0; trial < 20000; trial++) {
    int ncols = 1 + rng () % 200;      // variables
    int nrows = 1 + rng () % 120;      // xor rows
    std::vector<RefRow> rows (nrows);
    int words = (ncols + 63) / 64;
    // 50%: build a guaranteed-consistent system from a planted solution;
    // 50%: fully random (may be inconsistent).
    bool planted = (rng () & 1);
    std::vector<int> sol (ncols);
    for (int i = 0; i < ncols; i++)
      sol[i] = rng () & 1;
    for (int r = 0; r < nrows; r++) {
      rows[r].bits.assign (words, 0);
      int rhs = 0;
      for (int c = 0; c < ncols; c++)
        if (rng () % 4 == 0) { // sparse-ish
          rows[r].bits[c / 64] |= (1ULL << (c % 64));
          rhs ^= sol[c];
        }
      rows[r].rhs = planted ? rhs : (int) (rng () & 1);
    }
    bool a = ref_consistent (rows, ncols);
    bool b = packed_consistent (rows, ncols);
    total++;
    if (a)
      sat++;
    else
      unsat++;
    if (a != b) {
      mismatches++;
      if (mismatches <= 5)
        std::cerr << "MISMATCH trial=" << trial << " ref=" << a
                  << " packed=" << b << " ncols=" << ncols
                  << " nrows=" << nrows << "\n";
    }
    if (planted && !a) {
      // planted systems must be consistent
      std::cerr << "BUG: planted system reported inconsistent (trial "
                << trial << ")\n";
      mismatches++;
    }
  }
  std::cout << "trials=" << total << " consistent=" << sat
            << " inconsistent=" << unsat << " mismatches=" << mismatches
            << "\n";
  return mismatches ? 1 : 0;
}
