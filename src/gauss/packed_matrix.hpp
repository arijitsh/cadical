/******************************************
Gauss-Jordan matrix core for CaDiCaL, ported from CryptoMiniSat.

Copyright (C) 2009-2020 Authors of CryptoMiniSat, see AUTHORS file
Copyright (c) 2012  Cheng-Shen Han
Copyright (c) 2012  Jie-Hong Roland Jiang

For more information, see "When Boolean Satisfiability Meets Gaussian
Elimination in a Simplex Way." by Cheng-Shen Han and Jie-Hong Roland Jiang
in CAV (Computer Aided Verification), 2012: 410-426.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction ... (MIT license, see CryptoMiniSat).
***********************************************/

// This header contains ONLY the self-contained GF(2) bit-matrix arithmetic
// (PackedRow / PackedMatrix).  The solver-coupled propagation helpers
// (find_watchVar / propGauss / get_reason) live in the EGaussian port and are
// added in a later phase.  Keeping the arithmetic separate lets it be unit
// tested independently of CaDiCaL's internals.

#ifndef CADICAL_GAUSS_PACKED_MATRIX_HPP
#define CADICAL_GAUSS_PACKED_MATRIX_HPP

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>

namespace CaDiCaL {
namespace Gauss {

// A single row of the GF(2) matrix.  The backing store layout is exactly the
// CryptoMiniSat one: the int64 word immediately BEFORE 'mp' holds the row's
// right-hand side (rhs, 0/1) and words mp[0..size-1] hold the column bits,
// 64 columns per word.  PackedRow is a non-owning view; PackedMatrix owns the
// contiguous allocation.
class PackedRow {
public:
  PackedRow () = delete;

  PackedRow &operator= (const PackedRow &b) {
    // start from -1, because that is where the rhs lives
    for (int i = -1; i < size; i++)
      *(mp + i) = *(b.mp + i);
    return *this;
  }

  PackedRow &operator^= (const PackedRow &b) {
    for (int i = -1; i < size; i++)
      *(mp + i) ^= *(b.mp + i);
    return *this;
  }

  void xor_in (const PackedRow &b) {
    rhs_internal ^= b.rhs_internal;
    for (int i = 0; i < size; i++)
      *(mp + i) ^= *(b.mp + i);
  }

  const int64_t &rhs () const { return rhs_internal; }
  int64_t &rhs () { return rhs_internal; }

  bool isZero () const {
    for (int i = 0; i < size; i++)
      if (mp[i])
        return false;
    return true;
  }

  void setZero () { memset (mp, 0, sizeof (int64_t) * size); }
  void setOne () { memset (mp, 0xff, sizeof (int64_t) * size); }

  void clearBit (const uint32_t i) { mp[i / 64] &= ~(1LL << (i % 64)); }
  void setBit (const uint32_t i) { mp[i / 64] |= (1LL << (i % 64)); }
  // Toggle a bit: building an XOR row by toggling makes a variable that occurs
  // an even number of times (e.g. after equivalence substitution) cancel.
  void flipBit (const uint32_t i) { mp[i / 64] ^= (1LL << (i % 64)); }

  void invert_rhs (const bool b = true) { rhs_internal ^= (int) b; }

  void swapBoth (PackedRow b) {
    int64_t *__restrict mp1 = mp - 1;
    int64_t *__restrict mp2 = b.mp - 1;
    uint32_t i = size + 1;
    while (i != 0) {
      std::swap (*mp1, *mp2);
      mp1++;
      mp2++;
      i--;
    }
  }

  bool operator[] (const uint32_t i) const {
    return (mp[i / 64] >> (i % 64)) & 1;
  }

  uint32_t popcnt () const {
    uint32_t ret = 0;
    for (int i = 0; i < size; i++)
      ret += __builtin_popcountll ((uint64_t) mp[i]);
    return ret;
  }

  uint32_t popcnt_at_least_2 () const {
    uint32_t ret = 0;
    for (int i = 0; i < size && ret < 2; i++)
      ret += __builtin_popcountll ((uint64_t) mp[i]);
    return ret;
  }

  int get_size () const { return size; }

private:
  friend class PackedMatrix;
  friend std::ostream &operator<< (std::ostream &os, const PackedRow &m);

  PackedRow (const uint32_t _size, int64_t *const _mp)
      : mp (_mp + 1), rhs_internal (*_mp), size (_size) {}

  int64_t *__restrict const mp;
  int64_t &rhs_internal;
  const int size;
};

inline std::ostream &operator<< (std::ostream &os, const PackedRow &m) {
  for (int i = 0; i < m.size * 64; i++)
    os << (int) m[i];
  os << " -- rhs: " << m.rhs ();
  return os;
}

// Owns a contiguous block of rows.  Column count is rounded up to a multiple
// of 64; each row occupies (numCols+1) int64 words (the extra leading word is
// the rhs).
class PackedMatrix {
public:
  PackedMatrix () : mp (nullptr), numRows (0), numCols (0) {}

  ~PackedMatrix () { free (mp); }

  void resize (const uint32_t num_rows, uint32_t num_cols) {
    num_cols = num_cols / 64 + (bool) (num_cols % 64);
    if (numRows * (numCols + 1) < (int) num_rows * ((int) num_cols + 1)) {
      size_t size = sizeof (int64_t) * num_rows * (num_cols + 1);
      free (mp);
      int ret = posix_memalign ((void **) &mp, 16, size);
      (void) ret;
    }
    numRows = num_rows;
    numCols = num_cols;
  }

  void resizeNumRows (const uint32_t num_rows) { numRows = num_rows; }

  PackedMatrix &operator= (const PackedMatrix &b) {
    if (numRows * (numCols + 1) < b.numRows * (b.numCols + 1)) {
      size_t size = sizeof (int64_t) * b.numRows * (b.numCols + 1);
      free (mp);
      int ret = posix_memalign ((void **) &mp, 16, size);
      (void) ret;
    }
    numRows = b.numRows;
    numCols = b.numCols;
    memcpy (mp, b.mp, sizeof (int64_t) * numRows * (numCols + 1));
    return *this;
  }

  PackedRow operator[] (const uint32_t i) {
    return PackedRow (numCols, mp + i * (numCols + 1));
  }

  PackedRow operator[] (const uint32_t i) const {
    return PackedRow (numCols, mp + i * (numCols + 1));
  }

  uint32_t getSize () const { return numRows; }
  int get_num_col_words () const { return numCols; }

private:
  int64_t *mp;
  int numRows;
  int numCols;
};

} // namespace Gauss
} // namespace CaDiCaL

#endif // CADICAL_GAUSS_PACKED_MATRIX_HPP
