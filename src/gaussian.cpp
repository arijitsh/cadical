#include "gaussian.hpp"

#include <algorithm>
#include <iostream>

namespace CaDiCaL {

void Gaussian::add_clause(const Xor &clause) {
  equations.push_back(clause);
}

std::optional<std::vector<bool>> Gaussian::eliminate() const {
  std::cout << "[gaussian] running elimination on " << equations.size()
            << " equations" << std::endl;
  // Determine the maximum variable index
  int max_var = 0;
  for (const auto &eq : equations)
    for (int v : eq.vars)
      if (v > max_var)
        max_var = v;

  std::vector<std::vector<int>> mat;
  std::vector<int> rhs;
  for (const auto &eq : equations) {
    std::vector<int> row(max_var, 0);
    for (int v : eq.vars)
      if (v >= 1 && v <= max_var)
        row[v - 1] ^= 1;
    mat.push_back(row);
    rhs.push_back(eq.rhs ? 1 : 0);
  }

  int rows = mat.size();
  int cols = max_var;
  int r = 0;
  for (int c = 0; c < cols && r < rows; ++c) {
    int pivot = -1;
    for (int i = r; i < rows; ++i) {
      if (mat[i][c]) {
        pivot = i;
        break;
      }
    }
    if (pivot == -1)
      continue;
    if (pivot != r) {
      std::swap(mat[pivot], mat[r]);
      std::swap(rhs[pivot], rhs[r]);
    }
    for (int i = 0; i < rows; ++i) {
      if (i != r && mat[i][c]) {
        for (int j = c; j < cols; ++j)
          mat[i][j] ^= mat[r][j];
        rhs[i] ^= rhs[r];
      }
    }
    ++r;
  }

  // Check for inconsistency
  for (int i = r; i < rows; ++i) {
    bool all_zero = true;
    for (int c = 0; c < cols; ++c) {
      if (mat[i][c]) {
        all_zero = false;
        break;
      }
    }
    if (all_zero && rhs[i])
      return {};
  }

  std::vector<bool> solution(max_var, false);
  for (int i = r - 1; i >= 0; --i) {
    int pivot_col = -1;
    for (int c = 0; c < cols; ++c)
      if (mat[i][c]) {
        pivot_col = c;
        break;
      }
    if (pivot_col == -1)
      continue;
    int sum = rhs[i];
    for (int c = pivot_col + 1; c < cols; ++c)
      if (mat[i][c] && solution[c])
        sum ^= 1;
    solution[pivot_col] = sum;
  }

  return solution;
}

} // namespace CaDiCaL

