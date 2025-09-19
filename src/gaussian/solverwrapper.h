#ifndef GAUSSWRAPPER_H
#define GAUSSWRAPPER_H

#include <vector>

namespace CaDiCaL {
struct Internal;
}

namespace CMSat {

// Very small Gaussian elimination wrapper used for tests.  Full
// functionality is implemented in 'solverwrapper.cpp'.
class Solver {
  CaDiCaL::Internal *internal;

  struct XorClause {
    std::vector<unsigned> vars;
    bool rhs;
  };

  std::vector<XorClause> xors;

public:
  explicit Solver(CaDiCaL::Internal *i = nullptr);
  unsigned gqhead = 0;

  bool init_gauss();

  enum class gauss_ret { g_cont, g_nothing, g_false };

  gauss_ret run_gauss();
  void cancel();
  void run_top_level_gauss();
  bool add_xor_clause(const std::vector<unsigned> &vars, bool rhs);
};

} // namespace CMSat

#endif // GAUSSWRAPPER_H
