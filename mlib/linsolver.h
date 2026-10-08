// Optional linear solver selection, independent of FAST's matrix storage.
#ifndef FAST_LINSOLVER_H
#define FAST_LINSOLVER_H

#include <istream>
#include <string>

class cSysMatrix;
class cVector;

namespace LinearSolver
{
  void ReadBackend(std::istream &input);
  void ReadOptions(std::istream &input);
  void ReadOptionsFile(std::istream &input);
  void Initialize(const std::string &model);
  bool UsesPetsc();
  int Solve(cSysMatrix &matrix, cVector &rhs, cVector &solution,
            int maxiter, double tolerance);

  // Finalizes PETSc on both success and exception paths, only if we own it.
  class Runtime
  {
  public:
    Runtime(int argc, char **argv);
    ~Runtime();
    Runtime(const Runtime &) = delete;
    Runtime &operator=(const Runtime &) = delete;
  };
}

#endif
