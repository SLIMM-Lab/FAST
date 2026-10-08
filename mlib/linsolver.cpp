#include "linsolver.h"
#include "sysmat.h"
#include "vec.h"

#include <cmath>
#include <stdexcept>
#include <vector>

#ifdef FAST_WITH_PETSC
#include <petscksp.h>
#if defined(PETSC_USE_COMPLEX) || !defined(PETSC_USE_REAL_DOUBLE)
#error "FAST requires PETSc with real, double-precision scalars"
#endif
#endif

namespace
{
  bool usePetsc = false;
  std::string options, optionsFile;
  int argumentCount = 0;
  char **arguments = nullptr;
#ifdef FAST_WITH_PETSC
  bool ownsPetsc = false;
#endif

  std::string ReadQuoted(std::istream &input)
  {
    char quote;
    std::string value;
    if (!(input >> quote) || quote != '\'' ||
        !std::getline(input, value, '\'') || input.eof() ||
        value.find("\n%") != std::string::npos)
      throw std::runtime_error("Linear solver settings require a single-quoted string");
    return value;
  }

#ifdef FAST_WITH_PETSC
  void Check(PetscErrorCode error)
  {
    if (error)
    {
      const char *message = nullptr;
      PetscErrorMessage(error, &message, nullptr);
      throw std::runtime_error("PETSc error " + std::to_string(error) + ": " +
                               (message ? message : "check solver options"));
    }
  }

  struct SolveObjects
  {
    Mat matrix = nullptr;
    Vec rhs = nullptr, solution = nullptr;
    KSP solver = nullptr;
    ~SolveObjects()
    {
      KSPDestroy(&solver);
      VecDestroy(&solution);
      VecDestroy(&rhs);
      MatDestroy(&matrix);
    }
  };
#endif
}

namespace LinearSolver
{
  Runtime::Runtime(int argc, char **argv)
  {
    argumentCount = argc;
    arguments = argv;
  }

  Runtime::~Runtime()
  {
#ifdef FAST_WITH_PETSC
    if (ownsPetsc) PetscFinalize();
#endif
  }

  void ReadBackend(std::istream &input)
  {
    const std::string backend = ReadQuoted(input);
    if (backend != "native" && backend != "petsc")
      throw std::runtime_error("Unknown linear solver '" + backend +
                               "': choose 'native' or 'petsc'");
    usePetsc = backend == "petsc";
  }

  void ReadOptions(std::istream &input) { options = ReadQuoted(input); }
  void ReadOptionsFile(std::istream &input) { optionsFile = ReadQuoted(input); }
  bool UsesPetsc() { return usePetsc; }

  void Initialize(const std::string &model)
  {
    if (!usePetsc)
    {
      if (!options.empty() || !optionsFile.empty())
        throw std::runtime_error("PETSc options require LINEAR.SOLVER 'petsc'");
      return;
    }
#ifndef FAST_WITH_PETSC
    throw std::runtime_error("PETSc support is disabled; rebuild with USE_PETSC=1");
#else
    PetscBool initialized;
    Check(PetscInitialized(&initialized));
    if (!initialized)
    {
      // Initialize without argv; apply the command line last for precedence.
      Check(PetscInitializeNoArguments());
      ownsPetsc = true;
      Check(PetscPushErrorHandler(PetscReturnErrorHandler, nullptr));
    }
    int size;
    Check(MPI_Comm_size(PETSC_COMM_WORLD, &size));
    if (size != 1)
      throw std::runtime_error("FAST assembly is serial; run PETSc with one MPI rank");
    if (!optionsFile.empty())
    {
      std::string path = optionsFile;
      const auto slash = model.find_last_of("/\\");
      if (path[0] != '/' && slash != std::string::npos)
        path = model.substr(0, slash + 1) + path;
      Check(PetscOptionsInsertFile(PETSC_COMM_SELF, nullptr, path.c_str(), PETSC_TRUE));
    }
    Check(PetscOptionsInsertString(nullptr, options.c_str()));
    Check(PetscOptionsInsert(nullptr, &argumentCount, &arguments, nullptr));
#endif
  }

  int Solve(cSysMatrix &matrix, cVector &rhs, cVector &solution,
            int maxiter, double tolerance)
  {
#ifndef FAST_WITH_PETSC
    throw std::runtime_error("PETSc support is disabled; rebuild with USE_PETSC=1");
#else
    const int n = matrix.Dim();
    if (rhs.Dim() != n || solution.Dim() != n)
      throw std::runtime_error("PETSc linear system dimensions do not match");
    if (n == 0) return 1;

    SolveObjects data;
    std::vector<PetscInt> counts(n, 0);
    matrix.VisitEntries([&](int i, int, double value) {
      if (value != 0.0) ++counts[i];
    });
    Check(MatCreateSeqAIJ(PETSC_COMM_SELF, n, n, 0, counts.data(), &data.matrix));
    matrix.VisitEntries([&](int i, int j, double value) {
      if (value != 0.0) Check(MatSetValue(data.matrix, i, j, value, INSERT_VALUES));
    });
    Check(MatAssemblyBegin(data.matrix, MAT_FINAL_ASSEMBLY));
    Check(MatAssemblyEnd(data.matrix, MAT_FINAL_ASSEMBLY));
    Check(MatCreateVecs(data.matrix, &data.solution, &data.rhs));
    Check(VecSet(data.solution, 0.0));
    PetscScalar *values;
    Check(VecGetArray(data.rhs, &values));
    for (int i = 0; i < n; ++i) values[i] = rhs[i];
    Check(VecRestoreArray(data.rhs, &values));

    Check(KSPCreate(PETSC_COMM_SELF, &data.solver));
    Check(KSPSetOperators(data.solver, data.matrix, data.matrix));
    Check(KSPSetType(data.solver, KSPGMRES));
    PC pc;
    Check(KSPGetPC(data.solver, &pc));
    Check(PCSetType(pc, PCILU));
    Check(KSPSetTolerances(data.solver, tolerance, PETSC_DEFAULT, PETSC_DEFAULT,
                           maxiter > 0 ? maxiter : PETSC_DEFAULT));
    Check(KSPSetFromOptions(data.solver));
    // Always copy the RHS before solving, including Solve(v, v) in MPC paths.
    Check(KSPSolve(data.solver, data.rhs, data.solution));
    KSPConvergedReason reason;
    Check(KSPGetConvergedReason(data.solver, &reason));
    if (reason <= 0)
      throw std::runtime_error("PETSc linear solver did not converge (reason " +
                               std::to_string(reason) + ")");
    const PetscScalar *answer;
    Check(VecGetArrayRead(data.solution, &answer));
    bool finite = true;
    for (int i = 0; i < n; ++i) finite = finite && std::isfinite(answer[i]);
    if (!finite)
    {
      Check(VecRestoreArrayRead(data.solution, &answer));
      throw std::runtime_error("PETSc returned a non-finite solution");
    }
    for (int i = 0; i < n; ++i) solution[i] = answer[i];
    Check(VecRestoreArrayRead(data.solution, &answer));
    return 1;
#endif
  }
}
