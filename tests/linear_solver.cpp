#include "linsolver.h"
#include "sysmat.h"
#include "vec.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>

static void Require(bool condition, const char *message)
{
  if (!condition) throw std::runtime_error(message);
}

static void CheckSolution(cVector &x)
{
  for (int i = 0; i < 3; ++i)
    Require(std::isfinite(x[i]) && std::abs(x[i] - (i + 1)) < 1.e-8,
            "Incorrect linear solution");
}

int main(int argc, char **argv)
{
  try
  {
    LinearSolver::Runtime runtime(argc, argv);
    const bool petsc = argc > 1 && std::string(argv[1]).find("petsc") == 0;
    const bool spdOnly = argc > 1 && std::string(argv[1]) == "petsc-spd";
    std::istringstream backend(petsc ? "'petsc'" : "'native'");
    LinearSolver::ReadBackend(backend);
    LinearSolver::Initialize("test");
    int profile[] = {0, 0, 1};
    const eSMType types[] = {SYM_SKYLINE, UNSYM_SKYLINE, SYM_SPARSE,
                            UNSYM_SPARSE, EIGEN_UNSYM_SPARSE};
    for (eSMType type : types)
    {
      std::unique_ptr<cSysMatrix> matrix(cSysMatrix::CreateMatrix(type, 3, profile));
      matrix->SetParam(100, 1.e-12);
      cVector rhs(3), solution(3);
      for (int repeat = 0; repeat < 2; ++repeat)
      {
        matrix->Zero();
        // Reassembly must replace, rather than accumulate, previous values.
        const double scale = repeat + 1;
        const double entries[3][3] = {{4,1,0},{1,3,1},{0,1,2}};
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) matrix->Add(i, j, scale*entries[i][j]);
        rhs[0] = 6*scale; rhs[1] = 10*scale; rhs[2] = 8*scale;
        Require(matrix->Solve(rhs, solution) != 0, "Linear solve failed");
        CheckSolution(solution);
        // MPC and eigensolver callers use the same vector for both arguments.
        Require(matrix->Solve(rhs, rhs) != 0, "Aliased RHS solve failed");
        CheckSolution(rhs);
        if (petsc)
          Require(std::abs(matrix->Get(0,0) - 4*scale) < 1.e-12,
                  "PETSc changed the assembled matrix");
      }
      cMatrix constraints(1,3), lagrangian(1,1);
      constraints[0][0] = 1; constraints[0][1] = 0; constraints[0][2] = 0;
      cVector target(1), multipliers(1), correction(3), forces(3);
      target[0] = 0.5;
      Require(matrix->SolveMPC(constraints, target, solution, lagrangian,
                               multipliers, correction, forces) != 0,
              "MPC solve failed");
      Require(std::abs(solution[0] - correction[0] - 0.5) < 1.e-8,
              "MPC constraint was not enforced");
      rhs.Zero();
      Require(matrix->Solve(rhs, solution) != 0, "Zero RHS solve failed");
      Require(solution.Length() < 1.e-12, "Zero RHS gave a nonzero solution");
      std::cout << "PASS: storage " << type << ", reassembly, aliased RHS, MPC and zero RHS\n";
    }

    if (petsc)
    {
      // Actual nonsymmetric values exercise skyline's separate upper storage.
      for (eSMType type : {UNSYM_SKYLINE, UNSYM_SPARSE})
      {
        if (spdOnly) continue; // CG requires a symmetric positive definite operator.
        std::unique_ptr<cSysMatrix> matrix(cSysMatrix::CreateMatrix(type, 3, profile));
        matrix->Zero(); matrix->SetParam(100, 1.e-12);
        const double entries[3][3] = {{4,2,0},{1,3,1},{0,-1,2}};
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j) matrix->Add(i, j, entries[i][j]);
        cVector rhs(3), solution(3);
        rhs[0] = 8; rhs[1] = 10; rhs[2] = 4;
        Require(matrix->Solve(rhs, solution) != 0, "Nonsymmetric solve failed");
        CheckSolution(solution);
      }
      if (!spdOnly) std::cout << "PASS: nonsymmetric matrices\n";

      // Positive original diagonals can still yield a negative Crout pivot.
      int fullProfile[] = {0,0};
      for (eSMType type : {SYM_SKYLINE, UNSYM_SKYLINE})
      {
        std::unique_ptr<cSysMatrix> matrix(cSysMatrix::CreateMatrix(type, 2, fullProfile));
        matrix->Zero();
        matrix->Add(0,0,1); matrix->Add(1,1,1);
        matrix->Add(1,0,2); matrix->Add(0,1,2);
        Require(matrix->GetNgtPvt() == 1, "Incorrect stability pivot count");
        Require(matrix->Get(1,1) == 1, "Pivot diagnostic changed the matrix");
      }
      std::cout << "PASS: stability pivots preserved\n";
    }
    return 0;
  }
  catch (const std::exception &error)
  {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
