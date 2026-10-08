// -------------------------------------------------------------------------
// sysmat.h - Definition of classes to handle square sparse matrices of
//            large linear systems.
// -------------------------------------------------------------------------
// Remarks:
// -------------------------------------------------------------------------
// cSysMatrix => class for large scale square matrices
//   cSymSkylMatrix => symmetric skyline matrices (lower triangular)
//   cUnsymSkylMatrix => general skyline matrices
//   cSprsMatrix => general sparse matrix
//     cSymSprsMatrix => symmetric sparse matrix (lower triangular)
//
// ------------------------------------------------------------------------
// 
// int SolveMPC(cMatrix &P, cVector &b, cVector &uf, cMatrix &A,
//              cVector &lm, cVector &ul, cVector &fl)
//
//   P  - constraint matrix                                          (in)
//   b  - constant term in constraint equation                       (in)
//   A  - Lagrangian matrix                                          (out)
//   lm - Lagrangian multipliers                                     (out)
//   uf - unconstrained displacements                                (out)
//   ul - displacements to enforce the constraints                   (out)
//   fl - forces to enforce the constraints                          (out)
// 
// Enforces the multi-point (equality) constraints [P]{u} = {b} by using
// the Lagrange Multipliers Method.
// For linear problems the final displacements should be computed as
// {u} = {uf} - {ul}.
// For nonlinear problems using the Newton-Rapshon Method the above
// relation is valid, but {u}, {uf} and {ul} are incremental rather than
// total  displacements.
// ----------------------------------------------------------------------

#ifndef _SYSMAT_H
#define _SYSMAT_H

#include <functional>

// -------------------------------------------------------------------------
// SysMat types:
//
typedef enum
{
  SYM_SKYLINE,
  UNSYM_SKYLINE,
  SYM_SPARSE,
  UNSYM_SPARSE,
  EIGEN_UNSYM_SPARSE
} eSMType;

// -------------------------------------------------------------------------
// Forward declarations:
//
class cVector;
class cMatrix;
class cSprsMatrix;

// -------------------------------------------------------------------------
// Abstract base class to handle sparse matrices:
//
class cSysMatrix
{
 friend class cVector;  // Allow cVector access the private data

 protected:
  int dim;        // Matrix dimension (rows/columns)
  int nelm;       // Number of stored elements
  int maxit;      // Max iterations for linear solver
  double tol;     // Solver tolerance

 protected:
  int IsSup(int i, int j) { return(i < j); }
  int IsInf(int i, int j) { return(i > j); }

 public:
  static cSysMatrix *CreateMatrix(eSMType, int, int *);

                   cSysMatrix(int n = 1) { dim = n; nelm = 0;
                                           maxit = n; tol = 1.0e-05; }
  virtual         ~cSysMatrix(void) { }
          int      GetNumElm (void) { return nelm; }
          int      Dim       (void) { return dim; }
          void     SetParam  (int k, double t) { maxit = k; tol = t; }
          int      SolveMPC  (cMatrix&, cVector&, cVector&, cMatrix&, 
                              cVector&, cVector&, cVector&);
          int      SolveMPC  (cSprsMatrix &, cVector&, cVector&, cMatrix&,  
                              cVector&, cVector&, cVector&);
  virtual int      Symmetric (void) { return 0; }
  virtual double **Val       (void) { return 0; }
  virtual void     Zero      (void) { }
  virtual int      Solve     (cVector &, cVector &) { return 0; }
  virtual int      GetSmlPvt (void) { return 0; }
  virtual int      GetLgtPvt (void) { return 0; }
  virtual int      GetNgtPvt (void) { return 0; }
  virtual double   Get       (int, int) = 0;
  virtual void     Add       (int, int, double) = 0;
  virtual void     Print     (void) = 0;
  virtual void     MultVect  (double *, double *) = 0;
  // Visit the full, unfactorized matrix, including both symmetric halves.
  virtual void     VisitEntries(const std::function<void(int, int, double)> &) = 0;
  virtual void     AddMat    (double, cSysMatrix *) { }
};

// -------------------------------------------------------------------------
// Symmetric skyline matrices:
//
class cSymSkylMatrix : public cSysMatrix
{
 protected:
  int SolveNative(cVector &, cVector &);
  int    dec;     // Decomposition flag
  int    *skl;    // Matrix skyline
  double **val;   // Stored values

 public:
  void VisitEntries(const std::function<void(int, int, double)> &);
           cSymSkylMatrix(int, int *);
          ~cSymSkylMatrix(void);
  double&  operator()    (int i, int j) { return val[i][j]; }
  double **Val           (void) { return val; }
  int      Symmetric     (void) { return 1; }
  void     Zero          (void);
  int      Solve         (cVector &, cVector &);
  int      GetSmlPvt     (void);
  int      GetLgtPvt     (void);
  int      GetNgtPvt     (void);
  double   Get           (int, int);
  void     Add           (int, int, double);
  void     Print         (void);
  void     MultVect      (double *, double *);
  void     AddMat        (double, cSysMatrix *);
};

// -------------------------------------------------------------------------
// Unsymmetric matrices with symmetric skyline:
//
class cUnsymSkylMatrix : public cSysMatrix
{
 protected:
  int SolveNative(cVector &, cVector &);
  int    dec;     // Decomposition flag
  int    *skl;    // Matrix skyline
  double **L;     // Lower matrix
  double **U;     // Upper matrix

 public:
  void VisitEntries(const std::function<void(int, int, double)> &);
           cUnsymSkylMatrix(int, int *);
          ~cUnsymSkylMatrix(void);
  void     Zero          (void);
  int      Solve         (cVector &, cVector &);
  int      GetSmlPvt     (void);
  int      GetLgtPvt     (void);
  int      GetNgtPvt     (void);
  double   Get           (int, int);
  void     Add           (int, int, double);
  void     Print         (void);
  void     MultVect      (double *, double *);
  void     AddMat        (double, cSysMatrix *);
};

// -------------------------------------------------------------------------
// Sparse matrix element:
//
typedef struct _spelm SpElm;

struct _spelm
{
  int    col;
  double val;
  SpElm *nxt;
};

// -------------------------------------------------------------------------
// Class to handle nonsymentric sparse matrices:
//
class cSprsMatrix : public cSysMatrix
{
 protected:
  SpElm **row;   // Array of rows (list of nonzero elements)

  int *nc;

 public:
  void VisitEntries(const std::function<void(int, int, double)> &);
                  cSprsMatrix(int n = 1);
  virtual        ~cSprsMatrix(void);
  virtual void    Zero       (void);
  virtual int     Solve      (cVector &, cVector &);
  virtual double  Get        (int, int);
  virtual void    Add        (int, int, double);
  virtual void    Print      (void);
  virtual void    MultVect   (double *, double *);
  virtual void    MultVect   (cVector &, cVector &);  // For MPC solver 
  virtual void    MultTVect  (cVector &, cVector &);  // For MPC solver
  virtual void    MultMat    (cMatrix &, cMatrix &);  // For MPC solver
  virtual void    GetRow     (int, cVector &);        // For MPC solver
};

// -------------------------------------------------------------------------
// Class to handle symmetric sparse matrices:
//
class cSymSprsMatrix : public cSprsMatrix
{
 public:
          cSymSprsMatrix(int n = 1);
         ~cSymSprsMatrix(void) { }
  int     Symmetric     (void) { return 1; }
  int     Solve         (cVector &, cVector &);
  double  Get           (int, int);
  void    Add           (int, int, double);
  void    Print         (void);
  void    MultVect      (double *, double *);
};

// -------------------------------------------------------------------------
// Class to handle nonsymentric sparse matrices with Eigen lib:
//
class cSprsMatEigen : public cSprsMatrix
{
 public:
                  cSprsMatEigen(int n = 1) : cSprsMatrix(n) { }
  virtual        ~cSprsMatEigen(void) { }
  virtual int     Solve      (cVector &, cVector &);
};

// -------------------------------------------------------------------------
// Generic Solvers:
//
int PCGSolver(int maxit, double tol, cSysMatrix *K, cVector &f,
              cVector &u);

int FstEigPair(int maxit, double tol, cSysMatrix *K, double *, cVector &v);

#endif
