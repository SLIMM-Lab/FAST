// -------------------------------------------------------------------------
// ctrl.h - file containing the definition of the cControl class.
// -------------------------------------------------------------------------
//
// The cControl class is responsible for the control of the analysis
// process. It also performs several auxiliary tasks, as the assembly of
// global vectors and matrices, required by the analysis process. The
// typical problems considered are the linear static problems, linearized
// buckling, nonlinear path-following, linear and nonlinear dynamics.
//
// Control
// |-- LinearStatic
// |-- QuasiStatic
// |-- LoadControl (Newton-Raphson)
// |-- Incremental (Euler)
// |-- Secant
// |-- Path-Following
// |   |-- DisplacementControl
// |   |-- ArcLength
// |       |-- InitialOrthogonal (Riks)
// |       |-- UpdatedOrthogonal (Ramm)
// |       |-- ConsistentlyLinearized
// |       |-- Cilindrical (Crisfield)
// |-- LinearStability
// |-- NonlinearStability
// |-- Newmark
// |   |-- IncrementalNewmark
// |   |-- NonLinearNewmark
// |   |-- HHT
// |       |-- NonLinearHHT
// |   |-- GeneralizedAlpha 
// |       |-- NonlinearGeneralizedAlpha
// |-- Modal
// |-- HeatTransfer
//
// -------------------------------------------------------------------------
// Static public methods:
// -------------------------------------------------------------------------
//
// static void ReadAlg(void)
//
// This method reads the chosen algorithm from the input file, creates the
// algorithm object, and stores a pointer to it in the variable Algorithm.
// -------------------------------------------------------------------------
//
// static void ReadMaxStep(void)
//
// This method reads the maximum number of steps of a nonlinear analysis.
// -------------------------------------------------------------------------
//
// static void ReadMaxIter(void)
//
// This method reads the maximum number of iterations for convergence of
// the analysis algorithm.
// -------------------------------------------------------------------------
//
// static void ReadPrintStep(void)
//
// This method reads the steps used result printing, where 1 means that all
// steps are printed, 2 means that only even steps are printed and so on.
// -------------------------------------------------------------------------
//
// static void ReadPropLoad(void)
//
// This method reads a flag (0/1) indicating if the load is proportional.
// The use of this flag in nonlinear static analysis avoids the computation
// of the applied loads a each new analysis step.
// -------------------------------------------------------------------------
//
// static void ReadFactor(void)
//
// This method reads the incremental step factor.
// -------------------------------------------------------------------------
//
// static void ReadTimeStep(void)
//
// This method reads the time step used in the analysis of time dependent
// problems.
// -------------------------------------------------------------------------
//
// static void ReadTol(void)
//
// This method reads the tolerance require for the convergence of the
// iterations of the analysis algorithm.
// -------------------------------------------------------------------------
//
// static void ReadTimeStep(void)
//
// This method reads the time step used in the analysis of time dependent
// problems.
// -------------------------------------------------------------------------
//
// static void ReadPrintReactions(void)
//
// This method reads the flag controlling the printing of support reactions.
// -------------------------------------------------------------------------
//
// static void ReadPrintEqvMesh(void)
//
// This method reads the flag controlling the printing of equivalent mesh.
// -------------------------------------------------------------------------
//
// static void ReadPatStrSmt(void)
//
// This method reads the flag controlling stresses smooth at patch output
// points.
// -------------------------------------------------------------------------
//
// static void ReadLamLocStr(void)
//
// This method reads the flag controlling whether the output of laminated
// element stresses are calculated in laminate system.
// -------------------------------------------------------------------------
//
// static cControl *GetControl(void)
//
// This method returns a pointer to the control object.
// -------------------------------------------------------------------------
//
// static eCtrlType GetType(void)
//
// This method returns the type of the chosen algorithm.
// -------------------------------------------------------------------------
//
// static double GetStepFactor(void)
//
// This method returns the load step (increment) factor.
// -------------------------------------------------------------------------
//
// static double GetTotFactor(void)
//
// This method returns the total (accumulated) load factor.
// -------------------------------------------------------------------------
//
// static double GetTimeStep(void)
//
// This method returns the current time step.
// -------------------------------------------------------------------------
//
// static double GetTotTime(void)
//
// This method returns the total (accumulated) time.
// -------------------------------------------------------------------------
//
// static void Init(void)
//
// This method performs some tasks that should be completed prior to the
// beginning of the different solution procedures, as the numbering of the
// nodal dofs, the computation of the number of equations and of the
// profile of the stiffness matrix.
// -------------------------------------------------------------------------
//
// static void End(void)
//
// This method is responsible for the tasks that should be performed after
// the end of the numerical analysis, as the printing of the %END tag at
// the output file and the destruction of the objects created by the
// program.
// -------------------------------------------------------------------------
//
// static void ExtLoadVector(cVector &f)
//
//  f      -  external load vector                                   (out)
//
// This method computes the external load vector of the FE model. Since the
// time is not given it should be used only for static problems with
// proportional loads.
// -------------------------------------------------------------------------
//
// static void ExtLoadVector(double t, cVector &f)
//
//  t      -  current time                                            (in)
//  f      -  external load vector                                   (out)
//
// This method computes the external load vector of the FE model at the
// given time.
// -------------------------------------------------------------------------
//
// static int IntForceVector(cVector &g)
//
//  g      -  internal force vector                                  (out)
//
// This method computes the internal force vector of the FE model. It loops
// through the elements, gets the internal force vector of each element and
// add this vector to the global one. It returns (1) if the computation
// was successfully completed and (0) otherwise.
// -------------------------------------------------------------------------
//
// static void SuppReactions(cMatrix &R)
//
//  R      -  nodal reactions                                        (out)
//
// This method computes the nodal reactions at nodes with supports, spring
// supports and prescribed displacements. R = g - f is a matrix with nnode
// rows and 6 columns.
// -------------------------------------------------------------------------
//
// static void StiffnessMatrix(cSysMatrix *K)
//
//  K      -  global stiffness matrix                                (out)
//
// This method computes the stiffness matrix of the FE model. It loops
// through the elements, gets the stiffness matrix of each element
// and add this matrix to the global one.
// -------------------------------------------------------------------------
//
// static void MassMatrix(cSysMatrix *M)
//
//  M      -  global mass matrix                                     (out)
//
// This method computes the mass matrix of the FE model. It loops
// through the elements, gets the mass matrix of each element
// and add this matrix to the global one.
// -------------------------------------------------------------------------
//
// static void GeomStiffMatrix(cSysMatrix *G)
//
//  G      -  global geometric stiffness matrix                      (out)
//
// This method computes the geometric stiffness matrix of the FE model.
// It loops through the elements, gets the geometric stiffness matrix of
// each element and add this matrix to the global one.
// -------------------------------------------------------------------------
//
// static void AssignDispl(double *u)
//
//  u      -  global displacement vector                              (in)
//
// This method makes the results obtained by the global analysis algorithm
// accessible to the elements. To this aim, it extracts the displacements
// of each node from the global displacement vector and assigns zero to
// the displacements corresponding to fixed dofs.
// -------------------------------------------------------------------------
//
// static void PrintResult(void)
//
// This method prints in the output file the results computed in each step
// of the chosen analysis algorithm.
// -------------------------------------------------------------------------
//
// virtual void UpdateItera(cVector &du)
//
// This method loops through the mesh to update the internal variables of
// each finite element. It should be called after the correction of each
// nonlinear analysis iteration.
// -------------------------------------------------------------------------
//
// static void UpdateState(void)
//
// This method loops through the mesh to update the internal variables of
// each finite element. It should be called after the convergence of each
// nonlinear analysis step.
//
// -------------------------------------------------------------------------
// Pure virtual public methods:
// -------------------------------------------------------------------------
//
// virtual void Solver(void)
//
// This method implements the steps required to perform an specific
// analysis type: linear static, nonlinear static, linear dynamic, etc.
// It is responsible to obtain the required global matrices and vectors,
// compute the unknows of the problem, and print the computed results.
// -------------------------------------------------------------------------

#ifndef _CTRL_H
#define _CTRL_H

#include "sysmat.h"
#include <vector>

// -------------------------------------------------------------------------
// Forward declarations:
//
class cVector;
class cMatrix;
class cNode;

// -------------------------------------------------------------------------
// Micro model types:
//
typedef struct
{
  cNode *nodep; // plus
  cNode *nodem; // minus
} sMicroPBC;

typedef std::vector<cNode*> vMicroVtx;
typedef std::vector<sMicroPBC> vMicroPBC;

// -------------------------------------------------------------------------
// Control types:
//
typedef enum
{
  LINEAR_STATIC,   // Linear static analysis
  QUASISTATIC,     // Linear quasi-static analyis (e.g. viscoelasticity)
  LOAD_CONTROL,    // Nonlinear analysis - load control and NR iterations
  INCREMENTAL,     // Nonlinear analysis - pure incremental (Euler) method
  SECANT,          // Nonlinear analysis - secant method
  PATH_FOLLOWING,  // Nonlinear analysis - path-following methods
  LINEAR_STABILITY,// Computation of buckling loads and modes
  NONLIN_STABILITY,// Computation of nonlinear critical points
  NEWMARK,         // Linear dynamic analysis by the Newmark's method
  NEWMARK_INCR,    // Linear dyn. anl. in the incremental form
  NEWMARK_NL,      // Nonlinear dynamic analysis - Newmark Method
  HHT,             // Linear dynamic anl. - Hilber-Hughes-Taylor method
  HHT_NL,          // Nonlinear dyn. anl. - Hilber-Hughes-Taylor method
  GEN_ALPHA,       // Linear dynamic anl. - Generalized-Alpha method
  GEN_ALPHA_NL,    // Nonlinear dynamic anl. - Generalized-Alpha method
  MODAL,           // Computation of vibration modes and frequencies
  HEAT_TRANSFER    // Linear heat transfer analysis
} eCtrlType;

// -------------------------------------------------------------------------
// Definition of cControl base class:
//
class cControl
{
 protected:
  static bool      PrintReac;   // Print support reactions
  static eSMType   SysMatType;  // Type of system matrix
  static int       NumEq;       // Number of system equations
  static int      *Profile;     // Profile of the stiffness matrix
  static int       MaxStep;     // Maximum number of steps
  static int       MaxIter;     // Maximum number of iterations
  static int       TargetIter;  // Target number of iterations for the next step
  static int       PrevIter;    // Number of iterations of last step
  static int       CurrStep;    // Current step
  static int       CurrIter;    // Current iteration
  static int       PrintStep;   // Number of steps between result printing
  static int       PrintedSteps;// Auxiliar variable.
  static int       PropLoad;    // Flag for proportional loading (N-R)
  static int       OutPrec;     // Output precision.
  static double    StepFactor;  // Given load increment
  static double    TotFactor;   // Current load factor
  static double    TimeStep;    // Current time step
  static double    TotTime;     // Total time (current time)
  static double    Tol;         // Tolerance for convergence
  static cControl *Algorithm;   // Chosen control algorithm
  static eCtrlType Type;        // Type of chosen algorithm
  static long      CasePos[2];  // Position in the output file
  static bool      PrintEqvMsh; // Flag for print IGA equivalent mesh.
  static bool      PrintAvrStr; // Flag for print average nodal stresses.
  static bool      StrExtFlag;  // Stress extrapolation flag.
  static bool      PatStrSmt;   // Patch stress smoothing flag.
  static bool      LamLocStr;   // Output laminate local stress flag. 
  static int       EigAlgType;  // Eigen value problem algorithm type.
  static int       NumModes;    // Number of Eigen values/vectors computed.
  static int       ErrIntOrd;   // Integration order used on L2 error evaluation.
  // Temporary - Micromechanical analysis (Evandro)
  static int       MicroModel;  // Number of RVE dimensions (2D or 3D)
  static double    MicroVol;    // RVE volume
  static double    MicroRegFac; // Regularization factor
  static vMicroVtx MicroVecVtx; // RVE vertices
  static vMicroPBC MicroVecPBC; // PBC node  pairs
  static cVector   MacroStrain; // Macroscopic strains

 private:
  static  void      InitDofs(void);
  static  void      CompEqs(void);
  static  void      CompProf(void);
  static  void      PrintIntPntStress(void);
  static  void      PrintNodalStress(void);
  static  void      PrintAvrNodalStress(void);
  static  void      PrintAvrStress(void);
  static  void      PrintPatDispl(void);
  static  void      PrintPatStr(void);
  static  void      PrintEqvMshData(void);
  static  void      PrintEqvMshDispl(void);
  static  void      PrintEqvMshNodalStress(void); 
  static  void      PrintStrErrNormScalVis(int);

 protected:
  static  void      GetPrintIdx(int, int *, int, int *, int *);
  virtual int       Convergence(double, double);
  virtual int       Convergence(cVector &, cVector &);
  virtual int       Convergence(double, double, double);
  virtual void      CalcMacroStrLag(cVector &, cMatrix &);
  virtual void      MicroQMatrix2D(cMatrix &);
  virtual void      MicroQMatrix3D(cMatrix &);
  virtual int       MicroNumPBC(void);
  virtual void      MicroGetPBC(cMatrix &, cVector &);
  virtual void      MicroGetPBC(cSprsMatrix &, cVector &);

 public:
  static  void      ReadAlg(void);
  static  void      ReadMatrixType(void);
  static  void      ReadLinearSolver(void);
  static  void      ReadPetscOptions(void);
  static  void      ReadPetscOptionsFile(void);
  static  void      ReadMaxStep(void);
  static  void      ReadMaxIter(void);
  static  void      ReadTargetIter(void);
  static  void      ReadPrintStep(void);
  static  void      ReadPropLoad(void);
  static  void      ReadOutPrec(void);
  static  void      ReadFactor(void);
  static  void      ReadTimeStep(void);
  static  void      ReadTol(void);
  static  void      ReadPrintReactions(void);
  static  void      ReadPrintAverageStress(void);
  static  void      ReadPrintEqvMesh(void);
  static  void      ReadStrExtFlag(void);
  static  void      ReadPatStrSmt(void);
  static  void      ReadLamLocStr(void);
  static  void      ReadEigAlgType(void);
  static  void      ReadNumModes(void);
  static  void      ReadMicroModel(void);
  static  void      ReadMicroRegFac(void);
  static  void      ReadMicroPair(void);
  static  void      ReadMicroVtx(void);
  static  void      ReadMacroStrain(void);
  static  void      SetTargetIter(int nti) { TargetIter = nti; }
  static  cControl *GetControl(void) { return Algorithm; }
  static  eCtrlType GetType(void) { return Type; }
  static  int       GetTargetIter(void) { return TargetIter; }
  static  int       GetOutPrec(void) { return OutPrec; }
  static  double    GetStepFactor(void) { return StepFactor; }
  static  double    GetTotFactor(void) { return TotFactor; }
  static  double    GetTimeStep(void) { return TimeStep; }
  static  double    GetTotTime(void) { return TotTime; }
  static  double    GetIgaEqvMesh(void) { return PrintEqvMsh; }
  static  bool      GetStrExtFlag(void) { return StrExtFlag; }
  static  bool      GetPatStrSmt(void) { return PatStrSmt; }
  static  bool      GetLamLocStr(void) { return LamLocStr; }
  static  void      Init(void);
  static  void      End(void);
  static  void      ExtLoadVector(double, cVector &);
  static  void      ExtLoadVector(cVector &);
  static  int       IntForceVector(cVector &);
  static  void      SuppReactions(cMatrix &);
  static  void      StiffnessMatrix(cSysMatrix *);
  static  void      MassMatrix(cSysMatrix *);
  static  void      GeomStiffMatrix(cSysMatrix *);
  static  void      AssignDispl(cVector &);
  static  void      PrintResult(void);
  static  void      PrintStrainEnergy(void);
  static  void      PrintErrNorm(void);
  static  void      PrintStrErrNorm(int);
  static  void      PrintErrNorm2(void);
  static  void      UpdateItera(cVector &);
  static  void      UpdateState(void);
                    cControl(void);
  virtual          ~cControl(void);
  virtual void      Solver(void) = 0;
  static  void      RegisterDestroyFunc(void (*) (void));

  static  int       ExistStrainLoads    ( void );
  static  int       StrainLoadVector    ( cVector &, cVector & );
  static  void      SetTotFactor   (double tf) { TotFactor = tf;}
};

#endif
