// -------------------------------------------------------------------------
// ctrl.cpp - implementation of control class.
// -------------------------------------------------------------------------
// Created:      16-Nov-2000     Evandro Parente Junior
//
// Modified:     27-Nov-2000     Evandro Parente Junior
//               Implementation of the Newton-Rapson control class.
//
// Modified:     03-Oct-2001     Evandro Parente Junior
//               Use of vector/matrix classes.
//
// Modified:     11-Jul-2011     Rafael Fernandes da Silva
//               Implementation of vibration analysis.
//
// Modified:     08-Aug-2011     Evandro Parente Junior
//               Implementation of nodal constraints (master-slave).
//
// Modified:     12-Aug-2011     Evandro Parente Junior
//               Implementation of spring (elastic) supports.
//
// Modified:     22-Out-2011     Iuri Barcelos Rocha
//               Implementation of arc-length methods.
//
// Modified:     22-Nov-2011     Luiz Antonio Taumaturgo Mororo
//               Implementation of Nonlinear Newmark method.
//
// Modified:     04-Dec-2011     Iuri Barcelos Rocha
//               Implementation of linear stability analysis.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
//
// Modified:     12-Jul-2015     Elias Saraiva Barroso
//               Implementation of the stiffness matrix assembly
//               parallelization.
//
// Modified:     12-Jul-2015     Elias Saraiva Barroso
//               Implementation of nodal stress without extrapolation.
//
// Modified:     13-Oct-2015     Evandro Parente Junior
//               Evaluation and printing of support reactions.
//
// Modified:     07-Dec-2017     Bergson da Silva Matias
//               Implementation of Generalized-Alpha method.
//
// Modified:     07-Feb-2018     Elias Saraiva Barroso
//               Implementation of IGA post-processing + Destroy Stack
// -------------------------------------------------------------------------

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <stack>

using namespace std;

#include "ctrl.h"
#include "ctrllins.h"
#include "ctrlnr.h"
#include "ctrlinc.h"
#include "ctrlsec.h"
#include "ctrldsp.h"
#include "ctrlarclen.h"
#include "ctrlqstc.h"
#include "ctrlnewmark.h"
#include "ctrlhht.h"
#include "ctrlgenalpha.h"
#include "ctrlmodal.h"
#include "ctrllinstab.h"
#include "ctrlnlstab.h"
#include "ctrlheat.h"
#include "node.h"
#include "spring.h"
#include "element.h"
#include "elmdsh.h"
#include "shape.h"
#include "shpiga.h"
#include "load.h"
#include "prescdof.h"
#include "vec.h"
#include "mat.h"
#include "patch.h"
#include "sysmat.h"
#include "linsolver.h"
#include "matvec.h"
#include "eig.h"
#include "utl.h"
#include "gblvar.h"
#include "gbldef.h"
#include "utl.h"
#include "field.h"
#include "intpoint.h"

#ifdef _OMP_
#include "omp.h"
#endif

// -------------------------------------------------------------------------
// Static variables:
//
bool      cControl :: PrintReac   = false;
bool      cControl :: PrintAvrStr = false;
eSMType   cControl :: SysMatType  = SYM_SKYLINE;
int       cControl :: NumEq       = 0;
int*      cControl :: Profile     = 0;
int       cControl :: MaxStep     = 1;
int       cControl :: MaxIter     = 0;
int       cControl :: CurrStep    = 1;
int       cControl :: CurrIter    = 1;
int       cControl :: PrevIter    = 1;
int       cControl :: TargetIter  = 1;
int       cControl :: PrintStep   = 1;
int       cControl :: PropLoad    = 1;
int       cControl :: OutPrec     = 5;
int       cControl :: ErrIntOrd   = 15;
double    cControl :: TotFactor   = 1.0;
double    cControl :: StepFactor  = 1.0;
double    cControl :: TotTime     = 0.0;
double    cControl :: TimeStep    = 1.0;
double    cControl :: Tol         = 1.0e-05;
cControl* cControl :: Algorithm   = 0;
eCtrlType cControl :: Type;
long      cControl :: CasePos[2]  = {0,0};
bool      cControl :: PrintEqvMsh = false;
bool      cControl :: StrExtFlag  = true;
bool      cControl :: PatStrSmt   = true;
bool      cControl :: LamLocStr   = true;
int       cControl :: EigAlgType  = 0;
int       cControl :: NumModes    = 5;
int       cControl :: PrintedSteps = 0;
static stack<void (*) (void)> DestroyFunc;
//static FILE *graf = 0;

int       cControl :: MicroModel  = 0;      // Evandro
double    cControl :: MicroVol    = 0.0;    // Evandro
double    cControl :: MicroRegFac = 1.0e-8; // Evandro
vMicroVtx cControl :: MicroVecVtx;          // Evandro
vMicroPBC cControl :: MicroVecPBC;          // Evandro
cVector   cControl :: MacroStrain;          // Evandro

// -------------------------------------------------------------------------
// Public methods:
//

// ================================ ReadAlg ================================

void cControl :: ReadAlg(void)
{
  // Read the algorithm label.

  char label[100];
  if (!Utl::ReadString(in, label))
  {
    cout << "Error in the input of the algorithm label!\n";
    exit(0);
  }

  // Create the appropriate algorithm.

  if (string(label) == "Linear" || string(label) == "linear")
    Algorithm = new cLinStat( );
  else if (string(label) == "LoadControl" ||  string(label) == "lcm" ||
           string(label) == "NewtonRaphson" || string(label) == "snr")
    Algorithm = new cCtrlNR( );
  else if (string(label) == "Incremental" || string(label) == "Euler")
    Algorithm = new cCtrlIncremental( );
   else if (string(label) == "Secant" || string(label) == "secant")
    Algorithm = new cCtrlSecant( );
  else if (string(label) == "DisplacementControl" || string(label) == "dcm")
    Algorithm = new cCtrlDispl( );
  else if (string(label) == "QuasiStatic" || string(label) == "quasistatic")
    Algorithm = new cCtrlQuasiStatic( );
  else if (string(label) == "ModalAnalysis")
    Algorithm = new cCtrlModal( );
  else if (string(label) == "LinearStability")
    Algorithm = new cCtrlLinStab( );
  else if (string(label) == "NonlinearStability")
    Algorithm = new cCtrlNlStab( );
  else if (string(label) == "Newmark" || string(label) == "newmark")
    Algorithm = new cCtrlNewmark( );
  else if (string(label) == "NewmarkIncremental")
    Algorithm = new cCtrlIncNmk( );
  else if (string(label) == "NewmarkNonlinear" || string(label) == "newmarknl")
    Algorithm = new cCtrlNonlinNewmark( );
  else if (string(label) == "HHT" || string(label) == "hht")
    Algorithm = new cCtrlHHT( );
  else if (string(label) == "HHTNonlinear" || string(label) == "hhtnl")
    Algorithm = new cCtrlNonlinHHT( );
  else if (string(label) == "GenAlpha" || string(label) == "genalpha")
    Algorithm = new cCtrlGenAlpha( );
  else if (string(label) == "GenAlphaNonlinear" || string(label) == "genalphanl")
    Algorithm = new cCtrlNonlinGenAlpha( );
  else if (string(label) == "InitialOrthoArcLen" || string(label) == "ioal")
    Algorithm = new cIniOrtArcLength( );
  else if (string(label) == "UpdatedOrthoArcLen" || string(label) == "uoal")
    Algorithm = new cUpOrtArcLength( );
  else if (string(label) == "LinearizedArcLen" || string(label) == "clinal")
    Algorithm = new cConLinArcLength( );
  else if (string(label) == "CylindricalArcLen" || string(label) == "clal")
    Algorithm = new cCilArcLength( );
  else if (string(label) == "HeatTransfer" || string(label) == "heattransfer")
    Algorithm = new cHeatTransf( );
  else
  {
    cout << "Unknown algorithm: " << label << "\n";
    exit(0);
  }
}

// ============================= ReadMatrixType ============================

void cControl :: ReadMatrixType(void)
{
  // Read the matrix type.

  char label[100];
  if (!Utl::ReadString(in, label))
  {
    cout << "Error in the input of the algorithm label!\n";
    exit(0);
  }

  // Store the appropriate matrix type.

  if (string(label) == "skyline")
    SysMatType = SYM_SKYLINE;
  else if (string(label) == "unsym_skyline")
    SysMatType = UNSYM_SKYLINE;
  else if (string(label) == "sparse")
    SysMatType = SYM_SPARSE;
  else if (string(label) == "unsym_sparse")
    SysMatType = UNSYM_SPARSE;
  else if (string(label) == "eigen_unsym_sparse")
    SysMatType = EIGEN_UNSYM_SPARSE;
  else  //default
    SysMatType = SYM_SKYLINE;
}

void cControl::ReadLinearSolver() { LinearSolver::ReadBackend(in); }
void cControl::ReadPetscOptions() { LinearSolver::ReadOptions(in); }
void cControl::ReadPetscOptionsFile() { LinearSolver::ReadOptionsFile(in); }

// ============================== ReadMaxStep ==============================

void cControl :: ReadMaxStep(void)
{
  if (!(in >> MaxStep))
  {
    cout << "Error in the input of the number of steps!\n";
    exit(0);
  }
}

// ============================== ReadMaxIter ==============================

void cControl :: ReadMaxIter(void)
{
  if (!(in >> MaxIter))
  {
    cout << "Error in the input of the number of iterations!\n";
    exit(0);
  }
}

// ============================= ReadTargetIter ============================

void cControl :: ReadTargetIter(void)
{
  if (!(in >> TargetIter))
  {
    cout << "Error in the input of target iterations!\n";
    exit(0);
  }
}

// ============================= ReadPrintStep =============================

void cControl :: ReadPrintStep(void)
{
  if (!(in >> PrintStep))
  {
    cout << "Error in the input of the number of print steps!\n";
    exit(0);
  }
}

// ============================= ReadPropLoad ==============================

void cControl :: ReadPropLoad(void)
{
  if (!(in >> PropLoad))
  {
    cout << "Error in the input of the proportional loading flag!\n";
    exit(0);
  }
}

// ============================= ReadPropLoad ==============================

void cControl :: ReadOutPrec(void)
{
  if (!(in >> OutPrec) || (OutPrec < 0))
  {
    cout << "Error in the input of the post-processing decimal precision!\n";
    exit(0);
  }
}

// ============================== ReadFactor ===============================

void cControl :: ReadFactor(void)
{
  if (!(in >> StepFactor))
  {
    cout << "Error in the input of the step factor!\n";
    exit(0);
  }
}

// ============================= ReadTimeStep ==============================

void cControl :: ReadTimeStep(void)
{
  if (!(in >> TimeStep))
  {
    cout << "Error in the input of the time step!\n";
    exit(0);
  }
}

// ================================ ReadTol ================================

void cControl :: ReadTol(void)
{
  if (!(in >> Tol))
  {
    cout << "Error in the input of the algorithm tolerance!\n";
    exit(0);
  }
}

// =========================== ReadPrintReactions ==========================

void cControl :: ReadPrintReactions(void)
{
  int flag;
  if (!(in >> flag))
  {
    cout << "Error in the input of the flag for print reactions!\n";
    exit(0);
  }
  if (flag) PrintReac = true;
}

// =========================== ReadPrintAverageStress ======================

void cControl :: ReadPrintAverageStress(void)
{
  int flag;
  if (!(in >> flag))
  {
    cout << "Error in the input of the flag for print average nodal stress!\n";
    exit(0);
  }
  if (flag) PrintAvrStr = true;
}

// ========================== ReadPrintEqvMesh =============================

void cControl :: ReadPrintEqvMesh(void)
{
  if (!(in >> PrintEqvMsh))
  {
    cout << "Error in the input of the IGA equivalent mesh flag!\n";
    exit(0);
  }
}

// ========================== ReadStrExtFlag ===============================

void cControl :: ReadStrExtFlag(void)
{
  if (!(in >> StrExtFlag))
  {
    cout << "Error in the input of the stress extrapolation flag!\n";
    exit(0);
  }
}

// ========================== ReadPatStrSmt ================================

void cControl :: ReadPatStrSmt(void)
{
  if (!(in >> PatStrSmt))
  {
    cout << "Error in the input of the stress extrapolation flag!\n";
    exit(0);
  }
}

// ========================== ReadLamLocStr ================================

void cControl :: ReadLamLocStr(void)
{
  if (!(in >> LamLocStr))
  {
    cout << "Error in the input of the laminate local stress flag!\n";
    exit(0);
  }
}

// ========================== ReadStrExtFlag ===============================

void cControl :: ReadEigAlgType(void)
{
  eEigAlgType type;
  if (!(in >> type))
  {
    cout << "Error in the input of the eigenproblem type!\n";
    exit(0);
  }

  EigAlgType = static_cast<int> (type);
}

// ============================== ReadNumModes =============================

void cControl :: ReadNumModes(void)
{
  if (!(in >> NumModes))
  {
    cout << "Error in the input of the number of eigenmodes!\n";
    exit(0);
  }
}

// ============================ ReadMicroModel =============================

void cControl :: ReadMicroModel(void)
{
  if (!(in >> MicroModel) || (MicroModel < 2) || (MicroModel > 3) || 
      !(in >> MicroVol))
  {
    cout << "Error in the input of the dimension of the micromechanical model!\n";
    exit(0);
  }
//  cout << "MicroModel = " << MicroModel << "D\n";
//  cout << "MicroVol = " << scientific << MicroVol << endl;
}

// ============================ ReadMicroRegFac ============================

void cControl :: ReadMicroRegFac(void)
{
  if (!(in >> MicroRegFac) || (MicroRegFac < 0.0))
  {
    cout << "Error in the input of the regularization factor of the micromechanical model!\n";
    exit(0);
  }
//  cout << "MicroRegFac = " << scientific << MicroRegFac << endl;
}

// ============================= ReadMicroPair =============================

void cControl :: ReadMicroPair(void)
{
  int npair;
  if (!(in >> npair))
  {
    cout << "Error in the input of the number of RVE node pairs!\n";
    exit(0);
  }
//  cout << "Micromodel node pairs\n";
//  cout << npair << "\n";
  MicroVecPBC.resize(npair);

  int plus,minus;
  for (int i = 0; i < npair; i++) 
  {
    in >> plus >> minus;
//    cout << plus << "  " << minus << endl;
    MicroVecPBC[i].nodep = cNode :: FindNode(plus);
    MicroVecPBC[i].nodem = cNode :: FindNode(minus);
    if (!MicroVecPBC[i].nodep || !MicroVecPBC[i].nodem)
    {
      cout << "Error in the input of the RVE node pairs (invalid node)!\n";
      exit(0);
    }
  }
}

// ============================= ReadMicroVtx ==============================

void cControl :: ReadMicroVtx(void)
{
  int nvtx = 4;                   // 2D
  if (MicroModel == 3) nvtx = 8;  // 3D
  MicroVecVtx.resize(nvtx);

 // cout << "MicroVtx = ";
  int vtx;
  for (int i = 0; i < nvtx; i++) 
  {
    in >> vtx;
//    cout << vtx << "  ";
    MicroVecVtx[i] = cNode :: FindNode(vtx);
    if (!MicroVecVtx[i])
    {
      cout << "Error in the input of the RVE vertices (invalid node)!\n";
      exit(0);
    }
  }
//  cout << endl;
}

// ============================ ReadMacroStrain ============================

void cControl :: ReadMacroStrain(void)
{
  // Resize the strain vector.

  int nstr = 3;                   // 2D
  if (MicroModel == 3) nstr = 6;  // 3D
  MacroStrain.Resize(nstr);
  
  // Read macro strains.

  for (int i = 0; i < nstr; i++) 
  {
    in >> MacroStrain[i];
  }
//  cout << "{MacroStrain} = ";
//  MacroStrain.Print( );
//  cout << endl;
}

// ================================= Init ==================================

void cControl :: Init(void)
{
  // Build isogeometric equivalent mesh.

  if (PrintEqvMsh)
    PrintEqvMshData( );

  // Evaluate degenerated shell elements directors.
  cElmDegShell :: CompDegShlAxes( );

  // If necessary create the default algorithm.

  if (!Algorithm) Algorithm = new cLinStat( );

  // Mark the active dofs.

  InitDofs( );

  // Compute the node equations.

  CompEqs( );
  if (Feedback) cout << "\tNumber of equations = " << NumEq << endl;

  // Compute the matrix profile.

  CompProf( );
  if (SysMatType == SYM_SKYLINE || SysMatType == UNSYM_SKYLINE)
  {
    int mbw  = 0;   // Max bandwith
    int nlow = 0;   // Number of stored elements below the diagonal
    for (int i = 0; i < NumEq; i++)
    {
      int bw = i - Profile[i];
      nlow += bw;
      mbw   = MAX(mbw, bw);
    }

    if (Feedback)
    {
      if (SysMatType == SYM_SKYLINE)
      {
        cout << "\tSymmetric stiffness matrix\n";
        cout << "\tNumber of stored elements (skyline format) = ";
	cout << nlow + NumEq << endl;
      }
      else
      {
        cout << "\tUnsymmetric stiffness matrix\n";
        cout << "\tNumber of stored elements (skyline format) = ";
	cout << 2*nlow + NumEq << endl;
      }
      cout << "\tMaximum bandwith = " << mbw + 1 << "\n";
      cout << "\tAverage bandwith = " << fixed << setprecision(2);
      cout << double(nlow)/double(NumEq) + 1.0;
      cout << endl;
    }
  }
}

// ================================= End ===================================

void cControl :: End(void)
{
  // Print the end mark and the step labels.

  out << "\n%END\n";
  if (CasePos[0] != 0)
    out.seekp(CasePos[0], ios_base::beg); // fseek(pos, CasePos, SEEK_SET);

  out << "1   " << PrintedSteps << "\n";
  for (int i = 1; i <= PrintedSteps; i++)
  {
    out << setw(2) << i << "  'step" << i << "'\n";
  }

  if (PrintEqvMsh)
  {
    outem << "\n%END\n";
    if (CasePos[1] != 0)
      outem.seekp(CasePos[1], ios_base::beg); // fseek(pos, CasePos, SEEK_SET);
    
    outem << "1   " << PrintedSteps << "\n";
    for (int i = 1; i <= PrintedSteps; i++)
      outem << setw(2) << i << "  'step" << i << "'\n";
  }

//  fclosepgraf);

  // Release the class data.

  delete Algorithm;
  delete []Profile;

  // Destroy the model.

  cElement :: Destroy( );
  cNode :: Destroy( );

  while (DestroyFunc.size( ) > 0)
  {
    DestroyFunc.top( )( );
    DestroyFunc.pop( );
  }
}

// ============================= ExtLoadVector =============================

void cControl :: ExtLoadVector(cVector &f)
{
#if 0
  // Initialize the load vector.

  f.Zero( );

  // Compute and add the contribution of each element.

  int maxdof = cLoad :: GetMaxDof( );
  cVector felm(maxdof);
  for (cLoad *load = cLoad::GetHead( ); load != 0; load = load->GetNext( ))
  {
    load->ExtForce(1.0, felm);
    load->AddGlobVec(felm, f);
  }

  cout << "cControl :: ExtLoadVector (t = 1.0)\n";
#endif

  ExtLoadVector(1.0, f);
}

// ============================= ExtLoadVector =============================

void cControl :: ExtLoadVector(double t, cVector &f)
{
  // Initialize the load vector.

  f.Zero( );

  // Compute and add the contribution of each element.

  int maxdof = cLoad :: GetMaxDof( );
  cVector felm(maxdof);
  for (cLoad *load = cLoad::GetHead( ); load != 0; load = load->GetNext( ))
  {
    load->ExtForce(t, felm);
    load->AddGlobVec(felm, f);
  }

  if (ExistStrainLoads( ) && Type == LINEAR_STATIC)
  {
    cVector ft(NumEq);
    cVector fa(NumEq);     // aux. vector
    ft.Zero( );

    IntForceVector(fa);
    StrainLoadVector(fa, ft);

    f += ft;
  }

  // Modify the load vector to consider the prescribed dofs.

  int npd = cPrescDof :: GetNumPrescDofs( );
  if (npd > 0)
  {
    // Get the prescribed dofs and values.

    int *prdof = new int [npd];
    cVector prval(npd);
    cPrescDof :: GetPrescDofs(prdof);
    cPrescDof :: GetPrescVals(t, prval);

    // Add the penalty force.

    for (int i = 0; i < npd; i++)
      f[prdof[i]-1] += PENALTY_STIFFNESS*prval[i];

    delete [] prdof;
  }
}

// ============================ IntForceVector =============================

int cControl :: IntForceVector(cVector &g)
{
  int i;

  // Initialize the internal force vector.

  g.Zero( );

  // Compute and add the contribution of spring (elastic) supports.

  int nsd = cSpring :: GetMaxDof( );
  cVector gspr(nsd);
  int nss = cSpringSupp :: GetNumSpring( );
  for (i = 0; i < nss; i++)
  {
    // Get the current spring.

    cSpringSupp *spring = cSpringSupp :: GetSpring(i);

    // Evaluate the element internal force vector.

    spring->IntForce(gspr);

    // Add the element vector to the global vector.

    spring->AddGlobVec(gspr, g);
  }

  // Compute and add the contribution of spring (elastic) connections.

  int nsdc = cSpring :: GetMaxDof( );
  cVector gsprc(nsdc);
  int nssc = cSpringConn :: GetNumSpring( );
  for (i = 0; i < nssc; i++)
  {
    // Get the current spring.

    cSpringConn *spring = cSpringConn :: GetSpring(i);

    // Evaluate the element internal force vector.

    spring->IntForce(gsprc);

    // Add the element vector to the global vector.

    spring->AddGlobVec(gsprc, g);
  }

  // Alloc element internal force (max size).

  int maxdof = cElement :: GetMaxDof( );
  cVector gelm(maxdof);

  // Compute and add the contribution of each element.

  int code = 1;
  int nelm = cElement :: GetNumElm( );
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);

    // Evaluate the element internal force vector.

    code = elm->IntForce(gelm);
    if (!code) break;

    // Add the element vector to the global vector.

    elm->AddGlobVec(gelm, g);
  }

  // Modify the internal force vector to consider the prescribed dofs.

  int npd = cPrescDof :: GetNumPrescDofs( );
  if (npd > 0)
  {
    // Get the prescribed dofs and current nodal values.

    int* prdof = new int[npd];
    cVector crval(npd);
    cPrescDof :: GetPrescDofs(prdof);
    cPrescDof :: CurrNodalVals(crval);

    // Add the penalty force.

    for (i = 0; i < npd; i++)
      g[prdof[i]-1] += PENALTY_STIFFNESS*crval[i];

    delete[] prdof;
  }

  return(code);
}

// ============================ ExistStrainLoads ===========================

int cControl :: ExistStrainLoads(void)
{
  return(cElement :: ExistThermalLoad( ));
}

// ============================ StrainLoadVector ===========================

int cControl :: StrainLoadVector(cVector &g, cVector &dg)
{
  // Perturb the current factor and compute the forward internal force

  double delta  = 1.0e-06;
  double factor = cControl :: GetTotFactor( );
  cControl :: SetTotFactor(factor + delta);
  if (!IntForceVector(dg)) return(0);

  // Compute -d{g}/dfac by finite differences and restore the current factor

  dg  = g - dg;
  dg /= delta;
  cControl :: SetTotFactor(factor);

  return(1);
}

// ============================ SuppReactions ==============================

void cControl :: SuppReactions(cMatrix &R)
{
  // Initialize the support reactions.

  R.Zero( );

  // Add internal forces.

  int maxdof = cElement :: GetMaxDof( );
  cVector gelm(maxdof);
  int nelm = cElement :: GetNumElm( );
  for (int i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);

    // Compute the element internal forces and add to the nodal reactions.

    if (elm->HasSupp( ))
    {
      elm->IntForce(gelm);
      elm->AddReactions(gelm, R);
    }
  }

  // Add the external forces (negative sign).

  maxdof = cLoad :: GetMaxDof( );
  cVector felm(maxdof);
  for (cLoad *load = cLoad::GetHead( ); load != 0; load = load->GetNext( ))
  {
    if (load->HasSupp( ))
    {
      load->ExtForce(TotFactor, felm);
      load->AddReactions(felm, R);
    }
  }
}

// ============================ StiffnessMatrix ============================

void cControl :: StiffnessMatrix(cSysMatrix *K)
{
  int i;

  // Initialize the global stiffness matrix.

  K->Zero( );

  // Compute and add the contribution of spring (elastic) supports.

  int nsd = cSpring :: GetMaxDof( );
  cMatrix Kspr(nsd, nsd);
  int nss = cSpringSupp :: GetNumSpring( );
  for (i = 0; i < nss; i++)
  {
    // Get the current spring.

    cSpringSupp *spring = cSpringSupp :: GetSpring(i);

    // Evaluate the spring stiffness matrix.

    spring->StiffMat(Kspr);

    // Add the element spring to the global matrix.

    spring->AddGlobMat(Kspr, K);
  }

  // Compute and add the contribution of spring (elastic) connections.

  int nsdc = cSpring :: GetMaxDof( );
  cMatrix Ksprc(nsdc, nsdc);
  int nssc = cSpringConn :: GetNumSpring( );
  for (i = 0; i < nssc; i++)
  {
    // Get the current spring.

    cSpringConn *spring = cSpringConn :: GetSpring(i);

    // Evaluate the spring stiffness matrix.

    spring->StiffMat(Ksprc);

    // Add the element spring to the global matrix.

    spring->AddGlobMat(Ksprc, K);
  }

  // Alloc the element stiffness matrix (max size).

  int maxdof = cElement :: GetMaxDof( );
  cMatrix Kelm(maxdof, maxdof);

  // Compute and add the contribution of each element.

  int nelm = cElement :: GetNumElm( );
  cElement* elm = 0;
 
  #pragma omp parallel for firstprivate(elm,Kelm)
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.
    elm = cElement :: GetElm(i);

    // Evaluate the element stiffness matrix.
    elm->StiffMat(Kelm);

    // Add the element matrix to the global matrix.
    #pragma omp critical
    elm->AddGlobMat(Kelm, K);
  }

  // Modify the stiffness matrix to consider the prescribed dofs.

  int npd = cPrescDof :: GetNumPrescDofs( );
  if (npd > 0)
  {
    // Get array of prescribed dofs.

    int* prdof = new int[npd];
    cPrescDof :: GetPrescDofs(prdof);

    // Add the penalty stiffness to the matrix diagonal.

    for (i = 0; i < npd; i++)
      K->Add(prdof[i]-1, prdof[i]-1, PENALTY_STIFFNESS);

    delete[] prdof;
  }
}

// ============================ MassMatrix ============================

void cControl :: MassMatrix(cSysMatrix *M)
{
  int i,j;

  // Initialize the global stiffness matrix.

  M->Zero( );

  // Add the contribution of the nodal masses.

  int nnm = cNode :: GetNumMass( );
  int dof[6];
  for (i = 0; i < nnm; i++)
  {
    sNodeMass mass = cNode :: GetMass(i);

    for (j = 0; j < 6; j++) dof[j] = mass.node->GetDof(j);

    if (dof[0]) M->Add(dof[0]-1, dof[0]-1, mass.mx);
    if (dof[1]) M->Add(dof[1]-1, dof[1]-1, mass.my);
    if (dof[2]) M->Add(dof[2]-1, dof[2]-1, mass.mz);
    if (dof[3]) M->Add(dof[3]-1, dof[3]-1, mass.mrx);
    if (dof[4]) M->Add(dof[4]-1, dof[4]-1, mass.mry);
    if (dof[5]) M->Add(dof[5]-1, dof[5]-1, mass.mrz);
  }

  // Alloc the mass matrix (max size).

  int maxdof = cElement :: GetMaxDof( );
  cMatrix Melm(maxdof, maxdof);

  // Compute and add the contribution of each element.

  int nelm = cElement :: GetNumElm( );
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);

    // Evaluate the element mass matrix.

    elm->MassMat(Melm);

    // Add the element matrix to the global matrix.

    elm->AddGlobMat(Melm, M);
  }
}

// ============================ GeomStiffMatrix ============================

void cControl :: GeomStiffMatrix(cSysMatrix *G)
{
  int i;

  // Initialize the global geometric stiffness matrix

  G->Zero();

  // Alloc the element geometric stiffness matrix (max size).

  int maxdof = cElement :: GetMaxDof( );
  cMatrix Gelm(maxdof, maxdof);

  // Compute and add the contribution of each element.

  int nelm = cElement :: GetNumElm( );
  #pragma omp parallel for firstprivate(Gelm)
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);

    // Evaluate the element stiffness matrix.

    elm->GeomStiff(Gelm);

    // Add the element matrix to the global matrix.

    #pragma omp critical
    elm->AddGlobMat(Gelm, G);
  }

}

// ============================== PrintEqvMshData ==========================

void cControl :: PrintEqvMshData( )
{
  if (Feedback) 
    cout << "\n\tBuilding IGA equivalent mesh ............" << endl;
  
  // Process equivalent node list.
  cShapeIGA :: LoadEqvMsh( ); 

  // Open the equivalent mesh output file stream.
  string aux = fname + "_EqvMsh.pos";
  outem.open(aux.c_str( ));
  if (!outem.is_open())
  {
    cout << "Invalid output file." << "\n";
    exit(1);
  }

  outem << "%HEADER\n";
  outem << "Equivalent mesh output file generated by FAST.\n\n";

  // Print equivalent mesh.
  cShapeIGA :: PrintEqvMshData(outem);
}

// ============================== AssignDispl ==============================

void cControl :: AssignDispl(cVector &u)
{
  int nnode = cNode :: GetNumNode( );
  for (int i = 0; i < nnode; i++)
  {
    // Get the current node.
    cNode *node = cNode :: GetNode(i);

    // Get the nodal displacements from the global vector.

    for (int j = 0; j < 6; j++)
    {
      int dof = node->GetDof(j);
      if (dof)
        node->SetDispl(j, u[dof-1]);
      else
        node->SetDispl(j, 0.0);
    }
  }
}

// ============================== PrintResult ==============================

void cControl :: PrintResult(void)
{
  // Begin printing.

  Printing = true;

  // Print the result header.

  if (++PrintedSteps == 1)
  {
    out << "%RESULT\n";
    out << "1\n";
    out << "1  'result1'\n\n";

    out << "%RESULT.CASE" << endl; // Do not change endl to \n (use of tellp) 
    CasePos[0] = out.tellp( );        // Store current position
    for (int i = 0; i <= 1.1*MaxStep; i++)
      out << "             \n";

    if (PrintEqvMsh)
    {
      outem << "%RESULT\n";
      outem << "1\n";
      outem << "1  'result1'\n\n";
      
      outem << "%RESULT.CASE" << endl; // Do not change endl to \n (use of tellp) 
      CasePos[1] = outem.tellp( );        // Store current position
      for (int i = 0; i <= 1.1*MaxStep; i++)
        outem << "             \n";
    }
  }

  // Print the step header.

  out << "%RESULT.CASE.STEP\n";
  out << PrintedSteps << "\n\n";

  out << "%RESULT.CASE.STEP.FACTOR\n";
  out << scientific << setprecision(OutPrec) << TotFactor << "\n\n";

  if (PrintEqvMsh)
  {
    outem << "%RESULT.CASE.STEP\n";
    outem << PrintedSteps << "\n\n";
  
    outem << "%RESULT.CASE.STEP.FACTOR\n";
    outem << scientific << setprecision(OutPrec) << TotFactor << "\n\n";
  }

  // Print the nodal displacements.

  out << "%RESULT.CASE.STEP.NODAL.DISPLACEMENT\n";
  int nnode = cNode :: GetNumNode( );
  out << nnode << "  'displacement'\n";
//  double u1,u2;

  double dmax[3], dmin[3];
  for(int i = 0; i < 3; ++i) 
    dmax[i] = dmin[i] = cNode :: GetNode(0)->GetDispl(i);

  for (int i = 0; i < nnode; i++)
  {
    // Get the current node.

    cNode *node = cNode :: GetNode(i);
    out << left << setw(5) << node->GetLabel( ) << " ";

    // Get the nodal displacements.

    out << showpos << scientific << setprecision(OutPrec);
    for (int j = 0; j < 6; j++) out << node->GetDispl(j) << "  ";
    out << "\n" << resetiosflags(ios::showpos | ios::scientific);

    for(int i = 0; i < 3; ++i) dmax[i] = max(dmax[i],node->GetDispl(i));
    for(int i = 0; i < 3; ++i) dmin[i] = min(dmax[i],node->GetDispl(i));
    
    int nodeid = 0;
    //nodeid = 43;  // Scordelis-Lo 
    //nodeid = 52;  // Cantilever Beam Q8
    //nodeid = 21;  // Slit Annular Plate Q8 (A = 1, B = 21)
    //nodeid = 289;  // Hemispherical Shell Q8 (A = 17, B = 289)
    //nodeid = 284; // Pull-out Open Cylinder Q8 (A = 426, B = 337, C = 284)
    //nodeid = 288; // Pinched Cylinder Q8
    //nodeid = 2;   // Cylindrical Panel/Roof Q8 and StripFooting
    //nodeid = 27;   // Thermal Buckling Q8
    //nodeid = 4;  // Pinched Cylinder (A = 1, B = 4)
    if (node->GetLabel( ) == nodeid)
    {	    
      cout << "Node " << nodeid << "  ";
      cout << "u = "  << showpos << scientific << setprecision(5) << node->GetDispl(0) << "  ";
      cout << "v = "  << showpos << scientific << setprecision(5) << node->GetDispl(1) << "  ";
      cout << "w = "  << showpos << scientific << setprecision(5) << node->GetDispl(2) << "  ";
      cout << "LF = " << noshowpos << scientific << setprecision(OutPrec) << TotFactor << "\n";
    }
  }

  //cout << "max(u) = " << showpos << scientific << setprecision(5) << dmax[0] << endl;
  //cout << "max(v) = " << showpos << scientific << setprecision(5) << dmax[1] << endl;
  //cout << "max(w) = " << showpos << scientific << setprecision(5) << dmax[2] << endl;

  //cout << "min(u) = " << showpos << scientific << setprecision(5) << dmin[0] << endl;
  //cout << "min(v) = " << showpos << scientific << setprecision(5) << dmin[1] << endl;
  //cout << "min(w) = " << showpos << scientific << setprecision(5) << dmin[2] << endl;

  out << "\n";
 
  // Print support reactions.

  if (PrintReac)
  {
    // Evaluate the number of reaction nodes.

    int nrn = 0;  // Number nodes with reactions
    for (int i = 0; i < nnode; i++)
    {
      cNode *node = cNode::GetNode(i);
      if (node->HasSupp( )) nrn++;
    }

    // Evaluate and print the nodal reactions.

    if (nrn > 0)
    {
      out << "%RESULT.CASE.STEP.NODAL.REACTIONS\n";
      out << nrn << "\n";

      cMatrix R(nnode, 6);
      SuppReactions(R);
      for (int i = 0; i < nnode; i++)
      {
        cNode *node = cNode::GetNode(i);
        if (node->HasSupp( ))
        {
          out << left << setw(4) << fixed << node->GetLabel( ) << " ";
          out << scientific << right << setprecision(5);
          for (int j = 0; j < 6; j++) out << setw(13) << R(i, j) << " ";
          out << "\n";
        }
      }
      out << "\n";
    }
  }

  // Print stress and energy norms.

  // PrintStrainEnergy( );
  if (cField :: GetInstance( )->hasStressField( ))
  {  
    PrintStrErrNorm(0);
    if (cField :: GetInstance( )->hasStrainField( ))
      PrintStrErrNorm(1);

    //PrintStrErrNormScalVis(0);
  }

  // Print IGA Patch displacement points

  PrintPatDispl( );

  // Print stresses (in global or local system depending on the material).
  
  PrintIntPntStress( );
  PrintNodalStress( );

  // Print average stresses.

  if (PrintAvrStr)
  {
    PrintAvrNodalStress( );
    if (MicroModel) PrintAvrStress( );   // Micromechanical analysis 
  }

  // Print IGA Patch Stress points.

  PrintPatStr( );

  // Print IGA equivament nodal displacements and stresses.
  
  if (PrintEqvMsh)
  {
    PrintEqvMshDispl( );
    PrintEqvMshNodalStress( );
  }
 
  // End printing.

  Printing = false;
}

// ============================== UpdateItera ==============================

void cControl :: UpdateItera(cVector &du)
{
  int nelm = cElement :: GetNumElm( );
  for (int i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);

    // Update the element internal variables.

    elm->UpdateItera(du);
  }
}

// ============================== UpdateState ==============================

void cControl :: UpdateState(void)
{
  int nelm = cElement :: GetNumElm( );
  for (int i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);

    // Update the element internal variables.

    elm->UpdateState( );
  }
}

// =============================== cControl ================================

cControl :: cControl(void)
{
}

// =============================== ~cControl ===============================

cControl :: ~cControl(void)
{
}


// -------------------------------------------------------------------------
// Protected methods:
//

// ============================== Convergence ==============================
//
// This method checks the convergence of the equilibrium iterations using
// the criteria |r| < |q|*Tol. It returns (1) if the convergence criteria
// is satisfied and (0) otherwise.
//
//   lenq - length (|q|) of the reference load vector                  (in)
//   r    - residual (out-of-balance) force vector                     (in)
//
int cControl :: Convergence(double lenq, double lenr)
{
  double nrm = MAX(1.0, lenq);
  double err = lenr/nrm;

  if (Feedback)
  {
    cout << "  Iter: " << right << setw(3) << CurrIter;
    cout << "  Err: "  << scientific  << setprecision(4) << err;
    cout << "  LF: "   << fixed << showpos << TotFactor << endl << noshowpos;
  }

  if (err < Tol) return(1);
  return(0);
}

// ============================== Convergence ==============================
//
// This method checks the convergence of the equilibrium iterations using
// the criteria |r| < |q|*Tol. It returns (1) if the convergence criteria
// is satisfied and (0) otherwise.
//
//   lenq - length (|q|) of the reference load vector                  (in)
//   r    - residual (out-of-balance) force vector                     (in)
//
int cControl :: Convergence(cVector &q, cVector &r)
{
  int npd = cPrescDof :: GetNumPrescDofs( );
  if (npd == 0)
    return(Convergence(q.Length( ), r.Length( )));

  // Exclude the prescribed dofs.

  cVector qaux = q;
  cVector raux = r;
  int* prdof = new int[npd];
  cPrescDof :: GetPrescDofs(prdof);
  for (int i = 0; i < npd; i++)
  {
    int dof = prdof[i] - 1;
    qaux[dof] = 0.0;
    raux[dof] = 0.0;
  }

  delete[] prdof;

  return(Convergence(qaux.Length( ), raux.Length( )));
}

// ============================== Convergence ==============================
//
// This method checks the convergence of the equilibrium iterations using
// the criteria |r| < |q|*Tol. It returns (1) if the convergence criteria
// is satisfied and (0) otherwise.
//
//   lenq - length (|q|) of the reference load vector                  (in)
//   r    - residual (out-of-balance) force vector                     (in)
//
int cControl :: Convergence(double lenq, double lenr, double time)
{
  double nrm = MAX(1.0, lenq);
  double err = lenr/nrm;

  if (Feedback)
  {
    cout << "Iter: "  << right << setw(3) << CurrIter;
    cout << "  Err: " << scientific  << setprecision(4) << err;
    cout << "  Time: " << scientific << time << endl;
  }

  if (err < Tol) return(1);
  return(0);
}

// =========================== CalcMacroStrLag =============================

void cControl :: CalcMacroStrLag(cVector &lm, cMatrix &A)
{
  cout << "-----------------------------------------------------\n";
  cout << "MicroModel = " << MicroModel << "D\n";
  if (MicroVecVtx.size() == 0)
  {
    cout << "Invalid number of vertex nodes.\n";
    cout << "-----------------------------------------------------\n";
    return;
  }

  // Set the number of stress components.

  int nstr = 3;                   // 2D
  if (MicroModel == 3) nstr = 6;  // 3D

  // Assembly [Q].

  int nc = lm.Dim( ); // Number of constraints
  cMatrix Q(nstr, nc);
  if (MicroModel == 2)
    MicroQMatrix2D(Q);
  else
    MicroQMatrix3D(Q);

  // Evaluate the macroscopic stresses {sigM}.

  cVector sigM(nstr);
  sigM = (-1.0/MicroVol)*Q*lm;
  cout << "{strM} = " << scientific << setprecision(6);
  sigM.Print( );
  cout << "\n";

  // Evaluate the macroscopic constitutive matrix [CM].

  cMatrix CM(nstr, nstr);
  cMatrix D(nc, nstr);     // [D] = ([A]^-1)[Q]
  cVector v(nc);
  for (int i = 0; i < nstr; i++)  // [A][D] = [Q] 
  {
    for (int j = 0; j < nc; j++) v[j] = Q[i][j]; 
    A.SolveLU(v, v);
    for (int j = 0; j < nc; j++) D[j][i] = v[j];
  }
  CM = (1.0/MicroVol)*Q*D;
  cout << "[CM] = \n";
  CM.Print( );
  cout << "-----------------------------------------------------\n";
}
      
// ============================ MicroQMatrix2D =============================

void CalcH2D(double x, double y, cMatrix &H)
{
  H[0][0] = x;
  H[0][1] = 0.0;

  H[1][0] = 0.0;
  H[1][1] = y;

  H[2][0] = y/2.0;
  H[2][1] = x/2.0;
}

void cControl :: MicroQMatrix2D(cMatrix &Q)
{
  // Initialization.

  int ndsp = 2;   // Number of displacements
  int nstr = 3;   // Number of stresses
  int nvtx = 4;   // Number of vertices
  int j = 0;      // Counter
  Q.Zero( );

  // Vertex nodes.

  cMatrix H(nstr, ndsp);
  for (int i = 1; i < nvtx; i++) 
  {
    sNodeCoord coord = MicroVecVtx[i]->GetCoord( );
    CalcH2D(coord.x, coord.y, H);
//    cout << "Node = " << MicroVecVtx[i]->GetLabel( ) << "\n";
//    cout << "[H] = \n";
//    H.Print( );
    Q.Set(0, ndsp*j++, H);
  }
//  cout << endl;

  // Edge nodes.

  int npair = MicroVecPBC.size( );
  for (int i = 0; i < npair; i++) 
  {
    sNodeCoord p = MicroVecPBC[i].nodep->GetCoord( ); // Plus node
    sNodeCoord m = MicroVecPBC[i].nodem->GetCoord( ); // Minus node
    CalcH2D(p.x - m.x, p.y - m.y, H);
//    cout << "Plus = " << MicroVecPBC[i].nodep->GetLabel( ) << "  ";
//    cout << "Minus = " << MicroVecPBC[i].nodem->GetLabel( ) << "\n";
//    cout << "[H] = \n";
//    H.Print( );
    Q.Set(0, ndsp*j++, H);
  }
//  cout << endl;
}

// ============================ MicroQMatrix3D =============================

void CalcH3D(double x, double y, double z, cMatrix &H)
{
  H[0][0] = x;
  H[0][1] = 0.0;
  H[0][2] = 0.0;

  H[1][0] = 0.0;
  H[1][1] = y;
  H[1][2] = 0.0;

  H[2][0] = 0.0;
  H[2][1] = 0.0; 
  H[2][2] = z;

  H[3][0] = y/2.0;
  H[3][1] = x/2.0;
  H[3][2] = 0.0;

  H[4][0] = z/2.0;
  H[4][1] = 0.0;
  H[4][2] = x/2.0;

  H[5][0] = 0.0;
  H[5][1] = z/2.0;
  H[5][2] = y/2.0;
}

void cControl :: MicroQMatrix3D(cMatrix &Q)
{
  // Initialization.

  int ndsp = 3;   // Number of displacements
  int nstr = 6;   // Number of stresses
  int nvtx = 8;   // Number of vertices
  int j = 0;      // Counter
  Q.Zero( );

  // Vertex nodes.

  cMatrix H(nstr, ndsp);
  for (int i = 1; i < nvtx; i++) 
  {
    sNodeCoord coord = MicroVecVtx[i]->GetCoord( );
    CalcH3D(coord.x, coord.y, coord.z, H);
//    cout << "Node = " << MicroVecVtx[i]->GetLabel( ) << "\n";
//    cout << "[H] = \n";
//    H.Print( );
    Q.Set(0, ndsp*j++, H);
  }
//  cout << endl;

  // Edge and face nodes.

  int npair = MicroVecPBC.size( );
//  cout << "npair = " << npair << "\n";
  for (int i = 0; i < npair; i++) 
  {
    sNodeCoord p = MicroVecPBC[i].nodep->GetCoord( ); // Plus node
    sNodeCoord m = MicroVecPBC[i].nodem->GetCoord( ); // Minus node
    CalcH3D(p.x - m.x, p.y - m.y, p.z - m.z, H);
//    cout << "Plus = " << MicroVecPBC[i].nodep->GetLabel( ) << "  ";
//    cout << "Minus = " << MicroVecPBC[i].nodem->GetLabel( ) << "\n";
//    cout << "[H] = \n";
//    H.Print( );
    Q.Set(0, ndsp*j++, H);
  }
//  cout << endl;
}
 
// ============================== MicroNumPBC ==============================

int cControl :: MicroNumPBC(void)
{
  int ndsp = 2;                   // 2D
  if (MicroModel == 3) ndsp = 3;  // 3D

  int nvtx = MicroVecVtx.size( );
  int npair = MicroVecPBC.size( );
  int nc = ndsp*(nvtx - 1 + npair);
  //cout << "nc = " << nc << endl;
  return(nc);
}

// ============================== MicroGetPBC ==============================

void cControl :: MicroGetPBC(cMatrix &P, cVector &b)
{
  //cout << "MicroModel = " << MicroModel << "D\n";
  //cout << "Regularization Factor = " << scientific << MicroRegFac << endl;
  //cout << "{StrainM} = " << scientific;
  //MacroStrain.Print( );
  //cout << endl;

  // Initialization.

  int ndsp = 2;          // Number of displacements
  int nstr = 3;          // Number of stresses
  if (MicroModel == 3)
  {
    ndsp = 3;
    nstr = 6;
  }
  b.Zero( );
  P.Zero( );
  
  // Vertex nodes.

  int nvtx = MicroVecVtx.size( );
  //cout << "Number of vertex PBC constraints = " << ndsp*(nvtx - 1) << "\n";
  int k = 0;
  cVector u(ndsp);
  cMatrix H(nstr, ndsp);
  for (int i = 1; i < nvtx; i++) // Exclude the first vertex
  {
    //cout << "Node = " << setw(3) << MicroVecVtx[i]->GetLabel( ) << "  ";
    sNodeCoord coord = MicroVecVtx[i]->GetCoord( );
    if (MicroModel == 2)
      CalcH2D(coord.x, coord.y, H);
    else
      CalcH3D(coord.x, coord.y, coord.z, H);
    u = t(H)*MacroStrain;
    //cout << "{b} = " << scientific;
    //u.Print( );
    for (int j = 0; j < ndsp; j++)
    {
      int dof = MicroVecVtx[i]->GetDof(j) - 1;
      P[k][dof] = 1.0;
      b[k] = u[j];
      k++;
    }
  }
  //cout << endl;

  // Edge and face nodes.

  int npair = MicroVecPBC.size( );
  //cout << "Number of edge/face PBC constraints = " << ndsp*npair << "\n";
  for (int i = 0; i < npair; i++) 
  {
    //cout << "Plus = " << setw(3) << MicroVecPBC[i].nodep->GetLabel( ) << "  ";
    //cout << "Minus = " << setw(3) << MicroVecPBC[i].nodem->GetLabel( ) << "  ";
    sNodeCoord p = MicroVecPBC[i].nodep->GetCoord( ); // Plus node
    sNodeCoord m = MicroVecPBC[i].nodem->GetCoord( ); // Minus node
    if (MicroModel == 2)
      CalcH2D(p.x - m.x, p.y - m.y, H);
    else
      CalcH3D(p.x - m.x, p.y - m.y, p.z - m.z, H);
    u = t(H)*MacroStrain;
    //cout << "{b} = " << scientific;
    //u.Print( );
    for (int j = 0; j < ndsp; j++)
    {
      int dofp = MicroVecPBC[i].nodep->GetDof(j) - 1;
      int dofm = MicroVecPBC[i].nodem->GetDof(j) - 1;
      P[k][dofp] =  1.0;
      P[k][dofm] = -1.0;
      b[k] = u[j];
      k++;
    }
  }
  //cout << "\nNumber of PBC constraints = " << k << endl;
}

// ============================== MicroGetPBC ==============================

void cControl :: MicroGetPBC(cSprsMatrix &P, cVector &b)
{
  //cout << "MicroModel = " << MicroModel << "D\n";
  //cout << "Regularization Factor = " << scientific << MicroRegFac << endl;
  //cout << "{StrainM} = " << scientific;
  //MacroStrain.Print( );
  //cout << endl;

  // Initialization.

  int ndsp = 2;          // Number of displacements
  int nstr = 3;          // Number of stresses
  if (MicroModel == 3)
  {
    ndsp = 3;
    nstr = 6;
  }
  b.Zero( );
  
  // Vertex nodes.

  int nvtx = MicroVecVtx.size( );
  //cout << "Number of vertex PBC constraints = " << ndsp*(nvtx - 1) << "\n";
  int k = 0;
  cVector u(ndsp);
  cMatrix H(nstr, ndsp);
  for (int i = 1; i < nvtx; i++) // Exclude the first vertex
  {
    //cout << "Node = " << setw(3) << MicroVecVtx[i]->GetLabel( ) << "  ";
    sNodeCoord coord = MicroVecVtx[i]->GetCoord( );
    if (MicroModel == 2)
      CalcH2D(coord.x, coord.y, H);
    else
      CalcH3D(coord.x, coord.y, coord.z, H);
    u = t(H)*MacroStrain;
    //cout << "{b} = " << scientific;
    //u.Print( );
    for (int j = 0; j < ndsp; j++)
    {
      int dof = MicroVecVtx[i]->GetDof(j) - 1;
      P.Add(k, dof, 1.0);
      b[k] = u[j];
      k++;
    }
  }
  //cout << endl;

  // Edge and face nodes.

  int npair = MicroVecPBC.size( );
  //cout << "Number of edge/face PBC constraints = " << ndsp*npair << "\n";
  for (int i = 0; i < npair; i++) 
  {
    //cout << "Plus = " << setw(3) << MicroVecPBC[i].nodep->GetLabel( ) << "  ";
    //cout << "Minus = " << setw(3) << MicroVecPBC[i].nodem->GetLabel( ) << "  ";
    sNodeCoord p = MicroVecPBC[i].nodep->GetCoord( ); // Plus node
    sNodeCoord m = MicroVecPBC[i].nodem->GetCoord( ); // Minus node
    if (MicroModel == 2)
      CalcH2D(p.x - m.x, p.y - m.y, H);
    else
      CalcH3D(p.x - m.x, p.y - m.y, p.z - m.z, H);
    u = t(H)*MacroStrain;
    //cout << "{b} = " << scientific;
    //u.Print( );
    for (int j = 0; j < ndsp; j++)
    {
      int dofp = MicroVecPBC[i].nodep->GetDof(j) - 1;
      int dofm = MicroVecPBC[i].nodem->GetDof(j) - 1;
      P.Add(k, dofp,  1.0);
      P.Add(k, dofm, -1.0);
      b[k] = u[j];
      k++;
    }
  }
  //cout << "\nNumber of PBC constraints = " << k << endl;
}

// ============================ PrintAvrStress =============================

void cControl :: PrintAvrStress(void)
{
  // Initialize data.

  double elmvol;
  double totvol = 0.0;
  int ngsl = cElement :: GetNumSclLab( );
  cVector elmstr(ngsl);
  cVector totstr(ngsl);
  totstr.Zero( );

  // Compute the mesh volume and stress integral. 

  int nelm = cElement :: GetNumElm( );
  for (int i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);

    // Add the element volume and stress integral.

    elm->CalcStrInt(elmvol, elmstr);
    totvol += elmvol;
    totstr += elmstr;
  }
  //cout << "totvol = " << totvol << "\n";

  // Get the dimensions of the micromechanical model (RVE).
 
  int nnode = cNode :: GetNumNode( );
  double Lx = 0.0;
  double Ly = 0.0;
  double Lz = 0.0;
  for (int i = 0; i < nnode; i++)
  {
    // Get the current node coords.

    cNode *node = cNode::GetNode(i);
    sNodeCoord p = node->GetCoord( );

    // Update the dimensions.

    if (p.x > Lx) Lx = p.x;
    if (p.y > Ly) Ly = p.y;
    if (p.z > Lz) Lz = p.z;
  }
  //cout << "Lx = " << Lx << "\n";
  //cout << "Ly = " << Ly << "\n";
  //cout << "Lz = " << Lz << "\n";

  // Compute and print the average stresses.

  double RVEvol = Lx*Ly;              // 2D
  if (MicroModel == 3) RVEvol *= Lz;  // 3D
  cout << "\nRVEvol = " << RVEvol << "\n";
  totstr /= RVEvol;
  cout << "{AvrStr} = " << scientific << setprecision(6);
  totstr.Print( );

  // Initialize the support reactions.

  int maxdof = cElement :: GetMaxDof( );
  cVector gelm(maxdof);
  cMatrix R(nnode, maxdof);
  R.Zero( );

  // Add internal forces.

  for (int i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);

    // Compute the element internal forces and add to the nodal reactions.

    if (1) // (elm->HasSupp( )) // Evandro: temporary for microscale testing
    {
      elm->IntForce(gelm);
      elm->AddReactions(gelm, R);
    }
  }
  //cout << "[R] = \n";
  //R.Print( );
  //cout << "\n";

  // Evalute the macroscopic stresses using the nodal forces.

  if (MicroModel == 2)  // 2D
  {
    // Evaluate the resulting forces at right and top edges.
 
    int ndsp = 2;   // Number of displacements
    int nstr = 3;   // Number of stresses
    double Tol = 1.0e-6*min(Lx, Ly);
    cVector gr(ndsp);  // Right
    cVector gt(ndsp);  // Top
    gr.Zero( ); 
    gt.Zero( ); 
    for (int i = 0; i < nnode; i++)
    {
      // Get the current node.
 
      cNode *node = cNode::GetNode(i);
      sNodeCoord p = node->GetCoord( );
 
      // Add to the resulting force of right nodes.
 
      if (abs(p.x - Lx) <= Tol)
        for (int j = 0; j < ndsp; j++) gr[j] += R[i][j];
 
      // Add to the resulting force of top nodes.
 
      if (abs(p.y - Ly) <= Tol)
        for (int j = 0; j < ndsp; j++) gt[j] += R[i][j];
    }
    //cout << "{gr} = ";
    //gr.Print( );
    //cout << "\n";
    //cout << "{gt} = ";
    //gt.Print( );
    //cout << "\n";
 
    // Evaluate the macroscopic (homogenized) strain.
 
    cMatrix Hr(nstr, ndsp);  // Right
    cMatrix Ht(nstr, ndsp);  // Top
    CalcH2D(Lx,  0.0, Hr);
    CalcH2D(0.0, Ly,  Ht);
    cVector strM(nstr);
    strM  = Hr*gr;
    strM += Ht*gt;
    strM /= RVEvol;
    cout << "{strM_H} = ";
    strM.Print( );
    cout << "\n";
  }
  else  // 3D
  {
    // Evaluate the resulting forces at front, right and top faces.
 
    int ndsp = 3;   // Number of displacements
    int nstr = 6;   // Number of stresses
    double Tol = 1.0e-6*min(Lx, min(Ly, Lz));
    cVector gf(ndsp);  // Front
    cVector gr(ndsp);  // Right
    cVector gt(ndsp);  // Top
    gf.Zero( ); 
    gr.Zero( ); 
    gt.Zero( ); 
    for (int i = 0; i < nnode; i++)
    {
      // Get the current node.
 
      cNode *node = cNode::GetNode(i);
      sNodeCoord p = node->GetCoord( );
 
      // Add to the resulting force of front nodes.
 
      if (abs(p.x - Lx) <= Tol)
        for (int j = 0; j < ndsp; j++) gf[j] += R[i][j];
 
      // Add to the resulting force of right nodes.
 
      if (abs(p.y - Ly) <= Tol)
        for (int j = 0; j < ndsp; j++) gr[j] += R[i][j];
 
      // Add to the resulting force of top nodes.
 
      if (abs(p.z - Lz) <= Tol)
        for (int j = 0; j < ndsp; j++) gt[j] += R[i][j];
    }
 
    // Evaluate the macroscopic (homogenized) strain.
 
    cMatrix Hf(nstr, ndsp);  // Front
    cMatrix Hr(nstr, ndsp);  // Right
    cMatrix Ht(nstr, ndsp);  // Top
    CalcH3D(Lx,  0.0, 0.0, Hf);
    CalcH3D(0.0, Ly,  0.0, Hr);
    CalcH3D(0.0, 0.0, Lz,  Ht);
    cVector strM(nstr);
    strM  = Hf*gf;
    strM += Hr*gr;
    strM += Ht*gt;
    strM /= RVEvol;
    cout << "{strM_H} = ";
    strM.Print( );
    cout << "\n";
  }
}

// -------------------------------------------------------------------------
// Private methods:
//

// =============================== InitDofs ================================
//
// This method activate the free dofs of each node of the FE mesh. This task
// is performed in two steps:
//
// [1] Loops through the elements to activate the appropriate dofs of the
//     connected nodes.
// [2] Deactivate the fixed dofs (given in the input file) which were
//     activated in the previous step.
//
// Proceding this way it is not necessary to explicitly fix, in the input
// file, the nodal dofs that are not used by mesh elements, in order to
// avoid the singularity of the stiffness matrix.
//
void cControl :: InitDofs(void)
{
  int i,j;

  // Activate the nodal dofs according to the connected elements.

  int nelm = cElement :: GetNumElm( );
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);

    // Activate the dofs of element nodes.

    elm->ActivateDof( );
  }

  // Apply the boundary conditions.

  int nsupp = cNode :: GetNumSupp( );
  int dof[8];
  for (i = 0; i < nsupp; i++)
  {
    // Get the current support.

    sNodeSupp supp = cNode :: GetSupp(i);

    // Get the nodal dofs.

    for (j = 0; j < 6; j++) dof[j] = supp.node->GetDof(j);

    // Deactivate the fixed dofs.

    if (supp.dx) dof[0] = 0;
    if (supp.dy) dof[1] = 0;
    if (supp.dz) dof[2] = 0;
    if (supp.rx) dof[3] = 0;
    if (supp.ry) dof[4] = 0;
    if (supp.rz) dof[5] = 0;

    // Set the nodal dofs.

    for (j = 0; j < 6; j++) supp.node->SetDof(j, dof[j]);
  }
}

// =============================== CompEqs =================================
//
// This method sets the equation number of each active dof and compute the
// number of equations of the FE model (number of active dofs).
//
void cControl :: CompEqs(void)
{
  // Compute the nodal equations.

  NumEq = 0;
  int nnode = cNode :: GetNumNode( );
  for (int i = 0; i < nnode; i++)
  {
    // Get the current node.

    cNode *node = cNode :: GetNodeInSolverOrder(i);

    // Nodal constraints.

    sNodeConstr *constraint = cNode :: FindConstraint(node->GetLabel( ));

    // Compute the number of equations and set the nodal dofs.

    if (!constraint) // Node is not linked to other node
    {
      for (int j = 0; j < 8; j++)
      {
        if (node->GetDof(j))        // Free DOF
        {
          NumEq++;
          node->SetDof(j, NumEq);
        }
      }
    }
    else             // Node is linked to a master node
    {
      for (int j = 0; j < 8; j++)
      {
        if (constraint->doflink[j]) // Linked DOF
        {
          node->SetDof(j, constraint->master->GetDof(j));
        }
        else                        // Independent DOF
        {
          if (node->GetDof(j))      // Free DOF
          {
            NumEq++;
            node->SetDof(j, NumEq);
          }
        }
      }
    }
  }
}

// =============================== CompProf ================================
//
// This method computes the profile of the global stiffness matrix. This
// vector indicates the index of the first nonzero element of each line
// of the stiffness matrix.
//
void cControl :: CompProf(void)
{
  int i,j;

  // Alloc and initialize the profile vector.

  Profile = new int[NumEq];
  for (i = 0; i < NumEq; i++) Profile[i] = NumEq + 1;

  // Compute the matrix profile.

  int nelm = cElement :: GetNumElm( );
  int maxdof = cElement :: GetMaxDof( );
  int *elmdof = new int[maxdof];
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);

    // Get the element dofs.

    int ndof = elm->GetNumElmDof( );
    elm->GetElmDofs(elmdof);

    // Compute the lowest element active equation number.

    int min = NumEq + 1;
    for (j = 0; j < ndof; j++)
    {
      if ((elmdof[j]) && (elmdof[j] < min)) min = elmdof[j];
    }

    // Update the system profile.

    for (j = 0; j < ndof; j++)
    {
      if ((elmdof[j]) && (Profile[elmdof[j]-1] >= min))
        Profile[elmdof[j]-1] = min - 1;
    }
  }

  // Release memory.

  delete []elmdof;
}

// ============================== GetPrintIdx ==============================
//
// This method computes the indices of the element response labels in a
// given vector of global labels. It should be used to print the element
// responses in the appropriate position of the output file.
//
//   nel      - number of element labels                               (in)
//   elmlabel - vector of element labels                               (in)
//   ngl      - number of global labels                                (in)
//   glblabel - vector of global labels                                (in)
//   idx      - indices of elm labels in the array of glb labels      (out)
//
void cControl :: GetPrintIdx(int nel, int *elmlabel, int ngl, int *glblabel,
                             int *idx)
{
  for (int i = 0; i < nel; i++)
  {
    for (int l = 0; l < ngl; l++)
    {
      if (elmlabel[i] == glblabel[l])
      {
        idx[i] = l;
        break;
      }
    }
  }
}

// =========================== PrintIntPntStress ===========================

void cControl :: PrintIntPntStress(void)
{
  int i,j,k;

  // Check if there are integration points.

  int maxp = cElement :: GetMaxIntPnt( );
  if (maxp == 0) return;

  // Print the integration point scalar labels.

  out << "%RESULT.CASE.STEP.ELEMENT.GAUSS.SCALAR\n";
  int  ngsl = cElement :: GetNumSclLab( );
  int *vgsl = new int[ngsl];
  cElement :: GetVecSclLab(vgsl);
  out << ngsl << "\n";
  for (i = 0; i < ngsl; i++) out << "'" << VEC_SCL_LAB[vgsl[i]] << "'  ";
  out << "\n\n";

  // Print the integration point scalar data.

  out << "%RESULT.CASE.STEP.ELEMENT.GAUSS.SCALAR.DATA\n";
  int nelm = cElement :: GetNumElm( );
  out << nelm << "\n";

  int *vel = new int[ngsl];
  int *idx = new int[ngsl];
  cVector prn(ngsl);
  cMatrix stress(maxp,ngsl);
  cMatrix strain(maxp,ngsl);
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);
    out << elm->GetLabel( ) << "\n";
    int npt = elm->GetNumIntPnt( );
    out << npt << "\n";
    if (npt == 0) continue;

    // Get element stress labels and the respective indices.

    int nel = elm->GetNumStrCmp( );
    elm->GetStrLabels(vel);
    GetPrintIdx(nel, vel, ngsl, vgsl, idx);

    // Get integration point stresses.

    elm->IntPntStress(strain, stress);

    // Print stresses in the correct order.

    out << showpos << scientific << setprecision(OutPrec);
    
    for (j = 0; j < npt; j++)
    {
      prn.Zero( );
      for (k = 0; k < nel; k++) prn[idx[k]] = stress[j][k];
      for (k = 0; k < ngsl; k++) out << prn[k] << " ";
      out << "\n";
    }

    out << resetiosflags(ios::scientific | ios::showpos);
  }
  out << "\n";

  // Release memory.

  delete []vgsl;
  delete []vel;
  delete []idx;
}

// ============================ PrintNodalStress ===========================

void cControl :: PrintNodalStress(void)
{

  int i,j,k;

  // Print the integration point scalar labels.

  out << "%RESULT.CASE.STEP.ELEMENT.NODAL.SCALAR\n";
  int  maxen = cElement :: GetMaxNode( );
  int  ngsl  = cElement :: GetNumNodSclLab( );
  int *vgsl  = new int[ngsl];
  cElement :: GetVecNodSclLab(vgsl);
  out << ngsl << "\n";
  for (i = 0; i < ngsl; i++) out << "'" << VEC_SCL_LAB[vgsl[i]] << "'  ";
  out << "\n\n";

  // Print the integration point scalar data.

  out << "%RESULT.CASE.STEP.ELEMENT.NODAL.SCALAR.DATA\n";
  int nelm = cElement :: GetNumElm( );
  out << nelm << "\n";

  int *vel = new int[ngsl];
  int *idx = new int[ngsl];
  cVector prn(ngsl);
  cMatrix stress(maxen,ngsl);
//  cMatrix strain(maxen,ngsl);
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);
    out << elm->GetLabel( ) << "\n";

    // Get element stress labels and the respective indices.

    int nel = elm->GetNumNodStrCmp( );
    elm->GetNodStrLabels(vel);
    GetPrintIdx(nel, vel, ngsl, vgsl, idx);

    // Get nodal stresses.

    elm->NodalStress(stress);

    // Print stresses in the correct order.

    out << showpos << scientific << setprecision(OutPrec);
    
    int nen = elm->GetShape( )->GetNumNode( );
    for (j = 0; j < nen; j++)
    {
      prn.Zero( );
      for (k = 0; k < nel; k++) prn[idx[k]] = stress[j][k];
      for (k = 0; k < ngsl; k++) out << prn[k] << " ";
      out << "\n";
    }

    out << resetiosflags(ios::scientific | ios::showpos);
  }
  out << "\n";

  // Release memory.

  delete []vgsl;
  delete []vel;
  delete []idx;

/*
  int i,j,k;

  // Print the element nodal scalar labels.

  out << "%RESULT.CASE.STEP.ELEMENT.NODAL.SCALAR\n";
  int  ngsl = cElement :: GetNumSclLab( );
  int *vgsl = new int[ngsl];
  cElement :: GetVecSclLab(vgsl);
  out << ngsl << "\n";
  for (i = 0; i < ngsl; i++) out << "'" << VEC_SCL_LAB[vgsl[i]] << "'  ";
  out << "\n\n";

  // Print the element nodal scalar data.

  out << "%RESULT.CASE.STEP.ELEMENT.NODAL.SCALAR.DATA\n";
  int nelm = cElement :: GetNumElm( );
  out << nelm << "\n";

  int *vel = new int[ngsl];
  int *idx = new int[ngsl];
  cVector prn(ngsl);

  for (i = 0; i < nelm; i++)
  {
    // Get the current element.

    cElement *elm = cElement :: GetElm(i);
    out << elm->GetLabel( ) << "\n";

    // Get element stress labels and the respective indices.

    int nel = elm->GetNumStrCmp( );
    elm->GetStrLabels(vel);
    GetPrintIdx(nel, vel, ngsl, vgsl, idx);

    // Get element nodal stresses.

    int nen = cElement :: GetMaxNode( );
    cMatrix str(nen, nel);
    elm->NodalStress(str);

    // Print element nodal stresses in the correct order.

    out << showpos << scientific;
    for (j = 0; j < nen; j++)
    {
      prn.Zero( );
      for (k = 0; k < nel; k++) prn[idx[k]] = str[j][k];
      for (k = 0; k < ngsl; k++) out << setw(5) << prn[k] << " ";
      out << "\n";
    }
    out << resetiosflags(ios::showpos | ios::scientific);
  }
  out << "\n";

  // Release memory.

  delete []vgsl;
  delete []vel;
  delete []idx;
 */
}

// ========================== PrintAvrNodalStress ==========================

void cControl :: PrintAvrNodalStress( )
{
  int i,j,k,nl;

  // Evaluate average stresses.
  int     nnode  = cNode :: GetNumNode( );
  int     ngsl   = cElement :: GetNumNodSclLab( );
  cMatrix ncount(nnode,ngsl), nstr(nnode,ngsl);
  ncount.Zero( );
  nstr.Zero( );

  // Print nodal scalar labels.

  out << "%RESULT.CASE.STEP.NODAL.SCALAR\n";
  int  maxen = cElement :: GetMaxNode( );
  int *vgsl  = new int[ngsl];
  cElement :: GetVecNodSclLab(vgsl);
  out << ngsl << "\n";
  for (i = 0; i < ngsl; i++) out << "'" << VEC_SCL_LAB[vgsl[i]] << "'  ";
  out << "\n\n";

  // Evaluate average nodal scalar data.
  int nelm = cElement :: GetNumElm( );

  int *vel = new int[ngsl];
  int *idx = new int[ngsl];
  cMatrix stress(maxen,ngsl);
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.
    cElement *elm = cElement :: GetElm(i);

    // Get element stress labels and the respective indices.
    int nel = elm->GetNumNodStrCmp( );
    elm->GetNodStrLabels(vel);
    GetPrintIdx(nel, vel, ngsl, vgsl, idx);

    // Get nodal stresses.
    elm->NodalStress(stress);

    // Print stresses in the correct order.
    int nen = elm->GetShape( )->GetNumNode( );
    for (j = 0; j < nen; j++)
    {
      nl = elm->GetShape( )->GetNode(j)->GetLabel( );
      for (k = 0; k < nel; k++) 
      {
        nstr(nl-1,idx[k])   += stress[j][k];
        ncount(nl-1,idx[k]) += 1;
      }
    }
  }

  // Average computed results.
  for(i = 0; i < nnode; ++i)
    for(j = 0; j < ngsl; ++j)
      if (ncount(i,j)) nstr(i,j) /= ncount(i,j);

  // Print nodal scalar scalar data.
  out << "%RESULT.CASE.STEP.NODAL.SCALAR.DATA\n";
  out << nnode << "\n";

  //double maxSxx = -1e15; // temporario.
  for (i = 0; i < nnode; i++)
  {
    // Get the current node.
    out << left << setw(5) << i+1 << " ";

    // Get the nodal displacements.
    out << showpos << scientific << setprecision(OutPrec);
    for (j = 0; j < ngsl; j++) out << nstr(i,j) << "  ";
    out << "\n" << resetiosflags(ios::showpos | ios::scientific);
    
   // maxSxx = max(maxSxx,nstr(i,0));
  }
  out << "\n";
  //cout << scientific << setprecision(6) << endl;
  //cout << "Max Stress XX = " << maxSxx << endl; // Remover gambiarra depois!!

  // Release memory.

  delete []vgsl;
  delete []vel;
  delete []idx;
}

// ============================ PrintStrainEnergy ==========================

void cControl :: PrintStrainEnergy(void)
{
  // Alloc the element stiffness matrix (max size).
  int maxdof = cElement :: GetMaxDof( );
  cMatrix Kelm(maxdof, maxdof);
  cVector u(maxdof);
  cVector aux(maxdof);

  // Compute each element strain energy.
  int nelm = cElement :: GetNumElm( );
  cElement* elm = 0;
  double U = 0.0;
  
  for (int i = 0; i < nelm; i++)
  {
    // Get the current element.
    elm = cElement :: GetElm(i);

    // Evaluate the element stiffness matrix.
    elm->StiffMat(Kelm);

    // Get element displacements.
    elm->NodalDispl(u);
    aux = Kelm * u;

    for(int n = 0; n < maxdof; ++n)
      U += u[n] * aux[n]; 
  }
//  cout << scientific << setprecision(10) << endl;
//  cout << "U = " << U << endl; 
}

// ============================ PrintErrNorm ===============================

void cControl :: PrintErrNorm( )
{
  // Auxiliary variables.
  int code    = 0;   // 0 - displacements, 1 - stresses, 2 - Energy.
  double norm = 0.0; // L2 error norm.
  double vol  = 0.0; // model volum.
  int i,j,k;        // Counters.
  sNodeCoord cart;  // Cartesian coordinates of integration points.

  // Check if there are integration points.

  int maxp = cElement :: GetMaxIntPnt( );
  if (maxp == 0) return;

  // Alloc memory.
  int maxnode       = cElement :: GetMaxNode( );
  int  ngsl         = cElement :: GetNumSclLab( );
  double *N         = new double[maxnode];
  sNodeCoord *coord = new sNodeCoord[maxnode];

  // Evaluate L2 error norm in each element.
  cMatrix stress(maxp,ngsl);
  cMatrix strain(maxp,ngsl);
  int nelm = cElement :: GetNumElm( );
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.
    cElement *elm = cElement :: GetElm(i);

    // Get integration point stresses.
    stress.Zero( );
    strain.Zero( );
    elm->IntPntStress(strain,stress);
 
    // Evaluate error function.
    sNatCoord pnt;
    double    w;
    int       nip = elm->GetNumIntPnt( );
    cVector error(nip);
    error.Zero( );
    for(j = 0; j < nip; ++j)
    {
      // Get analytical solution.
      elm->GetIntPnt(j,pnt,w);
      elm->GetShape( )->Evaluate(pnt,N,coord,cart);

      double sigT[9];
      cField :: GetInstance( )->Sig(cart.x,cart.y,cart.z,sigT);
      sigT[2] = sigT[1];
      sigT[1] = sigT[4];

      for(k = 0; k < ngsl; ++k)
        error[j] += pow(stress[j][k]-sigT[k],2.0);       
        //error[j] += strain[j][k]*stress[j][k];       
    }

    // Integrate error function.
    norm += elm->IntFunc(error,vol);
  }
  norm = sqrt(norm);

  cout << "Stress error norm: " << norm << " " << norm/vol << ", vol " << vol  << endl;
  out << "\n";
}



static void getTrianQuadSampl(const int &level, int &npt, sNatCoord *&coord,
int &npb, int *&bind, int &nt, int *&tind)
{
  // Compute tesselation triangles.
  double delta = 1.0/(level);

  int cont1 = 0;
  int cont2 = 0;
  npt = (level+1) * (level+2) / 2;
  coord = new sNatCoord[npt];
  npb  = level*3; 
  bind = new int [npb];

  for(int j = 0; j < level+1; ++j)
    for(int i = 0; i < level-j+1; ++i, ++cont1)
    {
      coord[cont1].r = i*delta;
      coord[cont1].s = j*delta;
      if (j == 0 || i == (level-j))
        bind[cont2++] = cont1+1;  
      else if (i == 0)
        bind[npb-j] = cont1+1;
    }

  // Store triangulation index.
  nt    = level * level;
  tind  = new int [nt*3];
  cont1 = 0;
  cont2 = 0;
  for(int j = 0; j < level; ++j,++cont2)
    for(int i = 0; i < level-j; ++i,++cont2)
    {
      // First triangle.
      //tind[cont1++] = cont2+1;
      //tind[cont1++] = cont2+2;
      //tind[cont1++] = cont2+level-j+2;
      tind[cont1++] = ((2*level+3-j) * j)/2     + i + 1;
      tind[cont1++] = ((2*level+3-j) * j)/2     + (i+1) + 1;
      tind[cont1++] = ((2*level+2-j) * (j+1))/2 + i + 1;

      // Second triangle.
      if (i < level-j-1)
      {
        // Evaluate second triangle vertices.
        //tind[cont1++] = cont2+level-j+2;
        //tind[cont1++] = cont2+2;
        //tind[cont1++] = cont2+level-j+3;
        tind[cont1++] = ((2*level+2-j) * (j+1))/2 + i + 1;
        tind[cont1++] = ((2*level+3-j) * j)/2     + (i+1) + 1;
	tind[cont1++] = ((2*level+2-j) * (j+1))/2 + (i+1) + 1;
      }
    }
}

void cControl :: PrintStrErrNormScalVis(int code)
{
  int nelm = cElement :: GetNumElm( );
  int maxpnt = 10000;
  int maxelmpnt = round(double(maxpnt)/nelm);
  int level = std::max(sqrt(maxelmpnt),1.0);

  // Evaluate number of loads.
  if (1)
  {
    int nload = 0;
    for (cLoad *load = cLoad::GetHead( ); load != 0; load = load->GetNext( ))
      nload++;

    cout << "num load " << nload << endl;
    int cl = round(log10(nload)/log10(2));
    level  = pow(2.0,8-cl);
    level  = (level < 1) ? 20 : level;
    cout << "level " << level << endl;
  }


  string filename = fname + ".visdat"; 
  ofstream outvis(filename.c_str( ));
 
  // Get element sample points.
  int npt;
  int npb;
  int *bind;
  int nt;
  int *tind;
  sNatCoord *ncoords;
  getTrianQuadSampl(level,npt,ncoords,npb,bind,nt,tind);
  
  int i,j,k;
  int ord            = ErrIntOrd;
  int NumIntPnt      = 0;
  int maxnode        = cElement :: GetMaxNode( );
  double norm        = 0.0;
  double vol         = 0.0;
  double *N          = new double[maxnode];
  eTopType ttype;
  sNodeCoord *coord  = new sNodeCoord[maxnode];
  sNatCoord   *dN    = new sNatCoord[maxnode]; 

  // Loop over each element.
  int ngsl = cElement :: GetNumSclLab( );
  cMatrix stress(npt,ngsl);
  cMatrix strain(npt,ngsl);

  outvis << "%REGION.NUMBER\n" << nelm << "\n\n";
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.
    cElement *elm = cElement :: GetElm(i);

    outvis << "%REGION\n"; 
    outvis << "1 " << i+1 << " " << npt << "\n";

    // Write element vertex coords.
    for(int j = 0; j < npt; ++j)
    {
      sNodeCoord c;
      elm->GetShape( )->Evaluate(ncoords[j],N,coord,c);
      outvis << c.x << " " << c.y << " " << c.z << "\n";
    }
    
    // Write boundary information.
    outvis << "1 'p' " << npb << endl;
    
    // Write boundary vertex incidence.
    for(int j = 0; j < npb; ++j)
      outvis << bind[j] << " ";
    outvis << "\n";

    // Write element triangulation.
    outvis << "'i' " << nt << " ";
    for(int j = 0; j < nt*3; ++j)
      outvis << tind[j] << " ";
    outvis << "\n";

    // Write element results.
    //outvis << "5 'Jacobian Triangle Shape' 'Stress Norm Error' 'STRESS_XX' 'STRESS_YY' 'STRESS_XY'\n\n";;
    outvis << "5 'Jacobian Triangle Shape' 'Stress Norm Error' 'DISPLACEMENT_MAGNITUDE' 'STRESS_YY' 'STRESS_XY'\n\n";;

    // Evaluate integration point strains and stresses.
    stress.Zero( );
    strain.Zero( );
    elm->PntStress(ncoords,npt,strain,stress,true);

    // Write element results.
    outvis << "%VERTEX.RESULT\n5 "; 
    outvis << "'Stress Norm Error' 1 " << i+1 << " ";

    // Evaluate error function.
    for(j = 0; j < npt; ++j)
    {
      // Get analytical solution.
      sNodeCoord pos;
      elm->GetShape( )->Evaluate(ncoords[j],N,coord,pos);
    
      double error     = 0.0;
      double FieldNorm = 0.0;
      double NumNorm   = 0.0;
      double strs[6];
      cField :: GetInstance( )->Stress(pos.x,pos.y,pos.z,strs);
      strs[2] = strs[3];
    
      if (code == 0)
      {
        for(k = 0; k < ngsl; ++k)
	{
	  NumNorm   += pow(stress[j][k]-strs[k],2.0);
          FieldNorm += strs[k]*strs[k];
	}
	error = sqrt(NumNorm) / sqrt(FieldNorm);
      }
      else
      {
        double strn[6];
        cField :: GetInstance( )->Strain(pos.x,pos.y,pos.z,strn);
        strn[2] = strn[3];
    
        for(k = 0; k < ngsl; ++k)
          error += (strn[k]-strain[j][k]) * (strs[k]-stress[j][k]);
      }
      outvis << error << " "; 
    }

    // Print Jacobian triangle shape.
    outvis << "\n'Jacobian Triangle Shape' 1 " << i+1 << " ";
    for(j = 0; j < npt; ++j)
    {
      double metric = 0.0;

      // Evaluate Shape function derivatives.
      elm->GetShape( )->DrvShpRST(ncoords[j],dN);

      // Evaluate Jacobian matrix.
      int nn = elm->GetShape( )->GetNumNode( );
      elm->GetShape( )->NodalCoord(coord);

      double J[4] = {0,0,0,0};
      for (int k = 0; k < nn; k++)
      {
        J[0] += dN[k].r*coord[k].x;
        J[1] += dN[k].r*coord[k].y;
        J[2] += dN[k].s*coord[k].x;
        J[3] += dN[k].s*coord[k].y;
      } 

      // Check if jacobian is valid.
      if ((J[0]*J[3] - J[1]*J[2]) < 1e-8)
        metric = 0.0;
      else
      {
        // Compute metric tensor.
        double m11, m22, m12, det;
        m11 = J[0] * J[0] + J[1] * J[1];
        m22 = J[2] * J[2] + J[3] * J[3];
        m12 = J[0] * J[2] + J[1] * J[3];
        
        det = m11 * m22 - m12 * m12;
        metric  = sqrt(3.0 * det) / (m11 + m22 - m12);
      }
      outvis << metric << " "; 
    }

    // Print Stress XX.
    //outvis << "\n'STRESS_XX' 1 " << i+1 << " ";
    //for(j = 0; j < npt; ++j) outvis << stress[j][0] << " ";

    // Print displacement magnitude.
    outvis << "\n'DISPLACEMENT_MAGNITUDE' 1 " << i+1 << " ";
    for(j = 0; j < npt; ++j)
    {
      elm->GetShape( )->ShpFunc(ncoords[j],N);

      double u = 0.0, v = 0.0;
      int nn = elm->GetShape( )->GetNumNode( );
      for(int k = 0; k < nn; ++k)
      {
        u += elm->GetShape( )->GetNode(k)->GetDispl(0) * N[k];	      
        v += elm->GetShape( )->GetNode(k)->GetDispl(1) * N[k];	      
      }      
      outvis << sqrt(u*u+v*v) << " ";
    }

    // Print Stress YY.
    outvis << "\n'STRESS_YY' 1 " << i+1 << " ";
    for(j = 0; j < npt; ++j) outvis << stress[j][1] << " ";

    // Print Stress YY.
    outvis << "\n'STRESS_XY' 1 " << i+1 << " ";
    for(j = 0; j < npt; ++j) outvis << stress[j][2] << " ";

    outvis << "\n\n";
  }
  outvis.close( );

  // Release memory.
  delete [] coord;
  delete [] dN;
  delete [] ncoords;
  delete [] bind;
  delete [] tind;
}

void cControl :: PrintStrErrNorm(int code)
{
  int i,j,k;
  int ord            = ErrIntOrd;
  int NumIntPnt      = 0;
  int maxnode        = cElement :: GetMaxNode( );
  double norm        = 0.0;
  double vol         = 0.0;
  double *N          = new double[maxnode];
  eTopType ttype;
  cIntPoint *IntPnt  = 0;
  sNatCoord *ncoords = 0;
  sNodeCoord *coord  = new sNodeCoord[maxnode];

  // Loop over each element.
  int nelm = cElement :: GetNumElm( );
  int ngsl = cElement :: GetNumSclLab( );
  cMatrix stress;
  cMatrix strain;
  for (i = 0; i < nelm; i++)
  {
    // Get the current element.
    cElement *elm = cElement :: GetElm(i);

    // Create quadrature points.
    if (i == 0 || ttype != elm->GetShape( )->GetTopologyType( ))
    {
      // Release memory of old integration points.
      if (ncoords) delete [] ncoords;
      if (IntPnt)  delete [] IntPnt;

      int order[3]      = {ord,ord,ord};
      ttype             = elm->GetShape( )->GetTopologyType( );
      eQuadType qtype   = static_cast<eQuadType>(elm->GetIntType( ));
      eTopType  stype   = elm->GetShape( )->GetTopologyType( );
                 IntPnt = cIntPoint::CreateIntPoints(stype,qtype,order,&NumIntPnt);
      
	  // Store coordinates.
      ncoords = new sNatCoord [NumIntPnt];
      for(int i = 0; i < NumIntPnt; ++i) 
        ncoords[i] = IntPnt[i].GetCoord( ); 

      // Resize matrices.
      if (stress.NRow( ) < NumIntPnt)
      {
        stress.Resize(NumIntPnt,ngsl);
        strain.Resize(NumIntPnt,ngsl);
      }
    }

    // Evaluate integration point strains and stresses.
    stress.Zero( );
    strain.Zero( );
    elm->PntStress(ncoords,NumIntPnt,strain,stress,true);
 
    // Evaluate error function.
    cVector error(NumIntPnt);
    error.Zero( );
    for(j = 0; j < NumIntPnt; ++j)
    {
      // Get analytical solution.
      sNodeCoord pos;
      elm->GetShape( )->Evaluate(IntPnt[j].GetCoord( ),N,coord,pos);

      double strs[6];
      cField :: GetInstance( )->Stress(pos.x,pos.y,pos.z,strs);
      strs[2] = strs[3];

      if (code == 0)
        for(k = 0; k < ngsl; ++k)
          error[j] += pow(stress[j][k]-strs[k],2.0);
      else
      {
        double strn[6];
        cField :: GetInstance( )->Strain(pos.x,pos.y,pos.z,strn);
        strn[2] = strn[3];

        for(k = 0; k < ngsl; ++k)
          error[j] += (strn[k]-strain[j][k]) * (strs[k]-stress[j][k]);
      }
    }

    // Integrate error function.
    norm += elm->IntFunc(NumIntPnt,IntPnt,error,vol);
  }
  norm = sqrt(norm);

  out << scientific << setprecision(6);
  cout << scientific << setprecision(6);

  cout << ((code == 0) ? "L2 Stress Norm: " : "L2 Energy Norm: ") << norm << endl; 
  out  << ((code == 0) ? "%L2.STRESS.NORM\n" : "%L2.ENERGY.NORM\n") << norm << "\n\n";

  out << resetiosflags(ios::scientific | ios::showpos);
  cout << resetiosflags(ios::scientific | ios::showpos);

  // Release memory.
  delete [] IntPnt;
  delete [] ncoords;
}

// =========================== PrintPatDispl ===============================

void cControl :: PrintPatDispl( )
{
  if (cShapeIGA :: GetPatPntVec( ).empty( ))
    return;

  // Auxiliary data.
  sNatCoord       p;
  cPatch         *pat;
  cShape         *shp;
  const sPatPnt  *pnt;
  int             elmid;
  double          Displ[6];
  double          ElmRefCoord[3] = {0.0,0.0,0.0};

  // Alloc memory for shape functions.
  int maxnode = cElement :: GetMaxNode( );
  double *N   = new double [maxnode];
  
  // Printing output tag.
  out << "%RESULT.PATCH.DISPLACEMENT\n";
  out << cShapeIGA :: GetPatPntVec( ).size( ) << "\n";
 
  // Loop over all patch displacement points.
  for(unsigned int i = 0; i < cShapeIGA :: GetPatPntVec( ).size( ); ++i)
  {
    // Get patch pointer.
    pnt = &(cShapeIGA :: GetPatPntVec( )[i]);
    pat = cPatch :: FindPatch(pnt->patid);

    // Get shape containing this point.
    elmid = pat->getSpanID(pnt->parcoord[0],pnt->parcoord[1],pnt->parcoord[2]);

    // Test if exist shape in that patch.
    const PatElmStdMap &Map = cShapeIGA :: GetPatElmMap( );
    bool skip = false;

    if (Map.count(pat) == 0)
      skip = true;
    else if (Map.find(pat)->second[elmid] == 0)
      skip = true;

    if (skip)
    {
      // Print displacements.
      out << left << setw(4) << i+1 << " ";
      out << showpos << scientific << setprecision(OutPrec);
      for(int j = 0; j < 5; ++j)
        out << 0.0 << " ";
      out << 0.0 << "\n";
      out << resetiosflags(ios::showpos | ios::scientific);
 
      continue;
    }

    // Get element shape.
    shp = Map.find(pat)->second[elmid]->GetShape( );

    // Get shape span limits.
    double SpanLim[3][2];
    pat->getSpanLimits(elmid,SpanLim);

    // Compute parametric coordinates in element reference [-1,1].
    for(int pvar = 0; pvar < pat->getNumParVar( ); ++pvar)
    {
      ElmRefCoord[pvar]  = pnt->parcoord[pvar] - SpanLim[pvar][0];
      ElmRefCoord[pvar] *= 2.0/(SpanLim[pvar][1] - SpanLim[pvar][0]);
      ElmRefCoord[pvar] -= 1.0;
    }

    // Setup parametric coords values.
    p.r = ElmRefCoord[0];
    p.s = ElmRefCoord[1];
    p.t = ElmRefCoord[2];

    // Evaluate shape functions.
    shp->ShpFunc(p,N);

    // Evaluate displacements.
    for(int j = 0; j < 6; ++j)
      Displ[j]  = 0;

    for(int j = 0; j < shp->GetNumNode( ); ++j)
      for(int k = 0; k < 6; ++k)
        Displ[k] += shp->GetNode(j)->GetDispl(k) * N[j];

    // Evandro - Temporario!
    cout << "IGA_DISP-LOADFAC  ";
    cout << showpos << scientific << setprecision(OutPrec);
    for(int j = 0; j < 6; ++j)
      cout << Displ[j] << " ";
    cout << scientific << setprecision(OutPrec) << TotFactor << "\n\n";
    cout << resetiosflags(ios::showpos | ios::scientific);

    // Print displacements.
    out << left << setw(4) << i+1 << " ";
    out << showpos << scientific << setprecision(OutPrec);
    for(int j = 0; j < 5; ++j)
      out << Displ[j] << " ";
    out << Displ[5] << "\n";
    out << resetiosflags(ios::showpos | ios::scientific);
  }
  out << "\n";

  // Release memory.
  delete [] N;
}

// ============================ PrintIGAStress =============================

void cControl :: PrintPatStr( )
{
  if (cShapeIGA :: GetPatPntVec( ).size( ) == 0)
    return;

  // Auxiliary data.
  sNatCoord      p;
  cPatch        *pat;
  cElement      *elm;
  const sPatPnt *pnt;
  int            elmid;
  int            nel = 0;
  int            i,k,id,pvar;

  double ElmRefCoord[3] = {0.0,0.0,0.0};

  // Print the patch scalar labels.

  out << "%RESULT.CASE.STEP.PATCH.SCALAR\n";
  int  ngsl = cElement :: GetNumSclLab( );
  int *vgsl = new int[ngsl];
  cElement :: GetVecSclLab(vgsl);
  out << ngsl << "\n";
  for (int i = 0; i < ngsl; i++) out << "'" << VEC_SCL_LAB[vgsl[i]] << "'  ";
  out << "\n\n";

  // Print the patch scalar data.

  out << "%RESULT.CASE.STEP.PATCH.SCALAR.DATA\n";
  int npt = cShapeIGA :: GetPatPntVec( ).size( );
  out << npt << "\n";

  int *vel = new int[ngsl];
  int *idx = new int[ngsl];
  cVector prn(ngsl); 
  cMatrix stress(1,ngsl), strsum(1,ngsl);
  cMatrix strain(1,ngsl);

  int numelm;
  int inds[8];

  for (i = 0; i < npt; i++)
  {
    // Get patch pointer.
    pnt = &(cShapeIGA :: GetPatPntVec( )[i]);
    pat = cPatch :: FindPatch(pnt->patid);

    // Get shape containing this point.
    pat->getSpanID(pnt->parcoord[0],pnt->parcoord[1],pnt->parcoord[2],numelm,inds);
    out << left << setw(4) << i+1;

    // Clear stress sum vector.
    if (PatStrSmt)
    {
      strsum.Zero( );
      out << "\n"; 
    }
    else 
      out << "   " << numelm << "\n";

    // Loop over each element.
    for(id = 0; id < numelm; ++id)
    {
      elmid = inds[id];
    
      // Test if the target element exists.
      const PatElmStdMap &Map = cShapeIGA :: GetPatElmMap( );
      if (Map.count(pat) == 0)
      {
        cout << " 'Null Element' \n";
        continue;
      }
      else if (Map.find(pat)->second[elmid] == 0)
      {
        cout << " 'Null Element' \n";
        continue;
      }
      elm   = Map.find(pat)->second[elmid];
      
      // Get shape span limits.
      double SpanLim[3][2];
      pat->getSpanLimits(elmid,SpanLim);
      
      // Compute parametric coordinates in element reference [-1,1].
      for(pvar = 0; pvar < pat->getNumParVar( ); ++pvar)
      {
        ElmRefCoord[pvar]  = pnt->parcoord[pvar] - SpanLim[pvar][0];
        ElmRefCoord[pvar] *= 2.0/(SpanLim[pvar][1] - SpanLim[pvar][0]);
        ElmRefCoord[pvar] -= 1.0;
      }
      
      // Setup parametric coords values.
      p.r = ElmRefCoord[0];
      p.s = ElmRefCoord[1];
      p.t = ElmRefCoord[2];
      
      // Get element stress labels and the respective indices.
      nel = elm->GetNumStrCmp( );
      elm->GetStrLabels(vel);
      cControl :: GetControl( )->GetPrintIdx(nel, vel, ngsl, vgsl, idx);
      
      // Get element stresses at the requested parametric point.
      elm->PntStress(&p,1,strain,stress);
      
      // Accumulate stresses.
      if (PatStrSmt)
        strsum += stress;
      else
      {
        // Print patch parametric point stress.
        out << left << setw(4) << elm->GetLabel( ) << " ";
        out << showpos << scientific;
        prn.Zero( );
        for (k = 0; k < nel; k++) prn[idx[k]] = stress[0][k];
        for (k = 0; k < ngsl; k++) out << setw(5) << prn[k] << " ";
        out << "\n";
        out << resetiosflags(ios::showpos | ios::scientific);
      }
    }

    if (PatStrSmt)
    {
      // Print the smoothed results.
      out << showpos << scientific;
      prn.Zero( );
      for (k = 0; k < nel; k++) prn[idx[k]] = strsum[0][k];
      for (k = 0; k < ngsl; k++) out << setw(5) << prn[k]/numelm << " ";
      out << "\n";
      out << resetiosflags(ios::showpos | ios::scientific);
    }
  }
  out << "\n";

  // Release memory.

  delete []vgsl;
  delete []vel;
  delete []idx;

}

void cControl :: PrintEqvMshDispl( )
{
  cShapeIGA :: PrintEqvMshDispl(outem); 
}

// ============================ PrintEqvMshNodalStress =====================

void cControl :: PrintEqvMshNodalStress( )
{
  int i,j,k,e;

  // Print the integration point scalar labels.

  outem << "%RESULT.CASE.STEP.ELEMENT.NODAL.SCALAR\n";
  int  maxen  = 20; // Eqv Elem with most nodes is brick20.
       maxen *= max(cShapeIGA :: GetEqvMshDisc(0)+1,1); // Multiplied by user
       maxen *= max(cShapeIGA :: GetEqvMshDisc(1)+1,1); // defined discretization.
       maxen *= max(cShapeIGA :: GetEqvMshDisc(2)+1,1);
  int  ngsl  = cElement :: GetNumSclLab( );
  int *vgsl  = new int[ngsl];
  cElement :: GetVecSclLab(vgsl);
  outem << ngsl << "\n";
  for (i = 0; i < ngsl; i++) outem << "'" << VEC_SCL_LAB[vgsl[i]] << "'  ";
  outem << "\n\n";

  int GblEqvElmID = 0;

  // Print the element nodal scalar data.
  const ElmStdVec &elmvec = cShapeIGA :: GetParentElmVec( );

  outem << "%RESULT.CASE.STEP.ELEMENT.NODAL.SCALAR.DATA\n";
  int nelm = elmvec.size( );
  outem << cShapeIGA :: GetNumEqvElm( ) << "\n";

  int *vel = new int[ngsl];
  int *idx = new int[ngsl];
  cVector prn(ngsl);
  cMatrix stress(maxen,ngsl);
  cMatrix strain(maxen,ngsl);
  for (i = 0; i < nelm; i++)
  {
    // Get element stress labels and the respective indices.

    int nel = elmvec[i]->GetNumStrCmp( );
    elmvec[i]->GetStrLabels(vel);
    GetPrintIdx(nel, vel, ngsl, vgsl, idx);

    // Get equivalent element type.
    eTopType t = elmvec[i]->GetShape( )->GetTopologyType( );

    // Get
    sEqvElmQuad *Quad = &(cShapeIGA :: GetEqvElmPnt(t));
    
    // Get all equivalent element stresses.  

    elmvec[i]->PntStress(&(Quad->coords[0]),Quad->NumPnt,strain,stress);

    // Print stresses in the correct order.

    int ne  = Quad->NumElm;
    for(e = 0; e < ne; ++e)
    {
      int nen = Quad->NumElmPnt;

      outem << ++GblEqvElmID << "\n";
      outem << showpos << scientific << setprecision(OutPrec);
      for (j = 0; j < nen; j++)
      {
        prn.Zero( );
        for (k = 0; k < nel; k++) prn[idx[k]] = stress[j][k];
        for (k = 0; k < ngsl; k++) outem << prn[k] << " ";
        outem << "\n";
      }
      outem << noshowpos;

      outem << resetiosflags(ios::scientific | ios::showpos);
    }
  }
  outem << "\n";

  // Release memory.

  delete []vgsl;
  delete []vel;
  delete []idx;
}


// ============================ RegisterDestroyFunc ========================

void cControl :: RegisterDestroyFunc(void (*func) (void))
{
  // Register function pointer.
  DestroyFunc.push(func);
}

// ======================================================= End of file =====
