// -------------------------------------------------------------------------
// input.cpp - input of model data.
// -------------------------------------------------------------------------
// Created:      14-Nov-2000    Evandro Parente Junior
//
// Modified:     28-Nov-2000    Evandro Parente Junior
//               Input of control data.
//
// Modified:     06-Nov-2001    Evandro Parente Junior
//               Input of displacement control data.
//
// Modified:     06-Oct-2013    Evandro Parente Junior
//               Elimination of Thickness and redefinition of Section input.
//
// Modified:     08-Apr-2015     Elias Saraiva Barroso
//               Input/output operations using c++ streams.
//
// Modified:     11-Jun-2015     Elias Saraiva Barroso
//               Change the implementation of InpCopyData function: variable 
//               used to read each line (char line[BUFSIZ]) has been changed
//               to string object due to bugs when the line size is very big.
//
// -------------------------------------------------------------------------

#include <iostream>
#include <stdlib.h>
#include <string.h>

using namespace std;

#include "input.h"
#include "ctrl.h"
#include "ctrlpath.h"
#include "ctrldsp.h"
#include "ctrlarclen.h"
#include "ctrlnewmark.h"
#include "ctrlgenalpha.h"
#include "ctrlhht.h"
#include "node.h"
#include "element.h"
#include "elmdsh.h"
#include "material.h"
#include "section.h"
#include "secbar.h"
#include "spring.h"
#include "sprprop.h"
#include "load.h"
#include "prescdof.h"
#include "timefunc.h"
#include "heat.h"
#include "utl.h"
#include "gblvar.h"
#include "shpiga.h"
#include "shpshell.h"
#include "field.h"
#include "eig.h"

// -------------------------------------------------------------------------
// Local definitions:
//

// Structure storing a data label (see the Neutral File Documentation) and
// a pointer to the function which reads the associated data.
typedef struct
{
  const char *label;   // String with one data label of the neutral file.
  void (*func)(void);  // Function to read the data linked to the label.
} sInpData;

// -------------------------------------------------------------------------
// Local variables:
//

// List of the data labels that can be read by the program.
// These labels must be kept alphabetically sorted !!!
static sInpData TabKey[] = {
{ "BEAM.SECTION.ORIENTATION",            cElement :: ReadBeamOrt },
{ "CONTROL.POINT.WEIGHT",                cNode    :: ReadCtrlPntW },
{ "ELEMENT",                             cElement :: ReadNumElm },
{ "ELEMENT.AXISYM.L6_INF",               cElement :: ReadAxiSymL6Inf },
{ "ELEMENT.AXISYM.Q4",                   cElement :: ReadAxiSymQ4 },
{ "ELEMENT.AXISYM.Q4TL",                 cElement :: ReadAxiSymQ4TL },
{ "ELEMENT.AXISYM.Q8",                   cElement :: ReadAxiSymQ8 },
{ "ELEMENT.AXISYM.Q8TL",                 cElement :: ReadAxiSymQ8TL },
{ "ELEMENT.AXISYM.Q9",                   cElement :: ReadAxiSymQ9 },
{ "ELEMENT.AXISYM.Q9TL",                 cElement :: ReadAxiSymQ9TL },
{ "ELEMENT.AXISYM.T3",                   cElement :: ReadAxiSymT3 },
{ "ELEMENT.AXISYM.T3TL",                 cElement :: ReadAxiSymT3TL },
{ "ELEMENT.AXISYM.T6",                   cElement :: ReadAxiSymT6 },
{ "ELEMENT.AXISYM.T6TL",                 cElement :: ReadAxiSymT6TL },
{ "ELEMENT.BRICK20",                     cElement :: ReadBrick20 },
{ "ELEMENT.BRICK20TL",                   cElement :: ReadBrick20TL },
{ "ELEMENT.BRICK8",                      cElement :: ReadBrick8 },
{ "ELEMENT.BRICK8TL",                    cElement :: ReadBrick8TL },
{ "ELEMENT.BSPLINE.SOLID",               cElement :: ReadSolidBsp },
{ "ELEMENT.BSPLINE.SOLIDTL",             cElement :: ReadSolidBspTL },
{ "ELEMENT.DONNELLSHELL.BSPLINE.SURFACE",cElement :: ReadDonnellShellBsp },
{ "ELEMENT.DONNELLSHELL.BSPLINE.SURFACETL",cElement :: ReadDonnellShellBspTL },
{ "ELEMENT.DONNELLSHELL.Q8",             cElement :: ReadDonnellShellQ8 },
{ "ELEMENT.DONNELLSHELL.Q8TL",           cElement :: ReadDonnellShellQ8TL },
{ "ELEMENT.FORCES.OUTPUT.CONVENTION",    cElement :: ReadOutputConv },
{ "ELEMENT.FRAME3D",                     cElement :: ReadFrame3D },
{ "ELEMENT.FRAME3DCR",                   cElement :: ReadFrame3DCR },
{ "ELEMENT.FRAME3DCRTL",                 cElement :: ReadFrame3DCRTL },
{ "ELEMENT.FRAME3DTL",                   cElement :: ReadFrame3DTL },
{ "ELEMENT.GRID",                        cElement :: ReadGrid },     
//{ "ELEMENT.INTERFACE.BSPLINE.SURFACE",   cElement :: ReadInterfLineBsp },
{ "ELEMENT.HSDTPLATE.BSPLINE.SURFACE",   cElement :: ReadHSDTPlateBsp },
{ "ELEMENT.HSDTPLATE.BSPLINE.SURFACETL", cElement :: ReadHSDTPlateBspTL },	
{ "ELEMENT.INTERFACE.L2",                cElement :: ReadInterfL2 },
{ "ELEMENT.INTERFACE.L3",                cElement :: ReadInterfL3 },
{ "ELEMENT.INTERFACE.LINE.BSPLINE",      cElement :: ReadInterfLineBsp },
{ "ELEMENT.INTERFACE.Q4",                cElement :: ReadInterfQ4 },
{ "ELEMENT.INTERFACE.Q8",                cElement :: ReadInterfQ8 },
{ "ELEMENT.PLFRAME",                     cElement :: ReadPlFrame },
{ "ELEMENT.PLFRAMECR1",                  cElement :: ReadPlFrameCR1 },
{ "ELEMENT.PLFRAMECR2",                  cElement :: ReadPlFrameCR2 },
{ "ELEMENT.PLFRAMETL1",                  cElement :: ReadPlFrameTL1 },
{ "ELEMENT.PLFRAMETL2",                  cElement :: ReadPlFrameTL2 },
{ "ELEMENT.PLSTRAIN.BEZIER.SURFACE",     cElement :: ReadPlStrainBezSurf },
{ "ELEMENT.PLSTRAIN.BEZIER.SURFACETL",   cElement :: ReadPlStrainBezSurfTL },
{ "ELEMENT.PLSTRAIN.BEZIER.TRIANGLE",    cElement :: ReadPlStrainBezTri },
{ "ELEMENT.PLSTRAIN.BEZIER.TRIANGLETL",  cElement :: ReadPlStrainBezTriTL },
{ "ELEMENT.PLSTRAIN.BSPLINE.SURFACE",    cElement :: ReadPlStrainBsp },
{ "ELEMENT.PLSTRAIN.BSPLINE.SURFACETL",  cElement :: ReadPlStrainBspTL },
{ "ELEMENT.PLSTRAIN.Q4",                 cElement :: ReadPlStrainQ4 },
{ "ELEMENT.PLSTRAIN.Q4TL",               cElement :: ReadPlStrainQ4TL },
{ "ELEMENT.PLSTRAIN.Q8",                 cElement :: ReadPlStrainQ8 },
{ "ELEMENT.PLSTRAIN.Q8TL",               cElement :: ReadPlStrainQ8TL },
{ "ELEMENT.PLSTRAIN.Q9",                 cElement :: ReadPlStrainQ9 },
{ "ELEMENT.PLSTRAIN.Q9TL",               cElement :: ReadPlStrainQ9TL },
{ "ELEMENT.PLSTRAIN.T3",                 cElement :: ReadPlStrainT3 },
{ "ELEMENT.PLSTRAIN.T3TL",               cElement :: ReadPlStrainT3TL },
{ "ELEMENT.PLSTRAIN.T6",                 cElement :: ReadPlStrainT6 },
{ "ELEMENT.PLSTRAIN.T6TL",               cElement :: ReadPlStrainT6TL },
{ "ELEMENT.PLSTRESS.BEZIER.SURFACE",     cElement :: ReadPlStressBezSurf },
{ "ELEMENT.PLSTRESS.BEZIER.SURFACETL",   cElement :: ReadPlStressBezSurfTL },
{ "ELEMENT.PLSTRESS.BEZIER.TRIANGLE",    cElement :: ReadPlStressBezTri },
{ "ELEMENT.PLSTRESS.BEZIER.TRIANGLETL",  cElement :: ReadPlStressBezTriTL },
{ "ELEMENT.PLSTRESS.BSPLINE.SURFACE",    cElement :: ReadPlStressBsp },
{ "ELEMENT.PLSTRESS.BSPLINE.SURFACETL",  cElement :: ReadPlStressBspTL },
{ "ELEMENT.PLSTRESS.Q4",                 cElement :: ReadPlStressQ4 },
{ "ELEMENT.PLSTRESS.Q4TL",               cElement :: ReadPlStressQ4TL },
{ "ELEMENT.PLSTRESS.Q8",                 cElement :: ReadPlStressQ8 },
{ "ELEMENT.PLSTRESS.Q8TL",               cElement :: ReadPlStressQ8TL },
{ "ELEMENT.PLSTRESS.Q9",                 cElement :: ReadPlStressQ9 },
{ "ELEMENT.PLSTRESS.Q9TL",               cElement :: ReadPlStressQ9TL },
{ "ELEMENT.PLSTRESS.T3",                 cElement :: ReadPlStressT3 },
{ "ELEMENT.PLSTRESS.T3TL",               cElement :: ReadPlStressT3TL },
{ "ELEMENT.PLSTRESS.T6",                 cElement :: ReadPlStressT6 },
{ "ELEMENT.PLSTRESS.T6TL",               cElement :: ReadPlStressT6TL },
{ "ELEMENT.PLTRUSS",                     cElement :: ReadPlTruss },
{ "ELEMENT.PLTRUSSCR",                   cElement :: ReadPlTrussCR },
{ "ELEMENT.PLTRUSSTL",                   cElement :: ReadPlTrussTL },
{ "ELEMENT.SHALLOWSHELL.BEZIER.SURFACE",  cElement :: ReadShallowShellBezSurf },
{ "ELEMENT.SHALLOWSHELL.BEZIER.SURFACETL", cElement :: ReadShallowShellBezSurfTL },
{ "ELEMENT.SHALLOWSHELL.BEZIER.TRIANGLE", cElement :: ReadShallowShellBezTrian },
{ "ELEMENT.SHALLOWSHELL.BEZIER.TRIANGLETL", cElement :: ReadShallowShellBezTrianTL },
{ "ELEMENT.SHALLOWSHELL.BSPLINE.SURFACE",   cElement :: ReadShallowShellBsp },
{ "ELEMENT.SHALLOWSHELL.BSPLINE.SURFACETL", cElement :: ReadShallowShellBspTL },
{ "ELEMENT.SHALLOWSHELL.Q8",             cElement :: ReadShallowShellQ8 },
{ "ELEMENT.SHALLOWSHELL.Q8TL",           cElement :: ReadShallowShellQ8TL },
{ "ELEMENT.SHALLOWSHELL.T6",             cElement :: ReadShallowShellT6 },
{ "ELEMENT.SHALLOWSHELL.T6TL",           cElement :: ReadShallowShellT6TL },
{ "ELEMENT.SHELL.BEZIER.SURFACE",        cElement :: ReadShellBez },
{ "ELEMENT.SHELL.BEZIER.SURFACETL",      cElement :: ReadShellBezTL },
{ "ELEMENT.SHELL.BEZIER.TRIANGLE",       cElement :: ReadShellBezTrian },
{ "ELEMENT.SHELL.BEZIER.TRIANGLETL",     cElement :: ReadShellBezTrianTL },
{ "ELEMENT.SHELL.BSPLINE.SURFACE",       cElement :: ReadShellBsp },
{ "ELEMENT.SHELL.BSPLINE.SURFACETL",     cElement :: ReadShellBspTL },
{ "ELEMENT.SHELL.Q8",                    cElement :: ReadShellQ8 },
{ "ELEMENT.SHELL.Q8TL",                  cElement :: ReadShellQ8TL },
{ "ELEMENT.SHELL.Q9",                    cElement :: ReadShellQ9 },
{ "ELEMENT.SHELL.REFERENCE.VECTOR",      cShapeShell :: ReadReferenceVector},
{ "ELEMENT.SHELL.T6",                    cElement :: ReadShellT6 },
{ "ELEMENT.SHELL.UPDATE.OPTION",         cElmDegShellTL :: ReadUpdateOption},	
{ "ELEMENT.TET10",                       cElement :: ReadTet10 },
{ "ELEMENT.TET10TL",                     cElement :: ReadTet10TL },
{ "ELEMENT.TET4",                        cElement :: ReadTet4 },
{ "ELEMENT.TET4TL",                      cElement :: ReadTet4TL },
{ "ELEMENT.THERMAL2D.BEZIER.SURFACE",    cElement :: ReadPlThermBezSurf },
{ "ELEMENT.THERMAL2D.BEZIER.TRIANGLE",   cElement :: ReadPlThermBezTri },
{ "ELEMENT.THERMAL2D.Q4",                cElement :: ReadPlThermQ4 },
{ "ELEMENT.THERMAL2D.Q8",                cElement :: ReadPlThermQ8 },
{ "ELEMENT.THERMAL2D.T3",                cElement :: ReadPlThermT3 },
{ "ELEMENT.THERMAL2D.T6",                cElement :: ReadPlThermT6 },
{ "ELEMENT.THICKPLT.Q4",                 cElement :: ReadThickPltQ4 },
{ "ELEMENT.THICKPLT.Q8",                 cElement :: ReadThickPltQ8 },
{ "ELEMENT.THICKPLT.T3",                 cElement :: ReadThickPltT3 },
{ "ELEMENT.THICKPLT.T6",                 cElement :: ReadThickPltT6 },
{ "ELEMENT.TRUSS3D",                     cElement :: ReadTruss3D },
{ "ELEMENT.TRUSS3DCR",                   cElement :: ReadTruss3DCR },
{ "ELEMENT.TRUSS3DTL",                   cElement :: ReadTruss3DTL },
{ "GEOSTATIC.STRESSES",                  cElement :: ReadGeoStat },
{ "GEOSTATIC.VERTICAL.AXIS",             cElement :: ReadGeoVertAxis },
{ "HEADER.ANALYSIS.ALGORITHM",           cControl :: ReadAlg },
{ "HEADER.ANALYSIS.ARC-LENGTH",          cCtrlArcLength :: ReadData },
{ "HEADER.ANALYSIS.CRITICAL.POINT",      cCtrlPath :: ReadStabData },
{ "HEADER.ANALYSIS.DISPLACEMENT.CONTROL",cCtrlDispl :: ReadData },
{ "HEADER.ANALYSIS.EIGEN.ALGORITHM",     cControl :: ReadEigAlgType },
{ "HEADER.ANALYSIS.EIGEN.NUMBER.MODES",  cControl :: ReadNumModes },
{ "HEADER.ANALYSIS.EIGEN.SHIFT.FACTOR",  cGenEigShfPower :: ReadShiftFactor },
{ "HEADER.ANALYSIS.GEN-ALPHA.PARAMETERS",cCtrlGenAlpha :: ReadGenAlphaParam },
{ "HEADER.ANALYSIS.HHT.PARAMETER",       cCtrlHHT :: ReadHHTParam },
{ "HEADER.ANALYSIS.LINEAR.SOLVER",       cControl :: ReadLinearSolver },
{ "HEADER.ANALYSIS.MATRIX.TYPE",         cControl :: ReadMatrixType },
{ "HEADER.ANALYSIS.MAXIMUM.ITERATIONS",  cControl :: ReadMaxIter },
{ "HEADER.ANALYSIS.MAXIMUM.STEPS",       cControl :: ReadMaxStep },
{ "HEADER.ANALYSIS.NEWMARK.PARAMETERS",  cCtrlNewmark :: ReadNewmarkParam },
{ "HEADER.ANALYSIS.PETSC.OPTIONS",       cControl :: ReadPetscOptions },
{ "HEADER.ANALYSIS.PETSC.OPTIONS.FILE",  cControl :: ReadPetscOptionsFile },
{ "HEADER.ANALYSIS.PRINT.STEPS",         cControl :: ReadPrintStep },
{ "HEADER.ANALYSIS.PROPORTIONAL.LOADING",cControl :: ReadPropLoad },
{ "HEADER.ANALYSIS.RAYLEIGH.DAMPING",    cCtrlNewmark :: ReadRayleighDamp },
{ "HEADER.ANALYSIS.STEP.FACTOR",         cControl :: ReadFactor },
{ "HEADER.ANALYSIS.TARGET.ITERATIONS",   cControl :: ReadTargetIter },
{ "HEADER.ANALYSIS.TIME.STEP",           cControl :: ReadTimeStep },
{ "HEADER.ANALYSIS.TOLERANCE",           cControl :: ReadTol },
{ "HEADER.EQUIVALENT.MESH.DIVISION",     cShapeIGA :: ReadEqvMshDiv },
{ "HEADER.OUTPUT.PATCH.STRESS.SMOOTH",   cControl :: ReadPatStrSmt },
{ "HEADER.OUTPUT.PRECISION",             cControl :: ReadOutPrec },
{ "HEADER.OUTPUT.STRESS.EXTRAPOLATION",  cControl :: ReadStrExtFlag },
{ "HEADER.PRINT.AVERAGE.NODAL.STRESS",   cControl :: ReadPrintAverageStress },
{ "HEADER.PRINT.EQUIVALENT.MESH",        cControl :: ReadPrintEqvMesh },
{ "HEADER.PRINT.SUPPORT.REACTIONS",      cControl :: ReadPrintReactions },
{ "INITIAL.ELEMENT.STRESSES",            cElement :: ReadInitialStress },
{ "INTEGRATION.ORDER",                   cElement :: ReadIntOrd },
{ "INTEGRATION.TYPE",                    cElement :: ReadIntType },
{ "LAMINATE.LOCAL.STRESS.OUTPUT",        cControl :: ReadLamLocStr },
{ "LOAD.CASE.AREA.IGA.UNIFORM",          cLoad :: ReadFaceUnifIgaLoads },
{ "LOAD.CASE.AREA.UNIFORM",              cLoad :: ReadFaceUnifLoads },
{ "LOAD.CASE.BODY.FORCES",               cLoad :: ReadBodyForces },
{ "LOAD.CASE.BODY.HEAT.SOURCE",          cHeat :: ReadBodySources },
{ "LOAD.CASE.BODY.HEAT.SOURCE.FIELD",    cHeat :: ReadFieldBodySources },
{ "LOAD.CASE.BODY.HEAT.SOURCE.GENERAL",  cHeat :: ReadBodyGenSources },
{ "LOAD.CASE.DONNELL.UNIFORM",           cLoad :: ReadDonnellUnifLoads },
{ "LOAD.CASE.ELEMENT.TEMPERATURE",        cElement :: ReadElemTemp },
{ "LOAD.CASE.FRAME3D.UNIFORM",           cLoad :: ReadFrame3DUnifLoads},
{ "LOAD.CASE.GRID.UNIFORM",              cLoad :: ReadGridUnifLoads},
{ "LOAD.CASE.LINE.FORCE.FIELD",          cLoad :: ReadLineFieldLoads },
{ "LOAD.CASE.LINE.FORCE.GENERAL",        cLoad :: ReadLineGenLoads },
{ "LOAD.CASE.LINE.FORCE.IGA.GENERAL",    cLoad :: ReadLineGenIgaLoads },
{ "LOAD.CASE.LINE.FORCE.IGA.UNIFORM",    cLoad :: ReadLineUnifIgaLoads },
{ "LOAD.CASE.LINE.FORCE.UNIFORM",        cLoad :: ReadLineUnifLoads },
{ "LOAD.CASE.LINE.HEAT.FLUX.FIELD",      cHeat :: ReadLineFieldFluxes },
{ "LOAD.CASE.LINE.HEAT.FLUX.GENERAL",    cHeat :: ReadLineGenFluxes },
{ "LOAD.CASE.LINE.HEAT.FLUX.UNIFORM",    cHeat :: ReadLineUnifFluxes },
{ "LOAD.CASE.LINE.MOMENT.IGA.UNIFORM",   cLoad :: ReadLineUnifIgaMoments },
{ "LOAD.CASE.LINE.MOMENT.UNIFORM",       cLoad :: ReadLineUnifMoments },
{ "LOAD.CASE.NODAL.FORCE",               cLoad :: ReadNodalLoads },
{ "LOAD.CASE.PARAMETRIC.CONCENTRATE",    cLoad :: ReadParConcForces },
{ "LOAD.CASE.PATCH.CONCENTRATE",         cLoad :: ReadPatConcForces },
{ "LOAD.CASE.PLFRAME.UNIFORM",           cLoad :: ReadPlFrameUnifLoads },
{ "LOAD.CASE.PRESCRIBED.DOF",            cPrescDof :: ReadPrescDofs },
{ "LOAD.CASE.PRESCRIBED.TEMPERATURE",    cPrescDof :: ReadPrescTemps },
{ "LOAD.CASE.SHELL.UNIFORM",             cLoad :: ReadShellUnifLoads },
{ "MATERIAL",                            cMaterial :: ReadNumMat },
{ "MATERIAL.CONCRETE.EUROCEB",           cMaterial :: ReadConcEuroCEB },
{ "MATERIAL.CONCRETE.LEE.FENVES",        cMaterial :: ReadLeeFenves },
{ "MATERIAL.CONCRETE.MAZARS",            cMaterial :: ReadDamageMazars },
{ "MATERIAL.CONCRETE.MU-MODEL",          cMaterial :: ReadDamageMuModel },
{ "MATERIAL.CONCRETE.NBRCEB",            cMaterial :: ReadConcNBRCEB },
{ "MATERIAL.CZM.MODEI.BILINEAR",         cMaterial :: ReadCZMModeIBilinear },
{ "MATERIAL.CZM.MODEI.EXPONENTIAL",      cMaterial :: ReadCZMModeIExponential },
{ "MATERIAL.DENSITY",                    cMaterial :: ReadDensity },
{ "MATERIAL.ELASTOPLASTIC",              cMaterial :: ReadPlastLinHard},
{ "MATERIAL.ELASTOPLASTIC.ISOEXPHARD",   cMaterial :: ReadPlastIsoExpHard},
{ "MATERIAL.ELASTOPLASTIC.ISOLINHARD",   cMaterial :: ReadPlastIsoLinHard},
//{ "MATERIAL.ELASTOPLASTIC.ISONONLINHARD",cMaterial :: ReadPlastIsoNonlinHard},
{ "MATERIAL.FGM",                        cMaterial :: ReadFGM },
{ "MATERIAL.FGM.TTO",                    cMaterial :: ReadFGMTTO },
{ "MATERIAL.ISOTROPIC",                  cMaterial :: ReadElastIso },
{ "MATERIAL.NONLINEAR.ELASTIC",          cMaterial :: ReadElastNonlin },
{ "MATERIAL.ORTHOTROPIC",                cMaterial :: ReadElastOrtho },
{ "MATERIAL.PROGFAIL.ENGELSTAD",         cMaterial :: ReadProgFailEngelstad },
{ "MATERIAL.PROGFAIL.HASHIN",            cMaterial :: ReadProgFailHashin },
{ "MATERIAL.PROGFAIL.TSAI",              cMaterial :: ReadProgFailTsai },
{ "MATERIAL.THERMAL.EXPANSION",          cMaterial :: ReadThermExpan },
{ "MATERIAL.THERMAL.ISOTROPIC",          cMaterial :: ReadThermalIso },
{ "MATERIAL.THERMAL.ORTHOTROPIC",        cMaterial :: ReadThermalOrtho },
{ "MATERIAL.VISCOELASTIC",               cMaterial :: ReadViscoElastic },
{ "MATHEMATICAL.EXPRESSION",             Utl       :: ReadExp },
{ "MICRO.MODEL",                         cControl :: ReadMicroModel },
{ "MICRO.MODEL.MACRO.STRAINS",           cControl :: ReadMacroStrain },
{ "MICRO.MODEL.NODE.PAIRS",              cControl :: ReadMicroPair },
{ "MICRO.MODEL.REGULARIZATION.FACTOR",   cControl :: ReadMicroRegFac },
{ "MICRO.MODEL.VERTEX.NODES",            cControl :: ReadMicroVtx },
{ "NODAL.MASS",                          cNode :: ReadNodeMass },
{ "NODE",                                cNode :: ReadNumNode },
{ "NODE.CONSTRAINT",                     cNode :: ReadNodeConstraint },
{ "NODE.COORD",                          cNode :: ReadCoord },
{ "NODE.DEFAULT.THICKNESS",              cNode :: ReadDefThk },
{ "NODE.MULTI-POINT.CONSTRAINT",         cNode :: ReadNodeMPC },	
{ "NODE.NORMAL",                         cNode :: ReadNodeNormal },
{ "NODE.SOLVER.ORDER",                   cNode :: ReadSolverOrder },
{ "NODE.SUPPORT",                        cNode :: ReadSupp },
{ "NODE.THICKNESS",                      cNode :: ReadThk },
{ "PATCH",                               ReadPatch },
{ "PATCH.OUTPUT.POINT",                  cShapeIGA :: ReadPatPnt },
{ "SECTION",                             cSection :: ReadNumSec },
{ "SECTION.BAR.CIRCLE",                  cSection :: ReadSecCirc },
{ "SECTION.BAR.GENERAL",                 cSection :: ReadSecGen },
{ "SECTION.BAR.INTEGRATION",             cSection :: ReadSecInteg },
{ "SECTION.BAR.RECTANGLE",               cSection :: ReadSecRect },
{ "SECTION.BAR.TUBE",                    cSection :: ReadSecTube },
{ "SECTION.FGM.2D",                      cSection :: ReadFGM2D },
{ "SECTION.FGM.3D",                      cSection :: ReadFGM3D },
{ "SECTION.FGM.SHELL",                   cSection :: ReadFGMShell },
{ "SECTION.FRAME3D.COUPLED",             cSection :: ReadSecB3DCoupled },
{ "SECTION.GENERAL.SHELL",               cSection :: ReadGeneralShell },
{ "SECTION.HOMOGENEOUS.ISOTROPIC.2D",    cSection :: ReadHomIso2D },
{ "SECTION.HOMOGENEOUS.ISOTROPIC.3D",    cSection :: ReadHomIso3D },
{ "SECTION.HOMOGENEOUS.ORTHOTROPIC.2D",  cSection :: ReadHomOrtho2D },
{ "SECTION.HOMOGENEOUS.ORTHOTROPIC.3D",  cSection :: ReadHomOrtho3D },
{ "SECTION.HSDT.FGM.PLATE",              cSection :: ReadHSDTFGMPlate },
{ "SECTION.HSDT.PLATE",                  cSection :: ReadHSDTPlate },
{ "SECTION.INTERFACE.2D",                cSection :: ReadInterface2D },
{ "SECTION.INTERFACE.3D",                cSection :: ReadInterface3D },
{ "SECTION.LAMINATED.2D",                cSection :: ReadLaminated2D },
{ "SECTION.LAMINATED.SHELL",             cSection :: ReadLaminatedShell },
{ "SECTION.LAMINATED.SOLID",             cSection :: ReadLaminatedSolid },
{ "SECTION.RC.RECTANGLE",                cSection :: ReadSecRCRect },
{ "SOLUTION.DISPLACEMENT.FIELD",         cField :: ReadDisplField     },
{ "SOLUTION.STRAIN.FIELD",               cField :: ReadStrainField    },
{ "SOLUTION.STRESS.FIELD",               cField :: ReadStressField    },
{ "SOLUTION.TEMPERATURE.FIELD",          cField :: ReadTempField      },
{ "SPRING.CONNECTION",                   cSpringConn :: ReadSpringConn },
{ "SPRING.PROPERTY",                     cSpringProp :: ReadSpringProp },
{ "SPRING.SUPPORT",                      cSpringSupp :: ReadSpringSupp },
{ "TIME.FUNCTION.HALFSINE",              cTimeFunc :: ReadTimeFuncHalfSine },
{ "TIME.FUNCTION.HARMONIC",              cTimeFunc :: ReadTimeFuncHarmonic },
{ "TIME.FUNCTION.TABLE",                 cTimeFunc :: ReadTimeFuncTable }
};

static int KeySiz = sizeof(TabKey[0]);
static int NumKey = sizeof(TabKey)/KeySiz;

// -------------------------------------------------------------------------
// Local functions:
//

// =============================== CompFunc ================================

static int CompFunc(const void *a, const void *b)
{
  const char *sa = (const char *)a;
  const char *sb = ((sInpData *)b)->label;

  return(strcmp(sa, sb));
}

// -------------------------------------------------------------------------
// Public functions:
//

// ============================== InpReadData ==============================

void InpReadData(void)
{
#ifdef _DEBUG_
  for (int i = 0; i < NumKey-1; i++)
  {
    if (strcmp(TabKey[i].label, TabKey[i+1].label) >= 0)
    {
      cout << "\nError: keywords not sorted in input.cpp!\n";
      cout << TabKey[i].label << "  -  " << TabKey[i+1].label << "\n\n";
      exit(0);
    }
  }
#endif

  // Read the input file.

  char label[BUFSIZ];
  sInpData *p;
  while (1)
  {
    // Read the next label.

    if (!Utl::NextLabel(in, label))
    {
      cout << "Input error or label END not found!\n";
      exit(0);
    }
//    cout << "label = |" << label << "|\n";

    // Check for the end of file.

    if (string(label) == "END") break;

    // Get the input function from the given label.

    p = (sInpData *)bsearch(label, TabKey, NumKey, KeySiz, CompFunc);

#ifdef _DEBUG_
    if (!p && Feedback) 
      cout << "Warning - ignored label:" << label << "\n";
#endif

    // Read the data from the input file.

    if (p && p->func) p->func( );
  }
}

// ============================== InpCopyData ==============================

void InpCopyData(void)
{
  // Scan the first line (ignore the label HEADER).

  string line;
  in.seekg(0, ios::beg);
  getline(in,line);

  // Copy the input file in the output file (before the %END).

  while(!in.eof())
  {
    getline(in,line);
    if (string(line) == "%END") break;
    out << line << "\n";
  }
}

// ============================================== Fim do arquivo ===========
