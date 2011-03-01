#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"

#include "FunctionBasis/ParticleOnDiskFunctionBasis.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHEChernInsulatorCorrelation" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PlotOptionGroup = new OptionGroup ("plot options");  
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += PlotOptionGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the vector file describing the state whose density has to be plotted");
  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density insted of density-density correlation", false);

  (*PlotOptionGroup) += new SingleStringOption ('\n', "output", "output file name (default output name replace the .vec extension of the input file with .rho.dat)", 0);
  (*PlotOptionGroup) += new SingleIntegerOption ('\n', "nbr-samplesx", "number of samples along the x direction", 100, true, 10);
  (*PlotOptionGroup) += new SingleIntegerOption ('\n', "nbr-samplesy", "number of samples along the y direction", 100, true, 10);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEChernInsulatorCorrelation -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("state") == 0)
    {
      cout << "FQHEChernInsulatorCorrelation requires an input state" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state")) == false)
    {
      cout << "can't find vector file " << Manager.GetString("state") << endl;
      return -1;      
    }

  int NbrParticles = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  int MomentumX = 0;
  int MomentumY = 0;
  int NbrSamplesX = Manager.GetInteger("nbr-samplesx");
  int NbrSamplesY = Manager.GetInteger("nbr-samplesy");
  bool DensityFlag = Manager.GetBoolean("density");
  bool Statistics = true;

  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(Manager.GetString("state"),
							  NbrParticles, NbrSiteX, NbrSiteY, MomentumX, MomentumY, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("state") << endl;
      return -1;
    }
  
  FermionOnSquareLatticeMomentumSpace* Space = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, MomentumX, MomentumY);
  ComplexVector ComplexState;
  if (ComplexState.ReadVector (Manager.GetString("eigenstate")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
      return -1;      
    }
  if (DensityFlag == false)
    {
//       ParticleOnSphereDensityDensityOperator Operator (Space, i);
//       PrecalculatedValues[i] = Operator.MatrixElement(ComplexState, ComplexState);
    }
  else
    {
//       ParticleOnSphereDensityOperator Operator (Space, i);
//       PrecalculatedValues[i] = Operator.MatrixElement(ComplexState, ComplexState);
    }  
  delete Space;
  ofstream File;
  File.precision(14);
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = 0;
      if (DensityFlag == false)
	{
	  TmpFileName = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "rhorho");
	}
      else
	{
	  TmpFileName = ReplaceExtensionToFileName(Manager.GetString("state"), "vec", "rho");
	}
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << Manager.GetString("state") << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  File.close();
  return 0;
}
