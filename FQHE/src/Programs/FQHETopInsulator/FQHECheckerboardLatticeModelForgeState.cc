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

#include "FunctionBasis/ParticleOnChernInsulatorSingleBandFunctionBasis.h"
#include "FunctionBasis/ParticleOnCheckerboardLatticeFunctionBasis.h"

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
  OptionManager Manager ("FQHECheckerboardLatticeModelForgeState" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the vector file describing the state whose density has to be plotted");
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0 - 0.5 * M_SQRT2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tpp", "second next nearest neighbor hoping amplitude", 0.5 * (M_SQRT2 - 1.0));
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECheckerboardLatticeModelForgeState -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("state") == 0)
    {
      cout << "FQHECheckerboardLatticeModelForgeState requires an input state" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("state")) == false)
    {
      cout << "can't find vector file " << Manager.GetString("state") << endl;
      return -1;      
    }

  int NbrParticles = 2;
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey"); 

  bool Statistics = true;
  
  // ParticleOnCheckerboardLatticeFunctionBasis Basis ();
  for (int i = 0; i < NbrSitesX; ++i)
    {
      for (int j = 0; j < NbrSitesY; ++j)
	{
	  FermionOnSquareLatticeMomentumSpace* Space = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSitesX, NbrSitesY, i, j);
	  ComplexVector State (Space->GetHilbertSpaceDimension(), true);
	  for (int kx1 = 0; kx1 < NbrSitesX; ++kx1)
	    for (int kx2 = 0; kx2 < NbrSitesX; ++kx2)
	      {
		if (((kx1 + kx2) % NbrSitesX) == i)
		  {
		    for (int ky1 = 0; ky1 < NbrSitesY; ++ky1)
		      for (int ky2 = 0; ky2 < NbrSitesY; ++ky2)
			{
			  if (((ky1 + ky2) % NbrSitesY) == j)
			    {
			      
			    }  
			}
		  }
	      }
	  char* EigenstateOutputFile = new char [512];
	  if (Manager.GetDouble("mu-s") == 0.0)
	    sprintf (EigenstateOutputFile, "fermions_singleband_checkerboardlattice_twoparticles_n_%d_x_%d_y_%d_t1_%f_t2_%f_gx_%f_gy_%f_kx_%d_ky_%d.0.vec", NbrParticles, NbrSitesX, NbrSitesY, 
		     Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), i, j);
	  else
	    sprintf (EigenstateOutputFile, "fermions_singleband_checkerboardlattice_twoparticles_n_%d_x_%d_y_%d_t1_%f_t2_%f_gx_%f_gy_%f_mus_%f_kx_%d_ky_%d.0.vec", NbrParticles, NbrSitesX, NbrSitesY, 
		     Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Manager.GetDouble("mu-s"), i, j);
	  State.WriteVector(EigenstateOutputFile);
	  delete[] EigenstateOutputFile;
	  delete Space;
	}
    }

  return 0;
}
