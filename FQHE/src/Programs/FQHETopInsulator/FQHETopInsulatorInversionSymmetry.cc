#include "Vector/ComplexVector.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETopInsulatorInversionSymmetry", "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\0', "input-state", "name of the file corresponding to the state on which inversion symmetry will be applied");
  (*OutputGroup) += new SingleStringOption ('o', "output-state", "use this output file name instead of the one that can be deduced from the input state file name");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETopInsulatorInversionSymmetry -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-state") == 0)
    {
      cout << "error, a input state file should be provided. See man page for option syntax or type FQHETopInsulatorInversionSymmetry -h" << endl;
      return -1;
    }
  if (IsFile(Manager.GetString("input-state")) == false)
    {
      cout << "can't open file " << Manager.GetString("input-state") << endl;
      return -1;
    }

  int TotalKx = 0;
  int TotalKy = 0;
  int NbrParticles = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  bool Statistics = true;
  double Mass = 0.0;
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
							  NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy, Mass, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
      return -1;
    }

  int InvertTotalKx = (NbrSiteX - TotalKx) % NbrSiteX;
  int InvertTotalKy = (NbrSiteY - TotalKy) % NbrSiteY;

  ComplexVector InputState;
  ComplexVector OutputState;
  if (InputState.ReadVector (Manager.GetString("input-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
      return -1;
    }

  if (Statistics == true)
    {
      FermionOnSquareLatticeMomentumSpace* InputSpace = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
      FermionOnSquareLatticeMomentumSpace* OutputSpace = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, InvertTotalKx, InvertTotalKy);
      cout << "inversion symmetry from (N=" << NbrParticles << ",Nx=" << NbrSiteX << ",Ny=" << NbrSiteY << ",kx=" << TotalKx << ",kx=" << TotalKy << ") to (N=" 
	   << NbrParticles << ",Nx=" << NbrSiteX << ",Ny=" << NbrSiteY << ",kx=" << InvertTotalKx << ",kx=" << InvertTotalKy << ")" << endl; 
      OutputState = OutputSpace->InversionSymmetry(InputState, InputSpace);
      delete InputSpace;
      delete OutputSpace;
    }
  else
    {

    }
  char* OutputFileName;
  if (Manager.GetString("output-state") == 0)
    {
      OutputFileName = new char [32 + strlen(Manager.GetString("input-state"))];
      char* PrefixPosition = strstr (Manager.GetString("input-state"), "kx_");
      char* SuffixPosition = strstr (Manager.GetString("input-state"), "ky_");
      SuffixPosition += 3;
      while (((*SuffixPosition) >= '0') && ((*SuffixPosition) <= '9') && ((*SuffixPosition) != '\0'))
	++SuffixPosition;
      (*PrefixPosition) = '\0';
      sprintf (OutputFileName, "%skx_%d_ky_%d%s", Manager.GetString("input-state"), InvertTotalKx, InvertTotalKy, SuffixPosition);
      (*PrefixPosition) = 'k';
    }
  else
    {
      OutputFileName = new char [1 + strlen(Manager.GetString("output-state"))];
      strcpy (OutputFileName, Manager.GetString("output-state"));
    }
  if (OutputState.WriteVector (OutputFileName) == false)
    {
      cout << "can't write vector file " << OutputFileName << endl;
      return -1;
    }
  delete[] OutputFileName;
  return 0;
}
