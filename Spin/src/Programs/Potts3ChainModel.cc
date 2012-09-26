#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/Potts3ChainHamiltonian.h"

#include "HilbertSpace/Potts3Chain.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"


#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("Potts3ChainModel" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 10);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 0);
  (*SystemGroup) += new  SingleDoubleOption ('u', "potential", "potnetial term coupling to neighboring sites", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('f', "flip", "on-site flip term", 0.0);
  (*SystemGroup) += new  BooleanOption  ('\n', "periodic", "use periodic boundary conditions");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type Potts3ChainModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSpins = Manager.GetInteger("nbr-spin");

  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  if (Manager.GetBoolean("periodic") == false)
    {
      sprintf (OutputFileName, "potts3_openchain_n_%d", NbrSpins);
      sprintf (CommentLine, " open potts 3 chain with %d sites \n# Sz", NbrSpins);
    }
  else
    {
      sprintf (OutputFileName, "potts3_closechain_n_%d", NbrSpins);
      sprintf (CommentLine, " close potts 3 chain with %d sites \n# Sz", NbrSpins);
    }
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);
  double* PotentialTerms;
  if (Manager.GetBoolean("periodic") == false)
    {
      PotentialTerms = new double [NbrSpins - 1];
      PotentialTerms[0] = Manager.GetDouble("potential");
      for (int i = 1; i < (NbrSpins - 1); ++i)
	PotentialTerms[i] = PotentialTerms[0];
    }
  else
    {
      PotentialTerms = new double [NbrSpins];
      PotentialTerms[0] = Manager.GetDouble("potential");
      for (int i = 1; i < NbrSpins; ++i)
	PotentialTerms[i] = PotentialTerms[0];
    }
  double* FlipTerms = new double [NbrSpins];
  FlipTerms[0] = Manager.GetDouble("flip");
  for (int i = 1; i < NbrSpins; ++i)
    FlipTerms[i] = FlipTerms[0];

  int MaxSzValue = 2;
  int InitalSzValue = 0;
  if (Manager.GetInteger("initial-sz") >= 0)
    {
      InitalSzValue = Manager.GetInteger("initial-sz") % 3;
    }
  if (Manager.GetInteger("nbr-sz") > 0)
    {
      MaxSzValue = InitalSzValue + Manager.GetInteger("nbr-sz");
      if (MaxSzValue >= 3)
	MaxSzValue = 2;
    }
  bool FirstRun = true;
  for (; InitalSzValue <= MaxSzValue; ++InitalSzValue)
    {
      AbstractSpinChain* Chain = new Potts3Chain (NbrSpins, InitalSzValue, 1000000);
      
//       cout << "--------------------------------------" << endl;
//       cout << "Sz = " << InitalSzValue << endl;
//       for (int i = 0; i < Chain->GetHilbertSpaceDimension(); ++i)
// 	{
// 	  cout << i << " : ";
// 	  Chain->PrintState(cout, i) << endl;
// 	}

      Potts3ChainHamiltonian Hamiltonian (Chain, NbrSpins, PotentialTerms, FlipTerms, Manager.GetBoolean("periodic"));
      char* TmpSzString = new char[64];
      sprintf (TmpSzString, "%d", InitalSzValue);
      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
      sprintf (TmpEigenstateString, "%s_sz_%d", OutputFileName, InitalSzValue);
      GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
 			       FirstRun, TmpEigenstateString);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      FirstRun = false;
      delete Chain;
      delete[] TmpSzString;
    }
  return 0;
}
