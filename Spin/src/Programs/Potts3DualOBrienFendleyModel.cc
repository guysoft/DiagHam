#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/Potts3ChainDualOBrienFendleyHamiltonian.h"

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
  OptionManager Manager ("Potts3DualOBrienFendleyModel" , "0.01");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 10);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-q", "initial q sector that has to computed (can be either 0, 1 or 2)", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-q", "number of q value to evaluate (0 for all q sectors)", 0);
  (*SystemGroup) += new  BooleanOption  ('\n', "periodic", "use periodic boundary conditions");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new SingleStringOption  ('\n', "export-hamiltonian", "export the hamiltonian in a column formatted ASCII file");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type Potts3DualOBrienFendleyModel -h" << endl;
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
      sprintf (OutputFileName, "potts3_opendualobrienfendley_n_%d", NbrSpins);
      sprintf (CommentLine, " Potts 3 dual O'Brien-Fendley model with %d sites and open boundary conditions \n#\n# Q", NbrSpins);
    }
  else
    {
      sprintf (OutputFileName, "potts3_periodicdualobrienfendley_n_%d", NbrSpins);
      sprintf (CommentLine, " Potts 3 dual O'Brien-Fendley model with %d sites and periodic boundary conditions\n#\n# Q", NbrSpins);
    }

  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);

  int MaxQValue = 2;
  int InitialQValue = 0;
  if (Manager.GetInteger("initial-q") >= 0)
    {
      InitialQValue = Manager.GetInteger("initial-q") % 3;
    }
  if (Manager.GetInteger("nbr-q") > 0)
    {
      MaxQValue = InitialQValue + Manager.GetInteger("nbr-q") - 1;
      if (MaxQValue >= 3)
	MaxQValue = 2;
    }
  bool FirstRun = true;

  for (; InitialQValue <= MaxQValue; ++InitialQValue)
    {
      Potts3Chain* Chain = new Potts3Chain (NbrSpins, InitialQValue, 1000000);      
      Potts3ChainDualOBrienFendleyHamiltonian* Hamiltonian = 0;
      Hamiltonian = new Potts3ChainDualOBrienFendleyHamiltonian (Chain, NbrSpins, Manager.GetBoolean("periodic"), 
								 ((long) Manager.GetInteger("memory")) << 20);
      char* TmpQString = new char[64];
      sprintf (TmpQString, "%d", InitialQValue);
      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
      sprintf (TmpEigenstateString, "%s_q_%d", OutputFileName, InitialQValue);
      GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpQString, CommentLine, 0.0,  FullOutputFileName,
			       FirstRun, TmpEigenstateString);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      FirstRun = false;
      delete Chain;
      delete[] TmpQString;
      delete Hamiltonian;
    }
  return 0;
}
