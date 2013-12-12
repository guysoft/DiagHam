#include "Hamiltonian/SpinChainZ2InteractingHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainFull.h"

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
  OptionManager Manager ("Z2InteractingChain" , "0.01");
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

  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 8);
  (*SystemGroup) += new  SingleDoubleOption('j', "j-value", "", 1.0);
  (*SystemGroup) += new  SingleDoubleOption('f', "f-value", "", 1.0);
  (*SystemGroup) += new  SingleDoubleOption('v', "v-value", "", 1.0);
  (*SystemGroup) += new  SingleIntegerOption ('b', "boundary-conditions", "boundary conditions (0 for open, 1 for periodic, -1 for antiperiodic)", 0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type Z2InteractingChain -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSpins = Manager.GetInteger("nbr-spin");
  double JValue = Manager.GetDouble("j-value");
  double FValue = Manager.GetDouble("f-value");
  double VValue = Manager.GetDouble("v-value");
  int BValue = Manager.GetInteger("boundary-conditions");
  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  sprintf (OutputFileName, "z2interactingchain_n_%d_j_%.6f_f_%.6f_v_%.6f_b_%d", NbrSpins, JValue, FValue, VValue, BValue);
  sprintf (CommentLine, " Z2 interacting chain with %d sites, J=%.6f, F= %.6f, V=%.6f and boundary conditions B=%d\n# ", NbrSpins, JValue, FValue, VValue, BValue);

  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);

  bool FirstRun = true;
  Spin1_2Chain* Chain = new Spin1_2ChainFull (NbrSpins, 1000000);

  if (Chain->GetHilbertSpaceDimension() > 0)
    {
      SpinChainZ2InteractingHamiltonian Hamiltonian (Chain, NbrSpins, JValue, FValue, VValue, (double) BValue);
      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
      sprintf (TmpEigenstateString, "%s", OutputFileName);
      char TmpEntry = '\0';
      GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, &TmpEntry, CommentLine, 0.0,  FullOutputFileName,
			       FirstRun, TmpEigenstateString);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      FirstRun = false;
      delete[] TmpEigenstateString;
    }
  delete Chain;
  return 0;
}
