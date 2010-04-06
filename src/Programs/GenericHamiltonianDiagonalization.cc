#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/FileBasedHamiltonian.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"

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
  OptionManager Manager ("GenericHamiltonianDiagonalization" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleStringOption ('\0', "hamiltonian", "text file where the hamiltonian matrix elements are stored");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericHamiltonianDiagonalization -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("hamiltonian") == 0)
    {
      cout << "no hamiltonian provided" << endl; 
      return -1;
    }
  if (IsFile(Manager.GetString("hamiltonian")) == false)
    {
      cout << "can't open hamiltonian " << Manager.GetString("hamiltonian") << endl; 
      return -1;
    }

  FileBasedHamiltonian Hamiltonian (Manager.GetString("hamiltonian"), 0, false, true);

  char* CommentLine = new char [strlen(Manager.GetString("hamiltonian")) + 256];
  sprintf (CommentLine, "eigenvalues of %s", Manager.GetString("hamiltonian"));

  GenericRealMainTask Task(&Manager, Hamiltonian.GetHilbertSpace(), &Lanczos, &Hamiltonian, " ", CommentLine, 0.0, "toto");
  MainTaskOperation TaskOperation (&Task);
  TaskOperation.ApplyOperation(Architecture.GetArchitecture());

  return 0;
}
