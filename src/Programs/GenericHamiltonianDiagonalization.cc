#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/FileBasedHamiltonian.h"
#include "Hamiltonian/FileBasedHermitianHamiltonian.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"
#include "MainTask/GenericComplexMainTask.h"

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

  (*SystemGroup) += new  SingleStringOption ('\0', "hamiltonian", "text file where the hamiltonian matrix elements are stored");
  (*SystemGroup) += new  BooleanOption ('\n', "fortran", "assume indices are 1-based instead of 0-based (i.e. fortran index convention)");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "skip-lines", "skip the first n-tf lines of the input file", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "data-column", "index of the column that contains the matrix elements (or their real part)", 0);
  (*SystemGroup) += new BooleanOption  ('c', "complex", "indicate that the Hamiltonian is complex");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "prefix to use for output file names", "dummy");
  (*OutputGroup) += new SingleStringOption ('\n', "eigenstate-file", "prefix to use for the eigenstate output file names", "dummy");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
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

  AbstractHamiltonian* Hamiltonian = 0;
  if (Manager.GetBoolean("complex") == false)
    {
      Hamiltonian  = new FileBasedHamiltonian(Manager.GetString("hamiltonian"), Manager.GetInteger("data-column"), false, Manager.GetBoolean("fortran"), Manager.GetInteger("skip-lines"));
    }
  else
    {
      Lanczos.SetComplexAlgorithms();
      Hamiltonian  = new FileBasedHermitianHamiltonian(Manager.GetString("hamiltonian"), Manager.GetInteger("data-column"), false, Manager.GetBoolean("fortran"), Manager.GetInteger("skip-lines"));
    }

  Architecture.GetArchitecture()->SetDimension(Hamiltonian->GetHilbertSpaceDimension());	

  char* CommentLine = new char [strlen(Manager.GetString("hamiltonian")) + 256];
  sprintf (CommentLine, "eigenvalues of %s\n #", Manager.GetString("hamiltonian"));

  char* EigenvectorFileName = Manager.GetString("eigenstate-file");
  
  if (Manager.GetBoolean("complex") == false)
    {
      GenericRealMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, " ", CommentLine, 0.0,  Manager.GetString("output-file"), true, EigenvectorFileName);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
    }
  else
    {
      GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, " ", CommentLine, 0.0,  Manager.GetString("output-file"), true, EigenvectorFileName);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
    }

  return 0;
}
