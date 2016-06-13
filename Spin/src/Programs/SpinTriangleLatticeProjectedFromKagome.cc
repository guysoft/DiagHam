#include "Hamiltonian/TwoDimensionalTriangularLatticeWithPseudospinHamiltonian.h"

#include "HilbertSpace/Spin1_2ChainWithPseudospin.h"

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
  OptionManager Manager ("SpinTriangleLatticeProjectedFromKagome" , "0.01");
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

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 1);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new BooleanOption  ('\n', "cylinder", "use periodic boundary in the y direction only");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 0);
  (*SystemGroup) += new  SingleDoubleOption ('j', "j-value", "coupling constant value", 1.0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericPeriodicSpinChain -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int SpinValue = Manager.GetInteger("spin");
  int NbrSitesX = Manager.GetInteger("nbr-sitex");
  int NbrSitesY = Manager.GetInteger("nbr-sitey");
  int NbrSpins = NbrSitesX * NbrSitesY;
  double JValue =  Manager.GetDouble("j-value");
  
  
  char* OutputFileName = new char [512];
  if (Manager.GetBoolean("cylinder"))
    sprintf (OutputFileName, "spin_1_2_triangle_cylinder_pseudospin_x_%d_y_%d_j_%.6f", NbrSitesX, NbrSitesY, JValue);
  else
    sprintf (OutputFileName, "spin_1_2_triangle_pseudospin_x_%d_y_%d_j_%.6f", NbrSitesX, NbrSitesY, JValue);
  char* CommentLine = new char [512];
  sprintf (CommentLine, "spin 1/2 system with boundary conditions on the triangle lattice, pseudospin 1/2 and %d sites in the x direction, %d sites in the y direction \n# Sz", NbrSitesX, NbrSitesY);
  
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);
  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
  
  
  Spin1_2ChainWithPseudospin* Space = 0;
  
  int MaxSzValue = NbrSpins * SpinValue;
  int InitalSzValue = -MaxSzValue;
  if (Manager.GetInteger("initial-sz") > 1)
    {
      InitalSzValue += (Manager.GetInteger("initial-sz") & ~1);
    }
  if (Manager.GetInteger("nbr-sz") > 0)
    {
      MaxSzValue = InitalSzValue + ((Manager.GetInteger("nbr-sz") - 1) * 2);
    }
  bool FirstRun = true;
  for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
    {
      Space = new Spin1_2ChainWithPseudospin(NbrSpins, InitalSzValue, 1000000);
      cout << "2Sz = " << InitalSzValue << endl; 
      cout << (Space->GetHilbertSpaceDimension()) << endl;
      if (Space->GetHilbertSpaceDimension() > 0)
      {
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	  TwoDimensionalTriangularLatticeWithPseudospinHamiltonian* Hamiltonian = 0;
	  Hamiltonian = new TwoDimensionalTriangularLatticeWithPseudospinHamiltonian(Space, NbrSitesX, NbrSitesY, JValue, (!Manager.GetBoolean("cylinder")));
	  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	  sprintf (TmpEigenstateString, "%s_sz_%d", OutputFileName, InitalSzValue);
	  char* TmpSzString = new char[64];
	  sprintf (TmpSzString, "%d", InitalSzValue);
	
	  GenericRealMainTask Task(&Manager, Space, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
				   FirstRun, TmpEigenstateString);
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  FirstRun = false;
	  delete Hamiltonian;
	  delete[] TmpSzString;
	  delete[] TmpEigenstateString;
      }
    }

  delete[] FullOutputFileName;
  delete[] OutputFileName;
  delete[] TmpEigenstateString;
  delete[] CommentLine;
  return 0;
}