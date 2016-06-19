#include "Hamiltonian/TwoDimensionalKagomeLatticeHamiltonian.h"

#include "HilbertSpace/AbstractSpinChain.h"
#include "HilbertSpace/Spin1_2ChainNew.h"

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
  OptionManager Manager ("SpinKagomeLattice" , "0.01");
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
  (*SystemGroup) += new BooleanOption  ('\n', "force-negativesz", "compute negative Sz sectors");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 0);
  (*SystemGroup) += new  SingleDoubleOption ('j', "j-value", "coupling constant value for nearest neighbors", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('a', "anisotropy", "anisotropy between up and down triangles", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "easy-plane", "easy plane anisotropy", 1.0);
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
  int NbrSpins = NbrSitesX * NbrSitesY * 3;
  double JValue =  Manager.GetDouble("j-value");
  double JDownValue = JValue * Manager.GetDouble("anisotropy");
  double JEasyPlane = JValue * Manager.GetDouble("easy-plane");
  double JDownEasyPlane = JEasyPlane * Manager.GetDouble("anisotropy");
  
  if (Manager.GetDouble("easy-plane") != 1.0)
    cout << "Warning: easy-plane anisotropy is not tested in this code" << endl;
  
  char* ParametersName = new char[256];
  if (Manager.GetDouble("anisotropy") == 1.0)
  {
    if (Manager.GetDouble("easy-plane") == 1.0)
      sprintf(ParametersName, "heisenberg");
    else
      sprintf(ParametersName, "jx_%.6f_jy_%.6f_jz_%.6f", JEasyPlane, JEasyPlane, JValue);
  }
  else
  {
    if (Manager.GetDouble("easy-plane") == 1.0)
      sprintf(ParametersName, "anisotropy_jup_%.6f_jdown_%.6f", JValue, JDownValue);
    else
      sprintf(ParametersName, "jxup_%.6f_jyup_%.6f_jzup_%.6f_jxdown_%.6f_jydown_%.6f_jzdown_%.6f", JEasyPlane, JEasyPlane, JValue, JDownEasyPlane, JDownEasyPlane, JDownValue);
  }
  
  char* OutputFileName = new char [512];
  if (Manager.GetBoolean("cylinder"))
    sprintf (OutputFileName, "spin_1_2_kagome_cylinder_x_%d_y_%d_%s", NbrSitesX, NbrSitesY, ParametersName);
  else
    sprintf (OutputFileName, "spin_1_2_kagome_x_%d_y_%d_%s", NbrSitesX, NbrSitesY, ParametersName);
  char* CommentLine = new char [512];
  sprintf (CommentLine, "spin 1/2 system with boundary conditions on the kagome lattice and %d sites in the x direction, %d sites in the y direction \n# Sz", NbrSitesX, NbrSitesY);
  
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);
  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
  
  
  AbstractSpinChain* Space = 0;
  
  int MaxSzValue = NbrSpins * SpinValue;
  int InitalSzValue = NbrSpins & 1;
  if (Manager.GetBoolean("force-negativesz"))
    InitalSzValue = -MaxSzValue;
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
      Space = new Spin1_2ChainNew (NbrSpins, InitalSzValue, 1000000);
      cout << "2Sz = " << InitalSzValue << endl; 
      cout << (Space->GetHilbertSpaceDimension()) << endl;
      if (Space->GetHilbertSpaceDimension() > 0)
      {
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	  TwoDimensionalKagomeLatticeHamiltonian* Hamiltonian = 0;
	  Hamiltonian = new TwoDimensionalKagomeLatticeHamiltonian(Space, NbrSitesX, NbrSitesY, JValue, JDownValue, JEasyPlane, JDownEasyPlane, (!Manager.GetBoolean("cylinder")));
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

  delete[] ParametersName;
  delete[] FullOutputFileName;
  delete[] OutputFileName;
  delete[] TmpEigenstateString;
  delete[] CommentLine;
  return 0;
}