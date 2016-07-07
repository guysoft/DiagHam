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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nx1", "first coordinate of the first spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ny1", "second coordinate of the first spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nx2", "first coordinate of the second spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ny2", "second coordinate of the second spanning vector of the tilted lattice", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "real-offset", "second coordinate in real space of the second spanning vector of the real space lattice (0 if lattice is untilted)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "cylinder", "use periodic boundary in the y direction only");
  (*SystemGroup) += new BooleanOption  ('\n', "force-negativesz", "compute negative Sz sectors");
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
  
  int nx1 = Manager.GetInteger("nx1");
  int ny1 = Manager.GetInteger("ny1");
  int nx2 = Manager.GetInteger("nx2");
  int ny2 = Manager.GetInteger("ny2");
  
  int OffsetReal = Manager.GetInteger("real-offset");
  bool TiltedFlag = true;
  if ( ((nx1 == 0) && (ny1 == 0)) || ((nx2 == 0) && (ny2 == 0)) )
    TiltedFlag = false;
  else
    {
      if ((nx1*ny2 - nx2*ny1) != NbrSitesX * NbrSitesY)
	{
	  cout << "Boundary conditions define a lattice that has a number of sites different from NbrSiteX * NbrSiteY - should have (nx1*ny2 - nx2*ny1) = NbrSiteX * NbrSiteY " << endl;
	  return 0;
	}
      
      if ((((OffsetReal*ny2 + nx2) % NbrSitesX) != 0 || ((nx1 + OffsetReal*ny1) % NbrSitesX != 0)))
      {
	  cout << "Tilted lattice not properly defined. Should have ((offset*ny2 + nx2) % NbrSiteX) = 0 and ((nx1 + offset*ny1) % NbrSiteX = 0) to verify momentum conservation" << endl;
	  return 0;
      }
	
	
      cout << "Using tilted boundary conditions" << endl;
    }
    
    
  
  
  char* OutputFileName = new char [512];
  if (Manager.GetBoolean("cylinder"))
  {
    if (TiltedFlag)
      cout << "Error: tilted boundary conditions not supported for cylinder geometry" << endl;
    else
      sprintf (OutputFileName, "spin_1_2_triangle_cylinder_pseudospin_x_%d_y_%d_j_%.6f", NbrSitesX, NbrSitesY, JValue);
  }
  else
  {
    if (TiltedFlag)
      sprintf (OutputFileName, "spin_1_2_triangle_pseudospin_x_%d_y_%d_nx1_%d_ny1_%d_nx2_%d_ny2_%d_off_%d_j_%.6f", NbrSitesX, NbrSitesY, nx1, ny1, nx2, ny2, OffsetReal, JValue);
    else
      sprintf (OutputFileName, "spin_1_2_triangle_pseudospin_x_%d_y_%d_j_%.6f", NbrSitesX, NbrSitesY, JValue);
  }
  char* CommentLine = new char [512];
  sprintf (CommentLine, "spin 1/2 system with boundary conditions on the triangle lattice, pseudospin 1/2 and %d sites in the x direction, %d sites in the y direction \n# Sz", NbrSitesX, NbrSitesY);
  
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);
  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
  
  
  Spin1_2ChainWithPseudospin* Space = 0;
  
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
      Space = new Spin1_2ChainWithPseudospin(NbrSpins, InitalSzValue, 1000000);
      cout << "2Sz = " << InitalSzValue << endl; 
      cout << (Space->GetHilbertSpaceDimension()) << endl;
      if (Space->GetHilbertSpaceDimension() > 0)
      {
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	  TwoDimensionalTriangularLatticeWithPseudospinHamiltonian* Hamiltonian = 0;
	  Hamiltonian = new TwoDimensionalTriangularLatticeWithPseudospinHamiltonian(Space, NbrSitesX, NbrSitesY, JValue, (!Manager.GetBoolean("cylinder")), OffsetReal);
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