#include "Hamiltonian/ShastrySutherlandHamiltonian.h"
#include "Hamiltonian/ShastrySutherlandAnd2DTranslationHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "HilbertSpace/Spin1_2ChainNewAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainMirrorSymmetry.h"
#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1_2ChainFullAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainFullInversionAnd2DTranslation.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"
#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

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
  OptionManager Manager ("ShastrySutherlandModel" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of unit cells along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of unit cells along the y direction", 3);
  (*SystemGroup) += new SingleDoubleOption ('\n', "j-value", "Heisenberg coupling constant between neighboring sites", 1.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "jp-value", "dimer coupling constant", 2.0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);
  (*SystemGroup) += new SingleStringOption ('\n', "selected-points", "provide a two column ascii file that indicates which momentum sectors have to be computed");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-momentum", "disable momentum quantum numbers even if the system is translation invariant");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-inversion", "disable the inversion symmetry quantum number");
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
      cout << "see man page for option syntax or type ShastrySutherlandModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSitesX = Manager.GetInteger("nbr-sitex");
  int NbrSitesY = Manager.GetInteger("nbr-sitey");
  int NbrSpins = 4 * NbrSitesX * NbrSitesY;
  
  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  char* BoundaryName = new char [16];
  sprintf (OutputFileName, "spin_1_2_shastrysutherland_n_%d_x_%d_y_%d", NbrSpins, NbrSitesX, NbrSitesY);
  sprintf (CommentLine, " ising with %d unit cells in the x direction, %d unit cells in the y direction \n#", NbrSitesX, NbrSitesY);

  char* OutputParameterFileName = new char [256];
  sprintf (OutputParameterFileName + strlen(OutputParameterFileName), "_j_%.6f_jp_%.6f", Manager.GetDouble("j-value"), Manager.GetDouble("jp-value"));
  if (Manager.GetBoolean("disable-momentum") == false)
    Lanczos.SetComplexAlgorithms();
  
  char* FullOutputFileName = new char [strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
  sprintf (FullOutputFileName, "%s_%s.dat", OutputFileName, OutputParameterFileName);

  int MaxSzValue = NbrSpins;
  int InitalSzValue = 0;
  if (Manager.GetInteger("initial-sz") > 1)
    {
      InitalSzValue += (Manager.GetInteger("initial-sz") & ~1);
    }
  if (Manager.GetInteger("nbr-sz") > 0)
    {
      MaxSzValue = InitalSzValue + ((Manager.GetInteger("nbr-sz") - 1) * 2);
    }

  if (Manager.GetBoolean("disable-momentum") == false)
    {
      int NbrMomenta = 0;
      int* XMomenta = 0;
      int* YMomenta = 0;
      if (Manager.GetString("selected-points") == 0)
	{
	  int MinXMomentum = Manager.GetInteger("only-kx");
	  int MaxXMomentum = MinXMomentum;
	  if (MinXMomentum < 0)
	    {
	      MinXMomentum = 0;
	      MaxXMomentum = NbrSitesX - 1;
	    }
	  int MinYMomentum = Manager.GetInteger("only-ky");
	  int MaxYMomentum = MinYMomentum;
	  if (MinYMomentum < 0)
	    {
	      MinYMomentum = 0;
	      MaxYMomentum = NbrSitesY - 1;
	    }
	  NbrMomenta = (MaxXMomentum - MinXMomentum + 1) * (MaxYMomentum - MinYMomentum + 1);
	  XMomenta = new int[NbrMomenta];
	  YMomenta = new int[NbrMomenta];
	  int TmpIndex = 0;
	  for (int i = MinXMomentum; i <= MaxXMomentum; ++i)
	    {
	      for (int j = MinYMomentum; j <= MaxYMomentum; ++j)
		{
		  XMomenta[TmpIndex] = i;
		  YMomenta[TmpIndex] = j;
		  ++TmpIndex;
		}
	    }
	}
      else
	{
	  MultiColumnASCIIFile MomentumFile;
	  if (MomentumFile.Parse(Manager.GetString("selected-points")) == false)
	    {
	      MomentumFile.DumpErrors(cout);
	      return -1;
	    }
	  NbrMomenta = MomentumFile.GetNbrLines();
	  XMomenta = MomentumFile.GetAsIntegerArray(0);
	  YMomenta = MomentumFile.GetAsIntegerArray(1);
	}
      bool FirstRun = true;
      for (int TotalSz = InitalSzValue; TotalSz <= MaxSzValue; TotalSz += 2)
	{ 
	  for (int MomentumSector = 0; MomentumSector < NbrMomenta; ++MomentumSector)
	    {
// 	  int TmpInvXMomenta = (NbrSitesX - XMomenta[MomentumSector]) % NbrSitesX;
// 	  int TmpInvYMomenta = (NbrSitesY - YMomenta[MomentumSector]) % NbrSitesY;
// 	  int MaxInversionSector = -1;
// 	  if ((TmpInvXMomenta == XMomenta[MomentumSector]) && (TmpInvYMomenta == YMomenta[MomentumSector]) && (Manager.GetBoolean("disable-inversion") == false))
// 	    {
// 	      MaxInversionSector = 1;
// 	    }
// 	  for (int InversionSector = -1; InversionSector <= MaxInversionSector; InversionSector += 2)
// 	    {
// 	      cout << "-------------------------------------------" << endl;
 	      cout << "sz=" << TotalSz << "  kx=" << XMomenta[MomentumSector] << "  ky=" << YMomenta[MomentumSector] << endl; 
 	      AbstractSpinChain* Chain = 0;
// 	      if (MaxInversionSector == -1)
	      Chain = new Spin1_2ChainNewAnd2DTranslation (NbrSpins, TotalSz, XMomenta[MomentumSector], NbrSitesX, YMomenta[MomentumSector], NbrSitesY);
// 	      else
// 		Chain = new Spin1_2ChainFullInversionAnd2DTranslation (InversionSector, XMomenta[MomentumSector], NbrSitesX, YMomenta[MomentumSector], NbrSitesY);
 	      if (Chain->GetHilbertSpaceDimension() > 0)
 		{
 		  Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
 		  ShastrySutherlandAnd2DTranslationHamiltonian* Hamiltonian = 0;
 		  Hamiltonian = new ShastrySutherlandAnd2DTranslationHamiltonian(Chain, XMomenta[MomentumSector], NbrSitesX, 
										 YMomenta[MomentumSector], NbrSitesY, 
										 Manager.GetDouble("j-value"), Manager.GetDouble("jp-value"));
		  
 		  char* TmpEigenstateString = new char[strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
// 		  if (MaxInversionSector == -1)
 		    sprintf (TmpEigenstateString, "%s_%s_kx_%d_ky_%d_sz_%d", OutputFileName, OutputParameterFileName, XMomenta[MomentumSector], 
			     YMomenta[MomentumSector], TotalSz);
// 		  else
// 		    sprintf (TmpEigenstateString, "%s_%s_kx_%d_ky_%d_i_%d", OutputFileName, OutputParameterFileName, XMomenta[MomentumSector], YMomenta[MomentumSector], InversionSector);
		  
 		  char* TmpString = new char[32];
// 		  if (Manager.GetBoolean("disable-inversion") == true)
// 		    {
		  sprintf (TmpString, "%d %d %d", TotalSz, XMomenta[MomentumSector], YMomenta[MomentumSector]);
// 		    }
// 		  else
// 		    {
// 		      if (MaxInversionSector == -1)
// 			{
// 			  sprintf (TmpString, "%d %d 0", XMomenta[MomentumSector], YMomenta[MomentumSector]);
// 			}
// 		      else
// 			{
// 			  sprintf (TmpString, "%d %d %d", XMomenta[MomentumSector], YMomenta[MomentumSector], InversionSector);
// 			}
// 		    }
 		  GenericComplexMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
 					      FirstRun, TmpEigenstateString);
 		  MainTaskOperation TaskOperation (&Task);
 		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
 		  FirstRun = false;
 		  delete Hamiltonian;
 		  delete[] TmpString;
 		  delete[] TmpEigenstateString;
 		}
 	      delete Chain;
 	    }
	}      
    }
  else
    {
      bool FirstRun = true;
      for (int TotalSz = InitalSzValue; TotalSz <= MaxSzValue; TotalSz += 2)
	{ 
	  AbstractSpinChain* Chain = new Spin1_2ChainNew (NbrSpins, TotalSz, 10000000);
	  if (Chain->GetHilbertSpaceDimension() > 0)
	    {
	      Architecture.GetArchitecture()->SetDimension(Chain->GetHilbertSpaceDimension());	
	      ShastrySutherlandHamiltonian* Hamiltonian = 0;
	      Hamiltonian = new ShastrySutherlandHamiltonian(Chain, 2 * NbrSitesX, 2 * NbrSitesY, Manager.GetDouble("j-value"), Manager.GetDouble("jp-value"));
	      char* TmpEigenstateString = new char[strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
	      sprintf (TmpEigenstateString, "%s_%s_sz_%d", OutputFileName, OutputParameterFileName, TotalSz);
	      char* TmpString = new char[64];
	      sprintf (TmpString, "%d ", TotalSz);
	      GenericRealMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpString, CommentLine, 0.0,  FullOutputFileName,
				       FirstRun, TmpEigenstateString);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      FirstRun = false;
	      delete Hamiltonian;
	      delete[] TmpString;
	      delete[] TmpEigenstateString;
	    }
	  delete Chain;
	}
    }

  delete[] OutputFileName;
  delete[] CommentLine;
  delete[] FullOutputFileName;
  return 0;
}
