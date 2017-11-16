#include "Hamiltonian/SpinChainAKLTHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin1ChainWithSzSymmetry.h"
#include "HilbertSpace/Spin1ChainWithInversionSymmetry.h"
#include "HilbertSpace/Spin1ChainWithSzInversionSymmetries.h"

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
  OptionManager Manager ("SpinChainAKLT" , "0.01");
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

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 2);
  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 8);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "use-periodic", "use periodic boundary conditions");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-szsymmetry", "disable the Sz<->-Sz symmetry");
  (*SystemGroup) += new  BooleanOption ('\n', "disable-inversionsymmetry", "disable the inversion symmetry");
  (*SystemGroup) += new  SingleDoubleOption ('\n', "additional-quadratic", "coefficient in front of the additional quadratic term (0 being the pure AKLT hamiltonian)", 0.0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinChainAKLT -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int SpinValue = Manager.GetInteger("spin");
  int NbrSpins = Manager.GetInteger("nbr-spin");

  char* OutputFileName = new char [512];
  char* CommentLine = new char [512];
  char* BoundaryName = new char [16];
  if (Manager.GetBoolean("use-periodic") == false)
    sprintf (BoundaryName, "open");
  else
    sprintf (BoundaryName, "closed");
  if ((SpinValue & 1) == 0)
    {
      if (Manager.GetDouble("additional-quadratic") != 0.0)
	{
	  sprintf (OutputFileName, "spin_%d_%saklt_quadratic_%.6f_n_%d", (SpinValue / 2), BoundaryName, 
		   Manager.GetDouble("additional-quadratic"), NbrSpins);
	}
      else
	{
	  sprintf (OutputFileName, "spin_%d_%saklt_n_%d", (SpinValue / 2), BoundaryName, NbrSpins);
	}
      if (Manager.GetBoolean("disable-szsymmetry") == false)
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz SzSym InvSym ", BoundaryName, (SpinValue / 2), NbrSpins);
	    }
	  else
	    {
	      sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz SzSym ", BoundaryName, (SpinValue / 2), NbrSpins);
	    }
	}
      else
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz InvSym ", BoundaryName, (SpinValue / 2), NbrSpins);
	    }
	  else
	    {
	      sprintf (CommentLine, " %s spin %d chain with %d sites \n# 2Sz ", BoundaryName, (SpinValue / 2), NbrSpins);
	    }
	}
    }
  else
    {
      if (Manager.GetDouble("additional-quadratic") != 0.0)
	{
	  sprintf (OutputFileName, "spin_%d_2_%saklt_quadratic_%.6f_n_%d", SpinValue, BoundaryName, 
		   Manager.GetDouble("additional-quadratic"), NbrSpins);
	}
      else
	{
	  sprintf (OutputFileName, "spin_%d_2_%saklt_n_%d", SpinValue, BoundaryName, NbrSpins);
	}
      if (Manager.GetBoolean("disable-szsymmetry") == false)
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz SzSym InvSym ", BoundaryName, SpinValue, NbrSpins);
	    }
	  else
	    {
	      sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz SzSym ", BoundaryName, SpinValue, NbrSpins);
	    }
	}
      else
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz InvSym ", BoundaryName, SpinValue, NbrSpins);
	    }
	  else
	    {
	      sprintf (CommentLine, " %s spin %d/2 chain with %d sites \n# 2Sz ", BoundaryName, SpinValue, NbrSpins);
	    }
	}
    }
  char* FullOutputFileName = new char [strlen(OutputFileName)+ 16];
  sprintf (FullOutputFileName, "%s.dat", OutputFileName);

  int MaxSzValue = NbrSpins * SpinValue;
  int InitalSzValue = MaxSzValue & 1;
  if (Manager.GetInteger("initial-sz") != 0)
    {
      InitalSzValue += (Manager.GetInteger("initial-sz") & ~1);
    }
  if (Manager.GetInteger("nbr-sz") > 0)
    {
      MaxSzValue = InitalSzValue + ((Manager.GetInteger("nbr-sz") - 1) * 2);
    }
  bool FirstRun = true;

  if ((InitalSzValue == 0) && (Manager.GetBoolean("disable-szsymmetry") == false))
    {
      for (int SzSymmetrySector = -1; SzSymmetrySector <= 1; SzSymmetrySector += 2)
	{
	  if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	    {
	      for (int InversionSymmetrySector = -1; InversionSymmetrySector <= 1; InversionSymmetrySector += 2)
		{
		  AbstractSpinChain* Chain = 0;
		  switch (SpinValue)
		    {
		    case 2 :
		      Chain = new Spin1ChainWithSzInversionSymmetries (NbrSpins, InversionSymmetrySector, SzSymmetrySector, InitalSzValue, 1000000);
		      break;
		    default :
		      {
			if ((SpinValue & 1) == 0)
			  cout << "spin " << (SpinValue / 2) << " are not available" << endl;
			else 
			  cout << "spin " << SpinValue << "/2 are not available" << endl;
			return -1;
		      }
		    }
		  
		  if (Chain->GetHilbertSpaceDimension() > 0)
		    {
		      cout << "Sz = " << InitalSzValue << ", Sz<->-Sz sector=" << SzSymmetrySector << ", inversion sector=" << InversionSymmetrySector << endl; 
		      SpinChainAKLTHamiltonian Hamiltonian (Chain, NbrSpins, 1.0 + 3.0 * Manager.GetDouble("additional-quadratic"), Manager.GetBoolean("use-periodic"));
		      char* TmpSzString = new char[64];
		      sprintf (TmpSzString, "%d %d %d", InitalSzValue, SzSymmetrySector, InversionSymmetrySector);
		      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		      sprintf (TmpEigenstateString, "%s_sz_%d_invsym_%d_szsym_%d", OutputFileName, InitalSzValue, InversionSymmetrySector, SzSymmetrySector);
		      GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
					       FirstRun, TmpEigenstateString);
		      MainTaskOperation TaskOperation (&Task);
		      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		      FirstRun = false;
		      delete[] TmpSzString;
		    }
		  delete Chain;
		}
	    }
	  else
	    {
	      AbstractSpinChain* Chain = 0;
	      switch (SpinValue)
		{
		case 2 :
		  Chain = new Spin1ChainWithSzSymmetry (NbrSpins, SzSymmetrySector, InitalSzValue, 1000000);
		  break;
		default :
		  {
		    if ((SpinValue & 1) == 0)
		  cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	      
	      if (Chain->GetHilbertSpaceDimension() > 0)
		{
		  cout << "Sz = " << InitalSzValue << ", Sz<->-Sz sector=" << SzSymmetrySector << endl; 
		  SpinChainAKLTHamiltonian Hamiltonian (Chain, NbrSpins, 1.0 + 3.0 * Manager.GetDouble("additional-quadratic"), Manager.GetBoolean("use-periodic"));
		  char* TmpSzString = new char[64];
		  sprintf (TmpSzString, "%d %d", InitalSzValue, SzSymmetrySector);
		  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		  sprintf (TmpEigenstateString, "%s_sz_%d_szsym_%d", OutputFileName, InitalSzValue, SzSymmetrySector);
		  GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
					   FirstRun, TmpEigenstateString);
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  FirstRun = false;
		  delete[] TmpSzString;
		}
	      delete Chain;
	    }
	}
      InitalSzValue += 2;
    }
  for (; InitalSzValue <= MaxSzValue; InitalSzValue +=2)
    {
      if (Manager.GetBoolean("disable-inversionsymmetry") == false)
	{
	  for (int InversionSymmetrySector = -1; InversionSymmetrySector <= 1; InversionSymmetrySector += 2)
	    {
	      AbstractSpinChain* Chain = 0;
	      switch (SpinValue)
		{
		case 2 :
		  Chain = new Spin1ChainWithInversionSymmetry (NbrSpins, InversionSymmetrySector, InitalSzValue, 1000000);
	      break;
		default :
		  {
		    if ((SpinValue & 1) == 0)
		      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	      
	      if (Chain->GetHilbertSpaceDimension() > 0)
		{
		  cout << "Sz = " << InitalSzValue << ", inversion sector=" << InversionSymmetrySector << endl; 
		  SpinChainAKLTHamiltonian Hamiltonian (Chain, NbrSpins, 1.0 + 3.0 * Manager.GetDouble("additional-quadratic"), Manager.GetBoolean("use-periodic"));
		  char* TmpSzString = new char[64];
		  if (Manager.GetBoolean("disable-szsymmetry") == false)
		    {
		      sprintf (TmpSzString, "%d 0 %d", InitalSzValue, InversionSymmetrySector);
		    }
		  else
		    {
		      sprintf (TmpSzString, "%d %d", InitalSzValue, InversionSymmetrySector);
		    }
		  char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
		  sprintf (TmpEigenstateString, "%s_sz_%d_invsym_%d", OutputFileName, InitalSzValue, InversionSymmetrySector);
		  GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
					   FirstRun, TmpEigenstateString);
		  MainTaskOperation TaskOperation (&Task);
		  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		  FirstRun = false;
		  delete[] TmpSzString;
		}
	      delete Chain;
	    }
	}
      else
	{
	  AbstractSpinChain* Chain = 0;
	  switch (SpinValue)
	    {
	    case 1 :
	      Chain = new Spin1_2Chain (NbrSpins, InitalSzValue, 1000000);
	      break;
	    case 2 :
	      Chain = new Spin1Chain (NbrSpins, InitalSzValue, 1000000);
	      break;
	    default :
	      {
		if ((SpinValue & 1) == 0)
		  cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		else 
		  cout << "spin " << SpinValue << "/2 are not available" << endl;
		return -1;
	      }
	    }
	  
	  if (Chain->GetHilbertSpaceDimension() > 0)
	    {
	      cout << "Sz = " << InitalSzValue << endl; 
	      SpinChainAKLTHamiltonian Hamiltonian (Chain, NbrSpins, 1.0 + 3.0 * Manager.GetDouble("additional-quadratic"), Manager.GetBoolean("use-periodic"));
	      char* TmpSzString = new char[64];
	      if (Manager.GetBoolean("disable-szsymmetry") == false)
		{
		  sprintf (TmpSzString, "%d 0", InitalSzValue);
		}
	      else
		{
		  sprintf (TmpSzString, "%d", InitalSzValue);
		}
	      char* TmpEigenstateString = new char[strlen(OutputFileName) + 64];
	      sprintf (TmpEigenstateString, "%s_sz_%d", OutputFileName, InitalSzValue);
	      GenericRealMainTask Task(&Manager, Chain, &Lanczos, &Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
				       FirstRun, TmpEigenstateString);
	      MainTaskOperation TaskOperation (&Task);
	      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	      FirstRun = false;
	      delete[] TmpSzString;
	    }
	  delete Chain;
	}
    }
  return 0;
}
