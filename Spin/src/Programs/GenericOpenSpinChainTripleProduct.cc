#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Hamiltonian/SpinChainTripleProductHamiltonian.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1Chain.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

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
  OptionManager Manager ("GenericOpenSpinChainTripleProduct" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleIntegerOption ('s', "spin", "twice the spin value", 1);
  (*SystemGroup) += new  SingleIntegerOption ('p', "nbr-spin", "number of spins", 10);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "initial-sz", "twice the initial sz sector that has to computed", 0);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "nbr-sz", "number of sz value to evaluate (0 for all sz sectors)", 0);
  (*SystemGroup) += new  SingleDoubleOption ('j', "j-value", "isotropic coupling constant value", 1.0);
  (*SystemGroup) += new  SingleDoubleOption ('z', "djz-value", "delta compare to the coupling constant value along z", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('c', "chi-value", "coupling constant  in front of the triple product S_i . (S_{i+1} x S_{i+2})", 0.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "hz-value", "amplitude of the Zeeman term along the z axis", 0.0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericOpenSpinChainTripleProduct -h" << endl;
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
  if ((SpinValue & 1) == 0)
    {
      sprintf (OutputFileName, "spin_%d_openchain_n_%d", (SpinValue / 2), NbrSpins);
      sprintf (CommentLine, " open spin %d chain with %d sites \n# 2Sz ", (SpinValue / 2), NbrSpins);
    }
  else
    {
      sprintf (OutputFileName, "spin_%d_2_openchain_n_%d", SpinValue, NbrSpins);
      sprintf (CommentLine, " open spin %d/2 chain with %d sites \n# 2Sz", SpinValue, NbrSpins);
    }
  char* OutputParameterFileName = new char [256];
  if (Manager.GetDouble("djz-value") == 0)
    {      
      if (Manager.GetDouble("hz-value") == 0)
	{
	  sprintf (OutputParameterFileName, "j_%f_chi_%f", Manager.GetDouble("j-value"), Manager.GetDouble("chi-value"));
	}
      else
	{
	  sprintf (OutputParameterFileName, "j_%f_chi_%f_hz_%f", Manager.GetDouble("j-value"), Manager.GetDouble("chi-value"), Manager.GetDouble("hz-value"));
	}
    }
  else
    {
      if (Manager.GetDouble("hz-value") == 0)
	{
	  sprintf (OutputParameterFileName, "j_%f_djz_%f_chi_%f", Manager.GetDouble("j-value"), Manager.GetDouble("djz-value"), Manager.GetDouble("chi-value"));
	}
      else
	{
	  sprintf (OutputParameterFileName, "j_%f_djz_%f_chi_%f_hz_%f", Manager.GetDouble("j-value"), Manager.GetDouble("djz-value"), Manager.GetDouble("chi-value"), Manager.GetDouble("hz-value"));
	}
    }
    
  char* FullOutputFileName = new char [strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
  sprintf (FullOutputFileName, "%s_%s.dat", OutputFileName, OutputParameterFileName);
  double* JValues = new double [NbrSpins - 1];
  JValues[0] = Manager.GetDouble("j-value");
  for (int i = 1; i < (NbrSpins - 1); ++i)
    JValues[i] = JValues[0];
  double* JzValues = new double [NbrSpins - 1];
  double TmpDeltaJz = Manager.GetDouble("djz-value");
  for (int i = 0; i < (NbrSpins - 1); ++i)
    JzValues[i] = JValues[i] + TmpDeltaJz;
  double* ChiValues = new double [NbrSpins - 1];
  ChiValues[0] = Manager.GetDouble("chi-value");
  for (int i = 1; i < (NbrSpins - 1); ++i)
    ChiValues[i] = ChiValues[0];
  double* HzValues = 0;
  if (Manager.GetDouble("hz-value") != 0)
    {
      HzValues = new double [NbrSpins];
      HzValues[0] = Manager.GetDouble("hz-value");
      for (int i = 1; i < NbrSpins; ++i)
	HzValues[i] = HzValues[0];
    }

  int MaxSzValue = NbrSpins * SpinValue;
  int InitalSzValue = MaxSzValue & 1;
  if  (Manager.GetDouble("hz-value") != 0)
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
      
      SpinChainTripleProductHamiltonian* Hamiltonian = 0;
      if (HzValues == 0)
	Hamiltonian = new SpinChainTripleProductHamiltonian(Chain, NbrSpins, JValues, JzValues, ChiValues);
      else
	Hamiltonian = new SpinChainTripleProductHamiltonian(Chain, NbrSpins, JValues, JzValues, ChiValues, HzValues);
      char* TmpSzString = new char[64];
      char* TmpEigenstateString = new char[strlen(OutputFileName) + strlen(OutputParameterFileName) + 64];
      sprintf (TmpSzString, "%d", InitalSzValue);
      sprintf (TmpEigenstateString, "%s_%s_sz_%d", OutputFileName, OutputParameterFileName, InitalSzValue);
      GenericComplexMainTask Task(&Manager, Chain, &Lanczos, Hamiltonian, TmpSzString, CommentLine, 0.0,  FullOutputFileName,
				  FirstRun, TmpEigenstateString);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      FirstRun = false;
      delete Hamiltonian;
      delete Chain;
      delete[] TmpSzString;
      delete[] TmpEigenstateString;
    }
  delete[] OutputFileName;
  delete[] CommentLine;
  delete[] JValues;
  delete[] FullOutputFileName;
  delete[] JzValues;
  return 0;
}
