#include "Vector/ComplexVector.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "HilbertSpace/Spin1_2ChainNew.h"
#include "HilbertSpace/Spin1_2ChainWithPseudospin.h"



#include "Options/Options.h"
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <cstring> 

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("SpinProjectKagomeOntoTriangleLattice" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption('i', "input-state", "name of the file containing the state defined in the translation invariant basis");
  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinSystemConvertFromTranslationInvariantBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSites = 0;
  int* XMomentum = 0;
  int* YMomentum = 0;
  int* InversionSector = 0;
  int* SzValue = 0;
  int* SzParitySector = 0;
  bool TotalSpinConservedFlag = true;
  bool InversionFlag = false;
  bool SzSymmetryFlag = false;
  int XPeriodicity = 0;
  int YPeriodicity = 0;
  bool Statistics = true;
  int NbrInputStates = 0;
  int SpinValue = 1;

  if (Manager.GetString("input-state") == 0)
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type SpinSystemConvertFromTranslationInvariantBasis -h" << endl;
      return -1;
    }

  if (Manager.GetString("input-state") != 0)
    {
      NbrInputStates = 1;
      SzValue = new int [1];
      SzParitySector = new int[1];
      SzParitySector[0] = 0;
      if (IsFile(Manager.GetString("input-state")) == false)
	{
	  cout << "can't open file " << Manager.GetString("input-state") << endl;
	  return -1;
	}      
      if (SpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue[0], SpinValue) == false)
	{
	  TotalSpinConservedFlag = false;
	  if (SpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SpinValue, XMomentum[0], XPeriodicity, 
								    YMomentum[0], YPeriodicity) == false)
	    {
	      cout << "error while retrieving system parameters from file name " <<Manager.GetString("input-state")  << endl;
	      return -1;
	    }
	  
	}
      else
	{
	  if (SzValue[0] == 0)
	  {
	    SpinFindSystemInfoFromVectorFileName((Manager.GetString("input-state")), NbrSites, SzValue[0], SpinValue);
	    cout << NbrSites << " " << SzValue[0] << " " << SpinValue << " " << endl;
	    
	  }
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      NbrInputStates = DegenerateFile.GetNbrLines();
      SzValue = new int [NbrInputStates];
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  // 	if ((FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0,i), NbrParticles, NbrSites, XMomentum[i], YMomentum[i],  XPeriodicity, YPeriodicity, Statistics, GutzwillerFlag) == false) && (FTIHubbardModelWith1DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0,i), NbrParticles, NbrSites, XMomentum[i], XPeriodicity, Statistics, GutzwillerFlag) == false))
	  // 	  {
	  // 	    cout << "error while retrieving system parameters from file name " << DegenerateFile(0, i) << endl;
	  // 	    return -1;
	  // 	  }
	  // 	  TotalSpinConservedFlag = FTIHubbardModelWithSzFindSystemInfoFromVectorFileName (DegenerateFile(0, i), NbrParticles, NbrSites, SzValue[i], Statistics, GutzwillerFlag);
	}
      
    }
 
  cout << "Nbr sites=" << NbrSites << " Sz =" << SzValue[0] << endl;
 

  ComplexVector* InputStates = new ComplexVector[NbrInputStates];
  char** InputStateNames = new char*[NbrInputStates];
  if (Manager.GetString("input-state") != 0)
    {
      InputStateNames[0] = new char [strlen(Manager.GetString("input-state")) + 1];
      strcpy (InputStateNames[0], Manager.GetString("input-state"));
      if (InputStates[0].ReadVector (Manager.GetString("input-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("input-state") << endl;
	  return -1;      
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      if (InputStates[0].ReadVector (DegenerateFile(0, 0)) == false)
	{
	  cout << "can't open vector file " << DegenerateFile(0, 0) << endl;
	  return -1;      
	}	  
      InputStateNames[0] = new char [strlen(DegenerateFile(0, 0)) + 1];
      strcpy (InputStateNames[0], DegenerateFile(0, 0));
      for (int i = 1; i < NbrInputStates; ++i)
	{
	  if (InputStates[i].ReadVector (DegenerateFile(0, i)) == false)
	    {
	      cout << "can't open vector file " << DegenerateFile(0, i) << endl;
	      return -1;      
	    }	  
	  if (InputStates[0].GetVectorDimension() != InputStates[i].GetVectorDimension())
	    {
	      cout << "error, " << DegenerateFile(0, 0) << " and " <<  DegenerateFile(0, i) << " don't have the same  dimension (" << InputStates[0].GetVectorDimension() << " and " << InputStates[i].GetVectorDimension()<< ")" << endl;
	      return -1;
	    }
	  InputStateNames[i] = new char [strlen(DegenerateFile(0, i)) + 1];
	  strcpy (InputStateNames[i], DegenerateFile(0, i));
	}
    }


    Spin1_2ChainNew* InputSpace = new Spin1_2ChainNew (NbrSites, SzValue[0], 1000000);
//   AbstractSpinChain** InputSpace = new AbstractSpinChain* [NbrInputStates];
//   for (int i = 0; i < NbrInputStates; ++i)
//     {
//       if (TotalSpinConservedFlag == false)
//       {
// 	cout << "spin not conserved not available" << endl;
// 	return -1;	      
//       }
//       else
//       {
// 	switch (SpinValue)
// 	{
// 	  case 1:
// 	    InputSpace[i] = new Spin1_2ChainNew (NbrSites, SzValue[0], 1000000);
// 	    break;
// 	  default :
// 	  {
// 	    if ((SpinValue & 1) == 0)
// 	      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
// 	    else 
// 	      cout << "spin " << SpinValue << "/2 are not available" << endl;
// 	    return -1;
// 	  }
// 	}
//       }
//     }
  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (InputSpace->GetHilbertSpaceDimension() != InputStates[i].GetVectorDimension())
	{
	  cout << "error, dimension mismatch for vector " << i << " (" << InputStates[i].GetVectorDimension() << ", should be " << InputSpace->GetHilbertSpaceDimension() << ")" << endl;
	  return -1;
	}
    }
  
  Spin1_2ChainWithPseudospin* OutputSpace = new Spin1_2ChainWithPseudospin ((NbrSites / 3), SzValue[0], 1000000);
//   AbstractSpinChain* OutputSpace = new AbstractSpinChain;
//   for (int i = 0; i < NbrInputStates; ++i)
//     {
//       if (TotalSpinConservedFlag == false)
// 	{
// 	}
//       else
// 	{
// 	  switch (SpinValue)
// 	    {
// 	    case 1 :
// 	      OutputSpace[i] = new Spin1_2ChainWithPseudospin (NbrSites, SzValue[i], 1000000);
// 	      break;
// 	    default :
// 	      {
// 		if ((SpinValue & 1) == 0)
// 		  cout << "spin " << (SpinValue / 2) << " are not available" << endl;
// 		else 
// 		  cout << "spin " << SpinValue << "/2 are not available" << endl;
// 		return -1;
// 	      }
// 	    }
// 	}
//     }
 
  for (int i = 0; i < NbrInputStates; ++i)
    {      
      char* VectorOutputName = new char[strlen(InputStateNames[i]) + 64];
      char* OutputName = new char[strlen(InputStateNames[i]) + 64];
      sprintf(VectorOutputName, "%s.proj", InputStateNames[i]);
      sprintf(OutputName, "%s.proj.norm", InputStateNames[i]);
      RealVector TmpVector = OutputSpace->ProjectToEffectiveSubspaceThreeToOne(InputStates[i], InputSpace);
      
      ofstream File;
      File.open(OutputName, ios::binary | ios::out);
      File.precision(14);
      File << "#Norm projected Vector" << endl << TmpVector.Norm();
      cout << "Projected vector norm = " << TmpVector.Norm() << endl;
      File.close();
      
      TmpVector.Normalize();
      if (TmpVector.WriteVector(VectorOutputName) == false)
	{
	  cout << "error, can't write vector " << VectorOutputName << endl;
	}

      delete[] VectorOutputName;
      delete[] OutputName;
    }
}
