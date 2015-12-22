#include "Vector/ComplexVector.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "GeneralTools/StringTools.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1_2ChainFullAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainFullInversionAnd2DTranslation.h"


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
  OptionManager Manager ("SpinSystemConvertFromTranslationInvariantBasis" , "0.01");
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
  bool TotalSpinConservedFlag = false;
  bool InversionFlag = false;
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
      XMomentum = new int[1];
      YMomentum = new int[1];
      InversionSector = new int[1];
      InversionSector[0] = 0;
      SzValue = new int [1];
      if (IsFile(Manager.GetString("input-state")) == false)
	{
	  cout << "can't open file " << Manager.GetString("input-state") << endl;
	  return -1;
	}      
      if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue[0], SpinValue, XMomentum[0], XPeriodicity,
								YMomentum[0], YPeriodicity) == false)
	{
	  TotalSpinConservedFlag = false;
	  if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SpinValue, XMomentum[0], XPeriodicity, 
								    YMomentum[0], YPeriodicity) == false)
	    {
	      cout << "error while retrieving system parameters from file name " <<Manager.GetString("input-state")  << endl;
	      return -1;
	    }
	  InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SpinValue, XMomentum[0], XPeriodicity, 
											 YMomentum[0], YPeriodicity, InversionSector[0]);
	  if (InversionFlag == false)
	    cout << "error" << endl;
	  else
	    cout << "OK" << " " << InversionSector[0] << endl;
	}
      else
	{
	  InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue[0], SpinValue, XMomentum[0], XPeriodicity,
											 YMomentum[0], YPeriodicity, InversionSector[0]);
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
      XMomentum = new int [NbrInputStates];
      YMomentum = new int [NbrInputStates];
      InversionSector = new int[NbrInputStates];
      SzValue = new int [NbrInputStates];
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  InversionSector[i] = 0;
	  // 	if ((FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0,i), NbrParticles, NbrSites, XMomentum[i], YMomentum[i],  XPeriodicity, YPeriodicity, Statistics, GutzwillerFlag) == false) && (FTIHubbardModelWith1DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0,i), NbrParticles, NbrSites, XMomentum[i], XPeriodicity, Statistics, GutzwillerFlag) == false))
	  // 	  {
	  // 	    cout << "error while retrieving system parameters from file name " << DegenerateFile(0, i) << endl;
	  // 	    return -1;
	  // 	  }
	  // 	  TotalSpinConservedFlag = FTIHubbardModelWithSzFindSystemInfoFromVectorFileName (DegenerateFile(0, i), NbrParticles, NbrSites, SzValue[i], Statistics, GutzwillerFlag);
	}
      
    }
  
  if (YPeriodicity == 0)
    {
      cout << "Convert from kx basis to real space basis" << endl;
      for (int i = 0; i < NbrInputStates; ++i)
	cout << " Nbr sites=" << NbrSites << " Kx=" << XMomentum[i] << endl;
    }
  else
    {
      cout << "Convert from (kx,ky) basis to real space basis" << endl;
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  if ((InversionFlag == false) || (InversionSector[i] == 0))
	    {
	      cout << " Nbr sites=" << NbrSites << "=" << XPeriodicity << "x" << YPeriodicity << " Kx=" << XMomentum[i] << " Ky=" << YMomentum[i] << endl;
	    }
	  else
	    {
	      cout << " Nbr sites=" << NbrSites << "=" << XPeriodicity << "x" << YPeriodicity << " Kx=" << XMomentum[i] << " Ky=" << YMomentum[i] << " I=" << InversionSector[i] << endl;
	    }
	}
    }

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


  AbstractSpinChain** InputSpace = new AbstractSpinChain* [NbrInputStates];
  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (YPeriodicity == 0)
	{
	  if (TotalSpinConservedFlag == false)
	    {
	      switch (SpinValue)
		{
		default :
		  {
		    if ((SpinValue & 1) == 0)
		      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	    }
	  else
	    {
	      switch (SpinValue)
		{
		default :
		  {
		    if ((SpinValue & 1) == 0)
		      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	    }
	}
      else
	{
	  if (TotalSpinConservedFlag == false)
	    {
	      if ((InversionFlag == false) || (InversionSector[i] == 0))
		{
		  switch (SpinValue)
		    {
		    case 1 :
		      InputSpace[i] = new Spin1_2ChainFullAnd2DTranslation (XMomentum[i], XPeriodicity, YMomentum[i], YPeriodicity);
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
		}
	      else
		{
		  switch (SpinValue)
		    {
		    case 1 :
		      InputSpace[i] = new Spin1_2ChainFullInversionAnd2DTranslation (InversionSector[i], XMomentum[i], XPeriodicity, YMomentum[i], YPeriodicity);
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
		}
	    }
	  else
	    {
	      switch (SpinValue)
		{
		default :
		  {
		    if ((SpinValue & 1) == 0)
		      cout << "spin " << (SpinValue / 2) << " are not available" << endl;
		    else 
		      cout << "spin " << SpinValue << "/2 are not available" << endl;
		    return -1;
		  }
		}
	    }
	}
    }
  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (InputSpace[i]->GetHilbertSpaceDimension() != InputStates[i].GetVectorDimension())
	{
	  cout << "error, dimension mismatch for vector " << i << " (" << InputStates[i].GetVectorDimension() << ", should be " << InputSpace[i]->GetHilbertSpaceDimension() << ")" << endl;
	  return -1;
	}
    }
  
  AbstractSpinChain** OutputSpace = new AbstractSpinChain* [NbrInputStates];
  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (TotalSpinConservedFlag == false)
	{
	  switch (SpinValue)
	    {
	    case 1 :
	      OutputSpace[i] = new Spin1_2ChainFull (NbrSites);
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
	}
      else
	{
	  switch (SpinValue)
	    {
	    case 1 :
	      OutputSpace[i] = new Spin1_2Chain (NbrSites, SzValue[i], 1000000);
	      break;
	    case 2 :
	      OutputSpace[i] = new Spin1Chain (NbrSites, SzValue[i], 1000000);
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
	}
    }

  char* XPeriodicityString = new char[256];
  char* XMomentumString = new char[256];
  
  
  for (int i = 0; i < NbrInputStates; ++i)
    {
      if (YPeriodicity == 0)
	{
	  sprintf (XPeriodicityString, "_x_%d", XPeriodicity);
	  sprintf (XMomentumString, "_kx_%d", XMomentum[i]);
	}
      else
	{
	  sprintf (XPeriodicityString, "_x_%d_y_%d", XPeriodicity, YPeriodicity);
	  sprintf (XMomentumString, "_kx_%d_ky_%d", XMomentum[i], YMomentum[i]);
	}
      
      char* VectorOutputName = ReplaceString(InputStateNames[i], XMomentumString, "");
      ComplexVector TmpVector = InputSpace[i]->ConvertFromKxKyBasis(InputStates[i], OutputSpace[i]);
      char* Extension = new char[256];
      sprintf (Extension, "%d.vec", i);
      VectorOutputName = ReplaceExtensionToFileName(VectorOutputName, "vec", Extension);
      if (TmpVector.WriteVector(VectorOutputName) == false)
	{
	  cout << "error, can't write vector " << VectorOutputName << endl;
	}
      delete[] Extension;
      delete[] VectorOutputName;
    }
}
