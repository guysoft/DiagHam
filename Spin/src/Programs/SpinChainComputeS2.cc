#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Operator/SpinS2Operator.h"
#include "Operator/SpinWith1DTranslationS2Operator.h"

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1Chain.h"
#include "HilbertSpace/Spin1_2ChainFull.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"
#include "HilbertSpace/Spin1_2ChainFullAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainFullInversionAnd2DTranslation.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "Options/Options.h"

#include "GeneralTools/Endian.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;




int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("SpinChainComputeS2" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "state whose total S^2 has to be computed");  
  (*SystemGroup) += new BooleanOption  ('c', "complex", "consider complex wave function");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinChainComputeS2 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("input-state") == 0)
    {
      cout << "error, an eigenstate file should be provided. See man page for option syntax or type SpinChainComputeS2 -h" << endl;
      return -1;
    }

  int SpinValue = 0;
  int NbrSpins = 0;
  int TotalSz = 0;
  bool SzFlag = true;
  bool Momentum2DFlag = false;
  bool Momentum1DFlag = false;
  bool InversionFlag = false;  
  int XMomentum = 0;
  int XPeriodicity = 0;
  int YMomentum = 0;
  int YPeriodicity = 0;
  int InversionSector = 0;

  if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSpins, TotalSz, SpinValue, XMomentum, XPeriodicity,
							    YMomentum, YPeriodicity) == false)
    {
      if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSpins, SpinValue, XMomentum, XPeriodicity, 
								YMomentum, YPeriodicity) == false)
	{
	  if (SpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSpins, TotalSz, SpinValue, XMomentum) == false)
	    {
	      if (SpinFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSpins, TotalSz, SpinValue) == false)
		{
		  SzFlag = false;
		  if (SpinFindSystemInfoFromFileName(Manager.GetString("input-state"), NbrSpins, SpinValue) == false)
		    {
		      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
		      return -1;
		    }
		}
	    }
	  else
	    {
	      XPeriodicity = NbrSpins;
	      Momentum1DFlag = true;	      
	    }
	}
      else
	{
	  Momentum2DFlag = true;
	  SzFlag = false;
	  InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSpins, SpinValue, XMomentum, XPeriodicity, 
											 YMomentum, YPeriodicity, InversionSector);
	}
    }
  else
    {
      InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSpins, TotalSz, SpinValue, XMomentum, XPeriodicity,
										     YMomentum, YPeriodicity, InversionSector);
      SzFlag = true;
      Momentum2DFlag = true;
    }
  if ((Momentum2DFlag == false) && (Momentum1DFlag == false))
    {
      if (SzFlag == true)
	cout << "N=" << NbrSpins << " Sz=" <<  TotalSz << " 2s=" << SpinValue << endl;
      else
	cout << "N=" << NbrSpins << " 2s=" << SpinValue << endl;
    }
  else
    {
      if (Momentum1DFlag == false)
	{
	  if (SzFlag == true)
	    {
	      if ((InversionFlag == false) || (InversionSector == 0))
		{
		  cout << "N=" << NbrSpins << "=" << XPeriodicity << "x" << YPeriodicity << " Sz=" <<  TotalSz << " 2s=" << SpinValue << endl;
		}
	      else
		{
		  cout << "N=" << NbrSpins << "=" << XPeriodicity << "x" << YPeriodicity << " Sz=" <<  TotalSz << " 2s=" << SpinValue << " I=" << InversionSector << endl;
		}
	    }
	  else
	    {
	      if ((InversionFlag == false) || (InversionSector == 0))
		{
		  cout << "N=" << NbrSpins << "=" << XPeriodicity << "x" << YPeriodicity << " 2s=" << SpinValue << endl;
		}
	      else
		{
		  cout << "N=" << NbrSpins << "=" << XPeriodicity << "x" << YPeriodicity << " 2s=" << SpinValue << " I=" << InversionSector << endl;
		}
	    }
	}
      else
	{
	  if (SzFlag == true)
	    cout << "N=" << NbrSpins << " Sz=" <<  TotalSz << " 2s=" << SpinValue << " kx=" << XMomentum << endl;
	  else
	    cout << "N=" << NbrSpins << " 2s=" << SpinValue << " kx=" << XMomentum  << endl;
	}
   }

  if ((Momentum2DFlag == false) && (Momentum1DFlag == false))
    {
      AbstractSpinChain* Space;
      
      if (SzFlag == true)
	{
	  switch (SpinValue)
	    {
	    case 1 :
	      Space = new Spin1_2Chain (NbrSpins, TotalSz, 1000000);
	      break;
	    case 2 :
	      Space = new Spin1Chain (NbrSpins, TotalSz, 1000000);
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
	      Space = new Spin1_2ChainFull (NbrSpins);
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
      SpinS2Operator TmpOperator(Space, NbrSpins);
      if (Manager.GetBoolean("complex") == false)
	{
	  RealVector TmpState;
	  if (TmpState.ReadVector (Manager.GetString("input-state")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
	      return -1;      
	    }
	  
	  Complex TmpS2 = TmpOperator.MatrixElement(TmpState, TmpState);
	  cout << "S^2 = " << TmpS2.Re << endl; 
	  cout << "2S = " << (sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0) << " " << round(sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0) << endl; 
	}
      else
	{
	  ComplexVector TmpState;
	  if (TmpState.ReadVector (Manager.GetString("input-state")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
	      return -1;      
	    }
	  
	  Complex TmpS2 = TmpOperator.MatrixElement(TmpState, TmpState);
	  cout << "<S^2>=" << TmpS2.Re << " <S>=" << (0.5 * (sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0)) << endl;
	  cout << "round(<2S>)=" << round(sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0) << endl; 
	}
    }
  else
    {      
      if (Momentum1DFlag == false)
	{
	  AbstractSpinChain* TmpSpace;
	  if ((InversionFlag == false) || (InversionSector == 0))
	    {
	      if (SzFlag == true)
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
		    case 1 :
		      TmpSpace = new Spin1_2ChainFullAnd2DTranslation (XMomentum, XPeriodicity, YMomentum, YPeriodicity);
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
	      if (SzFlag == true)
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
		    case 1 :
		      TmpSpace = new Spin1_2ChainFullInversionAnd2DTranslation (InversionSector, XMomentum, XPeriodicity, YMomentum, YPeriodicity);
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
	}
      else
	{
	  AbstractSpinChainWithTranslations* Space = 0;
	  if (SzFlag == true)
	    {
	      switch (SpinValue)
		{
		case 1 :
		  Space = new Spin1_2ChainWithTranslations (NbrSpins, XMomentum, 1, TotalSz, 1000000, 1000000);
		  break;
		case 2 :
		  Space = new Spin1ChainWithTranslations (NbrSpins, XMomentum, TotalSz, 1000000, 1000000);
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
	    }
	  SpinWith1DTranslationS2Operator TmpOperator(Space, NbrSpins);
	  ComplexVector TmpState;
	  if (TmpState.ReadVector (Manager.GetString("input-state")) == false)
	    {
	      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
	      return -1;      
	    }
	  
	  Complex TmpS2 = TmpOperator.MatrixElement(TmpState, TmpState);
	  cout << "<S^2>=" << TmpS2.Re << " <S>=" << (0.5 * (sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0)) << endl;
	  cout << "round(<2S>)=" << round(sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0) << endl; 
	}
    }
  
  return 0;
}

