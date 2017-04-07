#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "HilbertSpace/Spin1_2ChainWithPseudospin.h"
#include "HilbertSpace/Spin1_2ChainWithPseudospinAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation.h"

// #include "Operator/BondEnergySpinPseudospinOperator.h"
#include "Operator/SpinWithPseudospin2DTranslationSpinSpinCorrelationOperator.h"
#include "Operator/SpinWith2DTranslationBondBondCorrelationOperator.h"


#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


// get a linearized position index from the 2d coordinates
//
// xPosition = position along the x direction
// yPosition = position along the y direction
// return value = linearized index

int GetLinearizedIndex(int xPosition, int yPosition, int nbrSitesX, int nbrSitesY);

int main(int argc, char** argv)
{
  OptionManager Manager ("SpinComputeBondEnergyTriangleLatticeProjectedFromKagome" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  
  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "name of the file corresponding to the state |Psi_R> (the order parameter being <Psi_L|c^+c^+|Psi_R>");   (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
//   (*SystemGroup) += new BooleanOption ('\n', "only-cc", "compute only the parameters c^+_sigma c^+_sigma' instead of their linear combinations");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardSuperconductorOrderParameter -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSites = 0;
  int XMomentum = 0;
  int YMomentum = 0;
  int InversionSector = 0;
  int SzValue = 0;
  int SzParitySector = 0;
  bool TotalSpinConservedFlag = true;
  bool InversionFlag = false;
  bool SzSymmetryFlag = false;
  int XPeriodicity = 0;
  int YPeriodicity = 0;
  bool Statistics = true;
  int NbrInputStates = 0;
  int SpinValue = 1;
  int Offset = 0;
  
  if (Manager.GetString("input-state") == 0)
    {
      cout << "error, an input file has to be provided. See man page for option syntax or type SpinSystemConvertFromTranslationInvariantBasis -h" << endl;
      return -1;
    }

  else
    {
      InversionSector = 0;
      SzParitySector = 0;
      if (IsFile(Manager.GetString("input-state")) == false)
	{
	  cout << "can't open file " << Manager.GetString("input-state") << endl;
	  return -1;
	}   
	
      SpinFindSystemInfoFromVectorFileName (Manager.GetString("input-state"), NbrSites, SzValue, SpinValue, XMomentum, InversionSector, SzParitySector, Offset);
      
//       (char* filename, int& nbrSpins, int& sz, int& spin, int& momentum, int& inversion, int& szSymmetry, int& offset)
      cout << Offset << endl;
      if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue, SpinValue, XMomentum, XPeriodicity,
								YMomentum, YPeriodicity) == false)
	{
	  TotalSpinConservedFlag = false;
	  if (SpinWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SpinValue, XMomentum, XPeriodicity, 
								    YMomentum, YPeriodicity) == false)
	    {
	      cout << "error while retrieving system parameters from file name " <<Manager.GetString("input-state")  << endl;
	      return -1;
	    }
	  InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SpinValue, XMomentum, XPeriodicity, 
											 YMomentum, YPeriodicity, InversionSector);
	  if (InversionFlag == false)
	    cout << "error" << endl;
	  else
	    cout << "OK" << " " << InversionSector << endl;
	}
      else
	{
	  InversionFlag = SpinWith2DTranslationInversionFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrSites, SzValue, SpinValue, XMomentum, XPeriodicity,
											 YMomentum, YPeriodicity, InversionSector);
	  if (SzValue == 0)
	  {
	    SpinFindSystemInfoFromVectorFileName((Manager.GetString("input-state")), NbrSites, SzValue, SpinValue, InversionSector, SzParitySector);
	    cout << NbrSites << " " << SzValue << " " << SpinValue << " " << XMomentum << " " << XPeriodicity << " " << YMomentum << " " <<  YPeriodicity << " " << SzParitySector << endl;
	    if (SzParitySector != 0)
	      SzSymmetryFlag = true;
	  }
	}
    }
   
      
  ComplexVector State;
  if (Manager.GetString("input-state") != 0)
    {
      if (State.ReadVector (Manager.GetString("input-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("input-state") << endl;
	  return -1;      
	}
    }


  cout << (State.GetVectorDimension()) << endl;
  Spin1_2ChainWithPseudospinAnd2DTranslation* Space = 0;
//   cout << NbrSites << " " << SzValue << " " << SzParitySector << " " << XMomentum << " " << XPeriodicity << " " << YMomentum << " " << YPeriodicity << endl;
  Space = new Spin1_2ChainWithPseudospinSzSymmetryAnd2DTranslation(NbrSites, SzValue, SzParitySector, XMomentum, XPeriodicity, YMomentum, YPeriodicity, 1000000);

  if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
    {
      cout << "error, " << Manager.GetString("input-state")  << " has a wrong dimension (" << State.GetVectorDimension() << ", should be " << Space->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  
 
  ofstream File;
  char* OutputFileName;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      File.open(OutputFileName, ios::binary | ios::out);
    }
  else
    {
      if (Manager.GetString("input-state") != 0)
	{
	  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("input-state"), "vec", "spinspin.dat");
	  if (OutputFileName == 0)
	    {
	      cout << "no vec extension was find in " << Manager.GetString("input-state") << " file name" << endl;
	      return 0;
	    }
	  File.open(OutputFileName, ios::binary | ios::out);
	}
    }
  File.precision(14);
  cout.precision(14);
  
  File << "# i j <S_{0,0} S_{i,j}}>" << endl;
  AbstractOperator* Operator;
  Complex SpinSpinCorrelations;
  
  for (int i = 0; i < XPeriodicity; ++i)
  {
    for (int j = 0; j < YPeriodicity; ++j)
    {
      SpinSpinCorrelations = 0.0;
      for (int nx = 0; nx < XPeriodicity; ++nx)
      {
	for (int ny = 0; ny < YPeriodicity; ++ny)
	{
	  int TmpIndex1 = GetLinearizedIndex(nx, ny, XPeriodicity, YPeriodicity);
	  int TmpIndex2 = GetLinearizedIndex(nx + i, ny + j, XPeriodicity, YPeriodicity);
	  
	  Operator = new SpinWithPseudospin2DTranslationSpinSpinCorrelationOperator(Space, XMomentum, XPeriodicity, YMomentum, YPeriodicity, TmpIndex1, TmpIndex2);
	  OperatorMatrixElementOperation Operation(Operator, State, State, State.GetVectorDimension());
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  SpinSpinCorrelations += Operation.GetScalar();
// 	cout << NeighborSpinSpinCorrelation << " " ;
	  delete Operator;
	}
      }
      File << i << " " << j << " " << (SpinSpinCorrelations.Re / (XPeriodicity * YPeriodicity)) << endl;
    }
  }

  File.close();
  
  
  
  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("input-state"), "vec", "bondbond.dat");
  if (OutputFileName == 0)
  {
    cout << "no vec extension was find in " << Manager.GetString("input-state") << " file name" << endl;
    return 0;
  }
  File.open(OutputFileName, ios::binary | ios::out);
  File.precision(14);
  cout.precision(14);
  
  File << "# l_0 i j l_{i,j} <S_{0,0} S_{1,0} S_{i,j}S_{i,j,l}}>" << endl;
  int* TmpIndex2b = new int[3];
  Complex* BondBondCorrelations = new Complex[3];
  Complex TmpOnSite;
//   for (int m = 0; m < 3; ++m)
//   {
  for (int i = 0; i < XPeriodicity; ++i)
  {
    for (int j = 0; j < YPeriodicity; ++j)
    {
      TmpOnSite = 0.0;
      for (int l = 0; l < 3; ++l)
	BondBondCorrelations[l] = 0.0;
      for (int nx = 0; nx < XPeriodicity; ++nx)
      {
	for (int ny = 0; ny < YPeriodicity; ++ny)
	{
	  int TmpIndex1a = GetLinearizedIndex(nx, ny, XPeriodicity, YPeriodicity);
	  int TmpIndex1b = GetLinearizedIndex(nx + Offset, ny + 1, XPeriodicity, YPeriodicity);
	  
	  int TmpIndex2a = GetLinearizedIndex(nx + i, ny + j, XPeriodicity, YPeriodicity);
	  TmpIndex2b[0] = GetLinearizedIndex(nx + i + 1, ny + j, XPeriodicity, YPeriodicity);
	  TmpIndex2b[1] = GetLinearizedIndex(nx + i + Offset, ny + j + 1, XPeriodicity, YPeriodicity);
	  TmpIndex2b[2] = GetLinearizedIndex(nx + i + 1 - Offset, ny + j - 1, XPeriodicity, YPeriodicity);
	  
	
	  Operator = new SpinWith2DTranslationBondBondCorrelationOperator(Space, XMomentum, XPeriodicity, YMomentum, YPeriodicity, TmpIndex1a, TmpIndex1b, TmpIndex1a, TmpIndex1b);
	  OperatorMatrixElementOperation Operation(Operator, State, State, State.GetVectorDimension());
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  TmpOnSite += Operation.GetScalar();
// 	cout << NeighborSpinSpinCorrelation << " " ;
	  delete Operator;
	    
	  
	  for (int l = 0; l < 3; ++l)
	  {
	    Operator = new SpinWith2DTranslationBondBondCorrelationOperator(Space, XMomentum, XPeriodicity, YMomentum, YPeriodicity, TmpIndex1a, TmpIndex1b, TmpIndex2a, TmpIndex2b[l]);
	    OperatorMatrixElementOperation Operation(Operator, State, State, State.GetVectorDimension());
	    Operation.ApplyOperation(Architecture.GetArchitecture());
	    BondBondCorrelations[l] += Operation.GetScalar();
// 	cout << NeighborSpinSpinCorrelation << " " ;
	    delete Operator;
	  }
	}
      }
      cout << (TmpOnSite / (XPeriodicity * YPeriodicity)) << endl;
      for (int l = 0; l < 3; ++l)
	File << i << " " << j << " " << l << " " << (BondBondCorrelations[l].Re / (XPeriodicity * YPeriodicity)) << endl;
    }
  }

  File.close();
	  
	  
//   ofstream FileFourierTransform;
//   if ((RightMomentumFlag == true) && (LeftMomentumFlag == true))
//     {
//       char* TmpFileName = ReplaceExtensionToFileName(OutputFileName, "dat", "fourier.dat");
//       if (TmpFileName == 0)
// 	{
// 	  cout << "no dat extension was find in " << OutputFileName << " file name" << endl;
// 	  return 0;
// 	}
//       FileFourierTransform.open(TmpFileName, ios::binary | ios::out);
//       FileFourierTransform.precision(14);
//     }

 
//   if ((RightMomentumFlag == true) && (LeftMomentumFlag == true))
//     {
//       FileFourierTransform.close();
//     }
  return 0;
}


// get a linearized position index from the 2d coordinates
//
// xPosition = position along the x direction
// yPosition = position along the y direction
// return value = linearized index

int GetLinearizedIndex(int xPosition, int yPosition, int nbrSitesX, int nbrSitesY)
{
  if (xPosition < 0)
    xPosition += nbrSitesX;
  if (xPosition >= nbrSitesX)
    xPosition -= nbrSitesX;
  if (yPosition < 0)
    yPosition += nbrSitesY;
  if (yPosition >= nbrSitesY)
    yPosition -= nbrSitesY;
  return ((xPosition * nbrSitesY) + yPosition);
}