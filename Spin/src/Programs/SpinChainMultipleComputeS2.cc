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
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1_2ChainWithTranslationsAndSzInversionSymmetries.h"
#include "HilbertSpace/Spin1ChainWithTranslations.h"
#include "HilbertSpace/Spin1_2ChainFullAnd2DTranslation.h"
#include "HilbertSpace/Spin1_2ChainFullInversionAnd2DTranslation.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndInversionSymmetry.h"
#include "HilbertSpace/Spin1ChainWithTranslationsAndSzInversionSymmetries.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/OperatorMultipleMatrixElementOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "Options/Options.h"

#include "GeneralTools/Endian.h"
#include "GeneralTools/StringTools.h"

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
  OptionManager Manager ("SpinChainMultipleComputeS2" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "multiple-states", "provide as a matrix a series of states whose S^2 computed");  
  (*SystemGroup) += new SingleStringOption  ('\n', "spectrum", "provide the spectrum to automatically look for degenerate states");  
  (*SystemGroup) += new BooleanOption  ('c', "complex", "consider complex wave function");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type SpinChainMultipleComputeS2 -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("multiple-states") == 0)
    {
      cout << "error, an eigenstate file should be provided. See man page for option syntax or type SpinChainMultipleComputeS2 -h" << endl;
      return -1;
    }

  int SpinValue = 0;
  int NbrSpins = 0;
  int TotalSz = 0;
  bool SzFlag = true;
  bool Momentum1DFlag = false;
  bool InversionFlag = false;  
  bool SzSymmetryFlag = false;  
  int XMomentum = 0;
  int XPeriodicity = 0;
  int InversionSector = 0;
  int SzSymmetrySector = 0;
  double Error = 1e-12;
 
  if (SpinFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, TotalSz, SpinValue, XMomentum, InversionSector, SzSymmetrySector) == false)
    {
      if (SpinFindSystemInfoFromVectorFileName(Manager.GetString("multiple-states"), NbrSpins, TotalSz, SpinValue) == false)
	{
	  SzFlag = false;
	  if (SpinFindSystemInfoFromFileName(Manager.GetString("multiple-states"), NbrSpins, SpinValue) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << Manager.GetString("multiple-states") << endl;
	      return -1;
	    }
	}
    }
  else
    {
      XPeriodicity = NbrSpins;
      Momentum1DFlag = true;	      
      if (InversionSector != 0)
	InversionFlag = true;
      if (SzSymmetrySector != 0)
	SzSymmetryFlag = true;
    }
  
  if (Momentum1DFlag == false)
    {
      if (SzFlag == true)
	cout << "N=" << NbrSpins << " Sz=" <<  TotalSz << " 2s=" << SpinValue << endl;
      else
	cout << "N=" << NbrSpins << " 2s=" << SpinValue << endl;
    }
  else
    {
      if (SzFlag == true)
	cout << "N=" << NbrSpins << " Sz=" <<  TotalSz << " 2s=" << SpinValue << " kx=" << XMomentum;
      else
	cout << "N=" << NbrSpins << " 2s=" << SpinValue << " kx=" << XMomentum;
      if ((InversionFlag == true) && (InversionSector != 0))
	{
	  cout << " inversion=" << InversionSector;
	}
      if ((SzSymmetryFlag == true) && (SzSymmetrySector != 0))
	{
	  cout << " Sz<->-Sz=" << SzSymmetrySector;
	}
      cout << endl;
    }
  
  RealMatrix RealEigenstates;
  ComplexMatrix ComplexEigenstates;
  int NbrStates = 0;
  if (Manager.GetBoolean("complex") == false) 
    {
      if (RealEigenstates.ReadMatrix(Manager.GetString("multiple-states")) == false)
	{
	  cout << "can't read " << Manager.GetString("multiple-states") << endl;
	}
      NbrStates = RealEigenstates.GetNbrColumn();
    }
  else
    {
      if (ComplexEigenstates.ReadMatrix(Manager.GetString("multiple-states")) == false)
	{
	  cout << "can't read " << Manager.GetString("multiple-states") << endl;
	}
      NbrStates = ComplexEigenstates.GetNbrColumn();
    }
  
  int* Degeneracies = 0;
  int NbrDegeneracyIndices = 0;
  double* Spectrum = 0;
  if (Manager.GetString("spectrum") == 0)
    {
      int NbrDegeneracyIndices = NbrStates;
      Degeneracies = new int [NbrDegeneracyIndices];
      for (int i = 0; i < NbrDegeneracyIndices; ++i)
	Degeneracies[i] = 1;
    }
  else
    {
      MultiColumnASCIIFile SpectrumFile;
      if (SpectrumFile.Parse(Manager.GetString("spectrum")) == false)
	{
	  SpectrumFile.DumpErrors(cout);
	  return -1;
	}      
      Spectrum = new double[NbrStates];
      int TotalNbrEnergies = SpectrumFile.GetNbrLines();
      if (Momentum1DFlag == false)
	{
	  int* TmpSzValues = SpectrumFile.GetAsIntegerArray(0);
	  double* TmpEnergies = SpectrumFile.GetAsDoubleArray(1);
	  int TmpIndex = 0; 
	  while ((TmpIndex < TotalNbrEnergies) && (TmpSzValues[TmpIndex] != TotalSz))
	    ++TmpIndex;
	  if (TmpIndex == TotalNbrEnergies)
	    {
	      cout << "error, the spectrum has no eigenvalues corresponding to the required quantum numbers" << endl;
	      return 0;
	    }
	  int TmpIndex2 = 0;
	  while ((TmpIndex < TotalNbrEnergies) && (TmpIndex2 < NbrStates) && (TmpSzValues[TmpIndex] == TotalSz))
	    {
	      Spectrum[TmpIndex2] = TmpEnergies[TmpIndex];
	      ++TmpIndex2;
	      ++TmpIndex;
	    }
	  if (TmpIndex2 < NbrStates)
	    {
	      cout << "error, the spectrum has less eigenvalues corresponding to the required quantum numbers (" 
		   << TmpIndex2 << " vs " << NbrStates << ")" << endl;
	      return 0;
	    }
	}
      else
	{
	  if (InversionFlag == true)
	    {
	      if (SzSymmetryFlag == true)
		{
		  int* TmpSzValues = SpectrumFile.GetAsIntegerArray(0);
		  int* TmpKValues = SpectrumFile.GetAsIntegerArray(1);
		  int* TmpInvValues = SpectrumFile.GetAsIntegerArray(3);
		  int* TmpSzSymValues = SpectrumFile.GetAsIntegerArray(2);
		  double* TmpEnergies = SpectrumFile.GetAsDoubleArray(4);
		  int TmpIndex = 0; 
		  while ((TmpIndex < TotalNbrEnergies) && 
			 ((TmpSzValues[TmpIndex] != TotalSz) || (TmpKValues[TmpIndex] != XMomentum) || 
			  (TmpInvValues[TmpIndex] != InversionSector) || (TmpSzSymValues[TmpIndex] != SzSymmetrySector)))
		    ++TmpIndex;
		  if (TmpIndex == TotalNbrEnergies)
		    {
		      cout << "error, the spectrum has no eigenvalues corresponding to the required quantum numbers" << endl;
		      return 0;
		    }
		  int TmpIndex2 = 0;
		  while ((TmpIndex < TotalNbrEnergies) && (TmpIndex2 < NbrStates) && 
			 ((TmpSzValues[TmpIndex] == TotalSz) && (TmpKValues[TmpIndex] == XMomentum) && 
			  (TmpInvValues[TmpIndex] == InversionSector) && (TmpSzSymValues[TmpIndex] == SzSymmetrySector)))
		    {
		      Spectrum[TmpIndex2] = TmpEnergies[TmpIndex];
		      ++TmpIndex2;
		      ++TmpIndex;
		    }
		  if (TmpIndex2 < NbrStates)
		    {
		      cout << "error, the spectrum has less eigenvalues corresponding to the required quantum numbers (" 
			   << TmpIndex2 << " vs " << NbrStates << ")" << endl;
		      return 0;
		    }
		}
	      else
		{
		  int* TmpSzValues = SpectrumFile.GetAsIntegerArray(0);
		  int* TmpKValues = SpectrumFile.GetAsIntegerArray(1);
		  int* TmpInvValues = SpectrumFile.GetAsIntegerArray(2);
		  double* TmpEnergies = SpectrumFile.GetAsDoubleArray(3);
		  int TmpIndex = 0; 
		  while ((TmpIndex < TotalNbrEnergies) && 
			 ((TmpSzValues[TmpIndex] != TotalSz) || (TmpKValues[TmpIndex] != XMomentum) || 
			  (TmpInvValues[TmpIndex] != InversionSector)))
		    ++TmpIndex;
		  if (TmpIndex == TotalNbrEnergies)
		    {
		      cout << "error, the spectrum has no eigenvalues corresponding to the required quantum numbers" << endl;
		      return 0;
		    }
		  int TmpIndex2 = 0;
		  while ((TmpIndex < TotalNbrEnergies) && (TmpIndex2 < NbrStates) && 
			 ((TmpSzValues[TmpIndex] == TotalSz) && (TmpKValues[TmpIndex] == XMomentum) && 
			  (TmpInvValues[TmpIndex] == InversionSector)))
		    {
		      Spectrum[TmpIndex2] = TmpEnergies[TmpIndex];
		      ++TmpIndex2;
		      ++TmpIndex;
		    }
		  if (TmpIndex2 < NbrStates)
		    {
		      cout << "error, the spectrum has less eigenvalues corresponding to the required quantum numbers (" 
			   << TmpIndex2 << " vs " << NbrStates << ")" << endl;
		      return 0;
		    }
		}
	    }
	  else
	    {
	      if (SzSymmetryFlag == true)
		{
		  int* TmpSzValues = SpectrumFile.GetAsIntegerArray(0);
		  int* TmpKValues = SpectrumFile.GetAsIntegerArray(1);
		  int* TmpSzSymValues = SpectrumFile.GetAsIntegerArray(2);
		  double* TmpEnergies = SpectrumFile.GetAsDoubleArray(4);
		  int TmpIndex = 0; 
		  while ((TmpIndex < TotalNbrEnergies) && 
			 ((TmpSzValues[TmpIndex] != TotalSz) || (TmpKValues[TmpIndex] != XMomentum) || 
			  (TmpSzSymValues[TmpIndex] != SzSymmetrySector)))
		    ++TmpIndex;
		  if (TmpIndex == TotalNbrEnergies)
		    {
		      cout << "error, the spectrum has no eigenvalues corresponding to the required quantum numbers" << endl;
		      return 0;
		    }
		  int TmpIndex2 = 0;
		  while ((TmpIndex < TotalNbrEnergies) && (TmpIndex2 < NbrStates) && 
			 ((TmpSzValues[TmpIndex] == TotalSz) && (TmpKValues[TmpIndex] == XMomentum) && 
			  (TmpSzSymValues[TmpIndex] == SzSymmetrySector)))
		    {
		      Spectrum[TmpIndex2] = TmpEnergies[TmpIndex];
		      ++TmpIndex2;
		      ++TmpIndex;
		    }
		  if (TmpIndex2 < NbrStates)
		    {
		      cout << "error, the spectrum has less eigenvalues corresponding to the required quantum numbers (" 
			   << TmpIndex2 << " vs " << NbrStates << ")" << endl;
		      return 0;
		    }
		}
	      else
		{
		  int* TmpSzValues = SpectrumFile.GetAsIntegerArray(0);
		  int* TmpKValues = SpectrumFile.GetAsIntegerArray(1);
		  double* TmpEnergies = SpectrumFile.GetAsDoubleArray(2);
		  int TmpIndex = 0; 
		  while ((TmpIndex < TotalNbrEnergies) && 
			 ((TmpSzValues[TmpIndex] != TotalSz) || (TmpKValues[TmpIndex] != XMomentum)))
		    ++TmpIndex;
		  if (TmpIndex == TotalNbrEnergies)
		    {
		      cout << "error, the spectrum has no eigenvalues corresponding to the required quantum numbers" << endl;
		      return 0;
		    }
		  int TmpIndex2 = 0;
		  while ((TmpIndex < TotalNbrEnergies) && (TmpIndex2 < NbrStates) && 
			 ((TmpSzValues[TmpIndex] == TotalSz) && (TmpKValues[TmpIndex] == XMomentum)))
		    {
		      Spectrum[TmpIndex2] = TmpEnergies[TmpIndex];
		      ++TmpIndex2;
		      ++TmpIndex;
		    }
		  if (TmpIndex2 < NbrStates)
		    {
		      cout << "error, the spectrum has less eigenvalues corresponding to the required quantum numbers (" 
			   << TmpIndex2 << " vs " << NbrStates << ")" << endl;
		      return 0;
		    }
		}
	    }
	}
      NbrDegeneracyIndices = 0;
      Degeneracies = new int [NbrStates];
//       for (int i = 0; i < NbrStates; ++i)
// 	{
// 	  Degeneracies[i] = 1;
// 	  ++NbrDegeneracyIndices;
// 	}
      int TmpIndex = 1;
      while (TmpIndex < NbrStates)
	{
	  int TmpDegeneracy = 1;
	  while ((TmpIndex < NbrStates) && ((fabs(Spectrum[TmpIndex] - Spectrum[TmpIndex - 1]) < Error) || 
					    (fabs(Spectrum[TmpIndex] - Spectrum[TmpIndex - 1]) < (Error * fabs(Spectrum[TmpIndex])))))
	    {
	      ++TmpIndex;
	      ++TmpDegeneracy;
	    }
	  Degeneracies[NbrDegeneracyIndices] = TmpDegeneracy;
	  ++NbrDegeneracyIndices;
	  ++TmpIndex;	  
	}
      int TmpSum = 0;
      for (int i = 0; i < NbrDegeneracyIndices; ++i)
	{
	  TmpSum += Degeneracies[i];
	}
      if (TmpSum != NbrStates)
	{
	  Degeneracies[NbrDegeneracyIndices] = 1;
	  ++NbrDegeneracyIndices;
	}
    }
  
  double* S2Values = new double[NbrStates];
  
  if (Momentum1DFlag == false)
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
// 	  Complex TmpS2 = TmpOperator.MatrixElement(TmpState, TmpState);
// 	  cout << "S^2 = " << TmpS2.Re << endl; 
// 	  cout << "2S = " << (sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0) << " " << round(sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0) << endl; 
	}
      else
	{
	  // 	  Complex TmpS2 = TmpOperator.MatrixElement(TmpState, TmpState);
// 	  cout << "<S^2>=" << TmpS2.Re << " <S>=" << (0.5 * (sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0)) << endl;
// 	  cout << "round(<2S>)=" << round(sqrt((4.0 * TmpS2.Re) + 1.0) - 1.0) << endl; 
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
	      {
		if (InversionFlag == true)
		  {
		    if (SzSymmetryFlag == true)
		      {
			Space = new Spin1_2ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, XMomentum, 1, InversionSector, SzSymmetrySector, TotalSz, 1000000, 1000000);
		      }
		    else
		      {
			Space = new Spin1_2ChainWithTranslationsAndInversionSymmetry (NbrSpins, XMomentum, 1, InversionSector, TotalSz, 1000000, 1000000);
		      }
		  }
		else
		  {
		    if (SzSymmetryFlag == true)
		      {
			Space = new Spin1_2ChainWithTranslationsAndSzSymmetry (NbrSpins, XMomentum, 1, SzSymmetrySector, TotalSz, 1000000, 1000000);
		      }
		    else
		      {
			Space = new Spin1_2ChainWithTranslations (NbrSpins, XMomentum, 1, TotalSz, 1000000, 1000000);
		      }
		  }
	      }
	      break;
	    case 2 :
	      {
		if (InversionFlag == true)
		  {
		    if (SzSymmetryFlag == true)
		      {
			Space = new Spin1ChainWithTranslationsAndSzInversionSymmetries (NbrSpins, XMomentum, InversionSector, SzSymmetrySector, TotalSz);
		      }
		    else
		      {
			Space = new Spin1ChainWithTranslationsAndInversionSymmetry (NbrSpins, XMomentum, InversionSector, TotalSz);
		      }
		  }
		else
		  {
		    if (SzSymmetryFlag == true)
		      {
			Space = new Spin1ChainWithTranslationsAndSzSymmetry (NbrSpins, XMomentum, SzSymmetrySector, TotalSz);
		      }
		    else
		      {
			Space = new Spin1ChainWithTranslations (NbrSpins, XMomentum, TotalSz);
		      }
		  }
	      }
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
      
      if (RealEigenstates.GetNbrColumn() != 0)
	{
	  if (RealEigenstates.GetNbrRow() != Space->GetHilbertSpaceDimension())
	    {
	      cout << "dimension mismatch between the eigenstates (" << RealEigenstates.GetNbrRow() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
	      return 0;
	    }
	}
      else
	{
	  if (ComplexEigenstates.GetNbrRow() != Space->GetHilbertSpaceDimension())
	    {
	      cout << "dimension mismatch between the eigenstates (" << ComplexEigenstates.GetNbrRow() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
	      return 0;
	    }
	}
      
      SpinWith1DTranslationS2Operator TmpOperator(Space, NbrSpins);
      int NbrNonDegenerateStates = 0;
      for (int k = 0; k < NbrDegeneracyIndices; ++k)
	{
	  if (Degeneracies[k] == 1)
	    {
	      ++NbrNonDegenerateStates;
	    }
	}
      RealVector* TmpNonDegenerateRealEigenstates = new RealVector[NbrNonDegenerateStates];
      ComplexVector* TmpNonDegenerateComplexEigenstates = new ComplexVector[NbrNonDegenerateStates];
      Complex* TmpNonDegenerateS2Values = new Complex[NbrNonDegenerateStates];
      int TmpIndex = 0;
      NbrNonDegenerateStates = 0;
      if (Manager.GetBoolean("complex") == false)
	{
	  for (int k = 0; k < NbrDegeneracyIndices; ++k)
	    {
	      if (Degeneracies[k] == 1)
		{
		  TmpNonDegenerateRealEigenstates[NbrNonDegenerateStates] = RealEigenstates[TmpIndex];
		  ++NbrNonDegenerateStates;
		}
	      TmpIndex += Degeneracies[k];	  
	    }
	}
      else
	{
	  for (int k = 0; k < NbrDegeneracyIndices; ++k)
	    {
	      if (Degeneracies[k] == 1)
		{
		  TmpNonDegenerateComplexEigenstates[NbrNonDegenerateStates] = RealEigenstates[TmpIndex];
		  ++NbrNonDegenerateStates;
		}
	      TmpIndex += Degeneracies[k];	  
	    }
	}
      if (Manager.GetBoolean("complex") == false)
	{
	  OperatorMultipleMatrixElementOperation TmpOperation (&TmpOperator, TmpNonDegenerateRealEigenstates, TmpNonDegenerateRealEigenstates, NbrNonDegenerateStates);
	  TmpOperation.ApplyOperation(Architecture.GetArchitecture());
	  for (int k = 0; k < NbrNonDegenerateStates; ++k)
	    {
	      TmpNonDegenerateS2Values[k] = TmpOperation.GetMatrixElement(k);
	    }
	}
      else
	{
	  OperatorMultipleMatrixElementOperation TmpOperation (&TmpOperator, TmpNonDegenerateComplexEigenstates, TmpNonDegenerateComplexEigenstates, NbrNonDegenerateStates);
	  TmpOperation.ApplyOperation(Architecture.GetArchitecture());
	  for (int k = 0; k < NbrNonDegenerateStates; ++k)
	    {
	      TmpNonDegenerateS2Values[k] = TmpOperation.GetMatrixElement(k);
	    }
	}

      TmpIndex = 0;
      NbrNonDegenerateStates = 0;
      for (int k = 0; k < NbrDegeneracyIndices; ++k)
	{
	  if (Degeneracies[k] == 1)
	    {
	      S2Values[TmpIndex] = TmpNonDegenerateS2Values[NbrNonDegenerateStates].Re;
	      ++NbrNonDegenerateStates;
	    }
	  else
	    {
	      cout << "degeneracy at " << TmpIndex << "(" << Degeneracies[k] << " states)" << endl;
	      HermitianMatrix S2Matrix(Degeneracies[k], true);
	      RealVector TmpRealVector;
	      ComplexVector TmpComplexVector;
	      if (Manager.GetBoolean("complex") == false)
		{
		  RealVector TmpVector(Space->GetHilbertSpaceDimension());
		  for (int i = 0; i < Degeneracies[k]; ++i)
		    {		 
		      TmpOperator.Multiply(RealEigenstates[TmpIndex + i], TmpVector);
		      for (int j = i; j < Degeneracies[k]; ++j)
			{
			  Complex TmpS2 = RealEigenstates[TmpIndex + j] * TmpVector;
			  S2Matrix.SetMatrixElement(j, i, TmpS2);
			}
		    }
		}
	      else
		{
		  ComplexVector TmpVector(Space->GetHilbertSpaceDimension());
		  for (int i = 0; i < Degeneracies[k]; ++i)
		    {		      
		      TmpOperator.Multiply(ComplexEigenstates[TmpIndex + i], TmpVector);
		      for (int j = i; j < Degeneracies[k]; ++j)
			{
			  Complex TmpS2= ComplexEigenstates[TmpIndex + j] * TmpVector;
			  S2Matrix.SetMatrixElement(j, i, TmpS2);
			}
		    }
		}	    
	      RealDiagonalMatrix TmpS2Eigenvalues(Degeneracies[k]);
	      S2Matrix.LapackDiagonalize(TmpS2Eigenvalues);
	      for (int i = 0; i < Degeneracies[k]; ++i)
		{
		  S2Values[TmpIndex + i] = TmpS2Eigenvalues[i];
		}
	    }
	  TmpIndex += Degeneracies[k];
	}
    }
  if (Spectrum == 0)
    {
      for (int i = 0; i < NbrStates; ++i)
	{
	  double TmpS2 = S2Values[i];
	  cout << i << " : <S^2>=" << TmpS2 << " <S>=" << (0.5 * (sqrt((4.0 * TmpS2) + 1.0) - 1.0)) <<  " round(<2S>)=" <<  round(sqrt((4.0 * TmpS2) + 1.0) - 1.0) << endl; 
	}
    }
  else
    {
      char* OutputFileName = new char[strlen(Manager.GetString("multiple-states")) + 32];
      sprintf (OutputFileName, "%s.s", Manager.GetString("multiple-states"));
      ofstream File;
      File.open(OutputFileName, ios::out);
      File.precision(14);
      char* LinePrefix = new char[512];
      if (Momentum1DFlag == false)
	{
	  sprintf (LinePrefix, "%d ", TotalSz);
	}
      else
	{
	  if (InversionFlag == true)
	    {
	      if (SzSymmetryFlag == true)
		{
		  sprintf (LinePrefix, "%d %d %d %d ", TotalSz, XMomentum, SzSymmetrySector, InversionSector);
		}
	      else
		{
		  sprintf (LinePrefix, "%d %d %d ", TotalSz, XMomentum, InversionSector);
		}
	    }
	  else
	    {
	      if (SzSymmetryFlag == true)
		{
		  sprintf (LinePrefix, "%d %d %d ", TotalSz, XMomentum, SzSymmetrySector);
		}
	      else
		{
		  sprintf (LinePrefix, "%d %d ", TotalSz, XMomentum);
		}
	    }
	}
      for (int i = 0; i < NbrStates; ++i)
	{
	  double TmpS2 = S2Values[i];
	  File << LinePrefix << Spectrum[i] << " " << TmpS2 << " " << (0.5 * (sqrt((4.0 * TmpS2) + 1.0) - 1.0)) <<  " " <<  round(sqrt((4.0 * TmpS2) + 1.0) - 1.0) << endl; 
	}
      File.close();
    }
  return 0;
}

