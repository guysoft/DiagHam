////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of bosons on disk using the Haldane basis             //
//                            for system size such that                       //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 27/01/2011                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "config.h"
#include "HilbertSpace/BosonOnDiskHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/ArrayTools.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h" 
#include "GeneralTools/Endian.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHESphereJackGeneratorOperation.h"

#include <math.h>
#include <fstream>
#include <sys/time.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::fstream;

using std::ios;


// default constructor
//

BosonOnDiskHaldaneHugeBasisShort::BosonOnDiskHaldaneHugeBasisShort ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson
// maxFileSize = maximum file size (in MBytes)
// referenceState = array that describes the reference state to start from
// memory = amount of memory granted for precalculations
// referenceState = array that describes the reference state to start from
// symmetricFlag = indicate if a symmetric basis has to be used (only available if totalLz = 0)
// fullDimension = provide the full (i.e. without squeezing) Hilbert space dimension (0 if it has to be computed)

BosonOnDiskHaldaneHugeBasisShort::BosonOnDiskHaldaneHugeBasisShort (int nbrBosons, int totalLz, int lzMax, unsigned long maxFileSize, int* referenceState, unsigned long memory, bool symmetricFlag, long fullDimension)
{
  this->TargetSpace = this;
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  int ShiftedLzMax = this->LzMax + nbrBosons - 1;
  int* TmpReferenceState = new int[ShiftedLzMax + 1];
  int TmpIndex = 0;
  for (int i = 0; i <= ShiftedLzMax; ++i)
    TmpReferenceState[i] = 0;   
  for (int i = 0; i <= this->LzMax; ++i)
    {
      for (int j = 0; j < referenceState[i]; ++j)
	{
	  TmpReferenceState[TmpIndex] = 1;
	  ++TmpIndex;
	}
      ++TmpIndex;
    }
  this->FermionBasis = 0;
  this->FermionHugeBasis = new FermionOnSphereHaldaneHugeBasis(nbrBosons, totalLz, this->LzMax + nbrBosons - 1, maxFileSize, TmpReferenceState, memory);
  delete[] TmpReferenceState;
  this->HilbertSpaceDimension = this->FermionHugeBasis->GetHilbertSpaceDimension();
  this->LargeHilbertSpaceDimension = this->FermionHugeBasis->GetLargeHilbertSpaceDimension();

  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->TemporaryMonomial = new unsigned long [this->NbrBosons];
  this->TemporaryMonomial2 = new unsigned long [this->NbrBosons];

  this->Flag.Initialize();
  int TmpLzMax = this->LzMax;


  if (this->ShiftedTotalLz < TmpLzMax)
    {
      TmpLzMax = this->ShiftedTotalLz;	  
    }
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;

}

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memoryHilbert = amount of memory granted to store the Hilbert space (in Mbytes)

BosonOnDiskHaldaneHugeBasisShort::BosonOnDiskHaldaneHugeBasisShort (char* fileName, unsigned long memoryHilbert)
{
  timeval TotalStartingTime;
  gettimeofday (&(TotalStartingTime), 0);

  this->FermionBasis = 0;
  this->FermionHugeBasis = new FermionOnSphereHaldaneHugeBasis(fileName, memoryHilbert);
  this->HilbertSpaceDimension = this->FermionHugeBasis->GetHilbertSpaceDimension();
  this->LargeHilbertSpaceDimension = this->FermionHugeBasis->GetLargeHilbertSpaceDimension();
  
  timeval TotalEndingTime;
  gettimeofday (&(TotalEndingTime), 0);
  double Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
  cout << "Hilbert space factorization done in " << Dt << "s" <<endl;
 
//  gettimeofday (&(TotalStartingTime), 0);

//   ifstream FileHilbert;
//   FileHilbert.open(this->FermionHugeBasis->HilbertSpaceFileName, ios::binary | ios::in);
//   FileHilbert.seekg (this->FermionHugeBasis->FileHeaderSize, ios::beg);
//   unsigned long CurrentPartition;
//   for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
//     {
//       ReadLittleEndian(FileHilbert, CurrentPartition);
//       if (i != this->FermionHugeBasis->FindStateIndexFactorized(CurrentPartition))
// //      if (CurrentPartition != this->FermionHugeBasis->GetStateFactorized(i))
// 	cout << i << " " << this->FermionHugeBasis->FindStateIndexFactorized(CurrentPartition) << " " << hex << CurrentPartition << " " << this->FermionHugeBasis->GetStateFactorized(i) << dec << endl;
//     }
//   FileHilbert.close();

//   long Test = 5;
//   long TmpShift = 470;
//   long TmpStep = 5;
//   unsigned long* TmpStateArray = new unsigned long [Test];
//   long* TmpIndexArray = new long[Test];
//   for (long i = 0; i < Test; ++i)
//     TmpStateArray[i] = this->FermionHugeBasis->GetStateFactorized((i * TmpStep) + TmpShift);
//   this->FermionHugeBasis->FindMultipleStateIndexFactorized(TmpStateArray, Test, TmpIndexArray);
//   for (long i = 0; i < Test; ++i)
//     if (((i * TmpStep) + TmpShift) != TmpIndexArray[i])
//       cout << ((i * TmpStep) + TmpShift) << " : " << TmpIndexArray[i] << " " << endl;
 
//   gettimeofday (&(TotalEndingTime), 0);
//   Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
//     ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
//   cout << "Hilbert space consistency check done in " << Dt << "s" <<endl;
//  exit(0) ;

  this->NbrBosons = this->FermionHugeBasis->NbrFermions;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = this->FermionHugeBasis->TotalLz;
  this->LzMax = this->FermionHugeBasis->LzMax - this->NbrBosons + 1;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;

  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->TemporaryMonomial = new unsigned long [this->NbrBosons];
  this->TemporaryMonomial2 = new unsigned long [this->NbrBosons];
  this->Flag.Initialize();
  int TmpLzMax = this->LzMax;

  if (this->ShiftedTotalLz < TmpLzMax)
    {
      TmpLzMax = this->ShiftedTotalLz;	  
    }
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;
}


// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnDiskHaldaneHugeBasisShort::BosonOnDiskHaldaneHugeBasisShort(const BosonOnDiskHaldaneHugeBasisShort& bosons)
{
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->FermionBasis = 0;
  this->FermionHugeBasis = (FermionOnSphereHaldaneHugeBasis*) bosons.FermionHugeBasis->Clone();
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];  
  this->TemporaryMonomial = new unsigned long [this->NbrBosons];
  this->TemporaryMonomial2 = new unsigned long [this->NbrBosons];
}

// destructor
//

BosonOnDiskHaldaneHugeBasisShort::~BosonOnDiskHaldaneHugeBasisShort ()
{
  if (this->FermionHugeBasis != 0)
    delete this->FermionHugeBasis;
  delete[] this->TemporaryMonomial;
  delete[] this->TemporaryMonomial2;
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnDiskHaldaneHugeBasisShort& BosonOnDiskHaldaneHugeBasisShort::operator = (const BosonOnDiskHaldaneHugeBasisShort& bosons)
{
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
    }
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->FermionBasis = 0;
  this->FermionHugeBasis = (FermionOnSphereHaldaneHugeBasis*) bosons.FermionHugeBasis->Clone();
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->TemporaryMonomial = new unsigned long [this->NbrBosons];
  this->TemporaryMonomial2 = new unsigned long [this->NbrBosons];

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnDiskHaldaneHugeBasisShort::Clone()
{
  return new BosonOnDiskHaldaneHugeBasisShort(*this);
}


// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& BosonOnSphereHaldaneHugeBasisShort::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  double Factor = 1.0 / state[reference];
  state[reference] = 1.0;
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  unsigned long TmpState = 0x0ul;
  if (this->FermionHugeBasis->CheckDiskStorage() == false)
    {
      TmpState = this->FermionHugeBasis->StateDescription[reference];
    }
  else
    {
      TmpState = this->FermionHugeBasis->GetStateFactorized(reference);
    }
   while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
  this->ConvertToMonomial(TmpState, TmpLzMax, this->TemporaryMonomial2);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
    }
  FactorialCoefficient ReferenceFactorial;
  FactorialCoefficient Factorial;
  this->FermionToBoson(TmpState, TmpLzMax, 
		       this->TemporaryState, this->TemporaryStateLzMax);
  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
    if (this->TemporaryState[k] > 1)
      ReferenceFactorial.FactorialDivide(this->TemporaryState[k]);
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    if (i != reference)
      {
	if (this->FermionHugeBasis->CheckDiskStorage() == false)
	  {
	    TmpState = this->FermionHugeBasis->StateDescription[i];
	  }
	else
	  {
	    TmpState = this->FermionHugeBasis->GetStateFactorized(i);
	  }
	while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	  --TmpLzMax;
	this->ConvertToMonomial(TmpState, TmpLzMax, this->TemporaryMonomial);
	int Index1 = 0;
	int Index2 = 0;
	double Coefficient = Factor;
	while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	  {
	    while ((Index1 < this->NbrBosons) && (this->TemporaryMonomial2[Index1] > this->TemporaryMonomial[Index2]))
	      {
		Coefficient *= InvSqrtCoefficients[this->TemporaryMonomial2[Index1]];
		++Index1;
	      }
	    while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (this->TemporaryMonomial2[Index1] == this->TemporaryMonomial[Index2]))
	      {
		++Index1;
		++Index2;
	      }
	    while ((Index2 < this->NbrBosons) && (this->TemporaryMonomial2[Index1] < this->TemporaryMonomial[Index2]))
	      {
		Coefficient *= SqrtCoefficients[this->TemporaryMonomial[Index2]];
		++Index2;
	      }	  
	  }
	while (Index1 < this->NbrBosons)
	  {
	    Coefficient *= InvSqrtCoefficients[this->TemporaryMonomial2[Index1]];
	    ++Index1;
	  }
	while (Index2 < this->NbrBosons)
	  {
	    Coefficient *= SqrtCoefficients[this->TemporaryMonomial2[Index2]];
	    ++Index2;
	  }
	if (symmetryFactor == true)
	  {
	    Factorial = ReferenceFactorial;
	    this->FermionToBoson(TmpState, TmpLzMax, 
				 this->TemporaryState, this->TemporaryStateLzMax);
	    for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	      if (this->TemporaryState[k] > 1)
		Factorial.FactorialMultiply(this->TemporaryState[k]);
	    Coefficient *= sqrt(Factorial.GetNumericalValue());
	  }
	state[i] *= Coefficient;
      if ((i & 0x3fffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
      }
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// symmetryFactor = if true also add the symmetry factors
// return value = converted state

RealVector& BosonOnSphereHaldaneHugeBasisShort::ConvertFromUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  double Factor = 1.0;
  int TmpLzMax = this->FermionHugeBasis->LzMax;
  unsigned long TmpState = 0x0ul;
  if (this->FermionHugeBasis->CheckDiskStorage() == false)
    {
      TmpState = this->FermionHugeBasis->StateDescription[reference];
    }
  else
    {
      TmpState = this->FermionHugeBasis->GetStateFactorized(reference);
    }
   while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
    --TmpLzMax;
   this->ConvertToMonomial(TmpState, TmpLzMax, this->TemporaryMonomial2);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      SqrtCoefficients[k] = 1.0 / InvSqrtCoefficients[k];
    }
  FactorialCoefficient ReferenceFactorial;
  FactorialCoefficient Factorial;
  this->FermionToBoson(TmpState, TmpLzMax, 
		       this->TemporaryState, this->TemporaryStateLzMax);
  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
    if (this->TemporaryState[k] > 1)
      ReferenceFactorial.FactorialMultiply(this->TemporaryState[k]);
  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
    if (i != reference)
      {
	if (this->FermionHugeBasis->CheckDiskStorage() == false)
	  {
	    TmpState = this->FermionHugeBasis->StateDescription[i];
	  }
	else
	  {
	    TmpState = this->FermionHugeBasis->GetStateFactorized(i);
	  }
	while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	  --TmpLzMax;
	this->ConvertToMonomial(TmpState, TmpLzMax, this->TemporaryMonomial);
	int Index1 = 0;
	int Index2 = 0;
	double Coefficient = Factor;
	if (symmetryFactor == true)
	  {
	    while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	      {
		while ((Index1 < this->NbrBosons) && (this->TemporaryMonomial2[Index1] > this->TemporaryMonomial[Index2]))
		  {
		    Coefficient *= InvSqrtCoefficients[this->TemporaryMonomial2[Index1]];
		    ++Index1;
		  }
		while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (this->TemporaryMonomial2[Index1] == this->TemporaryMonomial[Index2]))
		  {
		    ++Index1;
		    ++Index2;
		  }
		while ((Index2 < this->NbrBosons) && (this->TemporaryMonomial2[Index1] < this->TemporaryMonomial[Index2]))
		  {
		    Coefficient *= SqrtCoefficients[this->TemporaryMonomial[Index2]];
		    ++Index2;
		  }	  
	      }
	    while (Index1 < this->NbrBosons)
	      {
		Coefficient *= InvSqrtCoefficients[this->TemporaryMonomial2[Index1]];
		++Index1;
	      }
	    while (Index2 < this->NbrBosons)
	      {
		Coefficient *= SqrtCoefficients[this->TemporaryMonomial2[Index2]];
		++Index2;
	      }
	    Factorial = ReferenceFactorial;
	    this->FermionToBoson(TmpState, TmpLzMax, 
				 this->TemporaryState, this->TemporaryStateLzMax);
	    for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	      if (this->TemporaryState[k] > 1)
		Factorial.FactorialDivide(this->TemporaryState[k]);
	    Coefficient *= sqrt(Factorial.GetNumericalValue());
	  }
	else
	  {
	    Factorial = ReferenceFactorial;
	    this->FermionToBoson(TmpState, TmpLzMax, 
				 this->TemporaryState, this->TemporaryStateLzMax);
	    for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	      if (this->TemporaryState[k] > 1)
		Factorial.FactorialDivide(this->TemporaryState[k]);
	    Coefficient *= sqrt(Factorial.GetNumericalValue());
	  }
	state[i] *= Coefficient;
	if ((i & 0x3fffl) == 0l)
	  {
	    cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	    cout.flush();
	  }
      }
  state /= state.Norm();
  return state;
}

