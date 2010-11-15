////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of bosons on sphere  using the Haldane basis           //
//                            for system size such that                       //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 08/07/2008                      //
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
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/ArrayTools.h"
#include "Polynomial/RationalPolynomial.h"

#include <math.h>


using std::cout;
using std::endl;


// default constructor
//

BosonOnSphereHaldaneBasisShort::BosonOnSphereHaldaneBasisShort ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson
// referenceState = array that describes the reference state to start from

BosonOnSphereHaldaneBasisShort::BosonOnSphereHaldaneBasisShort (int nbrBosons, int& totalLz, int lzMax, int* referenceState)
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
  this->TotalLz = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      for (int j = 0; j < referenceState[i]; ++j)
	{
	  TmpReferenceState[TmpIndex] = 1;
	  ++TmpIndex;
	  this->TotalLz += i;
	}
      ++TmpIndex;
    }
  this->TotalLz = (2 * this->TotalLz) - (this->NbrBosons *  this->LzMax);
  this->FermionBasis = new FermionOnSphereHaldaneBasis(nbrBosons, totalLz, this->LzMax + nbrBosons - 1, TmpReferenceState);
  totalLz = this->TotalLz;
  delete[] TmpReferenceState;
  this->HilbertSpaceDimension = this->FermionBasis->GetHilbertSpaceDimension();
  this->LargeHilbertSpaceDimension = this->FermionBasis->GetLargeHilbertSpaceDimension();

  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
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
// memory = amount of memory granted for precalculations

BosonOnSphereHaldaneBasisShort::BosonOnSphereHaldaneBasisShort (char* fileName)
{
  this->FermionBasis = new FermionOnSphereHaldaneBasis(fileName);
  this->HilbertSpaceDimension = this->FermionBasis->GetHilbertSpaceDimension();
  this->LargeHilbertSpaceDimension = this->FermionBasis->GetLargeHilbertSpaceDimension();

  this->NbrBosons = this->FermionBasis->NbrFermions;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = this->FermionBasis->TotalLz;
  this->LzMax = this->FermionBasis->LzMax - this->NbrBosons + 1;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;

  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
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

BosonOnSphereHaldaneBasisShort::BosonOnSphereHaldaneBasisShort(const BosonOnSphereHaldaneBasisShort& bosons)
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
  this->FermionBasis = (FermionOnSphereHaldaneBasis*) bosons.FermionBasis->Clone();
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];  
}

// destructor
//

BosonOnSphereHaldaneBasisShort::~BosonOnSphereHaldaneBasisShort ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereHaldaneBasisShort& BosonOnSphereHaldaneBasisShort::operator = (const BosonOnSphereHaldaneBasisShort& bosons)
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
  this->FermionBasis = (FermionOnSphereHaldaneBasis*) bosons.FermionBasis->Clone();
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];

  return *this;
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool BosonOnSphereHaldaneBasisShort::WriteHilbertSpace (char* fileName)
{
  return this->FermionBasis->WriteHilbertSpace(fileName);
}

// convert a given state from Haldane basis to the usual n-body basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector BosonOnSphereHaldaneBasisShort::ConvertToNbodyBasis(RealVector& state, BosonOnSphereShort& nbodyBasis)
{
  RealVector TmpVector (nbodyBasis.GetHilbertSpaceDimension(), true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      TmpVector[nbodyBasis.FermionBasis->FindStateIndex(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i])] = state[i];
    }
  return TmpVector;
}

// convert a given state from the usual n-body basis to the Haldane basis
//
// state = reference on the vector to convert
// nbodyBasis = reference on the nbody-basis to use
// return value = converted vector

RealVector BosonOnSphereHaldaneBasisShort::ConvertFromNbodyBasis(RealVector& state, BosonOnSphereShort& nbodyBasis)
{
  RealVector TmpVector (this->HilbertSpaceDimension, true);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    TmpVector[i] = state[nbodyBasis.FermionBasis->FindStateIndex(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i])];
  TmpVector /= TmpVector.Norm();
  return TmpVector;
}

// create the Jack polynomial decomposition corresponding to the root partition
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& BosonOnSphereHaldaneBasisShort::GenerateJackPolynomial(RealVector& jack, double alpha)
{
  jack[0] = 1.0;
  double InvAlpha =  2.0 / alpha;

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));
  int ReducedNbrBosons = this->NbrBosons - 1;

  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (jack[i] == 0.0)
	{
	  double Rho = 0.0;
	  unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
	  this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
	  for (int j = 0; j < this->NbrBosons; ++j)
	    Rho += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));
	  double Coefficient = 0.0;
	  for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
	    for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
	      {
		double Diff = (double) (TmpMonomial[j1] - TmpMonomial[j2]);
		unsigned int Max = TmpMonomial[j2];
		unsigned long TmpState = 0x0ul;
		int Tmpj1 = j1;
		int Tmpj2 = j2;
		for (int l = 0; l < this->NbrBosons; ++l)
		  TmpMonomial2[l] = TmpMonomial[l];	    
		for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
		  {
		    ++TmpMonomial2[Tmpj1];
		    --TmpMonomial2[Tmpj2];
		    Diff += 2.0;
		    while ((Tmpj1 > 0) && (TmpMonomial2[Tmpj1] > TmpMonomial2[Tmpj1 - 1]))
		      {
			unsigned long Tmp = TmpMonomial2[Tmpj1 - 1];
			TmpMonomial2[Tmpj1 - 1] = TmpMonomial2[Tmpj1];
			TmpMonomial2[Tmpj1] = Tmp;
			--Tmpj1;
		      }
		    while ((Tmpj2 < ReducedNbrBosons) && (TmpMonomial2[Tmpj2] < TmpMonomial2[Tmpj2 + 1]))
		      {
			unsigned long Tmp = TmpMonomial2[Tmpj2 + 1];
			TmpMonomial2[Tmpj2 + 1] = TmpMonomial2[Tmpj2];
			TmpMonomial2[Tmpj2] = Tmp;
			++Tmpj2;
		      }
		    TmpState = this->ConvertFromMonomial(TmpMonomial2);
		    if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
		      {
			long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
			if (TmpIndex < this->HilbertSpaceDimension)
			  Coefficient += Diff * jack[TmpIndex];
		      }
		  }
	      }
	  jack[i] = Coefficient * InvAlpha / (RhoRoot - Rho);
	}
      else
	{
	  cout << "toto" << endl;
	}
      if ((i & 0xffl) == 0l)
	{
	  //	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  //	  cout.flush();
	}
    }
  delete[] TmpMonomial;
  cout << endl;

  return jack;
}

// create the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alphaNumerator = numerator of the Jack polynomial alpha coefficient
// alphaDenominator = numerator of the Jack polynomial alpha coefficient
// symbolicDepth = use symbolic calculation to solve singular values if non zero, if greater than zero it will use symbolic calculation up to a given depth (below that depth it will rely on numerical calculation),
//                 -1 if the symbolic calculation has to be done up to the point the singular value problem has been solved
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RationalVector& BosonOnSphereHaldaneBasisShort::GenerateJackPolynomial(RationalVector& jack, long alphaNumerator, long alphaDenominator, int symbolicDepth)
{
  jack[0] = 1l;
  Rational InvAlpha (2l * alphaDenominator, alphaNumerator);

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];

  Rational RhoRoot = 0l;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      RhoRoot += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
    }
  int ReducedNbrBosons = this->NbrBosons - 1;

  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (jack[i].Num() == 0l)
	{
	  Rational Rho = 0l;
	  unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
	  this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
	  for (int j = 0; j < this->NbrBosons; ++j)
	    Rho += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
	  if (Rho == RhoRoot)
	    {
	      if (symbolicDepth == 0)
		{
		  cout << "warning : singular value detected at position " << i << ", skipping the rest of the calculation" << endl;
		  return jack;
		}
	      int Depth = 1;
	      bool FullSymbolicFlag = false;
	      bool SolvedFlag = false;
	      while (((Depth <= symbolicDepth) || (symbolicDepth < 0)) && (SolvedFlag == false) && (FullSymbolicFlag == false))
		{
		  RationalPolynomial TmpNumerator;
		  RationalPolynomial TmpDenominator;
		  FullSymbolicFlag = this->GenerateSingleJackPolynomialCoefficient(jack, i, TmpNumerator, TmpDenominator, Depth);
		  Rational Tmp = TmpNumerator.PolynomialEvaluate(InvAlpha);
		  if (Tmp.Num() == 0l)
		    {
		      TmpNumerator.MonomialDivision(InvAlpha);
		      Tmp = TmpNumerator.PolynomialEvaluate(InvAlpha);
		      Tmp /= TmpDenominator.PolynomialEvaluate(InvAlpha);
		      jack[i] = (Tmp * InvAlpha) / (RhoRoot - Rho);
		      SolvedFlag = true;
		    }
		  ++Depth;
		}
	      if (SolvedFlag == false)
		{
		  if (FullSymbolicFlag == true)
		    cout << "warning : singular value detected at position " << i << ", skipping the rest of the calculation" << endl;
		  return jack;

		}
	    }
	  else
	    {
	      Rational Coefficient = 0;
	      for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
		for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
		  {
		    long Diff = (long) (TmpMonomial[j1] - TmpMonomial[j2]);
		    unsigned int Max = TmpMonomial[j2];
		    unsigned long TmpState = 0x0ul;
		    int Tmpj1 = j1;
		    int Tmpj2 = j2;
		    for (int l = 0; l < this->NbrBosons; ++l)
		      TmpMonomial2[l] = TmpMonomial[l];	    
		    for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
		      {
			++TmpMonomial2[Tmpj1];
			--TmpMonomial2[Tmpj2];
			Diff += 2l;
			while ((Tmpj1 > 0) && (TmpMonomial2[Tmpj1] > TmpMonomial2[Tmpj1 - 1]))
			  {
			    unsigned long Tmp = TmpMonomial2[Tmpj1 - 1];
			    TmpMonomial2[Tmpj1 - 1] = TmpMonomial2[Tmpj1];
			    TmpMonomial2[Tmpj1] = Tmp;
			    --Tmpj1;
			  }
			while ((Tmpj2 < ReducedNbrBosons) && (TmpMonomial2[Tmpj2] < TmpMonomial2[Tmpj2 + 1]))
			  {
			    unsigned long Tmp = TmpMonomial2[Tmpj2 + 1];
			    TmpMonomial2[Tmpj2 + 1] = TmpMonomial2[Tmpj2];
			    TmpMonomial2[Tmpj2] = Tmp;
			    ++Tmpj2;
			  }
			TmpState = this->ConvertFromMonomial(TmpMonomial2);
			if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
			  {
			    long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
			    if (TmpIndex < this->HilbertSpaceDimension)
			      Coefficient += Diff * jack[TmpIndex];
			  }
		      }
		  }
	      jack[i] = (Coefficient * InvAlpha) / (RhoRoot - Rho);
	    }
	}
      if ((i & 0xffffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  delete[] TmpMonomial;
  cout << endl;
  return jack;
}

// compute a single coefficient of the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur and using (partial symbolic calculation)
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// index = index of the component to compute
// numerator = reference on the polynomial where the numerator has to be stored
// denominator = reference on the polynomial where the denominator has to be stored
// depth = depth in the recurrence describing up to which point the symbolic calculation has to be performed 
// return value = true if a fully symbolic calculation has been performed

bool BosonOnSphereHaldaneBasisShort::GenerateSingleJackPolynomialCoefficient(RationalVector& jack, long index, RationalPolynomial& numerator, RationalPolynomial& denominator, int depth)
{
//   Rational InvAlpha (2l * alphaDenominator, alphaNumerator);

//   unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
//   unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];

//   Rational RhoRoot = 0l;
//   unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
//   this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
//   for (int j = 0; j < this->NbrBosons; ++j)
//     {
//       RhoRoot += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
//     }
//   int ReducedNbrBosons = this->NbrBosons - 1;

//   for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
//     {
//       if (jack[i].Num() == 0l)
// 	{
// 	  Rational Rho = 0l;
// 	  unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
// 	  this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
// 	  for (int j = 0; j < this->NbrBosons; ++j)
// 	    Rho += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
// 	  if (Rho == RhoRoot)
// 	    {
// 	      cout << "warning : singular value detected at position " << i << ", skipping the rest of the calculation" << endl;
// 	      return jack;
// 	    }
// 	  Rational Coefficient = 0;
// 	  for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
// 	    for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
// 	      {
// 		long Diff = (long) (TmpMonomial[j1] - TmpMonomial[j2]);
// 		unsigned int Max = TmpMonomial[j2];
// 		unsigned long TmpState = 0x0ul;
// 		int Tmpj1 = j1;
// 		int Tmpj2 = j2;
// 		for (int l = 0; l < this->NbrBosons; ++l)
// 		  TmpMonomial2[l] = TmpMonomial[l];	    
// 		for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
// 		  {
// 		    ++TmpMonomial2[Tmpj1];
// 		    --TmpMonomial2[Tmpj2];
// 		    Diff += 2l;
// 		    while ((Tmpj1 > 0) && (TmpMonomial2[Tmpj1] > TmpMonomial2[Tmpj1 - 1]))
// 		      {
// 			unsigned long Tmp = TmpMonomial2[Tmpj1 - 1];
// 			TmpMonomial2[Tmpj1 - 1] = TmpMonomial2[Tmpj1];
// 			TmpMonomial2[Tmpj1] = Tmp;
// 			--Tmpj1;
// 		      }
// 		    while ((Tmpj2 < ReducedNbrBosons) && (TmpMonomial2[Tmpj2] < TmpMonomial2[Tmpj2 + 1]))
// 		      {
// 			unsigned long Tmp = TmpMonomial2[Tmpj2 + 1];
// 			TmpMonomial2[Tmpj2 + 1] = TmpMonomial2[Tmpj2];
// 			TmpMonomial2[Tmpj2] = Tmp;
// 			++Tmpj2;
// 		      }
// 		    TmpState = this->ConvertFromMonomial(TmpMonomial2);
// 		    if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
// 		      {
// 			long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
// 			if (TmpIndex < this->HilbertSpaceDimension)
// 			  Coefficient += Diff * jack[TmpIndex];
// 		      }
// 		  }
// 	      }
// 	  jack[i] = (Coefficient * InvAlpha) / (RhoRoot - Rho);
// 	}
//       if ((i & 0xffffl) == 0l)
// 	{
// 	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
// 	  cout.flush();
// 	}
//     }
//   delete[] TmpMonomial;
//   cout << endl;
  return false;
}


// create the Jack polynomial decomposition corresponding to the root partition assuming the resulting state is invariant under the Lz<->-Lz symmetry
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// return value = decomposition of the corresponding Jack polynomial on the unnormalized basis

RealVector& BosonOnSphereHaldaneBasisShort::GenerateSymmetrizedJackPolynomial(RealVector& jack, double alpha)
{
  jack[0] = 1.0;
  double InvAlpha =  2.0 / alpha;

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));
  int ReducedNbrBosons = this->NbrBosons - 1;
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    if (jack[i] == 0.0)
      {
	double Rho = 0.0;
	unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
	this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
	for (int j = 0; j < this->NbrBosons; ++j)
	  Rho += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));
	double Coefficient = 0.0;
	for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
	  for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
	    {
	      double Diff = (double) (TmpMonomial[j1] - TmpMonomial[j2]);
	      unsigned int Max = TmpMonomial[j2];
	      unsigned long TmpState = 0x0ul;
	      int Tmpj1 = j1;
	      int Tmpj2 = j2;
	      for (int l = 0; l < this->NbrBosons; ++l)
		TmpMonomial2[l] = TmpMonomial[l];	    
	      for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
		{
		  ++TmpMonomial2[Tmpj1];
		  --TmpMonomial2[Tmpj2];
		  Diff += 2.0;
		  while ((Tmpj1 > 0) && (TmpMonomial2[Tmpj1] > TmpMonomial2[Tmpj1 - 1]))
		    {
		      unsigned long Tmp = TmpMonomial2[Tmpj1 - 1];
		      TmpMonomial2[Tmpj1 - 1] = TmpMonomial2[Tmpj1];
		      TmpMonomial2[Tmpj1] = Tmp;
		      --Tmpj1;
		    }
		  while ((Tmpj2 < ReducedNbrBosons) && (TmpMonomial2[Tmpj2] < TmpMonomial2[Tmpj2 + 1]))
		    {
		      unsigned long Tmp = TmpMonomial2[Tmpj2 + 1];
		      TmpMonomial2[Tmpj2 + 1] = TmpMonomial2[Tmpj2];
		      TmpMonomial2[Tmpj2] = Tmp;
		      ++Tmpj2;
		    }
		  TmpState = this->ConvertFromMonomial(TmpMonomial2);
		  if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
		    {
		      long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
		      if (TmpIndex < this->HilbertSpaceDimension)
			Coefficient += Diff * jack[TmpIndex];
		    }
		}
	    }
	
	long TmpIndex = this->FermionBasis->FindStateIndex(((FermionOnSphereHaldaneBasis*) this->FermionBasis)->GetSymmetricState(CurrentPartition), (this->LzMax - TmpMonomial[ReducedNbrBosons]) + ReducedNbrBosons);
	Coefficient *= InvAlpha;
	Coefficient /= (RhoRoot - Rho);
	if (i < TmpIndex)
	  jack[TmpIndex] = Coefficient;
	jack[i] = Coefficient;
	if ((i & 0xffl) == 0l)
	  {
	    cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100) / this->LargeHilbertSpaceDimension) << "%)           \r";
	    cout.flush();
	  }
      }
  delete[] TmpMonomial;
  cout << endl;

  return jack;
}

// check partitions that may lead to singular coefficient in a given Jack polynomial decomposition
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// error = error when comparing two rho values
// return value = vector with non-zero component being rho factor of possible singular coefficients

RealVector& BosonOnSphereHaldaneBasisShort::CheckPossibleSingularCoefficientsInJackPolynomial(RealVector& jack, double alpha, double error)
{
  double InvAlpha =  2.0 / alpha;

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));
  int ReducedNbrBosons = this->NbrBosons - 1;

  jack[0] = RhoRoot;
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      double Rho = 0.0;
      unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
      for (int j = 0; j < this->NbrBosons; ++j)
	Rho += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));
      if ((fabs(RhoRoot - Rho) < error) || (fabs(RhoRoot - Rho) < (error * fabs(RhoRoot))))
	jack[i] = Rho;
      else
	jack[i] = 0.0;
      }
  delete[] TmpMonomial;

  return jack;
}
  
// check partitions that may lead to singular coefficient in a given Jack polynomial decomposition, assuming only rational numbers occur
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alphaNumerator = numerator of the Jack polynomial alpha coefficient
// alphaDenominator = numerator of the Jack polynomial alpha coefficient
// return value = vector with non-zero component being rho factor of possible singular coefficients

RationalVector& BosonOnSphereHaldaneBasisShort::CheckPossibleSingularCoefficientsInJackPolynomial(RationalVector& jack, long alphaNumerator, long alphaDenominator)
{
  return jack;
}
  
// check partitions that may lead to singular coefficient in a given Jack polynomial decomposition
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// alpha = value of the Jack polynomial alpha coefficient
// error = error when comparing two rho values
// return value = vector with non-zero component being rho factor of possible singular coefficients

void BosonOnSphereHaldaneBasisShort::CheckMaximumConnectedStateInJackPolynomial()
{
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  int ReducedNbrBosons = this->NbrBosons - 1;
  long LargestDistance = 0;
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
      long MinIndex = this->LargeHilbertSpaceDimension;
      for (int j1 = 0; j1 < ReducedNbrBosons; ++j1)
	for (int j2 = j1 + 1; j2 < this->NbrBosons; ++j2)
	  {
	    unsigned int Max = TmpMonomial[j2];
	    unsigned long TmpState = 0x0ul;
	    int Tmpj1 = j1;
	    int Tmpj2 = j2;
	    for (int l = 0; l < this->NbrBosons; ++l)
	      TmpMonomial2[l] = TmpMonomial[l];	    
	    for (unsigned int k = 1; (k <= Max) && (TmpState < MaxRoot); ++k)
	      {
		++TmpMonomial2[Tmpj1];
		--TmpMonomial2[Tmpj2];
		while ((Tmpj1 > 0) && (TmpMonomial2[Tmpj1] > TmpMonomial2[Tmpj1 - 1]))
		  {
		    unsigned long Tmp = TmpMonomial2[Tmpj1 - 1];
		    TmpMonomial2[Tmpj1 - 1] = TmpMonomial2[Tmpj1];
		    TmpMonomial2[Tmpj1] = Tmp;
		    --Tmpj1;
		  }
		while ((Tmpj2 < ReducedNbrBosons) && (TmpMonomial2[Tmpj2] < TmpMonomial2[Tmpj2 + 1]))
		  {
		    unsigned long Tmp = TmpMonomial2[Tmpj2 + 1];
		    TmpMonomial2[Tmpj2 + 1] = TmpMonomial2[Tmpj2];
		    TmpMonomial2[Tmpj2] = Tmp;
		    ++Tmpj2;
		  }
		TmpState = this->ConvertFromMonomial(TmpMonomial2);
		if ((TmpState <= MaxRoot) && (TmpState > CurrentPartition))
		  {
		    long TmpIndex = this->FermionBasis->FindStateIndex(TmpState, TmpMonomial2[0] + ReducedNbrBosons);
		    if ((TmpIndex < this->HilbertSpaceDimension) && (TmpIndex < MinIndex))
		      {
			MinIndex = TmpIndex;
		      }
		  }
	      }
	  }
      //      this->PrintState(cout, i) << " : " << MinIndex << endl;
      if ((MinIndex < this->LargeHilbertSpaceDimension) && (LargestDistance < (i - MinIndex)))
	{
	  LargestDistance = (i - MinIndex);
	}
    }
  cout << "largest distance = " << LargestDistance << "   (" << ((LargestDistance * 100.0) / ((double) this->LargeHilbertSpaceDimension))<< "%,  dim = " << this->LargeHilbertSpaceDimension << ")" << endl; 
}
  
