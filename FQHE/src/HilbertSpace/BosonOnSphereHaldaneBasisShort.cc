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

  RationalPolynomial* TmpNumerators = new RationalPolynomial[this->LargeHilbertSpaceDimension];
  RationalPolynomial* TmpDenominators = new RationalPolynomial[this->LargeHilbertSpaceDimension];		  
  Rational* Roots = new Rational[this->LargeHilbertSpaceDimension];

  this->GenerateSingleJackPolynomialCoefficient(jack, 0, TmpNumerators, TmpDenominators, Roots);
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->GenerateSingleJackPolynomialCoefficient(jack, i, TmpNumerators, TmpDenominators, Roots);
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
	      int Depth = -1;//-1;
	      bool FullSymbolicFlag = false;
	      bool SolvedFlag = false;
	      //	      while (((Depth <= symbolicDepth) || (symbolicDepth < 0)) && (SolvedFlag == false) && (FullSymbolicFlag == false))
		{
// 		  RationalPolynomial TmpNumerator;
// 		  RationalPolynomial TmpDenominator;		  
// 		  FullSymbolicFlag = this->GenerateSingleJackPolynomialCoefficient(jack, i, TmpNumerator, TmpDenominator, Depth);
		  FullSymbolicFlag = this->GenerateSingleJackPolynomialCoefficient(jack, i, TmpNumerators, TmpDenominators, Roots);
		  Rational Tmp = TmpNumerators[i].PolynomialEvaluate(InvAlpha);
		  cout << "--------------------------------" << endl
		       << "result = " << endl;
		  cout << TmpNumerators[i] << endl;
		  cout << TmpDenominators[i] << endl;
		  cout << Tmp.Num() << endl;
		  if (Tmp.Num() == 0l)
		    {
		      cout << TmpNumerators[i] << endl;
		      cout << TmpDenominators[i] << endl;
		      TmpNumerators[i].MonomialDivision(InvAlpha);
		      TmpDenominators[i].MonomialDivision(InvAlpha);
		      Tmp = TmpNumerators[i].PolynomialEvaluate(InvAlpha);
		      Tmp /= TmpDenominators[i].PolynomialEvaluate(InvAlpha);
		      jack[i] = Tmp;
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
	      Rational TmpR = TmpNumerators[i].PolynomialEvaluate(InvAlpha) / TmpDenominators[i].PolynomialEvaluate(InvAlpha);
	      cout << jack[i] << " " << TmpR << endl;
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
  if (index == 0l)
    {
      cout << "contribution from root" << endl;
      numerator = RationalPolynomial(0);
      denominator = RationalPolynomial(0);
      numerator[0] = 1l;
      denominator[0] = 1l;
      cout << numerator << endl;
      cout << denominator << endl;
      return true;
    }

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];
  
  Rational RhoRootInvAlphaCoef = 0l;
  Rational RhoRootConstCoef = 0l;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      RhoRootInvAlphaCoef -= TmpMonomial[j] * ((long) j);
      RhoRootConstCoef += TmpMonomial[j] * (TmpMonomial[j] - 1l);
    }
  int ReducedNbrBosons = this->NbrBosons - 1;
  
  Rational RhoInvAlphaCoef = 0l;
  Rational RhoConstCoef = 0l;

  unsigned long CurrentPartition = this->FermionBasis->StateDescription[index];
  this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[index], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      RhoInvAlphaCoef -= TmpMonomial[j] * ((long) j);
      RhoConstCoef += TmpMonomial[j] * (TmpMonomial[j] - 1l);
    }
  
  int Pos = 0;
  long* ConnectedIndices = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients  = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
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
		  {
		    ConnectedIndices[Pos] = TmpIndex;
		    ConnectedCoefficients[Pos] = Diff;
		    ++Pos;
		  }
	      }
	  }
      }

  int NbrConnected = 1l;
  if (Pos > 1)
    {
      SortArrayDownOrdering<long>(ConnectedIndices, ConnectedCoefficients, Pos);
//       for (int i = 0; i < Pos; ++i)
// 	cout << ConnectedIndices[i] << " (" << ConnectedCoefficients[i] << ") ";
//       cout << endl;
      int TmpIndex = 1;
      while (TmpIndex < Pos)
	{
	  while ((TmpIndex < Pos) && (ConnectedIndices[TmpIndex] == ConnectedIndices[TmpIndex - 1]))
	    ++TmpIndex;
	  if (TmpIndex < Pos)
	    ++NbrConnected;
	  ++TmpIndex;
	}
      cout << "NbrConnected=" << NbrConnected << endl;
      long* TmpConnectedIndices = new long [NbrConnected];
      long* TmpConnectedCoefficients  = new long [NbrConnected];
      TmpConnectedIndices[0] = ConnectedIndices[0];
      TmpConnectedCoefficients[0] = ConnectedCoefficients[0];
      TmpIndex = 1;
      NbrConnected = 1;
      while (TmpIndex < Pos)
	{
	  while ((TmpIndex < Pos) && (ConnectedIndices[TmpIndex] == ConnectedIndices[TmpIndex - 1]))
	    {
	      TmpConnectedCoefficients[NbrConnected - 1] += ConnectedCoefficients[TmpIndex];
	      ++TmpIndex;
	    }
	  if (TmpIndex < Pos)
	    {
	      TmpConnectedIndices[NbrConnected] = ConnectedIndices[TmpIndex];
	      TmpConnectedCoefficients[NbrConnected] = ConnectedCoefficients[TmpIndex];	   
	      ++NbrConnected;
	    }
	  ++TmpIndex;
	}
      delete[] ConnectedIndices;
      delete[] ConnectedCoefficients;
      ConnectedIndices = TmpConnectedIndices;
      ConnectedCoefficients = TmpConnectedCoefficients;
      for (int i = 0; i < NbrConnected; ++i)
 	cout << ConnectedIndices[i] << " (" << ConnectedCoefficients[i] << ") " << endl;
      cout << endl;
    }
  bool SymbolicFlag = false;
  if (depth == 0)
    {
      numerator = RationalPolynomial(1);
      denominator = RationalPolynomial(1);
      numerator[0] = 0l;
      numerator[1] = 0l;
      Rational Tmp2 = 1l / (RhoRootInvAlphaCoef - RhoInvAlphaCoef);
      for (int i = 0; i < NbrConnected; ++i)
	{
	  numerator[1] += (jack[ConnectedIndices[i]] * ConnectedCoefficients[i]) * Tmp2;
	  cout << "computing non symbolic contribution to " << index << " from component " << ConnectedIndices[i] << " at depth " << depth << ", no symbolic calculation " << endl;
	}
      denominator[0] = (RhoRootConstCoef - RhoConstCoef) * Tmp2;
      denominator[1] = 1l;
    }
  else
    {
      SymbolicFlag = true;
      RationalPolynomial* TmpNumerators = new RationalPolynomial[Pos];
      RationalPolynomial* TmpDenominators = new RationalPolynomial[Pos];
      cout << "entering calculation of component " <<  index << "(connected to " << NbrConnected << " components) " << endl;
      for (int i = 0; i < NbrConnected; ++i)
        {
	  bool TmpFlag = this->GenerateSingleJackPolynomialCoefficient(jack, ConnectedIndices[i], TmpNumerators[i], TmpDenominators[i], depth - 1);
	  cout << "computing contribution to " << index << " component " << ConnectedIndices[i] << " at depth " << depth << endl;
	  cout << "numerator = " << TmpNumerators[i] << endl;
	  cout << "denominator = " << TmpDenominators[i] << endl;
	  if (TmpFlag == false)
	    SymbolicFlag = TmpFlag;
	}
      denominator = RationalPolynomial(1);
      Rational Tmp2 = (RhoRootInvAlphaCoef - RhoInvAlphaCoef);
      denominator[0] = (RhoRootConstCoef - RhoConstCoef);// * Tmp2;
      denominator[1] = Tmp2;
      for (int i = 0; i < NbrConnected; ++i)
        {
	  denominator *= TmpDenominators[i];
	  for (int j= 0; j < NbrConnected; ++j)
	    {
	      if (i != j)
		{
		  TmpNumerators[i] *= TmpDenominators[j];
		}
	    }
	  TmpNumerators[i] *= ConnectedCoefficients[i];
	}
      numerator = RationalPolynomial();
      for (int i = 0; i < NbrConnected; ++i)
        {
	  numerator += TmpNumerators[i];
	  cout << "sum=" << i << endl
	       << numerator << endl
	       << TmpNumerators[i] << endl;
	}
      numerator.ShiftPowers(1);
      //      numerator *= Tmp2;
    }
  delete[] TmpMonomial;
  delete[] TmpMonomial2;
  delete[] ConnectedIndices;
  delete[] ConnectedCoefficients;
  return SymbolicFlag;
}

// compute a single coefficient of the Jack polynomial decomposition corresponding to the root partition, assuming only rational numbers occur and using (partial symbolic calculation)
//
// jack = vector where the ecomposition of the corresponding Jack polynomial on the unnormalized basis will be stored
// index = index of the component to compute
// numerator = reference on the polynomial where the numerator has to be stored
// denominator = reference on the polynomial where the denominator has to be stored
// depth = depth in the recurrence describing up to which point the symbolic calculation has to be performed 
// return value = true if a fully symbolic calculation has been performed

bool BosonOnSphereHaldaneBasisShort::GenerateSingleJackPolynomialCoefficient(RationalVector& jack, long index, RationalPolynomial* numerators, RationalPolynomial* denominators, Rational* roots)
{
  if (numerators[index].Defined())
    return true;
  if (index == 0l)
    {
      if (!numerators[0l].Defined())
	{
	  numerators[0l] = RationalPolynomial(0);
	  denominators[0l] = RationalPolynomial(0);
	  numerators[0l][0] = 1l;
	  denominators[0l][0] = 1l;
	  roots[0l] = Rational(0l, 0l);
	  return true;
	}
    }

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial2 = new unsigned long [this->NbrBosons];
  
  Rational RhoRootInvAlphaCoef = 0l;
  Rational RhoRootConstCoef = 0l;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      RhoRootInvAlphaCoef -= TmpMonomial[j] * ((long) j);
      RhoRootConstCoef += TmpMonomial[j] * (TmpMonomial[j] - 1l);
    }
  int ReducedNbrBosons = this->NbrBosons - 1;
  
  Rational RhoInvAlphaCoef = 0l;
  Rational RhoConstCoef = 0l;

  unsigned long CurrentPartition = this->FermionBasis->StateDescription[index];
  this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[index], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      RhoInvAlphaCoef -= TmpMonomial[j] * ((long) j);
      RhoConstCoef += TmpMonomial[j] * (TmpMonomial[j] - 1l);
    }
  
  int Pos = 0;
  long* ConnectedIndices = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
  long* ConnectedCoefficients  = new long [((this->NbrBosons * ReducedNbrBosons) >> 1) * (this->LzMax + 1)];
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
		  {
		    ConnectedIndices[Pos] = TmpIndex;
		    ConnectedCoefficients[Pos] = Diff;
		    ++Pos;
		  }
	      }
	  }
      }

  int NbrConnected = 1l;
  if (Pos > 1)
    {
      SortArrayDownOrdering<long>(ConnectedIndices, ConnectedCoefficients, Pos);
      int TmpIndex = 1;
      while (TmpIndex < Pos)
	{
	  while ((TmpIndex < Pos) && (ConnectedIndices[TmpIndex] == ConnectedIndices[TmpIndex - 1]))
	    ++TmpIndex;
	  if (TmpIndex < Pos)
	    ++NbrConnected;
	  ++TmpIndex;
	}
      cout << "NbrConnected=" << NbrConnected << endl;
      long* TmpConnectedIndices = new long [NbrConnected];
      long* TmpConnectedCoefficients  = new long [NbrConnected];
      TmpConnectedIndices[0] = ConnectedIndices[0];
      TmpConnectedCoefficients[0] = ConnectedCoefficients[0];
      TmpIndex = 1;
      NbrConnected = 1;
      while (TmpIndex < Pos)
	{
	  while ((TmpIndex < Pos) && (ConnectedIndices[TmpIndex] == ConnectedIndices[TmpIndex - 1]))
	    {
	      TmpConnectedCoefficients[NbrConnected - 1] += ConnectedCoefficients[TmpIndex];
	      ++TmpIndex;
	    }
	  if (TmpIndex < Pos)
	    {
	      TmpConnectedIndices[NbrConnected] = ConnectedIndices[TmpIndex];
	      TmpConnectedCoefficients[NbrConnected] = ConnectedCoefficients[TmpIndex];	   
	      ++NbrConnected;
	    }
	  ++TmpIndex;
	}
      delete[] ConnectedIndices;
      delete[] ConnectedCoefficients;
      ConnectedIndices = TmpConnectedIndices;
      ConnectedCoefficients = TmpConnectedCoefficients;
    }
  bool SymbolicFlag = false;

  SymbolicFlag = true;
  cout << "entering calculation of component " <<  index << "(connected to " << NbrConnected << " components) " << endl;
  for (int i = 0; i < NbrConnected; ++i)
    {
      if (!numerators[ConnectedIndices[i]].Defined())
        {
	  this->GenerateSingleJackPolynomialCoefficient(jack, ConnectedIndices[i], numerators, denominators, roots);
	}
    }

  denominators[index] = RationalPolynomial(1);
  Rational Tmp2 = (RhoRootInvAlphaCoef - RhoInvAlphaCoef);
  denominators[index][0] = (RhoRootConstCoef - RhoConstCoef);// * Tmp2;
  denominators[index][1] = Tmp2;
  numerators[index] = RationalPolynomial();
  for (int i = 0; i < NbrConnected; ++i)
    {
      RationalPolynomial TmpNumerator = numerators[ConnectedIndices[i]];
      for (int j= 0; j < NbrConnected; ++j)
	{
	  if (i != j)
	    {
	      TmpNumerator *= denominators[ConnectedIndices[j]];
	    }
	}
      TmpNumerator *= ConnectedCoefficients[i];
      numerators[index] += TmpNumerator;
      denominators[index] *= denominators[ConnectedIndices[i]];
      cout << numerators[index] << endl;
      for (int i = 1l; i < index; ++i)
	{
	  if (roots[i].Den() != 0l)
	    {
	      cout << roots[i] << endl;
	      Rational Tmp4 = numerators[index].PolynomialEvaluate(roots[i]);
	      if (Tmp4.Num() == 0l)
		{
		  cout << "root found" << endl;
		  numerators[index].MonomialDivision(roots[i]);
		  denominators[index].MonomialDivision(roots[i]);
		}
	    }
	}
    }
  numerators[index].ShiftPowers(1);

//   for (int i = 1l; i < index; ++i)
//     {
//       if (roots[i].Den() != 0l)
// 	{
// 	  Rational Tmp4 = numerators[index].PolynomialEvaluate(roots[i]);
// 	  if (Tmp4.Num() == 0l)
// 	    {
// 	      numerators[index].MonomialDivision(roots[i]);
// 	      denominators[index].MonomialDivision(roots[i]);
// 	    }
// 	}
//     }

  if (Tmp2.Num() != 0l)
    {
      denominators[index] /= Tmp2;
      numerators[index] /= Tmp2;
      Rational Tmp3 = RhoRootConstCoef - RhoConstCoef;
      Tmp3 /= -Tmp2;
      Rational Tmp4 = numerators[index].PolynomialEvaluate(Tmp3);
      if (Tmp4.Num() == 0l)
	{
	  numerators[index].MonomialDivision(Tmp3);
	  denominators[index].MonomialDivision(Tmp3);
	}
      roots[index] = Tmp3;
    }
  else
    {
      Tmp2 = (RhoRootConstCoef - RhoConstCoef);
      if (Tmp2.Num() != 0)
	{
	  denominators[index] /= Tmp2;
	  numerators[index] /= Tmp2;
	}
      roots[index] = Rational(0l, 0l);
    }

  cout << "component " << index << " : " << endl;
  cout << "num  = " << numerators[index] << endl;
  cout << "den  = " << denominators[index] << endl;

  delete[] TmpMonomial;
  delete[] TmpMonomial2;
  delete[] ConnectedIndices;
  delete[] ConnectedCoefficients;
  return SymbolicFlag;
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

  double RhoRoot = 0.0;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += TmpMonomial[j] * (TmpMonomial[j] - 1.0 - InvAlpha * ((double) j));

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
  Rational InvAlpha (2l * alphaDenominator, alphaNumerator);

  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];

  Rational RhoRoot = 0;
  unsigned long MaxRoot = this->FermionBasis->StateDescription[0];
  this->ConvertToMonomial(MaxRoot, this->FermionBasis->StateLzMax[0], TmpMonomial);
  for (int j = 0; j < this->NbrBosons; ++j)
    RhoRoot += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));

  jack[0] = RhoRoot;
  for (long i = 1; i < this->LargeHilbertSpaceDimension; ++i)
    {
      Rational Rho = 0l;
      unsigned long CurrentPartition = this->FermionBasis->StateDescription[i];
      this->ConvertToMonomial(CurrentPartition, this->FermionBasis->StateLzMax[i], TmpMonomial);
      for (int j = 0; j < this->NbrBosons; ++j)
	Rho += TmpMonomial[j] * ((TmpMonomial[j] - 1l) - InvAlpha * ((long) j));
      if (RhoRoot == Rho)
	jack[i] = Rho;
      else
	jack[i] = 0l;
      }
  delete[] TmpMonomial;

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
  
