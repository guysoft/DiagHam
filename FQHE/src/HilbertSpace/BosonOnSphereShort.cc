////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on sphere for system size such that            //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 10/10/2007                      //
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
#include "HilbertSpace/BosonOnSphereShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"

#include <math.h>


using std::cout;
using std::endl;


// default constructor
//

BosonOnSphereShort::BosonOnSphereShort ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson

BosonOnSphereShort::BosonOnSphereShort (int nbrBosons, int totalLz, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  this->FermionBasis = new FermionOnSphere(nbrBosons, totalLz, lzMax + nbrBosons - 1);
  this->HilbertSpaceDimension = this->FermionBasis->HilbertSpaceDimension;

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

BosonOnSphereShort::BosonOnSphereShort(const BosonOnSphereShort& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];  
}

// destructor
//

BosonOnSphereShort::~BosonOnSphereShort ()
{
  if (this->FermionBasis != 0)
    delete this->FermionBasis;
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
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereShort& BosonOnSphereShort::operator = (const BosonOnSphereShort& bosons)
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
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereShort::Clone()
{
  return new BosonOnSphereShort(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereShort::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereShort::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereShort::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShort::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  if ((n1 > this->TemporaryStateLzMax) || (n2 > this->TemporaryStateLzMax) || (this->TemporaryState[n1] == 0) || (this->TemporaryState[n2] == 0) || ((n1 == n2) && (this->TemporaryState[n1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = this->TemporaryState[n2];
  --this->TemporaryState[n2];
  coefficient *= this->TemporaryState[n1];
  --this->TemporaryState[n1];
  ++this->TemporaryState[m2];
  coefficient *= this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShort::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  --nbrIndices;
  for (int i = 0; i <= nbrIndices; ++i)
    {
      if ((n[i] > this->TemporaryStateLzMax) || (this->TemporaryState[n[i]] == 0))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  coefficient = 1.0;
  for (int i = nbrIndices; i >= 0; --i)
    {
      if (this->TemporaryState[n[i]] == 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension; 	    
	}
      coefficient *= (double) this->TemporaryState[n[i]];
      --this->TemporaryState[n[i]];
    }
  for (int i = nbrIndices; i >= 0; --i)
    {
      ++this->TemporaryState[m[i]];
      coefficient *= (double) this->TemporaryState[m[i]];
    }
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnSphereShort::AA (int index, int n1, int n2)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
  if ((n1 > this->ProdATemporaryStateLzMax) || (n2 > this->ProdATemporaryStateLzMax) || 
      (this->ProdATemporaryState[n1] == 0) || (this->ProdATemporaryState[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryState[n1] == 1)))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryState[n2];
  --this->ProdATemporaryState[n2];
  Coefficient *= this->ProdATemporaryState[n1];
  --this->ProdATemporaryState[n1];
  for (int i = this->ProdATemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0ul;
  return sqrt(Coefficient);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphereShort::ProdA (int index, int* n, int nbrIndices)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
  int TmpCoefficient = 1;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (n[i] > this->ProdATemporaryStateLzMax)
	return 0.0;
      unsigned long& Tmp = this->ProdATemporaryState[n[i]];
      if (Tmp == 0)
	return 0.0;
      TmpCoefficient *= Tmp;
      --Tmp;
    }
  for (int i = this->ProdATemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0;
  return sqrt((double) TmpCoefficient);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShort::AdAd (int m1, int m2, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  ++this->TemporaryState[m2];
  coefficient = this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShort::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  int TmpCoefficient = 1;
  for (int i = 0; i < nbrIndices; ++i)
    TmpCoefficient *= ++this->TemporaryState[m[i]];
  coefficient = sqrt((double) TmpCoefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereShort::AdA (int index, int m)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  if (this->TemporaryStateLzMax < m)  
    return 0.0;
  return (double) (this->TemporaryState[m]);  
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereShort::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[state], this->FermionBasis->StateLzMax[state], this->TemporaryState, this->TemporaryStateLzMax);
  int i = 0;
  for (; i <= this->TemporaryStateLzMax; ++i)
    Str << this->TemporaryState[i] << " ";
  for (; i <= this->LzMax; ++i)
    Str << "0 ";
  Str << "   lzmax = " << this->TemporaryStateLzMax;
  return Str;
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  BosonOnSphereShort::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState)
{  
  RealSymmetricMatrix TmpDensityMatrix;
  return TmpDensityMatrix;	  
//   if (subsytemSize <= 0)
//     {
//       if ((lzSector == 0) && (nbrBosonSector == 0))
// 	{
// 	  RealSymmetricMatrix TmpDensityMatrix(1);
// 	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
// 	  return TmpDensityMatrix;
// 	}
//       else
// 	{
// 	  RealSymmetricMatrix TmpDensityMatrix;
// 	  return TmpDensityMatrix;	  
// 	}
//     }
//   if (subsytemSize > this->LzMax)
//     {
//       if ((lzSector == this->TotalLz) && (nbrBosonSector == this->NbrBosons))
// 	{
// 	  RealSymmetricMatrix TmpDensityMatrix(this->HilbertSpaceDimension);
// 	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	    for (int j = i; j < this->HilbertSpaceDimension; ++j)
// 	      TmpDensityMatrix.SetMatrixElement(i, j, groundState[i] * groundState[j]);
// 	}
//       else
// 	{
// 	  RealSymmetricMatrix TmpDensityMatrix;
// 	  return TmpDensityMatrix;	  
// 	}
//     }
//   if (subsytemSize == 1)
//     {
//       if (lzSector == 0)
// 	{
// 	  double TmpValue = 0;
// 	  for (int MinIndex = 0; MinIndex < this->HilbertSpaceDimension; ++MinIndex)
// 	    if (this->StateDescription[MinIndex][0] == nbrBosonSector)
// 	       TmpValue += groundState[MinIndex] * groundState[MinIndex];
// 	  RealSymmetricMatrix TmpDensityMatrix(1);
// 	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
// 	  return TmpDensityMatrix;
// 	}
//       else
// 	{
// 	  RealSymmetricMatrix TmpDensityMatrix;
// 	  return TmpDensityMatrix;	  
// 	}      
//     }
//   if (nbrBosonSector == 0)
//     {
//       if (lzSector == 0)
// 	{
// 	  double TmpValue = 0;
// 	  int MinIndex = 0;
// 	  while ((MinIndex < this->HilbertSpaceDimension) && (this->StateLzMax[MinIndex] >= subsytemSize))
// 	    {
// 	      int TmpPos = 0;
// 	      int* TmpState = this->StateDescription[MinIndex];
// 	      while ((TmpPos < subsytemSize) && (TmpState[TmpPos] == 0))
// 		++TmpPos;
// 	      if (TmpPos == subsytemSize)
// 		{
// 		  TmpValue += groundState[MinIndex] * groundState[MinIndex];
// 		}
// 	      ++MinIndex;
// 	    }
// 	  RealSymmetricMatrix TmpDensityMatrix(1);
// 	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
// 	  return TmpDensityMatrix;
// 	}
//       else
// 	{
// 	  RealSymmetricMatrix TmpDensityMatrix;
// 	  return TmpDensityMatrix;	  
// 	}
//     }

//   int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
//   int ShiftedLzSector = (lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1;
//   int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
//   if (ShiftedLzComplementarySector < 0)
//     {
//       RealSymmetricMatrix TmpDensityMatrix;
//       return TmpDensityMatrix;	  
//     }
//   int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
//   int MinIndex = 0;
//   int MaxIndex = this->HilbertSpaceDimension - 1;
//   if (nbrBosonSector == 1)
//     {
//       double TmpValue = 0.0;
//       int TmpLzMax = this->StateLzMax[MinIndex];
//       while ((MinIndex <= MaxIndex) && (subsytemSize <= TmpLzMax))
// 	{
// 	  int* TmpStateDescription = this->StateDescription[MinIndex];
// 	  if (TmpStateDescription[ShiftedLzSector] == 1)
// 	    {	      
// 	      int TmpPos = 0;
// 	      int TmpNbrBosons = 0;
// 	      while (TmpPos < subsytemSize)
// 		TmpNbrBosons += TmpStateDescription[TmpPos++];
// 	      if (TmpNbrBosons == 1)
// 		TmpValue += groundState[MinIndex] * groundState[MinIndex];	    
// 	    }
// 	  ++MinIndex;
// 	  if (MinIndex <= MaxIndex)
// 	    TmpLzMax = this->StateLzMax[MinIndex];
// 	}
//       RealSymmetricMatrix TmpDensityMatrix(1, true);
//       TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
//       return TmpDensityMatrix;
//     }
//   if (NbrBosonsComplementarySector == 0)
//     {
//       if (ShiftedLzComplementarySector != 0)
// 	{
// 	  RealSymmetricMatrix TmpDensityMatrix;
// 	  return TmpDensityMatrix;	  
// 	}
//       BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
//       cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
//       RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
//       MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
//       double TmpValue;
//       for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
// 	{
// 	  TmpValue = groundState[MinIndex + i];
// 	  for (int j = i; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
// 	    TmpDensityMatrix.SetMatrixElement(i, j, TmpValue * groundState[MinIndex + j]);
// 	}
//       return TmpDensityMatrix;
//     }


//   int TmpNbrBosons;
//   int TmpTotalLz;
//   int TmpIndex;
//   BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  
//   cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
//   int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
//   long TmpNbrNonZeroElements = 0;
//   int TmpComplementarySubsystemLzMax = this->StateLzMax[MinIndex];
//   while ((MinIndex <= MaxIndex) && (TmpComplementarySubsystemLzMax >= subsytemSize))
//     {
//       int* TmpComplementarySubsystem = this->StateDescription[MinIndex];
//       TmpIndex = MinIndex + 1;
//       int TmpPos = TmpComplementarySubsystemLzMax;	  
//       int TmpLzMax = TmpPos;
//       while ((TmpIndex <= MaxIndex) && (TmpPos == TmpLzMax))
// 	{
// 	  if (TmpLzMax == this->StateLzMax[TmpIndex])
// 	    {
// 	      TmpPos = subsytemSize;
// 	      while ((TmpPos <= TmpLzMax) && (this->StateDescription[TmpIndex][TmpPos] == TmpComplementarySubsystem[TmpPos]))
// 		++TmpPos;
// 	      if (TmpPos > TmpLzMax)
// 		{
// 		  ++TmpIndex;
// 		  --TmpPos;
// 		}
// 	      else
// 		{
// 		  TmpPos = -1;
// 		}
// 	    }
// 	  else
// 	    {
// 	      TmpPos = -1;
// 	    }
// 	}
//       TmpNbrBosons = 0;
//       TmpTotalLz = 0;
//       TmpPos = subsytemSize;	  
//       while (TmpPos <= TmpComplementarySubsystemLzMax)
// 	{
// 	  TmpNbrBosons += TmpComplementarySubsystem[TmpPos];
// 	  TmpTotalLz += TmpComplementarySubsystem[TmpPos] * TmpPos;
// 	  ++TmpPos;
// 	}
//       if ((TmpNbrBosons == NbrBosonsComplementarySector) && (ShiftedLzComplementarySector == TmpTotalLz))
// 	{
// 	  int Pos = 0;
// 	  for (int i = MinIndex; i < TmpIndex; ++i)
// 	    {
// 	      int* TmpState = this->StateDescription[i];
// 	      int TmpLzMax = subsytemSize - 1;
// 	      while (TmpState[TmpLzMax] == 0) 
// 		--TmpLzMax;
// 	      TmpStatePosition[Pos] = TmpDestinationHilbertSpace.FindStateIndex(TmpState, TmpLzMax);
// 	      ++Pos;
// 	    }
// 	  int Pos2;
// 	  Pos = 0;
// 	  int Pos3;
// 	  double TmpValue;
// 	  for (int i = MinIndex; i < TmpIndex; ++i)
// 	    {
// 	      Pos2 = 0;
// 	      Pos3 = TmpStatePosition[Pos];
// 	      TmpValue = groundState[i];
// 	      for (int j = MinIndex; j < TmpIndex; ++j)
// 		{
// 		  if (Pos3 <=  TmpStatePosition[Pos2])
// 		    {
// 		      TmpDensityMatrix.AddToMatrixElement(Pos3, TmpStatePosition[Pos2], TmpValue * groundState[j]);
// 		      ++TmpNbrNonZeroElements;
// 		    }
// 		  ++Pos2;
// 		}
// 	      ++Pos;
// 	    }
// 	}
//       MinIndex = TmpIndex;
//       if (MinIndex <= MaxIndex)
// 	TmpComplementarySubsystemLzMax = StateLzMax[MinIndex];
//     }
//   delete[] TmpStatePosition;
//   if (TmpNbrNonZeroElements > 0)	
//     return TmpDensityMatrix;
//   else
//     {
//       RealSymmetricMatrix TmpDensityMatrixZero;
//       return TmpDensityMatrixZero;
//     }
}
