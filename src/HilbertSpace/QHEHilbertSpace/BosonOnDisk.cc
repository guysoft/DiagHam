////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class of bosons on disk                         //
//                                                                            //
//                        last modification : 04/06/2002                      //
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
#include "HilbertSpace/QHEHilbertSpace/BosonOnDisk.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SpinQuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"

#include <math.h>
#include <iostream>


using std::cout;
using std::endl;


// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value

BosonOnDisk::BosonOnDisk (int nbrBosons, int totalLz)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->TotalLz, this->TotalLz);
  this->Flag.Initialize();
  this->StateDescription = new int* [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrBosons, this->TotalLz, this->TotalLz, this->TotalLz, 0);
  this->KeyMultiplicationTable = new int [this->TotalLz + 1];
  this->GenerateLookUpTable(0);
  this->TemporaryState = new int [this->TotalLz + 1];
  this->ProdATemporaryState = new int [this->TotalLz + 1];
#ifdef __DEBUG__
  int UsedMemory = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    UsedMemory += (this->StateLzMax[i] + 1) * sizeof(int) + sizeof(int*);
  UsedMemory += (this->TotalLz + 1) * sizeof(int);
  UsedMemory += this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += ((this->TotalLz + 1) * this->IncNbrBosons) * sizeof(int);
  UsedMemory += this->HilbertSpaceDimension * sizeof(int);
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnDisk::BosonOnDisk(const BosonOnDisk& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->LzMaxPosition = bosons.LzMaxPosition;
  this->Keys =  bosons.Keys;
  this->Flag = bosons.Flag;
  this->TemporaryState = new int [this->TotalLz + 1];
  this->ProdATemporaryState = new int [this->TotalLz + 1];
  this->KeyMultiplicationTable = bosons.KeyMultiplicationTable;
}

// destructor
//

BosonOnDisk::~BosonOnDisk ()
{
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->Keys;
      delete[] this->LzMaxPosition;
      delete[] this->KeyMultiplicationTable;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnDisk& BosonOnDisk::operator = (const BosonOnDisk& bosons)
{
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
    }
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->Flag = bosons.Flag;
  this->TemporaryState = new int [this->TotalLz + 1];
  this->ProdATemporaryState = new int [this->TotalLz + 1];
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnDisk::Clone()
{
  return new BosonOnDisk(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnDisk::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnDisk::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnDisk::ExtractSubspace (AbstractQuantumNumber& q, 
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

int BosonOnDisk::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int LzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  if ((n1 > LzMax) || (n2 > LzMax) || (State[n1] == 0) || (State[n2] == 0) || ((n1 == n2) && (State[n1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = LzMax;
  if (LzMax < m1)
    NewLzMax = m1;
  if (LzMax < m2)
    NewLzMax = m2;
  int i = 0;
  for (; i <= LzMax; ++i)
    this->TemporaryState[i] = State[i];
  for (; i <= NewLzMax; ++i)
    this->TemporaryState[i] = 0;
  coefficient = this->TemporaryState[n2];
  --this->TemporaryState[n2];
  coefficient *= this->TemporaryState[n1];
  --this->TemporaryState[n1];
  ++this->TemporaryState[m2];
  coefficient *= this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  int DestIndex = this->FindStateIndex(this->TemporaryState, NewLzMax);
  return DestIndex;
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnDisk::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  int CurrentLzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  --nbrIndices;
  for (int i = 0; i <= nbrIndices; ++i)
    {
      if ((n[i] > CurrentLzMax) || (State[n[i]] == 0))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  int NewLzMax = CurrentLzMax;
 
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->TemporaryState[i] = State[i];
  for (; i <= this->TotalLz; ++i)
    this->TemporaryState[i] = 0;
  coefficient = 1.0;
  for (i = nbrIndices; i >= 0; --i)
    {
      if (this->TemporaryState[n[i]] == 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension; 	    
	}
      coefficient *= (double) this->TemporaryState[n[i]];
      --this->TemporaryState[n[i]];
    }
  for (i = nbrIndices; i >= 0; --i)
    {
      ++this->TemporaryState[m[i]];
      coefficient *= (double) this->TemporaryState[m[i]];
      if (m[i] > NewLzMax)
	NewLzMax = m[i];
    }
  coefficient = sqrt(coefficient);
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  int DestIndex = this->FindStateIndex(this->TemporaryState, NewLzMax);
  return DestIndex;
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnDisk::AA (int index, int n1, int n2)
{
  int CurrentLzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax) || (State[n1] == 0) || (State[n2] == 0) || ((n1 == n2) && (State[n1] == 1)))
    {
      return 0.0;
    }
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->ProdATemporaryState[i] = State[i];
  double Coefficient = this->ProdATemporaryState[n2];
  --this->ProdATemporaryState[n2];
  Coefficient *= this->ProdATemporaryState[n1];
  --this->ProdATemporaryState[n1];
  for (i = CurrentLzMax + 1; i <= this->TotalLz; ++i)
    this->ProdATemporaryState[i] = 0;
  return Coefficient;
}
// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnDisk::ProdA (int index, int* n, int nbrIndices)
{
  int CurrentLzMax = this->StateLzMax[index];
  int* State = this->StateDescription[index];
  int i = 0;
  for (; i <= CurrentLzMax; ++i)
    this->ProdATemporaryState[i] = State[i];
  int TmpCoefficient = 1;
  for (i = nbrIndices - 1; i >= 0; --i)
    {
      if (n[i] > CurrentLzMax)
        {
          return 0.0;
        }
      int& Tmp = this->ProdATemporaryState[n[i]];
      if (Tmp == 0)
        {
          return 0.0;
        }
      TmpCoefficient *= Tmp;
      --Tmp;
    }
  for (i = CurrentLzMax + 1; i <= this->TotalLz; ++i)
    this->ProdATemporaryState[i] = 0;
  return sqrt((double) TmpCoefficient);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnDisk::AdAd (int m1, int m2, double& coefficient)
{
  int i = 0;
  for (; i <= this->TotalLz; ++i)
    {
      this->TemporaryState[i] = this->ProdATemporaryState[i];
    }
  ++this->TemporaryState[m2];
  coefficient *= this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  int NewLzMax = this->TotalLz;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnDisk::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  int i = 0;
  for (; i <= this->TotalLz; ++i)
    {
      this->TemporaryState[i] = this->ProdATemporaryState[i];
    }
  int TmpCoefficient = 1;
  for (i = 0; i < nbrIndices; ++i)
    TmpCoefficient *= ++this->TemporaryState[m[i]];
  coefficient = sqrt((double) TmpCoefficient);
  int NewLzMax = this->TotalLz;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  return this->FindStateIndex(this->TemporaryState, NewLzMax);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnDisk::AdA (int index, int m)
{
  if (this->StateLzMax[index] < m)
    return 0.0;
  else
    return ((double) this->StateDescription[index][m]);
}

// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

int BosonOnDisk::FindStateIndex(int* stateDescription, int lzmax)
{
  int TmpKey = this->GenerateKey(stateDescription, lzmax);
  int i;
  int* TmpStateDescription;
  int Start = this->LzMaxPosition[lzmax * this->IncNbrBosons + stateDescription[lzmax]];
  int* TmpKey2 = &(this->Keys[Start]);
  for (; Start < this->HilbertSpaceDimension; ++Start)    
    {
      if ((*TmpKey2) == TmpKey)
	{
	  i = 0;
	  TmpStateDescription = this->StateDescription[Start];
	  while (i <= lzmax)
	    {
	      if (stateDescription[i] != TmpStateDescription[i])
		i = lzmax + 2;
	      ++i;
	    }
	  if (i == (lzmax + 1))
	    return Start;
	}
      ++TmpKey2;
    }
  return this->HilbertSpaceDimension;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnDisk::PrintState (ostream& Str, int state)
{
  int* TmpState = this->StateDescription[state];
  int Max = this->StateLzMax[state];
  int i = 0;
  for (; i <= Max; ++i)
    Str << TmpState[i] << " ";
  for (; i <= this->TotalLz; ++i)
    Str << "0 ";
  Str << " key = " << this->Keys[state] << " lzmax position = " << this->LzMaxPosition[Max * (this->NbrBosons + 1) + TmpState[Max]]
      << " position = " << FindStateIndex(TmpState, Max);
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson in the state
// currentLzMax = momentum maximum value for bosons that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int BosonOnDisk::GenerateStates(int nbrBosons, int lzMax, int currentLzMax, int totalLz, int pos)
{
  if ((nbrBosons == 0) || ((nbrBosons * currentLzMax) < totalLz) || (pos == this->HilbertSpaceDimension))
    {
      return pos;
    }
  if ((nbrBosons * currentLzMax) == totalLz)
    {
      this->StateDescription[pos] = new int [lzMax + 1];
      int* TmpState = this->StateDescription[pos];
      for (int i = 0; i <= lzMax; ++i)
	TmpState[i] = 0;
      TmpState[currentLzMax] = nbrBosons;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }
  if ((currentLzMax == 0) || (totalLz == 0))
    {
      this->StateDescription[pos] = new int [lzMax + 1];
      int* TmpState = this->StateDescription[pos];
      for (int i = 1; i <= lzMax; ++i)
	TmpState[i] = 0;
      TmpState[0] = nbrBosons;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }

  int TmpTotalLz = totalLz / currentLzMax;
  int TmpNbrBosons = nbrBosons - TmpTotalLz;
  TmpTotalLz = totalLz - TmpTotalLz * currentLzMax;
  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpPos = pos;
  while (TmpNbrBosons < nbrBosons)
    {
      TmpPos = this->GenerateStates(TmpNbrBosons, lzMax, ReducedCurrentLzMax, TmpTotalLz, pos);
      for (int i = pos; i < TmpPos; i++)
	this->StateDescription[i][currentLzMax] = nbrBosons - TmpNbrBosons;
      ++TmpNbrBosons;
      pos = TmpPos;
      TmpTotalLz += currentLzMax;
    }
  if (lzMax == currentLzMax)
    return this->GenerateStates(nbrBosons, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, pos);
  else
    return this->GenerateStates(nbrBosons, lzMax, ReducedCurrentLzMax, totalLz, pos);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnDisk::GenerateLookUpTable(int memeory)
{
  this->Keys = new int [this->HilbertSpaceDimension];
  for (int i = 0; i <= this->TotalLz; ++i)
    this->KeyMultiplicationTable[i] = i * i * this->IncNbrBosons;

  this->LzMaxPosition = new int [(this->TotalLz + 1) * this->IncNbrBosons];
  int CurrentLzMax = this->TotalLz;
  int CurrentNbrLzMax = 1;
  this->LzMaxPosition[CurrentLzMax * (this->NbrBosons +1) + CurrentNbrLzMax] = 0; 
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->Keys[i] = this->GenerateKey(this->StateDescription[i], this->StateLzMax[i]);
      if (CurrentLzMax != this->StateLzMax[i])
	{
	  CurrentLzMax = this->StateLzMax[i];
	  CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	  this->LzMaxPosition[CurrentLzMax * this->IncNbrBosons + CurrentNbrLzMax] = i; 
	}
      else
	if (this->StateDescription[i][CurrentLzMax] != CurrentNbrLzMax)
	  {
	    CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	    this->LzMaxPosition[CurrentLzMax * this->IncNbrBosons + CurrentNbrLzMax] = i; 
	  }
    }
}

// generate look-up table associated to current Hilbert space
// 
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = key associated to the state

int BosonOnDisk::GenerateKey(int* stateDescription, int lzmax)
{
  int Key = 0;
  for (int i = 0; i <= lzmax; ++i)
    {
      Key += this->KeyMultiplicationTable[i] * stateDescription[i];
    }
  return Key;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension

int BosonOnDisk::EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  if ((nbrBosons == 0) || ((nbrBosons * lzMax) < totalLz))
    return 0;
  if (((nbrBosons * lzMax) == totalLz) || (lzMax == 0) || (totalLz == 0))
    return 1;
  int TmpDim = 0;
  while (totalLz >= 0)
    {
      TmpDim += this->EvaluateHilbertSpaceDimension(nbrBosons, lzMax - 1, totalLz);
      --nbrBosons;
      totalLz -= lzMax;
    }
  return TmpDim;
}

void BosonOnDisk::ErasthothenesSlieve(int maxNumber, int nbrWantedPrime, int factor)
{
   maxNumber &= ~((int) 0x1);
  int Size = (int) sqrt((double) maxNumber);
  int NbrPrime = 1;
  int* PrimeNumbers = new int [Size];
  int NbrKeyFactor = 1;
  this->KeyMultiplicationTable[0] = 0;
  PrimeNumbers[0] = 3;
  int Lim;
  int Lim2 = factor;
  int j;  
  for (int i = 5; i < maxNumber; i +=2)
    {
      Lim = (int) sqrt((double) i);
      for (j = 0; j < NbrPrime; ++j)
	if ((j > Lim) || ((i % PrimeNumbers[j]) == 0))
	  {
	    j = NbrPrime + 1;
	  }
      if (j == NbrPrime)
	{
	  if(i < Size)
	    {
	      PrimeNumbers[NbrPrime] = i;
	      ++NbrPrime;
	    }
	  if (i > Lim2)
	    {
	      this->KeyMultiplicationTable[NbrKeyFactor] = i;
	      ++NbrKeyFactor;
	      Lim2 += i * factor;
	      if (NbrKeyFactor == nbrWantedPrime)
		i = maxNumber;
	    }
	}
    }
  for (int i = 0; i <= nbrWantedPrime; ++i)
    cout << this->KeyMultiplicationTable[i] << " ";
  cout << endl;
  delete[] PrimeNumbers;
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex BosonOnDisk::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
					     int firstComponent, int nbrComponent)
{
  Complex Value;
  Complex Tmp;
  ComplexMatrix Perm(this->NbrBosons, this->NbrBosons);
  ComplexMatrix Functions(this->TotalLz + 1, this->NbrBosons);
  RealVector TmpCoordinates(2);
  int* Indices = new int [this->NbrBosons];
  int Pos;
  int Lz;
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      TmpCoordinates[0] = position[j << 1];
      TmpCoordinates[1] = position[1 + (j << 1)];
      for (int i = 0; i <= this->TotalLz; ++i)
	{
	  basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	  Functions[j].Re(i) = Tmp.Re;
	  Functions[j].Im(i) = Tmp.Im;
	}
    }
  double* Factors = new double [this->NbrBosons + 1];
  Factors[0] = 1.0;
  Factors[1] = 1.0;
  for (int i = 2; i <= this->NbrBosons; ++i)
    Factors[i] = Factors[i - 1] / sqrt((double) i);
  double TmpFactor;
  int* ChangeBitSign;
  int* ChangeBit;
  int TmpStateDescription;
  Perm.EvaluateFastPermanentPrecalculationArray(ChangeBit, ChangeBitSign);
//  cout << Functions << endl << endl;
  int LastComponent = firstComponent + nbrComponent;
  for (int k = firstComponent; k < LastComponent; ++k)
    {
      Pos = 0;
      Lz = 0;
      TmpFactor = state[k] * Factors[this->NbrBosons];
      while (Pos < this->NbrBosons)
	{
	  TmpStateDescription = this->StateDescription[k][Lz];
	  if (TmpStateDescription != 0)
	    {
	      TmpFactor *= Factors[TmpStateDescription];
	      for (int j = 0; j < TmpStateDescription; ++j)
		{
		  Indices[Pos] = Lz;
		  ++Pos;
		}
	    }
	  ++Lz;
	}
      for (int i = 0; i < this->NbrBosons; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i];	  
	  for (int j = 0; j < this->NbrBosons; ++j)
	    {
	      Perm[i].Re(j) = TmpColum2.Re(Indices[j]);
	      Perm[i].Im(j) = TmpColum2.Im(Indices[j]);
	    }
	}
//      Complex Truc = Perm.FastPermanent(ChangeBit, ChangeBitSign);
//      cout << Truc<< " " << TmpFactor <<  endl;
//      Value += Truc * TmpFactor;
//      cout << Perm.FastPermanent(ChangeBit, ChangeBitSign) << " " << TmpFactor << " " << (Perm.FastPermanent(ChangeBit, ChangeBitSign) * TmpFactor) << endl;
//      cout << Perm << endl;
      Value += Perm.FastPermanent(ChangeBit, ChangeBitSign) * TmpFactor;
    }
  delete[] ChangeBitSign;
  delete[] ChangeBit;
  delete[] Factors;
  delete[] Indices;
  return Value;
}

// evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = index of the coordinate that will be changed during the next time iteration
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex BosonOnDisk::EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
							      int nextCoordinates, int firstComponent, int nbrComponent)
{
  double* Factors = new double [this->NbrBosons + 1];
  Factors[0] = 1.0;
  Factors[1] = 1.0;
  for (int i = 2; i <= this->NbrBosons; ++i)
    Factors[i] = Factors[i - 1] / sqrt((double) i);
  double TmpFactor;
  Complex Value;
  Complex Tmp;
  Complex TmpPerm;
  int* Indices = new int [this->NbrBosons];
  int Pos;
  int Lz;
  int TmpStateDescription;
  int LastComponent = firstComponent + nbrComponent;
  if ((*(this->KeptCoordinates)) == -1)
    {
      ComplexMatrix Perm(this->NbrBosons, this->NbrBosons);
      ComplexMatrix Functions(this->TotalLz + 1, this->NbrBosons);
      RealVector TmpCoordinates(2);
      for (int j = 0; j < this->NbrBosons; ++j)
	{
	  TmpCoordinates[0] = position[j << 1];
	  TmpCoordinates[1] = position[1 + (j << 1)];
	  for (int i = 0; i <= this->TotalLz; ++i)
	    {
	      basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	      Functions[j].Re(i) = Tmp.Re;
	      Functions[j].Im(i) = Tmp.Im;
	    }
	}
      int* ChangeBitSign;
      int* ChangeBit;
      int TmpStateDescription;
      Perm.EvaluateFastPermanentPrecalculationArray(ChangeBit, ChangeBitSign, true);
      int LastComponent = firstComponent + nbrComponent;
      for (int k = firstComponent; k < LastComponent; ++k)
	{
	  Pos = 0;
	  Lz = 0;
	  TmpFactor = state[k] * Factors[this->NbrBosons];
	  while (Pos < this->NbrBosons)
	    {
	      TmpStateDescription = this->StateDescription[k][Lz];
	      if (TmpStateDescription != 0)
		{
		  TmpFactor *= Factors[TmpStateDescription];
		  for (int j = 0; j < TmpStateDescription; ++j)
		    {
		      Indices[Pos] = Lz;
		      ++Pos;
		    }
		}
	      ++Lz;
	    }
	  for (int i = 0; i < this->NbrBosons; ++i)
	    {
	      ComplexVector& TmpColum2 = Functions[i];	  
	      for (int j = 0; j < this->NbrBosons; ++j)
		{
		  Perm[i].Re(j) = TmpColum2.Re(Indices[j]);
		  Perm[i].Im(j) = TmpColum2.Im(Indices[j]);
		}
	    }
	  if (this->Minors[k] == 0)
	    {
	      this->Minors[k] = new Complex [this->NbrBosons];
	    }
	  Perm.FastPermanentMinorDevelopment(ChangeBit, ChangeBitSign, nextCoordinates, this->Minors[k]);
	  TmpPerm = 0.0;
	  for (int i = 0; i < this->NbrBosons; ++i)
	    TmpPerm += this->Minors[k][i] * Complex (Perm[nextCoordinates].Re(i), 
						     Perm[nextCoordinates].Im(i));
	  Value += TmpPerm * TmpFactor;
	}
      delete[] ChangeBitSign;
      delete[] ChangeBit;
      (*(this->KeptCoordinates)) = nextCoordinates;
    }
  else
    {
      Complex* Functions = new Complex[this->TotalLz + 1];
      RealVector TmpCoordinates(2);
      TmpCoordinates[0] = position[(*(this->KeptCoordinates)) << 1];
      TmpCoordinates[1] = position[1 + ((*(this->KeptCoordinates)) << 1)];
      for (int i = 0; i <= this->TotalLz; ++i)
	{
	  basis.GetFunctionValue(TmpCoordinates, Functions[i], i);
	}
      for (int k = firstComponent; k < LastComponent; ++k)
	{
	  Pos = 0;
	  Lz = 0;
	  TmpFactor = Factors[this->NbrBosons] * state[k];
	  while (Pos < this->NbrBosons)
	    {
	      TmpStateDescription = this->StateDescription[k][Lz];
	      if (TmpStateDescription != 0)
		{
		  TmpFactor *= Factors[TmpStateDescription];
		  for (int j = 0; j < TmpStateDescription; ++j)
		    {
		      Indices[Pos] = Lz;
		      ++Pos;
		    }
		}
	      ++Lz;
	    }
	  Complex* TmpMinors = this->Minors[k];
	  TmpPerm = 0.0;
	  for (int i = 0; i < this->NbrBosons; ++i)
	    TmpPerm += TmpMinors[i] * Functions[Indices[i]];
	  Value += TmpPerm * TmpFactor;
	}
      delete[] Functions;
      (*(this->KeptCoordinates)) = -1;
    }
  delete[] Factors;
  delete[] Indices;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void BosonOnDisk::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
  if ((timeCoherence == true) && (this->Minors == 0))
    {
      this->Minors = new Complex* [this->HilbertSpaceDimension];
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	this->Minors[i] = 0;
    }
}


