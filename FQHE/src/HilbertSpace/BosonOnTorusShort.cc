////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of bosons on torus                         //
//                                                                            //
//                        last modification : 03/09/2002                      //
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
#include "HilbertSpace/BosonOnTorusShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "MathTools/IntegerAlgebraTools.h"

#include <math.h>


using std::cout;
using std::endl;


// basic constructor
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson

BosonOnTorusShort::BosonOnTorusShort (int nbrBosons, int maxMomentum)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->MomentumConstraintFlag = false;

  this->TemporaryState = new unsigned long [this->MaxMomentum + 1];
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->MaxMomentum);
  cout << "dim = " << this->HilbertSpaceDimension << endl;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateMaxMomentum = new int [this->HilbertSpaceDimension];
  cout << this->GenerateStates(this->NbrBosons, this->MaxMomentum - 1, this->MaxMomentum - 1, 0) << endl;
  this->GenerateLookUpTable(0);
#ifdef __DEBUG__
  unsigned long UsedMemory = 0;
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * ((unsigned long) this->LookUpTableMemorySize) * sizeof(int);
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

// constructor with a constraint of the total momentum of states
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// momentumConstraint = index of the momentum orbit

BosonOnTorusShort::BosonOnTorusShort (int nbrBosons, int maxMomentum, int momentumConstraint)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrLzValue = this->MaxMomentum + 1;
  this->MomentumConstraint = momentumConstraint;
  this->MomentumConstraintFlag = true;
  this->GCDMaxMomentum = FindGCD(this->NbrBosons, this->MaxMomentum);
  this->TemporaryState = new unsigned long [this->MaxMomentum + 1];

  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->MaxMomentum);
  cout << "dim = " << this->HilbertSpaceDimension << endl;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateMaxMomentum = new int [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->MaxMomentum - 1, this->MaxMomentum - 1, 0, 0);
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  unsigned long UsedMemory = 0;
  UsedMemory += ((unsigned long) this->HilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * ((unsigned long) this->LookUpTableMemorySize) * sizeof(int);
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

BosonOnTorusShort::BosonOnTorusShort(const BosonOnTorusShort& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrLzValue = bosons.NbrLzValue;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;
  this->TemporaryState = new unsigned long [this->MaxMomentum + 1];
  this->Flag = bosons.Flag;
}

// destructor
//

BosonOnTorusShort::~BosonOnTorusShort ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;
    }
  delete[] this->TemporaryState;
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnTorusShort& BosonOnTorusShort::operator = (const BosonOnTorusShort& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;
      delete[] this->TemporaryState;
    }
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrLzValue = bosons.NbrLzValue;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;
  this->TemporaryState = new unsigned long [this->MaxMomentum + 1];
  this->Flag = bosons.Flag;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnTorusShort::Clone()
{
  return new BosonOnTorusShort(*this);
}

// get momemtum value of a given state
//
// index = state index
// return value = state momentum

int BosonOnTorusShort::GetMomentumValue(int index)
{
  this->FermionToBoson(this->StateDescription[index], this->StateMaxMomentum[index]  + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateMaxMomentum);
  int Momentum = 0;
  for (int i = 0; i <= this->TemporaryStateMaxMomentum; ++i)
    {
      Momentum += (this->TemporaryState[i] * i);
    }
  return (Momentum % this->MaxMomentum);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnTorusShort::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->MomentumConstraintFlag == false)
    {
      for (int i = 0; i < this->MaxMomentum; ++i)
	TmpList += new PeriodicMomentumQuantumNumber (i, this->MaxMomentum);
    }
  else
    TmpList += new PeriodicMomentumQuantumNumber (this->MomentumConstraint, this->MaxMomentum);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnTorusShort::GetQuantumNumber (int index)
{
  if (this->MomentumConstraintFlag == false)
    {
      return  new PeriodicMomentumQuantumNumber (this->GetMomentumValue(index), this->MaxMomentum);
    }
  else
    return new PeriodicMomentumQuantumNumber (this->MomentumConstraint, this->MaxMomentum);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnTorusShort::ExtractSubspace (AbstractQuantumNumber& q, 
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

int BosonOnTorusShort::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int CurrentLzMax = this->StateMaxMomentum[index];
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }    
  this->FermionToBoson(this->StateDescription[index], CurrentLzMax + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateMaxMomentum);
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax) || (this->TemporaryState[n1] == 0) || (this->TemporaryState[n2] == 0) || ((n1 == n2) && (this->TemporaryState[n1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = CurrentLzMax;
  if (NewLzMax < m1)
    NewLzMax = m1;
  if (NewLzMax < m2)
    NewLzMax = m2;
  for (int i = CurrentLzMax + 1; i <= NewLzMax; ++i)
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
  while (this->TemporaryState[NewLzMax] == 0x0ul)
    --NewLzMax;
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateMaxMomentum), NewLzMax + this->NbrBosons - 1);
}

// return matrix representation of the annihilation operator a_i
//
// i = operator index
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& BosonOnTorusShort::A (int i, Matrix& M)
{
  return M;
}

// return matrix representation of the creation operator a^+_i
//
// i = operator index
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& BosonOnTorusShort::Ad (int i, Matrix& M)
{
  return M;
}

// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

int BosonOnTorusShort::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    return PosMin;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnTorusShort::PrintState (ostream& Str, int state)
{
  int Max = this->StateMaxMomentum[state];
  this->FermionToBoson(this->StateDescription[state], Max + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateMaxMomentum);
  int i = 0;
  for (; i <= Max; ++i)
    Str << this->TemporaryState[i] << " ";
  for (; i < this->MaxMomentum; ++i)
    Str << "0 ";
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int BosonOnTorusShort::GenerateStates(int nbrBosons, int maxMomentum, int currentMaxMomentum, int pos)
{
  if ((nbrBosons == 0) || (currentMaxMomentum < 0))
    {
      return pos;
    }
  if (currentMaxMomentum == 0)
    {
      this->StateDescription[pos] = 0ul;      
      this->StateDescription[pos] = ((1ul << (nbrBosons + 1)) - 1ul);
      this->StateMaxMomentum[pos] = maxMomentum;
      return pos + 1;
    }

  this->StateDescription[pos] = 0ul;
  this->StateMaxMomentum[pos] = maxMomentum;
  ++pos;

  int TmpNbrBosons = 1;
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  int TmpPos = pos;
  while (TmpNbrBosons < nbrBosons)
    {
      unsigned long Mask = ((0x1ul << (nbrBosons - TmpNbrBosons)) - 0x1ul) << (currentMaxMomentum + TmpNbrBosons - 1);
      TmpPos = this->GenerateStates(TmpNbrBosons, maxMomentum, ReducedCurrentMaxMomentum, pos);
      for (; pos <TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
      ++TmpNbrBosons;
    }
  if (maxMomentum == currentMaxMomentum)
    return this->GenerateStates(nbrBosons, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, pos);
  else
    return this->GenerateStates(nbrBosons, maxMomentum, ReducedCurrentMaxMomentum, pos);
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson in the state
// currentMaxMomentum = momentum maximum value for bosons that are still to be placed
// pos = position in StateDescription array where to store states
// currentMomentum = current value of the momentum
// return value = position from which new states have to be stored

int BosonOnTorusShort::GenerateStates(int nbrBosons, int maxMomentum, int currentMaxMomentum, int pos, int currentMomentum)
{
  if (nbrBosons == 0)
    {
      if ((currentMomentum % this->MaxMomentum) == this->MomentumConstraint)
	{
	  this->StateDescription[pos] = 0x0ul;
	  this->StateMaxMomentum[pos] = maxMomentum;
	  return pos + 1;
	}
      else
	{
	  return pos;
	}
    }
  if (currentMaxMomentum == 0)
    {
      if ((currentMomentum % this->MaxMomentum) == this->MomentumConstraint)
	{
	  this->StateDescription[pos] = ((0x1ul << nbrBosons) - 0x1ul) << currentMaxMomentum;
	  this->StateMaxMomentum[pos] = maxMomentum;
	  return pos + 1;
	}
      else
	return pos;
    }

  int TmpNbrBosons = 0;
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  int TmpPos = pos;
  while (TmpNbrBosons < nbrBosons)
    {
      TmpPos = this->GenerateStates(TmpNbrBosons, maxMomentum, ReducedCurrentMaxMomentum, pos, currentMomentum + (nbrBosons - TmpNbrBosons) * currentMaxMomentum);
      unsigned long Mask = ((0x1ul << (nbrBosons - TmpNbrBosons)) - 0x1ul) << (currentMaxMomentum + TmpNbrBosons - 1);
      for (; pos <TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
      ++TmpNbrBosons;
    }
  if (maxMomentum == currentMaxMomentum)
    return this->GenerateStates(nbrBosons, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, pos, currentMomentum);
  else
    return this->GenerateStates(nbrBosons, maxMomentum, ReducedCurrentMaxMomentum, pos, currentMomentum);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnTorusShort::GenerateLookUpTable(int memory)
{
  // evaluate look-up table size
  this->NbrLzValue =  this->MaxMomentum + this->NbrBosons;
  memory /= (sizeof(int*) * this->NbrLzValue);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrLzValue)
    this->MaximumLookUpShift = this->NbrLzValue;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrLzValue];
  this->LookUpTableShift = new int [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];
  int CurrentLzMax = this->MaxMomentum + this->NbrBosons - 1;
  while ((this->StateDescription[0] >> CurrentLzMax) == 0x0ul)
    --CurrentLzMax;
  int* TmpLookUpTable = this->LookUpTable[CurrentLzMax];
  if (CurrentLzMax < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentLzMax] = 0;
  else
    this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentLzMax];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = this->StateDescription[i];
      int TmpLzMax = CurrentLzMax;
      while ((TmpState >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
     if (CurrentLzMax != TmpLzMax)
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
 	  CurrentLzMax = TmpLzMax;
	  TmpLookUpTable = this->LookUpTable[CurrentLzMax];
	  if (CurrentLzMax < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentLzMax] = 0;
	  else
	    this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentLzMax];
	  TmpLookUpTableValue = TmpState >> CurrentShift;
	  CurrentLookUpTableValue = this->LookUpTableMemorySize;
	  while (CurrentLookUpTableValue > TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = TmpState >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  while (CurrentLookUpTableValue > 0)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// return value = Hilbert space dimension

int BosonOnTorusShort::EvaluateHilbertSpaceDimension(int nbrBosons, int maxMomentum)
{
  FactorialCoefficient Dimension; 
  Dimension.PartialFactorialMultiply(maxMomentum, maxMomentum + nbrBosons - 1); 
  Dimension.FactorialDivide(nbrBosons);
  return Dimension.GetIntegerValue();
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson in the state
// currentMaxMomentum = momentum maximum value for bosons that are still to be placed
// currentMomentum = current value of the momentum
// return value = Hilbert space dimension

long BosonOnTorusShort::EvaluateHilbertSpaceDimension(int nbrBosons, int maxMomentum, int currentMaxMomentum, int currentMomentum)
{
  if (nbrBosons == 0)
    {
      if ((currentMomentum % this->MaxMomentum) == this->MomentumConstraint)
	{
	  return 1l;
	}
      else
	{
	  return 0l;
	}
    }
  if (currentMaxMomentum == 0)
    {
      if ((currentMomentum % this->MaxMomentum) == this->MomentumConstraint)
	{
	  return 1l;
	}
      else
	return 0l;
    }

  int TmpNbrBosons = 0;
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  long TmpNbrStates = 0l;
  while (TmpNbrBosons < nbrBosons)
    {
     TmpNbrStates += this->GenerateStates(TmpNbrBosons, maxMomentum, ReducedCurrentMaxMomentum, currentMomentum + (nbrBosons - TmpNbrBosons) * currentMaxMomentum);
     ++TmpNbrBosons;
    }
  if (maxMomentum == currentMaxMomentum)
    return TmpNbrStates + this->GenerateStates(nbrBosons, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, currentMomentum);
  else
    return TmpNbrStates + this->GenerateStates(nbrBosons, maxMomentum, ReducedCurrentMaxMomentum, currentMomentum);
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// kySector = Ky sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnTorusShort::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int kySector, RealVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((kySector == 0) && (nbrBosonSector == 0))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  if (subsytemSize > this->MaxMomentum)
    {
      if ((kySector == this->MomentumConstraint) && (nbrBosonSector == this->NbrBosons))
	{
	  RealSymmetricMatrix TmpDensityMatrix(this->HilbertSpaceDimension);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    for (int j = i; j < this->HilbertSpaceDimension; ++j)
	      TmpDensityMatrix.SetMatrixElement(i, j, groundState[i] * groundState[j]);
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  long TmpNbrNonZeroElements = 0;
  BosonOnTorusShort TmpDestinationHilbertSpace(nbrBosonSector, subsytemSize - 1, kySector);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  BosonOnTorusShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, this->MaxMomentum - subsytemSize, 0);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      //      for (int i = 0; i <= TmpHilbertSpace.StateMaxMomentum[MinIndex]; ++i)
      //	this->TemporaryState[i] = TmpHilbertSpace.StateDescription[MinIndex][i];
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{

// 	  unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
// 	  int TmpLzMax = this->FermionBasis->LzMax;
// 	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
// 	    --TmpLzMax;
// 	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
// 	  if (TmpPos != this->HilbertSpaceDimension)
// 	    {
// 	      TmpStatePosition[Pos] = TmpPos;
// 	      TmpStatePosition2[Pos] = j;
// 	      ++Pos;
// 	    }
	}

      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]]);
	    }
	}
    }
  
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}
