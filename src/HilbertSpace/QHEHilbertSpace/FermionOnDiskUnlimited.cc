////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of fermions on disk with no restriction on the         //
//               number of reachable states or the number of fermions         //
//                                                                            //
//                        last modification : 03/03/2004                      //
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
#include "HilbertSpace/QHEHilbertSpace/FermionOnDiskUnlimited.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"

#include <math.h>


using std::cout;
using std::endl;
using std::hex;
using std::dec;


// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = momentum total value

FermionOnDiskUnlimited::FermionOnDiskUnlimited (int nbrFermions, int totalLz)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = this->TotalLz - (((this->NbrFermions - 1) * (this->NbrFermions - 2)) / 2);
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz);
  this->Flag.Initialize();
  this->StateDescription = new FermionOnSphereLongState [this->HilbertSpaceDimension];
  this->TemporaryStateReducedNbrState = FermionOnSphereLongStateGetReducedNbrState(this->NbrLzValue);
  this->TemporaryState.Resize(this->TemporaryStateReducedNbrState);
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  this->ReducedNbrState = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrFermions, this->LzMax, this->LzMax, this->TotalLz, 0);
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    UsedMemory += ((this->ReducedNbrState[i] + 1) * sizeof(unsigned long));
  UsedMemory += this->HilbertSpaceDimension * (3 * sizeof(int));
  UsedMemory += this->NbrLzValue * sizeof(int);
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      if (this->NbrStateInLookUpTable[i] != 0)
	{
	  for (unsigned long j = 0; j <= this->HashKeyMask; ++j)
	    if (this->NbrStateInLookUpTable[i][j] > 0)      
	      UsedMemory += this->NbrStateInLookUpTable[i][j] * sizeof(int);
	}
    }
  UsedMemory += this->NbrLzValue * sizeof(int*);
  UsedMemory +=  (1 << this->MaximumSignLookUp) * sizeof(double);
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
// fermions = reference on the hilbert space to copy to copy

FermionOnDiskUnlimited::FermionOnDiskUnlimited(const FermionOnDiskUnlimited& fermions)
{
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  this->LookUpTable = fermions.LookUpTable;
  this->NbrStateInLookUpTable = fermions.NbrStateInLookUpTable;
  this->HashKeyMask = fermions.HashKeyMask;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->TemporaryStateReducedNbrState = FermionOnSphereLongStateGetReducedNbrState(this->NbrLzValue);
  this->TemporaryState.Resize(this->TemporaryStateReducedNbrState);
}

// destructor
//

FermionOnDiskUnlimited::~FermionOnDiskUnlimited ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->ReducedNbrState;
      delete[] this->SignLookUpTable;
      for (int i = 0; i < this->NbrLzValue; ++i)
	{
	  if (this->NbrStateInLookUpTable[i] != 0)
	    {
	      for (unsigned long j = 0; j <= this->HashKeyMask; ++j)
		if (this->NbrStateInLookUpTable[i][j] > 0)
		  delete[] this->LookUpTable[i][j];
	      delete[] this->LookUpTable[i];
	      delete[] this->NbrStateInLookUpTable[i];
	    }
	}
      delete[] this->LookUpTable;
      delete[] this->NbrStateInLookUpTable;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnDiskUnlimited& FermionOnDiskUnlimited::operator = (const FermionOnDiskUnlimited& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->ReducedNbrState;
      for (int i = 0; i < this->NbrLzValue; ++i)
	{
	  if (this->NbrStateInLookUpTable[i] != 0)
	    {
	      for (unsigned long j = 0; j <= this->HashKeyMask; ++j)
		if (this->NbrStateInLookUpTable[i][j] > 0)
		  delete[] this->LookUpTable[i][j];
	      delete[] this->LookUpTable[i];
	      delete[] this->NbrStateInLookUpTable[i];
	    }
	}
      delete[] this->LookUpTable;
      delete[] this->NbrStateInLookUpTable;
   }
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  this->LookUpTable = fermions.LookUpTable;
  this->NbrStateInLookUpTable = fermions.NbrStateInLookUpTable;
  this->HashKeyMask = fermions.HashKeyMask;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->TemporaryStateReducedNbrState = FermionOnSphereLongStateGetReducedNbrState(this->NbrLzValue);
  this->TemporaryState.Resize(this->TemporaryStateReducedNbrState);
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnDiskUnlimited::Clone()
{
  return new FermionOnDiskUnlimited(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnDiskUnlimited::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnDiskUnlimited::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnDiskUnlimited::ExtractSubspace (AbstractQuantumNumber& q, 
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

int FermionOnDiskUnlimited::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateLzMax = this->StateLzMax[index];
  if ((n1 > StateLzMax) || (n2 > StateLzMax) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
//  cout << "-----------" << m1 << " " << m2 << " " << n1 << " " << n2 << endl;
  int TmpReducedNbrState = this->ReducedNbrState[index];
  this->TemporaryState.EmptyState(this->TemporaryStateReducedNbrState);
  this->TemporaryState.Assign(this->StateDescription[index], TmpReducedNbrState);
  if ((this->TemporaryState.GetOccupation(n2) == 0) || (this->TemporaryState.GetOccupation(n1) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  this->TemporaryState.GetPermutationSign(n2, TmpReducedNbrState, this->SignLookUpTable, coefficient);
  this->TemporaryState.DecrementOccupation(n2);
  this->TemporaryState.PrintState(cout, TmpReducedNbrState, FermionOnSphereLongStateGetRemainderNbrState(StateLzMax + 1)) << endl;
  cout << coefficient << endl;
  this->TemporaryState.GetPermutationSign(n1, TmpReducedNbrState, this->SignLookUpTable, coefficient);
  this->TemporaryState.DecrementOccupation(n1);
  this->TemporaryState.PrintState(cout, TmpReducedNbrState, FermionOnSphereLongStateGetRemainderNbrState(StateLzMax + 1)) << endl;
  cout << coefficient << endl;

  StateLzMax = this->TemporaryState.GetHighestIndex(TmpReducedNbrState);
  TmpReducedNbrState = FermionOnSphereLongStateGetReducedNbrState(StateLzMax + 1);

  if (this->TemporaryState.GetOccupation(m2) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > StateLzMax)
    {
      StateLzMax = m2;
      TmpReducedNbrState = FermionOnSphereLongStateGetReducedNbrState(StateLzMax + 1);
   }  
  else
    {
      this->TemporaryState.GetPermutationSign(m2, TmpReducedNbrState, this->SignLookUpTable, coefficient);
    }
  this->TemporaryState.PrintState(cout, TmpReducedNbrState, FermionOnSphereLongStateGetRemainderNbrState(StateLzMax + 1)) << endl;
  cout << coefficient << endl;
  this->TemporaryState.IncrementOccupation(m2);
  if (this->TemporaryState.GetOccupation(m1) != 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > StateLzMax)
    {
      StateLzMax = m1;
      TmpReducedNbrState = FermionOnSphereLongStateGetReducedNbrState(StateLzMax + 1);
    }  
  else
    {
      this->TemporaryState.GetPermutationSign(m1, TmpReducedNbrState, this->SignLookUpTable, coefficient);
    }
  this->TemporaryState.IncrementOccupation(m1);
  this->TemporaryState.PrintState(cout, TmpReducedNbrState, FermionOnSphereLongStateGetRemainderNbrState(StateLzMax + 1)) << endl;
  cout << coefficient << endl;
/*  this->StateDescription[index].PrintState(cout,  this->ReducedNbrState[index], FermionOnSphereLongStateGetRemainderNbrState(this->StateLzMax[index] + 1)) << endl;
  this->TemporaryState.PrintState(cout, TmpReducedNbrState, FermionOnSphereLongStateGetRemainderNbrState(StateLzMax + 1)) << endl;
  cout << TmpReducedNbrState << " " << StateLzMax << endl;*/
  return this->FindStateIndex(this->TemporaryState, TmpReducedNbrState, StateLzMax);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnDiskUnlimited::AdA (int index, int m)
{
  if (this->StateDescription[index].GetOccupation(m) == 0)
    return 0.0;
  else
    return 1.0;
}

// find state index
//
// stateDescription = reference on the state description
// reducedNbrState =reduced number of state (aka the number of unsigned long per state) minus 1
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnDiskUnlimited::FindStateIndex(FermionOnSphereLongState& stateDescription, int reducedNbrState, int lzmax)
{
  int TmpKey = stateDescription.GetHashKey(reducedNbrState, this->HashKeyMask);
  int Max = this->NbrStateInLookUpTable[lzmax][TmpKey] - 1;
  int* TmpLookUpTable = this->LookUpTable[lzmax][TmpKey];
  int Min = 0;
  int Pos = 0;
  while ((Max - Min) > 1)
    {
      Pos = (Max + Min) >> 1;
      if (stateDescription.Greater (this->StateDescription[TmpLookUpTable[Pos]], reducedNbrState))
	Max = Pos;
      else
	Min = Pos;
    }
//  cout << TmpKey << endl;
//  cout << Min << " " << Max << " " << TmpLookUpTable[Min] << endl;
  if (stateDescription.Different(this->StateDescription[TmpLookUpTable[Min]], reducedNbrState))
    return TmpLookUpTable[Max];
  else
    return TmpLookUpTable[Min];
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnDiskUnlimited::PrintState (ostream& Str, int state)
{
  this->StateDescription[state].PrintState(Str, this->ReducedNbrState[state], FermionOnSphereLongStateGetRemainderNbrState(this->StateLzMax[state] + 1));
  Str << " position = " << this->FindStateIndex(this->StateDescription[state], this->ReducedNbrState[state], this->StateLzMax[state])  << "   ";
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// currentLzMax = momentum maximum value for fermions that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int FermionOnDiskUnlimited::GenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int pos)
{
  if ((nbrFermions == 0) || (totalLz < 0) || (currentLzMax < (nbrFermions - 1)))
    return pos;
  int LzTotalMax = ((2 * currentLzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (totalLz > LzTotalMax)
    return pos;
  if ((nbrFermions == 1) && (currentLzMax >= totalLz))
    {
      this->ReducedNbrState[pos] = FermionOnSphereLongStateGetReducedNbrState(lzMax + 1);
      this->StateDescription[pos].Resize(this->ReducedNbrState[pos]);
      this->StateLzMax[pos] = lzMax;
      this->StateDescription[pos].SetOccupation(totalLz);
      return pos + 1;
    }
  if (LzTotalMax == totalLz)
    {
      this->ReducedNbrState[pos] = FermionOnSphereLongStateGetReducedNbrState(lzMax + 1);
      this->StateDescription[pos].Resize(this->ReducedNbrState[pos]);
      this->StateLzMax[pos] = lzMax;
      for (int i = currentLzMax - nbrFermions + 1; i <= currentLzMax; ++i)
	this->StateDescription[pos].SetOccupation(i);
      return pos + 1;
    }

  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpPos = this->GenerateStates(nbrFermions - 1, lzMax, ReducedCurrentLzMax, totalLz - currentLzMax, pos);
  for (int i = pos; i < TmpPos; i++)
    this->StateDescription[i].SetOccupation(currentLzMax);
  if (lzMax == currentLzMax)
    return this->GenerateStates(nbrFermions, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, TmpPos);
  else
    return this->GenerateStates(nbrFermions, lzMax, ReducedCurrentLzMax, totalLz, TmpPos);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnDiskUnlimited::GenerateLookUpTable(int memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int) * this->NbrLzValue);
  this->HashKeyMask = (unsigned long) 0x1;
  while (memory > 0)
    {
      memory >>= 1;
      this->HashKeyMask <<= 1;
      this->HashKeyMask |= (unsigned long) 0x1;
    }

//  cout <<  this->HashKeyMask << " " << this->HilbertSpaceDimension << endl;
  // construct  look-up tables for searching states
  this->LookUpTable = new int** [this->NbrLzValue];
  int CurrentLzMax = this->StateLzMax[0];
  this->LookUpTable[CurrentLzMax] = new int* [this->HashKeyMask + (unsigned long) 1];
  this->NbrStateInLookUpTable = new int* [this->NbrLzValue];
  for (int i = 0; i < this->LzMax; ++i)
    this->NbrStateInLookUpTable[i] = 0;
  this->NbrStateInLookUpTable[CurrentLzMax] = new int [this->HashKeyMask + (unsigned long) 1];
  for (unsigned long j = 0; j <= this->HashKeyMask; ++j)
    this->NbrStateInLookUpTable[CurrentLzMax][j] = 0;
  unsigned long* TmpHashTable = new unsigned long [this->HilbertSpaceDimension];
  unsigned long TmpKey = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
//      cout << i << ": ";
//      this->StateDescription[i].PrintState(cout, this->ReducedNbrState[i], FermionOnSphereLongStateGetRemainderNbrState(this->StateLzMax[i] + 1));
//      cout << "(" << this->ReducedNbrState[i] << ")"<< endl;
      if (CurrentLzMax != this->StateLzMax[i])
	{
	  CurrentLzMax = this->StateLzMax[i];
	  this->LookUpTable[CurrentLzMax] = new int* [this->HashKeyMask + (unsigned long) 1];
	  this->NbrStateInLookUpTable[CurrentLzMax] = new int [this->HashKeyMask + (unsigned long) 1];
	  TmpKey = this->StateDescription[i].GetHashKey(this->ReducedNbrState[i], this->HashKeyMask);
	  TmpHashTable[i] = TmpKey;
	  for (unsigned long j = 0; j <= this->HashKeyMask; ++j)
	    this->NbrStateInLookUpTable[CurrentLzMax][j] = 0;
	  ++this->NbrStateInLookUpTable[CurrentLzMax][TmpKey];	  
	}
      else
	{
	  TmpKey = this->StateDescription[i].GetHashKey(this->ReducedNbrState[i], this->HashKeyMask);
	  TmpHashTable[i] = TmpKey;
	  ++this->NbrStateInLookUpTable[CurrentLzMax][TmpKey];
	}
    }
  for (int i = 0; i < this->NbrLzValue; ++i)
    if (this->NbrStateInLookUpTable[i] != 0)
      for (unsigned long j = 0; j <= this->HashKeyMask; ++j)
	if (this->NbrStateInLookUpTable[i][j] != 0)
	  {
	    this->LookUpTable[i][j] = new int [this->NbrStateInLookUpTable[i][j]];
//	    cout << i << " " << j << " " << this->NbrStateInLookUpTable[i][j] << endl;
	    this->NbrStateInLookUpTable[i][j] = 0;
	  }
	else
	  {
	    this->LookUpTable[i][j] = 0;
	  }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
//      cout << i << " " << this->StateLzMax[i] << " " << TmpHashTable[i] << " " << this->NbrStateInLookUpTable[this->StateLzMax[i]][TmpHashTable[i]] << endl;
      this->LookUpTable[this->StateLzMax[i]][TmpHashTable[i]][this->NbrStateInLookUpTable[this->StateLzMax[i]]
							      [TmpHashTable[i]]++] = i;
    }
  delete[] TmpHashTable;
  
  // look-up tables for evaluating sign when applying creation/annihilation operators
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double [Size];
  int Count;
  int TmpNbr;
  for (int j = 0; j < Size; ++j)
    {
      Count = 0;
      TmpNbr = j;
      while (TmpNbr != 0)
	{
	  if (TmpNbr & 0x1)
	    ++Count;
	  TmpNbr >>= 1;
	}
      if (Count & 1)
	this->SignLookUpTable[j] = -1.0;
      else
	this->SignLookUpTable[j] = 1.0;
    }
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

int FermionOnDiskUnlimited::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (lzMax < (nbrFermions - 1)))
    return 0;
  int LzTotalMax = ((2 * lzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (LzTotalMax < totalLz)
    return 0;
  if ((nbrFermions == 1) && (lzMax >= totalLz))
    return 1;
  if (LzTotalMax == totalLz)
    return 1;
  return  (this->EvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax)
	   + this->EvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz));
}
