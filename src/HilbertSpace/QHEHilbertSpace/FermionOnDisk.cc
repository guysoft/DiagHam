////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of fermions on disk with restriction on the            //
//              number of reachable states or the number of fermions          //
//                                                                            //
//                        last modification : 01/02/2004                      //
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
#include "HilbertSpace/QHEHilbertSpace/FermionOnDisk.h"
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

FermionOnDisk::FermionOnDisk (int nbrFermions, int totalLz)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = this->TotalLz - (((this->NbrFermions - 1) * (this->NbrFermions - 2)) / 2);
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz);
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrFermions, this->LzMax, this->LzMax, this->TotalLz, 0);
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
  UsedMemory += this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
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

FermionOnDisk::FermionOnDisk(const FermionOnDisk& fermions)
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
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTableMask = fermions.LookUpTableMask;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
}

// destructor
//

FermionOnDisk::~FermionOnDisk ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->LookUpTableShift;
      delete[] this->LookUpTableMask;
      for (int i = 0; i < this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnDisk& FermionOnDisk::operator = (const FermionOnDisk& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->LookUpTableShift;
      delete[] this->LookUpTableMask;
      for (int i = 0; i < this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
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
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTableMask = fermions.LookUpTableMask;
  this->LookUpTable = fermions.LookUpTable;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnDisk::Clone()
{
  return new FermionOnDisk(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnDisk::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnDisk::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnDisk::ExtractSubspace (AbstractQuantumNumber& q, 
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

int FermionOnDisk::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int StateLzMax = this->StateLzMax[index];
  unsigned long State = this->StateDescription[index];
  if ((n1 > StateLzMax) || (n2 > StateLzMax) || ((State & (((unsigned long) (0x1)) << n1)) == 0) 
      || ((State & (((unsigned long) (0x1)) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = StateLzMax;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) (0x1)) << n2);
  if (NewLzMax == n2)
    while ((TmpState >> NewLzMax) == 0)
      --NewLzMax;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) (0x1)) << n1);
  if (NewLzMax == n1)
    while ((TmpState >> NewLzMax) == 0)
      --NewLzMax;
  if ((TmpState & (((unsigned long) (0x1)) << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewLzMax)
    {
      NewLzMax = m2;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m2);
  if ((TmpState & (((unsigned long) (0x1)) << m1))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewLzMax)
    {
      NewLzMax = m1;
    }
  else
    {      
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnDisk::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  int StateLzMax = this->StateLzMax[index];
  unsigned long State = this->StateDescription[index];
  --nbrIndices;
  for (int i = 0; i < nbrIndices; ++i)
    {
      if ((n[i] > StateLzMax) || ((State & (((unsigned long) (0x1)) << n[i])) == 0))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      for (int j = i + 1; j <= nbrIndices; ++j)
	if ((n[i] == n[j]) || (m[i] == m[j]))
	  {
	    coefficient = 0.0;
	    return this->HilbertSpaceDimension; 	    
	  }
    }
  int NewLzMax = StateLzMax;
  unsigned long TmpState = State;
  int Index;
  coefficient = 1.0;
  for (int i = nbrIndices; i >= 0; --i)
    {
      Index = n[i];
      coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
      coefficient *= this->SignLookUpTable[(TmpState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      TmpState &= ~(((unsigned long) (0x1)) << Index);
      if (NewLzMax == Index)
	while ((TmpState >> NewLzMax) == 0)
	  --NewLzMax;
    }
  for (int i = nbrIndices; i >= 0; --i)
    {
      Index = m[i];
      if ((TmpState & (((unsigned long) (0x1)) << Index))!= 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
	}
      TmpState |= (((unsigned long) (0x1)) << Index);
    }
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnDisk::AdA (int index, int m)
{
  if ((this->StateLzMax[index] & (((unsigned long) (0x1)) << m)) == 0)
    return 0.0;
  else
    return 1.0;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnDisk::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  int Pos = this->LookUpTable[lzmax][(stateDescription >> this->LookUpTableShift[lzmax]) & this->LookUpTableMask[lzmax]];
  while (this->StateDescription[Pos] != stateDescription)
    ++Pos;
  return Pos;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnDisk::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  for (int i = 0; i < this->NbrLzValue; ++i)
    Str << ((TmpState >> i) & ((unsigned long) 0x1)) << " ";
  Str << " position = " << this->FindStateIndex(TmpState, this->StateLzMax[state]) << "   " << TmpState;
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

int FermionOnDisk::GenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int pos)
{
  if ((nbrFermions == 0) || (totalLz < 0) || (currentLzMax < (nbrFermions - 1)))
    return pos;
  int LzTotalMax = ((2 * currentLzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (totalLz > LzTotalMax)
    return pos;
  if ((nbrFermions == 1) && (currentLzMax >= totalLz))
    {
      this->StateDescription[pos] = ((unsigned long) 0x1) << totalLz;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }
  if (LzTotalMax == totalLz)
    {
      unsigned long Mask = 0;
      for (int i = currentLzMax - nbrFermions + 1; i <= currentLzMax; ++i)
	Mask |= (((unsigned long) 1) << i);
      this->StateDescription[pos] = Mask;
      this->StateLzMax[pos] = lzMax;
      return pos + 1;
    }

  int ReducedCurrentLzMax = currentLzMax - 1;
  int TmpPos = this->GenerateStates(nbrFermions - 1, lzMax, ReducedCurrentLzMax, totalLz - currentLzMax, pos);
  unsigned long Mask = ((unsigned long) 1) << currentLzMax;
  for (int i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (lzMax == currentLzMax)
    return this->GenerateStates(nbrFermions, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, TmpPos);
  else
    return this->GenerateStates(nbrFermions, lzMax, ReducedCurrentLzMax, totalLz, TmpPos);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnDisk::GenerateLookUpTable(int memory)
{
  // evaluate look-up table size
  memory /= (4 * this->NbrLzValue);
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
  this->LookUpTableMask = new unsigned long [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize];
  int CurrentLzMax = this->StateLzMax[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentLzMax];
  if (CurrentLzMax < this->MaximumLookUpShift)
    {
      this->LookUpTableShift[CurrentLzMax] = 0;
#ifdef __64_BITS__
      this->LookUpTableMask[CurrentLzMax] = 0xffffffffffffffff;
#else
      this->LookUpTableMask[CurrentLzMax] = 0xffffffff;
#endif
    }
  else
    {
      this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
#ifdef __64_BITS__
      this->LookUpTableMask[CurrentLzMax] = 0x7fffffffffffffff >> (this->LookUpTableShift[CurrentLzMax] -1);
#else
      this->LookUpTableMask[CurrentLzMax] = 0x7fffffff >> (this->LookUpTableShift[CurrentLzMax] - 1);
#endif
    }
  int CurrentShift = this->LookUpTableShift[CurrentLzMax];
  unsigned long CurrentMask = this->LookUpTableMask[CurrentLzMax];
  unsigned long CurrentLookUpTableValue = (this->StateDescription[0] >> CurrentShift) & CurrentMask;
  unsigned long TmpLookUpTableValue;
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentLzMax != this->StateLzMax[i])
	{
 	  CurrentLzMax = this->StateLzMax[i];
	  TmpLookUpTable = this->LookUpTable[CurrentLzMax];
	  if (CurrentLzMax < this->MaximumLookUpShift)
	    {
	      this->LookUpTableShift[CurrentLzMax] = 0;
#ifdef __64_BITS__
	      this->LookUpTableMask[CurrentLzMax] = 0xffffffffffffffff;
#else
	      this->LookUpTableMask[CurrentLzMax] = 0xffffffff;
#endif
	    }
	  else
	    {
	      this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
#ifdef __64_BITS__
	      this->LookUpTableMask[CurrentLzMax] = 0x7fffffffffffffff >> (this->LookUpTableShift[CurrentLzMax] -1);
#else
	      this->LookUpTableMask[CurrentLzMax] = 0x7fffffff >> (this->LookUpTableShift[CurrentLzMax] - 1);
#endif
	    }
	  CurrentShift = this->LookUpTableShift[CurrentLzMax];
	  CurrentMask = this->LookUpTableMask[CurrentLzMax];
	  CurrentLookUpTableValue = (this->StateDescription[i] >> CurrentShift) & CurrentMask;
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = (this->StateDescription[i] >> CurrentShift) & CurrentMask;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      CurrentLookUpTableValue = TmpLookUpTableValue;
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  
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
#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#endif
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

int FermionOnDisk::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
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
