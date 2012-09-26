////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                            class of spin 1 chain                           //
//                                                                            //
//                        last modification : 04/04/2001                      //
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


#include "HilbertSpace/Potts3Chain.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include <iostream>


using std::cout;
using std::endl;


#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880
#endif


// default constructor
//

Potts3Chain::Potts3Chain () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableSize = 0;
  this->LookUpTableMask = 0;
  this->LookUpPosition = 0;
  this->HilbertSpaceDimension = 0;
  this->ChainDescription = 0;
  this->ChainLength = 0;
  this->Sz = 0;
  this->FixedQuantumNumberFlag = false;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// memorySize = memory size in bytes allowed for look-up table

Potts3Chain::Potts3Chain (int chainLength, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedQuantumNumberFlag = false;
  this->LargeHilbertSpaceDimension = 3l;
  for (int i = 1; i < chainLength; i++)
    this->LargeHilbertSpaceDimension *= 3l;
  

  this->LookUpPosition = 0;
  this->LookUpTableSize = 4;
  memorySize >>= 2;
  this->LookUpTableMask = 0xfffffffc;
  while ((this->LookUpPosition <= this->ChainLength) && (memorySize >=  4))
    {
      this->LookUpTableMask <<= 2;
      this->LookUpTableSize <<= 2;
      memorySize >>= 2;
      this->LookUpPosition++;
    }
  this->LookUpTableMask = ~this->LookUpTableMask;
  this->LookUpTable = new int [this->LookUpTableSize];

  this->LargeHilbertSpaceDimension = this->GenerateStates (this->ChainLength - 1, this->Sz, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
}

// constructor for complete Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table

Potts3Chain::Potts3Chain (int chainLength, int sz, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Sz = sz % 3;
  this->FixedQuantumNumberFlag = true;

  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->ChainLength - 1, this->Sz);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->LookUpPosition = 0;
  this->LookUpTableSize = 4;
  memorySize >>= 1;
  this->LookUpTableMask = 0xfffffffc;
  while ((this->LookUpPosition < this->ChainLength) && (memorySize >=  4))
    {
      this->LookUpTableMask <<= 2;
      this->LookUpTableSize <<= 2;
      memorySize >>= 2;
      this->LookUpPosition++;
    }
  this->LookUpTableMask = ~this->LookUpTableMask;
  this->LookUpTable = new int [this->LookUpTableSize];

  this->ChainDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->LargeHilbertSpaceDimension = this->GenerateStates (this->ChainLength - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
}

// constructor from pre-constructed datas
//
// largehilbertSpaceDimension = Hilbert space dimension
// chainDescription = array describing states
// chainLength = number of spin 1
// sz = twice the value of total Sz component
// fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
// lookUpTable = look-up table
// lookUpTableSize = look-Up table size
// lookUpTablePosition = last position described by the look-Up table
// lookUpTableMask = look-Up table mask

Potts3Chain::Potts3Chain (long largeHilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
			  int sz, bool fixedQuantumNumberFlag, int* lookUpTable, int lookUpTableSize, 
			  int lookUpPosition, unsigned long lookUpTableMask)
{
  this->Flag.Initialize();
  this->LookUpTable = lookUpTable;
  this->LookUpTableMask = lookUpTableMask;
  this->LookUpPosition = lookUpPosition;
  this->LookUpTableSize = lookUpTableSize;
  this->ChainDescription = chainDescription;
  this->Sz = sz;
  this->FixedQuantumNumberFlag = fixedQuantumNumberFlag;
  this->ChainLength = chainLength;
  this->LargeHilbertSpaceDimension = largeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
}
  
// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Potts3Chain::Potts3Chain (const Potts3Chain& chain) 
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpPosition = chain.LookUpPosition;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->ChainDescription = chain.ChainDescription;
      this->Sz = chain.Sz;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableSize = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
      this->LargeHilbertSpaceDimension = 0l;
    }
}

// destructor
//

Potts3Chain::~Potts3Chain () 
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->ChainDescription;
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Potts3Chain& Potts3Chain::operator = (const Potts3Chain& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->ChainDescription;
      delete[] this->LookUpTable;
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainDescription = chain.ChainDescription;
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->LookUpPosition = chain.LookUpPosition;
      this->Sz = chain.Sz;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
      this->LargeHilbertSpaceDimension = chain.LargeHilbertSpaceDimension;
   }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableSize = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
      this->LargeHilbertSpaceDimension = 0l;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Potts3Chain::Clone()
{
  return new Potts3Chain (*this);
}

// evaluate Hilbert space dimension
//
// currentSite = current site to occupy
// currentSzValue = state current Sz value 
// return value = Hilbert space dimension

long Potts3Chain::EvaluateHilbertSpaceDimension(int currentSite, int currentSzValue)
{
  if (currentSite < 0)
    {
      if ((currentSzValue %3) == this->Sz)
	return 1l;
      else
	return 0l;
    }
  long TmpDimension = this->EvaluateHilbertSpaceDimension(currentSite - 1, currentSzValue + 2);
  TmpDimension += this->EvaluateHilbertSpaceDimension(currentSite - 1, currentSzValue + 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(currentSite - 1, currentSzValue);
  return TmpDimension;
}


// generate all states with no constraint on total Sz
//
// currentSite = current site to occupy
// currentPosition = current position of the state that has to be considered
// return value = number of generated states

long Potts3Chain::GenerateStates(int currentSite, long currentPosition) 
{
  if (currentSite < 0)
    {
      this->ChainDescription[currentPosition] = 0x0l;
      return (currentPosition + 1l);
    }
  long TmpPosition = this->GenerateStates(currentSite - 1, currentPosition);
  unsigned long TmpMask = 0x2ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    this->ChainDescription[currentPosition] |= TmpMask;
  TmpPosition = this->GenerateStates(currentSite - 1, currentPosition);
  TmpMask = 0x1ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    this->ChainDescription[currentPosition] |= TmpMask;
  return this->GenerateStates(currentSite - 1, currentPosition);
}

// generate all states corresponding to a given total Sz
//
// currentSite = current site to occupy
// currentSzValue = state current Sz value 
// currentPosition = current position of the state that has to be considered
// return value = number of generated states

long Potts3Chain::GenerateStates(int currentSite, int currentSzValue, long currentPosition) 
{
  if (currentSite < 0)
    {
      if ((currentSzValue %3) == this->Sz)
	{
	  this->ChainDescription[currentPosition] = 0x0l;
	  return (currentPosition + 1l);
	}
      else
	return currentPosition;
    }
  long TmpPosition = this->GenerateStates(currentSite - 1, currentSzValue + 2, currentPosition);
  unsigned long TmpMask = 0x2ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    {
      this->ChainDescription[currentPosition] |= TmpMask;
    }
  TmpPosition = this->GenerateStates(currentSite - 1, currentSzValue + 1, currentPosition);
  TmpMask = 0x1ul << (currentSite << 1);
  for (; currentPosition < TmpPosition; ++currentPosition)
    {
      this->ChainDescription[currentPosition] |= TmpMask;
    }
  return this->GenerateStates(currentSite - 1, currentSzValue, currentPosition);
}


// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Potts3Chain::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->FixedQuantumNumberFlag == true)
    {
      TmpList += new SzQuantumNumber (this->Sz);
    }
  else
    {
      for (int i = 0; i < 3; i++)
	{
	  TmpList += new SzQuantumNumber (i);
	}
    }
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* Potts3Chain::GetQuantumNumber (int index)
{ 
  return new SzQuantumNumber (this->TotalSz(index));
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Potts3Chain::TotalSz (int index)
{
  if (this->FixedQuantumNumberFlag == true)
    return this->Sz;
  unsigned long State = this->ChainDescription[index];
  unsigned long TmpSz = 0l;
  for (int i = 0; i < this->ChainLength; ++i)
    {
      TmpSz += (State & 0x3ul);
      State >>= 2;
    }
  return (((int) TmpSz) % 3);
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Potts3Chain::SmiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long tmpState2 = tmpState;
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << j);
  switch (tmpState2)
    {
    case 0x1ul:
      tmpState |= 0x2ul << j;
      break;
    case 0x0ul:
      tmpState |= 0x1ul << j;
      break;
    }	  
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x2ul:
      tmpState |= 0x1ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x2ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Potts3Chain::SpiSpj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long tmpState2 = tmpState;
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << j);
  switch (tmpState2)
    {
    case 0x1ul:
      tmpState |= 0x2ul << j;
      break;
    case 0x0ul:
      tmpState |= 0x1ul << j;
      break;
    }	  
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x1ul:
      tmpState |= 0x2ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x1ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Potts3Chain::SmiSmj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long tmpState2 = tmpState;
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << j);
  switch (tmpState2)
    {
    case 0x2ul:
      tmpState |= 0x1ul << j;
      break;
    case 0x0ul:
      tmpState |= 0x2ul << j;
      break;
    }	  
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x2ul:
      tmpState |= 0x1ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x2ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Potts3Chain::SpiSzj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long tmpState2 = tmpState;
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  if (tmpState2 == 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = (double) tmpState2;
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x1ul:
      tmpState |= 0x2ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x1ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Potts3Chain::SmiSzj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long tmpState2 = tmpState;
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  if (tmpState2 == 0x0ul)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = (double) tmpState2;
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x2ul:
      tmpState |= 0x1ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x2ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Potts3Chain::Spi (int i, int state, double& coefficient)
{
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x1ul:
      tmpState |= 0x2ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x1ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Potts3Chain::Smi (int i, int state, double& coefficient)
{
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  tmpState2 = tmpState;
  tmpState2 >>= i;
  tmpState2 &= 0x3ul;
  tmpState &= ~(0x3ul << i);
  switch (tmpState2)
    {
    case 0x2ul:
      tmpState |= 0x1ul << i;      
      break;
    case 0x0ul:
      tmpState |= 0x2ul << i;      
      break;
    }	  
  return this->FindStateIndex(tmpState);
}

// translate a state assuming the system have periodic boundary conditions (increasing the site index)
//
// nbrTranslations = number of translations to apply
// state = index of the state to translate 
// return value = index of resulting state

int Potts3Chain::TranslateState (int nbrTranslations, int state)
{
  unsigned long TmpState = this->ChainDescription[state];
  TmpState = (((TmpState & (0x3ul << ((this->ChainLength - nbrTranslations) - 1ul)) << 1) << nbrTranslations)
	      | (TmpState >> ((this->ChainLength - nbrTranslations) << 1)));
  return this->FindStateIndex(TmpState);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* Potts3Chain::ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter)
{
  if (q.GetQuantumNumberType() != AbstractQuantumNumber::Sz)
    return new Potts3Chain();
  int TmpSz = ((SzQuantumNumber&) q).GetSz();
  long LargeHilbertSubspaceDimension = 0l;
  int* TmpConvArray = new int [this->LargeHilbertSpaceDimension];
  for (int i = 0; i < this->LargeHilbertSpaceDimension; i++)
    {
      if (this->TotalSz(i) == TmpSz)
	{
	  TmpConvArray[LargeHilbertSubspaceDimension] = i;
	  LargeHilbertSubspaceDimension++;	  
	}
    }
  int* ConvArray = new int [LargeHilbertSubspaceDimension];
  unsigned long* SubspaceDescription = new unsigned long [LargeHilbertSubspaceDimension];
  int* SubspaceLookUpTable = new int [this->LookUpTableSize];
  unsigned long TestMask = this->ChainDescription[TmpConvArray[0]] & this->LookUpTableMask;
  SubspaceLookUpTable[TestMask] = 0;
  SubspaceDescription[0] = this->ChainDescription[TmpConvArray[0]];
  ConvArray[0] = TmpConvArray[0];
  for (long i = 1; i < LargeHilbertSubspaceDimension; i++)
    {
      if ((this->ChainDescription[TmpConvArray[i]] & this->LookUpTableMask) != TestMask)
	{
	  TestMask = this->ChainDescription[TmpConvArray[i]] & this->LookUpTableMask;
	  SubspaceLookUpTable[TestMask] = i;
	}
      SubspaceDescription[i] = this->ChainDescription[TmpConvArray[i]];
      ConvArray[i] = TmpConvArray[i];
    }
  converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, (int) LargeHilbertSubspaceDimension, ConvArray);
  return (AbstractSpinChain*) new Potts3Chain (LargeHilbertSubspaceDimension, SubspaceDescription, this->ChainLength,
					       TmpSz, true, SubspaceLookUpTable, this->LookUpTableSize, 
					       this->LookUpPosition, this->LookUpTableMask);
}

// find state index
//
// state = state description
// return value = corresponding index

int Potts3Chain::FindStateIndex(unsigned long state)
{
  int index = 0;//this->LookUpTable[state & this->LookUpTableMask];
  unsigned long* tmpState = &(this->ChainDescription[index]);
  while ((index < this->HilbertSpaceDimension) && (state != *(tmpState++)))
    ++index;
  return index;   
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Potts3Chain::PrintState (ostream& Str, int state)
{
  unsigned long StateDescription = this->ChainDescription[state];  
  for (int j = 0; j < this->ChainLength; ++j)
    {
      Str << (StateDescription & 0x3ul) << " ";
      StateDescription >>= 2;
    }
  return Str;
}
