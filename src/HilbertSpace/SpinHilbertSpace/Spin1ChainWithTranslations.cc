////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of spin 1 chain with translations                  //
//                                                                            //
//                        last modification : 15/10/2003                      //
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


#include "HilbertSpace/SpinHilbertSpace/Spin1ChainWithTranslations.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "GeneralTools/ArrayTools.h"
#include <iostream>


using std::cout;
using std::endl;


#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880
#endif


// default constructor
//

Spin1ChainWithTranslations::Spin1ChainWithTranslations () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableShift = 0;
  this->HilbertSpaceDimension = 0;
  this->ChainDescription = 0;
  this->ChainLength = 0;
  this->Momentum = 0;
  this->ComplementaryStateShift = 0;
  this->Sz = 0;
  this->FixedQuantumNumberFlag = false;
}

// constructor for Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// memorySize = memory size in bytes allowed for look-up table

Spin1ChainWithTranslations::Spin1ChainWithTranslations (int chainLength, int momentum, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedQuantumNumberFlag = false;
  this->ComplementaryStateShift = (this->ChainLength - 1) << 1;

  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;

  this->HilbertSpaceDimension = this->GenerateStates (300000);
}

// constructor for Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table

Spin1ChainWithTranslations::Spin1ChainWithTranslations (int chainLength, int momentum, int sz, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Sz = sz;
  this->FixedQuantumNumberFlag = true;
  this->ComplementaryStateShift = (this->ChainLength - 1) << 1;

  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;

  this->ChainDescription = new unsigned long [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates (300000, this->Sz);
}

// constructor from pre-constructed datas
//
// hilbertSpaceDimension = Hilbert space dimension
// chainDescription = array describing states
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
// lookUpTable = look-up table
// lookUpTableShift = shift to apply to a state to obtain an index to the look-up table 
// complementaryStateShift = shift to apply to move the spin from one end to the other one

Spin1ChainWithTranslations::Spin1ChainWithTranslations (int hilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
							int momentum, int sz, bool fixedQuantumNumberFlag, long* lookUpTable, int lookUpTableShift, 
							int complementaryStateShift)
{
  this->Flag.Initialize();
  this->ComplementaryStateShift = complementaryStateShift;
  this->LookUpTable = lookUpTable;
  this->LookUpTableShift = lookUpTableShift;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->ChainDescription = chainDescription;
  this->Momentum = momentum;
  this->Sz = sz;
  this->FixedQuantumNumberFlag = fixedQuantumNumberFlag;
  this->ChainLength = chainLength;
}
  
// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1ChainWithTranslations::Spin1ChainWithTranslations (const Spin1ChainWithTranslations& chain) 
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->ChainDescription = chain.ChainDescription;
      this->Sz = chain.Sz;
      this->Momentum = chain.Momentum;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->ComplementaryStateShift = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
    }
}

// destructor
//

Spin1ChainWithTranslations::~Spin1ChainWithTranslations () 
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

Spin1ChainWithTranslations& Spin1ChainWithTranslations::operator = (const Spin1ChainWithTranslations& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->ChainDescription;
      delete[] this->LookUpTable;
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->ChainDescription = chain.ChainDescription;
      this->Sz = chain.Sz;
      this->Momentum = chain.Momentum;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->ComplementaryStateShift = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescription = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1ChainWithTranslations::Clone()
{
  return new Spin1ChainWithTranslations (*this);
}

// return Hilbert space dimension
//
// return value = Hilbert space dimension

int Spin1ChainWithTranslations::GetHilbertSpaceDimension()
{
  return this->HilbertSpaceDimension;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Spin1ChainWithTranslations::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->FixedQuantumNumberFlag == true)
    {
      TmpList += new SzQuantumNumber (this->Sz);
    }
  else
    {
      int TmpSz = - 2 * this->ChainLength;
      for (int i = 0; i <= (2 * this->ChainLength); i++)
	{
	  TmpList += new SzQuantumNumber (TmpSz);
	  TmpSz += 2;
	}
    }
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* Spin1ChainWithTranslations::GetQuantumNumber (int index)
{ 
  return new SzQuantumNumber (this->TotalSz(index));
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Spin1ChainWithTranslations::TotalSz (int index)
{
  if (this->FixedQuantumNumberFlag == true)
    return this->Sz;
  unsigned long State = this->ChainDescription[index];
  int TmpSz = 0;
  unsigned long TmpState;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpState = State & 0x00000003;
      switch (TmpState)
	{
	case 0x00000003:
	  TmpSz += 2;
	  break;
	case 0x00000000:
	  TmpSz -= 2;
	  break;
	}
      State >>= 2;
    }
  return TmpSz;
}

// return value of twice spin projection on (Oz) for a given state
//
// stateDescription = state to which the spin projection has to be evaluated
// return value = twice spin projection on (Oz)

int Spin1ChainWithTranslations::GetTotalSz (unsigned long stateDescription)
{
  int TmpSz = 0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      switch (stateDescription & 0x3)
	{
	case 0x00000003:
	  TmpSz += 2;
	  break;
	case 0x00000000:
	  TmpSz -= 2;
	  break;
	}
      stateDescription >>= 2;
    }
  return TmpSz;
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

double Spin1ChainWithTranslations::SziSzj (int i, int j, int state)
{  
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long tmpState2 = (tmpState >> (j << 1)) & 0x00000003;
  tmpState >>= (i << 1);
  tmpState &= 0x00000003;
  if ((tmpState == 0x00000002) || (tmpState2 == 0x00000002))
    return 0.0;
  if (tmpState == tmpState2)
    return 1.0;
  else
    return -1.0;
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslations::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{  
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x00000003;
  switch (tmpState)
    {
    case 0x00000003:
      coefficient = M_SQRT2;
      State &= ~(0x00000001 << i);
      break;
    case 0x00000002:
      coefficient = M_SQRT2;
      State&= ~(0x00000002 << i);
      break;
    case 0x00000000:
      return this->HilbertSpaceDimension;
    }	  
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x00000003;
  switch (tmpState2)
    {
    case 0x00000003:
      return this->HilbertSpaceDimension;
    case 0x00000002:
      coefficient *= M_SQRT2;
      return this->FindStateIndex(this->FindCanonicalForm(State | (0x00000001 << j), nbrTranslation));
    case 0x00000000:
      coefficient *= M_SQRT2;
      return this->FindStateIndex(this->FindCanonicalForm(State | (0x00000002 << j), nbrTranslation));
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslations::SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x00000003;
  switch (tmpState)
    {
    case 0x00000003:
      return this->HilbertSpaceDimension;
      break;
    case 0x00000002:
      coefficient = M_SQRT2;
      State |= (0x00000001 << i);
      break;
    case 0x00000000:
      coefficient = M_SQRT2;
      State |= (0x00000002 << i);
      break;
    }	  
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x00000003;
  switch (tmpState2)
    {
    case 0x00000003:
      return this->HilbertSpaceDimension;
    case 0x00000002:
      coefficient *= M_SQRT2;
      return this->FindStateIndex(this->FindCanonicalForm(State | (0x00000001 << j), nbrTranslation));
    case 0x00000000:
      coefficient *= M_SQRT2;
      return this->FindStateIndex(this->FindCanonicalForm(State | (0x00000002 << j), nbrTranslation));
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslations::SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x00000003;
  switch (tmpState)
    {
    case 0x00000003:
      coefficient = M_SQRT2;
      State &= ~(0x00000001 << i);
      break;
    case 0x00000002:
      coefficient = M_SQRT2;
      State&= ~(0x00000002 << i);
      break;
    case 0x00000000:
      return this->HilbertSpaceDimension;
    }	  
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x00000003;
  switch (tmpState2)
    {
    case 0x00000003:
      {
	coefficient *= M_SQRT2;
	return this->FindStateIndex(this->FindCanonicalForm(State & ~(0x00000001 << j), nbrTranslation));
      }
      break;
    case 0x00000002:
      {
	coefficient *= M_SQRT2;
	return this->FindStateIndex(this->FindCanonicalForm(State & ~(0x00000002 << j), nbrTranslation));
      }
      break;
    case 0x00000000:
      {
	coefficient = 0;
	return this->HilbertSpaceDimension;
      }
      break;
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i Sz_j operator on a given state
//
// i = position of S+ operator
// j = position of Sz operator
// state = index of the state to be applied on S+_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslations::SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = (tmpState >> (j << 1)) & 0x00000003;
  if (tmpState2 == 0x00000002)
    {
      coefficient = 0;
      return this->HilbertSpaceDimension;
    }
  if (tmpState2 == 0x00000003)
    {
      coefficient = 1.0;
    }
  else
    {
      coefficient = -1.0;
    }
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x00000003;
  switch (tmpState)
    {
    case 0x00000003:
      {
	coefficient = 0;
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x00000002:
      {
	coefficient *= M_SQRT2;
	return this->FindStateIndex(this->FindCanonicalForm(State | (0x00000001 << i), nbrTranslation));
      }
      break;
    case 0x00000000:
      {
	coefficient *= M_SQRT2;
	return this->FindStateIndex(this->FindCanonicalForm(State | (0x00000002 << i), nbrTranslation));
      }
      break;
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1ChainWithTranslations::SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long tmpState = this->ChainDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = (tmpState >> (j << 1)) & 0x00000003;
  if (tmpState2 == 0x00000002)
    {
      coefficient = 0;
      return this->HilbertSpaceDimension;
    }
  if (tmpState2 == 0x00000003)
    {
      coefficient = 1.0;
    }
  else
    {
      coefficient = -1.0;
    }
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x00000003;
  switch (tmpState)
    {
    case 0x00000003:
      {
	coefficient *= M_SQRT2;
	return this->FindStateIndex(this->FindCanonicalForm(State & ~(0x00000001 << i), nbrTranslation));
      }
      break;
    case 0x00000002:
      {
	coefficient *= M_SQRT2;
	return this->FindStateIndex(this->FindCanonicalForm(State & ~(0x00000001 << i), nbrTranslation));
      }
      break;
    case 0x00000000:
      {
	coefficient *= 0;
	return this->HilbertSpaceDimension;
      }
      break;
    }	  
  return this->HilbertSpaceDimension;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* Spin1ChainWithTranslations::ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter)
{
/*  if (q.GetQuantumNumberType() != AbstractQuantumNumber::Sz)
    return new Spin1ChainWithTranslations();
  int TmpSz = ((SzQuantumNumber&) q).GetSz();
  int HilbertSubspaceDimension = 0;
  int* TmpConvArray = new int [this->HilbertSpaceDimension];
  for (int i = 0; i < this->HilbertSpaceDimension; i++)
    {
      if (this->TotalSz(i) == TmpSz)
	{
	  TmpConvArray[HilbertSubspaceDimension] = i;
	  HilbertSubspaceDimension++;	  
	}
    }
  int* ConvArray = new int [HilbertSubspaceDimension];
  unsigned long* SubspaceDescription = new unsigned long [HilbertSubspaceDimension];
  int* SubspaceLookUpTable = new int [this->LookUpTableSize];
  unsigned long TestMask = this->ChainDescription[TmpConvArray[0]] & this->LookUpTableMask;
  SubspaceLookUpTable[TestMask] = 0;
  SubspaceDescription[0] = this->ChainDescription[TmpConvArray[0]];
  ConvArray[0] = TmpConvArray[0];
//  unsigned long TestMask = this->LookUpTableMask;
  for (int i = 1; i < HilbertSubspaceDimension; i++)
    {
      if ((this->ChainDescription[TmpConvArray[i]] & this->LookUpTableMask) != TestMask)
	{
	  TestMask = this->ChainDescription[TmpConvArray[i]] & this->LookUpTableMask;
	  SubspaceLookUpTable[TestMask] = i;
	}
      SubspaceDescription[i] = this->ChainDescription[TmpConvArray[i]];
      ConvArray[i] = TmpConvArray[i];
    }
  converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, HilbertSubspaceDimension, ConvArray);
  return new Spin1ChainWithTranslations (HilbertSubspaceDimension, SubspaceDescription, this->ChainLength,
					 TmpSz, true, SubspaceLookUpTable, this->LookUpTableSize, 
					 this->LookUpPosition, this->LookUpTableMask);*/
  return 0;
}

// find the canonical form of a state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// return value = canonical form of the state

unsigned long Spin1ChainWithTranslations::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState = stateDescription;
  int index = 1;  
  while (index < this->ChainLength)
    {
      stateDescription = (stateDescription >> 2) | ((stateDescription & 0x3) << this->ComplementaryStateShift);
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = index;
	}
      ++index;
    }
  return CanonicalState;
}

// find how many translations are needed to obtain the same state
//
// stateDescription = unsigned integer describing the state
// return value = number of translation needed to obtain the same state

int Spin1ChainWithTranslations::FindNumberTranslation(unsigned long stateDescription)
{
  unsigned long TmpState = (stateDescription >> 2) | ((stateDescription & 0x3) << this->ComplementaryStateShift);
  int index = 1;  
  while (TmpState != stateDescription)
    {
      TmpState = (TmpState >> 2) | ((TmpState & 0x3) << this->ComplementaryStateShift);
      ++index;
    }
  return index;
}

// find state index
//
// state = state description
// return value = corresponding index

int Spin1ChainWithTranslations::FindStateIndex(unsigned long state)
{
  unsigned long MidPos = state >> this->LookUpTableShift;
//  unsigned long LowPos = this->LookUpTable[MidPos];
//  unsigned long HighPos = this->LookUpTable[MidPos + 1];
  unsigned long LowPos = 0;
  unsigned long HighPos = this->HilbertSpaceDimension - 1;
  while (HighPos != LowPos)
    {
      MidPos = (HighPos + LowPos) >> 1;
      if (this->ChainDescription[MidPos] >= state)
	HighPos = MidPos;
      else
	LowPos = MidPos;
    }
  return HighPos;   
}



// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Spin1ChainWithTranslations::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long tmp;
  unsigned long StateDescription = this->ChainDescription[state];  
//  Str << this->FindStateIndex(StateDescription) << " : ";
  Str << state << " : ";
  for (int j = 0; j < this->ChainLength; j++)
    {
      tmp = StateDescription & 0x3;
      if (tmp == 0)
	Str << "-1 ";
      else
	if (tmp == 0x2)
	  Str << "0 ";
	else
	  Str << "1 ";
      StateDescription >>= 2;
    }
  return Str;
}


// generate all states with no constraint on total Sz
//
// memorySlice = maximum amount of memory (in unsigned long unit) that can be allocated to partially evalauted the states
// return value = number of generated states

int Spin1ChainWithTranslations::GenerateStates(long memorySlice) 
{
  long MaximumNbrState = 3;
  for (int i = 1; i < this->ChainLength; ++i)
    {
      MaximumNbrState *= 3;
    }

  bool* CompatibilityWithMomentum = new bool [this->ChainLength + 1];
  for (int i = 1; i <= this->ChainLength; ++i)
    if (((i * this->Momentum) % this->ChainLength) == 0)
      CompatibilityWithMomentum[i] = true;
    else
      CompatibilityWithMomentum[i] = false;

  List<unsigned long*> TmpGeneratedStateList;
  List<long> TmpNbrGeneratedStateList;
  unsigned long TmpState = 0;
  unsigned long TmpState2;
  int NbrTranslation = 0;
  int TmpHilbertSpaceDimension = 0;
  unsigned long* TmpGeneratedStates = new unsigned long [memorySlice];
  int Shift = 0;
  int TwiceChainLength = (this->ChainLength << 1) - 2;
  while (MaximumNbrState > 0)
    {
      //test each state
      long Pos = 0;
      while ((Pos < memorySlice) && (MaximumNbrState > 0))
	{
	  Shift = TwiceChainLength;	  
	  while ((Shift >= 0) && (((TmpState >> Shift) & 0x3) != 0x1))
	    {
	      Shift -= 2;
	    }
	  while (Shift >= 0)
	    {
	      ++TmpState;
	      Shift = TwiceChainLength;	  
	      while ((Shift >= 0) && (((TmpState >> Shift) & 0x3) != 0x1))
		{
		  Shift -= 2;
		}
	    }
	  TmpState2 = this->FindCanonicalForm(TmpState, NbrTranslation);
	  if ((NbrTranslation == 0))// && (CompatibilityWithMomentum[this->FindNumberTranslation(TmpState2)]))
	    {
	      TmpGeneratedStates[Pos] = TmpState2;
	      ++Pos;
	    }
	  ++TmpState;
	  --MaximumNbrState;
	}
      cout << "Pos = " << Pos << " MaximumNbrState " << MaximumNbrState << endl;
      // delete duplicate entries
      if (Pos > 0)
	{
/*	  SortArrayUpOrdering(TmpGeneratedStates, Pos);
	  unsigned long* TmpGeneratedStates2 = new unsigned long [Pos];
	  --Pos;
	  int Pos2 = 0;
	  for (int i = 0; i < Pos; ++i)
	    {	      
	      TmpGeneratedStates2[Pos2] = TmpGeneratedStates[i];
	      ++Pos2;
	      while ((i < Pos) && (TmpGeneratedStates[i] == TmpGeneratedStates[i + 1]))
		++i;
	    }
	  if (TmpGeneratedStates[Pos] != TmpGeneratedStates[Pos - 1])
	    {
	      TmpGeneratedStates2[Pos2] = TmpGeneratedStates[Pos];
	      ++Pos2;	      
	    }
	  TmpGeneratedStateList += TmpGeneratedStates2;
	  TmpNbrGeneratedStateList += Pos2;
	  TmpHilbertSpaceDimension += Pos2;*/
	  TmpGeneratedStateList += TmpGeneratedStates;
	  TmpNbrGeneratedStateList += Pos;
	  TmpHilbertSpaceDimension += Pos;
	}
    }
  delete[] CompatibilityWithMomentum;
  this->ChainDescription = SmartMergeArrayListIntoArray(TmpGeneratedStateList, TmpNbrGeneratedStateList);
//  SortArrayUpOrdering(this->ChainDescription, TmpHilbertSpaceDimension);

  // create the look-up table
  unsigned long Max = ((unsigned long) 1) << ((this->ChainLength << 1) - this->LookUpTableShift - 1);
  cout << "Max " << Max << endl;
  this->LookUpTable = new long [Max + 1];
  long LowPos;
  long MidPos;
  long HighPos;
  for (unsigned long i = 0; i < Max; ++i)
    {
      LowPos = 0;
      HighPos = TmpHilbertSpaceDimension - 1;
      while ((HighPos - LowPos) > 1)
	{
	  MidPos = (HighPos + LowPos) >> 1;
	  if (this->ChainDescription[MidPos] >= i)
	    HighPos = MidPos;
	  else
	    LowPos = MidPos;
	}      
      this->LookUpTable[i] = LowPos;
    }
  this->LookUpTable[Max] = this->ChainDescription[TmpHilbertSpaceDimension - 1];
  return TmpHilbertSpaceDimension;
}

// generate all states with constraint on total Sz
//
// sz = twice the sz value
// memorySlice = maximum amount of memory (in unsigned long unit) that can be allocated to partially evalauted the states
// return value = number of generated states

int Spin1ChainWithTranslations::GenerateStates(int sz, long memorySlice) 
{
  long MaximumNbrState = 3;
  for (int i = 1; i < this->ChainLength; ++i)
    {
      MaximumNbrState *= 3;
    }

  bool* CompatibilityWithMomentum = new bool [this->ChainLength];
  for (int i = 0; i < this->ChainLength; ++i)
    if (((i * this->Momentum) % this->ChainLength))
      CompatibilityWithMomentum[i] = true;
    else
      CompatibilityWithMomentum[i] = false;

  List<unsigned long*> TmpGeneratedStateList;
  List<long> TmpNbrGeneratedStateList;
  unsigned long TmpState = 0;
  unsigned long TmpState2;
  int NbrTranslation;
  int TmpHilbertSpaceDimension = 0;
  unsigned long* TmpGeneratedStates = new unsigned long [memorySlice];
  int Shift = 0;
  int TwiceChainLength = (this->ChainLength << 1) - 2;
  while (MaximumNbrState > 0)
    {
      //test each state
      long Pos = 0;
      while ((Pos < memorySlice) && (MaximumNbrState > 0))
	{
	  Shift = TwiceChainLength;	  
	  while ((Shift >= 0) && ((TmpState >> Shift) != 0x1))
	    {
	      Shift -= 2;
	    }
	  while (Shift < 0)
	    {
	      ++TmpState;
	      Shift = TwiceChainLength;	  
	      while ((Shift >= 0) && ((TmpState >> Shift) != 0x1))
		{
		  Shift -= 2;
		}
	    }
	  TmpState2 = this->FindCanonicalForm(TmpState, NbrTranslation);
	  if ((NbrTranslation == 0) && (CompatibilityWithMomentum[this->FindNumberTranslation(TmpState2)]))
	    {
	      TmpGeneratedStates[Pos] = TmpState2;
	      ++Pos;
	    }
	  ++TmpState;
	  --MaximumNbrState;
	}
      // delete duplicate entries
      if (Pos > 0)
	{
/*	  SortArrayUpOrdering(TmpGeneratedStates, Pos);
	  unsigned long* TmpGeneratedStates2 = new unsigned long [Pos];
	  --Pos;
	  int Pos2 = 0;
	  for (int i = 0; i < Pos; ++i)
	    {	      
	      TmpGeneratedStates2[Pos2] = TmpGeneratedStates[i];
	      ++Pos2;
	      while ((i < Pos) && (TmpGeneratedStates[i] == TmpGeneratedStates[i + 1]))
		++i;
	    }
	  if (TmpGeneratedStates[Pos] != TmpGeneratedStates[Pos - 1])
	    {
	      TmpGeneratedStates2[Pos2] = TmpGeneratedStates[Pos];
	      ++Pos2;	      
	    }
	  TmpGeneratedStateList += TmpGeneratedStates2;
	  TmpNbrGeneratedStateList += Pos2;
	  TmpHilbertSpaceDimension += Pos2;*/
	  TmpGeneratedStateList += TmpGeneratedStates;
	  TmpNbrGeneratedStateList += Pos;
	  TmpHilbertSpaceDimension += Pos;
	}
    }
  delete[] TmpGeneratedStates;
  delete[] CompatibilityWithMomentum;
  this->ChainDescription = SmartMergeArrayListIntoArray(TmpGeneratedStateList, TmpNbrGeneratedStateList);
//  SortArrayUpOrdering(this->ChainDescription, TmpHilbertSpaceDimension);

  // create the look-up table
  unsigned long Max = ((unsigned long) 1) << ((this->ChainLength << 1) - this->LookUpTableShift - 1);
  this->LookUpTable = new long [Max + 1];
  long LowPos;
  long MidPos;
  long HighPos;
  for (unsigned long i = 0; i < Max; ++i)
    {
      LowPos = 0;
      HighPos = TmpHilbertSpaceDimension - 1;
      while ((HighPos - LowPos) > 1)
	{
	  MidPos = (HighPos + LowPos) >> 1;
	  if (this->ChainDescription[MidPos] >= i)
	    HighPos = MidPos;
	  else
	    LowPos = MidPos;
	}      
      this->LookUpTable[i] = LowPos;
    }
  this->LookUpTable[Max] = this->ChainDescription[TmpHilbertSpaceDimension - 1];
  return TmpHilbertSpaceDimension;
}

