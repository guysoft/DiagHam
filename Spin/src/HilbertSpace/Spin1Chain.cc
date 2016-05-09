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


#include "HilbertSpace/Spin1Chain.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;


#define MAX_HILBERTSPACE_DIMENSION 327680
#ifndef M_SQRT2
#define M_SQRT2	1.41421356237309504880
#endif


// default constructor
//

Spin1Chain::Spin1Chain () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableSize = 0;
  this->LookUpTableMask = 0;
  this->LookUpPosition = 0;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->ChainLength = 0;
  this->Sz = 0;
  this->FixedQuantumNumberFlag = false;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// constructor for complete Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// memorySize = memory size in bytes allowed for look-up table

Spin1Chain::Spin1Chain (int chainLength, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedQuantumNumberFlag = false;
  this->HilbertSpaceDimension = 3;
  for (int i = 1; i < chainLength; i++)
    this->HilbertSpaceDimension *= 3;

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

  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->GenerateStates (0, 0, 0xffffffff);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// constructor for complete Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table

Spin1Chain::Spin1Chain (int chainLength, int sz, int memorySize) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Sz = sz;
  this->FixedQuantumNumberFlag = true;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->Sz, this->ChainLength);

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

  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->HilbertSpaceDimension = this->GenerateStates (0, 0, 0xffffffff, this->ChainLength * 2);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// constructor from pre-constructed datas
//
// hilbertSpaceDimension = Hilbert space dimension
// chainDescription = array describing states
// chainLength = number of spin 1
// sz = twice the value of total Sz component
// fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
// lookUpTable = look-up table
// lookUpTableSize = look-Up table size
// lookUpTablePosition = last position described by the look-Up table
// lookUpTableMask = look-Up table mask

Spin1Chain::Spin1Chain (int hilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
			int sz, bool fixedQuantumNumberFlag, int* lookUpTable, int lookUpTableSize, 
			int lookUpPosition, unsigned long lookUpTableMask)
{
  this->Flag.Initialize();
  this->LookUpTable = lookUpTable;
  this->LookUpTableMask = lookUpTableMask;
  this->LookUpPosition = lookUpPosition;
  this->LookUpTableSize = lookUpTableSize;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->StateDescription = chainDescription;
  this->Sz = sz;
  this->FixedQuantumNumberFlag = fixedQuantumNumberFlag;
  this->ChainLength = chainLength;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}
  
// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1Chain::Spin1Chain (const Spin1Chain& chain) 
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
      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableSize = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->HilbertSpaceDimension = 0;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

Spin1Chain::~Spin1Chain () 
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1Chain& Spin1Chain::operator = (const Spin1Chain& chain)
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTable;
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->StateDescription = chain.StateDescription;
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableMask = chain.LookUpTableMask;
      this->LookUpTableSize = chain.LookUpTableSize;
      this->LookUpPosition = chain.LookUpPosition;
      this->Sz = chain.Sz;
      this->FixedQuantumNumberFlag = chain.FixedQuantumNumberFlag;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableSize = 0;
      this->LookUpTableMask = 0;
      this->LookUpPosition = 0;
      this->HilbertSpaceDimension = 0;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Sz = 0;
      this->FixedQuantumNumberFlag = false;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1Chain::Clone()
{
  return new Spin1Chain (*this);
}

// evaluate Hilbert space dimension
//
// sz = twice the Sz value
// nbrSites = number of sites
// return value = Hilbert space dimension

long Spin1Chain::EvaluateHilbertSpaceDimension(int sz, int nbrSites)
{
  if (nbrSites == 0)
    {
      if (sz == 0)
	{
	  return 1l;	  
	}
      else
	{
	  return 0l;	  
	}
    }
  long TmpDimension = this->EvaluateHilbertSpaceDimension(sz - 2, nbrSites - 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(sz, nbrSites - 1);
  TmpDimension += this->EvaluateHilbertSpaceDimension(sz + 2, nbrSites - 1);
  return TmpDimension;
}

// generate all states with no constraint on total Sz
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// currentStateDescription = description of current state
// return value = number of generated states

int Spin1Chain::GenerateStates(int statePosition,int sitePosition, unsigned long currentStateDescription) 
{
  int NextSitePosition = sitePosition + 1;     
  int NbrGeneratedState = -statePosition;
  sitePosition &= 0x0000000f;
  sitePosition <<= 1;
  unsigned long mask;
  if (NextSitePosition != this->ChainLength)
    {
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription);
      
      mask = ~(0x1ul << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription & mask);
      
      mask = ~(0x3ul << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription & mask);
    }
  else
    {
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      this->StateDescription[statePosition] = currentStateDescription;
      statePosition++;
      
      mask = ~(0x1ul << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      this->StateDescription[statePosition] = (currentStateDescription & mask);
      statePosition++;
      
      mask = ~(0x3ul << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      this->StateDescription[statePosition] = (currentStateDescription & mask);
      statePosition++;
    }
  NbrGeneratedState += statePosition;
  return NbrGeneratedState;
}

// generate all states corresponding to a given total Sz
//
// statePosition = position for the new states
// sitePosition = site on chain where spin has to be changed
// currentStateDescription = description of current state
// currentSz = total Sz value of current state
// return value = number of generated states

int Spin1Chain::GenerateStates(int statePosition,int sitePosition, unsigned long currentStateDescription, int currentSz) 
{
  int MaxSz = (this->ChainLength - sitePosition) * 4;
  if ((currentSz - this->Sz) > MaxSz)
    return 0;
  int NextSitePosition = sitePosition + 1;     
  int NbrGeneratedState = -statePosition;
  sitePosition &= 0x0000000f;
  sitePosition <<= 1;
  unsigned long mask;
  if (NextSitePosition != this->ChainLength)
    {
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription, currentSz);
      
      mask = ~(0x1ul << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription & mask, currentSz - 2);
      
      mask = ~(0x3ul << sitePosition);
      if (NextSitePosition == this->LookUpPosition)
	this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
      statePosition += this->GenerateStates(statePosition, NextSitePosition, currentStateDescription & mask, currentSz - 4);
    }
  else
    {
      if (currentSz == this->Sz)
	{
	  if (NextSitePosition == this->LookUpPosition)
	    this->LookUpTable[currentStateDescription & this->LookUpTableMask] = statePosition;
	  this->StateDescription[statePosition] = currentStateDescription;
	  statePosition++;
	}
      else
	{
	  currentSz -= 2;
	  if (currentSz == this->Sz)
	    {
	      mask = ~(0x1ul << sitePosition);
	      if (NextSitePosition == this->LookUpPosition)
		this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
	      this->StateDescription[statePosition] = (currentStateDescription & mask);
	      statePosition++;
	    }
	  else
	    {
	      currentSz -= 2;
	      if (currentSz == this->Sz)
		{
		  mask = ~(0x3ul << sitePosition);
		  if (NextSitePosition == this->LookUpPosition)
		    this->LookUpTable[(currentStateDescription & mask) & this->LookUpTableMask] = statePosition;
		  this->StateDescription[statePosition] = (currentStateDescription & mask);
		  statePosition++;
		}
	    }
	}
    }
  NbrGeneratedState += statePosition;
  return NbrGeneratedState;
}

// return Hilbert space dimension
//
// return value = Hilbert space dimension

int Spin1Chain::GetHilbertSpaceDimension()
{
  return this->HilbertSpaceDimension;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Spin1Chain::GetQuantumNumbers ()
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

AbstractQuantumNumber* Spin1Chain::GetQuantumNumber (int index)
{ 
  return new SzQuantumNumber (this->TotalSz(index));
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Spin1Chain::TotalSz (int index)
{
  if (this->FixedQuantumNumberFlag == true)
    return this->Sz;
  unsigned long State = this->StateDescription[index];
  int TmpSz = 0;
  unsigned long TmpState;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpState = State & 0x3ul;
      switch (TmpState)
	{
	case 0x3ul:
	  TmpSz += 2;
	  break;
	case 0x0ul:
	  TmpSz -= 2;
	  break;
	}
      State >>= 2;
    }
  return TmpSz;
}

// return matrix representation of Sx
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& Spin1Chain::Sxi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  i <<= 1;
  unsigned long tmpState;
  unsigned long State;
  double Factor = M_SQRT2 * 0.5;
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      tmpState = this->StateDescription[j];
      State = tmpState;
      tmpState >>= i;
      tmpState &= 0x3ul;
      switch (tmpState)
	{
	case 0x3ul:
	  M(this->FindStateIndex((State & ~(0x1ul << i))), j) = Factor;
	  break;
	case 0x2ul:
	  M(this->FindStateIndex((State | (0x1ul << i))), j) = Factor;
	  M(this->FindStateIndex((State & ~(0x2ul << i))), j) = Factor;
	  break;
	case 0x0ul:
	  M(this->FindStateIndex((State | (0x2ul << i))), j) = Factor;
	  break;
	}
    }
  return M;
}

// return matrix representation of i * Sy
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& Spin1Chain::Syi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  i <<= 1;
  unsigned long tmpState;
  unsigned long State;
  double Factor = M_SQRT2 * 0.5;
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      tmpState = this->StateDescription[j];
      State = tmpState;
      tmpState >>= i;
      tmpState &= 0x3ul;
      switch (tmpState)
	{
	case 0x3ul:
	  M(this->FindStateIndex((State & ~(0x1ul << i))), j) = -Factor;
	  break;
	case 0x2ul:
	  M(this->FindStateIndex((State | (0x1ul << i))), j) = Factor;
	  M(this->FindStateIndex((State & ~(0x2ul << i))), j) = -Factor;
	  break;
	case 0x0ul:
	  M(this->FindStateIndex((State | (0x2ul << i))), j) = Factor;
	  break;
	}
    }
  return M;
}

// return matrix representation of Sz
//
// i = operator position
// M = matrix where representation has to be stored
// return value = corresponding matrix

Matrix& Spin1Chain::Szi (int i, Matrix& M)
{
  if ((M.GetNbrRow() != this->HilbertSpaceDimension) || (M.GetNbrColumn() != this->HilbertSpaceDimension))
    M.ResizeAndClean(this->HilbertSpaceDimension, this->HilbertSpaceDimension);
  i <<= 1;
  unsigned long tmpState;
  unsigned long State;
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      tmpState = this->StateDescription[j];
      State = tmpState;
      tmpState >>= i;
      tmpState &= 0x3ul;
      switch (tmpState)
	{
	case 0x3ul:
	  M(j, j) = 1.0;
	  break;
	case 0x0ul:
	  M(j, j) = -1.0;
	  break;
	}
    }
  return M;
}

// return index of resulting state from application of S+_i operator on a given state
//
// i = position of S+ operator
// state = index of the state to be applied on S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1Chain::Spi (int i, int state, double& coefficient)
{
  unsigned long State = (this->StateDescription[state] >> (i << 1));
  unsigned long tmpState = State;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x2ul:
      {
	coefficient = M_SQRT2;
	return this->FindStateIndex(State | (0x1ul << (i << 1)));
      }
      break;
    case 0x0ul:
      {
	coefficient = M_SQRT2;
	return this->FindStateIndex(State | (0x2ul << (i << 1)));
      }
      break;
    case 0x3ul:
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
      break;
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i operator on a given state
//
// i = position of S- operator
// state = index of the state to be applied on S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1Chain::Smi (int i, int state, double& coefficient)
{
  unsigned long State = (this->StateDescription[state] >> (i << 1));
  unsigned long tmpState = State;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x2ul:
      {
	coefficient = M_SQRT2;
	return this->FindStateIndex(State & ~(0x2ul << (i << 1)));
      }
      break;
    case 0x0ul:
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x3ul:
      {
	coefficient = M_SQRT2;
	return this->FindStateIndex(State & ~(0x1ul << (i << 1)));
      }
      break;
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of Sz_i operator on a given state
//
// i = position of Sz operator
// state = index of the state to be applied on Sz_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1Chain::Szi (int i, int state, double& coefficient)
{
  unsigned long tmpState = (this->StateDescription[state] >> (i << 1));
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x2ul:
      coefficient = 0.0;
      break;
    case 0x0ul:
      coefficient = -1.0;
      break;
    case 0x3ul:
      coefficient = 1.0;
      break;
    }
  return state;
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

double Spin1Chain::SziSzj (int i, int j, int state)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpState2 = (tmpState >> (j << 1)) & 0x3ul;
  tmpState >>= (i << 1);
  tmpState &= 0x3ul;
  if ((tmpState == 0x2ul) || (tmpState2 == 0x2ul))
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
// return value = index of resulting state

int Spin1Chain::SmiSpj (int i, int j, int state, double& coefficient)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      coefficient = M_SQRT2;
      State &= ~(0x1ul << i);
      break;
    case 0x2ul:
      coefficient = M_SQRT2;
      State&= ~(0x2ul << i);
      break;
    case 0x0ul:
      return this->HilbertSpaceDimension;
    }	  
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  switch (tmpState2)
    {
    case 0x3ul:
      return this->HilbertSpaceDimension;
    case 0x2ul:
      coefficient *= M_SQRT2;
      return this->FindStateIndex(State | (0x1ul << j));
    case 0x0ul:
      coefficient *= M_SQRT2;
      return this->FindStateIndex(State | (0x2ul << j));
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i S+_j operator on a given state
//
// i = position of first S+ operator
// j = position of second S+ operator
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1Chain::SpiSpj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      return this->HilbertSpaceDimension;
      break;
    case 0x2ul:
      coefficient = M_SQRT2;
      State |= (0x1ul << i);
      break;
    case 0x0ul:
      coefficient = M_SQRT2;
      State |= (0x2ul << i);
      break;
    }	  
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  switch (tmpState2)
    {
    case 0x3ul:
      return this->HilbertSpaceDimension;
    case 0x2ul:
      coefficient *= M_SQRT2;
      return this->FindStateIndex(State | (0x1ul << j));
    case 0x0ul:
      coefficient *= M_SQRT2;
      return this->FindStateIndex(State | (0x2ul << j));
    }	  
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_j operator on a given state
//
// i = position of first S- operator
// j = position of second S- operator
// state = index of the state to be applied on S-_i S-_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// return value = index of resulting state

int Spin1Chain::SmiSmj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      coefficient = M_SQRT2;
      State &= ~(0x1ul << i);
      break;
    case 0x2ul:
      coefficient = M_SQRT2;
      State&= ~(0x2ul << i);
      break;
    case 0x0ul:
      return this->HilbertSpaceDimension;
    }	  
  j <<= 1;
  tmpState2 >>= j;
  tmpState2 &= 0x3ul;
  switch (tmpState2)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT2;
	return this->FindStateIndex(State & ~(0x1ul << j));
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT2;
	return this->FindStateIndex(State & ~(0x2ul << j));
      }
      break;
    case 0x0ul:
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
// return value = index of resulting state

int Spin1Chain::SpiSzj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = (tmpState >> (j << 1)) & 0x3ul;
  if (tmpState2 == 0x2ul)
    {
      coefficient = 0;
      return this->HilbertSpaceDimension;
    }
  if (tmpState2 == 0x3ul)
    {
      coefficient = 1.0;
    }
  else
    {
      coefficient = -1.0;
    }
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      {
	coefficient = 0;
	return this->HilbertSpaceDimension;
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT2;
	return this->FindStateIndex(State | (0x1ul << i));
      }
      break;
    case 0x0ul:
      {
	coefficient *= M_SQRT2;
	return this->FindStateIndex(State | (0x2ul << i));
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
// return value = index of resulting state

int Spin1Chain::SmiSzj (int i, int j, int state, double& coefficient)
{
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = (tmpState >> (j << 1)) & 0x3ul;
  if (tmpState2 == 0x2ul)
    {
      coefficient = 0;
      return this->HilbertSpaceDimension;
    }
  if (tmpState2 == 0x3ul)
    {
      coefficient = 1.0;
    }
  else
    {
      coefficient = -1.0;
    }
  i <<= 1;
  tmpState >>= i;
  tmpState &= 0x3ul;
  switch (tmpState)
    {
    case 0x3ul:
      {
	coefficient *= M_SQRT2;
	return this->FindStateIndex(State & ~(0x1ul << i));
      }
      break;
    case 0x2ul:
      {
	coefficient *= M_SQRT2;
	return this->FindStateIndex(State & ~(0x1ul << i));
      }
      break;
    case 0x0ul:
      {
	coefficient *= 0;
	return this->HilbertSpaceDimension;
      }
      break;
    }	  
  return this->HilbertSpaceDimension;
}

// translate a state assuming the system have periodic boundary conditions (increasing the site index)
//
// nbrTranslations = number of translations to apply
// state = index of the state to translate 
// return value = index of resulting state

int Spin1Chain::TranslateState (int nbrTranslations, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  TmpState = (((TmpState & (0x3ul << ((this->ChainLength - nbrTranslations) - 1ul)) << 1) << nbrTranslations)
	      | (TmpState >> ((this->ChainLength - nbrTranslations) << 1)));
  return this->FindStateIndex(TmpState);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* Spin1Chain::ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter)
{
  if (q.GetQuantumNumberType() != AbstractQuantumNumber::Sz)
    return new Spin1Chain();
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
  unsigned long TestMask = this->StateDescription[TmpConvArray[0]] & this->LookUpTableMask;
  SubspaceLookUpTable[TestMask] = 0;
  SubspaceDescription[0] = this->StateDescription[TmpConvArray[0]];
  ConvArray[0] = TmpConvArray[0];
//  unsigned long TestMask = this->LookUpTableMask;
  for (int i = 1; i < HilbertSubspaceDimension; i++)
    {
      if ((this->StateDescription[TmpConvArray[i]] & this->LookUpTableMask) != TestMask)
	{
	  TestMask = this->StateDescription[TmpConvArray[i]] & this->LookUpTableMask;
	  SubspaceLookUpTable[TestMask] = i;
	}
      SubspaceDescription[i] = this->StateDescription[TmpConvArray[i]];
      ConvArray[i] = TmpConvArray[i];
    }
  converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, HilbertSubspaceDimension, ConvArray);
  return (AbstractSpinChain*) new Spin1Chain (HilbertSubspaceDimension, SubspaceDescription, this->ChainLength,
					      TmpSz, true, SubspaceLookUpTable, this->LookUpTableSize, 
					      this->LookUpPosition, this->LookUpTableMask);
}

// find state index
//
// state = state description
// return value = corresponding index

int Spin1Chain::FindStateIndex(unsigned long state)
{
  int index = this->LookUpTable[state & this->LookUpTableMask];
  unsigned long* tmpState = &(this->StateDescription[index]);
  while ((index < this->HilbertSpaceDimension) && (state != *(tmpState++)))
    index++;
  return index;   
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Spin1Chain::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long tmp;
  unsigned long StateDescription = this->StateDescription[state];  
  //Str << this->FindStateIndex(StateDescription) << " : ";
  for (int j = 0; j < this->ChainLength; j++)
    {
      tmp = StateDescription & 0x3ul;
      if (tmp == 0x0ul)
	Str << "-1 ";
      else
	if (tmp == 0x2ul)
	  Str << "0 ";
	else
	  Str << "1 ";
      StateDescription >>= 2;
    }
  return Str;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix Spin1Chain::EvaluatePartialDensityMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      if (szSector == 0)
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
  if (nbrSites == this->ChainLength)
    {
      if (szSector == this->Sz)
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
  Spin1Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin1Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = nbrSites * 2;
  unsigned long Mask = (0x1ul << Shift) - 0x1ul;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;
  long TmpNbrNonZeroElements = 0l;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = (TmpHilbertSpace.StateDescription[MinIndex] << Shift) & 0xfffffffful;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask);
	  int TmpPos = this->FindStateIndex(TmpState2);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]];
	      for (int k = 0; k < Pos; ++k)
		{
		  if (TmpStatePosition2[k] >= Pos2)
		    {
		      TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]]);
		    }
		}
	    }
	}
    }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  return TmpDensityMatrix;
}
	
// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix Spin1Chain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      if (szSector == 0)
	{
	  RealMatrix TmpEntanglementMatrix(1, 1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
      
    }
  if (nbrSites == this->ChainLength)
    {
      if (szSector == this->Sz)
	{
	  RealMatrix TmpEntanglementMatrix(1, 1);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}      
    }
  Spin1Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin1Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);
  RealMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);

  int Shift = 2 * nbrSites;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << Shift;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | TmpDestinationHilbertSpace.StateDescription[j];
	  int TmpPos = this->FindStateIndex(TmpState2);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	       TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos]);
	    }
	}
    }

   return TmpEntanglementMatrix;
}
	
// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix Spin1Chain::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrSites == 0)
    {
      if (szSector == 0)
	{
	  ComplexMatrix TmpEntanglementMatrix(1, 1);
          Complex Tmp(1.0, 0.0);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
      
    }
  if (nbrSites == this->ChainLength)
    {
      if (szSector == this->Sz)
	{
	  ComplexMatrix TmpEntanglementMatrix(1, 1);
          Complex Tmp(1.0, 0.0);
	  TmpEntanglementMatrix.SetMatrixElement(0, 0, Tmp);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}      
    }
  Spin1Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin1Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);

  ComplexMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int Shift = 2 * nbrSites;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;

  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << Shift;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState2 = TmpState | TmpDestinationHilbertSpace.StateDescription[j];
	  int TmpPos = this->FindStateIndex(TmpState2);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	       TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos]);
	    }
	}
    }

   return TmpEntanglementMatrix;
}
