////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of doubled spin 0 +1/2 chain with translations           //
//                                                                            //
//                        last modification : 21/01/2016                      //
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


#include "HilbertSpace/DoubledSpin0_1_2_chainWithTranslations.h"


#include <iostream>
#include <math.h>

using std::cout;
using std::endl;


// default constructor
//

DoubledSpin0_1_2_chainWithTranslations::DoubledSpin0_1_2_chainWithTranslations () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableShift = 0;
  this->HilbertSpaceDimension = 0;
  this->ChainDescriptionBra = 0;
  this->ChainDescriptionKet = 0;
  this->ChainLength = 0;
  this->Momentum = 0;
  this->ComplementaryStateShift = 0;
  this->DiffSz = 0;
  this->FixedSpinProjectionFlag = false;
  this->CompatibilityWithMomentum = 0;
  this->RescalingFactors = 0;
  this->NbrStateInOrbit = 0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}



// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

DoubledSpin0_1_2_chainWithTranslations::DoubledSpin0_1_2_chainWithTranslations (int chainLength, int diffSz, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = diffSz;
  this->FixedSpinProjectionFlag = true;
  this->ComplementaryStateShift = 2*(this->ChainLength - 1);
  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;
  
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, this->DiffSz);
  
  this->ChainDescriptionBra = new unsigned long [this->LargeHilbertSpaceDimension];
  this->ChainDescriptionKet = new unsigned long [this->LargeHilbertSpaceDimension];
  
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, this->DiffSz, 0l);
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2Spin!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  

  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memorySize);
    }
  this->RescalingFactors = 0;
}
 


// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

DoubledSpin0_1_2_chainWithTranslations::DoubledSpin0_1_2_chainWithTranslations (int chainLength, int momentum, int diffSz, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->DiffSz = diffSz;
  this->FixedSpinProjectionFlag = true;
  this->Momentum = momentum;
  this->ComplementaryStateShift = 2*(this->ChainLength - 1);
  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;
  
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, this->DiffSz);

  this->ChainDescriptionBra = new unsigned long [this->LargeHilbertSpaceDimension];
  this->ChainDescriptionKet = new unsigned long [this->LargeHilbertSpaceDimension];
  
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, this->DiffSz, 0l);
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2Spin!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = 0l;
  unsigned long TmpCanonicalStateBra;
  unsigned long TmpCanonicalStateKet;
  int NbrTranslation;
  int CurrentNbrStateInOrbit;
  unsigned long DicardFlag = ~0x0ul;
  unsigned long TmpStateBra, TmpStateKet;
  this->CreatePrecalculationTable();

  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      TmpStateBra = this->ChainDescriptionBra[i];
      TmpStateKet = this->ChainDescriptionKet[i];
      this->FindCanonicalForm(TmpStateBra, TmpStateKet,TmpCanonicalStateBra, TmpCanonicalStateKet,NbrTranslation);

      if ((TmpStateBra  == TmpCanonicalStateBra) && (TmpStateKet  == TmpCanonicalStateKet ) )
	{
	  CurrentNbrStateInOrbit = this->FindNumberTranslation(TmpCanonicalStateBra,TmpCanonicalStateKet);
	  if (this->CompatibilityWithMomentum[CurrentNbrStateInOrbit] == true)
	    {
	      ++this->LargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->ChainDescriptionBra[i] = DicardFlag;
	    }
	}
      else
	{
	  this->ChainDescriptionBra[i] = DicardFlag;
	}
    }
  
  unsigned long* TmpStateDescriptionBra = new unsigned long [this->LargeHilbertSpaceDimension];
  unsigned long* TmpStateDescriptionKet = new unsigned long [this->LargeHilbertSpaceDimension];
  this->NbrStateInOrbit = new int [this->LargeHilbertSpaceDimension];

  this->LargeHilbertSpaceDimension = 0l;
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      if ( this->ChainDescriptionBra[i] != DicardFlag)
	{
	  TmpStateDescriptionBra[this->LargeHilbertSpaceDimension] = this->ChainDescriptionBra[i];
	  TmpStateDescriptionKet[this->LargeHilbertSpaceDimension] = this->ChainDescriptionKet[i];
	  CurrentNbrStateInOrbit = this->FindNumberTranslation(this->ChainDescriptionBra[i],this->ChainDescriptionKet[i]);
	  this->NbrStateInOrbit[this->LargeHilbertSpaceDimension] = CurrentNbrStateInOrbit;
	  ++this->LargeHilbertSpaceDimension;
	}
    }
  
  delete[]  this->ChainDescriptionBra;
  delete[]  this->ChainDescriptionKet;
  this->ChainDescriptionBra = TmpStateDescriptionBra;
  this->ChainDescriptionKet = TmpStateDescriptionKet;
  
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;  
  this->LookUpTable =0;

  if (this->HilbertSpaceDimension > 0)
    this->GenerateLookUpTable(memorySize);
}


 
// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

DoubledSpin0_1_2_chainWithTranslations::DoubledSpin0_1_2_chainWithTranslations (const DoubledSpin0_1_2_chainWithTranslations & chain)
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->ChainDescriptionBra = chain.ChainDescriptionBra;
      this->ChainDescriptionKet = chain.ChainDescriptionKet;
      this->DiffSz = chain.DiffSz;
      this->Momentum = chain.Momentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->UniqueStateDescriptionBra = chain.UniqueStateDescriptionBra;
      this->UniqueStateDescriptionSubArraySizeBra = chain.UniqueStateDescriptionSubArraySizeBra;
      this->NbrUniqueStateDescriptionBra = chain.NbrUniqueStateDescriptionBra;
      this->FirstIndexUniqueStateDescriptionBra = chain.FirstIndexUniqueStateDescriptionBra;

    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->ComplementaryStateShift = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescriptionBra = 0;
      this->ChainDescriptionKet = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->DiffSz = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

DoubledSpin0_1_2_chainWithTranslations::~DoubledSpin0_1_2_chainWithTranslations () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

DoubledSpin0_1_2_chainWithTranslations & DoubledSpin0_1_2_chainWithTranslations::operator = (const DoubledSpin0_1_2_chainWithTranslations & chain)
{
  AbstractDoubledSpinChainWithTranslations::operator =(chain);
/*
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->LargeHilbertSpaceDimension > 0l)
	{
 	  delete[] this->ChainDescriptionBra;
	  delete[] this->ChainDescriptionKet;
	  delete[] this->UniqueStateDescriptionBra;
	  delete[] this->UniqueStateDescriptionSubArraySizeBra;
	  delete[] this->FirstIndexUniqueStateDescriptionBra;
	  delete[] this->CompatibilityWithMomentum;
	  for (int i = 1; i <= this->ChainLength; ++i)
	    {
	      delete[] this->RescalingFactors[i];
	    } 
	  delete[] this->RescalingFactors;
	  delete[] this->NbrStateInOrbit;
	}
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->ChainDescriptionBra = chain.ChainDescriptionBra;
      this->ChainDescriptionKet = chain.ChainDescriptionKet;
      this->DiffSz = chain.DiffSz;
      this->Momentum = chain.Momentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
   }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->ComplementaryStateShift = 0;
      this->HilbertSpaceDimension = 0;
      this->ChainDescriptionBra = 0;
      this->ChainDescriptionKet = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->DiffSz = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;*/

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* DoubledSpin0_1_2_chainWithTranslations::Clone()
{
  return new DoubledSpin0_1_2_chainWithTranslations (*this);
}

// return value of twice spin projection on (Oz) for a given state
//
// stateDescription = state to which the spin projection has to be evaluated
// return value = twice spin projection on (Oz)

inline int DoubledSpin0_1_2_chainWithTranslations::GetTotalSz (unsigned long stateDescriptionBra,unsigned long stateDescriptionKet)
{
  int TmpSz = 0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      switch (stateDescriptionBra & 0x3ul)
	{
	case 0x3:
	  TmpSz += 1;
	  break;
	case 0x0:
	  TmpSz -= 1;
	  break;
	}
      stateDescriptionBra >>= 2;
      switch (stateDescriptionKet & 0x3ul)
	{
	case 0x3:
	  TmpSz -= 1;
	  break;
	case 0x0:
	  TmpSz += 1;
	  break;
	}
      stateDescriptionKet >>= 2;
    }
  return TmpSz;
}

// find the canonical form of a state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// return value = canonical form of the state

inline void DoubledSpin0_1_2_chainWithTranslations::FindCanonicalForm(unsigned long stateDescriptionBra,unsigned long stateDescriptionKet,unsigned long & canonicalStateBra , unsigned long & canonicalStateKet, int& nbrTranslation)
{
  nbrTranslation = 0;
  canonicalStateBra = stateDescriptionBra;
  canonicalStateKet = stateDescriptionKet;
  int index = 1;  
  while (index < this->ChainLength)
    {
      stateDescriptionBra = (stateDescriptionBra >> 2) | ((stateDescriptionBra & 0x3ul) << this->ComplementaryStateShift);
      stateDescriptionKet = (stateDescriptionKet >> 2) | ((stateDescriptionKet & 0x3ul) << this->ComplementaryStateShift);
      if ((stateDescriptionBra < canonicalStateBra)||((stateDescriptionBra == canonicalStateBra)&&(stateDescriptionKet < canonicalStateKet))  )
	{
	  canonicalStateBra = stateDescriptionBra;
	  canonicalStateKet = stateDescriptionKet;
	  nbrTranslation = index;
	}
      ++index;
    }
}

// find the canonical form of a state and find how many translations are needed to obtain the same state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state
// return value = canonical form of the state

inline void DoubledSpin0_1_2_chainWithTranslations::FindCanonicalForm(unsigned long stateDescriptionBra,unsigned long stateDescriptionKet, unsigned long & canonicalStateBra , unsigned long & canonicalStateKet, int& nbrTranslation, int& nbrTranslationToIdentity)
{
  nbrTranslation = 0;
  nbrTranslationToIdentity = 1;
  canonicalStateBra = stateDescriptionBra;
  canonicalStateKet = stateDescriptionKet;
  unsigned long ReferenceStateBra = stateDescriptionBra;
  unsigned long ReferenceStateKet = stateDescriptionKet;

  stateDescriptionBra = (stateDescriptionBra >> 2) | ((stateDescriptionBra & 0x3ul) << this->ComplementaryStateShift);
  stateDescriptionKet = (stateDescriptionKet >> 2) | ((stateDescriptionKet & 0x3ul) << this->ComplementaryStateShift);

  while ((ReferenceStateBra != stateDescriptionBra) && (ReferenceStateKet != stateDescriptionKet) && (nbrTranslationToIdentity < this->ChainLength))
    {
      if ((stateDescriptionBra < canonicalStateBra)||((stateDescriptionBra == canonicalStateBra)&&(stateDescriptionKet < canonicalStateKet))  )
	{
	  canonicalStateBra = stateDescriptionBra;
	  canonicalStateKet = stateDescriptionKet;
	  nbrTranslation = nbrTranslationToIdentity;
	}
      
      stateDescriptionBra = (stateDescriptionBra >> 2) | ((stateDescriptionBra & 0x3ul) << this->ComplementaryStateShift);
      stateDescriptionKet = (stateDescriptionKet >> 2) | ((stateDescriptionKet & 0x3ul) << this->ComplementaryStateShift);
      ++nbrTranslationToIdentity;
    }
}

// find how many translations are needed to obtain the same state
//
// stateDescription = unsigned integer describing the state
// return value = number of translation needed to obtain the same state

inline int DoubledSpin0_1_2_chainWithTranslations::FindNumberTranslation(unsigned long stateDescriptionBra,unsigned long stateDescriptionKet)
{
  unsigned long TmpStateBra = (stateDescriptionBra >> 2) | ((stateDescriptionBra & 0x3ul) << this->ComplementaryStateShift);
  unsigned long TmpStateKet = (stateDescriptionKet >> 2) | ((stateDescriptionKet & 0x3ul) << this->ComplementaryStateShift);
  int index = 1;  
  while ((TmpStateBra != stateDescriptionBra)||(TmpStateKet != stateDescriptionKet ))
    {     
      TmpStateBra = (TmpStateBra >> 2) | ((TmpStateBra & 0x3ul) << this->ComplementaryStateShift);
      TmpStateKet = (TmpStateKet >> 2) | ((TmpStateKet & 0x3ul) << this->ComplementaryStateShift);
      ++index;
    }
  return index;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& DoubledSpin0_1_2_chainWithTranslations::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long tmpBra,tmpKet;
  unsigned long StateDescriptionBra = this->ChainDescriptionBra[state];  
  unsigned long StateDescriptionKet = this->ChainDescriptionKet[state];  
  Str << this->FindStateIndex(StateDescriptionBra,StateDescriptionKet) << " : ";
  for (int j = this->ChainLength; j >0; j--)
    {
      tmpBra = ((StateDescriptionBra >> (( j -1) << 1)) & 0x3ul);
      tmpKet =  ((StateDescriptionKet >> ((j  - 1)<< 1)) & 0x3ul);

      Str << "(";
      if (tmpBra == 0)
	Str << "d ";
      else
	if (tmpBra == 0x1ul)
	  Str << "0 ";
	else
	  Str << "u ";
      Str << ",";
      if (tmpKet == 0)
	Str << "d ";
      else
	if (tmpKet == 0x1ul)
	  Str << "0 ";
	else
	  Str << "u ";
      Str << ") ";
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// lengthBra = length of the chain to be decided for bra spins
// lengthBra = length of the chain to be decided for ket spins
// diffSz = difference of spin projection between bra and ket chain
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored
long DoubledSpin0_1_2_chainWithTranslations::GenerateStates(int lengthBra, int lengthKet, int diffSz, long pos)
{
  if (lengthKet == 0)
    {
      if ( diffSz == -1 )
	{
	  this->ChainDescriptionKet[pos] = 0x1ul<<1;
	  pos ++;
	  return pos;
	}
      if ( diffSz == 0)
	{
	  this->ChainDescriptionKet[pos] = 0x1ul;
	  pos ++;
	  return pos;
	}
      if ( diffSz == 1)
	{
	  this->ChainDescriptionKet[pos] = 0x0ul;
	  pos ++;
	  return pos;
	}
      return pos;
    }

 if(lengthBra == 0)
    { 
      long TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz-1, pos); 
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescriptionBra[pos] = (0x1ul << 1);
	}
      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz, pos); 
      for (; pos < TmpPos; ++pos)
	this->ChainDescriptionBra[pos] = 0x1ul;
      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz+1, pos); 
      for (; pos < TmpPos; ++pos)
	this->ChainDescriptionBra[pos] = 0x0ul;
      return pos;
    }

  long TmpPos;
  unsigned long MaskBra;
  unsigned long MaskKet;
  
  if(lengthBra > 0)
    { 
      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz-1, pos); 
      MaskBra = (((0x1ul << 1)) << ((lengthBra<<1)));
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescriptionBra[pos] |= MaskBra;
	}
      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz, pos); 
      MaskBra = (((0x1ul)) << ((lengthBra<<1)));
      for (; pos < TmpPos; ++pos)
	this->ChainDescriptionBra[pos] |= MaskBra;
      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz+1, pos); 
      MaskBra = (((0x0ul)) << ((lengthBra<<1)));
      for (; pos < TmpPos; ++pos)
	this->ChainDescriptionBra[pos] |= MaskBra;
      return pos;
    }
  
  if (lengthKet > 0)
    {
      TmpPos = this->GenerateStates(lengthBra,lengthKet-1, diffSz+1, pos); 
      MaskKet = (((0x1ul << 1)) << (lengthKet<<1));
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescriptionKet[pos] |= MaskKet;
	}
      TmpPos = this->GenerateStates(lengthBra,lengthKet-1, diffSz, pos); 
      MaskKet = ((0x1ul) << (lengthKet<<1));
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescriptionKet[pos] |= MaskKet;
	}
      TmpPos = this->GenerateStates(lengthBra,lengthKet-1, diffSz-1, pos); 
      MaskKet = ((0x0ul) << (lengthKet<<1));
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescriptionKet[pos] |= MaskKet;
	}
      return pos;
    }
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// nbrNUp = number of particles with quantum number up
// nbrNDown = number of particles with quantum number down
// return value = Hilbert space dimension

long DoubledSpin0_1_2_chainWithTranslations::ShiftedEvaluateHilbertSpaceDimension(int lengthBra, int lengthKet, int diffSz)
{
  if ((lengthBra < 0) || (lengthKet < 0))
    return 0;
  
  if ((lengthBra == 0) && (lengthKet == 0))
    {
      if (diffSz == 0) 
	{
	  return 3;
	}
      if (diffSz == 1) 
	{
	  return 2;
	}
      if (diffSz == -1) 
	{
	  return 2;
	}

      if (diffSz == 2) 
	{
	  return 1;
	}
      if (diffSz == -2) 
	{
	  return 1;
	}

    }  
  long Tmp=0;
  
  if (lengthBra == 0)
    {
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra,lengthKet-1, diffSz-1);   
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra,lengthKet-1, diffSz); 
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra,lengthKet-1, diffSz+1); 
      return Tmp;
    }

  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra-1,lengthKet, diffSz+1); 


  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(lengthBra-1,lengthKet, diffSz); 

  Tmp +=  this->ShiftedEvaluateHilbertSpaceDimension(lengthBra-1,lengthKet, diffSz-1);
  return Tmp;
}
