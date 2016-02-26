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


#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslations.h"


#include <iostream>
#include <math.h>

using std::cout;
using std::endl;


// default constructor
//

DoubledSpin0_1_2_ChainWithTranslations::DoubledSpin0_1_2_ChainWithTranslations () 
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
  this->ShiftNegativeDiffSz = 0;
  this->BraShiftNegativeSz = 0;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}



// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

DoubledSpin0_1_2_ChainWithTranslations::DoubledSpin0_1_2_ChainWithTranslations (int chainLength, int diffSz, int memorySize, int memorySlice) 
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
  this->ShiftNegativeDiffSz =   this->LargeHilbertSpaceDimension;
  if (this->DiffSz !=0 )
    this->LargeHilbertSpaceDimension += this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, -this->DiffSz);
  
  this->ChainDescriptionBra = new unsigned long [this->LargeHilbertSpaceDimension];
  this->ChainDescriptionKet = new unsigned long [this->LargeHilbertSpaceDimension];
  
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, this->DiffSz, 0l);
  if (this->DiffSz != 0)
    TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, -this->DiffSz, TmpHilbertSpaceDimension);
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in DoubledSpin0_1_2_ChainWithTranslations!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  
/*  for(int i = 0 ; i <this->HilbertSpaceDimension ;i++)
    cout <<i<<" "<<this->ChainDescriptionBra[i]<<" "<< this->ChainDescriptionKet[i]<<" "<<this->GetTotalSz (this->ChainDescriptionBra[i],  this->ChainDescriptionKet[i])<<endl;*/

  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memorySize);
    }
  this->RescalingFactors = 0;
//  cout <<"BraShiftNegativeSz = "<<BraShiftNegativeSz<<endl;
}
 


// constructor for Hilbert space corresponding to a given total spin projection Sz no contraint on Momentum
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

DoubledSpin0_1_2_ChainWithTranslations::DoubledSpin0_1_2_ChainWithTranslations (int chainLength, int momentum, int diffSz, int memorySize, int memorySlice) 
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
  this->ShiftNegativeDiffSz =   this->LargeHilbertSpaceDimension;
  if (this->DiffSz !=0 )
    this->LargeHilbertSpaceDimension += this->ShiftedEvaluateHilbertSpaceDimension(this->ChainLength-1, this->ChainLength-1, -this->DiffSz);

  this->ChainDescriptionBra = new unsigned long [this->LargeHilbertSpaceDimension];
  this->ChainDescriptionKet = new unsigned long [this->LargeHilbertSpaceDimension];

  
  long TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, this->DiffSz, 0l);
  if (this->DiffSz != 0)
    TmpHilbertSpaceDimension = GenerateStates(this->ChainLength-1, this->ChainLength-1, -this->DiffSz, TmpHilbertSpaceDimension);
  
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

DoubledSpin0_1_2_ChainWithTranslations::DoubledSpin0_1_2_ChainWithTranslations (const DoubledSpin0_1_2_ChainWithTranslations & chain)
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
      this->ShiftNegativeDiffSz = chain.ShiftNegativeDiffSz;
      this->BraShiftNegativeSz =  chain.BraShiftNegativeSz;
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
      this->ShiftNegativeDiffSz =0;
      this->BraShiftNegativeSz =0;
    }
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

DoubledSpin0_1_2_ChainWithTranslations::~DoubledSpin0_1_2_ChainWithTranslations () 
{
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

DoubledSpin0_1_2_ChainWithTranslations & DoubledSpin0_1_2_ChainWithTranslations::operator = (const DoubledSpin0_1_2_ChainWithTranslations & chain)
{
  AbstractDoubledSpinChainWithTranslations::operator =(chain);
  this->ShiftNegativeDiffSz = chain.ShiftNegativeDiffSz;
  this->BraShiftNegativeSz = chain.BraShiftNegativeSz;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* DoubledSpin0_1_2_ChainWithTranslations::Clone()
{
  return new DoubledSpin0_1_2_ChainWithTranslations (*this);
}

// return value of twice spin projection on (Oz) for a given state
//
// stateDescription = state to which the spin projection has to be evaluated
// return value = twice spin projection on (Oz)

inline int DoubledSpin0_1_2_ChainWithTranslations::GetTotalSz (unsigned long stateDescriptionBra, unsigned long stateDescriptionKet)
{
  int TmpSz = 0;
  int Sign;
  for (int i = 0; i < this->ChainLength; i++)
    {
      if (i % 2==0 ) 
	Sign = 1;
      else
	Sign = -1;
      
      switch (stateDescriptionBra & 0x3ul)
	{
	  
	case 0x2:
	  TmpSz += Sign;
	  break;
	case 0x0:
	  TmpSz -= Sign;
	  break;
	}
      stateDescriptionBra >>= 2;
      switch (stateDescriptionKet & 0x3ul)
	{
	case 0x2:
	  TmpSz -= Sign;
	  break;
	case 0x0:
	  TmpSz += Sign;
	  break;
	}
      stateDescriptionKet >>= 2;
    }
  return TmpSz;
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& DoubledSpin0_1_2_ChainWithTranslations::PrintState (ostream& Str, int state)
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
long DoubledSpin0_1_2_ChainWithTranslations::GenerateStates(int lengthBra, int lengthKet, int diffSz, long pos)
{

  if ((lengthKet == 0) && (lengthBra == 0))
    {

      if (diffSz == 0) 
	{
	  this->ChainDescriptionBra[pos] = 0x1ul<<1;
	  this->ChainDescriptionKet[pos] = 0x1ul<<1;
	  pos ++;
	  this->ChainDescriptionBra[pos] = 0x1ul;
	  this->ChainDescriptionKet[pos] = 0x1ul;
	  pos ++;
	  this->ChainDescriptionBra[pos] = 0x0ul;
	  this->ChainDescriptionKet[pos] = 0x0ul;
	  pos ++;
	  return pos;
	}
      if (diffSz == 1) 
	{
	  this->ChainDescriptionBra[pos] = 0x1ul<<1;
	  this->ChainDescriptionKet[pos] = 0x1ul;
	  pos ++;
	  this->ChainDescriptionBra[pos] = 0x1ul;
	  this->ChainDescriptionKet[pos] = 0x0ul;
	  pos ++;
	  return pos;
	}
      if (diffSz == -1) 
	{
	  this->ChainDescriptionBra[pos] = 0x1ul;
	  this->ChainDescriptionKet[pos] = 0x1ul<<1;
	  pos ++;
	  this->ChainDescriptionBra[pos] = 0x0ul;
	  this->ChainDescriptionKet[pos] = 0x1ul;
	  pos ++;
	  return pos;
	}

      if (diffSz == 2) 
	{
	  this->ChainDescriptionBra[pos] = 0x1ul<<1;
	  this->ChainDescriptionKet[pos] = 0x0ul;
	  pos ++;
	  return pos;
	}
      if (diffSz == -2) 
	{
	  this->ChainDescriptionBra[pos] = 0x0ul;
	  this->ChainDescriptionKet[pos] = 0x1ul<<1;
	  pos ++;
	  return pos;
	}
   
      return pos;
    }

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

      int Sign;
      if (lengthBra%2==0)
	Sign = -1;
      else
	{
	  Sign = 1;
	}

      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz+Sign, pos); 
      MaskBra = (((0x1ul << 1)) << ((lengthBra<<1)));
      for (; pos < TmpPos; ++pos)
	{
	  this->ChainDescriptionBra[pos] |= MaskBra;
	}
      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz, pos); 
      MaskBra = (((0x1ul)) << ((lengthBra<<1)));
      for (; pos < TmpPos; ++pos)
	this->ChainDescriptionBra[pos] |= MaskBra;
      TmpPos = this->GenerateStates(lengthBra-1,lengthKet, diffSz-Sign, pos); 
      MaskBra = (((0x0ul)) << ((lengthBra<<1)));
      for (; pos < TmpPos; ++pos)
	this->ChainDescriptionBra[pos] |= MaskBra;
      return pos;
    }
  
  if (lengthKet > 0)
    {
      
      int Sign;
      if (lengthKet%2==0)
	Sign = -1;
      else
	{
	  Sign = 1;
	}
      
      TmpPos = this->GenerateStates(lengthBra,lengthKet-1, diffSz-Sign, pos); 
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
      TmpPos = this->GenerateStates(lengthBra,lengthKet-1, diffSz+Sign, pos); 
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

long DoubledSpin0_1_2_ChainWithTranslations::ShiftedEvaluateHilbertSpaceDimension(int lengthBra, int lengthKet, int diffSz)
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
      return 0;
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



// find state index
//
// state = state description
// return value = corresponding index

int DoubledSpin0_1_2_ChainWithTranslations::FindStateIndex(unsigned long stateBra,unsigned long stateKet) 
{
  int PosMin = 0;
  int PosMax = this->BraShiftNegativeSz - 1;
  if( this->GetTotalSz (stateBra, stateKet) <= 0)
    {
      PosMin = this->BraShiftNegativeSz;
      PosMax = this->NbrUniqueStateDescriptionBra - 1;     
    }

  int PosMid = (PosMin + PosMax) >> 1; 
  unsigned long CurrentState = this->UniqueStateDescriptionBra[PosMid];
  while (( (PosMax  - PosMin ) > 1 ) && (CurrentState != stateBra))
    {
//      cout <<PosMin<<" "<<PosMax<<" "<<PosMid<<endl;
       if (CurrentState > stateBra)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->UniqueStateDescriptionBra[PosMid];
    }
  
  if (CurrentState != stateBra)
    PosMid = PosMax;

  if (this->UniqueStateDescriptionBra[PosMid] != stateBra)
    {
      return this->HilbertSpaceDimension;
    }

  PosMin = this->FirstIndexUniqueStateDescriptionBra[PosMid];
  PosMax = PosMin + this->UniqueStateDescriptionSubArraySizeBra[PosMid] - 1;
  PosMid = (PosMin + PosMax) >> 1;
  CurrentState = this->ChainDescriptionKet[PosMid];
  while ( ( (PosMax  - PosMin ) > 1) && (CurrentState != stateKet))
    {
      if (CurrentState > stateKet)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->ChainDescriptionKet[PosMid];
    }
  if (this->ChainDescriptionKet[PosMax] == stateKet)
    return PosMax;
  if (this->ChainDescriptionKet[PosMid] == stateKet)
    return PosMid;
  return this->HilbertSpaceDimension;
}



// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void  DoubledSpin0_1_2_ChainWithTranslations::GenerateLookUpTable(unsigned long memory)
{  
  long TmpUniquePartition = 1l;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->ChainDescriptionBra[i - 1] == this->ChainDescriptionBra[i]))
	{
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	++TmpUniquePartition;
    }
  
  this->NbrUniqueStateDescriptionBra = TmpUniquePartition;
  this->UniqueStateDescriptionBra = new unsigned long [this->NbrUniqueStateDescriptionBra];
  this->UniqueStateDescriptionSubArraySizeBra = new int [this->NbrUniqueStateDescriptionBra];
  this->FirstIndexUniqueStateDescriptionBra = new int [this->NbrUniqueStateDescriptionBra];
  TmpUniquePartition = 0l;
  this->UniqueStateDescriptionBra[0l] = this->ChainDescriptionBra[0l];
  this->UniqueStateDescriptionSubArraySizeBra[0] = 1;
  this->FirstIndexUniqueStateDescriptionBra[0] = 0;
  this->BraShiftNegativeSz = 0;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->ChainDescriptionBra[i - 1] == this->ChainDescriptionBra[i]))
	{
	  ++this->UniqueStateDescriptionSubArraySizeBra[TmpUniquePartition];
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	{
	  ++TmpUniquePartition;
	  this->UniqueStateDescriptionBra[TmpUniquePartition] = this->ChainDescriptionBra[i];
	  this->UniqueStateDescriptionSubArraySizeBra[TmpUniquePartition] = 1; 
	  this->FirstIndexUniqueStateDescriptionBra[TmpUniquePartition] = i;
	  if(this->ChainDescriptionBra[i - 1]<  this->ChainDescriptionBra[i])
	    {
	      this->BraShiftNegativeSz = TmpUniquePartition;
	    }
	}
    }
}
