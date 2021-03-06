////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of spin 1/2 chain with translastion invariance          //
//                                                                            //
//                        last modification : 29/01/2002                      //
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

#include "HilbertSpace/Spin1_2Chain.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "MathTools/Complex.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "GeneralTools/ArrayTools.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;
using std::hex;
using std::dec;


#define M_SQRT3 1.73205080756888
#define M_SQRT6 2.44948974278318


// default constructor
//

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations () 
{
  this->Flag.Initialize();
  this->LookUpTable = 0;
  this->LookUpTableShift = 0;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->ChainLength = 0;
  this->Momentum = 0;
  this->StateMask = 0x0ul;
  this->StateShift = 0;
  this->ComplementaryStateShift = 0;
  this->Sz = 0;
  this->FixedSpinProjectionFlag = false;
  this->CompatibilityWithMomentum = 0;
  this->RescalingFactors = 0;
  this->NbrStateInOrbit = 0;
  this->MaxXMomentum = 0;
}

// constructor for Hilbert space with no restriction on total spin projection Sz
//
// chainLength = number of spin 1
// momemtum = total momentum of each state
// translationStep = indicates the step for an elementary translation
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations (int chainLength, int momentum, int translationStep, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->FixedSpinProjectionFlag = false;
  this->ComplementaryStateShift = this->ChainLength - translationStep;
  this->StateMask = (0x1ul << translationStep) - 1ul;
  this->StateShift = translationStep;
  this->Momentum = momentum;
  this->MaxXMomentum = chainLength;

  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;

  this->CreatePrecalculationTable();

  cout << "warning : untested code" << endl;
  long TmpHilbertSpaceDimension = (1l <<  this->ChainLength );
  this->LargeHilbertSpaceDimension = 0l;
  unsigned long TmpState;
  unsigned long TmpState2;
  int NbrTranslation;
  int CurrentNbrStateInOrbit;
  for (long i = TmpHilbertSpaceDimension - 1l; i >= 0l; --i)
    {
      TmpState = (unsigned long) i;
      TmpState2 = this->FindCanonicalForm(TmpState, NbrTranslation);
      if (TmpState2 == TmpState)
	{
	  CurrentNbrStateInOrbit = this->FindNumberTranslation(TmpState2);
	  if (this->CompatibilityWithMomentum[CurrentNbrStateInOrbit] == true)
	    {
	      ++this->LargeHilbertSpaceDimension;
	    }
	}
    }
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->NbrStateInOrbit = new int [this->LargeHilbertSpaceDimension];
  
  this->LargeHilbertSpaceDimension = 0l;
  for (long i = TmpHilbertSpaceDimension - 1l; i >= 0l; --i)
    {
      TmpState = (unsigned long) i;
      TmpState2 = this->FindCanonicalForm(TmpState, NbrTranslation);
      if (TmpState2 == TmpState)
	{
	  CurrentNbrStateInOrbit = this->FindNumberTranslation(TmpState2);
	  if (this->CompatibilityWithMomentum[CurrentNbrStateInOrbit] == true)
	    {
	      this->StateDescription[this->LargeHilbertSpaceDimension] = TmpState2;
	      this->NbrStateInOrbit[this->LargeHilbertSpaceDimension] = CurrentNbrStateInOrbit;
	      ++this->LargeHilbertSpaceDimension;
	    }
	}
    }
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;  

  if (this->HilbertSpaceDimension > 0)
    this->CreateLookUpTable();
}

// constructor for Hilbert space corresponding to a given total spin projection Sz
//
// chainLength = number of spin 1/2
// momemtum = total momentum of each state
// translationStep = indicates the step for an elementary translation
// sz = twice the value of total Sz component
// memorySize = memory size in bytes allowed for look-up table
// memorySlice = maximum amount of memory that can be allocated to partially evalauted the states

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations (int chainLength, int momentum, int translationStep, int sz, int memorySize, int memorySlice) 
{
  this->Flag.Initialize();
  this->ChainLength = chainLength;
  this->Sz = sz;
  this->FixedSpinProjectionFlag = true;
  this->ComplementaryStateShift = this->ChainLength - translationStep;
  this->StateMask = (0x1ul << translationStep) - 1ul;
  this->StateShift = translationStep;
  this->Momentum = momentum;
  this->MaxXMomentum = chainLength;

  memorySize /= sizeof(long);
  this->LookUpTableShift = 1;
  while ((1 << this->LookUpTableShift) <= memorySize)
    ++this->LookUpTableShift;
  if (this->LookUpTableShift < (this->ChainLength << 1))
    this->LookUpTableShift = (this->ChainLength << 1) - this->LookUpTableShift + 1;
  else
    this->LookUpTableShift = 0;

  this->CreatePrecalculationTable();

  this->StateDescription = new unsigned long [this->EvaluateHilbertSpaceDimension(this->ChainLength, this->Sz)];
  long TmpHilbertSpaceDimension = this->GenerateStates(0l, this->ChainLength - 1, (this->ChainLength + this->Sz) >> 1);
  this->LargeHilbertSpaceDimension = 0l;
  unsigned long TmpState;
  unsigned long TmpState2;
  int NbrTranslation;
  int CurrentNbrStateInOrbit;
  unsigned long DicardFlag = ~0x0ul;
  for (long i = 0l; i < TmpHilbertSpaceDimension; ++i)
    {
      TmpState = this->StateDescription[i];
      TmpState2 = this->FindCanonicalForm(TmpState, NbrTranslation);
      if (TmpState2 == TmpState)
	{
	  CurrentNbrStateInOrbit = this->FindNumberTranslation(TmpState2);
	  if (this->CompatibilityWithMomentum[CurrentNbrStateInOrbit] == true)
	    {
	      ++this->LargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->StateDescription[i] = DicardFlag;
	    }
	}
      else
	{
	  this->StateDescription[i] = DicardFlag;
	}
    }
  
  unsigned long* TmpStateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->NbrStateInOrbit = new int [this->LargeHilbertSpaceDimension];
  
  this->LargeHilbertSpaceDimension = 0l;
  for (long i = TmpHilbertSpaceDimension - 1l; i >= 0l; --i)
    {
      if (this->StateDescription[i] != DicardFlag)
	{
	  TmpStateDescription[this->LargeHilbertSpaceDimension] = this->StateDescription[i];
	  CurrentNbrStateInOrbit = this->FindNumberTranslation(this->StateDescription[i]);
	  this->NbrStateInOrbit[this->LargeHilbertSpaceDimension] = CurrentNbrStateInOrbit;
	  ++this->LargeHilbertSpaceDimension;
	}
    }
  
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;  
  this->LookUpTable =0;
  if (this->HilbertSpaceDimension > 0)
    this->CreateLookUpTable();
}

// constructor from pre-constructed datas
//
// hilbertSpaceDimension = Hilbert space dimension
// chainDescription = array describing states
// chainLength = number of spin 1
// momemtum = total momentum of each state
// sz = twice the value of total Sz component
// fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
// lookUpTableShift = shift to apply to a state to obtain an index to the look-up table 
// complementaryStateShift = shift to apply to move the spin from one end to the other one

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations (int hilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
							    int momentum, int sz, bool fixedQuantumNumberFlag, int lookUpTableShift, 
							    int complementaryStateShift)
{
  this->Flag.Initialize();
  this->ComplementaryStateShift = complementaryStateShift;
  this->LookUpTableShift = lookUpTableShift;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->StateDescription = chainDescription;
  this->Momentum = momentum;
  this->MaxXMomentum = chainLength;
  this->Sz = sz;
  this->FixedSpinProjectionFlag = fixedQuantumNumberFlag;
  this->ChainLength = chainLength;
  this->CreatePrecalculationTable();
  this->CreateLookUpTable();
}
  
// copy constructor (without duplicating datas)
//
// chain = reference on chain to copy

Spin1_2ChainWithTranslations::Spin1_2ChainWithTranslations (const Spin1_2ChainWithTranslations& chain) 
{
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->Momentum = chain.Momentum;
      this->MaxXMomentum = chain.MaxXMomentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->StateMask = chain.StateMask;
      this->StateShift = chain.StateShift;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
    }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->ComplementaryStateShift = 0;
      this->HilbertSpaceDimension = 0;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
      this->StateMask = 0x0ul;
      this->StateShift = 0;
      this->ComplementaryStateShift = 0;
      this->MaxXMomentum = 0;
    }
}

// destructor
//

Spin1_2ChainWithTranslations::~Spin1_2ChainWithTranslations () 
{
  if ((this->ChainLength != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTable;
      delete[] this->CompatibilityWithMomentum;
      int TmpPeriodicity = this->ChainLength / this->StateShift;
      for (int i = 1; i <= TmpPeriodicity; ++i)
	{
	  delete[] this->RescalingFactors[i];
	} 
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
}

// assignement (without duplicating datas)
//
// chain = reference on chain to copy
// return value = reference on current chain

Spin1_2ChainWithTranslations& Spin1_2ChainWithTranslations::operator = (const Spin1_2ChainWithTranslations& chain)
{
  if ((this->ChainLength != 0) && (this->HilbertSpaceDimension > 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->LookUpTable;
      delete[] this->CompatibilityWithMomentum;
      int TmpPeriodicity = this->ChainLength / this->StateShift;
      for (int i = 1; i <= TmpPeriodicity; ++i)
	{
	  delete[] this->RescalingFactors[i];
	} 
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }  
  this->Flag = chain.Flag;
  if (chain.ChainLength != 0)
    {
      this->ChainLength = chain.ChainLength;
      this->HilbertSpaceDimension = chain.HilbertSpaceDimension;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
      this->LookUpTable = chain.LookUpTable;
      this->LookUpTableShift = chain.LookUpTableShift;
      this->StateDescription = chain.StateDescription;
      this->Sz = chain.Sz;
      this->Momentum = chain.Momentum;
      this->MaxXMomentum = chain.MaxXMomentum;
      this->FixedSpinProjectionFlag = chain.FixedSpinProjectionFlag;
      this->CompatibilityWithMomentum = chain.CompatibilityWithMomentum;
      this->RescalingFactors = chain.RescalingFactors;
      this->NbrStateInOrbit = chain.NbrStateInOrbit;
      this->StateMask = chain.StateMask;
      this->StateShift = chain.StateShift;
      this->ComplementaryStateShift = chain.ComplementaryStateShift;
   }
  else
    {
      this->LookUpTable = 0;
      this->LookUpTableShift = 0;
      this->ComplementaryStateShift = 0;
      this->HilbertSpaceDimension = 0;
      this->StateDescription = 0;
      this->ChainLength = 0;
      this->Momentum = 0;
      this->MaxXMomentum = 0;
      this->Sz = 0;
      this->FixedSpinProjectionFlag = false;
      this->CompatibilityWithMomentum = 0;
      this->RescalingFactors = 0;
      this->NbrStateInOrbit = 0;
      this->StateMask = 0x0ul;
      this->StateShift = 0;
      this->ComplementaryStateShift = 0;
    }
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* Spin1_2ChainWithTranslations::Clone()
{
  return new Spin1_2ChainWithTranslations (*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> Spin1_2ChainWithTranslations::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  if (this->FixedSpinProjectionFlag == true)
    {
      List<AbstractQuantumNumber*> TmpList2;
      TmpList2 += new PeriodicMomentumQuantumNumber (this->Momentum, this->ChainLength);
      TmpList2 += new SzQuantumNumber (this->Sz);
      TmpList += new VectorQuantumNumber (TmpList2);
    }
  else
    {
      int TmpSz = - 2 * this->ChainLength;
      for (int i = 0; i <= (2 * this->ChainLength); i++)
	{
	  List<AbstractQuantumNumber*> TmpList2;
	  TmpList2 += new PeriodicMomentumQuantumNumber (this->Momentum, this->ChainLength);
	  TmpList2 += new SzQuantumNumber (TmpSz);
	  TmpList += new VectorQuantumNumber (TmpList2);
	  TmpSz += 2;
	}
    }
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* Spin1_2ChainWithTranslations::GetQuantumNumber (int index)
{ 
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->Momentum, this->ChainLength);
  TmpList += new SzQuantumNumber (this->TotalSz(index));
  return new VectorQuantumNumber (TmpList);
}

// return value of twice spin projection on (Oz) for a given state
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

int Spin1_2ChainWithTranslations::TotalSz (int index)
{
  if (this->FixedSpinProjectionFlag == true)
    return this->Sz;
  unsigned long State = this->StateDescription[index];
  int TmpSz = 0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpSz += (State & 0x1ul) << 1;
      State >>= 1;
    }
  TmpSz -= this->ChainLength;
  return TmpSz;
}

// return value of the value of the sum of the square of spin projection on (Oz) 
//
// index = index of the state to test
// return value = twice spin projection on (Oz)

inline double Spin1_2ChainWithTranslations::TotalSzSz (int index)
{  
  return (((double) this->ChainLength) * 0.25);
}

// return value of twice spin projection on (Oz) for a given state
//
// stateDescription = state to which the spin projection has to be evaluated
// return value = twice spin projection on (Oz)

inline int Spin1_2ChainWithTranslations::GetTotalSz (unsigned long stateDescription)
{
  int TmpSz = 0;
  for (int i = 0; i < this->ChainLength; i++)
    {
      TmpSz += (stateDescription & 0x1ul);
      stateDescription >>= 1;
    }
  return ((2 * TmpSz) - this->ChainLength);
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

double Spin1_2ChainWithTranslations::SziSzj (int i, int j, int state)
{  
  unsigned long Mask = ((0x1ul << i) | (0x1ul << j));
  unsigned long tmpState = this->StateDescription[state] & Mask;
  if ((tmpState == 0x0ul) || (tmpState == Mask))
    return 0.25;
  else
    return -0.25;
}

// compute the parity (prod_i Sz_i) for a given state
//
// state = index of the state to be applied on Sz_i operator
// return value = 0 if prod_i Sz_i = 1, 1 if prod_i Sz_i = -1

unsigned long Spin1_2ChainWithTranslations::Parity (int state)
{
  unsigned long TmpState = this->StateDescription[state];
#ifdef __64_BITS__
  TmpState ^= TmpState >> 32;
#endif
  TmpState ^= TmpState >> 16;
  TmpState ^= TmpState >> 8;
  TmpState ^= TmpState >> 4;
  TmpState ^= TmpState >> 2;
  TmpState ^= TmpState >> 1;
  return (TmpState & 0x1ul);
  
}
// return index of resulting state from application of P_ij operator on a given state
//
// i = first position
// j = second position
// state = index of the state to be applied on P_ij operator
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::Pij (int i, int j, int state, double& coefficient, int& nbrTranslation)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long tmpMask = (0x1ul << i) | (0x1ul << j);
  unsigned long tmpState2 = tmpState & tmpMask;
  unsigned long tmpState3 = ~tmpState & tmpMask;
  if ((tmpState2 == 0x0ul) || (tmpState3 == 0x0ul))
    return this->HilbertSpaceDimension;
  else
    {
      tmpState &= ~tmpMask;
      tmpState |= tmpState3;
      coefficient = 1.0;
      return this->SymmetrizeResult(tmpState, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
//       tmpState = this->FindCanonicalForm((tmpState & ~tmpMask) | tmpState3, nbrTranslation, i);
//       if (this->CompatibilityWithMomentum[i] == false)
// 	return this->HilbertSpaceDimension;
//       j = this->FindStateIndex(tmpState);
//       coefficient *= this->RescalingFactors[this->NbrStateInOrbit[state]][i];
//       return j;
    }
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{  
  unsigned long tmpState = this->StateDescription[state];
  unsigned long State = tmpState;
  unsigned long tmpState2 = tmpState;
  tmpState >>= i;
  tmpState &= 0x1ul;
  if (i != j)
    {
      tmpState2 >>= j; 
      tmpState2 &= 0x1ul;
      tmpState2 <<= 1;
      tmpState2 |= tmpState;
      if (tmpState2 == 0x1ul)
	{
	  State |= (0x1ul << j);
	  State &= ~(0x1ul << i);
	  coefficient = 1.0;
	  return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
	}
      else
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  if (tmpState == 0)
    {
      coefficient = -0.25;
      return state;
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i operator on a given state (only valid if there is no constraint on total Sz)
//
// i = operator position
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::Spi (int i, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long State = this->StateDescription[state];
  unsigned long TmpMask = (0x1ul << i);
  if ((State & TmpMask) == 0)
    {
      State |= TmpMask;
      coefficient = 1.0;
      return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i operator on a given state (only valid if there is no constraint on total Sz)
//
// i = operator position
// state = index of the state to be applied on S+_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::Smi (int i, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long State = this->StateDescription[state];
  unsigned long TmpMask = (0x1ul << i);
  if ((State & TmpMask) != 0)
    {
      State &= ~TmpMask;
      coefficient = 1.0;
      return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
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

int Spin1_2ChainWithTranslations::SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  if (i == j)
    return this->HilbertSpaceDimension;
  unsigned long State = this->StateDescription[state];
  unsigned long TmpMask = (0x1ul << j) | (0x1ul << i);
  if ((State & TmpMask) == 0)
    {
      State |= TmpMask;
      coefficient = 1.0;
      return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
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

int Spin1_2ChainWithTranslations::SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  if (i == j)
    return this->HilbertSpaceDimension;
  unsigned long State = this->StateDescription[state];
  unsigned long TmpMask = (0x1ul << j) | (0x1ul << i);
  if ((State & TmpMask) == TmpMask)
    {
      State &= ~TmpMask;
      coefficient = 1.0;
      return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
    }
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S+_i S+_i operator on a given state
//
// i = position of first S+ operator
// state = index of the state to be applied on S+_i S+_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::SpiSpi (int i, int state, double& coefficient, int& nbrTranslation)
{
  return this->HilbertSpaceDimension;
}

// return index of resulting state from application of S-_i S-_i operator on a given state
//
// i = position of the S- operator
// state = index of the state to be applied on S-_i S-_i operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::SmiSmi (int i, int state, double& coefficient, int& nbrTranslation)
{
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

int Spin1_2ChainWithTranslations::SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long State = this->StateDescription[state];
  if (((State >> j) & (unsigned long) 0x1) == 0)
    coefficient = 0.5;
  else
    coefficient = -0.5;
  if ((State & (0x1ul << i)) != 0)
    return this->HilbertSpaceDimension;
  State |= 0x1ul << i;
  return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
}

// return index of resulting state from application of S-_i Sz_j operator on a given state
//
// i = position of S- operator
// j = position of Sz operator
// state = index of the state to be applied on S-_i Sz_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

int Spin1_2ChainWithTranslations::SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{
  unsigned long State = this->StateDescription[state];
  if (((State >> j) & (unsigned long) 0x1) == 0)
    coefficient = 0.5;
  else
    coefficient = -0.5;
  if ((State & (0x1ul << i)) == 0)
    return this->HilbertSpaceDimension;
  State &= ~(0x1ul << i);
  return this->SymmetrizeResult(State, this->NbrStateInOrbit[state], coefficient, nbrTranslation);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* Spin1_2ChainWithTranslations::ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter)
{
  if (((VectorQuantumNumber&) q).GetQuantumNumbers().GetNbrElement() != 2)
    return 0;
  if (((((VectorQuantumNumber&) q)[0]))->GetQuantumNumberType() != AbstractQuantumNumber::PeriodicMomentum)
    return 0;
  if (((((VectorQuantumNumber&) q)[1]))->GetQuantumNumberType() != AbstractQuantumNumber::Sz)
    return 0;
  if (this->Momentum != ((PeriodicMomentumQuantumNumber*) (((VectorQuantumNumber&) q)[0]))->GetMomentum())
    return 0;
  if (this->FixedSpinProjectionFlag == true)
    {
      if (this->Sz != ((SzQuantumNumber*) (((VectorQuantumNumber&) q)[1]))->GetSz())
	return 0;
      else
	return this;
    }
  int TmpSz = ((SzQuantumNumber*) (((VectorQuantumNumber&) q)[1]))->GetSz();
  if ((TmpSz < (-2 * this->ChainLength)) || (TmpSz > (2 * this->ChainLength)))
    return 0;
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
  SubspaceDescription[0] = this->StateDescription[TmpConvArray[0]];
  ConvArray[0] = TmpConvArray[0];
//  unsigned long TestMask = this->LookUpTableMask;
  for (int i = 1; i < HilbertSubspaceDimension; i++)
    {
      SubspaceDescription[i] = this->StateDescription[TmpConvArray[i]];
      ConvArray[i] = TmpConvArray[i];
    }
  converter = SubspaceSpaceConverter (this->HilbertSpaceDimension, HilbertSubspaceDimension, ConvArray);
  return new Spin1_2ChainWithTranslations (HilbertSubspaceDimension, SubspaceDescription, this->ChainLength,
					 this->Momentum, TmpSz, true, 
					 this->LookUpTableShift, this->ComplementaryStateShift);
}


// find state index
//
// state = state description
// return value = corresponding index

inline int Spin1_2ChainWithTranslations::FindStateIndex(unsigned long state)
{
  unsigned long MidPos = state >> this->LookUpTableShift;
  unsigned long LowPos = this->LookUpTable[MidPos];
  unsigned long HighPos = this->LookUpTable[MidPos + 1];
  while ((HighPos - LowPos) > 1)
    {
      MidPos = (HighPos + LowPos) >> 1;
      if (this->StateDescription[MidPos] >= state)
	HighPos = MidPos;
      else
	LowPos = MidPos;
    }

  if (this->StateDescription[LowPos] == state) 
    return LowPos;
  if (this->StateDescription[HighPos] == state) 
    return HighPos;   
  return this->HilbertSpaceDimension;
}



// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& Spin1_2ChainWithTranslations::PrintState (ostream& Str, int state)
{
  if (state >= this->HilbertSpaceDimension)    
    return Str;
  unsigned long Mask = 0x1ul;
  for (int k = 0; k < this->ChainLength; k++)    
    {
      if ((this->StateDescription[state] & Mask) == 0x0ul)
	Str << "- ";
      else
	Str << "+ ";
      Mask <<= 1;
    }
  Str << " " << hex << this->StateDescription[state] << dec;
  return Str;
}

// generate all states with constraint on total Sz
//
// position = current position in the Hilbert space basis
// sitePosition = largest sit position that has to be filled 
// currentNbrSpinUp = number of spin to put in each state
// return value = new current position in the Hilbert space basis

long Spin1_2ChainWithTranslations::GenerateStates(long position, int sitePosition, int currentNbrSpinUp) 
{
  if ((currentNbrSpinUp > (sitePosition + 1)) || (currentNbrSpinUp < 0) || (sitePosition < 0)) 
    return position;
  if (currentNbrSpinUp == (sitePosition + 1))
    {
      this->StateDescription[position] = (0x1ul << currentNbrSpinUp) - 0x1ul;
      return position + 1;
    }
  if (currentNbrSpinUp == 0)
    {
      this->StateDescription[position] = 0x0ul;
      return position + 1;
    }
  unsigned long Mask = 0x1ul << sitePosition;
  long TmpPosition = this->GenerateStates(position, sitePosition  - 1, currentNbrSpinUp - 1);
  for (; position < TmpPosition; ++position)
    this->StateDescription[position] |= Mask;
  return this->GenerateStates(position, sitePosition - 1, currentNbrSpinUp);
}

// create precalculation tables
//

void Spin1_2ChainWithTranslations::CreatePrecalculationTable()
{
  int TmpPeriodicity = this->ChainLength / this->StateShift;
  this->CompatibilityWithMomentum = new bool [TmpPeriodicity + 1];
  for (int i = 0; i <= TmpPeriodicity; ++i)
    if (((i * this->Momentum) % TmpPeriodicity) == 0)
      this->CompatibilityWithMomentum[i] = true;
    else
      this->CompatibilityWithMomentum[i] = false;

  this->RescalingFactors = new double* [TmpPeriodicity + 1];
  for (int i = 1; i <= TmpPeriodicity; ++i)
    {
      this->RescalingFactors[i] = new double [TmpPeriodicity + 1];
      for (int j = 1; j <= TmpPeriodicity; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }
}


// create look-up table used to speed up index search
//

void Spin1_2ChainWithTranslations::CreateLookUpTable()
{
  int TmpHilbertSpaceDimension = this->HilbertSpaceDimension;
  // create the look-up table
  unsigned long Max = ((unsigned long) 1) << ((this->ChainLength << 1) - this->LookUpTableShift + 1);
  this->LookUpTable = new long [Max + 1];
  long LowPos;
  long MidPos;
  long HighPos;
  unsigned long Max2 = (this->StateDescription[TmpHilbertSpaceDimension - 1]) >> this->LookUpTableShift;
  for (unsigned long i = 0; i <= Max2; ++i)
    {
      LowPos = 0;
      HighPos = TmpHilbertSpaceDimension - 1;
      while ((HighPos - LowPos) > 1)
	{
	  MidPos = (HighPos + LowPos) >> 1;
	  if (this->StateDescription[MidPos] >= (i << this->LookUpTableShift))
	    HighPos = MidPos;
	  else
	    LowPos = MidPos;
	}      
      this->LookUpTable[i] = LowPos;
    }
  --TmpHilbertSpaceDimension;
  for (unsigned long i = Max2 + 1; i <= Max; ++i)    
    this->LookUpTable[i] = TmpHilbertSpaceDimension;
  ++TmpHilbertSpaceDimension;
}

// evaluate Hilbert space dimension
//
// nbrSpins = number of spins
// sz = twice the z projection of the total momentum
// return value = Hilbert space dimension

long Spin1_2ChainWithTranslations::EvaluateHilbertSpaceDimension(int nbrSpins, int szMax)
{
   FactorialCoefficient Coef;
   Coef.SetToOne();
   Coef.FactorialMultiply(nbrSpins);
   Coef.FactorialDivide((nbrSpins + szMax) / 2);
   Coef.FactorialDivide((nbrSpins - szMax) / 2);
   return Coef.GetIntegerValue();
}

// evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
// 
// nbrSites = number of sites that are part of the A subsytem 
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix Spin1_2ChainWithTranslations::EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture)
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
  Spin1_2Chain TmpDestinationHilbertSpace(nbrSites, szSector, 1000000);
  Spin1_2Chain TmpHilbertSpace(this->ChainLength - nbrSites, this->Sz - szSector, 1000000);

  ComplexMatrix TmpEntanglementMatrix(TmpHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  
  int Shift = nbrSites;
  int MinIndex = 0;
  int MaxIndex = TmpHilbertSpace.HilbertSpaceDimension;
  int TmpNbrTranslation;
  int TmpNbrTranslationToIdentity;
  Complex* TmpPhases = new Complex [this->ChainLength];
  double Coef = 2.0 * M_PI * ((double) this->Momentum) / ((double) this->ChainLength);
  for (int i = 0; i < this->ChainLength; ++i)
    {
      TmpPhases[i] = Phase(Coef * ((double) i));
    }

  unsigned long Mask1 = (0x1ul << Shift) - 0x1ul;
  unsigned long Mask2 = (0x1ul << this->ChainLength) - 0x1ul;
  for (; MinIndex < MaxIndex; ++MinIndex)    
    {
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex] << Shift;
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
  {
    unsigned long TmpState2 = (TmpState | (TmpDestinationHilbertSpace.StateDescription[j] & Mask1)) & Mask2;
    double Coefficient = 1.0;
    int TmpPos = this->SymmetrizeResult(TmpState2, 1, Coefficient, TmpNbrTranslation);
    if (TmpPos != this->HilbertSpaceDimension)
      {
        TmpEntanglementMatrix.AddToMatrixElement(MinIndex, j, groundState[TmpPos] * TmpPhases[TmpNbrTranslation] / sqrt((double) this->NbrStateInOrbit[TmpPos]));
      }
  }
    }
  delete[] TmpPhases;
  return TmpEntanglementMatrix;
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition.
// 
// nbrSpinUp = number of spin up that belong to the subsytem 
// kSector = momentum of the subsystem
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix Spin1_2ChainWithTranslations::EvaluatePartialDensityMatrixParticlePartition (int nbrSpinUpSector, int kSector, RealVector& groundState)
{
  if (nbrSpinUpSector == 0)
    {
      if (kSector == 0)
	{
	  HermitianMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }

  int TmpTotalNbrSpinUp = (this->ChainLength + this->Sz) >> 1;

  if (nbrSpinUpSector == TmpTotalNbrSpinUp)
    {
      if (kSector == this->Momentum)
	{
	  HermitianMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }

  int ComplementaryNbrSpinUpSector = TmpTotalNbrSpinUp - nbrSpinUpSector;
  int ComplementaryKSector = this->Momentum - kSector;
  if (ComplementaryKSector < 0)
    ComplementaryKSector += (this->ChainLength / this->StateShift);
  BinomialCoefficients TmpBinomial (TmpTotalNbrSpinUp);
  double TmpInvBinomial = 1.0 / (TmpBinomial(TmpTotalNbrSpinUp, nbrSpinUpSector));

  Spin1_2ChainWithTranslations TmpDestinationHilbertSpace(this->ChainLength, kSector, this->StateShift, nbrSpinUpSector, 1 << 18, 1 << 18);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpNbrTranslations = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;
  Spin1_2ChainWithTranslations TmpHilbertSpace(this->ChainLength, ComplementaryKSector, this->StateShift, ComplementaryNbrSpinUpSector, 1 << 18, 1 << 18);
  TmpInvBinomial = sqrt(TmpInvBinomial);
  int TmpTotalTranslations = this->ChainLength / this->StateShift;
  Complex* TmpExponentialArray = new Complex[TmpTotalTranslations];
  for (int i = 0; i < TmpTotalTranslations; ++i)
    TmpExponentialArray[i] = Polar(1.0, 2.0 * M_PI * ((double) (i * this->Momentum)) / ((double) TmpTotalTranslations));

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpState = TmpHilbertSpace.StateDescription[MinIndex];
      for (int k= 0; k < TmpTotalTranslations; ++k)
	{
	  TmpState = (TmpState >> this->StateShift) | ((TmpState & this->StateMask) << this->ComplementaryStateShift);
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState2 = TmpDestinationHilbertSpace.StateDescription[j];
	      if ((TmpState & TmpState2) == 0x0ul)
		{
		  int TmpPos = this->FindStateIndex(TmpState | TmpState2);
		  if (TmpPos != this->HilbertSpaceDimension)
		    {
		      TmpStatePosition[Pos] = TmpPos;
		      TmpStatePosition2[Pos] = j;
		      TmpStateCoefficient[Pos] = TmpInvBinomial;
		      TmpNbrTranslations[Pos] = k;
		      ++Pos;
		    }
		}
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      Complex TmpValue = TmpExponentialArray[TmpTotalTranslations - TmpNbrTranslations[j]];
	      TmpValue *= groundState[TmpStatePosition[j]] * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], 
							TmpExponentialArray[TmpNbrTranslations[k]] * TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		  }
	    }
	}
    }
  delete[] TmpExponentialArray;
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  delete[] TmpStateCoefficient;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }  
}


// return the Bosonic Occupation of a given state in the basis
//
// index = index of the state in the basis
// finalState = reference on the array where the monomial representation has to be stored

void Spin1_2ChainWithTranslations::GetBosonicOccupation (unsigned int index, int * finalState)
{
  for (int i = 0; i < this->ChainLength; i++)
    {
      finalState[i] = (this->StateDescription[index] >> ((unsigned long) i) )& 0x1ul;
    }
}

// convert the state on the site to its binary representation
//
// state = state to be stored
// sitePosition = position on the chain of the state
// return integer that code the state

unsigned long Spin1_2ChainWithTranslations::EncodeSiteState(int physicalState, int sitePosition)
{
  return  physicalState << sitePosition;
}

