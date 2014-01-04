////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of Potts-3 chain with the translation symmetry            //
//                                                                            //
//                        last modification : 01/01/2014                      //
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


#ifndef POTTS3CHAINWITHTRANSLATIONS_H
#define POTTS3CHAINWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"

#include <iostream>


using std::ostream;


class HermitianMatrix;
class RealMatrix;
class Matrix;
class SubspaceSpaceConverter;
class AbstractQuantumNumber;


class Potts3ChainWithTranslations : public AbstractSpinChainWithTranslations
{

 protected:

  // number of sites on th chain
  int ChainLength;

  // momentum 
  int Momentum;

  // shift to apply to move the spin from one end to the other one
  int ComplementaryStateShift;

  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;
  // number of state in each orbit
  int* NbrStateInOrbit;

  // total quantum number for the chain (i.e. either 0, 1 or 2)
  int Sz;
  // flag to indicate if Sz is fixed 
  bool FixedQuantumNumberFlag;
  
  int* LookUpTable;
  unsigned long LookUpTableMask;
  int LookUpTableSize;
  int LookUpPosition;

 public:

  unsigned long* ChainDescription;

  // default constructor
  //
  Potts3ChainWithTranslations ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // momemtum = momentum sector
  // memorySize = memory size in bytes allowed for look-up table
  Potts3ChainWithTranslations (int chainLength, int momentum, int memorySize);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // sz = twice the value of total Sz component
  // momemtum = momentum sector
  // memorySize = memory size in bytes allowed for look-up table
  Potts3ChainWithTranslations (int chainLength, int sz, int momentum, int memorySize);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Potts3ChainWithTranslations (const Potts3ChainWithTranslations& chain);

  // destructor
  //
  ~Potts3ChainWithTranslations ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Potts3ChainWithTranslations& operator = (const Potts3ChainWithTranslations& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  AbstractQuantumNumber* GetQuantumNumber (int index);

  // get the momentum of each state in the current Hilbert space
  //
  // return value = momentum value
  int GetMomentum();

  // return value of spin projection on (Oz) for a given state
  //
  // index = index of the state to test
  // return value = spin projection on (Oz)
  int TotalSz (int index);

  // return value of the value of the sum of the square of spin projection on (Oz) 
  //
  // index = index of the state to test
  // return value = twice spin projection on (Oz)
  double TotalSzSz (int index);

  // return index of resulting state from application of S+_i operator on a given state (only valid if there is no constraint on total Sz)
  //
  // i = operator position
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int Spi (int i, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i operator on a given state (only valid if there is no constraint on total Sz)
  //
  // i = operator position
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int Smi (int i, int state, double& coefficient, int& nbrTranslation);
    
  // return index of resulting state from application of S+_i S+_j operator on a given state
  //
  // i = position of first S+ operator
  // j = position of second S+ operator
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i S-_j operator on a given state
  //
  // i = position of first S- operator
  // j = position of second S- operator
  // state = index of the state to be applied on S-_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S+_i S+_i operator on a given state
  //
  // i = position of first S+ operator
  // state = index of the state to be applied on S+_i S+_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SpiSpi (int i, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i S-_i operator on a given state
  //
  // i = position of the S- operator
  // state = index of the state to be applied on S-_i S-_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SmiSmi (int i, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S+_i Sz_j operator on a given state
  //
  // i = position of S+ operator
  // j = position of Sz operator
  // state = index of the state to be applied on S+_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i Sz_j operator on a given state
  //
  // i = position of S- operator
  // j = position of Sz operator
  // state = index of the state to be applied on S-_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of Sz_i operator on a given state
  //
  // i = position of Sz operator
  // state = index of the state to be applied on Sz_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // return value = index of resulting state
  int Szi (int i, int state, double& coefficient);

  // return eigenvalue of Sz_i Sz_j associated to a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = corresponding eigenvalue
  double SziSzj (int i, int j, int state);

  // return value the product of the tau_j operators (Q operator) on a given state
  //
  // index = index of the state 
  // return value = Q value (0 for 1, 1 for exp(i 2 \pi / 3), -1 1 for exp(i 2 \pi / 3)) 
  double QValue (int index);

  // translate a state assuming the system have periodic boundary conditions (increasing the site index)
  //
  // nbrTranslations = number of translations to apply
  // state = index of the state to translate 
  // return value = index of resulting state
  int TranslateState (int nbrTranslations, int state);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter);

  // find state index
  //
  // state = state description
  // return value = corresponding index
  int FindStateIndex(unsigned long state);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

 protected:
  
  // find the index of state and the number of translations required to get the canonical form
  //
  // inputState = index of the state that was used as the initial state
  // outputState = description of the state once the operator has been applied
  // nbrTranslation = number of translations required to get the canonical form
  // nbrTranslationToIdentity = number of translations required to map the state on itself
  // coefficient = reference on the numerical coefficient
  // return value = index of the output state canonical form
  int FindStateIndexAndTransaltion(int& inputState, unsigned long& outputState, int& nbrTranslation, 
				   int& nbrTranslationToIdentity, double& coefficient);

  // find the canonical form of a state
  //
  // state = state description
  // nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
  // return value = canonical form of the state
  unsigned long FindCanonicalForm(unsigned long state, int& nbrTranslation);

  // find the canonical form of a state and find how many translations are needed to obtain the same state
  //
  // stateDescription = state description
  // nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
  // nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state
  // return value = canonical form of the state
  unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslation, int& nbrTranslationToIdentity);

  // find how many translations are needed to obtain the same state
  //
  // stateDescription = unsigned integer describing the state
  // return value = number of translation needed to obtain the same state
  int FindNumberTranslation(unsigned long stateDescription);

  // evaluate Hilbert space dimension
  //
  // currentSite = current site to occupy
  // currentSzValue = state current Sz value 
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int currentSite, int currentSzValue);

  // generate all states with no constraint on total Sz
  //
  // currentSite = current site to occupy
  // currentPosition = current position of the state that has to be considered
  // return value = number of generated states
  long GenerateStates(int currentSite, long currentPosition);

  // generate all states corresponding to a given total Sz
  //
  // currentSite = current site to occupy
  // currentSzValue = state current Sz value 
  // currentPosition = current position of the state that has to be considered
  // return value = number of generated states
  long GenerateStates(int currentSite, int currentSzValue, long currentPosition);

  // generate all states corresponding to a given momnetum
  //
  void GenerateStates();

};

// return index of resulting state from application of Sz_i operator on a given state
//
// i = position of Sz operator
// state = index of the state to be applied on Sz_i operator
// coefficient = reference on double where numerical coefficient has to be stored (0 for 1.0, 1.0 for exp(i 2 \pi / 3), 2.0 for exp(i 4 \pi / 3)) 
// return value = index of resulting state

inline int Potts3ChainWithTranslations::Szi (int i, int state, double& coefficient)
{
  coefficient = (((double) ((this->ChainDescription[state] >> (i << 1)) & 0x1ul))
		 - ((double) ((this->ChainDescription[state] >> ((i << 1) + 1)) & 0x1ul)));
  return state;
}

// return eigenvalue of Sz_i Sz_j associated to a given state
//
// i = first position
// j = second position
// state = index of the state to consider
// return value = corresponding eigenvalue

inline double Potts3ChainWithTranslations::SziSzj (int i, int j, int state)
{  
  return ((double) (((this->ChainDescription[state] >> (i << 1)) & 0x3ul)
		    * ((this->ChainDescription[state] >> (j << 1)) & 0x3ul)));
}

// return index of resulting state from application of S-_i S+_j operator on a given state
//
// i = position of S- operator
// j = position of S+ operator
// state = index of the state to be applied on S-_i S+_j operator
// coefficient = reference on double where numerical coefficient has to be stored
// nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of resulting state

inline int Potts3ChainWithTranslations::SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation)
{  
  unsigned long TmpState = this->ChainDescription[state];
  unsigned long TmpState2 = TmpState;
  j <<= 1;
  TmpState2 >>= j;
  TmpState2 &= 0x3ul;
  TmpState &= ~(0x3ul << j);
  switch (TmpState2)
    {
    case 0x1ul:
      TmpState |= 0x2ul << j;
      break;
    case 0x0ul:
      TmpState |= 0x1ul << j;
      break;
    }	  
  i <<= 1;
  TmpState2 = TmpState;
  TmpState2 >>= i;
  TmpState2 &= 0x3ul;
  TmpState &= ~(0x3ul << i);
  switch (TmpState2)
    {
    case 0x2ul:
      TmpState |= 0x1ul << i;      
      break;
    case 0x0ul:
      TmpState |= 0x2ul << i;      
      break;
    }	  
  return this->FindStateIndexAndTransaltion(state, TmpState, nbrTranslation, i, coefficient);
}

// find the index of state and the number of translations required to get the canonical form
//
// inputState = index of the state that was used as the initial state
// outputState = description of the state once the operator has been applied
// nbrTranslation = number of translations required to get the canonical form
// nbrTranslationToIdentity = number of translations required to map the state on itself
// coefficient = reference on the numerical coefficient
// return value = index of the output state canonical form

inline int Potts3ChainWithTranslations::FindStateIndexAndTransaltion(int& inputState, unsigned long& outputState, 
								     int& nbrTranslation, int& nbrTranslationToIdentity, double& coefficient)
{
  outputState = this->FindCanonicalForm(outputState, nbrTranslation, nbrTranslationToIdentity);
  if (((nbrTranslationToIdentity * this->Momentum) % this->ChainLength) != 0)
    {
      return this->HilbertSpaceDimension;
    }
  coefficient = this->RescalingFactors[this->NbrStateInOrbit[inputState]][nbrTranslationToIdentity];
  return this->FindStateIndex(outputState);
}

// return value the product of the tau_j operators (Q operator) on a given state
//
// index = index of the state 
// return value = Q value (0 for 1, 1 for exp(i 2 \pi / 3), -1 1 for exp(i 2 \pi / 3)) 

inline double Potts3ChainWithTranslations::QValue (int index)
{
  unsigned long Tmp = 0x0ul;
  unsigned long Tmp2 = this->ChainDescription[index];
  for (int i = 0; i < this->ChainLength; ++i)
    {
      Tmp += Tmp2 & 0x3ul;
      Tmp2 >>= 2;
    }
  return ((double) (Tmp % 3));
}

// find the canonical form of a state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// return value = canonical form of the state

inline unsigned long Potts3ChainWithTranslations::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState = stateDescription;
  int index = 1;  
  while (index < this->ChainLength)
    {
      stateDescription = (stateDescription >> 2) | ((stateDescription & ((unsigned long) 0x3)) << this->ComplementaryStateShift);
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = index;
	}
      ++index;
    }
  return CanonicalState;
}

// find the canonical form of a state and find how many translations are needed to obtain the same state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state
// return value = canonical form of the state

inline unsigned long Potts3ChainWithTranslations::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslation, int& nbrTranslationToIdentity)
{
  nbrTranslation = 0;
  nbrTranslationToIdentity = 1;
  unsigned long CanonicalState = stateDescription;
  unsigned long ReferenceState = stateDescription;
  stateDescription = (stateDescription >> 2) | ((stateDescription & ((unsigned long) 0x3)) << this->ComplementaryStateShift);
  while ((ReferenceState != stateDescription) && (nbrTranslationToIdentity < this->ChainLength))
    {
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = nbrTranslationToIdentity;
	}
      stateDescription = (stateDescription >> 2) | ((stateDescription & ((unsigned long) 0x3)) << this->ComplementaryStateShift);
      ++nbrTranslationToIdentity;
    }
  return CanonicalState;
}

// find how many translations are needed to obtain the same state
//
// stateDescription = unsigned integer describing the state
// return value = number of translation needed to obtain the same state

inline int Potts3ChainWithTranslations::FindNumberTranslation(unsigned long stateDescription)
{
  unsigned long TmpState = (stateDescription >> 2) | ((stateDescription & ((unsigned long) 0x3)) << this->ComplementaryStateShift);
  int index = 1;  
  while (TmpState != stateDescription)
    {
      TmpState = (TmpState >> 2) | ((TmpState & ((unsigned long) 0x3)) << this->ComplementaryStateShift);
      ++index;
    }
  return index;
}

// get the momentum of each state in the current Hilbert space
//
// return value = momentum value

inline int Potts3ChainWithTranslations::GetMomentum()
{
  return this->Momentum;
}


#endif


