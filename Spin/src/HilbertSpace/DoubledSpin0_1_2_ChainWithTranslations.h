////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of doubled spin 0 +1/2 chain with translations              //
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


#ifndef DOUBLEDSPIN0_1_2_CHAINWITHTRANSLATIONS_H
#define DOUBLEDSPIN0_1_2_CHAINWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractDoubledSpinChainWithTranslations.h"


#include <iostream>


using std::ostream;

class DoubledSpin0_1_2_chainWithTranslations : public AbstractDoubledSpinChainWithTranslations
{
  
 protected:
 public:

  // default constructor
  //
  DoubledSpin0_1_2_chainWithTranslations ();

  // constructor for complete Hilbert space with no restriction on momentum
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // memorySize = memory size in bytes allowed for look-up table
  // memorySlice = maximum amount of memory that can be allocated to partially evalauted the states
  DoubledSpin0_1_2_chainWithTranslations (int chainLength,  int diffSz, int memorySize, int memorySlice);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  DoubledSpin0_1_2_chainWithTranslations (int chainLength, int momentum, int sz, int memorySize, int memorySlice);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  DoubledSpin0_1_2_chainWithTranslations (const DoubledSpin0_1_2_chainWithTranslations & chain);

  // destructor
  //
  ~DoubledSpin0_1_2_chainWithTranslations ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  DoubledSpin0_1_2_chainWithTranslations& operator = (const DoubledSpin0_1_2_chainWithTranslations & chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  // return the Bosonic Occupation of a given state in the basis
  //
  // index = index of the state in the basis
  // return value bosonic occupation 
  inline void GetBosonicOccupation (unsigned int index, int * finalStateBra,int * finalStateKet);

  // convert the state on the site to its binary representation
  //
  // state = state to be stored
  // sitePosition = position on the chain of the state
  // return integer that code the state
  inline unsigned long EncodeSiteStateBra(int physicalState, int sitePosition);

  // convert the state on the site to its binary representation
  //
  // state = state to be stored
  // sitePosition = position on the chain of the state
  // return integer that code the state
  inline unsigned long EncodeSiteStateKet(int physicalState, int sitePosition);

 protected:

  // return value of twice spin projection on (Oz) for a given state
  //
  // stateDescription = state to which the spin projection has to be evaluated
  // return value = twice spin projection on (Oz)
  int GetTotalSz (unsigned long stateDescriptionBra, unsigned long stateDescriptionKet);
  
  // find the canonical form of a state
  //
  // state = state description
  // nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
  // return value = canonical form of the state
  inline void FindCanonicalForm ( unsigned long stateDescriptionBra, unsigned long stateDescriptionKet, unsigned long & canonicalStateBra, unsigned long & canonicalStateKet, int& nbrTranslation);

  // find the canonical form of a state and find how many translations are needed to obtain the same state
  //
  // stateDescription = state description
  // nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
  // nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state
  // return value = canonical form of the state
  inline void FindCanonicalForm ( unsigned long stateDescriptionBra, unsigned long stateDescriptionKet, unsigned long & canonicalStateBra, unsigned long & canonicalStateKet, int& nbrTranslation, int& nbrTranslationToIdentity);
  
  // find how many translations are needed to obtain the same state
  //
  // stateDescription = unsigned integer describing the state
  // return value = number of translation needed to obtain the same state
  int FindNumberTranslation(unsigned long stateDescriptionBra, unsigned long stateDescriptionKet);

  // generate all states corresponding to the constraints
  // 
  // lengthBra = length of the chain to be decided for bra spins
  // lengthBra = length of the chain to be decided for ket spins
  // diffSz = difference of spin projection between bra and ket chain
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int lengthBra, int lengthKet, int diffSz, long pos);
  
  long ShiftedEvaluateHilbertSpaceDimension(int lengthBra, int lengthKet, int diffSz);

};

// return the Bosonic Occupation of a given state in the basis
//
// index = index of the state in the basis
// finalState = reference on the array where the monomial representation has to be stored

inline void DoubledSpin0_1_2_chainWithTranslations::GetBosonicOccupation (unsigned int index, int * finalStateBra,int * finalStateKet)
{
  for (unsigned long i = 0; i < this->ChainLength; i++)
    {
      finalStateBra[i] = ((this->ChainDescriptionBra[index] >> ((unsigned long) 2*i) )& 0x3ul);
      finalStateKet[i] = ((this->ChainDescriptionKet[index] >> ((unsigned long) 2*i) )& 0x3ul);
    }
}

// convert the state on the site to its binary representation
//
// state = state to be stored
// sitePosition = position on the chain of the state
// return integer that code the state

inline unsigned long DoubledSpin0_1_2_chainWithTranslations::EncodeSiteStateBra(int physicalState, int sitePosition)
{
  return  physicalState << (2*sitePosition);
}

// convert the state on the site to its binary representation
//
// state = state to be stored
// sitePosition = position on the chain of the state
// return integer that code the state

inline unsigned long DoubledSpin0_1_2_chainWithTranslations::EncodeSiteStateKet(int physicalState, int sitePosition)
{
  return  physicalState << (2*sitePosition);
}


#endif


