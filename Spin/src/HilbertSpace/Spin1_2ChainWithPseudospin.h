////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of spin 1/2 chain with a fixed Sz value                //
//                                                                            //
//                        last modification : 29/06/2015                      //
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


#ifndef SPIN1_2CHAINWITHPSEUDOSPIN_H
#define SPIN1_2CHAINWITHPSEUDOSPIN_H


#include "config.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::ostream;


class Spin1_2ChainWithPseudospin : public Spin1_2ChainNew
{

 protected:

 public:


  // default constructor
  //
  Spin1_2ChainWithPseudospin ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin 1/2
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Spin1_2ChainWithPseudospin (int chainLength, int sz, int memorySize);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1_2ChainWithPseudospin (const Spin1_2ChainWithPseudospin& chain);

  // destructor
  //
  ~Spin1_2ChainWithPseudospin ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1_2ChainWithPseudospin& operator = (const Spin1_2ChainWithPseudospin& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // return value of spin projection on (Oz) for a given state
  //
  // index = index of the state to test
  // return value = spin projection on (Oz)
  virtual int TotalSz (int index);

  // compute the parity (prod_i Sz_i) for a given state
  //
  // state = index of the state to be applied on Sz_i operator
  // return value = total Sz value
  virtual unsigned long GetParity (int state);
  
  // translate a state assuming the system have periodic boundary
  // conditions (increasing the site index)
  //
  // nbrTranslations = number of translations to apply
  // state = index of the state to translate 
  // return value = index of resulting state
  virtual int TranslateState (int nbrTranslations, int state);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

 protected:

  // find state index
  //
  // stateDescription = state description
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription);

  // evaluate Hilbert space dimension
  //
  // sz = twice the Sz value
  // nbrSites = number of sites
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int sz, int nbrSites);

  // generate all states
  // 
  // nbrSpinUp = number of spin up
  // currentPosition = current position to consider in the chain
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  long GenerateStates(int nbrSpinUp, int currentPosition, long pos);

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  // stateMask = an optional mask to apply to each state to focus on the relevant bits
#ifdef __64_BITS__
  virtual void GenerateLookUpTable(unsigned long memory, unsigned long stateMask = 0xfffffffffffffffful);
#else
  virtual void GenerateLookUpTable(unsigned long memory, unsigned long stateMask = 0xfffffffful);
#endif

  // return the Bosonic Occupation of a given state in the basis
  //
  // index = index of the state in the basis
  // finalState = reference on the array where the monomial representation has to be stored
  virtual void GetBosonicOccupation (unsigned int index, int * finalState);
  
};

#endif


