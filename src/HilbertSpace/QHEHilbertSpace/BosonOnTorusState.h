////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of state of boson on a torus                   //
//                                                                            //
//                        last modification : 14/10/2003                      //
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


#ifndef BOSONONTORUSSTATE_H
#define BOSONONTORUSSTATE_H


#include "config.h"

#include <iostream>


using std::ostream;


// get the reduced number of state (aka the number of unsigned long per state)
// 
// nbrState = number of states in which bosons can lie
// return value = reduced number of state
int GetReducedNbrState (int nbrState);


class BosonOnTorusState
{

 private:

  // array describing the boson state
  unsigned long* StateDescription;

 public:

  // default constructor
  // 
  BosonOnTorusState();
  
  // basic constructor
  // 
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state)
  BosonOnTorusState(int reducedNbrState);
  
  // copy constructor
  // 
  // state = reference on the state to copy
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state)
  BosonOnTorusState(BosonOnTorusState& state, const int& reducedNbrState);

  // destructor
  // 
  ~BosonOnTorusState();

  // assign a state to the current one
  //
  // state = reference on the state to assign
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state)
  // return value = reference on  the current state
  BosonOnTorusState& Assign(BosonOnTorusState& state, int& reducedNbrState);

  // get hash key associated to the state
  //
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state)
  // keyMask = mask giving how many bits have to be kept
  // return value = hash key
  unsigned long GetHashKey (int& reducedNbrState, unsigned long& keyMask);

  // get the highest state index for which the state is occupied
  //
  // reducedNbrState = reduced number of state (aka the number of unsigned long per state)
  // return value = highest state index
  int GetHighestIndex (int& reducedNbrState);

  // set occupation of a state 
  //
  // stateIndex = index of the state whose occupation has to be set 
  // nbrBosons = number of bosons in the given state
  void SetOccupation (const int& stateIndex, const int& nbrBosons);

  // get occupation of a state 
  //
  // stateIndex = index of the state whose occupation has to be set 
  // return value = number of bosons in the given state
  int GetOccupation (const int& stateIndex);

  // increment occupation of a state 
  //
  // stateIndex = index of the state whose occupation has to be set 
  // tmpShift1 = refence on temporary variable used to evaluated a bit shift
  // tmpShift2 = refence on temporary variable used to evaluated a bit shift
  void IncrementOccupation (const int& stateIndex, int& tmpShift1, int& tmpShift2);
  
  // decrement occupation of a state (without testing if the state es empty)
  //
  // stateIndex = index of the state whose occupation has to be set 
  // tmpShift1 = refence on temporary variable used to evaluated a bit shift
  // tmpShift2 = refence on temporary variable used to evaluated a bit shift
  void DecrementOccupation (const int& stateIndex, int& tmpShift1, int& tmpShift2);
  
  // test if the state is empty an if it is not, decrement its occupation
  //
  // stateIndex = index of the state whose occupation has to be set 
  // tmpShift1 = refence on temporary variable used to evaluated a bit shift
  // tmpShift2 = refence on temporary variable used to evaluated a bit shift
  // return value = false if the state is empty
  inline bool TestAndDecrementOccupation (const int& stateIndex, int& tmpShift1, int& tmpShift2);
  
  // swap two states
  //
  // state = reference on the state to swap with the current one
  void SwapStates (BosonOnTorusState& state);

  // test if the current state is identical to another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the two states are identical
  bool Equal (BosonOnTorusState& state, int& reducedNbrState);

  // test if the current state is different to another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the two states are different
  bool Different (BosonOnTorusState& state, int& reducedNbrState);
  
  // test if the current state is greater than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater than the other state
  bool Greater (BosonOnTorusState& state, int& reducedNbrState);

  // test if the current state is greater or equal than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater or equal than the other state
  bool GreaterOrEqual (BosonOnTorusState& state, int& reducedNbrState);

  // test if the current state is greater than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater than the other state
  bool Lower (BosonOnTorusState& state, int& reducedNbrState);

  // test if the current state is greater or equal than another state
  //
  // state = reference on the state to compare with
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // return value = true if the current state is greater or equal than the other state
  bool LowerOrEqual (BosonOnTorusState& state, int& reducedNbrState);

  // put the state in a canonical form
  // 
  // 
  //  void PutInCanonicalForm(int& nbrTranslation);

  // shift a state to the left
  //
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
  // nbrTranslation = magnitude of the translation to apply
  void LeftShiftState(int& reducedNbrState, int& nbrStateRemainder, int nbrTranslation);

  // print a given state
  //
  // str = reference on current output stream 
  // reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
  // nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
  // return value = reference on current output stream 
  ostream& PrintState (ostream& str, int& reducedNbrState, int& nbrStateRemainder);

};

// get the reduced number of state (aka the number of unsigned long per state)
// 
// nbrState = number of states in which bosons can lie
// return value = reduced number of state
inline int GetReducedNbrState (int nbrState)
{
#ifdef __64_BITS__
  if (nbrState & ((unsigned long) 0x7))
    return (nbrState >> 3);
  else
    return ((nbrState >> 3) + 1);
#else
  if (nbrState & ((unsigned long) 0x3))
    return (nbrState >> 2);
  else
    return ((nbrState >> 2) + 1);
#endif
}

// assign a state to the current one
//
// state = reference on the state to assign
// reducedNbrState = reduced number of state (aka the number of unsigned long per state)
// return value = reference on  the current state

inline BosonOnTorusState& BosonOnTorusState::Assign(BosonOnTorusState& state, int& reducedNbrState)
{
  for (int i = 0; i < reducedNbrState; ++i)
    this->StateDescription[i] = state.StateDescription[i];
  return *this;
}

// get hash key associated to the state
//
// reducedNbrState = reduced number of state (aka the number of unsigned long per state)
// keyMask = mask giving how many bits have to be kept
// return value = hash key

inline unsigned long BosonOnTorusState::GetHashKey (int& reducedNbrState, unsigned long& keyMask)
{
  unsigned long Key = this->StateDescription[0];
  for (int i = 1; i < reducedNbrState; ++i)
    Key += this->StateDescription[i];
  return (Key & keyMask);
}

// get the highest state index for which the state is occupied
//
// reducedNbrState = reduced number of state (aka the number of unsigned long per state)
// return value = highest state index

inline int BosonOnTorusState::GetHighestIndex (int& reducedNbrState)
{
  int Index = 0;
  for (int i = reducedNbrState - 1; ((i >= 0) && (Index == 0)); --i)
    {
#ifdef __64_BITS__
      if ((this->StateDescription[i] & ((unsigned long) 0xff00000000000000)) != 0)
	Index = (i << 3) + 7;
      else
	if ((this->StateDescription[i] & ((unsigned long) 0xff000000000000)) != 0)
	  Index = (i << 3) + 6;
	else
	  if ((this->StateDescription[i] & ((unsigned long) 0xff0000000000)) != 0)
	    Index = (i << 3) + 5;
	  else
	    if ((this->StateDescription[i] & ((unsigned long) 0xff00000000)) != 0)
	      Index = (i << 3) + 4;
	    else
	      if ((this->StateDescription[i] & ((unsigned long) 0xff000000)) != 0)
		Index = (i << 3) + 3;
	      else
		if ((this->StateDescription[i] & ((unsigned long) 0xff0000)) != 0)
		  Index = (i << 3) + 2;
		else
		  if ((this->StateDescription[i] & ((unsigned long) 0xff00)) != 0)
		    Index = (i << 3) + 1;
		  else
		    if ((this->StateDescription[i] & ((unsigned long) 0xff)) != 0)
		      Index = (i << 3);
#else
      if ((this->StateDescription[i] & ((unsigned long) 0xff000000)) != 0)
	Index = (i << 2) + 3;
      else
	if ((this->StateDescription[i] & ((unsigned long) 0xff0000)) != 0)
	  Index = (i << 2) + 2;
	else
	  if ((this->StateDescription[i] & ((unsigned long) 0xff00)) != 0)
	    Index = (i << 2) + 1;
	  else
	    if ((this->StateDescription[i] & ((unsigned long) 0xff)) != 0)
	      Index = (i << 2);
#endif
    }
  return Index;  
}

// set occupation of a state 
//
// stateIndex = index of the state whose occupation has to be set 
// nbrBosons = number of bosons in the given state

inline void BosonOnTorusState::SetOccupation (const int& stateIndex, const int& nbrBosons)
{
#ifdef __64_BITS__
  this->StateDescription[stateIndex >> 3] &= ~(((unsigned long) 0xff) << ((stateIndex & ((unsigned long) 0x7)) << 3));
  this->StateDescription[stateIndex >> 3] |= nbrBosons << ((stateIndex & ((unsigned long) 0x7)) << 3);
#else
  this->StateDescription[stateIndex >> 2] &= ~(((unsigned long) 0xff) << ((stateIndex & ((unsigned long) 0x3)) << 3));
  this->StateDescription[stateIndex >> 2] |= nbrBosons << ((stateIndex & ((unsigned long) 0x3)) << 3);
#endif
}

// get occupation of a state 
//
// stateIndex = index of the state whose occupation has to be set 
// return value = number of bosons in the given state

inline int BosonOnTorusState::GetOccupation (const int& stateIndex)
{
#ifdef __64_BITS__
  return ((this->StateDescription[stateIndex >> 3] >> ((stateIndex & ((unsigned long) 0x7)) << 3)) & ((unsigned long) 0xff));
#else
  return ((this->StateDescription[stateIndex >> 2] >> ((stateIndex & ((unsigned long) 0x3)) << 3)) & ((unsigned long) 0xff));
#endif
}

// increment occupation of a state 
//
// stateIndex = index of the state whose occupation has to be set 
// tmpShift1 = refence on temporary variable used to evaluated a bit shift
// tmpShift2 = refence on temporary variable used to evaluated a bit shift

inline void BosonOnTorusState::IncrementOccupation (const int& stateIndex, int& tmpShift1, int& tmpShift2)
{
#ifdef __64_BITS__
  tmpShift1 = stateIndex >> 3;
  tmpShift2 = (stateIndex & ((unsigned long) 0x7)) << 3;
#else
  tmpShift1 = stateIndex >> 2;
  tmpShift2 = (stateIndex & ((unsigned long) 0x3)) << 3;
#endif
  this->StateDescription[tmpShift1] &= ~(((unsigned long) 0xff) << tmpShift2);
  this->StateDescription[tmpShift1] |= (((this->StateDescription[tmpShift1] >> tmpShift2) & ((unsigned long) 0xff)) + 1) << tmpShift2;
}

// decrement occupation of a state (without testing if the state es empty)
//
// stateIndex = index of the state whose occupation has to be set 
// tmpShift1 = refence on temporary variable used to evaluated a bit shift
// tmpShift2 = refence on temporary variable used to evaluated a bit shift

inline void BosonOnTorusState::DecrementOccupation (const int& stateIndex, int& tmpShift1, int& tmpShift2)
{
#ifdef __64_BITS__
  tmpShift1 = stateIndex >> 3;
  tmpShift2 = (stateIndex & 0x7) << 3;
#else
  tmpShift1 = stateIndex >> 2;
  tmpShift2 = (stateIndex & ((unsigned long) 0x3)) << 3;
#endif
  this->StateDescription[tmpShift1] &= ~(((unsigned long) 0xff) << tmpShift2);
  this->StateDescription[tmpShift1] |= (((this->StateDescription[tmpShift1] >> tmpShift2) & ((unsigned long) 0xff)) - 1) << tmpShift2;
}

// test if the state is empty an if it is not, decrement its occupation
//
// stateIndex = index of the state whose occupation has to be set 
// tmpShift1 = refence on temporary variable used to evaluated a bit shift
// tmpShift2 = refence on temporary variable used to evaluated a bit shift
// return value = false if the state is empty

inline bool BosonOnTorusState::TestAndDecrementOccupation (const int& stateIndex, int& tmpShift1, int& tmpShift2)
{
#ifdef __64_BITS__
  tmpShift1 = stateIndex >> 3;
  tmpShift2 = (stateIndex & ((unsigned long) 0x7)) << 3;
#else
  tmpShift1 = stateIndex >> 2;
  tmpShift2 = (stateIndex & ((unsigned long) 0x3)) << 3;
#endif
  if (this->StateDescription[tmpShift1] >> tmpShift2)
    {
      this->StateDescription[tmpShift1] &= ~(((unsigned long) 0xff) << tmpShift2);
      this->StateDescription[tmpShift1] |= (((this->StateDescription[tmpShift1] >> tmpShift2) & ((unsigned long) 0xff)) - 1) << tmpShift2;
      return true;
    }
  else
    return false;
}

// swap two states
//
// state = reference on the state to swap with the current one

inline void BosonOnTorusState::SwapStates (BosonOnTorusState& state)
{
  unsigned long* TmpStateDescription = this->StateDescription;
  this->StateDescription = state.StateDescription;
  state.StateDescription = TmpStateDescription;
}

// test if the current state is identical to another state
//
// state = reference on the state to compare with
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the two states are identical

inline bool BosonOnTorusState::Equal (BosonOnTorusState& state, int& reducedNbrState)
{
  while (reducedNbrState >= 0)
    if (state.StateDescription[reducedNbrState] != this->StateDescription[reducedNbrState])
      return false;
    else
      --reducedNbrState;
  return true;
}

// test if the current state is different to another state
//
// state = reference on the state to compare with
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the two states are different

inline bool BosonOnTorusState::Different (BosonOnTorusState& state, int& reducedNbrState)
{
  while (reducedNbrState >= 0)
    if (state.StateDescription[reducedNbrState] != this->StateDescription[reducedNbrState])
      return true;
    else
      --reducedNbrState;
  return false;
}

// test if the current state is greater than another state
//
// state = reference on the state to compare with
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the current state is greater than the other state

inline bool BosonOnTorusState::Greater (BosonOnTorusState& state, int& reducedNbrState)
{
  while (reducedNbrState >= 0)
    if (state.StateDescription[reducedNbrState] >= this->StateDescription[reducedNbrState])
      return false;
    else
      --reducedNbrState;
  return true;
}


// test if the current state is greater or equal than another state
//
// state = reference on the state to compare with
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the current state is greater or equal than the other state

inline bool BosonOnTorusState::GreaterOrEqual (BosonOnTorusState& state, int& reducedNbrState)
{
  while (reducedNbrState >= 0)
    if (state.StateDescription[reducedNbrState] > this->StateDescription[reducedNbrState])
      return false;
    else
      --reducedNbrState;
  return true;
}

// test if the current state is lower than another state
//
// state = reference on the state to compare with
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the current state is lower than the other state

inline bool BosonOnTorusState::Lower (BosonOnTorusState& state, int& reducedNbrState)
{
  while (reducedNbrState >= 0)
    if (state.StateDescription[reducedNbrState] <= this->StateDescription[reducedNbrState])
      return false;
    else
      --reducedNbrState;
  return true;
}


// test if the current state is lower or equal than another state
//
// state = reference on the state to compare with
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// return value = true if the current state is lower or equal than the other state

inline bool BosonOnTorusState::LowerOrEqual (BosonOnTorusState& state, int& reducedNbrState)
{
  while (reducedNbrState >= 0)
    if (state.StateDescription[reducedNbrState] < this->StateDescription[reducedNbrState])
      return false;
    else
      --reducedNbrState;
  return true;
}

// put the state in a canonical form
// 
// 

//inline void BosonOnTorusState::PutInCanonicalForm(BosonOnTorusState& tmpState, int& nbrTranslation)
//{
  //  if ()
//}

// shift a state to the left
//
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
// nbrTranslation = magnitude of the translation to apply

inline void BosonOnTorusState::LeftShiftState(int& reducedNbrState, int& nbrStateRemainder, int nbrTranslation)
{
#ifdef __64_BITS__
  unsigned long Remainder;
  // WARNING : this part of the code is not speed optimal but it doesn't need too much temporary variables and is fast enough
  // if the translation is lower than 15 (true for most of the system size reachable with a 64 bits architecture)

#else
  unsigned long Remainder;
  // WARNING : this part of the code is not speed optimal but it doesn't need too much temporary variables and is fast enough
  // if the translation is lower than 8 (true for most of the system size reachable with a 32 bits architecture)
  for (int j = 0; j < (nbrTranslation >> 2); ++j)
    {
      Remainder = this->StateDescription[0];
      for (int i = 0; i < reducedNbrState; ++i)
	{
	  this->StateDescription[i] = this->StateDescription[i + 1];
	}      
      this->StateDescription[reducedNbrState - 1] |= Remainder << (nbrStateRemainder << 3);
      this->StateDescription[reducedNbrState] = Remainder >> ((4 - nbrStateRemainder) << 3);
    }
  switch (nbrTranslation & 0x3)
    {
    case 1:
      {	
	Remainder = this->StateDescription[0] & ((unsigned long) 0xff);
	for (int i = 0; i < reducedNbrState; ++i)
	  {
	    this->StateDescription[i] >>= 8;
	    this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xff)) << 24;
	  }
	this->StateDescription[reducedNbrState] >>= 8;
	this->StateDescription[reducedNbrState] |= Remainder << ((nbrStateRemainder - 1) << 3);
      }
      break;
    case 2:
      {	
	if (nbrStateRemainder >= 2)
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffff);
	    for (int i = 0; i < reducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 16;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffff)) << 16;
	      }
	    this->StateDescription[reducedNbrState] >>= 16;
	    this->StateDescription[reducedNbrState] |= Remainder << ((nbrStateRemainder - 2) << 3);
	  }
	else
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffff);
	    int ReducedNbrState = reducedNbrState - 1;
	    for (int i = 0; i < ReducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 16;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffff)) << 16;
	      }
	    this->StateDescription[ReducedNbrState] >>= 16;
	    this->StateDescription[ReducedNbrState] |= this->StateDescription[reducedNbrState] << 16;
	    this->StateDescription[ReducedNbrState] |= (Remainder & ((unsigned long) 0xff)) << 24;
	    this->StateDescription[reducedNbrState] = Remainder >> 8;
	  }
      }
      break;
    case 3:
      {	
	if (nbrStateRemainder >= 3)
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffff);
	    for (int i = 0; i < reducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 24;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffff)) << 8;
	      }
	    this->StateDescription[reducedNbrState] >>= 24;
	    this->StateDescription[reducedNbrState] |= Remainder << ((nbrStateRemainder - 3) << 3);
	  }
	else
	  {
	    Remainder = this->StateDescription[0] & ((unsigned long) 0xffffff);
	    int ReducedNbrState = reducedNbrState - 1;
	    for (int i = 0; i < ReducedNbrState; ++i)
	      {
		this->StateDescription[i] >>= 24;
		this->StateDescription[i] |= (this->StateDescription[i + 1]  & ((unsigned long) 0xffffff)) << 8;
	      }
	    this->StateDescription[ReducedNbrState] >>= 24;
	    this->StateDescription[ReducedNbrState] |= this->StateDescription[reducedNbrState] << 8;
	    this->StateDescription[ReducedNbrState] |= (Remainder << ((nbrStateRemainder + 1) << 3));
	    this->StateDescription[reducedNbrState] = Remainder >> ((3 - nbrStateRemainder) << 3);
	  }
      }
      break;
    }
#endif
}

// print a given state
//
// str = reference on current output stream 
// reducedNbrState = reference on the reduced number of state (aka the number of unsigned long per state) minus 1
// nbrStateRemainder = number of the state in the last unsigned long array describing the whole state
// return value = reference on current output stream 

inline ostream& BosonOnTorusState::PrintState (ostream& str, int& reducedNbrState, int& nbrStateRemainder)
{
  for (int i = 0; i < reducedNbrState; ++i)
    {
#ifdef __64_BITS__
      for (int j = 0; j < 8; ++j)
#else
      for (int j = 0; j < 4; ++j)
#endif
	str << ((this->StateDescription[i] >> (j << 3)) & ((unsigned long) 0xff)) << " ";
    }
  for (int i = 0; i < nbrStateRemainder; ++i)
    str << ((this->StateDescription[reducedNbrState] >> (i << 3)) & ((unsigned long) 0xff)) << " ";
  return str;
}


#endif
