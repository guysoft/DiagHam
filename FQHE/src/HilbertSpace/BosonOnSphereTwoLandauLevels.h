////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                     class of bosons on sphere including two                //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 08/09/2009                      //
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


#ifndef BOSONONSPHERETWOLANDAULEVELS_H
#define BOSONONSPHERETWOLANDAULEVELS_H


#include "config.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"

#include <iostream>



class BosonOnSphereTwoLandauLevels :  public BosonOnSphereWithSpin
{

 protected:

  // maximum Lz value reached by a boson with a spin up
  int LzMaxUp;
  // maximum Lz value reached by a boson with a spin down
  int LzMaxDown;
  // shift to apply on the spin up part
  int LzShiftUp;
  // shift to apply on the spin down part
  int LzShiftDown;
  // sum of LzShiftUp and LzShiftDown
  int LzTotalShift;

 public:

  // default constructor
  //
  BosonOnSphereTwoLandauLevels();

  // basic constructor with no contraint on the number of particles per spin component
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMaxUp = twice the maximum Lz value reached by a boson with a spin up
  // lzMaxDown = twice the maximum Lz value reached by a boson with a spin down
  // memory = amount of memory granted for precalculations
  BosonOnSphereTwoLandauLevels (int nbrBosons, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereTwoLandauLevels(const BosonOnSphereTwoLandauLevels& bosons);

  // destructor
  //
  ~BosonOnSphereTwoLandauLevels ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereTwoLandauLevels& operator = (const BosonOnSphereTwoLandauLevels& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // create an SU(2) state from two U(1) state
  //
  // upState = vector describing the up spin part of the output state
  // upStateSpace = reference on the Hilbert space associated to the up spin part
  // downState = vector describing the down spin part of the output state
  // downStateSpace = reference on the Hilbert space associated to the down spin part  
  // return value = resluting SU(2) state
  virtual RealVector ForgeSU2FromU1(RealVector& upState, BosonOnSphere& upStateSpace, RealVector& downState, BosonOnSphere& downStateSpace);

 protected:

  // convert a bosonic state into its fermionic counterpart
  //
  // initialState = reference on the array where initialbosonic  state is stored
  // initialStateLzMax = reference on the initial bosonic state maximum Lz value
  // return value = corresponding fermionic state
  unsigned long BosonToFermion(unsigned long*& initialState, int& initialStateLzMax);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialState = initial fermionic state
  // initialStateLzMax = initial fermionic state maximum Lz value
  // finalState = reference on the array where the bosonic state has to be stored
  // finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored
  void FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState, int& finalStateLzMax);

  // evaluate Hilbert space dimension without constraint on the number of particles per level
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // lzMaxUp = momentum maximum value for a boson
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int lzMax, int totalLz, long pos);

};

// convert a bosonic state into its fermionic counterpart
//
// initialStateUp = reference on the array where the initial up bosonic state is stored
// initialStateLzMaxUp = reference on the initial up bosonic state maximum Lz value
// initialStateDown = reference on the array where the initial down bosonic state is stored
// initialStateLzMaxDown = reference on the initial down bosonic state maximum Lz value
// return value = corresponding fermionic state

inline unsigned long BosonOnSphereTwoLandauLevels::BosonToFermion(unsigned long*& initialStateUp, int& initialStateLzMaxUp, unsigned long*& initialStateDown, int& initialStateLzMaxDown)
{
  unsigned long TmpState = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= initialStateLzMaxUp; ++i)
    {
      TmpState |= ((1ul << initialStateUp[i]) - 1ul) << Shift;
      Shift += initialStateUp[i];
      Shift += 2;
    }
  for (int i = 0; i <= initialStateLzMaxDown; ++i)
    {
      TmpState |= ((1ul << initialStateDown[i]) - 1ul) << Shift;
      Shift += initialStateDown[i];
      Shift += 2;
    }
  return TmpState;
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateLzMax = initial fermionic state maximum Lz value
// finalState = reference on the array where the bosonic state has to be stored
// finalStateLzMax = reference on the integer where the bosonic state maximum Lz value has to be stored

inline void BosonOnSphereTwoLandauLevels::FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState, int& finalStateLzMax)
{
  finalStateLzMax = 0;
  while (initialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialState - 1ul) ^ (~initialState);
      TmpState &= ~(TmpState >> 1);
//      cout << hex << initialState << "  " << TmpState << dec << endl;
#ifdef  __64_BITS__
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0f0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ffff0000ul) != 0) << 4;      
      TmpPower |= ((TmpState & 0xffffffff00000000ul) != 0) << 5;      
#else
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ul) != 0) << 4;      
#endif
//      cout << TmpPower << endl;
      finalState[finalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState >>= TmpPower;
      ++finalStateLzMax;
      initialStateLzMax -= TmpPower;
    }
  --finalStateLzMax;
}

#endif


