////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//              class of fermions on sphere with SU(3) spin including         //
//                           the Z3 discrete symmetry                         //
//                                                                            //
//                        last modification : 09/02/2008                      //
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


#ifndef FERMIONONSPHEREWITHSU3SPINZ3SYMMETRY_H
#define FERMIONONSPHEREWITHSU3SPINZ3SYMMETRY_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSU3SpinTzSymmetry.h"

#include <iostream>


class FermionOnSphereWithSU3SpinZ3Symmetry :  public FermionOnSphereWithSU3SpinTzSymmetry
{


 public:

  // default constructor
  // 
  FermionOnSphereWithSU3SpinZ3Symmetry();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalTz = twice the total Tz value
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSU3SpinZ3Symmetry (int nbrFermions, int totalLz, int lzMax, int totalTz, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSU3SpinZ3Symmetry(const FermionOnSphereWithSU3SpinZ3Symmetry& fermions);

  // constructor from a binary file that describes the Hilbert space
  //
  // fileName = name of the binary file
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSU3SpinZ3Symmetry (char* fileName, unsigned long memory);

  // destructor
  //
  ~FermionOnSphereWithSU3SpinZ3Symmetry ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSU3SpinZ3Symmetry& operator = (const FermionOnSphereWithSU3SpinZ3Symmetry& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

 protected:

  // get canonical expression of a given state
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state
  virtual unsigned long GetCanonicalState (unsigned long initialState);

  // get symmetry of a given state 
  //
  // initialState = referennce state whose symmetry has to be computed
  // return value = corresponding symmetry bit
  virtual unsigned long GetStateSymmetry (unsigned long& initialState);

  // compute the parity of the number of spin singlet aka the number of particle pairs (no more than two) per Lz value
  //
  // initialState = reference on the state whose parity has to be evaluated
  // return value = corresponding parity bit
  virtual unsigned long GetSignedCanonicalState (unsigned long& initialState);

  // compute the parity of the number of spin singlet aka the number of particle pairs (no more than two) per Lz value
  //
  // initialState = reference on the state whose parity has to be evaluated
  // return value = corresponding parity bit
  unsigned long GetStateSingletYParity(unsigned long& initialState);

  // factorized code that is used to symmetrize result of the AdxAdy operations
  //
  // state = reference on the state that has been produced with the AdxAdy operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state
  virtual int SymmetrizeAdAdResult(unsigned long& state, double& coefficient);

};

// get canonical expression of a given state
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state

inline unsigned long FermionOnSphereWithSU3SpinZ3Symmetry::GetCanonicalState (unsigned long initialState)
{
  unsigned long TmpState =  (((initialState & FERMION_SPHERE_SU3_1_MASK) << 2) |
			     ((initialState & FERMION_SPHERE_SU3_23_MASK) >> 1));
  if (TmpState < initialState)
    initialState = TmpState;
  TmpState =  (((TmpState & FERMION_SPHERE_SU3_1_MASK) << 2) |
	       ((TmpState & FERMION_SPHERE_SU3_23_MASK) >> 1));
  if (TmpState < initialState)
    return TmpState;
  else
    return initialState;
}

// get symmetry of a given state 
//
// initialState = reference on the state whose symmetry has to be computed
// return value = corresponding symmetry bit

inline unsigned long FermionOnSphereWithSU3SpinZ3Symmetry::GetStateSymmetry (unsigned long& initialState)
{
  unsigned long TmpState =  (((initialState & FERMION_SPHERE_SU3_1_MASK) << 2) |
			     ((initialState & FERMION_SPHERE_SU3_23_MASK) >> 1));
  if (TmpState == initialState)
    return FERMION_SPHERE_SU3_Z3_SYMMETRIC_BIT;
  else
    return 0ul;
}

// get canonical expression of a given state and its symmetry
//
// initialState = reference on the state that has to be converted to its canonical expression
// return value = corresponding symmetry bit

inline unsigned long FermionOnSphereWithSU3SpinZ3Symmetry::GetSignedCanonicalState (unsigned long& initialState)
{  
  unsigned long TmpState =  (((initialState & FERMION_SPHERE_SU3_1_MASK) << 2) |
			     ((initialState & FERMION_SPHERE_SU3_23_MASK) >> 1));
  if (TmpState == initialState)
    return FERMION_SPHERE_SU3_Z3_SYMMETRIC_BIT;
  if (TmpState < initialState)
    initialState = TmpState;
  TmpState =  (((TmpState & FERMION_SPHERE_SU3_1_MASK) << 2) |
	       ((TmpState & FERMION_SPHERE_SU3_23_MASK) >> 1));
  if (TmpState < initialState)
    {
      initialState = TmpState;
      return 0x0ul;
    }
  else
    {
      return 0x0ul;
    }
}

// compute the parity of the number of spin singlet aka the number of particle pairs (no more than two) per Lz value
//
// initialState = reference on the state whose parity has to be evaluated
// return value = corresponding parity bit

inline unsigned long FermionOnSphereWithSU3SpinZ3Symmetry::GetStateSingletYParity(unsigned long& initialState)
{
  unsigned long TmpState = initialState;
  TmpState ^= initialState >> 1;
  TmpState ^= initialState >> 2;
  TmpState ^= ((initialState >> 2) | (initialState >> 1) | initialState);
  TmpState &= FERMION_SPHERE_SU3_1_MASK;
#ifdef __64_BITS__
  TmpState ^= (TmpState >> 32);
#endif
  TmpState ^= (TmpState >> 16);
  TmpState ^= (TmpState >> 8);
  TmpState ^= (TmpState >> 4);
  TmpState ^= (TmpState >> 2);
  TmpState ^= (TmpState >> 1);
  return ((TmpState & 1) << FERMION_SPHERE_SU3_Z3SINGLETPARITY_SHIFT);
}

// factorized code that is used to symmetrize result of the AdxAdy operations
//
// state = reference on the state that has been produced with the AdxAdy operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int FermionOnSphereWithSU3SpinZ3Symmetry::SymmetrizeAdAdResult(unsigned long& state, double& coefficient)
{
  unsigned long TmpState2 = state;
  if ((this->GetSignedCanonicalState(state) & FERMION_SPHERE_SU3_Z3_SYMMETRIC_BIT) != 0x0ul)
    {
      if ((this->ProdASignature & FERMION_SPHERE_SU3_Z3_SYMMETRIC_BIT) == 0x0ul)
	coefficient *= sqrt(3.0);
      int NewLzMax = 2 + (this->LzMax * 3);
      while ((state >> NewLzMax) == 0x0ul)
	--NewLzMax;
      return this->FindStateIndex(state, NewLzMax);
    }
  int NewLzMax = 2 + (this->LzMax * 3);
  while ((state >> NewLzMax) == 0x0ul)
    --NewLzMax;
  if (TmpState2 != state)
    {
      coefficient *= (1.0 - 2.0 * ((double) ((this->GetStateSingletYParity(TmpState2) >> FERMION_SPHERE_SU3_Z3SINGLETPARITY_SHIFT) & 0x1ul)));
    }
  if ((this->ProdASignature & FERMION_SPHERE_SU3_Z3_SYMMETRIC_BIT) != 0x0ul)
    coefficient *= 1.0 / sqrt(3.0);
  return this->FindStateIndex(state, NewLzMax);
}
#endif


