////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//              class of fermions on sphere with SU(3) spin including         //
//                          Tz<->-Tz discrete symmetry                        //
//                                                                            //
//                        last modification : 04/02/2008                      //
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


#ifndef FERMIONONSPHEREWITHSU3SPIN_H
#define FERMIONONSPHEREWITHSU3SPIN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSU3Spin.h"

#include <iostream>


class FermionOnSphere;


class FermionOnSphereWithSU3SpinTzSymmetricBasis :  public ParticleOnSphereWithSU3Spin
{

  // number of fermions
  int NbrFermions;
  // number of fermions plus 1
  int IncNbrFermions;
  // momentum total value
  int TotalLz;
  // maximum Lz value reached by a fermion
  int LzMax;
  // number of Lz values in a state
  int NbrLzValue;
  // twice the total Tz value
  int TotalTz;
  // three time the total Y value
  int TotalY;

  // array describing each state
  unsigned long* StateDescription;
  // array giving maximum Lz value reached for a fermion in a given state
  int* StateHighestBit;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given lzmax sector
  unsigned long LookUpTableMemorySize;
  // shift used in each lzmax sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  int** LookUpTable;

  // a table containing ranging from 0 to 2^MaximumSignLookUp - 1
  double* SignLookUpTable;
  // a table containing the mask on the bits to keep for each shift that is requested by sign evaluation
  unsigned long* SignLookUpTableMask;
  // number to evalute size of SignLookUpTable
  int MaximumSignLookUp;

  // temporary state used when applying ProdA operator
  unsigned long ProdATemporaryState;
  // Lz maximum value associated to temporary state used when applying ProdA operator
  int ProdALzMax;
  // symmetry of the temporary state
  int ProdASymmetry;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalY = three time the total Y value
  // minusTzParity = select the  Tz <-> -Tz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSU3SpinTzSymmetricBasis (int nbrFermions, int totalLz, int lzMax, int totalY, bool minusTzParity, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSU3SpinTzSymmetricBasis(const FermionOnSphereWithSU3SpinTzSymmetricBasis& fermions);

  // destructor
  //
  ~FermionOnSphereWithSU3SpinTzSymmetricBasis ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSU3SpinTzSymmetricBasis& operator = (const FermionOnSphereWithSU3SpinTzSymmetricBasis& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // apply a^+_m_1 a_m_1 operator to a given state (only state 1 Tz=+1/2, Y=+1/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_1 a_m_1
  double Ad1A1 (int index, int m);

  // apply a^+_m_2 a_m_2 operator to a given state (only state 2 Tz=-1/2, Y=+1/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_2 a_m_2
  double Ad2A2 (int index, int m);

  // apply a^+_m_3 a_m_3 operator to a given state (only state 3 Tz=0, Y=-2/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_3 a_m_3
  double Ad3A3 (int index, int m);

  // apply a^+_m1_1 a^+_m2_1 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1Ad1 (int m1, int m2, double& coefficient);

  // apply a^+_m1_1 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1Ad2 (int m1, int m2, double& coefficient);

  // apply a^+_m1_1 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1Ad3 (int m1, int m2, double& coefficient);

  // apply a^+_m1_2 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad2Ad2 (int m1, int m2, double& coefficient);

  // apply a^+_m1_2 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad2Ad3 (int m1, int m2, double& coefficient);

  // apply a^+_m1_3 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad3Ad3 (int m1, int m2, double& coefficient);

  protected:

  // get canonical expression of a given state
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state
  virtual unsigned long GetCanonicalState (unsigned long initialState);

  // get symmetry of a given state 
  //
  // initialState = referennce state whose symmetry has to be computed
  virtual void GetStateSymmetry (unsigned long& initialState);

  // get canonical expression of a given state and its symmetry
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state (with symmetry bit)
  virtual unsigned long GetSignedCanonicalState (unsigned long initialState);

};

// get canonical expression of a given state
//
// initialState = state that has to be converted to its canonical expression
// return value = corresponding canonical state

inline unsigned long FermionOnSphereWithSpinSzSymmetry::GetCanonicalState (unsigned long initialState)
{
  unsigned long TmpMask =  ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpMask |= TmpMask << 1;
  if ((initialState ^ TmpMask) < initialState)
    return (initialState ^ TmpMask);
  else
    return initialState;
}

// get symmetry of a given state 
//
// initialState = reference on the state whose symmetry has to be computed
// return value = corresponding symmetry bit

unsigned long void FermionOnSphereWithSpinSzSymmetry::GetStateSymmetry (unsigned long& initialState)
{
  unsigned long TmpMask =  ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU3_TZ_MASK;
  TmpMask |= TmpMask << 1;
  if ((initialState ^ TmpMask) != initialState)    
    return 0ul;
  else
    return FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT;
}

// get canonical expression of a given state and its symmetry
//
// initialState = reference on the state that has to be converted to its canonical expression
// return value = corresponding symmetry bit

inline unsigned int FermionOnSphereWithSpinSzSymmetry::GetSignedCanonicalState (unsigned long& initialState)
{
  unsigned long TmpMask =  ((initialState >> 1) ^ initialState) & FERMION_SPHERE_SU2_TZ_MASK;
  TmpMask |= TmpMask << 1;
  if ((initialState ^ TmpMask) < initialState)
    {
      initialState ^= TmpMask;
      return 0x0;
    }
  else
    if ((initialState ^ TmpMask) != initialState)
      return FERMION_SPHERE_SU3_TZ_SYMMETRIC_BIT;
    else
      return 0x0u;
}

#endif


