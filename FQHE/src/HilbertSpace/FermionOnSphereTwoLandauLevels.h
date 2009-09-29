////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere including two                //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 19/05/2009                      //
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


#ifndef FERMIONONSPHERETWOLANDAULEVELS_H
#define FERMIONONSPHERETWOLANDAULEVELS_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"

#include <iostream>



class FermionOnSphereTwoLandauLevels :  public FermionOnSphereWithSpin
{

 protected:

  // maximum Lz value reached by a fermion with a spin up
  int LzMaxUp;
  // maximum Lz value reached by a fermion with a spin down
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
  FermionOnSphereTwoLandauLevels();

  // basic constructor with no contraint on the number of particles per spin component
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMaxUp = twice the maximum Lz value reached by a fermion with a spin up
  // lzMaxDown = twice the maximum Lz value reached by a fermion with a spin down
  // memory = amount of memory granted for precalculations
  FermionOnSphereTwoLandauLevels (int nbrFermions, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereTwoLandauLevels(const FermionOnSphereTwoLandauLevels& fermions);

  // destructor
  //
  ~FermionOnSphereTwoLandauLevels ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereTwoLandauLevels& operator = (const FermionOnSphereTwoLandauLevels& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // apply a^+_m_d a_m_d operator to a given state (only spin down)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AddAd (int index, int m);

  // apply a^+_m_u a_m_u operator to a given state  (only spin up)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AduAu (int index, int m);

  // apply a^+_m_u a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAu (int index, int m, int n, double& coefficient);

  // apply a^+_m_d a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAd (int index, int m, int n, double& coefficient);

  // apply a^+_m_u a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAd (int index, int m, int n, double& coefficient);

  // apply a^+_m_d a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAu (int index, int m, int n, double& coefficient);

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin up)
  // return value =  multiplicative factor 
  virtual double AuAu (int index, int n1, int n2);

  // apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin down)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double AdAd (int index, int n1, int n2);

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double AuAd (int index, int n1, int n2);

  // apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdu (int m1, int m2, double& coefficient);

  // apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAdd (int m1, int m2, double& coefficient);

  // apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdd (int m1, int m2, double& coefficient);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int* spinIndices, int nbrIndices);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // spinIndices = integer that gives the spin indices associated to each annihilation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int spinIndices, int nbrIndices);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // spinIndices = integer that gives the spin indices associated to each creation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int spinIndices, int nbrIndices, double& coefficient);

  // create an SU(2) state from two U(1) state
  //
  // upState = vector describing the up spin part of the output state
  // upStateSpace = reference on the Hilbert space associated to the up spin part
  // downState = vector describing the down spin part of the output state
  // downStateSpace = reference on the Hilbert space associated to the down spin part  
  // return value = resluting SU(2) state
  virtual RealVector ForgeSU2FromU1(RealVector& upState, FermionOnSphere& upStateSpace, RealVector& downState, FermionOnSphere& downStateSpace);

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermionsUp = number of fermions with spin up
  // nbrFermionsDown = number of fermions with spin down
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension      
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermionsUp, int nbrFermionsDown, int lzMax, int totalLz);

  // evaluate Hilbert space dimension without constraint on the number of particles per level
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateFullHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

  // generate all states corresponding to the constraints
  // 
  // nbrFermionsUp = number of fermions with spin up
  // nbrFermionsDown = number of fermions with spin down
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermionsUp, int nbrFermionsDown, int lzMax, int totalLz, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMaxUp = momentum maximum value for a fermion
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateFullStates(int nbrFermions, int lzMax, int totalLz, long pos);

};

#endif


