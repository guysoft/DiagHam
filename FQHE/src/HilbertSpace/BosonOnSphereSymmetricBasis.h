////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of bosons on sphere using                     //
//                            the Lz <-> -Lz symmetry                         //
//                                                                            //
//                        last modification : 07/10/2007                      //
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


#ifndef BOSONONSPHERESYMMETRICBASIS_H
#define BOSONONSPHERESYMMETRICBASIS_H


#include "config.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


#define BOSON_SPHERE_SYMMETRIC_BIT  0x80000000ul
#define BOSON_SPHERE_SYMMETRIC_MASK 0x7ffffffful


class BosonOnSphereSymmetricBasis :  public BosonOnSphere
{

 protected:

  // shift to apply to a state before inverting its expression
  int InvertShift;
  // shift to apply to a state after inverting its expression
  int InvertUnshift;

  // temporary state and its LzMax value used during symmetrization operations
  int* TemporarySymmetrizedState;
  int TemporarySymmetrizedStateLzMax;
 
 public:

  // default constructor
  //
  BosonOnSphereSymmetricBasis ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // lzMax = maximum Lz value reached by a boson
  BosonOnSphereSymmetricBasis (int nbrBosons, int lzMax);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereSymmetricBasis(const BosonOnSphereSymmetricBasis& bosons);

  // destructor
  //
  virtual ~BosonOnSphereSymmetricBasis ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereSymmetricBasis& operator = (const BosonOnSphereSymmetricBasis& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  virtual AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter);

  // convert a given state from Lz-symmetric basis to the usual n-body basis
  //
  // state = reference on the vector to convert
  // nbodyBasis = reference on the nbody-basis to use
  // return value = converted vector
  RealVector ConvertToNbodyBasis(RealVector& state, BsosnOnSphere& nbodyBasis);

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
  //
  // index = index of the state on which the operator has to be applied
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient);

  // apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AA (int index, int n1, int n2);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int nbrIndices);

  // apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAd (int m1, int m2, double& coefficient);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int nbrIndices, double& coefficient);

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdA (int index, int m);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // evaluate wave function in real space using a given basis and only for a given range of components
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
					int firstComponent, int nbrComponent);

  // evaluate wave function in real space using a given basis and only for a given range of components, using time coherence
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // nextCoordinates = index of the coordinate that will be changed during the next time iteration
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunctionWithTimeCoherence (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
							 int nextCoordinates, int firstComponent, int nbrComponent);

  // initialize evaluation of wave function in real space using a given basis and only for a given range of components and
  //
  // timeCoherence = true if time coherence has to be used
  void InitializeWaveFunctionEvaluation (bool timeCoherence = false);

 protected:

  // find state index
  //
  // stateDescription = array describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = corresponding index
   virtual int FindStateIndex(int* stateDescription, int lzmax);

  // get canonical expression of a given state
  //
  // initialState = state that has to be converted to its canonical expression
  // return value = corresponding canonical state
  unsigned long GetCanonicalState (unsigned long initialState);

  // get symmetry of a given state 
  //
  // initialState = referennce state whose symmetry has to be computed
  void GetStateSymmetry (unsigned long& initialState);

  // get canonical expression of a given state and its symmetry
  //
  // initialState = array describing the state that has to be converted to its canonical expression
  // lzmax = maximum Lz value reached by a boson in initialState
  void GetSignedCanonicalState (int* initialState, int initialStateLzMax);

};

// get canonical expression of a given state and its symmetry
//
// initialState = array describing the state that has to be converted to its canonical expression
// lzmax = maximum Lz value reached by a boson in initialState

inline void BosonOnSphereSymmetricBasis::GetSignedCanonicalState (int* initialState, int initialStateLzMax)
{
  this->TemporarySymmetrizedStateLzMax = 0;
  while (initialState[this->TemporarySymmetrizedStateLzMax] == 0)
    ++this->TemporarySymmetrizedStateLzMax;
  for (int i = this->TemporarySymmetrizedStateLzMax; i <= initialStateLzMax; ++i)
    this->TemporarySymmetrizedState[this->LzMax - i] = initialState[i];
  for (int i = ; i >= 0; ++i)
    this->TemporarySymmetrizedState[this->LzMax - i] =0;
  this->TemporarySymmetrizedStateLzMax = this->LzMax - this->TemporarySymmetrizedStateLzMax;
  if (this->TemporarySymmetrizedStateLzMa != initialStateLzMax)
    this->TemporarySymmetrizedStateLzMax |= BOSON_SPHERE_SYMMETRIC_BIT;
  else
    {
    }
}

#endif


