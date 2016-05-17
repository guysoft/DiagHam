////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                        Class author Cecile Repellin                        //
//                                                                            //
//                                                                            //
//              class of quasiholes on sphere with spin and pairing           //
//      (i.e. Sz conservation but no conservation of the particle number)     //
//                                                                            //
//                        last modification : 05/05/2016                      //
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


#ifndef QUASIHOLEONSPHEREWITHSPINANDPAIRING_H
#define QUASIHOLEONSPHEREWITHSPINANDPAIRING_H


#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "Matrix/SparseRealMatrix.h"

#include <iostream>


class QuasiholeOnSphereWithSpinAndPairing :  public FermionOnSphereWithSpin
{


 protected:
   // first index of the (k, r) exclusion principle in one layer
   int KValue;
   // second index of the (k, r) exclusion principle in one layer
   int RValue;
   // fermion factor (0 if boson, 1 if fermion)
   int FermionFactor;
   // number of admissible (NbrParticles, Lz) in one layer
   int NbrQuasiholeEntriesSingleLayer;
   // array containing the dimension of each quasihole subspace for a single layer
   int* NbrQuasiholesPerNPerLzSingleLayer;
   // array containing the first single layer index for each value of N and Lz
   int* SingleLayerIndices;
   // array containing the number of fermions in the up layer for each Hilbert space index
   int* NbrFermionUpFullSpace;
   // array containing the value of Lz in the up layer for each Hilbert space index
   int* LzValueUpFullSpace;
   // Hamiltonian matrix elements for quasiholes without spin
   SparseRealMatrix AnnihilationElementsOneLayer;
   // array of a^\dagger_m a_m matrix elements for quasiholes without spin
   SparseRealMatrix* AdAElementsOneLayer;
   // maximal number of fermions per layer
   int NbrFermionsUpMax;


 public:

  // default constructor
  //
  QuasiholeOnSphereWithSpinAndPairing();

  // basic constructor
  // 
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSpin = twice the total spin value
  // memory = amount of memory granted for precalculations
  QuasiholeOnSphereWithSpinAndPairing (int kExclusionPrinciple, int rExclusionPrinciple, int totalLz, int lzMax, int totalSpin, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // quasiholes = reference on the hilbert space to copy to copy
  QuasiholeOnSphereWithSpinAndPairing(const QuasiholeOnSphereWithSpinAndPairing& quasiholes);

  // destructor
  //
  ~QuasiholeOnSphereWithSpinAndPairing ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  QuasiholeOnSphereWithSpinAndPairing& operator = (const QuasiholeOnSphereWithSpinAndPairing& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();
  
  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);
  
  // apply a_u_m a_d_m to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for destruction operators
  // leftIndices = reference to an array containing the indices of the resulting states
  // interactionElements = reference to an array containing the matrix elements 
  // return value = number of left states that are connected to the initial state
  virtual int AuAd (int index, int m, int*& leftIndices, double*& interactionElements);
  
  // apply a^\dagger_u_m a^\dagger_d_m to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for creation operators
  // leftIndices = reference to an array containing the indices of the resulting states
  // interactionElements = reference to an array containing the matrix elements 
  // return value = number of left states that are connected to the initial state
  virtual int AduAdd (int index, int m, int*& leftIndices, double*& interactionElements);
  
  // apply a^\dagger_u_m a_u_m to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for destruction operators
  // leftIndices = reference to an array containing the indices of the resulting states
  // interactionElements = reference to an array containing the matrix elements 
  // return value = number of left states that are connected to the initial state
  virtual int AduAu (int index, int m, int*& leftIndices, double*& interactionElements);
  
  // apply a^\dagger_d_m a_d_m to a given state
  //
  // index = index of the state on which the operator has to be applied
  // m = index for destruction operators
  // leftIndices = reference to an array containing the indices of the resulting states
  // interactionElements = reference to an array containing the matrix elements 
  // return value = number of left states that are connected to the initial state
  virtual int AddAd (int index, int m, int*& leftIndices, double*& interactionElements);
  
  // get the number of particles in a given state
  //
  // index =index of the state whose number of particles has to be returned
  // return value = number of particles
  virtual int GetTotalNumberOfParticles (int index);


 protected:

  // evaluate Hilbert space dimension
  //
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // totalSpin = twice the total spin
  // return value = Hilbert space dimension      
  virtual long EvaluateHilbertSpaceDimension();


  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion in the state
  // totalLz = momentum total value
  // totalSpin = twice the total spin
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates();
  
  
  // get the linear index corresponding to a set of number of fermions and momentum
  //
  // NbrParticles = number of fermions
  // totalLz = value of the angular momentum
  // return value = linearized index
  virtual int GetLinearIndexSingleLayer(int nbrParticles, int totalLz);
  
  // get the maximal value of Lz in one layer for a given number of particles in this layer
  //
  // nbrParticles = number of particles
  // return value = maximal Lz
  virtual int GetMaximalLzSingleLayer(int nbrParticles);
  
  // find state index from the value of each number in one layer
  //
  // nbrParticlesUp = number of particles with up spin
  // lzValueUp = value of the angular momentum for up spins
  // alpha = index of state with up spin
  // beta = index of state with down spin
  // return value = corresponding index, -1 if an error occured
  virtual int FindStateIndex(int nbrParticlesUp, int lzValueUp, int alpha, int beta);
  
  // find the values of beta in each layer for a given state
  //
  // index = state index
  // alpha = reference on the value of alpha (up spin)
  // beta = reference on the value of beta (down spsin)
  virtual void FindBetaIndices(int index, int& alpha, int& beta);
  
};


  // get the linear index corresponding to a set of number of fermions and momentum
  //
  // NbrParticles = number of fermions
  // totalLz = value of the angular momentum
  // return value = linearized index
  inline int QuasiholeOnSphereWithSpinAndPairing::GetLinearIndexSingleLayer(int nbrParticles, int totalLz)
  {
    int MaxTotalLz = this->GetMaximalLzSingleLayer(nbrParticles);

    int TmpIndex = nbrParticles + this->LzMax * (nbrParticles - 1) * nbrParticles / 2 - (this->RValue + this->KValue * this->FermionFactor) * nbrParticles * (nbrParticles - 1) * (nbrParticles - 2) / 3 + (totalLz + MaxTotalLz) / 2;
    return TmpIndex;
  }
  
  
  // get the maximal value of Lz in one layer for a given number of particles in this layer
  //
  // nbrParticles = number of particles
  // return value = maximal Lz
  inline int QuasiholeOnSphereWithSpinAndPairing::GetMaximalLzSingleLayer(int nbrParticles)
  {
    int maxTotalLz = nbrParticles * this->LzMax - (((this->RValue + (this->KValue * this->FermionFactor)) * nbrParticles * (nbrParticles - 1)));
    return maxTotalLz;
        
  }

  
  // get the number of particles in a given state
  //
  // index =index of the state whose number of particles has to be returned
  // return value = number of particles
  inline int QuasiholeOnSphereWithSpinAndPairing::GetTotalNumberOfParticles (int index)
  {
    return (2 * this->NbrFermionUpFullSpace[index] - this->TotalSpin);
  }
#endif