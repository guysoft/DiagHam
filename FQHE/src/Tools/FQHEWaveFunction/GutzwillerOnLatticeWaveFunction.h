////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                       Copyright (C) 2009 Gunnar Möller                     //
//                                                                            //
//                                                                            //
//           class for calculation of a Gutzwiller state on the lattice       //
//                                                                            //
//                        last modification : 27/10/2009                      //
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


#ifndef GUTZWILLERWAVEFUNCTION_H
#define GUTZWILLERWAVEFUNCTION_H

class AbstractQHEOnLatticeHamiltonian;
class AbstractArchitecture;

#include "Vector/ComplexVector.h"
#include "HilbertSpace/ParticleOnLattice.h"

class GutzwillerOnLatticeWaveFunction
{
 protected:

  // target state, internal use
  ComplexVector TargetVector;

  // number of particles in many-body state
  int NbrParticles;

  // number of sites on lattice
  int NbrSites;

  // flag for hardcore condition
  bool HardCoreFlag;

  // dimension of Space
  int Dim;

  // number of variational parameters
  int NbrVariationalParameters;
  
  // vector with variational parameters
  // order: {(a^0_i),(a^1_i),...}, i=0...N_s
  ComplexVector VariationalParameters;

  // target Hilbert space
  ParticleOnLattice *Space;

  // Hamiltonian
  AbstractQHEOnLatticeHamiltonian *Hamiltonian;

  // Architecture
  AbstractArchitecture *Architecture;

 public:
  // constructor
  // nbrParticles = particles in the condensate (should match space)
  // hardCore = flag indicating whether double occupations may occur or whether hardcore particles are present
  // space = target-space of many-body state
  // variationalParameters = initial set of trial parameters
  GutzwillerOnLatticeWaveFunction(int nbrParticles, bool hardCore, ParticleOnLattice *space, ComplexVector *variationalParameters=NULL);

  // destructor
  ~GutzwillerOnLatticeWaveFunction();

  // get Many-Body state
  // return = resultingState
  ComplexVector & GetGutzwillerWaveFunction();

  // set trial parameters
  void SetVariationalParameters(ComplexVector &variationalParameters);

  // define a Hamiltonian to enable immediate evaluation of the energy
  void SetHamiltonian(AbstractQHEOnLatticeHamiltonian *hamiltonian);
  
  // define an architecture to enable multi-processor operations
  void SetArchitecture(AbstractArchitecture *architecture);

 protected:
 
  // main recursion to calculate State \prod_i[\sum_k \sum a^k_i a^\dagger_i)^k] |state>
  // nextQ = value quantum number in next operator to be applied
  // nbrBosons = number of bosons already in state
  // state = state to be acted upon
  // prefactor = previous coefficients applied to state
  // in last stage of recursion, writes to this->TargetVector
  void Product (int nextQ, int nbrBosons, unsigned long state, Complex prefactor);


};


#endif
