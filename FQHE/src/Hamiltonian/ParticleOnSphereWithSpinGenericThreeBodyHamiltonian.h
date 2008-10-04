////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                     spin and generic 3-body interaction                    //
//                                                                            //
//                        last modification : 26/08/2008                      //
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


#ifndef PARTICLEONSPHEREGENERICTHREEBODYHAMILTONIAN_H
#define PARTICLEONSPHEREGENERICTHREEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "Hamiltonian/AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian.h"

#include <iostream>


using std::ostream;


class ClebschGordanCoefficients;


class ParticleOnSphereWithSpinGenericThreeBodyHamiltonian : public AbstractQHEOnSphereWithSpinNBodyInteractionHamiltonian
{

 protected:

  
  // array with the three-body pseudo-potentials between spin up - spin up, sorted with respect to the relative angular momentum, 
  // taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  double* ThreeBodyPseudoPotentialUpUp;
  // nuber of elements in the ThreeBodyPseudoPotential array
  int NbrThreeBodyPseudoPotentialUpUp;
  // maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential
  int MaxRelativeAngularMomentumUpUp;

  // array with the three-body pseudo-potentials between spin down - spin down, sorted with respect to the relative angular momentum, 
  // taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  double* ThreeBodyPseudoPotentialDownDown;
  // nuber of elements in the ThreeBodyPseudoPotential array
  int NbrThreeBodyPseudoPotentialDownDown;
  // maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential
  int MaxRelativeAngularMomentumDownDown;

  // array with the three-body pseudo-potentials between spin up - spin down, sorted with respect to the relative angular momentum, 
  // taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  double* ThreeBodyPseudoPotentialUpDown;
  // nuber of elements in the ThreeBodyPseudoPotential array
  int NbrThreeBodyPseudoPotentialUpDown;
  // maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential
  int MaxRelativeAngularMomentumUpDown;

  // array with the pseudo-potentials (ordered such that the last element corresponds to the delta interaction)
  double* PseudoPotential;

 public:

  // default constructor
  //
  ParticleOnSphereWithSpinGenericThreeBodyHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // threeBodyPseudoPotential = array with the three-body pseudo-potentials sorted with respect to the relative angular momentum, 
  //                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  // maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereWithSpinGenericThreeBodyHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, 
						      double* threeBodyPseudoPotential, int maxRelativeAngularMomentum, double l2Factor, 
						      AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
						      char* precalculationFileName = 0);

  // constructor from datas with a fully-defined two body interaction
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // threeBodyPseudoPotential = array with the three-body pseudo-potentials sorted with respect to the relative angular momentum, 
  //                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
  // maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in ThreeBodyPseudoPotential
  // l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
  // pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereWithSpinGenericThreeBodyHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int lzmax, 
						      double* threeBodyPseudoPotential, int maxRelativeAngularMomentum,
						      double l2Factor, double* pseudoPotential, 
						      AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
						      char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereWithSpinGenericThreeBodyHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();


 protected:
 
  // compute all projector coefficient associated to a given relative angular momentum between 3 particles
  //
  // relativeMomentum = value of twice the relative angular momentum between the 3 particles
  // degeneracyIndex = optional degeneracy index for relative angular momentum greater than 5 for bosons (8 for fermions)
  // indices = array that contains all possible sets of indices (size of the array is 3 * nbrIndexSets)
  // nbrIndexSets = number of sets
  double* ComputeProjectorCoefficients(int relativeMomentum, int degeneracyIndex, int* indices, int nbrIndexSets);

  // evaluate all interaction factors
  //   
  void EvaluateInteractionFactors();


};

#endif
