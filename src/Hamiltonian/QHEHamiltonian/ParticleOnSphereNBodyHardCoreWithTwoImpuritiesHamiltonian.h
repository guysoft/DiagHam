////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//  n-body hard core interaction and two localized impurities at the poles    //
//                                                                            //
//                        last modification : 04/05/2006                      //
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


#ifndef PARTICLEONSPHERENBODYHARDCOREWITHTWOIMPURITIESHAMILTONIAN_H
#define PARTICLEONSPHERENBODYHARDCOREWITHTWOIMPURITIESHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/QHEHilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereNBodyHardCoreHamiltonian.h"

#include <iostream>


using std::ostream;


class ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian : public ParticleOnSphereNBodyHardCoreHamiltonian
{

 private:

  // number of particle that interact simultaneously through the hard core interaction
  int NbrNbody;
  
  // weight of the different n-body interaction terms with respect to each other
  double* NBodyInteractionWeightFactors;

  // potential strength associted to the impurities
  double ImpurityPotential;

 public:

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // architecture = architecture to use for precalculation
  // nbrBody = number of particle that interact simultaneously through the hard core interaction
  // impurityPotential = potential strength associted to the impurities
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, int nbrBody, double impurityPotential,
							    AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
							    char* precalculationFileName = 0);

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // architecture = architecture to use for precalculation
  // maxNbrBody = maximum number of particle that interact simultaneously through the hard core interaction
  // nBodyFactors = weight of the different n-body interaction terms with respect to each other
  // impurityPotential = potential strength associted to the impurities
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
							 int maxNbrBody, double* nBodyFactors, double impurityPotential,
							 AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
							 char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();


 protected:
 
  // evaluate all interaction factors (including those arising from impurities)
  //     
  void EvaluateInteractionFactorsWithTwoImpurities();


};

#endif
