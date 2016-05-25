////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                      Class author: Cecile Repellin                         //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to Laughlin qh on a sphere with      //
//        SU(2) spin with opposite magnetic field for each species            //
//                            and pairing                                     //
//                                                                            //
//                        last modification : 11/05/2016                      //
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


#ifndef PARTICLEONSPHEREWITHSPINTIMEREVERSALSYMMETRYQUASIHOLEHAMILTONIANANDPAIRING_H
#define PARTICLEONSPHEREWITHSPINTIMEREVERSALSYMMETRYQUASIHOLEHAMILTONIANANDPAIRING_H


#include "config.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairing.h"
#include "Hamiltonian/AbstractQHEOnSphereHamiltonian.h"

#include <iostream>


using std::ostream;


class MathematicaOutput;
class AbstractArchitecture;


class ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing : public AbstractQHEOnSphereHamiltonian
{

   friend class QHEParticlePrecalculationOperation;

 protected:

  // array that contains all one-body interaction factors for particles with spin up
  double* OneBodyInteractionFactorsupup;
  // array that contains all one-body interaction factors for particles with spin down
  double* OneBodyInteractionFactorsdowndown;
  // array that contains all one-body interaction factors for the pairing term
  double* OneBodyInteractionFactorsPairing;
  // factor in front of the charging energy (i.e 1/(2C))
  double ChargingEnergy;
  // avearge number of particles in the system
  double AverageNumberParticles;
  // global shift to apply to the diagonal matrix elements
  double HamiltonianShift;
  // array where the states coupled to a given state can be temporarily stored
  int* TmpLeftIndices;
  // array where the interaction elements coupling a given state to a group of other states can be temporarily stored
  double* TmpInteractionElements;

 public:

  ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // lzmax = maximum Lz value reached by a particle in the state
  // architecture = architecture to use for precalculation
  // pseudoPotential = array with the pseudo-potentials (sorted such that the first element corresponds to the delta interaction)
  //                   first index refered to the spin sector (sorted as up-up, down-down, up-down)
  // onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
  // onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
  // onebodyPotentialPairing =  one-body pairing term (sorted from component on the lowest Lz state to component on the highest Lz state), on site, symmetric spin up / spin down
  // chargingEnergy = factor in front of the charging energy (i.e 1/(2C))
  // averageNumberParticles = average number of particles in the system
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing(QuasiholeOnSphereWithSpinAndPairing* particles, int lzmax,
									    double* onebodyPotentialUpUp, double* onebodyPotentialDownDown,
									    double* onebodyPotentialPairing, double chargingEnergy, double averageNumberParticles,
									    AbstractArchitecture* architecture, long memory = -1, 
									    bool onDiskCacheFlag = false, char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing();

  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);
  
  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  virtual AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  virtual int GetHilbertSpaceDimension ();

 protected:
 
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and add result to another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector at which result has to be added
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
				  int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual RealVector* LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						  int firstComponent, int nbrComponent);

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // nbrComponent  = number of components that has to be precalcualted
  // return value = number of non-zero matrix element
  virtual long PartialFastMultiplicationMemory(int firstComponent, int nbrComponent);

  // firstComponent = index of the first component that has to be precalcualted
  // nbrComponent  = number of components that has to be precalcualted
  virtual void PartialEnableFastMultiplication(int firstComponent, int nbrComponent);
  
};

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension
inline int ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::GetHilbertSpaceDimension ()
{
  return (this->Particles->GetHilbertSpaceDimension());
}
  
  
#endif
