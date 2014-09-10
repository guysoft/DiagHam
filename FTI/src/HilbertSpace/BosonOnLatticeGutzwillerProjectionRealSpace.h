////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of hardcore boson  on lattice in real space       //
//                                                                            //
//                        last modification : 10/09/2014                      //
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


#ifndef BOSONONLATTICEGUTZWILLERPROJECTIONREALSPACE_H
#define BOSONONLATTICEGUTZWILLERPROJECTIONREALSPACE_H

#include "config.h"
#include "HilbertSpace/FermionOnLatticeRealSpace.h"

#include <iostream>



class BosonOnLatticeGutzwillerProjectionRealSpace : public FermionOnLatticeRealSpace
{

 protected:

  // number of bosons
  int NbrBosons;
  int IncNbrBosons;

  // indices of the orbitals that are kept when performing an orbital cut
  int* KeptOrbitals;

 public:

  // default constructor
  // 
  BosonOnLatticeGutzwillerProjectionRealSpace ();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSite = number of sites
  // memory = amount of memory granted for precalculations
  BosonOnLatticeGutzwillerProjectionRealSpace (int nbrBosons, int nbrSite, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  BosonOnLatticeGutzwillerProjectionRealSpace(const BosonOnLatticeGutzwillerProjectionRealSpace& bosons);

  // destructor
  //
  ~BosonOnLatticeGutzwillerProjectionRealSpace ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnLatticeGutzwillerProjectionRealSpace & operator = (const BosonOnLatticeGutzwillerProjectionRealSpace& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // evaluate the orbital cut entanglement matrix. The entanglement matrix is only evaluated for fixed number of particles
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // keptOrbitals = array of orbitals that have to be kept, should be sorted from the smallest index to the largest index 
  // nbrKeptOrbitals = array of orbitals that have to be kept
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrParticleSector, int nbrKeptOrbitals, int* keptOrbitals, ComplexVector& groundState, AbstractArchitecture* architecture = 0);
  
  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
  // 
  // nbrParticleSector = number of particles that belong to the subsytem 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

 protected:

  // core part of the evaluation orbital cut entanglement matrix calculation
  // 
  // minIndex = first index to consider in source Hilbert space
  // nbrIndex = number of indices to consider in source Hilbert space
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  // return value = number of components that have been added to the density matrix
  virtual long EvaluatePartialEntanglementMatrixCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
						      ComplexVector& groundState, ComplexMatrix* entanglementMatrix);

};


#endif


