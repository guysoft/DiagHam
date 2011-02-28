////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                   class of quatum spin Hall restricted to two band         //
//                                                                            //
//                        last modification : 27/02/2011                      //
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


#ifndef PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDHAMILTONIAN_H
#define PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinChernInsulatorHamiltonian.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian : public ParticleOnLatticeWithSpinChernInsulatorHamiltonian
{

 protected:
  
  // szSymmetryBreaking = amplitude of the Sz symmetry breaking term
  double SzSymmetryBreaking;

 public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // bandParameter = band parameter
  // szSymmetryBreaking = amplitude of the Sz symmetry breaking term
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double bandParameter, double szSymmetryBreaking, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the part of the interaction coefficient coming from the two band truncated basis
  //
  // basisMatrices = array of basis matrices
  // momentumIndex1 = momentum index for the first particle
  // momentumIndex2 = momentum index for the second particle
  // bandIndex1 = band index of the first particle
  // bandIndex2 = band index of the second particle
  // return value = coefficient
  Complex ComputeBasisContribution(ComplexMatrix* basisMatrices, int momentumIndex1, int momentumIndex2, int bandIndex1, int bandIndex2);

};



#endif
