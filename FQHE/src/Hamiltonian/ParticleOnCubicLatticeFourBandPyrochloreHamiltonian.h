////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of 3d topological insulator based on the pryochlore       //
//                       model and restricted to four bands                   //
//                                                                            //
//                        last modification : 13/08/2012                      //
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


#ifndef PARTICLEONCUBICLATTICEFOURBANDPYROCHLOREHAMILTONIAN_H
#define PARTICLEONCUBICLATTICEFOURBANDPYROCHLOREHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian.h"
#include "Tools/FTITightBinding/Abstract3DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnCubicLatticeFourBandPyrochloreHamiltonian : public ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian
{

 protected:
 
  // number of sites in the z direction
  int NbrSiteZ;
  // number of sites in the direction perpendicular to X
  int NbrSiteYZ;

  // repulsive on-site potential strength between identical spins
  double UPotential;
  // repulsive on-site potential strength between opposite spins
  double VPotential;

  // pointer to the tight binding model
  Abstract3DTightBindingModel* TightBindingModel;

  // use flat band model
  bool FlatBand;


 public:

  // default constructor
  //
  ParticleOnCubicLatticeFourBandPyrochloreHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // uPotential = repulsive on-site potential strength between identical spins
  // vPotential = repulsive on-site potential strength between opposite spins
  // tightBindingModel = pointer to the tight binding model
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnCubicLatticeFourBandPyrochloreHamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, int nbrSiteZ, double uPotential, double vPotential,  Abstract3DTightBindingModel* tightBindingModel, 
						      bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnCubicLatticeFourBandPyrochloreHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();
  


};

#endif
