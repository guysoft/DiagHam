////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//      class of checkerboard lattice model with interacting particles        //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 03/04/2011                      //
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


#ifndef PARTICLEONLATTICEWITHSPINCHECKERBOARDLATTICESINGLEBANDHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINCHECKERBOARDLATTICESINGLEBANDHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeCheckerboardLatticeSingleBandHamiltonian : public ParticleOnLatticeChernInsulatorSingleBandHamiltonian
{

 protected:
  
  // hoping amplitude between neareast neighbor sites
  double NNHoping;
  // hoping amplitude between next neareast neighbor sites
  double NextNNHoping;
  // hoping amplitude between second next neareast neighbor sites
  double SecondNextNNHoping;
  // nearest neighbor density-density potential strength
  double UPotential;
  // boundary condition twisting angle along x
  double GammaX;
  // boundary condition twisting angle along y
  double GammaY;

  // use flat band model
  bool FlatBand;
  
 public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // t2p = hoping amplitude between second next neareast neighbor sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeCheckerboardLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential, double t1, double t2, double t2p, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeCheckerboardLatticeSingleBandHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body interaction between two sites A and B 
  //
  // kx1 = momentum along x for the A site
  // ky1 = momentum along y for the A site
  // kx2 = momentum along x for the B site
  // ky2 = momentum along y for the B site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementUpDown(int kx1, int ky1, int kx2, int ky2);

};



#endif
