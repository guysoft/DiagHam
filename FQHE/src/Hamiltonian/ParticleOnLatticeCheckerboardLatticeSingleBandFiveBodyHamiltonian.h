////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//       class of checkerboard lattice model with interacting particles       //
//         in the single band approximation and five body interaction         // 
//                                                                            //
//                        last modification : 07/08/2011                      //
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


#ifndef PARTICLEONLATTICEWITHSPINCHECKERBOARDLATTICESINGLEBANDFIVEBODYHAMILTONIAN_H
#define PARTICLEONLATTICEWITHSPINCHECKERBOARDLATTICESINGLEBANDFIVEBODYHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Hamiltonian/ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian.h"
#include "Vector/ComplexVector.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian : public ParticleOnLatticeCheckerboardLatticeSingleBandFourBodyHamiltonian
{

 protected:
  
  
 public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = strength of the repulsive five body neareast neighbor interaction
  // vPotential = strength of the repulsive two body neareast neighbor interaction
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // t2p = hoping amplitude between second next neareast neighbor sites
  // mus = sublattice staggered chemical potential 
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential, double vPotential, double t1, double t2, double t2p, double mus, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeCheckerboardLatticeSingleBandFiveBodyHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the five body interaction between one four A and three sites B 
  //
  // kx1 = momentum along x of the creation operator of the A site
  // ky1 = momentum along y of the creation operator of the A site
  // kx2 = momentum along x of the creation operator of the first B site
  // ky2 = momentum along y of the creation operator of the first B site
  // kx3 = momentum along x of the creation operator of the second B site
  // ky3 = momentum along y of the creation operator of the second B site
  // kx4 = momentum along x of the creation operator of the third B site
  // ky4 = momentum along y of the creation operator of the third B site
  // kx5 = momentum along x of the creation operator of the fourth B site
  // ky5 = momentum along y of the creation operator of the fourth B site
  // kx6 = momentum along x of the creation operator of the A site
  // ky6 = momentum along y of the creation operator of the A site
  // kx7 = momentum along x of the creation operator of the first B site
  // ky7 = momentum along y of the creation operator of the first B site
  // kx8 = momentum along x of the creation operator of the second B site
  // ky8 = momentum along y of the creation operator of the second B site
  // kx9 = momentum along x of the creation operator of the third B site
  // ky9 = momentum along y of the creation operator of the third B site
  // kx10 = momentum along x of the creation operator of the fourth B site
  // ky10 = momentum along y of the creation operator of the fourth B site
  // return value = corresponding matrix element
  Complex ComputeFiveBodyMatrixElementABBBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8, int kx9, int ky9, int kx10, int ky10);

  // compute the matrix element for the five body interaction between four sites A and one site B 
  //
  // kx1 = momentum along x of the creation operator of the B site
  // ky1 = momentum along y of the creation operator of the B site
  // kx2 = momentum along x of the creation operator of the first A site
  // ky2 = momentum along y of the creation operator of the first A site
  // kx3 = momentum along x of the creation operator of the second A site
  // ky3 = momentum along y of the creation operator of the second A site
  // kx4 = momentum along x of the creation operator of the third A site
  // ky4 = momentum along y of the creation operator of the third A site
  // kx5 = momentum along x of the creation operator of the fourth A site
  // ky5 = momentum along y of the creation operator of the fourth A site
  // kx6 = momentum along x of the creation operator of the B site
  // ky6 = momentum along y of the creation operator of the B site
  // kx7 = momentum along x of the creation operator of the first A site
  // ky7 = momentum along y of the creation operator of the first A site
  // kx8 = momentum along x of the creation operator of the second A site
  // ky8 = momentum along y of the creation operator of the second A site
  // kx9 = momentum along x of the creation operator of the third A site
  // ky9 = momentum along y of the creation operator of the third A site
  // kx10 = momentum along x of the creation operator of the fourth A site
  // ky10 = momentum along y of the creation operator of the fourth A site
  // return value = corresponding matrix element
  Complex ComputeFiveBodyMatrixElementBAAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6, int kx7, int ky7, int kx8, int ky8, int kx9, int ky9, int kx10, int ky10);


};



#endif
