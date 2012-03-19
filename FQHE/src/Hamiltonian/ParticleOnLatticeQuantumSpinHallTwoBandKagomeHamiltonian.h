////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                class of quatum spin Hall restricted to two bands           //
//                           using the kagome model                           //
//                                                                            //
//                        last modification : 17/03/2012                      //
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


#ifndef PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDKAGOMEHAMILTONIAN_H
#define PARTICLEONLATTICEQUANTUMSPINHALLTWOBANDKAGOMEHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnSquareLatticeTwoBandSimpleTIHamiltonian.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian : public ParticleOnSquareLatticeTwoBandSimpleTIHamiltonian
{

 protected:
  
  // hopping amplitude between neareast neighbor sites
  double NNHopping;
  // hopping amplitude between next neareast neighbor sites
  double NextNNHopping;
  // spin orbit coupling to neareast neighbor sites
  double NNSpinOrbit;
  // spin orbit coupling to next neareast neighbor sites
  double NextNNSpinOrbit;
  
  // mixing term coupling the two copies of the kagome lattice
  Complex MixingTerm;

  // four times the sublattice staggered chemical potential 
  double MuS;
  // boundary condition twisting angle along x
  double GammaX;
  // boundary condition twisting angle along y
  double GammaY;
  // nearest neighbor density-density potential strength
  double UPotential;
  // strength of the repulsive on site two body interaction between opposite spins
  double VPotential;
  // strength of the repulsive two body neareast neighbor interaction between opposite spins
  double WPotential;
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
  // vPotential = strength of the repulsive on site two body interaction between opposite spins
  // wPotential = strength of the repulsive two body neareast neighbor interaction between opposite spins
  // t1 = real part of the hopping amplitude between neareast neighbor sites
  // t2 = real part of the hopping amplitude between next neareast neighbor sites
  // lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
  // lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
  // mixingTermNorm = norm of the mixing term coupling the two copies of the kagome lattice
  // mixingTermArgv = argument of the mixing term coupling the two copies of the kagome lattice
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential, double vPotential, double wPotential, double t1, double t2, double lambda1, double lambda2, double mixingTermNorm, double mixingTermArg, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  //  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body interaction between two sites A and B  belonging to the same layer
  //
  // kx1 = momentum along x for the creation operator on A site with spin up
  // ky1 = momentum along y for the creation operator on A site with spin up
  // kx2 = momentum along x for the creation operator on B site with spin up
  // ky2 = momentum along y for the creation operator on B site with spin up
  // kx3 = momentum along x for the annihilation operator on A site with spin up
  // ky3 = momentum along y for the annihilation operator on A site with spin up
  // kx4 = momentum along x for the annihilation operator on B site with spin up
  // ky4 = momentum along y for the annihilation operator on B site with spin up
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two sites A and B with down spins
  //
  // kx1 = momentum along x for the creation operator on A site with spin down
  // ky1 = momentum along y for the creation operator on A site with spin down
  // kx2 = momentum along x for the creation operator on B site with spin down
  // ky2 = momentum along y for the creation operator on B site with spin down
  // kx3 = momentum along x for the annihilation operator on A site with spin down
  // ky3 = momentum along y for the annihilation operator on A site with spin down
  // kx4 = momentum along x for the annihilation operator on B site with spin down
  // ky4 = momentum along y for the annihilation operator on B site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementADownBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);
  
  // compute the matrix element for the two body interaction between two sites A and B with opposite spins
  //
  // kx1 = momentum along x for the creation operator on A site with spin down
  // ky1 = momentum along y for the creation operator on A site with spin down
  // kx2 = momentum along x for the creation operator on B site with spin up
  // ky2 = momentum along y for the creation operator on B site with spin up
  // kx3 = momentum along x for the annihilation operator on A site with spin down
  // ky3 = momentum along y for the annihilation operator on A site with spin down
  // kx4 = momentum along x for the annihilation operator on B site with spin up
  // ky4 = momentum along y for the annihilation operator on B site with spin up
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementADownBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two sites A and B with opposite spins
  //
  // kx1 = momentum along x for the creation operator on A site with spin up
  // ky1 = momentum along y for the creation operator on A site with spin up
  // kx2 = momentum along x for the creation operator on B site with spin down
  // ky2 = momentum along y for the creation operator on B site with spin down
  // kx3 = momentum along x for the annihilation operator on A site with spin up
  // ky3 = momentum along y for the annihilation operator on A site with spin up
  // kx4 = momentum along x for the annihilation operator on B site with spin down
  // ky4 = momentum along y for the annihilation operator on B site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementAUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two sites A with opposite spins 
  //
  // kx1 = momentum along x for the creation operator on A site with spin up
  // ky1 = momentum along y for the creation operator on A site with spin up
  // kx2 = momentum along x for the creation operator on A site with spin down
  // ky2 = momentum along y for the creation operator on A site with spin down
  // kx3 = momentum along x for the annihilation operator on A site with spin up
  // ky3 = momentum along y for the annihilation operator on A site with spin up
  // kx4 = momentum along x for the annihilation operator on A site with spin down
  // ky4 = momentum along y for the annihilation operator on A site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);

  // compute the matrix element for the two body interaction between two sites B with opposite spins 
  //
  // kx1 = momentum along x for the creation operator on B site with spin up
  // ky1 = momentum along y for the creation operator on B site with spin up
  // kx2 = momentum along x for the creation operator on B site with spin down
  // ky2 = momentum along y for the creation operator on B site with spin down
  // kx3 = momentum along x for the annihilation operator on B site with spin up
  // ky3 = momentum along y for the annihilation operator on B site with spin up
  // kx4 = momentum along x for the annihilation operator on B site with spin down
  // ky4 = momentum along y for the annihilation operator on B site with spin down
  // return value = corresponding matrix element
  virtual Complex ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4);
    
  // compute the matrix element for the two body interaction between two sites A and B  belonging to the same layer
  //
  // kx1 = momentum along x for the A site
  // ky1 = momentum along y for the A site
  // kx2 = momentum along x for the B site
  // ky2 = momentum along y for the B site
  // return value = corresponding matrix element
  //  Complex ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the two body interaction between two sites A with different layer indices 
  //
  // kx1 = momentum along x for the first A site
  // ky1 = momentum along y for the first A site
  // kx2 = momentum along x for the second A site
  // ky2 = momentum along y for the second A site
  // return value = corresponding matrix element
  //  Complex ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kx2, int ky2);

  // compute the matrix element for the two body interaction between two sites B with different layer indices 
  //
  // kx1 = momentum along x for the first B site
  // ky1 = momentum along y for the first B site
  // kx2 = momentum along x for the second B site
  // ky2 = momentum along y for the second B site
  // return value = corresponding matrix element
  //  Complex ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kx2, int ky2);
    

  // compute the one body transformation matrices and the optional one body band stucture contribution
  //
  // oneBodyBasis = array of one body transformation matrices
  virtual void ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis);


};


#endif
