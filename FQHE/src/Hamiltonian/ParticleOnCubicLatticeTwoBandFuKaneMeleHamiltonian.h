////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of 3d topological insulator based on the Fu-Kane-Mele     //
//                       model and restricted to two bands                    //
//                                                                            //
//                        last modification : 18/07/2011                      //
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


#ifndef PARTICLEONCUBICLATTICETWOBANDFUKANEMELECHECKERBOARDHAMILTONIAN_H
#define PARTICLEONCUBICLATTICETWOBANDFUKANEMELECHECKERBOARDHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian : public ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian
{

 protected:
 
  // number of sites in the z direction
  int NbrSiteZ;
  // number of sites in the direction perpendicular to X
  int NbrSiteYZ;
  
  // hoping amplitude between neareast neighbor sites
  double NNHoping;
  // hoping amplitude between next neareast neighbor sites
  double NextNNHoping;
  // hoping amplitude between second next neareast neighbor sites
  double SecondNextNNHoping;
  // boundary condition twisting angle along x
  double GammaX;
  // boundary condition twisting angle along y
  double GammaY;
  // nearest neighbor density-density potential strength
  double UPotential;
  // mixing term coupling the two copies of the checkerboard lattice
  Complex MixingTerm;

  // use flat band model
  bool FlatBand;

 public:

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // uPotential = strength of the repulsive two body neareast neighbor interaction
  // t1 = hoping amplitude between neareast neighbor sites
  // t2 = hoping amplitude between next neareast neighbor sites
  // t2p = hoping amplitude between second next neareast neighbor sites
  // mixingTermNorm = norm of the mixing term coupling the two copies of the checkerboard lattice
  // mixingTermArgv = argument of the mixing term coupling the two copies of the checkerboard lattice
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, int nbrSiteZ, double uPotential, double t1, double t2, double t2p, double mixingTermNorm, double mixingTermArg, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);

  // destructor
  //
  ~ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();

  // compute the matrix element for the two body interaction between two sites A and B  belonging to the same layer
  //
  // kx1 = momentum along x for the A site
  // ky1 = momentum along y for the A site
  // kz1 = momentum along z for the A site
  // kx2 = momentum along x for the B site
  // ky2 = momentum along y for the B site
  // kz2 = momentum along z for the B site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2);

  // compute the matrix element for the two body interaction between two sites A with different layer indices 
  //
  // kx1 = momentum along x for the first A site
  // ky1 = momentum along y for the first A site
  // kz1 = momentum along z for the first A site
  // kx2 = momentum along x for the second A site
  // ky2 = momentum along y for the second A site
  // kz2 = momentum along z for the second A site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2);

  // compute the matrix element for the two body interaction between two sites B with different layer indices 
  //
  // kx1 = momentum along x for the first B site
  // ky1 = momentum along y for the first B site
  // kz1 = momentum along z for the first B site
  // kx2 = momentum along x for the second B site
  // ky2 = momentum along y for the second B site
  // kz2 = momentum along z for the second B site
  // return value = corresponding matrix element
  Complex ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2);
    

  // compute the transformation basis contribution to the interaction matrix element
  // 
  // oneBodyBasis = array of transformation basis matrices
  // momentumIndex1 = compact momentum index of the first creation operator
  // momentumIndex2 = compact momentum index of the second creation operator
  // momentumIndex3 = compact momentum index of the first annihilation operator
  // momentumIndex4 = compact momentum index of the second annihiliation operator
  // energyIndex1 = energy index of the first creation operator
  // energyIndex2 = energy index of the second creation operator
  // energyIndex3 = energy index of the first annihilation operator
  // energyIndex4 = energy index of the second annihiliation operator
  // siteIndex1 = site index of the first creation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
  // siteIndex2 = site index of the second creation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
  // siteIndex3 = site index of the first annihilation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
  // siteIndex4 = site index of the second annihiliation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
  Complex ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
						int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
						int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
						int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4);

};

// compute the transformation basis contribution to the interaction matrix element
// 
// oneBodyBasis = array of transformation basis matrices
// momentumIndex1 = compact momentum index of the first creation operator
// momentumIndex2 = compact momentum index of the second creation operator
// momentumIndex3 = compact momentum index of the first annihilation operator
// momentumIndex4 = compact momentum index of the second annihiliation operator
// energyIndex1 = energy index of the first creation operator
// energyIndex2 = energy index of the second creation operator
// energyIndex3 = energy index of the first annihilation operator
// energyIndex4 = energy index of the second annihiliation operator
// siteIndex1 = site index of the first creation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
// siteIndex2 = site index of the second creation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
// siteIndex3 = site index of the first annihilation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)
// siteIndex4 = site index of the second annihiliation operator (0=Aup, 1=Bup, 2=Adown, 3=Bdown)

inline Complex ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ComputeTransfomationBasisContribution(ComplexMatrix* oneBodyBasis,
													 int momentumIndex1, int momentumIndex2, int momentumIndex3, int momentumIndex4, 
													 int energyIndex1, int energyIndex2, int energyIndex3, int energyIndex4,
													 int siteIndex1, int siteIndex2, int siteIndex3, int siteIndex4)
{
  return (Conj(oneBodyBasis[momentumIndex1][energyIndex1][siteIndex1]) * oneBodyBasis[momentumIndex3][energyIndex3][siteIndex3] * Conj(oneBodyBasis[momentumIndex2][energyIndex2][siteIndex2]) * oneBodyBasis[momentumIndex4][energyIndex4][siteIndex4]);
}


#endif
