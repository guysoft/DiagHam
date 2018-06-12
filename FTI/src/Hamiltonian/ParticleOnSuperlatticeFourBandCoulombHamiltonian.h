////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//            class of Moire superlattice with Coulomb interactions           //
//                         projected to four bands                            //
//                                                                            //
//                        last modification : 17/05/2018                      //
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


#ifndef PARTICLEONSUPERLATTICEFOURBANDCOULOMBHAMILTONIAN_H
#define PARTICLEONSUPERLATTICEFOURBANDCOULOMBHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFourBandHamiltonian.h"
#include "Matrix/ComplexMatrix.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


class ParticleOnSuperlatticeFourBandCoulombHamiltonian : public ParticleOnLatticeQuantumSpinHallFourBandHamiltonian
{

 protected:
 
  // numerical factor for momentum along x
  double KxFactor;
  // numerical factor for momentum along y
  double KyFactor;
  
  // boundary condition twisting angle along x
  double GammaX;
  // boundary condition twisting angle along y
  double GammaY;
  // nearest neighbor density-density potential strength
  double UPotential;
  double ScreeningDistance;
  // use flat band model
  bool FlatBand;
  // index of the band in the tight-binding model
  int BandIndex;
  
  // Tight Binding Model
  Abstract2DTightBindingModel* TightBindingModel;
  
  //coordinates of the reciprocal lattice vectors
  double Gx1;
  double Gy1;
  double Gx2;
  double Gy2;
  
  int TruncateQx;
  int TruncateQy;

 public:

  // default constructor
  //
  ParticleOnSuperlatticeFourBandCoulombHamiltonian();

  // constructor
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // uPotential = repulsive on-site potential strength between different orbitals
  // vPotential = repulsive on-site potential strength between opposite spins
  // mass = mass term of the simple TI model
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // flatBandFlag = use flat band model
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  ParticleOnSuperlatticeFourBandCoulombHamiltonian (ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, int truncateQx, int truncateQy, double uPotential, double screeningDistance, double gammaX, double gammaY,  Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory = -1);
  
  // destructor
  //
  ~ParticleOnSuperlatticeFourBandCoulombHamiltonian();
  

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();  
  
  // compute the amplitude of the Coulomb interaction
  //
  // Index1 = momentum of the first creation operator
  // Index2 = momentum of the second creation operator
  // Index3 = momentum of the first annihilation operator
  // Index4 = momentum of the second annihilation operator
  // valleyIndex1 = valley index of the first density operator
  // valleyIndex2 = valley index of the second density operator
  // return value = complex amplitude
  virtual Complex EvaluateInteractionCoefficient(int Index1, int Index2, int Index3, int Index4, int valleyIndex1, int valleyIndex2);
    
  // compute the form factor for the density operator 
  // 
  // kx = momentum along x of annihilation operator
  // ky = momentum along y of creation operator
  // qx = momentum transfer along x direction
  // qy = momentum transfer along y direction
  // valleyIndex = valley index of density operator
  virtual Complex ComputeDensityFormFactor(int kx, int ky, int qx, int qy, int valleyIndex);
  
  // evaluate V(q) for the Coulomb interaction
  //
  // kx = component of momentum along first Bravais vector
  // ky = component of momentum along second Bravais vector
  // return value = amplitude of V(q)
  virtual double GetVofQ (int kx, int ky);

};


#endif
