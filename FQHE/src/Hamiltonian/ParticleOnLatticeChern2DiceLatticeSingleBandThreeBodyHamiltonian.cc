//////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of ruby lattice model with interacting particles          //
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 25/10/2011                      //
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


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include "GeneralTools/StringTools.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;



// default constructor
//

ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// threeBodyPotential = strength of the repulsive three body neareast neighbor interaction
// uPotential = strength of the repulsive two body neareast neighbor interaction
// vPotential = strength of the repulsive two body neareast neighbor interaction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mus = sublattice chemical potential on A sites
// gammaX = boundary condition twisting angle along e_a (measured in units of 2pi)
// gammaY = boundary condition twisting angle along e_b (measured in units of 2pi)
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, double threeBodyPotential, double uPotential, double t, double epsilon, double lambda, double bfield1, double bfield3, double intraCoefficient, double interCoefficient, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrCellX;
  this->NbrSiteY = nbrCellY;
  this->LzMax = nbrCellX * nbrCellY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);

  this->HamiltonianShift = 0.0;
  this->TrHopping = 0.0;
  this->TiHopping = 0.0;
  this->T1rHopping = 0.0;
  this->T1iHopping = 0.0;
  this->T4Hopping = 0.0;

  this->THopping = t;
  this->Epsilon = epsilon;
  this->Lambda = lambda;
  this->BField1 = bfield1;
  this->BField3 = bfield3;
  
  this->MuS = 0.0;
  this->GammaX = gammaX;
  this->GammaY = gammaY;


  this->NBodyValue = 3;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->FlatBand = flatBandFlag;
  this->UPotential = threeBodyPotential;
  this->VPotential = uPotential;
  this->IntraCoefficient = intraCoefficient;
  this->InterCoefficient = interCoefficient;

  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;

  this->TwoBodyFlag = false;
  
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      cout << "fast = ";
      PrintMemorySize(cout, TmpMemory)<< endl;
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::~ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian()
{
}


// conventions adopted for matrix elements:
// modified with respect to Tang, Mei, Wen:
// unit cell is triangle standing on base. bottom left corner site A and bottom right corner B, tip is site C.


// compute the matrix element for the two body interaction between two sites A and B 
//
// k1a = creation momentum along x for the B site
// k1b = creation momentum along y for the B site
// k2a = annihilation momentum along x for the B site
// k2b = annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementAB(int k1a, int k1b, int k2a, int k2b)
{
  //Complex Tmp = (1+Phase(-this->KyFactor * ((double) (k2b - k1b))) + Phase(-this->KxFactor * ((double) (k2a - k1a))));
  return 0.0;
}

// compute the matrix element for the two body interaction between two sites A and C 
//
// k1a = creation momentum along x for the C site
// k1b = creation momentum along y for the C site
// k2a = annihilation momentum along x for the C site
// k2b = annihilation momentum along y for the C site
// return value = corresponding matrix element

/*Complex ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementAC(int k1a, int k1b, int k2a, int k2b)
{
 Complex Tmp = (1+Phase(-this->KyFactor * ((double) (k2b - k1b))) + Phase(-this->KxFactor * ((double) (k1a - k2a))+this->KyFactor * ((double) (k1b - k2b))));
 return Tmp;
 }*/

// compute the matrix element for the two body interaction between two sites B and C 
//
// k1a = creation momentum along x for the B site
// k1b = creation momentum along y for the B site
// k2a = creation momentum along x for the C site
// k2b = creation momentum along y for the C site
// k3a = annihilation momentum along x for the B site
// k3b = annihilation momentum along y for the B site
// k4a = annihilation momentum along x for the C site
// k4b = annihilation momentum along y for the C site
// return value = corresponding matrix element

/*
 Complex ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementBC(int k1a, int k1b, int k2a, int k2b, int k3a, int k3b, int k4a, int k4b)
{
  Complex Tmp = (1+Phase(-this->KxFactor * ((double) (k3a - k1a))) + Phase(-this->KxFactor * ((double) (k3a - k1a)) - this->KyFactor * ((double) (k1b - k3b))));
  return Tmp;
 }
*/


// compute the matrix element for on-site two body interaction involving A sites
//
// return value = corresponding matrix element
 
Complex ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementAA(int kx1, int ky1, int kx2, int ky2)
{
  return 0.0;
}

// compute the matrix element for on-site two body interaction involving B sites
//
// kx1 = first creation momentum along x for the B site
// ky1 = first creation momentum along y for the B site
// kx2 = second creation momentum along x for the B site
// ky2 = second creation momentum along y for the B site
// kx3 = first annihilation momentum along x for the B site
// ky3 = first annihilation momentum along y for the B site
// kx4 = second annihilation momentum along x for the B site
// ky4 = second annihilation momentum along y for the B site
// return value = corresponding matrix element
  /*
Complex ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementOnSiteBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return 0.0;
}
  */
// compute the matrix element for on-site two body interaction involving C sites
//
// kx1 = first creation momentum along x for the C site
// ky1 = first creation momentum along y for the C site
// kx2 = second creation momentum along x for the C site
// ky2 = second creation momentum along y for the C site
// kx3 = first annihilation momentum along x for the C site
// ky3 = first annihilation momentum along y for the C site
// kx4 = second annihilation momentum along x for the C site
// ky4 = second annihilation momentum along y for the C site
// return value = corresponding matrix element
   /*
Complex ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementOnSiteCC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return 0.0;
}

/* END TWO-BODY TERMS */



// compute the matrix element for the on-site three body interaction related to sites A
//
// kx1 = first creation momentum along x for the first A site
// ky1 = first creation momentum along y for the first A site
// kx2 = second creation momentum along x for the second A site
// ky2 = second creation momentum along y for the second A site
// kx3 = third creation momentum along x for the second A site
// ky3 = third creation momentum along y for the second A site
// kx4 = first annihilation momentum along x for the first A site
// ky4 = first annihilation momentum along y for the first A site
// kx5 = second annihilation momentum along x for the second A site
// ky5 = second annihilation momentum along y for the second A site
// kx6 = third annihilation momentum along x for the second A site
// ky6 = third annihilation momentum along y for the second A site
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the on-site three body interaction related to sites B
//
// kx1 = first creation momentum along x for the first B site
// ky1 = first creation momentum along y for the first B site
// kx2 = second creation momentum along x for the second B site
// ky2 = second creation momentum along y for the second B site
// kx3 = third creation momentum along x for the second B site
// ky3 = third creation momentum along y for the second B site
// kx4 = first annihilation momentum along x for the first B site
// ky4 = first annihilation momentum along y for the first B site
// kx5 = second annihilation momentum along x for the second B site
// ky5 = second annihilation momentum along y for the second B site
// kx6 = third annihilation momentum along x for the second B site
// ky6 = third annihilation momentum along y for the second B site
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteBBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the on-site three body interaction related to sites A
//
// kx1 = first creation momentum along x for the first C site
// ky1 = first creation momentum along y for the first C site
// kx2 = second creation momentum along x for the second C site
// ky2 = second creation momentum along y for the second C site
// kx3 = third creation momentum along x for the second C site
// ky3 = third creation momentum along y for the second C site
// kx4 = first annihilation momentum along x for the first C site
// ky4 = first annihilation momentum along y for the first C site
// kx5 = second annihilation momentum along x for the second C site
// ky5 = second annihilation momentum along y for the second C site
// kx6 = third annihilation momentum along x for the second C site
// ky6 = third annihilation momentum along y for the second C site
// return value = corresponding matrix element
/*
Complex ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteCCC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}
*/

void ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
  double KX;
  double KY;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  KX = this->KxFactor * (((double) kx) + this->GammaX);
	  KY = this->KyFactor * (((double) ky) + this->GammaY);
	  int Index = (kx * this->NbrSiteY) + ky;

          double angle = 2.0 * M_PI / 3.0;
          Complex GammaK = this->THopping * (1.0 + Phase(-KX) + Phase(-KY));
          Complex GammaP = this->Lambda * (1.0 + Phase(-KX + angle) + Phase(-KY - angle));
          Complex GammaM = this->Lambda * (1.0 + Phase(-KX - angle) + Phase(-KY + angle));

	  HermitianMatrix TmpOneBodyHamiltonian(6, true);
 	  TmpOneBodyHamiltonian.SetMatrixElement(0, 4, Conj(GammaK));
 	  TmpOneBodyHamiltonian.SetMatrixElement(1, 5, Conj(GammaK));
 	  TmpOneBodyHamiltonian.SetMatrixElement(0, 5, I()*Conj(GammaP));
 	  TmpOneBodyHamiltonian.SetMatrixElement(1, 4, I()*Conj(GammaM));

 	  TmpOneBodyHamiltonian.SetMatrixElement(2, 4, GammaK);
 	  TmpOneBodyHamiltonian.SetMatrixElement(3, 5, GammaK);
 	  TmpOneBodyHamiltonian.SetMatrixElement(2, 5, I()*GammaM);
 	  TmpOneBodyHamiltonian.SetMatrixElement(3, 4, I()*GammaP);

	  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, this->BField1);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -this->BField1);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 2, this->BField1);
	  TmpOneBodyHamiltonian.SetMatrixElement(3, 3, -this->BField1);
	  TmpOneBodyHamiltonian.SetMatrixElement(4, 4, this->Epsilon + this->BField3);
	  TmpOneBodyHamiltonian.SetMatrixElement(5, 5, this->Epsilon - this->BField3);

	  TmpOneBodyHamiltonian *= -1.0;

	  ComplexMatrix TmpMatrix(6, 6, true);
	  TmpMatrix.SetToIdentity();
	  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
	  oneBodyBasis[Index] = TmpMatrix;
	  if (this->FlatBand == false)
	    {
	      this->OneBodyInteractionFactors[Index] = 0.5 * TmpDiag(0, 0);
	  }
	  cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2) << " " 
	       << TmpDiag(3, 3) << " " << TmpDiag(4, 4) << " " << TmpDiag(5, 5) << endl;
	}
    }
}

// compute the matrix element for the three body interaction between one site B and two sites A
//
// kx1 = creation momentum along x for the first A site
// ky1 = creation momentum along y for the first A site
// kx2 = creation momentum along x for the second A site
// ky2 = creation momentum along y for the second A site
// kx3 = annihilation momentum along x for the first A site
// ky3 = annihilation momentum along y for the first A site
// kx4 = annihilation momentum along x for the second A site
// ky4 = annihilation momentum along y for the second A site
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementBAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  Complex Tmp = 0.0;
  return Tmp;
}

// compute the matrix element for the on-site three body interaction with two particles on a B site and one on a A site
//
// kx1 = creation momentum along x for the A site
// ky1 = creation momentum along y for the A site
// kx2 = first creation momentum along x for the B site
// ky2 = first creation momentum along y for the B site
// kx3 = second creation momentum along x for the B site
// ky3 = second creation momentum along y for the B site
// kx4 = annihilation momentum along x for the A site
// ky4 = annihilation momentum along y for the A site
// kx5 = first annihilation momentum along x for the B site
// ky5 = first annihilation momentum along y for the B site
// kx6 = second annihilation momentum along x for the B site
// ky6 = second annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  Complex Tmp = 0.0;
  return Tmp;
}

// compute the matrix element for the on-site three body interaction with two particles on a A site and one on a B site
//
// kx1 = first creation momentum along x for the A site
// ky1 = first creation momentum along y for the A site
// kx2 = second creation momentum along x for the A site
// ky2 = second creation momentum along y for the A site
// kx3 = creation momentum along x for the B site
// ky3 = creation momentum along y for the B site
// kx4 = first annihilation momentum along x for the A site
// ky4 = first annihilation momentum along y for the A site
// kx5 = second annihilation momentum along x for the A site
// ky5 = second annihilation momentum along y for the sA site
// kx6 = annihilation momentum along x for the B site
// ky6 = annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementOnSiteAAB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  Complex Tmp = 0.0;
  return Tmp;
}

Complex ParticleOnLatticeChern2DiceLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 0.0;
}
