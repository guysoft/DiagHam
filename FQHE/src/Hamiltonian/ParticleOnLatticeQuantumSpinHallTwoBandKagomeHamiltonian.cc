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


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = strength of the repulsive two body neareast neighbor interaction with identical spin
// vPotential = strength of the repulsive on site two body interaction with opposite spin
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

ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian::ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
														   int nbrSiteY, double uPotential, double vPotential, double wPotential, double t1, double t2, double lambda1, double lambda2, double mixingTermNorm, double mixingTermArg, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->HamiltonianShift = 0.0;
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->NNSpinOrbit = lambda1;
  this->NextNNSpinOrbit = lambda2;
  this->MixingTerm = mixingTermNorm * Phase(mixingTermArg);
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->WPotential = wPotential;

  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      if (TmpMemory < 1024)
	cout  << "fast = " <<  TmpMemory << "b ";
      else
	if (TmpMemory < (1 << 20))
	  cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	else
	  if (TmpMemory < (1 << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    {
	      cout  << "fast = " << (TmpMemory >> 30) << ".";
	      TmpMemory -= ((TmpMemory >> 30) << 30);
	      TmpMemory *= 100l;
	      TmpMemory >>= 30;
	      if (TmpMemory < 10l)
		cout << "0";
	      cout  << TmpMemory << " Gb ";
	    }
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian::~ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian()
{
}
  
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

Complex ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian::ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return (2.0 * (cos (0.5 * ((this->KxFactor * ((double) (kx4 - kx2))) - (this->KyFactor * ((double) (ky4 - ky2)))))
		 + cos (0.5 * ((this->KxFactor * ((double) (kx4 - kx2))) + (this->KyFactor * ((double) (ky4 - ky2)))))));
}

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

Complex ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian::ComputeTwoBodyMatrixElementADownBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return (2.0 * (cos (0.5 * ((this->KxFactor * ((double) (kx4 - kx2))) - (this->KyFactor * ((double) (ky4 - ky2)))))
		 + cos (0.5 * ((this->KxFactor * ((double) (kx4 - kx2))) + (this->KyFactor * ((double) (ky4 - ky2)))))));
}
  
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

Complex ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian::ComputeTwoBodyMatrixElementADownBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return (2.0 * (cos (0.5 * ((this->KxFactor * ((double) (kx3 - kx1))) - (this->KyFactor * ((double) (ky3 - ky1)))))
		 + cos (0.5 * ((this->KxFactor * ((double) (kx3 - kx1))) + (this->KyFactor * ((double) (ky3 - ky1)))))));
}

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

Complex ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian::ComputeTwoBodyMatrixElementAUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return (2.0 * (cos (0.5 * ((this->KxFactor * ((double) (kx4 - kx2))) - (this->KyFactor * ((double) (ky4 - ky2)))))
		 + cos (0.5 * ((this->KxFactor * ((double) (kx4 - kx2))) + (this->KyFactor * ((double) (ky4 - ky2)))))));
}

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

Complex ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian::ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 1.0;
  return Tmp;
}

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

Complex ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian::ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return  Phase (0.5 * ((this->KxFactor * ((double) (kx4 + kx3 - kx2 - kx1))) +
			(this->KyFactor * ((double) (ky4 + ky3 - ky2 - ky1)))));
}


// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = array of one body transformation matrices

void ParticleOnLatticeQuantumSpinHallTwoBandKagomeHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
  double KX, KY;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = (kx * this->NbrSiteY) + ky;
	oneBodyBasis[Index] = ComplexMatrix(6, 6, true);

	KX = 0.5 * this->KxFactor * (((double) kx) + this->GammaX);
	KY = 0.5 * this->KyFactor * (((double) ky) + this->GammaY);
	Complex HAB (-2.0 * this->NNHopping, -2.0 * this->NNSpinOrbit);
	HAB *= cos (KX);
	Complex HAC(-2.0 * this->NNHopping, 2.0 * this->NNSpinOrbit);
	HAC *= cos (KY);
	Complex HBC(-2.0 * this->NNHopping, -2.0 * this->NNSpinOrbit);
	HBC *= cos(KX - KY);

	Complex HAB2 (-2.0 * this->NextNNHopping, 2.0 * this->NextNNSpinOrbit);
	HAB2 *= cos (KX - 2.0 * KY);
	Complex HAC2 (-2.0 * this->NextNNHopping, -2.0 * this->NextNNSpinOrbit);
	HAC2 *= cos (2.0 * KX - KY);
	Complex HBC2 (-2.0 * this->NextNNHopping, 2.0  *  this->NextNNSpinOrbit);
	HBC2 *= cos (KX + KY);

	HAB += HAB2;
	HAC += HAC2;
	HBC += HBC2;
		
	HermitianMatrix TmpOneBodyHamiltonian(3, true);
	
	TmpOneBodyHamiltonian.SetMatrixElement(0, 1, HAB);
	TmpOneBodyHamiltonian.SetMatrixElement(0, 2, HAC);
	TmpOneBodyHamiltonian.SetMatrixElement(1, 2, HBC);
	ComplexMatrix TmpMatrix(3, 3, true);
	TmpMatrix.SetToIdentity();
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
	for (int i = 0; i < 3; ++i)
	  for (int j = 0; j < 3; ++j)
	    oneBodyBasis[Index][2 * i][j] = TmpMatrix[i][j];
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactorsupup[Index] = 0.5 * TmpDiag(0, 0);
	  }
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2) << endl;
		
	KX = 0.5 * this->KxFactor * (((double) -kx) + this->GammaX);
	KY = 0.5 * this->KyFactor * (((double) -ky) + this->GammaY);
	HAB = Complex(-2.0 * this->NNHopping, -2.0 * this->NNSpinOrbit);
	HAB *= cos (KX);
	HAC = Complex(-2.0 * this->NNHopping, 2.0 * this->NNSpinOrbit);
	HAC *= cos (KY);
	HBC = Complex(-2.0 * this->NNHopping, -2.0 * this->NNSpinOrbit);
	HBC *= cos(KX - KY);

	HAB2 = Complex(-2.0 * this->NextNNHopping, 2.0 * this->NextNNSpinOrbit);
	HAB2 *= cos (KX - 2.0 * KY);
	HAC2 = Complex(-2.0 * this->NextNNHopping, -2.0 * this->NextNNSpinOrbit);
	HAC2 *= cos (2.0 * KX - KY);
	HBC2 = Complex(-2.0 * this->NextNNHopping, 2.0 * this->NextNNSpinOrbit);
	HBC2 *= cos (KX + KY);

	HAB += HAB2;
	HAC += HAC2;
	HBC += HBC2;
		
	TmpOneBodyHamiltonian.ClearMatrix();	
	TmpOneBodyHamiltonian.SetMatrixElement(0, 1, Conj(HAB));
	TmpOneBodyHamiltonian.SetMatrixElement(0, 2, Conj(HAC));
	TmpOneBodyHamiltonian.SetMatrixElement(1, 2, Conj(HBC));
	TmpMatrix.SetToIdentity();
#ifdef __LAPACK__
	TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
	for (int i = 0; i < 3; ++i)
	  for (int j = 0; j < 3; ++j)
	    oneBodyBasis[Index][2 * i + 1][3 + j] = TmpMatrix[i][j];
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactorsdowndown[Index] = 0.5 * TmpDiag(0, 0);
	  }
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2) << endl;
      }
}
