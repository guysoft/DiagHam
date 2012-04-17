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
//                        last modification : 08/09/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeChern2TriangularLatticeSingleBandHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <cmath>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;
using std::sin;
using std::cos;


// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrCellX = number of sites in the x direction
// nbrCellY = number of sites in the y direction
// uPotential = strength of the repulsive two body neareast neighbor interaction
// vPotential = strength of the repulsive two body second nearest neighbor interaction
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

ParticleOnLatticeChern2TriangularLatticeSingleBandHamiltonian::ParticleOnLatticeChern2TriangularLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, double uPotential, double vPotential, double t1, double t2, double phi, double mus, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrCellX;
  this->NbrSiteY = nbrCellY;
  this->LzMax = nbrCellX * nbrCellY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->Phi = 2.0 * M_PI*phi;
  this->HamiltonianShift = 0.0;
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->NNSpinOrbit = 0;
  this->NextNNSpinOrbit = 0;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;

  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
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

ParticleOnLatticeChern2TriangularLatticeSingleBandHamiltonian::~ParticleOnLatticeChern2TriangularLatticeSingleBandHamiltonian()
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

Complex ParticleOnLatticeChern2TriangularLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementAB(int k1a, int k1b, int k2a, int k2b)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KxFactor * ((double) (k2a - k1a))));
  //Complex Tmp = Phase (0.5 * (this->KxFactor * ((double) (k2a - k1a))));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and C 
//
// k1a = creation momentum along x for the C site
// k1b = creation momentum along y for the C site
// k2a = annihilation momentum along x for the C site
// k2b = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2TriangularLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementAC(int k1a, int k1b, int k2a, int k2b)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KyFactor * ((double) (k2b - k1b))));
  //Complex Tmp = Phase (0.5 * (this->KyFactor * ((double) (k2b - k1b))));
  return Tmp;
}

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

Complex ParticleOnLatticeChern2TriangularLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementBC(int k1a, int k1b, int k2a, int k2b, int k3a, int k3b, int k4a, int k4b)
{
  Complex Tmp = 2.0 * cos (0.5 * ((this->KxFactor * ((double) (k3a - k1a))) + (this->KyFactor * ((double) (k4b - k2b)))));
  //Complex Tmp = Phase(0.5 * ((this->KxFactor * ((double) (k3a - k1a))) + (this->KyFactor * ((double) (k4b - k2b)))));
  return Tmp;
}


// compute the matrix element for on-site two body interaction involving A sites
//
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2TriangularLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteAA()
{
  return 1.0;
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

Complex ParticleOnLatticeChern2TriangularLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return (0.0 * Phase(-1.0*(((this->KxFactor * ((double) (kx1 + kx2 - kx3 -kx4))) +(this->KyFactor * ((double) (ky1 + ky2 - ky3 -ky4))))/3.0)));
}

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

Complex ParticleOnLatticeChern2TriangularLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteCC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return (0.0 * Phase(((this->KxFactor * ((double) (kx1 + kx2 - kx3 -kx4))) +(-2.0*this->KyFactor * ((double) (ky1 + ky2 - ky3 -ky4))))/3.0));
}

// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = array of one body transformation matrices

void ParticleOnLatticeChern2TriangularLatticeSingleBandHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
 double KA, KB;
  for (int ka = 0; ka < this->NbrSiteX; ++ka)
    for (int kb = 0; kb < this->NbrSiteY; ++kb)
      {
	KA = this->KxFactor * (((double) ka) + this->GammaX);
	KB = this->KyFactor * (((double) kb) + this->GammaY);
	int Index = (ka * this->NbrSiteY) + kb;
	
	Complex HAB = this->NNHopping * (1.0 + Phase(KA) + Phase(KB));
	Complex HAC = -1.0*this->NNHopping*(Phase(-2.0*this->Phi)+Phase(KB)+Phase(2.0*this->Phi+KB-KA));
	Complex HBC = this->NNHopping*(Phase(2.0*this->Phi) + Phase(KB - KA - 2.0*this->Phi) + Phase(-KA));

	Complex HAA = 2.0*this->NextNNHopping*(cos(KB+this->Phi)+cos(KA-this->Phi)+cos(KA - KB + this->Phi));
	Complex HBB = 2.0*this->NextNNHopping*(cos(KB-this->Phi)+cos(KA+this->Phi)+cos(KA - KB - this->Phi));
	Complex HCC = -2.0*this->NextNNHopping*(cos(KB)+cos(KA)+cos(KA - KB));
	
		
	HermitianMatrix TmpOneBobyHamiltonian(3, true);
	
	TmpOneBobyHamiltonian.SetMatrixElement(0, 1, HAB);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 2, HAC);
	TmpOneBobyHamiltonian.SetMatrixElement(1, 2, HBC);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 0, HAA);
	TmpOneBobyHamiltonian.SetMatrixElement(1, 1, HBB);
	TmpOneBobyHamiltonian.SetMatrixElement(2, 2, HCC);
	//cout << TmpOneBobyHamiltonian << endl;
	ComplexMatrix TmpMatrix(3, 3, true);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	TmpMatrix[2][2] = 1.0;
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBobyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBobyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
	//cout << TmpMatrix[0][0]<<" "<< TmpMatrix[0][1]<< " " << TmpMatrix[0][2]<<endl;
	oneBodyBasis[Index] = TmpMatrix;
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactors[Index] = 0.5*TmpDiag(0, 0);
	  }
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2) << endl;
      }
}
