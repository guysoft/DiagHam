 ////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//             class of dice lattice model with interacting particles         //
//                   in the single chern 2 band approximation                 // 
//                                                                            //
//                        last modification : 16/09/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian.h"
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
// uPotential = strength of the repulsive two body neareast neighbor interaction
// t1 = real part of the hoping amplitude between neareast neighbor sites
// t2 = real part of the hoping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hoping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hoping amplitude between next neareast neighbor sites
// mus = sublattice chemical potential on A sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, 
												       int nbrSiteY, double uPotential, double t1, double t2, double lambda1, double lambda2, double mus, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);

  this->HamiltonianShift = 0.0;

  this->T1Hoping = t1;
  this->T2Hoping = t2;
  this->Lambda1Hoping = lambda1;
  this->Lambda2Hoping = lambda2;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;

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

ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::~ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
  if (this->FlatBand == false)
    this->OneBodyInteractionFactors = new double [this->NbrSiteX * this->NbrSiteY];
  this->ComputeOneBodyMatrices(OneBodyBasis);
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;      
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = (kx1 * this->NbrSiteY) + ky1;
		int Index2 = (kx2 * this->NbrSiteY) + ky2;
		if (Index1 < Index2)
		  ++this->NbrSectorIndicesPerSum[(((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)];    
	      }
      this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  if (this->NbrSectorIndicesPerSum[i]  > 0)
	    {
	      this->SectorIndicesPerSum[i] = new int[2 * this->NbrSectorIndicesPerSum[i]];      
	      this->NbrSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = (kx1 * this->NbrSiteY) + ky1;
		int Index2 = (kx2 * this->NbrSiteY) + ky2;
		if (Index1 < Index2)
		  {
		    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }
      double FactorUAB = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUAC = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUBC = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));

      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteY;
		  int ky3 = Index3 % this->NbrSiteY;
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;

 		  this->InteractionFactors[i][Index] = FactorUAB * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteA(kx2, ky2, kx4, ky4);
 		  this->InteractionFactors[i][Index] -= FactorUAB * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteA(kx1, ky1, kx4, ky4);
 		  this->InteractionFactors[i][Index] -= FactorUAB * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteA(kx2, ky2, kx3, ky3);
 		  this->InteractionFactors[i][Index] += FactorUAB * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteA(kx1, ky1, kx3, ky3);

 		  this->InteractionFactors[i][Index] += FactorUAC * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteB(kx2, ky2, kx4, ky4);
 		  this->InteractionFactors[i][Index] -= FactorUAC * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteB(kx1, ky1, kx4, ky4);
 		  this->InteractionFactors[i][Index] -= FactorUAC * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteB(kx2, ky2, kx3, ky3);
 		  this->InteractionFactors[i][Index] += FactorUAC * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteB(kx1, ky1, kx3, ky3);

 		  this->InteractionFactors[i][Index] += FactorUBC * (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteC(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors[i][Index] -= FactorUBC * (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteC(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactors[i][Index] -= FactorUBC * (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteC(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactors[i][Index] += FactorUBC * (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteC(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

		  this->InteractionFactors[i][Index] *= -2.0;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the matrix element for the on-site two body interaction for site A 
//
// kx1 = creation momentum along x for the B site
// ky1 = creation momentum along y for the B site
// kx2 = annihilation momentum along x for the B site
// ky2 = annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteA(int kx1, int ky1, int kx2, int ky2)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KxFactor * ((double) (kx2 - kx1))));
  return Tmp;
}

// compute the matrix element for the on-site two body interaction for site B 
//
// kx1 = creation momentum along x for the C site
// ky1 = creation momentum along y for the C site
// kx2 = annihilation momentum along x for the C site
// ky2 = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteB(int kx1, int ky1, int kx2, int ky2)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KyFactor * ((double) (ky2 - ky1))));
  return Tmp;
}

// compute the matrix element for the on-site two body interaction for site C 
//
// kx1 = creation momentum along x for the B site
// ky1 = creation momentum along y for the B site
// kx2 = creation momentum along x for the C site
// ky2 = creation momentum along y for the C site
// kx3 = annihilation momentum along x for the B site
// ky3 = annihilation momentum along y for the B site
// kx4 = annihilation momentum along x for the C site
// ky4 = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 2.0 * cos (0.5 * ((this->KxFactor * ((double) (kx3 - kx1))) + (this->KyFactor * ((double) (ky4 - ky2)))));
  return Tmp;
}

// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = array of one body transformation matrices

void ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = (kx * this->NbrSiteY) + ky;
	Complex HAB (this->T1Hoping, -this->Lambda1Hoping);
	HAB *= cos (0.5 * this->KxFactor * ((double) kx));
	Complex HAC(this->T1Hoping, this->Lambda1Hoping);
	HAC *= cos (0.5 * this->KyFactor * ((double) ky));
	Complex HBC(this->T1Hoping, this->Lambda1Hoping);
	HBC *= cos (0.5 * (this->KxFactor * ((double) kx) - this->KyFactor * ((double) ky)));
	HermitianMatrix TmpOneBobyHamiltonian(3, true);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 1, HAB);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 2, HAC);
	TmpOneBobyHamiltonian.SetMatrixElement(1, 2, HBC);
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
	oneBodyBasis[Index] = TmpMatrix;	
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactors[Index] = TmpDiag(0, 0);
	  }
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2) << endl;
      }
}
