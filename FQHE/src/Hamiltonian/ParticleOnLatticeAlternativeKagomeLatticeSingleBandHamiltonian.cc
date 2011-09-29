////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Yang-Le Wu                            //
//                                                                            //
//            class of Kagome lattice model with interacting particles        //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 06/09/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian.h"
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
// uPotential = strength of the repulsive two body nearest neighbor interaction
// vPotential = strength of the repulsive two body second nearest neighbor interaction
// t1 = hoping amplitude between nearest neighbor sites
// t2 = hoping amplitude between next nearest neighbor sites
// l1 = Rashba coupling between nearest neighbor sites
// l2 = Rashba coupling between next nearest neighbor sites
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// bandIndex = index of the band that has to be partially filled
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential, double vPotential, double t1, double t2, double l1, double l2, double gammaX, double gammaY, int bandIndex, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->HamiltonianShift = 0.0;
  this->NNHoping = t1;
  this->NextNNHoping = t2;
  this->NNRashba = l1;
  this->NextNNRashba = l2;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  this->BandIndex = bandIndex;
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

ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::~ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
  if (this->FlatBand == false)
    this->OneBodyInteractionFactors = new double [this->NbrSiteX * this->NbrSiteY];
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
  {
    double x=2*M_PI*((double)kx + this->GammaX)/this->NbrSiteX;
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
        double y=2*M_PI*((double)ky + this->GammaY)/this->NbrSiteY;
	int Index = (kx * this->NbrSiteY) + ky;

        Complex nnBA = Complex(-this->NNHoping, -this->NNRashba) * (1 + Phase(x));
        Complex nnCA = Complex(-this->NNHoping, +this->NNRashba) * (1 + Phase(y));
        Complex nnCB = Complex(-this->NNHoping, -this->NNRashba) * (1 + Phase(y-x));
        Complex nnnBA = Complex(-this->NextNNHoping, +this->NextNNRashba) * (Phase(y) + Phase(x-y));
        Complex nnnCA = Complex(-this->NextNNHoping, -this->NextNNRashba) * (Phase(x) + Phase(y-x));
        Complex nnnCB = Complex(-this->NextNNHoping, +this->NextNNRashba) * (Phase(-x) + Phase(y));

	HermitianMatrix TmpOneBodyHamiltonian(3, true);
        TmpOneBodyHamiltonian.SetMatrixElement(1, 0, nnBA + nnnBA);
        TmpOneBodyHamiltonian.SetMatrixElement(2, 0, nnCA + nnnCA);
        TmpOneBodyHamiltonian.SetMatrixElement(2, 1, nnCB + nnnCB);
	ComplexMatrix TmpMatrix(3, 3, true);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	TmpMatrix[2][2] = 1.0;
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	OneBodyBasis[Index] = TmpMatrix;	
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactors[Index] = TmpDiag(0, 0);
	  }
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2) << "  e1=[" << TmpMatrix[0][0] << ", " << TmpMatrix[0][1] << ", " << TmpMatrix[0][2] << "]  e2=[" << TmpMatrix[1][0] << ", " << TmpMatrix[1][1] << ", " << TmpMatrix[1][2] << "]  e3=[" << TmpMatrix[2][0] << ", " << TmpMatrix[2][1] << ", " << TmpMatrix[2][2] << "]" << endl;
      }
  }
 
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
      double FactorU = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      if (this->FlatBand == false)
	FactorU *= this->UPotential;
      double FactorV = this->VPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1) // annihilation operators
	    {
	      int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2) // creation operators
		{
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteY;
		  int ky3 = Index3 % this->NbrSiteY;
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;
                  // the InteractionFactors is supposed to be the coefficients to   A+_3 A_1 A+_4 A_2
                  // tricky part: OneBodyBasis[Index] stores the result of LapackDiagonalize
                  // and its [0][_] elements are the COMPLEX CONJUGATE of wave functions < _ |lower band>. (See the end of HermitianMatrix.cc)

                  Complex sumU = 0.;
                  Complex sumV = 0.;

                  sumU += Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * this->ComputeTwoBodyMatrixElementNNAB(kx2, ky2, kx4, ky4);
                  sumU += Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeTwoBodyMatrixElementNNBC(kx2, ky2, kx4, ky4);
                  sumU += Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * this->ComputeTwoBodyMatrixElementNNCA(kx2, ky2, kx4, ky4);
                  sumU -= Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
                        * this->ComputeTwoBodyMatrixElementNNAB(kx2, ky2, kx3, ky3);
                  sumU -= Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
                        * this->ComputeTwoBodyMatrixElementNNBC(kx2, ky2, kx3, ky3);
                  sumU -= Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
                        * this->ComputeTwoBodyMatrixElementNNCA(kx2, ky2, kx3, ky3);
                  sumU -= Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * this->ComputeTwoBodyMatrixElementNNAB(kx1, ky1, kx4, ky4);
                  sumU -= Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeTwoBodyMatrixElementNNBC(kx1, ky1, kx4, ky4);
                  sumU -= Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * this->ComputeTwoBodyMatrixElementNNCA(kx1, ky1, kx4, ky4);
                  sumU += Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
                        * this->ComputeTwoBodyMatrixElementNNAB(kx1, ky1, kx3, ky3);
                  sumU += Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
                        * this->ComputeTwoBodyMatrixElementNNBC(kx1, ky1, kx3, ky3);
                  sumU += Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
                        * this->ComputeTwoBodyMatrixElementNNCA(kx1, ky1, kx3, ky3);





                  sumV += Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * this->ComputeTwoBodyMatrixElementNNNAB(kx2, ky2, kx4, ky4);
                  sumV += Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeTwoBodyMatrixElementNNNBC(kx2, ky2, kx4, ky4);
                  sumV += Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * this->ComputeTwoBodyMatrixElementNNNCA(kx2, ky2, kx4, ky4);
                  sumV -= Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
                        * this->ComputeTwoBodyMatrixElementNNNAB(kx2, ky2, kx3, ky3);
                  sumV -= Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
                        * this->ComputeTwoBodyMatrixElementNNNBC(kx2, ky2, kx3, ky3);
                  sumV -= Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
                        * this->ComputeTwoBodyMatrixElementNNNCA(kx2, ky2, kx3, ky3);
                  sumV -= Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * this->ComputeTwoBodyMatrixElementNNNAB(kx1, ky1, kx4, ky4);
                  sumV -= Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeTwoBodyMatrixElementNNNBC(kx1, ky1, kx4, ky4);
                  sumV -= Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * this->ComputeTwoBodyMatrixElementNNNCA(kx1, ky1, kx4, ky4);
                  sumV += Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index3][this->BandIndex][1]
                        * this->ComputeTwoBodyMatrixElementNNNAB(kx1, ky1, kx3, ky3);
                  sumV += Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index3][this->BandIndex][2]
                        * this->ComputeTwoBodyMatrixElementNNNBC(kx1, ky1, kx3, ky3);
                  sumV += Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index3][this->BandIndex][0]
                        * this->ComputeTwoBodyMatrixElementNNNCA(kx1, ky1, kx3, ky3);


		  this->InteractionFactors[i][Index] = -2.0 * (FactorU * sumU + FactorV * sumV);

		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the matrix element for the two body interaction between two sites A and B 
//
// kx2 = annihilation momentum along x for the B site
// ky2 = annihilation momentum along y for the B site
// kx4 = creation momentum along x for the B site
// ky4 = creation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementNNAB(int kx2, int ky2, int kx4, int ky4)
{
  double dx = 2.0 * M_PI * ((double)(kx2-kx4)) / this->NbrSiteX;
  double dy = 2.0 * M_PI * ((double)(ky2-ky4)) / this->NbrSiteY;
  Complex Tmp = 1 + Phase(-dx);
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B and C
//
// kx2 = annihilation momentum along x for the C site
// ky2 = annihilation momentum along y for the C site
// kx4 = creation momentum along x for the C site
// ky4 = creation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementNNBC(int kx2, int ky2, int kx4, int ky4)
{
  double dx = 2.0 * M_PI * ((double)(kx2-kx4)) / this->NbrSiteX;
  double dy = 2.0 * M_PI * ((double)(ky2-ky4)) / this->NbrSiteY;
  Complex Tmp = 1 + Phase(dx - dy);
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites C and A
//
// kx2 = annihilation momentum along x for the A site
// ky2 = annihilation momentum along y for the A site
// kx4 = creation momentum along x for the A site
// ky4 = creation momentum along y for the A site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementNNCA(int kx2, int ky2, int kx4, int ky4)
{
  double dx = 2.0 * M_PI * ((double)(kx2-kx4)) / this->NbrSiteX;
  double dy = 2.0 * M_PI * ((double)(ky2-ky4)) / this->NbrSiteY;
  Complex Tmp = 1 + Phase(dy);
  return Tmp;
}

// compute the matrix element for the two body interaction between two A sites
//
// kx2 = annihilation momentum along x for the second site
// ky2 = annihilation momentum along y for the second site
// kx4 = creation momentum along x for the second site
// ky4 = creation momentum along y for the second site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementNNNAB(int kx2, int ky2, int kx4, int ky4)
{
  double dx = 2.0 * M_PI * ((double)(kx2-kx4)) / this->NbrSiteX;
  double dy = 2.0 * M_PI * ((double)(ky2-ky4)) / this->NbrSiteY;
  Complex Tmp = Phase(-dy) + Phase(dy - dx);
  return Tmp;
}

// compute the matrix element for the two body interaction between two B sites
//
// kx2 = annihilation momentum along x for the second site
// ky2 = annihilation momentum along y for the second site
// kx4 = creation momentum along x for the second site
// ky4 = creation momentum along y for the second site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementNNNBC(int kx2, int ky2, int kx4, int ky4)
{
  double dx = 2.0 * M_PI * ((double)(kx2-kx4)) / this->NbrSiteX;
  double dy = 2.0 * M_PI * ((double)(ky2-ky4)) / this->NbrSiteY;
  Complex Tmp = Phase(dx) + Phase(-dy);
  return Tmp;
}

// compute the matrix element for the two body interaction between two C sites
//
// kx2 = annihilation momentum along x for the second site
// ky2 = annihilation momentum along y for the second site
// kx4 = creation momentum along x for the second site
// ky4 = creation momentum along y for the second site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementNNNCA(int kx2, int ky2, int kx4, int ky4)
{
  double dx = 2.0 * M_PI * ((double)(kx2-kx4)) / this->NbrSiteX;
  double dy = 2.0 * M_PI * ((double)(ky2-ky4)) / this->NbrSiteY;
  Complex Tmp = Phase(dx) + Phase(dy - dx);
  return Tmp;
}
