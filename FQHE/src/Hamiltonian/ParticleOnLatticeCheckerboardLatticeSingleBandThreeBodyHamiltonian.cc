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
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 13/07/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian.h"
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
// vPotential = strength of the repulsive three body neareast neighbor interaction
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// t2p = hoping amplitude between second next neareast neighbor sites
// mus = sublattice staggered chemical potential 
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, 
												       int nbrSiteY, double uPotential, double vPotential, double t1, double t2, double t2p, double mus, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->HamiltonianShift = 0.0;
  this->NNHoping = t1;
  this->NextNNHoping = t2;
  this->SecondNextNNHoping = t2p;
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
  this->TwoBodyFlag = false;
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

ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::~ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
  if (this->FlatBand == false)
    this->OneBodyInteractionFactors = new double [this->NbrSiteX * this->NbrSiteY];
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = (kx * this->NbrSiteY) + ky;
	Complex B1 = 4.0 * this->NNHoping * Complex (cos (1.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (1.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY)) * cos(M_PI * 0.25), 
						     sin (1.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * sin (1.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY)) * sin(M_PI * 0.25));
	double d1 = 4.0 * this->SecondNextNNHoping * cos (2.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (2.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY));
	double d3 =  this->MuS + (2.0 * this->NextNNHoping * (cos (2.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX))
							      - cos (2.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY))));
	HermitianMatrix TmpOneBobyHamiltonian(2, true);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 0, d1 + d3);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 1, B1);
	TmpOneBobyHamiltonian.SetMatrixElement(1, 1, d1 - d3);
	ComplexMatrix TmpMatrix(2, 2, true);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBobyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBobyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	OneBodyBasis[Index] = TmpMatrix;	
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactors[Index] = TmpDiag(0, 0);
	  }
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << "  e1=[" << TmpMatrix[0][0] << ", " << TmpMatrix[0][1] << "]  e2=[" << TmpMatrix[1][0] << ", " << TmpMatrix[1][1] << "]" << endl;
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

      if (this->TwoBodyFlag == false)
	{
	  for (int i = 0; i < this->NbrSectorSums; ++i)	  
	    this->NbrSectorIndicesPerSum[i] = 0;
	  this->InteractionFactors = new Complex* [this->NbrSectorSums];
	}

      this->NbrThreeBodySectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrThreeBodySectorIndicesPerSum = new int[this->NbrThreeBodySectorSums];
      for (int i = 0; i < this->NbrThreeBodySectorSums; ++i)
	this->NbrThreeBodySectorIndicesPerSum[i] = 0;      
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
	    for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
		  {
		    int Index1 = (kx1 * this->NbrSiteY) + ky1;
		    int Index2 = (kx2 * this->NbrSiteY) + ky2;
		    int Index3 = (kx3 * this->NbrSiteY) + ky3;
		    if ((Index1 < Index2) && (Index2 < Index3))
		      ++this->NbrThreeBodySectorIndicesPerSum[(((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY)];    
		  }
      this->ThreeBodySectorIndicesPerSum = new int* [this->NbrThreeBodySectorSums];
      for (int i = 0; i < this->NbrThreeBodySectorSums; ++i)
	{
	  if (this->NbrThreeBodySectorIndicesPerSum[i]  > 0)
	    {
	      this->ThreeBodySectorIndicesPerSum[i] = new int[3 * this->NbrThreeBodySectorIndicesPerSum[i]];      
	      this->NbrThreeBodySectorIndicesPerSum[i] = 0;
	    }
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
	    for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	      for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
		for (int ky3 = 0; ky3 < this->NbrSiteY; ++ky3) 
		  {
		    int Index1 = (kx1 * this->NbrSiteY) + ky1;
		    int Index2 = (kx2 * this->NbrSiteY) + ky2;
		    int Index3 = (kx3 * this->NbrSiteY) + ky3;
		    if ((Index1 < Index2) && (Index2 < Index3))
		      {
			int TmpSum = (((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY);
			this->ThreeBodySectorIndicesPerSum[TmpSum][this->NbrThreeBodySectorIndicesPerSum[TmpSum] * 3] = Index1;
			this->ThreeBodySectorIndicesPerSum[TmpSum][1 + (this->NbrThreeBodySectorIndicesPerSum[TmpSum] * 3)] = Index2;
			this->ThreeBodySectorIndicesPerSum[TmpSum][2 + (this->NbrThreeBodySectorIndicesPerSum[TmpSum] * 3)] = Index3;
			++this->NbrThreeBodySectorIndicesPerSum[TmpSum];    
		      }
		  }


      double FactorU = this->UPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      double FactorV = this->VPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      this->ThreeBodyInteractionFactors = new Complex* [this->NbrThreeBodySectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->ThreeBodyInteractionFactors[i] = new Complex[this->NbrThreeBodySectorIndicesPerSum[i] * this->NbrThreeBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrThreeBodySectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->ThreeBodySectorIndicesPerSum[i][j1 * 3];
	      int Index2 = this->ThreeBodySectorIndicesPerSum[i][(j1 * 3) + 1];
	      int Index3 = this->ThreeBodySectorIndicesPerSum[i][(j1 * 3) + 2];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrThreeBodySectorIndicesPerSum[i]; ++j2)
		{
		  int Index4 = this->ThreeBodySectorIndicesPerSum[i][j2 * 3];
		  int Index5 = this->ThreeBodySectorIndicesPerSum[i][(j2 * 3) + 1];
		  int Index6 = this->ThreeBodySectorIndicesPerSum[i][(j2 * 3) + 2];
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;
		  int kx5 = Index5 / this->NbrSiteY;
		  int ky5 = Index5 % this->NbrSiteY;
		  int kx6 = Index6 / this->NbrSiteY;
		  int ky6 = Index6 % this->NbrSiteY;


                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5);
		  
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4);
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0] * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5);

		  this->ThreeBodyInteractionFactors[i][Index] *= 2.0 * FactorU;

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
// kx1 = momentum along x for the A site
// ky1 = momentum along y for the A site
// kx2 = momentum along x for the B site
// ky2 = momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2)
{
  Complex Tmp = 2.0 * (cos (M_PI * ((((double) (kx2 - kx1)) / ((double) this->NbrSiteX)) - ((((double) (ky2 - ky1)) / ((double) this->NbrSiteY))))) 
		       + cos (M_PI * ((((double) (kx2 - kx1)) / ((double) this->NbrSiteX)) + ((((double) (ky2 - ky1)) / ((double) this->NbrSiteY))))));
  return Tmp;
}

// compute the matrix element for the three body interaction between two sites B and one site A
//
// kx1 = creation momentum along x for the first B site
// ky1 = creation momentum along y for the first B site
// kx2 = creation momentum along x for the second B site
// ky2 = creation momentum along y for the second B site
// kx3 = annihilation momentum along x for the first B site
// ky3 = annihilation momentum along y for the first B site
// kx4 = annihilation momentum along x for the second B site
// ky4 = annihilation momentum along y for the second B site
// return value = corresponding matrix element

Complex ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 2.0 * (cos (M_PI * ((((double) ((kx1 - kx3) - (kx2 - kx4))) / ((double) this->NbrSiteX)) + ((((double) ((ky1 - ky3) + (ky2 - ky4))) / ((double) this->NbrSiteY))))) 
		       + cos (M_PI * ((((double) ((kx1 - kx3) + (kx2 - kx4))) / ((double) this->NbrSiteX)) - ((((double) ((ky1 - ky3) - (ky2 - ky4))) / ((double) this->NbrSiteY))))));
  return Tmp;
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

Complex ParticleOnLatticeCheckerboardLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementBAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 0.0;
  return Tmp;
}
