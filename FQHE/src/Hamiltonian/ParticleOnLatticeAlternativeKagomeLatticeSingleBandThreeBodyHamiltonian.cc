////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                           class author: Yang-Le Wu                         //
//                                                                            //
//            class of Kagome lattice model with interacting particles        //
//         in the single band approximation with three body interaction       // 
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
#include "Hamiltonian/ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian.h"
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

// default constructor
//

ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian::ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = strength of the repulsive two body nearest neighbor interaction
// vPotential = strength of the repulsive two body second nearest neighbor interaction
// wPotential = strength of the repulsive three body nearest neighbor interaction
// sPotential = strength of the repulsive three body next-to-nearest neighbor interaction
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

ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian::ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, 
																		 double uPotential, double vPotential, double wPotential, double sPotential,
																		 double t1, double t2, double l1, double l2, double gammaX, double gammaY, int bandIndex, 
																		 bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NBodyValue = 3;
  this->ComputePhaseArray();

  this->HamiltonianShift = 0.0;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->NNHoping = t1;
  this->NextNNHoping = t2;
  this->NNRashba = l1;
  this->NextNNRashba = l2;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->BandIndex = bandIndex;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->WPotential = wPotential;
  this->SPotential = sPotential;
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

ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian::~ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian()
{
  delete[] this->XPhaseTable;
  delete[] this->XHalfPhaseTable;
  delete[] this->YPhaseTable;
  delete[] this->YHalfPhaseTable;
}

// evaluate all interaction factors
//   

void ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian::EvaluateInteractionFactors()
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
      this->SectorIndicesPerSum = new int* [this->NbrSectorSums];
      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      if (this->TwoBodyFlag == true)
      {
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
          double FactorU = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
          double FactorV = this->VPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
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

                      sumU += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0]
                            * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                            * this->ComputeTwoBodyMatrixElementNNAB(kx2, ky2, kx4, ky4);
                      sumU += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1]
                            * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2]
                            * this->ComputeTwoBodyMatrixElementNNBC(kx2, ky2, kx4, ky4);
                      sumU += Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2]
                            * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                            * this->ComputeTwoBodyMatrixElementNNCA(kx2, ky2, kx4, ky4);
                      sumU -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                            * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1]
                            * this->ComputeTwoBodyMatrixElementNNAB(kx2, ky2, kx3, ky3);
                      sumU -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                            * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2]
                            * this->ComputeTwoBodyMatrixElementNNBC(kx2, ky2, kx3, ky3);
                      sumU -= Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2]
                            * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0]
                            * this->ComputeTwoBodyMatrixElementNNCA(kx2, ky2, kx3, ky3);
                      sumU -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0]
                            * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                            * this->ComputeTwoBodyMatrixElementNNAB(kx1, ky1, kx4, ky4);
                      sumU -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1]
                            * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2]
                            * this->ComputeTwoBodyMatrixElementNNBC(kx1, ky1, kx4, ky4);
                      sumU -= Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2]
                            * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                            * this->ComputeTwoBodyMatrixElementNNCA(kx1, ky1, kx4, ky4);
                      sumU += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                            * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1]
                            * this->ComputeTwoBodyMatrixElementNNAB(kx1, ky1, kx3, ky3);
                      sumU += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                            * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2]
                            * this->ComputeTwoBodyMatrixElementNNBC(kx1, ky1, kx3, ky3);
                      sumU += Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2]
                            * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0]
                            * this->ComputeTwoBodyMatrixElementNNCA(kx1, ky1, kx3, ky3);



                      sumV += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0]
                            * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                            * this->ComputeTwoBodyMatrixElementNNNAB(kx2, ky2, kx4, ky4);
                      sumV += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1]
                            * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2]
                            * this->ComputeTwoBodyMatrixElementNNNBC(kx2, ky2, kx4, ky4);
                      sumV += Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2]
                            * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                            * this->ComputeTwoBodyMatrixElementNNNCA(kx2, ky2, kx4, ky4);
                      sumV -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                            * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1]
                            * this->ComputeTwoBodyMatrixElementNNNAB(kx2, ky2, kx3, ky3);
                      sumV -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                            * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2]
                            * this->ComputeTwoBodyMatrixElementNNNBC(kx2, ky2, kx3, ky3);
                      sumV -= Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2]
                            * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0]
                            * this->ComputeTwoBodyMatrixElementNNNCA(kx2, ky2, kx3, ky3);
                      sumV -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0]
                            * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                            * this->ComputeTwoBodyMatrixElementNNNAB(kx1, ky1, kx4, ky4);
                      sumV -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1]
                            * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2]
                            * this->ComputeTwoBodyMatrixElementNNNBC(kx1, ky1, kx4, ky4);
                      sumV -= Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2]
                            * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                            * this->ComputeTwoBodyMatrixElementNNNCA(kx1, ky1, kx4, ky4);
                      sumV += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                            * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1]
                            * this->ComputeTwoBodyMatrixElementNNNAB(kx1, ky1, kx3, ky3);
                      sumV += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                            * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2]
                            * this->ComputeTwoBodyMatrixElementNNNBC(kx1, ky1, kx3, ky3);
                      sumV += Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2]
                            * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0]
                            * this->ComputeTwoBodyMatrixElementNNNCA(kx1, ky1, kx3, ky3);


                      this->InteractionFactors[i][Index] = -2.0 * (FactorU * sumU + FactorV * sumV);


                      TotalNbrInteractionFactors++;
                      ++Index;
                  }
              }
          }
      }
      cout << "nbr 2-body interaction = " << TotalNbrInteractionFactors << endl;
      TotalNbrInteractionFactors=0;


      this->NbrNBodySectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrNBodySectorIndicesPerSum = new int[this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	this->NbrNBodySectorIndicesPerSum[i] = 0;      
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
		      ++this->NbrNBodySectorIndicesPerSum[(((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY)];    
		  }
      this->NBodySectorIndicesPerSum = new int* [this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  if (this->NbrNBodySectorIndicesPerSum[i]  > 0)
	    {
	      this->NBodySectorIndicesPerSum[i] = new int[this->NBodyValue * this->NbrNBodySectorIndicesPerSum[i]];      
	      this->NbrNBodySectorIndicesPerSum[i] = 0;
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
			this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum] * 3] = Index1;
			this->NBodySectorIndicesPerSum[TmpSum][1 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = Index2;
			this->NBodySectorIndicesPerSum[TmpSum][2 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = Index3;
			++this->NbrNBodySectorIndicesPerSum[TmpSum];    
		      }
		  }

      double FactorW = this->WPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), this->NBodyValue - 1);
      double FactorS = this->SPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), this->NBodyValue - 1);
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1) // annihilation operators
	    {
	      int Index1 = this->NBodySectorIndicesPerSum[i][j1 * 3];
	      int Index2 = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 1];
	      int Index3 = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 2];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      int kx3 = Index3 / this->NbrSiteY;
	      int ky3 = Index3 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2) // creation operators
		{
		  int Index4 = this->NBodySectorIndicesPerSum[i][j2 * 3];
		  int Index5 = this->NBodySectorIndicesPerSum[i][(j2 * 3) + 1];
		  int Index6 = this->NBodySectorIndicesPerSum[i][(j2 * 3) + 2];
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;
		  int kx5 = Index5 / this->NbrSiteY;
		  int ky5 = Index5 % this->NbrSiteY;
		  int kx6 = Index6 / this->NbrSiteY;
		  int ky6 = Index6 % this->NbrSiteY;

                  Complex sumW=0.;
                  Complex sumS=0.;

                  sumW += Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index5][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][2]) * OneBodyBasis[Index6][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6);
                        
                  sumW -= Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index6][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][2]) * OneBodyBasis[Index5][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5);
                        
                  sumW -= Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index5][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][2]) * OneBodyBasis[Index6][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6);
                        
                  sumW += Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index5][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index6][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4);
                        
                  sumW += Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index6][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][2]) * OneBodyBasis[Index5][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5);
                        
                  sumW -= Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index6][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index5][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4);
                        
                  sumW -= Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][1]) * OneBodyBasis[Index5][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index6][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6);
                        
                  sumW += Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][1]) * OneBodyBasis[Index6][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index5][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5);
                        
                  sumW += Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index5][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index6][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6);
                        
                  sumW -= Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index5][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][1]) * OneBodyBasis[Index6][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4);
                        
                  sumW -= Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index6][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index5][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5);
                        
                  sumW += Conj(OneBodyBasis[Index1][this->BandIndex][0]) * OneBodyBasis[Index6][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][1]) * OneBodyBasis[Index5][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4);
                        
                  sumW -= Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index5][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][2]) * OneBodyBasis[Index6][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6);
                        
                  sumW += Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index6][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][2]) * OneBodyBasis[Index5][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5);
                        
                  sumW += Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index5][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][2]) * OneBodyBasis[Index6][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6);
                        
                  sumW -= Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index5][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index6][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4);
                        
                  sumW -= Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index6][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][2]) * OneBodyBasis[Index5][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5);
                        
                  sumW += Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index6][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index5][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4);
                        
                  sumW += Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][1]) * OneBodyBasis[Index5][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index6][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6);
                        
                  sumW -= Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][1]) * OneBodyBasis[Index6][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index5][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5);
                        
                  sumW -= Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index5][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index6][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6);
                        
                  sumW += Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index5][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][1]) * OneBodyBasis[Index6][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4);
                        
                  sumW += Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index6][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index5][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5);
                        
                  sumW -= Conj(OneBodyBasis[Index2][this->BandIndex][0]) * OneBodyBasis[Index6][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index3][this->BandIndex][1]) * OneBodyBasis[Index5][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4);
                        
                  sumW += Conj(OneBodyBasis[Index3][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index5][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index6][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6);
                        
                  sumW -= Conj(OneBodyBasis[Index3][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index6][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index5][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5);
                        
                  sumW -= Conj(OneBodyBasis[Index3][this->BandIndex][0]) * OneBodyBasis[Index5][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index6][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6);
                        
                  sumW += Conj(OneBodyBasis[Index3][this->BandIndex][0]) * OneBodyBasis[Index5][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index6][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4);
                        
                  sumW += Conj(OneBodyBasis[Index3][this->BandIndex][0]) * OneBodyBasis[Index6][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index5][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5);
                        
                  sumW -= Conj(OneBodyBasis[Index3][this->BandIndex][0]) * OneBodyBasis[Index6][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][1]) * OneBodyBasis[Index5][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4);
                        
                  sumW -= Conj(OneBodyBasis[Index3][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index5][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index6][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6);
                        
                  sumW += Conj(OneBodyBasis[Index3][this->BandIndex][0]) * OneBodyBasis[Index4][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index6][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index5][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5);
                        
                  sumW += Conj(OneBodyBasis[Index3][this->BandIndex][0]) * OneBodyBasis[Index5][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index6][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6);
                        
                  sumW -= Conj(OneBodyBasis[Index3][this->BandIndex][0]) * OneBodyBasis[Index5][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index6][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4);
                        
                  sumW -= Conj(OneBodyBasis[Index3][this->BandIndex][0]) * OneBodyBasis[Index6][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index4][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index5][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5);
                        
                  sumW += Conj(OneBodyBasis[Index3][this->BandIndex][0]) * OneBodyBasis[Index6][this->BandIndex][0]
                        * Conj(OneBodyBasis[Index2][this->BandIndex][1]) * OneBodyBasis[Index5][this->BandIndex][1]
                        * Conj(OneBodyBasis[Index1][this->BandIndex][2]) * OneBodyBasis[Index4][this->BandIndex][2]
                        * this->ComputeThreeBodyMatrixElementNNABC(kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4);

                  this->NBodyInteractionFactors[i][Index] = 2.0 * (FactorW * sumW + FactorS * sumS);

		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
      cout << "nbr 3-body interaction = " << TotalNbrInteractionFactors << endl;
      cout << "====================================" << endl;
  }
}

// compute the matrix element for the two body interaction between two sites A and B 
//
// kx2 = annihilation momentum along x for the B site
// ky2 = annihilation momentum along y for the B site
// kx4 = creation momentum along x for the B site
// ky4 = creation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementNNAB(int kx2, int ky2, int kx4, int ky4)
{
  double dx = ((double)(kx2-kx4)) * this->KxFactor;
  double dy = ((double)(ky2-ky4)) * this->KyFactor;
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

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementNNBC(int kx2, int ky2, int kx4, int ky4)
{
  double dx = ((double)(kx2-kx4)) * this->KxFactor;
  double dy = ((double)(ky2-ky4)) * this->KyFactor;
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

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementNNCA(int kx2, int ky2, int kx4, int ky4)
{
  double dx = ((double)(kx2-kx4)) * this->KxFactor;
  double dy = ((double)(ky2-ky4)) * this->KyFactor;
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

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementNNNAB(int kx2, int ky2, int kx4, int ky4)
{
  double dx = ((double)(kx2-kx4)) * this->KxFactor;
  double dy = ((double)(ky2-ky4)) * this->KyFactor;
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

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementNNNBC(int kx2, int ky2, int kx4, int ky4)
{
  double dx = ((double)(kx2-kx4)) * this->KxFactor;
  double dy = ((double)(ky2-ky4)) * this->KyFactor;
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

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementNNNCA(int kx2, int ky2, int kx4, int ky4)
{
  double dx = ((double)(kx2-kx4)) * this->KxFactor;
  double dy = ((double)(ky2-ky4)) * this->KyFactor;
  Complex Tmp = Phase(dx) + Phase(dy - dx);
  return Tmp;
}

// compute the matrix element for the three body interaction between one site A and two sites B 
//
// kx2 = annihilation momentum along x for the first B site
// ky2 = annihilation momentum along y for the first B site
// kx3 = annihilation momentum along x for the second B site
// ky3 = annihilation momentum along y for the second B site
// kx5 = creation momentum along x for the first B site
// ky5 = creation momentum along y for the first B site
// kx6 = creation momentum along x for the second B site
// ky6 = creation momentum along y for the second B site
// return value = corresponding matrix element

Complex ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementNNABC(int kx2, int ky2, int kx3, int ky3, int kx5, int ky5, int kx6, int ky6)
{
    double dx2 = ((double)(kx2 - kx5)) * this->KxFactor;
    double dx3 = ((double)(kx3 - kx6)) * this->KxFactor;
    double dy2 = ((double)(ky2 - ky5)) * this->KyFactor;
    double dy3 = ((double)(ky3 - ky6)) * this->KyFactor;
    Complex Tmp = 1 + Phase(- dx2 - dy3);
    return Tmp;
}

// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = array of one body transformation matrices

void ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
  {
    double x=((double)kx + this->GammaX) * this->KxFactor;
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
        double y=((double)ky + this->GammaY) * this->KyFactor;
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
	oneBodyBasis[Index] = TmpMatrix;	
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactors[Index] = TmpDiag(0, 0);
	  }
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2) << "  e1=[" << TmpMatrix[0][0] << ", " << TmpMatrix[0][1] << ", " << TmpMatrix[0][2] << "]  e2=[" << TmpMatrix[1][0] << ", " << TmpMatrix[1][1] << ", " << TmpMatrix[1][2] << "]  e3=[" << TmpMatrix[2][0] << ", " << TmpMatrix[2][1] << ", " << TmpMatrix[2][2] << "]" << endl;
      }
  }
}
// compute all the phase precalculation arrays 
//

void ParticleOnLatticeAlternativeKagomeLatticeSingleBandThreeBodyHamiltonian::ComputePhaseArray()
{
  this->XPhaseTable = new Complex [2 * this->NBodyValue * this->NbrSiteX];
  this->XHalfPhaseTable = new Complex [2 * this->NBodyValue * this->NbrSiteX];
  this->XPhaseTableShift = this->NBodyValue * this->NbrSiteX;
  for (int i = -this->XPhaseTableShift; i < this->XPhaseTableShift; ++i)
    {
      this->XPhaseTable[this->XPhaseTableShift + i] = Phase(this->KxFactor * ((double) i));
      this->XHalfPhaseTable[this->XPhaseTableShift + i] = Phase(0.5 * this->KxFactor * ((double) i));
    }
  this->YPhaseTable = new Complex [2 * this->NBodyValue * this->NbrSiteY];
  this->YHalfPhaseTable = new Complex [2 * this->NBodyValue * this->NbrSiteY];
  this->YPhaseTableShift = this->NBodyValue * this->NbrSiteY;
  for (int i = -this->YPhaseTableShift; i < this->YPhaseTableShift; ++i)
    {
      this->YPhaseTable[this->YPhaseTableShift + i] = Phase(this->KyFactor * ((double) i));
      this->YHalfPhaseTable[this->YPhaseTableShift + i] = Phase(0.5 * this->KyFactor * ((double) i));
    }
}
