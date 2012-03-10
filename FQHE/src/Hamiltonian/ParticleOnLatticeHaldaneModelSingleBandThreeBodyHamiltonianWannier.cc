 ////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                           class author: Yang-Le Wu                         //
//                                                                            //
//               class of Haldane model with interacting particles            //
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 16/08/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier.h"
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

ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier()
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
// t3 = hoping amplitude between next to next nearest neighbor sites
// phi =  Haldane phase on nnn hopping
// mus = sublattice staggered chemical potential 
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, 
        double uPotential, double vPotential, double wPotential, double sPotential,
															 double t1, double t2, double t3, double phi, double mus, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
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
  this->NextNextNNHoping = t3;
  this->HaldanePhase = phi;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->WPotential = wPotential;
  this->SPotential = sPotential;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  this->TwoBodyFlag = true;
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

ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::~ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier()
{
  delete[] this->XPhaseTable;
  delete[] this->XHalfPhaseTable;
  delete[] this->YPhaseTable;
  delete[] this->YHalfPhaseTable;
}

// evaluate all interaction factors
//   

void ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
   
  if (this->FlatBand == false)
    this->OneBodyInteractionFactors = new double [this->NbrSiteX * this->NbrSiteY];
  this->ComputeOneBodyMatrices(OneBodyBasis);
	      	  cout << "End of ComputeOneBodyMatrices" <<endl; 
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
                      this->InteractionFactors[i][Index]  = FactorU * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                                                                    * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
                      this->InteractionFactors[i][Index] -= FactorU * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                                                                    * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);
                      this->InteractionFactors[i][Index] -= FactorU * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1])
                                                                    * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
                      this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1])
                                                                    * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);

                      this->InteractionFactors[i][Index] += FactorV * (  (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0])
                                                                       + (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]))
                                                                    * this->ComputeTwoBodyMatrixElementAA(kx2, ky2, kx4, ky4);
                      this->InteractionFactors[i][Index] -= FactorV * (  (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0])
                                                                       + (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]))
                                                                    * this->ComputeTwoBodyMatrixElementAA(kx1, ky1, kx4, ky4);
                      this->InteractionFactors[i][Index] -= FactorV * (  (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0])
                                                                       + (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1]))
                                                                    * this->ComputeTwoBodyMatrixElementAA(kx2, ky2, kx3, ky3);
                      this->InteractionFactors[i][Index] += FactorV * (  (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0])
                                                                       + (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1]))
                                                                    * this->ComputeTwoBodyMatrixElementAA(kx1, ky1, kx3, ky3);

                      this->InteractionFactors[i][Index] *= -2.0;

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


                  sumW += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6);
                  sumW -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5);
                  sumW -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6);
                  sumW += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4);
                  sumW += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5);
                  sumW -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4);
                  sumW -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6);
                  sumW += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5);
                  sumW += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6);
                  sumW -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4);
                  sumW -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5);
                  sumW += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4);
                  sumW -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6);
                  sumW += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5);
                  sumW += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6);
                  sumW -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4);
                  sumW -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5);
                  sumW += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4);
                  sumW += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6);
                  sumW -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5);
                  sumW -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6);
                  sumW += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4);
                  sumW += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5);
                  sumW -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4);
                  sumW += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6);
                  sumW -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5);
                  sumW -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6);
                  sumW += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4);
                  sumW += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5);
                  sumW -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4);
                  sumW -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6);
                  sumW += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5);
                  sumW += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6);
                  sumW -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4);
                  sumW -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5);
                  sumW += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                        * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4);

                  sumW += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6));
                  sumW -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5));
                  sumW -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6));
                  sumW += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4));
                  sumW += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5));
                  sumW -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4));
                  sumW -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6));
                  sumW += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5));
                  sumW += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6));
                  sumW -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4));
                  sumW -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5));
                  sumW += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4));
                  sumW -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6));
                  sumW += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5));
                  sumW += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6));
                  sumW -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4));
                  sumW -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5));
                  sumW += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4));
                  sumW += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6));
                  sumW -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5));
                  sumW -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6));
                  sumW += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4));
                  sumW += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5));
                  sumW -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4));
                  sumW += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6));
                  sumW -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5));
                  sumW -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6));
                  sumW += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4));
                  sumW += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5));
                  sumW -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4));
                  sumW -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6));
                  sumW += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5));
                  sumW += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6));
                  sumW -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4));
                  sumW -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5));
                  sumW += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                        * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                        * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                        * Conj(this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4));






                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6);

                  sumS -= (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5);

                  sumS -= (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5);

                  sumS -= (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4);

                  sumS -= (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6);

                  sumS -= (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4);

                  sumS -= (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5);

                  sumS += (  Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4);

                  sumS -= (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6);

                  sumS -= (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4);

                  sumS -= (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6);

                  sumS -= (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5);

                  sumS -= (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4);

                  sumS += (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5);

                  sumS -= (  Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6);

                  sumS -= (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5);

                  sumS -= (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5);

                  sumS -= (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4);

                  sumS -= (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6);

                  sumS -= (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4);

                  sumS -= (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5);

                  sumS += (  Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                           * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                           * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                          +  Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                           * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                           * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])
                        * this->ComputeThreeBodyMatrixElementAAA(kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4);



                  this->NBodyInteractionFactors[i][Index] = 0.;
                  this->NBodyInteractionFactors[i][Index] += 2.0 * FactorW * sumW;
                  this->NBodyInteractionFactors[i][Index] += 2.0 * FactorS * sumS;

		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
      cout << "nbr 3-body interaction = " << TotalNbrInteractionFactors << endl;
      cout << "====================================" << endl;
    }
  // Bosonic case
  else
    {

      // Two-Body
      this->NbrSectorSums = this->NbrSiteY;
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
      	this->NbrSectorIndicesPerSum[i] = 0;      
      for (int x1 = 0; x1 < this->NbrSiteX; ++x1)
      	for (int x2 = 0; x2 < this->NbrSiteX; ++x2)
      	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
      	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
      	      {
      		int Index1 = (x1 * this->NbrSiteY) + ky1;
      		int Index2 = (x2 * this->NbrSiteY) + ky2;
      		if (Index1 <= Index2)
      		  ++this->NbrSectorIndicesPerSum[((ky1 + ky2) % this->NbrSiteY)];    
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
      for (int x1 = 0; x1 < this->NbrSiteX; ++x1)
      	for (int x2 = 0; x2 < this->NbrSiteX; ++x2)
      	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
      	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
      	      {
      		int Index1 = (x1 * this->NbrSiteY) + ky1;
      		int Index2 = (x2 * this->NbrSiteY) + ky2;
      		if (Index1 <= Index2)
      		  {
      		    int TmpSum =  ((ky1 + ky2) % this->NbrSiteY);
      		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
      		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
      		    ++this->NbrSectorIndicesPerSum[TmpSum];    
      		  }
      	      }


      	
      	  double FactorU = this->UPotential*0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));

	  
      	  this->InteractionFactors = new Complex* [this->NbrSectorSums];
      	  for (int i = 0; i < this->NbrSectorSums; ++i)
      	    {
      	      this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
      	      int Index = 0;
      	      for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
      		{
      		  int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
      		  int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
      		  int x1 = Index1 / this->NbrSiteY;
      		  int ky1 = Index1 % this->NbrSiteY;
      		  int x2 = Index2 / this->NbrSiteY;
      		  int ky2 = Index2 % this->NbrSiteY;
      		  for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
      		    {
      		      int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
      		      int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
      		      int x3 = Index3 / this->NbrSiteY;
      		      int ky3 = Index3 % this->NbrSiteY;
      		      int x4 = Index4 / this->NbrSiteY;
      		      int ky4 = Index4 % this->NbrSiteY;

      		      Complex sumU=0.0;

		      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
			for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
			  for (int kx3 = 0; kx3 < this->NbrSiteX; ++kx3)
			    for (int kx4 = 0; kx4 < this->NbrSiteX; ++kx4) 
			      {
				if( (kx1+kx2-kx3-kx4) %  this->NbrSiteX == 0)
				  {

				int IndexBis1 = kx1 * this->NbrSiteY + ky1;
				int IndexBis2 = kx2 * this->NbrSiteY + ky2;
				int IndexBis3 = kx3 * this->NbrSiteY + ky3;
				int IndexBis4 = kx4 * this->NbrSiteY + ky4;

				Complex tmpsumU;
			      
      		      tmpsumU = (Conj(OneBodyBasis[IndexBis1][0][0]) * OneBodyBasis[IndexBis3][0][0] * Conj(OneBodyBasis[IndexBis2][0][0]) * OneBodyBasis[IndexBis4][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
      		      tmpsumU += (Conj(OneBodyBasis[IndexBis2][0][0]) * OneBodyBasis[IndexBis3][0][0] * Conj(OneBodyBasis[IndexBis1][0][0]) * OneBodyBasis[IndexBis4][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
      		      tmpsumU += (Conj(OneBodyBasis[IndexBis1][0][0]) * OneBodyBasis[IndexBis4][0][0] * Conj(OneBodyBasis[IndexBis2][0][0]) * OneBodyBasis[IndexBis3][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();
      		      tmpsumU += (Conj(OneBodyBasis[IndexBis2][0][0]) * OneBodyBasis[IndexBis4][0][0] * Conj(OneBodyBasis[IndexBis1][0][0]) * OneBodyBasis[IndexBis3][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteAA();

      		      tmpsumU += (Conj(OneBodyBasis[IndexBis1][0][1]) * OneBodyBasis[IndexBis3][0][1] * Conj(OneBodyBasis[IndexBis2][0][1]) * OneBodyBasis[IndexBis4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
      		      tmpsumU += (Conj(OneBodyBasis[IndexBis2][0][1]) * OneBodyBasis[IndexBis3][0][1] * Conj(OneBodyBasis[IndexBis1][0][1]) * OneBodyBasis[IndexBis4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
      		      tmpsumU += (Conj(OneBodyBasis[IndexBis1][0][1]) * OneBodyBasis[IndexBis4][0][1] * Conj(OneBodyBasis[IndexBis2][0][1]) * OneBodyBasis[IndexBis3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
      		      tmpsumU += (Conj(OneBodyBasis[IndexBis2][0][1]) * OneBodyBasis[IndexBis4][0][1] * Conj(OneBodyBasis[IndexBis1][0][1]) * OneBodyBasis[IndexBis3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteBB(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

		      tmpsumU *= Conj(Conj(WannierBasis[ky1][x1][kx1])*Conj(WannierBasis[ky2][x2][kx2])*WannierBasis[ky3][x3][kx3]*WannierBasis[ky4][x4][kx4]);

		      sumU += tmpsumU;

				  }
			      }


      		      if (Index3 == Index4)
      			{
      			  sumU *= 0.5;
      			}
      		      if (Index1 == Index2)
      			{
      			  sumU *= 0.5;
      			}

      		      this->InteractionFactors[i][Index] = 2.0 * FactorU * sumU ;

      		      TotalNbrInteractionFactors++;
      		      ++Index;

      		    }
      		}
      	    }
      	



      // Three Body
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
		    if ((Index1 <= Index2) && (Index2 <= Index3))
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
		    if ((Index1 <= Index2) && (Index2 <= Index3))
		      {
			int TmpSum = (((kx1 + kx2 + kx3) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2 + ky3) % this->NbrSiteY);
			this->NBodySectorIndicesPerSum[TmpSum][this->NbrNBodySectorIndicesPerSum[TmpSum] * 3] = Index1;
			this->NBodySectorIndicesPerSum[TmpSum][1 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = Index2;
			this->NBodySectorIndicesPerSum[TmpSum][2 + (this->NbrNBodySectorIndicesPerSum[TmpSum] * 3)] = Index3;
			++this->NbrNBodySectorIndicesPerSum[TmpSum];    
		      }
		  }


      int KxIn[3];
      int KyIn[3];
      int IndexIn[3];

      int** Permutations = 0; 
      double* PermutationSign = 0; 
      int NbrPermutations = this->ComputePermutations(Permutations, PermutationSign);

      int TmpLargestSector = 0;
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	if (this->NbrNBodySectorIndicesPerSum[i] > TmpLargestSector)
	  TmpLargestSector = this->NbrNBodySectorIndicesPerSum[i];

      Complex* TmpAAAIn = new Complex[TmpLargestSector];
      Complex* TmpAAAOut = new Complex[TmpLargestSector];
      Complex* TmpBBBIn = new Complex[TmpLargestSector];
      Complex* TmpBBBOut = new Complex[TmpLargestSector];

      // double FactorU = this->UPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      double ThreeBodyFactor = this->WPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      IndexIn[0] = this->NBodySectorIndicesPerSum[i][j1 * 3];
	      IndexIn[1] = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 1];
	      IndexIn[2] = this->NBodySectorIndicesPerSum[i][(j1 * 3) + 2];
	      KxIn[0] = IndexIn[0] / this->NbrSiteY;
	      KyIn[0] = IndexIn[0] % this->NbrSiteY;
	      KxIn[1] = IndexIn[1] / this->NbrSiteY;
	      KyIn[1] = IndexIn[1] % this->NbrSiteY;
	      KxIn[2] = IndexIn[2] / this->NbrSiteY;
	      KyIn[2] = IndexIn[2] % this->NbrSiteY;

	      Complex TmpAAAIn2 = 0.0;
	      Complex TmpAAAOut2 = 0.0;
	      Complex TmpBBBIn2 = 0.0;
	      Complex TmpBBBOut2 = 0.0;


	      for (int l1 = 0; l1 < NbrPermutations; ++l1)
		{
		  int* TmpPerm = Permutations[l1];
		  Complex SumTmpAAAIn2 = Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0]) ;
 		  Complex SumTmpBBBIn2 = Conj(OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][1]) ;
		  
 		  Complex SumTmpAAAOut2 = OneBodyBasis[IndexIn[TmpPerm[0]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][0] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][0] ;
 		  Complex SumTmpBBBOut2 = OneBodyBasis[IndexIn[TmpPerm[0]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[1]]][0][1] * OneBodyBasis[IndexIn[TmpPerm[2]]][0][1] ;

				  
		  if ((IndexIn[TmpPerm[0]] == IndexIn[TmpPerm[1]]) && (IndexIn[TmpPerm[1]] == IndexIn[TmpPerm[2]]))
		    {
		      SumTmpAAAIn2 *= 1.0/6.0;
		      SumTmpBBBIn2 *= 1.0/6.0;
		      SumTmpAAAOut2 *= 1.0/6.0;
		      SumTmpBBBOut2 *= 1.0/6.0;
		      
		    }
		  else
		    {
		      if( (IndexIn[TmpPerm[0]] == IndexIn[TmpPerm[1]])  || (IndexIn[TmpPerm[1]] == IndexIn[TmpPerm[2]])  || (IndexIn[TmpPerm[0]] == IndexIn[TmpPerm[2]]) )
			{
            	      	  SumTmpAAAIn2 *= 0.5;
			  SumTmpBBBIn2 *= 0.5;
			  SumTmpAAAOut2 *= 0.5;
			  SumTmpBBBOut2 *= 0.5;
			}
		    }
		  
 		  TmpAAAIn2 += SumTmpAAAIn2 ;
 		  TmpBBBIn2 += SumTmpBBBIn2 ;
		  
 		  TmpAAAOut2 += SumTmpAAAOut2 ;
 		  TmpBBBOut2 += SumTmpBBBOut2 ;


 		}

 	      TmpAAAIn[j1] =  TmpAAAIn2;
 	      TmpAAAOut[j1] = TmpAAAOut2;
 	      TmpBBBIn[j1] =  TmpBBBIn2;
 	      TmpBBBOut[j1] = TmpBBBOut2;


	    }

	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  this->NBodyInteractionFactors[i][Index] = 2.0 * ThreeBodyFactor * (TmpAAAIn[j1] * TmpAAAOut[j2]  + TmpBBBIn[j1] * TmpBBBOut[j2]);							    
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }

	}
      for (int i=0; i<NbrPermutations; ++i)
	delete [] Permutations[i];
      delete [] Permutations;
      delete [] PermutationSign;


      delete[] TmpAAAIn;
      delete[] TmpAAAOut;
      delete[] TmpBBBIn;
      delete[] TmpBBBOut;

    }
}


// compute the matrix element for the two body interaction between two sites A and B 
//
// kx1 = annihilation momentum along x for the B site
// ky1 = annihilation momentum along y for the B site
// kx2 = creation momentum along x for the B site
// ky2 = creation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2)
{
  double dx = ((double)(kx1-kx2)) * this->KxFactor;
  double dy = ((double)(ky1-ky2)) * this->KyFactor;
  Complex Tmp = 1 + Phase(dx + dy) + Phase(dy);
  return Tmp;
}

// compute the matrix element for the two body interaction between two A sites (or two B sites) 
//
// kx1 = annihilation momentum along x for the second site
// ky1 = annihilation momentum along y for the second site
// kx2 = creation momentum along x for the second site
// ky2 = creation momentum along y for the second site
// return value = corresponding matrix element

Complex ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::ComputeTwoBodyMatrixElementAA(int kx1, int ky1, int kx2, int ky2)
{
  double dx = ((double)(kx1-kx2)) * this->KxFactor;
  double dy = ((double)(ky1-ky2)) * this->KyFactor;
  Complex Tmp = Phase(dx) + Phase(dy) + Phase(dx + dy);
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

Complex ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::ComputeThreeBodyMatrixElementABB(int kx2, int ky2, int kx3, int ky3, int kx5, int ky5, int kx6, int ky6)
{
  //double dx2 = ((double)(kx2 - kx5)) * this->KxFactor;
  double dx3 = ((double)(kx3 - kx6)) * this->KxFactor;
  double dy2 = ((double)(ky2 - ky5)) * this->KyFactor;
  double dy3 = ((double)(ky3 - ky6)) * this->KyFactor;
  Complex Tmp = Phase(dy2 + dx3 + dy3) + Phase(dx3 + dy3) + Phase(dy3);
  return Tmp;
}

// compute the matrix element for the three body interaction between NNN sites 
//
// kx2 = annihilation momentum along x for the second site
// ky2 = annihilation momentum along y for the second site
// kx3 = annihilation momentum along x for the third site
// ky3 = annihilation momentum along y for the third site
// kx5 = creation momentum along x for the second site
// ky5 = creation momentum along y for the second site
// kx6 = creation momentum along x for the third site
// ky6 = creation momentum along y for the third site
// return value = corresponding matrix element

Complex ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::ComputeThreeBodyMatrixElementAAA(int kx2, int ky2, int kx3, int ky3, int kx5, int ky5, int kx6, int ky6)
{
    double dx2 = ((double)(kx2 - kx5)) * this->KxFactor;
    double dx3 = ((double)(kx3 - kx6)) * this->KxFactor;
    double dy2 = ((double)(ky2 - ky5)) * this->KyFactor;
    double dy3 = ((double)(ky3 - ky6)) * this->KyFactor;
    Complex Tmp = Phase(dx2 + dy2) * (Phase(dx3) + Phase(dy3));
    return Tmp;
}

// compute the matrix element for on-site two body interaction involving A sites
//
// return value = corresponding matrix element

Complex ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::ComputeTwoBodyMatrixElementOnSiteAA()
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

Complex ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::ComputeTwoBodyMatrixElementOnSiteBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return 1.0;
}

// compute the matrix element for the on-site three body interaction related to sites A
//
// kx1 = first annihilation momentum along x for the first A site
// ky1 = first annihilation momentum along y for the first A site
// kx2 = second annihilation momentum along x for the second A site
// ky2 = second annihilation momentum along y for the second A site
// kx3 = third annihilation momentum along x for the second A site
// ky3 = third annihilation momentum along y for the second A site
// kx4 = first creation momentum along x for the first A site
// ky4 = first creation momentum along y for the first A site
// kx5 = second creation momentum along x for the second A site
// ky5 = second creation momentum along y for the second A site
// kx6 = third creation momentum along x for the second A site
// ky6 = third creation momentum along y for the second A site
// return value = corresponding matrix element

Complex ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::ComputeThreeBodyMatrixElementOnSiteAAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
  return 1.0;
}

// compute the matrix element for the on-site three body interaction related to sites A
//
// kx1 = first annihilation momentum along x for the first A site
// ky1 = first annihilation momentum along y for the first A site
// kx2 = second annihilation momentum along x for the second A site
// ky2 = second annihilation momentum along y for the second A site
// kx3 = third annihilation momentum along x for the second A site
// ky3 = third annihilation momentum along y for the second A site
// kx4 = first creation momentum along x for the first A site
// ky4 = first creation momentum along y for the first A site
// kx5 = second creation momentum along x for the second A site
// ky5 = second creation momentum along y for the second A site
// kx6 = third creation momentum along x for the second A site
// ky6 = third creation momentum along y for the second A site
// return value = corresponding matrix element

Complex ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::ComputeThreeBodyMatrixElementOnSiteBBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4, int kx5, int ky5, int kx6, int ky6)
{
//   Complex Tmp = Phase (0.5 * ((((double) (kx6 + kx5 + kx4 - kx3 - kx2 - kx1)) * this->KxFactor)
// 			      + ((((double) (ky6 + ky5 + ky4 - ky3 - ky2 - ky1)) * this->KyFactor))));
//   return Tmp;

    //   double dx = ((double)(-kx1-kx2-kx3+kx4+kx5+kx6)) * this->KxFactor;
    //    double dy = ((double)(-ky1-ky2-ky3+ky4+ky5+ky6)) * this->KyFactor;
    //    Complex Tmp = Phase((-dx-2.0*dy)/3.0);
    //return Tmp;
  
  return 1.0;
}


// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = array of one body transformation matrices

void ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
  {
    double kxl=((double)kx + this->GammaX) * this->KxFactor;
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
        double kyl=((double)ky + this->GammaY) * this->KyFactor;
	int Index = (kx * this->NbrSiteY) + ky;

	// My Convention
        // Complex B1 = - this->NNHoping * Complex(1 + cos(kxl) + cos(kyl), - sin(kxl) - sin(kyl));
	// Complex B2 = - this->NextNextNNHoping * Complex(cos(kxl+kyl)+2*cos(kxl-kyl),-sin(kxl+kyl) );
        // double d0 = - 2.0 * this->NextNNHoping * cos(this->HaldanePhase) * (cos(kxl) + cos(kyl) + cos(kxl-kyl));
        // double d3 = - 2.0 * this->NextNNHoping * sin(this->HaldanePhase) * (sin(kxl) - sin(kyl) - sin(kxl-kyl)) + this->MuS;

        Complex B1 = - this->NNHoping * Complex(1 + cos(kxl+kyl) + cos(kyl), + sin(kxl+kyl) + sin(kyl));
	Complex B2 = - this->NextNextNNHoping * Complex(2* cos(kxl) + cos(kxl+2*kyl),  sin(kxl+2*kyl));
        double d0 = - 2.0 * this->NextNNHoping * cos(this->HaldanePhase) * (cos(kxl) + cos(kyl) + cos(kxl+kyl));
        double d3 = - 2.0 * this->NextNNHoping * sin(this->HaldanePhase) * (sin(kxl) + sin(kyl) - sin(kxl+kyl)) + this->MuS;

	HermitianMatrix TmpOneBodyHamiltonian(2, true);
	TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d0 + d3);
	TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B1+B2);
	TmpOneBodyHamiltonian.SetMatrixElement(1, 1, d0 - d3);
	ComplexMatrix TmpMatrix(2, 2, true);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	oneBodyBasis[Index] = TmpMatrix;	
	// if (this->FlatBand == false)
	//   {
	//     this->OneBodyInteractionFactors[Index] = 0.5*TmpDiag(0, 0);
	//     // The 1/2 factor is needed because this one body term is counted twice somewhere in ChernInsulatorSingleBandHamiltonian because it was written for spin half particles. This factor is also present for the Kagome Lattice
	//   }
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << "  e1=[" << TmpMatrix[0][0] << ", " << TmpMatrix[0][1] << "]  e2=[" << TmpMatrix[1][0] << ", " << TmpMatrix[1][1] << "]" << endl;
      }
  }


  // BerryConnection 
  this->BerryConnectionX = new double* [this->NbrSiteY];
  this->BerryConnectionXCumSum = new double* [this->NbrSiteY];
  this->Theta = new double [this->NbrSiteY];

  for (int ky = 0; ky < this->NbrSiteY; ++ky)
    {
      double TempTheta = 0;
      this->BerryConnectionX[ky] = new double[this->NbrSiteX];
      this->BerryConnectionXCumSum[ky] = new double[this->NbrSiteX];
      for (int kx = 0; kx < this->NbrSiteX; ++kx)
	{
	  int Index1=kx*this->NbrSiteY+ky;
	  int Index2;
	  if(kx!=this->NbrSiteX-1)
	    { 
	      Index2=(kx+1)*this->NbrSiteY+ky;
	    }
	  else
	    {
	      Index2=ky;
	    }
	  double ArgTemp = -Arg(oneBodyBasis[Index2][0][0]*Conj(oneBodyBasis[Index1][0][0])+oneBodyBasis[Index2][0][1]*Conj(oneBodyBasis[Index1][0][1]));
	  this->BerryConnectionX[ky][kx]= ArgTemp ;
	  this->BerryConnectionXCumSum[ky][kx]= TempTheta ;
	  TempTheta += ArgTemp;
	}
      if(ky!=0)
	{
	  if(TempTheta > this->Theta[ky-1])
	    TempTheta -= 2.0 * M_PI ;
	}
      this->Theta[ky]=TempTheta;
      cout << "ky = " << ky << endl;
      cout << "theta = " << TempTheta << endl;
    }

  //WannierBasis

  this->WannierBasis = new Complex** [this->NbrSiteY];
  for (int ky = 0; ky < this->NbrSiteY; ++ky)
    {
      double kyl=((double)ky + this->GammaY) * this->KyFactor;
      this->WannierBasis[ky] = new Complex* [this->NbrSiteX];
      for (int x = 0; x < this->NbrSiteX; ++x)
	{
	  this->WannierBasis[ky][x] = new Complex [this->NbrSiteX];
	  for (int kx = 0; kx < this->NbrSiteX; ++kx)
	    {
	      Complex Temp;
	      double kxl=((double)kx + this->GammaX) * this->KxFactor;
	      if(kx==0)
		{
		  Temp=1.0;
		}
	      else
		{
		  Temp = Phase(-this->BerryConnectionXCumSum[ky][kx] - kxl * ((double) x - this->Theta[ky] /(2.0*M_PI)));
		}

	      Temp *= 1.0 / sqrt(this->NbrSiteX);
	      this->WannierBasis[ky][x][kx] = Temp ;
	    }
	}
    }




}
// compute all the phase precalculation arrays 
//

void ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonianWannier::ComputePhaseArray()
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
