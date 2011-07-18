 ////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Yang-Le Wu                            //
//                                                                            //
//               class of Haldane model with interacting particles            //
//         in the single band approximation and three body interaction        // 
//                                                                            //
//                        last modification : 17/07/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian.h"
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
// vPotential = strength of the repulsive two body second neareast neighbor interaction
// wPotential = strength of the repulsive three body neareast neighbor interaction
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// phi =  Haldane phase on nnn hopping
// mus = sublattice staggered chemical potential 
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian::ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, int nbrSiteY, double uPotential, double vPotential, double wPotential, double t1, double t2, double phi, double mus, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->HamiltonianShift = 0.0;
  this->NNHoping = t1;
  this->NextNNHoping = t2;
  this->HaldanePhase = phi;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->WPotential = wPotential;
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

ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian::~ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian::EvaluateInteractionFactors()
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
        Complex B1 = this->NNHoping * Complex(1 + cos(x+y) + cos(y), + sin(x+y) + sin(y));
        double d0 = + 2.0 * this->NextNNHoping * cos(this->HaldanePhase) * (cos(x) + cos(y) + cos(x+y));
        double d3 = + 2.0 * this->NextNNHoping * sin(this->HaldanePhase) * (sin(x) + sin(y) - sin(x+y)) + this->MuS;
	HermitianMatrix TmpOneBobyHamiltonian(2, true);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 0, d0 + d3);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 1, B1);
	TmpOneBobyHamiltonian.SetMatrixElement(1, 1, d0 - d3);
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
  }
 
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


      double FactorW = this->WPotential * 0.5 / pow(((double) (this->NbrSiteX * this->NbrSiteY)), 2);
      if (this->FlatBand == false)
          FactorW *= this->WPotential;
      this->ThreeBodyInteractionFactors = new Complex* [this->NbrThreeBodySectorSums];
      for (int i = 0; i < this->NbrThreeBodySectorSums; ++i)
	{
	  this->ThreeBodyInteractionFactors[i] = new Complex[this->NbrThreeBodySectorIndicesPerSum[i] * this->NbrThreeBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrThreeBodySectorIndicesPerSum[i]; ++j1) // annihilation operators
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
	      for (int j2 = 0; j2 < this->NbrThreeBodySectorIndicesPerSum[i]; ++j2) // creation operators
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

                  this->ThreeBodyInteractionFactors[i][Index] = 0.;

                  // ABB
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * this->ComputeThreeBodyMatrixElementABB(kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4); 

                  // BAA
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx5, ky5, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx6, ky6, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx4, ky4, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx6, ky6, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx4, ky4, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx3, ky3, kx5, ky5, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx5, ky5, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx6, ky6, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx4, ky4, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx6, ky6, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx4, ky4, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx2, ky2, kx5, ky5, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx5, ky5, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx6, ky6, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx4, ky4, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx6, ky6, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx4, ky4, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx3, ky3, kx5, ky5, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx5, ky5, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx6, ky6, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx4, ky4, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx6, ky6, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx4, ky4, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index3][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx3, ky3, kx1, ky1, kx5, ky5, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx5, ky5, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx6, ky6, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx4, ky4, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx6, ky6, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx4, ky4, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx1, ky1, kx2, ky2, kx5, ky5, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx5, ky5, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index4][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx6, ky6, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx4, ky4, kx6, ky6); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index5][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index6][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx6, ky6, kx4, ky4); 
                  this->ThreeBodyInteractionFactors[i][Index] -= Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx4, ky4, kx5, ky5); 
                  this->ThreeBodyInteractionFactors[i][Index] += Conj(OneBodyBasis[Index3][0][1]) * OneBodyBasis[Index6][0][1]
                                                               * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index5][0][0]
                                                               * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]
                                                               * this->ComputeThreeBodyMatrixElementBAA(kx2, ky2, kx1, ky1, kx5, ky5, kx4, ky4); 
                  
                  this->ThreeBodyInteractionFactors[i][Index] *= + 2.0 * FactorW; // FIXME the plus sign is a work around

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
// kx1 = annihilation momentum along x for the B site
// ky1 = annihilation momentum along y for the B site
// kx2 = creation momentum along x for the B site
// ky2 = creation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementAB(int kx1, int ky1, int kx2, int ky2)
{
  double dx=((double)(kx1-kx2))/this->NbrSiteX;
  double dy=((double)(ky1-ky2))/this->NbrSiteY;
  Complex Tmp = 1 + Phase(2.0 * M_PI * (dx + dy)) + Phase(2.0 * M_PI * dy);
  return Tmp;
}

// compute the matrix element for the two body interaction between two A sites (or two B sites) 
//
// kx1 = annihilation momentum along x for the second site
// ky1 = annihilation momentum along y for the second site
// kx2 = creation momentum along x for the second site
// ky2 = creation momentum along y for the second site
// return value = corresponding matrix element

Complex ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian::ComputeTwoBodyMatrixElementAA(int kx1, int ky1, int kx2, int ky2)
{
  double dx=((double)(kx1-kx2))/this->NbrSiteX;
  double dy=((double)(ky1-ky2))/this->NbrSiteY;
  Complex Tmp = Phase(2.0 * M_PI * dx) + Phase(2.0 * M_PI * dy) + Phase(2.0 * M_PI * (+ dx + dy));
  return Tmp;
}

// compute the matrix element for the three body interaction between two sites B and one site A
//
// kx1 = annihilation momentum along x for the first B site
// ky1 = annihilation momentum along y for the first B site
// kx2 = annihilation momentum along x for the second B site
// ky2 = annihilation momentum along y for the second B site
// kx3 = creation momentum along x for the first B site
// ky3 = creation momentum along y for the first B site
// kx4 = creation momentum along x for the second B site
// ky4 = creation momentum along y for the second B site
// return value = corresponding matrix element

Complex ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementABB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  double dx13=((double)(kx1-kx3))/this->NbrSiteX;
  double dy13=((double)(ky1-ky3))/this->NbrSiteY;
  double dx24=((double)(kx2-kx4))/this->NbrSiteX;
  double dy24=((double)(ky2-ky4))/this->NbrSiteY;
  Complex Tmp = Phase(2.0 * M_PI * dy13) + Phase(2.0 * M_PI * (dy13 + dx24 + dy24)) + Phase(2.0 * M_PI * (dx13 + dy13));
  return Tmp;
}

// compute the matrix element for the three body interaction between one site B and two sites A
//
// kx1 = annihilation momentum along x for the first A site
// ky1 = annihilation momentum along y for the first A site
// kx2 = annihilation momentum along x for the second A site
// ky2 = annihilation momentum along y for the second A site
// kx3 = creation momentum along x for the first A site
// ky3 = creation momentum along y for the first A site
// kx4 = creation momentum along x for the second A site
// ky4 = creation momentum along y for the second A site
// return value = corresponding matrix element

Complex ParticleOnLatticeHaldaneModelSingleBandThreeBodyHamiltonian::ComputeThreeBodyMatrixElementBAA(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  double dx13=((double)(kx1-kx3))/this->NbrSiteX;
  double dy13=((double)(ky1-ky3))/this->NbrSiteY;
  double dx24=((double)(kx2-kx4))/this->NbrSiteX;
  double dy24=((double)(ky2-ky4))/this->NbrSiteY;
  Complex Tmp = Phase(2.0 * M_PI * (- dx13 - dy13 - dy24)) + Phase(2.0 * M_PI * (- dx13 - dy13)) + Phase(2.0 * M_PI * (-dy13));
  return Tmp;
}
