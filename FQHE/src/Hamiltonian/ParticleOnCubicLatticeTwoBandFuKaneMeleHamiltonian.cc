////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of 3d topological insulator based on the Fu-Kane-Mele     //
//                       model and restricted to two bands                    //
//                                                                            //
//                        last modification : 18/07/2011                      //
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
#include "Hamiltonian/ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian.h"
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
// nbrSiteZ = number of sites in the z direction
// uPotential = strength of the repulsive two body neareast neighbor interaction
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// t2p = hoping amplitude between second next neareast neighbor sites
// mixingTermNorm = norm of the mixing term coupling the two copies of the checkerboard lattice
// mixingTermArgv = argument of the mixing term coupling the two copies of the checkerboard lattice
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
												       int nbrSiteY, int nbrSiteZ, double uPotential, double t1, double t2, double t2p, double mixingTermNorm, double mixingTermArg, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteYZ = this->NbrSiteY * this->NbrSiteZ;
  this->LzMax = nbrSiteX * nbrSiteY * nbrSiteZ - 1;
  this->HamiltonianShift = 0.0;
  this->NNHoping = t1;
  this->NextNNHoping = t2;
  this->SecondNextNNHoping = t2p;
  this->MixingTerm = mixingTermNorm * Phase(mixingTermArg);
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
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

ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::~ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ];
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  if (this->FlatBand == false)
    {
      this->OneBodyInteractionFactorsupup = new double [this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ];
      this->OneBodyInteractionFactorsdowndown = new double [this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ];
    }

  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      for (int kz = 0; kz < this->NbrSiteZ; ++kz)
	{
	  HermitianMatrix TmpOneBodyHamiltonian(4, true);
	  int Index = (kx * this->NbrSiteY) + ky;
	  Complex B1 = 4.0 * this->NNHoping * Complex (cos (1.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (1.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY)) * cos(M_PI * 0.25), 
						       sin (1.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * sin (1.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY)) * sin(M_PI * 0.25));
	  double d1 = 4.0 * this->SecondNextNNHoping * cos (2.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (2.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY));
	  double d3 = 2.0 * this->NextNNHoping * (cos (2.0 * M_PI * (((double) kx) + this->GammaX) / ((double) this->NbrSiteX))
						  - cos (2.0 * M_PI * (((double) ky) + this->GammaY) / ((double) this->NbrSiteY)));
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d1 + d3);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B1);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 1, d1 - d3);
	  B1 = 4.0 * this->NNHoping * Complex (cos (1.0 * M_PI * (((double) -kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (1.0 * M_PI * (((double) -ky) + this->GammaY) / ((double) this->NbrSiteY)) * cos(M_PI * 0.25), 
					       sin (1.0 * M_PI * (((double) -kx) + this->GammaX) / ((double) this->NbrSiteX)) * sin (1.0 * M_PI * (((double) -ky) + this->GammaY) / ((double) this->NbrSiteY)) * sin(M_PI * 0.25));
	  d1 = 4.0 * this->SecondNextNNHoping * cos (2.0 * M_PI * (((double) -kx) + this->GammaX) / ((double) this->NbrSiteX)) * cos (2.0 * M_PI * (((double) -ky) + this->GammaY) / ((double) this->NbrSiteY));
	  d3 = 2.0 * this->NextNNHoping * (cos (2.0 * M_PI * (((double) -kx) + this->GammaX) / ((double) this->NbrSiteX))
					   - cos (2.0 * M_PI * (((double) -ky) + this->GammaY) / ((double) this->NbrSiteY)));
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 2, d1 + d3);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 3, Conj(B1));
	  TmpOneBodyHamiltonian.SetMatrixElement(3, 3, d1 - d3);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 3, - I() * this->MixingTerm);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 2, I() * this->MixingTerm);
	  //	TmpOneBodyHamiltonian.SetMatrixElement(1, 2, - I() * this->MixingTerm);
	  //	TmpOneBodyHamiltonian.SetMatrixElement(0, 3, I() * this->MixingTerm);
	  ComplexMatrix TmpMatrix(4, 4, true);
	  TmpMatrix[0][0] = 1.0;
	  TmpMatrix[1][1] = 1.0;
	  TmpMatrix[2][2] = 1.0;
	  TmpMatrix[3][3] = 1.0;
	  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	  OneBodyBasis[Index] = TmpMatrix;	
	  if (this->FlatBand == false)
	    {
	      this->OneBodyInteractionFactorsupup[Index] = TmpDiag(0, 0);
	      this->OneBodyInteractionFactorsdowndown[Index] = TmpDiag(1, 1);
	    }
	  cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2) << " " << TmpDiag(3, 3) << endl;
	}

 
  this->NbrInterSectorSums = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  this->NbrIntraSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    this->NbrIntraSectorIndicesPerSum[i] = 0;      

  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)      
	  for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
	    for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
	      ++this->NbrInterSectorIndicesPerSum[((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)) * this->NbrSiteZ + ((kz1 + kz2) % this->NbrSiteZ)];    
  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      if (this->NbrInterSectorIndicesPerSum[i] > 0)
	{
	  this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
	  this->NbrInterSectorIndicesPerSum[i] = 0;
	}
    }
  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)    
	  for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
	    for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
	      {
		int TmpSum = ((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)) * this->NbrSiteZ + ((kz1 + kz2) % this->NbrSiteZ);
		this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = ((kx1 * this->NbrSiteY) + ky1) * this->NbrSiteZ + kz1;
		this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = ((kx2 * this->NbrSiteY) + ky2) * this->NbrSiteZ + kz2;
		++this->NbrInterSectorIndicesPerSum[TmpSum];    
	      }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
		for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
		  {
		    int Index1 = ((kx1 * this->NbrSiteY) + ky1) * this->NbrSiteZ + kz1;
		    int Index2 = ((kx2 * this->NbrSiteY) + ky2) * this->NbrSiteZ + kz2;
		    if (Index1 < Index2)
		      ++this->NbrIntraSectorIndicesPerSum[((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)) * this->NbrSiteZ + ((kz1 + kz2) % this->NbrSiteZ)];    
		  }
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	    {
	      this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	      this->NbrIntraSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      for (int kz1 = 0; kz1 < this->NbrSiteZ; ++kz1)
		for (int kz2 = 0; kz2 < this->NbrSiteZ; ++kz2)      
		  {
		    int Index1 = ((kx1 * this->NbrSiteY) + ky1) * this->NbrSiteZ + kz1;
		    int Index2 = ((kx2 * this->NbrSiteY) + ky2) * this->NbrSiteZ + kz2;
		    if (Index1 < Index2)
		      {
			int TmpSum = ((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)) * this->NbrSiteZ + ((kz1 + kz2) % this->NbrSiteZ);
			this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
			this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
			++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		      }
		  }
      
      double Factor = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      if (this->FlatBand == false)
	Factor *= this->UPotential;

      this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupupupdown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupupdowndown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndownupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndowndowndown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndownupdown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupdownupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupdownupdown = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsupdowndowndown = new Complex* [this->NbrIntraSectorSums];

      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndowndowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsupupdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndownupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsupdowndowndown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsupdownupup[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteYZ;
	      int ky1 = Index1 % this->NbrSiteYZ;
	      int kz1 = ky1 % this->NbrSiteZ;
	      ky1 /= this->NbrSiteZ;
	      int kx2 = Index2 / this->NbrSiteYZ;
	      int ky2 = Index2 % this->NbrSiteYZ;
	      int kz2 = ky2 % this->NbrSiteZ;
	      ky2 /= this->NbrSiteZ;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteYZ;
		  int ky3 = Index3 % this->NbrSiteYZ;
		  int kz3 = ky3 % this->NbrSiteZ;
		  ky3 /= this->NbrSiteZ;
		  int kx4 = Index4 / this->NbrSiteYZ;
		  int ky4 = Index4 % this->NbrSiteYZ;
		  int kz4 = ky4 % this->NbrSiteZ;
		  ky4 /= this->NbrSiteZ;

		  //  upup upup coefficient
 		  this->InteractionFactorsupupupup[i][Index] = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx4, ky4, kz4);
 		  this->InteractionFactorsupupupup[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx4, ky4, kz4);
 		  this->InteractionFactorsupupupup[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx3, ky3, kz3);
 		  this->InteractionFactorsupupupup[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx3, ky3, kz3);

 		  this->InteractionFactorsupupupup[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx4, ky4, kz4);
 		  this->InteractionFactorsupupupup[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx4, ky4, kz4);
 		  this->InteractionFactorsupupupup[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx3, ky3, kz3);
 		  this->InteractionFactorsupupupup[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx3, ky3, kz3);

 		  this->InteractionFactorsupupupup[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx3, ky3, kz3);
 		  this->InteractionFactorsupupupup[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx3, ky3, kz3);;
 		  this->InteractionFactorsupupupup[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx3, ky3, kz3) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx4, ky4, kz4);
 		  this->InteractionFactorsupupupup[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx3, ky3, kz3) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx4, ky4, kz4);

 		  this->InteractionFactorsupupupup[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx3, ky3, kz3);
 		  this->InteractionFactorsupupupup[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx3, ky3, kz3);
 		  this->InteractionFactorsupupupup[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx3, ky3, kz3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx4, ky4, kz4);
 		  this->InteractionFactorsupupupup[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx3, ky3, kz3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx4, ky4, kz4);

		  this->InteractionFactorsupupupup[i][Index] *= -2.0 * Factor;

		  // downdown downdown coefficient
 		  this->InteractionFactorsdowndowndowndown[i][Index] = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx4, ky4, kz4);
 		  this->InteractionFactorsdowndowndowndown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx4, ky4, kz4);
 		  this->InteractionFactorsdowndowndowndown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx3, ky3, kz3);
 		  this->InteractionFactorsdowndowndowndown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx3, ky3, kz3);

 		  this->InteractionFactorsdowndowndowndown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx4, ky4, kz4);
 		  this->InteractionFactorsdowndowndowndown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx4, ky4, kz4);
 		  this->InteractionFactorsdowndowndowndown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx3, ky3, kz3);
 		  this->InteractionFactorsdowndowndowndown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx3, ky3, kz3);

 		  this->InteractionFactorsdowndowndowndown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx3, ky3, kz3);
 		  this->InteractionFactorsdowndowndowndown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx3, ky3, kz3);
 		  this->InteractionFactorsdowndowndowndown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx3, ky3, kz3) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx4, ky4, kz4);
 		  this->InteractionFactorsdowndowndowndown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx3, ky3, kz3) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx4, ky4, kz4);

 		  this->InteractionFactorsdowndowndowndown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx3, ky3, kz3);
 		  this->InteractionFactorsdowndowndowndown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx3, ky3, kz3);
 		  this->InteractionFactorsdowndowndowndown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx3, ky3, kz3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx4, ky4, kz4);
 		  this->InteractionFactorsdowndowndowndown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx3, ky3, kz3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx4, ky4, kz4);

		  this->InteractionFactorsdowndowndowndown[i][Index] *= -2.0 * Factor;



		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	  Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteYZ;
	      int ky1 = Index1 % this->NbrSiteYZ;
	      int kz1 = ky1 % this->NbrSiteZ;
	      ky1 /= this->NbrSiteZ;
	      int kx2 = Index2 / this->NbrSiteYZ;
	      int ky2 = Index2 % this->NbrSiteYZ;
	      int kz2 = ky2 % this->NbrSiteZ;
	      ky2 /= this->NbrSiteZ;
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteYZ;
		  int ky3 = Index3 % this->NbrSiteYZ;
		  int kz3 = ky3 % this->NbrSiteZ;
		  ky3 /= this->NbrSiteZ;
		  int kx4 = Index4 / this->NbrSiteYZ;
		  int ky4 = Index4 % this->NbrSiteYZ;
		  int kz4 = ky4 % this->NbrSiteZ;
		  ky4 /= this->NbrSiteZ;

		  // updown upup coefficient
		  this->InteractionFactorsupdownupup[i][Index] = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx4, ky4, kz4);
		  this->InteractionFactorsupdownupup[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx3, ky3, kz3);
		  this->InteractionFactorsupdownupup[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx4, ky4, kz4);
		  this->InteractionFactorsupdownupup[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx3, ky3, kz3);
		  this->InteractionFactorsupdownupup[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx3, ky3, kz3);
		  this->InteractionFactorsupdownupup[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx3, ky3, kz3) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx4, ky4, kz4);
		  this->InteractionFactorsupdownupup[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx3, ky3, kz3);
		  this->InteractionFactorsupdownupup[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 0, 0, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx3, ky3, kz3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx4, ky4, kz4);
		  this->InteractionFactorsupdownupup[i][Index] *= -2.0 * Factor;

		  // updown downdown coefficient
		  this->InteractionFactorsupdowndowndown[i][Index] = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx4, ky4, kz4);
		  this->InteractionFactorsupdowndowndown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx3, ky3, kz3);
		  this->InteractionFactorsupdowndowndown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx4, ky4, kz4);
		  this->InteractionFactorsupdowndowndown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx3, ky3, kz3);
		  this->InteractionFactorsupdowndowndown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx3, ky3, kz3);
		  this->InteractionFactorsupdowndowndown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx3, ky3, kz3) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx4, ky4, kz4);
		  this->InteractionFactorsupdowndowndown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx3, ky3, kz3);
		  this->InteractionFactorsupdowndowndown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 1, 1, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx3, ky3, kz3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx4, ky4, kz4);
		  this->InteractionFactorsupdowndowndown[i][Index] *= -2.0 * Factor;

		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}

      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdownupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  if (this->NbrIntraSectorIndicesPerSum[i] > 0)
	    {
	      this->InteractionFactorsupupupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	      this->InteractionFactorsdowndownupdown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	    }
	  else
	    {
	      this->InteractionFactorsupupupdown[i] = 0;
	      this->InteractionFactorsdowndownupdown[i] = 0;
	    }
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteYZ;
	      int ky1 = Index1 % this->NbrSiteYZ;
	      int kz1 = ky1 % this->NbrSiteZ;
	      ky1 /= this->NbrSiteZ;
	      int kx2 = Index2 / this->NbrSiteYZ;
	      int ky2 = Index2 % this->NbrSiteYZ;
	      int kz2 = ky2 % this->NbrSiteZ;
	      ky2 /= this->NbrSiteZ;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteYZ;
		  int ky3 = Index3 % this->NbrSiteYZ;
		  int kz3 = ky3 % this->NbrSiteZ;
		  ky3 /= this->NbrSiteZ;
		  int kx4 = Index4 / this->NbrSiteYZ;
		  int ky4 = Index4 % this->NbrSiteYZ;
		  int kz4 = ky4 % this->NbrSiteZ;
		  ky4 /= this->NbrSiteZ;

		  // upup updown coefficient
		  this->InteractionFactorsupupupdown[i][Index] = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx4, ky4, kz4);
		  this->InteractionFactorsupupupdown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx4, ky4, kz4);
		  this->InteractionFactorsupupupdown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx4, ky4, kz4);
		  this->InteractionFactorsupupupdown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx4, ky4, kz4);
		  this->InteractionFactorsupupupdown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx3, ky3, kz3);
		  this->InteractionFactorsupupupdown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx3, ky3, kz3);
		  this->InteractionFactorsupupupdown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx3, ky3, kz3);
		  this->InteractionFactorsupupupdown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx3, ky3, kz3);
		  this->InteractionFactorsupupupdown[i][Index] *= -2.0 * Factor;

		  // downdown updown coefficient
		  this->InteractionFactorsdowndownupdown[i][Index] = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx4, ky4, kz4);
		  this->InteractionFactorsdowndownupdown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx4, ky4, kz4);
		  this->InteractionFactorsdowndownupdown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx3, ky3, kz3);
		  this->InteractionFactorsdowndownupdown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx3, ky3, kz3);
		  this->InteractionFactorsdowndownupdown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 0, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx3, ky3, kz3);
		  this->InteractionFactorsdowndownupdown[i][Index] -= this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 0, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx3, ky3, kz3);
		  this->InteractionFactorsdowndownupdown[i][Index] *= -2.0 * Factor;


		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	  Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteYZ;
	      int ky1 = Index1 % this->NbrSiteYZ;
	      int kz1 = ky1 % this->NbrSiteZ;
	      ky1 /= this->NbrSiteZ;
	      int kx2 = Index2 / this->NbrSiteYZ;
	      int ky2 = Index2 % this->NbrSiteYZ;
	      int kz2 = ky2 % this->NbrSiteZ;
	      ky2 /= this->NbrSiteZ;
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteYZ;
		  int ky3 = Index3 % this->NbrSiteYZ;
		  int kz3 = ky3 % this->NbrSiteZ;
		  ky3 /= this->NbrSiteZ;
		  int kx4 = Index4 / this->NbrSiteYZ;
		  int ky4 = Index4 % this->NbrSiteYZ;
		  int kz4 = ky4 % this->NbrSiteZ;
		  ky4 /= this->NbrSiteZ;
 		  this->InteractionFactorsupdownupdown[i][Index] = this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx4, ky4, kz4);
 		  this->InteractionFactorsupdownupdown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kz2, kx4, ky4, kz4);
 		  this->InteractionFactorsupdownupdown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx3, ky3, kz3);
 		  this->InteractionFactorsupdownupdown[i][Index] += this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kz2, kx4, ky4, kz4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kz1, kx3, ky3, kz3);
		  this->InteractionFactorsupdownupdown[i][Index] *= -2.0 * Factor;
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}

    }

  delete[] OneBodyBasis;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the matrix element for the two body interaction between two sites A and B belonging to the same layer
//
// kx1 = momentum along x for the A site
// ky1 = momentum along y for the A site
// kx2 = momentum along x for the B site
// ky2 = momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2)
{
  Complex Tmp = 2.0 * (cos (M_PI * ((((double) (kx2 - kx1)) / ((double) this->NbrSiteX)) - ((((double) (ky2 - ky1)) / ((double) this->NbrSiteY))))) 
		       + cos (M_PI * ((((double) (kx2 - kx1)) / ((double) this->NbrSiteX)) + ((((double) (ky2 - ky1)) / ((double) this->NbrSiteY))))));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A with different layer indices 
//
// kx1 = momentum along x for the first A site
// ky1 = momentum along y for the first A site
// kx2 = momentum along x for the second A site
// ky2 = momentum along y for the second A site
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2)
{
  Complex Tmp = 1.0;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B with different layer indices 
//
// kx1 = momentum along x for the first B site
// ky1 = momentum along y for the first B site
// kx2 = momentum along x for the second B site
// ky2 = momentum along y for the second B site
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeTwoBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2)
{
  Complex Tmp = Phase (M_PI * ((((double) (kx2 - kx1)) / ((double) this->NbrSiteX)) + ((((double) (ky2 - ky1)) / ((double) this->NbrSiteY)))));
  return Tmp;
}
