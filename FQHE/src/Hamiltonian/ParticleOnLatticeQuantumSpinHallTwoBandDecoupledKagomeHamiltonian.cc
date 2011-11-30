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
//                      using two decoupled kagome models                     //
//                                                                            //
//                        last modification : 30/11/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian.h"
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
// vPotential = strength of the repulsive on site two body interaction between opposite spins
// wPotential = strength of the repulsive two body neareast neighbor interaction between opposite spins
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mus = sublattice staggered chemical potential 
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian::ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
																		 int nbrSiteY, double uPotential, double vPotential, double wPotential, double t1, double t2, double lambda1, double lambda2, double mus, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
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
  this->MuS = mus;
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
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->HermitianSymmetryFlag = true;
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

ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian::~ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  if (this->FlatBand == false)
    {
      this->OneBodyInteractionFactorsupup = new double [this->NbrSiteX * this->NbrSiteY];
      this->OneBodyInteractionFactorsdowndown = new double [this->NbrSiteX * this->NbrSiteY];
    }


  this->ComputeOneBodyMatrices(OneBodyBasis);

  this->NbrInterSectorSums = this->NbrSiteX * this->NbrSiteY;
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
	  ++this->NbrInterSectorIndicesPerSum[(((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)];    
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
	  {
	    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = (kx1 * this->NbrSiteY) + ky1;
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = (kx2 * this->NbrSiteY) + ky2;
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	  }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = (kx1 * this->NbrSiteY) + ky1;
		int Index2 = (kx2 * this->NbrSiteY) + ky2;
		if (Index1 < Index2)
		  ++this->NbrIntraSectorIndicesPerSum[(((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)];    
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
	      {
		int Index1 = (kx1 * this->NbrSiteY) + ky1;
		int Index2 = (kx2 * this->NbrSiteY) + ky2;
		if (Index1 < Index2)
		  {
		    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }

      double FactorAUpBUp = -2.0 * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorAUpCUp = FactorAUpBUp;
      double FactorBUpCUp = FactorAUpBUp;
      double FactorAUpADown = FactorAUpBUp * this->VPotential;
      double FactorBUpBDown = FactorAUpBUp * this->VPotential;
      double FactorCUpCDown = FactorAUpBUp * this->VPotential;
      double FactorAUpBDown = FactorAUpBUp * this->WPotential;
      if (this->FlatBand == false)
	FactorAUpBUp *= this->UPotential;
      this->InteractionFactorsupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndown = new Complex* [this->NbrIntraSectorSums];

      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndown[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->IntraSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->IntraSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteY;
		  int ky3 = Index3 % this->NbrSiteY;
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;

// upup upup coefficient
 		  this->InteractionFactorsupup[i][Index] = FactorAUpBUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpBUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpBUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsupup[i][Index] += FactorAUpBUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx3, ky3);

 		  this->InteractionFactorsupup[i][Index] += FactorAUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpCUp(kx2, ky2, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpCUp(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpCUp(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsupup[i][Index] += FactorAUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpCUp(kx1, ky1, kx3, ky3);

 		  this->InteractionFactorsupup[i][Index] += FactorBUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 1, 2, 1, 2) * this->ComputeTwoBodyMatrixElementBUpCUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorBUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 1, 2, 1, 2) * this->ComputeTwoBodyMatrixElementBUpCUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorBUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 1, 2, 1, 2) * this->ComputeTwoBodyMatrixElementBUpCUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactorsupup[i][Index] += FactorBUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 1, 2, 1, 2) * this->ComputeTwoBodyMatrixElementBUpCUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

 		  this->InteractionFactorsupup[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 0, 3, 0, 3);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 0, 3, 0, 3);
 		  this->InteractionFactorsupup[i][Index] -= FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 0, 3, 0, 3);
 		  this->InteractionFactorsupup[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 0, 3, 0, 3);

 		  this->InteractionFactorsupup[i][Index] += FactorBUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorBUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorBUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactorsupup[i][Index] += FactorBUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

 		  this->InteractionFactorsupup[i][Index] += FactorCUpCDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 0, 0, 0, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementCUpCDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorCUpCDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 0, 0, 0, 0, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementCUpCDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactorsupup[i][Index] -= FactorCUpCDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 0, 0, 0, 0, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementCUpCDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactorsupup[i][Index] += FactorCUpCDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 0, 0, 0, 0, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementCUpCDown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);



// downdown downdown coefficient

 		  this->InteractionFactorsdowndown[i][Index] = FactorAUpBUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpBUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpBUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsdowndown[i][Index] += FactorAUpBUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx3, ky3);

 		  this->InteractionFactorsdowndown[i][Index] += FactorAUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpCUp(kx2, ky2, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpCUp(kx1, ky1, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpCUp(kx2, ky2, kx3, ky3);
 		  this->InteractionFactorsdowndown[i][Index] += FactorAUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 2, 0, 2) * this->ComputeTwoBodyMatrixElementAUpCUp(kx1, ky1, kx3, ky3);

 		  this->InteractionFactorsdowndown[i][Index] += FactorBUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 1, 2, 1, 2) * this->ComputeTwoBodyMatrixElementBUpCUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorBUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 1, 2, 1, 2) * this->ComputeTwoBodyMatrixElementBUpCUp(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorBUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 1, 2, 1, 2) * this->ComputeTwoBodyMatrixElementBUpCUp(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactorsdowndown[i][Index] += FactorBUpCUp * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 1, 2, 1, 2) * this->ComputeTwoBodyMatrixElementBUpCUp(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

 		  this->InteractionFactorsdowndown[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 0, 3, 0, 3);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 0, 3, 0, 3);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 0, 3, 0, 3);
 		  this->InteractionFactorsdowndown[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 0, 3, 0, 3);

 		  this->InteractionFactorsdowndown[i][Index] += FactorBUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorBUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorBUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactorsdowndown[i][Index] += FactorBUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

 		  this->InteractionFactorsdowndown[i][Index] += FactorCUpCDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 1, 1, 1, 1, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementCUpCDown(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorCUpCDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index3, Index4, 1, 1, 1, 1, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementCUpCDown(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  this->InteractionFactorsdowndown[i][Index] -= FactorCUpCDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index4, Index3, 1, 1, 1, 1, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementCUpCDown(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  this->InteractionFactorsdowndown[i][Index] += FactorCUpCDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index2, Index1, Index4, Index3, 1, 1, 1, 1, 1, 4, 1, 4) * this->ComputeTwoBodyMatrixElementCUpCDown(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);

		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}	  
      
      this->InteractionFactorsupdown = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrInterSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->InterSectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->InterSectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteY;
		  int ky3 = Index3 % this->NbrSiteY;
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;

		  this->InteractionFactorsupdown[i][Index] = 0.0;

//  		  this->InteractionFactorsupdown[i][Index] = Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 1, 0, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);


//  		  this->InteractionFactorsupdown[i][Index] += Factor * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 2, 3, 2, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);

//  		  this->InteractionFactorsupdown[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 2, 0, 2);

//  		  this->InteractionFactorsupdown[i][Index] += FactorAUpADown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 1, 3, 1, 3) * this->ComputeTwoBodyMatrixElementBUpBDown(kx2, ky2, kx4, ky4) * this->ComputeTwoBodyMatrixElementBUpBDown(kx1, ky1, kx3, ky3);

//  		  this->InteractionFactorsupdown[i][Index] += FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 0, 3, 0, 3) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);


//  		  this->InteractionFactorsupdown[i][Index] += FactorAUpBDown * this->ComputeTransfomationBasisContribution(OneBodyBasis, Index1, Index2, Index3, Index4, 0, 1, 0, 1, 2, 1, 2, 1) * this->ComputeTwoBodyMatrixElementAUpBUp(kx2, ky2, kx4, ky4);

		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }

  delete[] OneBodyBasis;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the matrix element for the two body interaction between two sites A and B belonging to the same up layer
//
// kx1 = creation momentum along x for the B site
// ky1 = creation momentum along y for the B site
// kx2 = annihilation momentum along x for the B site
// ky2 = annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian::ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kx2, int ky2)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KxFactor * ((double) (kx2 - kx1))));
  return Tmp;
}


// compute the matrix element for the two body interaction between two sites A and C belonging to the same up layer
//
// kx1 = creation momentum along x for the C site
// ky1 = creation momentum along y for the C site
// kx2 = annihilation momentum along x for the C site
// ky2 = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian::ComputeTwoBodyMatrixElementAUpCUp(int kx1, int ky1, int kx2, int ky2)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KyFactor * ((double) (ky2 - ky1))));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B and C belonging to the same up layer
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

Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian::ComputeTwoBodyMatrixElementBUpCUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 2.0 * cos (0.5 * ((this->KxFactor * ((double) (kx3 - kx1))) + (this->KyFactor * ((double) (ky4 - ky2)))));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and B belonging to the same down layer
//
// kx1 = creation momentum along x for the B site
// ky1 = creation momentum along y for the B site
// kx2 = annihilation momentum along x for the B site
// ky2 = annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian::ComputeTwoBodyMatrixElementADownBDown(int kx1, int ky1, int kx2, int ky2)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kx2, ky2);
}

// compute the matrix element for the two body interaction between two sites A and C belonging to the same down layer
//
// kx1 = creation momentum along x for the C site
// ky1 = creation momentum along y for the C site
// kx2 = annihilation momentum along x for the C site
// ky2 = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian::ComputeTwoBodyMatrixElementADownCDown(int kx1, int ky1, int kx2, int ky2)
{
  return this->ComputeTwoBodyMatrixElementAUpCUp(kx1, ky1, kx2, ky2);
}

// compute the matrix element for the two body interaction between two sites B and C belonging to the same down layer
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

Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian::ComputeTwoBodyMatrixElementBDownCDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return this->ComputeTwoBodyMatrixElementBUpCUp(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
}

// compute the matrix element for the two body interaction between two sites A with different layer indices 
//
// kx1 = first creation momentum along x for the A site
// ky1 = first creation momentum along y for the A site
// kx2 = second creation momentum along x for the A site
// ky2 = second creation momentum along y for the A site
// kx3 = first annihilation momentum along x for the A site
// ky3 = first annihilation momentum along y for the A site
// kx4 = second annihilation momentum along x for the A site
// ky4 = second annihilation momentum along y for the A site

Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian::ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 0.0;
  return Tmp;
}
    
// compute the matrix element for the two body interaction between two sites B with different layer indices 
//
// kx1 = first creation momentum along x for the B site
// ky1 = first creation momentum along y for the B site
// kx2 = second creation momentum along x for the B site
// ky2 = second creation momentum along y for the B site
// kx3 = first annihilation momentum along x for the B site
// ky3 = first annihilation momentum along y for the B site
// kx4 = second annihilation momentum along x for the B site
// ky4 = second annihilation momentum along y for the B site

Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian::ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 0.0;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites C with different layer indices 
//
// kx1 = first creation momentum along x for the C site
// ky1 = first creation momentum along y for the C site
// kx2 = second creation momentum along x for the C site
// ky2 = second creation momentum along y for the C site
// kx3 = first annihilation momentum along x for the C site
// ky3 = first annihilation momentum along y for the C site
// kx4 = second annihilation momentum along x for the C site
// ky4 = second annihilation momentum along y for the C site

Complex ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian::ComputeTwoBodyMatrixElementCUpCDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 0.0;
  return Tmp;
}
    
// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = array of one body transformation matrices (the leftmost upper block for the spin up, the rightmsdt lower block for the spin down)

void ParticleOnLatticeQuantumSpinHallTwoBandDecoupledKagomeHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
  double KX, KY;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = (kx * this->NbrSiteY) + ky;

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
	    oneBodyBasis[Index][i][j] = TmpMatrix[i][j];
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
	    oneBodyBasis[Index][3 + i][3 + j] = TmpMatrix[i][j];
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactorsdowndown[Index] = 0.5 * TmpDiag(0, 0);
	  }
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2) << endl;
      }
}
