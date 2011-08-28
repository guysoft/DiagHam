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
//                       model and restricted to four bands                   //
//                                                                            //
//                        last modification : 26/08/2011                      //
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
#include "Hamiltonian/ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian.h"
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
// vPotential = strength of the repulsive two body on site interaction
// nnHopingDistortion111 = distortion of nearest neighbor hoping amplitude in the (111) direction
// spinOrbitCoupling = amplitude of the spin orbit coupling
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// gammaZ = boundary condition twisting angle along z
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian::ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int nbrSiteX, 
												       int nbrSiteY, int nbrSiteZ, double uPotential, double vPotential, double nnHopingDistortion111, double spinOrbitCoupling, double gammaX, double gammaY, double gammaZ, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteYZ = this->NbrSiteY * this->NbrSiteZ;
  this->LzMax = nbrSiteX * nbrSiteY * nbrSiteZ - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->KzFactor = 2.0 * M_PI / ((double) this->NbrSiteZ);
  this->HamiltonianShift = 0.0;
  this->NNHopingDistortion111 = nnHopingDistortion111;
  this->SpinOrbitCoupling = spinOrbitCoupling;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->GammaZ = gammaZ;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsupum = 0;
  this->OneBodyInteractionFactorsupdp = 0;
  this->OneBodyInteractionFactorsupdm = 0;
  this->OneBodyInteractionFactorsumum = 0;
  this->OneBodyInteractionFactorsumdp = 0;
  this->OneBodyInteractionFactorsumdm = 0;
  this->OneBodyInteractionFactorsdpdp = 0;
  this->OneBodyInteractionFactorsdpdm = 0;
  this->OneBodyInteractionFactorsdmdm = 0;
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

ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian::~ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  int NbrSites = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
  HermitianMatrix*OneBodyHamiltonian  = new HermitianMatrix [NbrSites];
  this->ComputeOneBodyHamiltonian(OneBodyHamiltonian);

  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  this->OneBodyInteractionFactorsupup = new double [NbrSites];
  this->OneBodyInteractionFactorsumum = new double [NbrSites];
  this->OneBodyInteractionFactorsdpdp = new double [NbrSites];
  this->OneBodyInteractionFactorsdmdm = new double [NbrSites];
  this->OneBodyInteractionFactorsupum = new Complex [NbrSites];
  this->OneBodyInteractionFactorsupdp = new Complex [NbrSites];
  this->OneBodyInteractionFactorsupdm = new Complex [NbrSites];
  this->OneBodyInteractionFactorsumdp = new Complex [NbrSites];
  this->OneBodyInteractionFactorsumdm = new Complex [NbrSites];
  this->OneBodyInteractionFactorsdpdm = new Complex [NbrSites];

  for (int i = 0; i < NbrSites; ++i)
    {
      double Tmp1;
      Complex Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(0, 0, Tmp1);      
      this->OneBodyInteractionFactorsupup[i] = Tmp1;
      OneBodyHamiltonian[i].GetMatrixElement(1, 1, Tmp1);      
      this->OneBodyInteractionFactorsumum[i] = Tmp1;
      OneBodyHamiltonian[i].GetMatrixElement(2, 2, Tmp1);      
      this->OneBodyInteractionFactorsdpdp[i] = Tmp1;
      OneBodyHamiltonian[i].GetMatrixElement(3, 3, Tmp1);      
      this->OneBodyInteractionFactorsdmdm[i] = Tmp1;

      OneBodyHamiltonian[i].GetMatrixElement(0, 1, Tmp2);      
      this->OneBodyInteractionFactorsupum[i] = Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(0, 2, Tmp2);      
      this->OneBodyInteractionFactorsupdp[i] = Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(0, 3, Tmp2);      
      this->OneBodyInteractionFactorsupdm[i] = Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(1, 2, Tmp2);      
      this->OneBodyInteractionFactorsumdp[i] = Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(1, 3, Tmp2);      
      this->OneBodyInteractionFactorsumdm[i] = Tmp2;
      OneBodyHamiltonian[i].GetMatrixElement(2, 3, Tmp2);      
      this->OneBodyInteractionFactorsdpdm[i] = Tmp2;
    }
  this->ComputeOneBodyHamiltonian(OneBodyHamiltonian);
 
  this->NbrInterSectorSums = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  this->NbrIntraSectorSums = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
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
      
      double Factor = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ));
      double FactorAUpADown = Factor * this->VPotential * 0.0;
      double FactorBUpBDown = Factor * this->VPotential * 0.0;
      if (this->FlatBand == false)
	Factor *= this->UPotential;
      Factor = 0.0;
      double FactorAUpBUp = Factor;
      double FactorADownBDown = Factor;
      double FactorAUpBDown = Factor;
      double FactorADownBUp = Factor;

      Complex Tmp;

      //  upup upup coefficient
      this->InteractionFactorsupupupup = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsumumumum = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdpdpdpdp = new Complex* [this->NbrIntraSectorSums];
      this->InteractionFactorsdmdmdmdm = new Complex* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupupupup[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsumumumum[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdpdpdpdp[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdmdmdmdm[i] = new Complex[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
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

                  this->InteractionFactorsupupupup[i][Index] = 0.0;
                  Tmp = this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp -=this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx4, ky4, kz4, kx3, ky3, kz3);
                  Tmp -= this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx3, ky3, kz3, kx4, ky4, kz4);
                  Tmp += this->ComputeTwoBodyMatrixElementAUpADown(kx2, ky2, kz2, kx1, ky1, kz1, kx4, ky4, kz4, kx3, ky3, kz3);
                  this->InteractionFactorsupupupup[i][Index] += -2.0 * FactorAUpADown * Tmp;
		  
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}



      //  updown updown coefficient
      this->InteractionFactorsupumupum = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsupdpupdp = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsupdmupdm = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsumdpumdp = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsumdmumdm = new Complex* [this->NbrInterSectorSums];
      this->InteractionFactorsdpdmdpdm = new Complex* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupumupum[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsupdpupdp[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsupdmupdm[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsumdpumdp[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsumdmumdm[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
	  this->InteractionFactorsdpdmdpdm[i] = new Complex[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
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
                  this->InteractionFactorsupumupum[i][Index] = 0.0;
                  Tmp = this->ComputeTwoBodyMatrixElementAUpADown(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
                  this->InteractionFactorsupumupum[i][Index] += -2.0 * FactorAUpADown * Tmp;

		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }

  delete[] OneBodyHamiltonian;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the matrix element for the two body interaction between two sites A and B with up spins
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kz1 = momentum along z for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on B site with spin up
// ky2 = momentum along y for the creation operator on B site with spin up
// kz2 = momentum along z for the creation operator on B site with spin up
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kz3 = momentum along z for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin up
// ky4 = momentum along y for the annihilation operator on B site with spin up
// kz4 = momentum along z for the annihilation operator on B site with spin up
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  Complex Tmp = 1.0 ;
  Tmp += Phase ((((double) (kx1 - kx3)) * this->KxFactor));
  Tmp += Phase ((((double) (kz1 - kz3)) * this->KzFactor));
  Tmp += Phase ((((double) (ky1 - ky3)) * this->KyFactor));
  Tmp *= Phase (0.25 * ((((double) (kx4 - kx2)) * this->KxFactor) + (((double) (ky4 - ky2)) * this->KyFactor)
			+ (((double) (kz4 - kz2)) * this->KzFactor)));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and B with down spins
//
// kx1 = momentum along x for the creation operator on A site with spin down
// ky1 = momentum along y for the creation operator on A site with spin down
// kz1 = momentum along z for the creation operator on A site with spin down
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kz2 = momentum along z for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin down
// ky3 = momentum along y for the annihilation operator on A site with spin down
// kz3 = momentum along z for the annihilation operator on A site with spin down
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// kz4 = momentum along z for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementADownBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
}

// compute the matrix element for the two body interaction between two sites A and B with opposite spins
//
// kx1 = momentum along x for the creation operator on A site with spin down
// ky1 = momentum along y for the creation operator on A site with spin down
// kz1 = momentum along z for the creation operator on A site with spin down
// kx2 = momentum along x for the creation operator on B site with spin up
// ky2 = momentum along y for the creation operator on B site with spin up
// kz2 = momentum along z for the creation operator on B site with spin up
// kx3 = momentum along x for the annihilation operator on A site with spin down
// ky3 = momentum along y for the annihilation operator on A site with spin down
// kz3 = momentum along z for the annihilation operator on A site with spin down
// kx4 = momentum along x for the annihilation operator on B site with spin up
// ky4 = momentum along y for the annihilation operator on B site with spin up
// kz4 = momentum along z for the annihilation operator on B site with spin up
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementADownBUp(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
}

// compute the matrix element for the two body interaction between two sites A and B with opposite spins
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kz1 = momentum along z for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kz2 = momentum along z for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kz3 = momentum along z for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// kz4 = momentum along z for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementAUpBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  return this->ComputeTwoBodyMatrixElementAUpBUp(kx1, ky1, kz1, kx2, ky2, kz2, kx3, ky3, kz3, kx4, ky4, kz4);
}

// compute the matrix element for the two body interaction between two sites A with opposite spins 
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kz1 = momentum along z for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on A site with spin down
// ky2 = momentum along y for the creation operator on A site with spin down
// kz2 = momentum along z for the creation operator on A site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kz3 = momentum along z for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on A site with spin down
// ky4 = momentum along y for the annihilation operator on A site with spin down
// kz4 = momentum along z for the annihilation operator on A site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  Complex Tmp = 1.0;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B with opposite spins 
//
// kx1 = momentum along x for the creation operator on B site with spin up
// ky1 = momentum along y for the creation operator on B site with spin up
// kz1 = momentum along z for the creation operator on B site with spin up
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kz2 = momentum along z for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on B site with spin up
// ky3 = momentum along y for the annihilation operator on B site with spin up
// kz3 = momentum along z for the annihilation operator on B site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// kz4 = momentum along z for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian::ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kz1, int kx2, int ky2, int kz2, int kx3, int ky3, int kz3, int kx4, int ky4, int kz4)
{
  Complex Tmp = Phase (0.25  * ((((double) (kx4 + kx3 - kx1 - kx2)) * this->KxFactor) + (((double) (ky4 + ky3 - ky1 - ky2)) * this->KyFactor) + (((double) (kz4 + kz3 - kz1 - kz2)) * this->KzFactor)));
  return Tmp;
}

// compute the one body hamiltonians related to the band stucture contribution
//
// oneBodyHamiltonians = array of one body hamiltonians

void ParticleOnCubicLatticeFourBandFuKaneMeleHamiltonian::ComputeOneBodyHamiltonian(HermitianMatrix* oneBodyHamiltonians)
{
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      for (int kz = 0; kz < this->NbrSiteZ; ++kz)
	{
	  HermitianMatrix TmpOneBodyHamiltonian(4, true);
	  int Index = ((kx * this->NbrSiteY) + ky) * this->NbrSiteZ + kz;
	  Complex B1 = 1.0 + this->NNHopingDistortion111 + Phase(((double) kx) * this->KxFactor) + Phase(((double) ky) * this->KyFactor) + Phase(((double) kz) * this->KzFactor);
	  Complex TmpPhase = Phase (-0.25 * ((((double) kx) * this->KxFactor) +(((double) ky) * this->KyFactor) + (((double) kz) * this->KzFactor)));
	  B1 *= TmpPhase;
	  double d3 = this->SpinOrbitCoupling * (sin ((((double) ky) * this->KyFactor))
						 - sin ((((double) kz) * this->KzFactor))
						 - sin ((((double) ky) * this->KyFactor) - (((double) kx) * this->KxFactor))
						 + sin ((((double) kz) * this->KzFactor) - (((double) kx) * this->KxFactor)));
	  double d4 = this->SpinOrbitCoupling * (sin ((((double) kz) * this->KzFactor))
						 - sin ((((double) kx) * this->KxFactor))
						 - sin ((((double) kz) * this->KzFactor) - (((double) ky) * this->KyFactor))
						 + sin ((((double) kx) * this->KxFactor) - (((double) ky) * this->KyFactor)));
	  double d5 = this->SpinOrbitCoupling * (sin ((((double) kx) * this->KxFactor))
						 - sin ((((double) ky) * this->KyFactor))
						 - sin ((((double) kx) * this->KxFactor) - (((double) kz) * this->KzFactor))
						 + sin ((((double) ky) * this->KyFactor) - (((double) kz) * this->KzFactor)));
	  Complex B2 = d3 - I() * d4;
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d5);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -d5);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 2, -d5);
	  TmpOneBodyHamiltonian.SetMatrixElement(3, 3, d5);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B1);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 3, B1);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 2, B2);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 3, -B2);

	  if (this->FlatBand == false)
	    {
	      oneBodyHamiltonians[Index] = TmpOneBodyHamiltonian;
	    }
	  else
	    {
	      oneBodyHamiltonians[Index].Copy(TmpOneBodyHamiltonian);
	      RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	      TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
	      TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif   	      
	      oneBodyHamiltonians[Index] /= fabs(TmpDiag(0, 0));
	    }
	}
}
