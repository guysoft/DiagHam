////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//            class of 3d topological insulator based on the pryochlore       //
//                       model and restricted to four bands                   //
//                                                                            //
//                        last modification : 13/08/2012                      //
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
#include "Hamiltonian/ParticleOnCubicLatticeFourBandPyrochloreHamiltonian.h"
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

ParticleOnCubicLatticeFourBandPyrochloreHamiltonian::ParticleOnCubicLatticeFourBandPyrochloreHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// nbrSiteZ = number of sites in the z direction
// uPotential = strength of the repulsive two body on site interactions
// vPotential = trength of the repulsive two body neareast neighbor interaction
// tightBindingModel = pointer to the tight binding model
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnCubicLatticeFourBandPyrochloreHamiltonian::ParticleOnCubicLatticeFourBandPyrochloreHamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int nbrSiteX, 
													 int nbrSiteY, int nbrSiteZ, double uPotential, double vPotential, 
													 Abstract3DTightBindingModel* tightBindingModel, 
													 bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteYZ = this->NbrSiteY * this->NbrSiteZ;
  this->LzMax = nbrSiteX * nbrSiteY * nbrSiteZ - 1;
  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;
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
  this->HermitianSymmetryFlag = true;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();

  int Dim = this->Particles->GetHilbertSpaceDimension();
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

ParticleOnCubicLatticeFourBandPyrochloreHamiltonian::~ParticleOnCubicLatticeFourBandPyrochloreHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnCubicLatticeFourBandPyrochloreHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  int NbrSites = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;

  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;

  this->InteractionFactorsupupupup = 0;
  this->InteractionFactorsupumupum = 0;
  this->InteractionFactorsupdpupdp = 0;
  this->InteractionFactorsupdmupdm = 0;
  this->InteractionFactorsumumumum = 0;
  this->InteractionFactorsumdpumdp = 0;
  this->InteractionFactorsumdmumdm = 0;
  this->InteractionFactorsdpdpdpdp = 0;
  this->InteractionFactorsdpdmdpdm = 0;
  this->InteractionFactorsdmdmdmdm = 0;

  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsumum = 0;
  this->OneBodyInteractionFactorsdpdp = 0;
  this->OneBodyInteractionFactorsdmdm = 0;
  this->OneBodyInteractionFactorsupum = 0;
  this->OneBodyInteractionFactorsupdp = 0;
  this->OneBodyInteractionFactorsupdm = 0;
  this->OneBodyInteractionFactorsumdp = 0;
  this->OneBodyInteractionFactorsumdm = 0;
  this->OneBodyInteractionFactorsdpdm = 0;

  ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];
  if (this->FlatBand == false)
    {
      this->OneBodyInteractionFactorsupup = new double [NbrSites];
      this->OneBodyInteractionFactorsumum = new double [NbrSites];
      this->OneBodyInteractionFactorsdpdp = new double [NbrSites];
      this->OneBodyInteractionFactorsdmdm = new double [NbrSites];
    }
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      for (int kz = 0; kz < this->NbrSiteZ; ++kz)
	{
	  int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky, kz);
	  if (this->FlatBand == false)
	    {
	      this->OneBodyInteractionFactorsupup[Index] = 0.5 * this->TightBindingModel->GetEnergy(0, Index);
	      this->OneBodyInteractionFactorsumum[Index] = 0.5 * this->TightBindingModel->GetEnergy(1, Index);
	      this->OneBodyInteractionFactorsdpdp[Index] = 0.5 * this->TightBindingModel->GetEnergy(2, Index);
	      this->OneBodyInteractionFactorsdmdm[Index] = 0.5 * this->TightBindingModel->GetEnergy(3, Index);
	    }
	  OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
	}


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
    }
  else
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
		    if (Index1 <= Index2)
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
		    if (Index1 <= Index2)
		      {
			int TmpSum = ((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY)) * this->NbrSiteZ + ((kz1 + kz2) % this->NbrSiteZ);
			this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
			this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
			++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		      }
		  }
      
      double Factor = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ));
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

