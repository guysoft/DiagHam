////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Cecile Repellin                       //
//                                                                            //
//            class of Moire superlattice with Coulomb interactions           //
//                         projected to four bands                            //
//                                                                            //
//                        last modification : 17/05/2018                      //
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
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Hamiltonian/ParticleOnSuperlatticeFourBandCoulombHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"
#include "GeneralTools/StringTools.h"


#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

ParticleOnSuperlatticeFourBandCoulombHamiltonian::ParticleOnSuperlatticeFourBandCoulombHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = strength of the repulsive two body on site interactions
// vPotential = trength of the repulsive two body neareast neighbor interaction
// kineticScale = global energy scale of the kinetic energy term (i.e t1 hopping term)
// mass = mass term of the simple TI model
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnSuperlatticeFourBandCoulombHamiltonian::ParticleOnSuperlatticeFourBandCoulombHamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int nbrSiteX, 
												       int nbrSiteY, int truncateQx, int truncateQy, double uPotential, double screeningDistance, double gammaX, double gammaY, Abstract2DTightBindingModel* tightBindingModel, 
												       bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;
  this->BandIndex = (this->TightBindingModel->GetNbrBands())/2 - 1;
  
  
  
  this->TruncateQx = truncateQx;
  this->TruncateQy = truncateQy;
  
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->ScreeningDistance = screeningDistance;
  
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
    
    double Factor = 0.5 / ((double) this->NbrSiteX * this->NbrSiteY);
}



// destructor
//

ParticleOnSuperlatticeFourBandCoulombHamiltonian::~ParticleOnSuperlatticeFourBandCoulombHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnSuperlatticeFourBandCoulombHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;    
  int NbrSites = this->TightBindingModel->GetNbrStatePerBand();
  if (this->FlatBand == false)
    {
      this->OneBodyInteractionFactorsupup = new double [NbrSites];
      this->OneBodyInteractionFactorsumum = new double [NbrSites];
      this->OneBodyInteractionFactorsdpdp = new double [NbrSites];
      this->OneBodyInteractionFactorsdmdm = new double [NbrSites];
    }
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  int IndexPlus = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
      int IndexMinus = this->TightBindingModel->GetLinearizedMomentumIndexSafe(-kx, -ky);
	  
      if (this->FlatBand == false)
	    {
	      this->OneBodyInteractionFactorsupup[IndexPlus] = this->TightBindingModel->GetEnergy(this->BandIndex, IndexPlus);
	      this->OneBodyInteractionFactorsumum[IndexPlus] = this->TightBindingModel->GetEnergy(this->BandIndex, IndexMinus);
	      this->OneBodyInteractionFactorsdpdp[IndexPlus] = this->TightBindingModel->GetEnergy(this->BandIndex, IndexPlus);
	      this->OneBodyInteractionFactorsdmdm[IndexPlus] = this->TightBindingModel->GetEnergy(this->BandIndex, IndexMinus);
	    }
	}


  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  
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
	  ++this->NbrInterSectorIndicesPerSum[((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY))];    
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
	    int TmpSum = ((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY));
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = ((kx1 * this->NbrSiteY) + ky1);
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = ((kx2 * this->NbrSiteY) + ky2);
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	  }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = ((kx1 * this->NbrSiteY) + ky1);
		int Index2 = ((kx2 * this->NbrSiteY) + ky2);
		if (Index1 < Index2)
		  ++this->NbrIntraSectorIndicesPerSum[((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY))];    
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
		int Index1 = ((kx1 * this->NbrSiteY) + ky1);
		int Index2 = ((kx2 * this->NbrSiteY) + ky2);
		if (Index1 < Index2)
		  {
		    int TmpSum = ((((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY));
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      
      
      ////
      
      
      double Factor = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      if (this->FlatBand == false)
            Factor *= this->UPotential* (2.0 / sqrt(3.0));

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
          
          this->InteractionFactorsupupupup[i][Index] = 0.0;
          this->InteractionFactorsumumumum[i][Index] = 0.0;
          this->InteractionFactorsdpdpdpdp[i][Index] = 0.0;
          this->InteractionFactorsdmdmdmdm[i][Index] = 0.0;          
          
          this->InteractionFactorsupupupup[i][Index] += Factor * this->EvaluateInteractionCoefficient(Index1, Index2, Index3, Index4, 1, 1);
          this->InteractionFactorsupupupup[i][Index] += Factor * this->EvaluateInteractionCoefficient(Index2, Index1, Index4, Index3, 1, 1);
          this->InteractionFactorsupupupup[i][Index] -= Factor * this->EvaluateInteractionCoefficient(Index2, Index1, Index3, Index4, 1, 1);
          this->InteractionFactorsupupupup[i][Index] -= Factor * this->EvaluateInteractionCoefficient(Index1, Index2, Index4, Index3, 1, 1);
          
          this->InteractionFactorsdpdpdpdp[i][Index] += Factor * this->EvaluateInteractionCoefficient(Index1, Index2, Index3, Index4, 1, 1);
          this->InteractionFactorsdpdpdpdp[i][Index] += Factor * this->EvaluateInteractionCoefficient(Index2, Index1, Index4, Index3, 1, 1);
          this->InteractionFactorsdpdpdpdp[i][Index] -= Factor * this->EvaluateInteractionCoefficient(Index2, Index1, Index3, Index4, 1, 1);
          this->InteractionFactorsdpdpdpdp[i][Index] -= Factor * this->EvaluateInteractionCoefficient(Index1, Index2, Index4, Index3, 1, 1);
          
          this->InteractionFactorsumumumum[i][Index] += Factor * this->EvaluateInteractionCoefficient(Index1, Index2, Index3, Index4, 0, 0);
          this->InteractionFactorsumumumum[i][Index] += Factor * this->EvaluateInteractionCoefficient(Index2, Index1, Index4, Index3, 0, 0);
          this->InteractionFactorsumumumum[i][Index] -= Factor * this->EvaluateInteractionCoefficient(Index2, Index1, Index3, Index4, 0, 0);
          this->InteractionFactorsumumumum[i][Index] -= Factor * this->EvaluateInteractionCoefficient(Index1, Index2, Index4, Index3, 0, 0);
          
          this->InteractionFactorsdmdmdmdm[i][Index] += Factor * this->EvaluateInteractionCoefficient(Index1, Index2, Index3, Index4, 0, 0);
          this->InteractionFactorsdmdmdmdm[i][Index] += Factor * this->EvaluateInteractionCoefficient(Index2, Index1, Index4, Index3, 0, 0);
          this->InteractionFactorsdmdmdmdm[i][Index] -= Factor * this->EvaluateInteractionCoefficient(Index2, Index1, Index3, Index4, 0, 0);
          this->InteractionFactorsdmdmdmdm[i][Index] -= Factor * this->EvaluateInteractionCoefficient(Index1, Index2, Index4, Index3, 0, 0);
          
          
          

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
                    
          double FactorInter = Factor;
          this->InteractionFactorsupumupum[i][Index] = -FactorInter * this->EvaluateInteractionCoefficient(Index1, Index2, Index4, Index3, 1, 0);
          this->InteractionFactorsupumupum[i][Index] -= FactorInter * this->EvaluateInteractionCoefficient(Index2, Index1, Index3, Index4, 0, 1);
          
          this->InteractionFactorsupdpupdp[i][Index] = -FactorInter * this->EvaluateInteractionCoefficient(Index1, Index2, Index4, Index3, 1, 1);
          this->InteractionFactorsupdpupdp[i][Index] -= FactorInter * this->EvaluateInteractionCoefficient(Index2, Index1, Index3, Index4, 1, 1);
          
          this->InteractionFactorsupdmupdm[i][Index] = -FactorInter * this->EvaluateInteractionCoefficient(Index1, Index2, Index4, Index3, 1, 0);
          this->InteractionFactorsupdmupdm[i][Index] -= FactorInter * this->EvaluateInteractionCoefficient(Index2, Index1, Index3, Index4, 0, 1);
          
          this->InteractionFactorsumdpumdp[i][Index] = -FactorInter * this->EvaluateInteractionCoefficient(Index1, Index2, Index4, Index3, 0, 1);  
          this->InteractionFactorsumdpumdp[i][Index] -= FactorInter * this->EvaluateInteractionCoefficient(Index2, Index1, Index3, Index4, 1, 0);          
          
          this->InteractionFactorsumdmumdm[i][Index] = -FactorInter * this->EvaluateInteractionCoefficient(Index1, Index2, Index4, Index3, 0, 0);
          this->InteractionFactorsumdmumdm[i][Index] -= FactorInter * this->EvaluateInteractionCoefficient(Index2, Index1, Index3, Index4, 0, 0);          
          
          this->InteractionFactorsdpdmdpdm[i][Index] = -FactorInter * this->EvaluateInteractionCoefficient(Index1, Index2, Index4, Index3, 1, 0);
          this->InteractionFactorsdpdmdpdm[i][Index] -= FactorInter * this->EvaluateInteractionCoefficient(Index2, Index1, Index3, Index4, 0, 1);
          

          
          ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }

  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}


// evaluate V(q) for the Coulomb interaction
//
// kx = component of momentum along first Bravais vector
// ky = component of momentum along second Bravais vector
// return value = amplitude of V(q)

double ParticleOnSuperlatticeFourBandCoulombHamiltonian::GetVofQ (int kx, int ky)
{
    if (kx == 0 && ky == 0)
        return 0.0;
    double Q = this->TightBindingModel->EvaluateNormQ(kx, ky);
    double TmpV = (1.0 - exp(-Q*this->ScreeningDistance)) / Q;    
//     cout << kx << " " << ky << " " << Q << " " << UPotential << " " << TmpV << endl;
    return TmpV;
    
}

/*
// compute the amplitude of the Coulomb interaction
//
// Index1 = momentum of the first creation operator
// Index2 = momentum of the second creation operator
// Index3 = momentum of the first annihilation operator
// Index4 = momentum of the second annihilation operator
// valleyIndex1 = valley index of the first density operator
// valleyIndex2 = valley index of the second density operator
// return value = complex amplitude

Complex ParticleOnSuperlatticeFourBandCoulombHamiltonian::EvaluateInteractionCoefficient(int Index1, int Index2, int Index3, int Index4, int valleyIndex1, int valleyIndex2)
{
    int kx1 = Index1 / this->NbrSiteY;
    int ky1 = Index1 % this->NbrSiteY;
	int kx2 = Index2 / this->NbrSiteY;
	int ky2 = Index2 % this->NbrSiteY;
    
    int kx3 = Index3 / this->NbrSiteY;
    int ky3 = Index3 % this->NbrSiteY;
	int kx4 = Index4 / this->NbrSiteY;
	int ky4 = Index4 % this->NbrSiteY;
    
    int Qx = (this->NbrSiteX + kx1 - kx4) % this->NbrSiteX;
    int Qy = (this->NbrSiteY + ky1 - ky4) % this->NbrSiteY;
        
    int TmpQx = Qx;
    int TmpQy = Qy;
    
    Complex Sum = 0.0;
    Complex Coefficient;
    
    int MaxBZX = this->TruncateQx;
    int MaxBZY = this->TruncateQy;
    if (Qx != 0)
        MaxBZX -= 1;
    if (Qy != 0)
        MaxBZY -=1;
        
    for (int BZIndexX = -this->TruncateQx; BZIndexX <= MaxBZX; ++BZIndexX)
    {
        for (int BZIndexY = -this->TruncateQy; BZIndexY <= MaxBZY; ++BZIndexY)
        {
            TmpQx = Qx + BZIndexX * this->NbrSiteX;
            TmpQy = Qy + BZIndexY * this->NbrSiteY;
            Coefficient = this->ComputeDensityFormFactor(kx4, ky4, TmpQx, TmpQy, valleyIndex1) * this->ComputeDensityFormFactor(kx3, ky3, -TmpQx, -TmpQy, valleyIndex2) * this->GetVofQ(TmpQx, TmpQy);
            
            Sum += Coefficient;
            
        }
    }    
  return Sum;    
}*/


// compute the form factor for the density operator 
// 
// kx = momentum along x of annihilation operator
// ky = momentum along y of creation operator
// qx = momentum transfer along x direction
// qy = momentum transfer along y direction
// valleyIndex = valley index of density operator

Complex ParticleOnSuperlatticeFourBandCoulombHamiltonian::ComputeDensityFormFactor(int kx, int ky, int qx, int qy, int valleyIndex)
{
  return this->TightBindingModel->ComputeDensityFormFactor(kx, ky, qx, qy, valleyIndex);
}


// compute the amplitude of the projected Coulomb interaction
//
// Index1 = momentum of the first creation operator
// Index2 = momentum of the second creation operator
// Index3 = momentum of the first annihilation operator
// Index4 = momentum of the second annihilation operator
// valleyIndex1 = valley index of the first density operator
// valleyIndex2 = valley index of the second density operator
// return value = complex amplitude

Complex ParticleOnSuperlatticeFourBandCoulombHamiltonian::EvaluateInteractionCoefficient(int Index1, int Index2, int Index3, int Index4, int valleyIndex1, int valleyIndex2)
{
    int kx1 = Index1 / this->NbrSiteY;
    int ky1 = Index1 % this->NbrSiteY;
	int kx2 = Index2 / this->NbrSiteY;
	int ky2 = Index2 % this->NbrSiteY;
    
    int kx3 = Index3 / this->NbrSiteY;
    int ky3 = Index3 % this->NbrSiteY;
	int kx4 = Index4 / this->NbrSiteY;
	int ky4 = Index4 % this->NbrSiteY;
    
    int Qx = (this->NbrSiteX + kx1 - kx4) % this->NbrSiteX;
    int Qy = (this->NbrSiteY + ky1 - ky4) % this->NbrSiteY;
    
    
    int TmpQx;
    int TmpQy = Qy;
    
    Complex Sum = 0.0;
    Complex Coefficient = 0.0;
    Complex Precision = 1.0;
    
    
    
    while (Norm(Precision) > MACHINE_PRECISION)
    {
//         TmpQx = Qx + ((-TmpQy / 2) % this->NbrSiteX);
        TmpQx = Qx;
        if ((TmpQy != 0) || (TmpQx != 0))
        {
            Coefficient = this->ComputeDensityFormFactor(kx4, ky4, TmpQx, TmpQy, valleyIndex1) * this->ComputeDensityFormFactor(kx3, ky3, -TmpQx, -TmpQy, valleyIndex2) * this->GetVofQ(TmpQx, TmpQy);
            Precision = Coefficient;
//             cout << TmpQx << " " << TmpQy << " " << Precision << endl;
        }
        else
        {
            Coefficient = 0.0; // yields non-zero terms only for non-singular interactions
            Precision = 1.0;
        }
        
        TmpQx += this->NbrSiteX; 
        while (Norm(Precision) > MACHINE_PRECISION)
        {
            Precision = this->ComputeDensityFormFactor(kx4, ky4, TmpQx, TmpQy, valleyIndex1) * this->ComputeDensityFormFactor(kx3, ky3, -TmpQx, -TmpQy, valleyIndex2) * this->GetVofQ(TmpQx, TmpQy);
            Coefficient += Precision;
//             cout << TmpQx << " " << TmpQy << " " << Precision << endl;
            TmpQx += this->NbrSiteX;
        }        
//         TmpQx = Qx + ((-TmpQy / 2) % this->NbrSiteX) - this->NbrSiteX;
        TmpQx = Qx + -this->NbrSiteX;
        Precision = 1.0;
        while (Norm(Precision) > MACHINE_PRECISION)
        {
            Precision = this->ComputeDensityFormFactor(kx4, ky4, TmpQx, TmpQy, valleyIndex1) * this->ComputeDensityFormFactor(kx3, ky3, -TmpQx, -TmpQy, valleyIndex2) * this->GetVofQ(TmpQx, TmpQy);
            Coefficient += Precision;
//             cout << TmpQx << " " << TmpQy << " " << Precision << endl;
            TmpQx -= this->NbrSiteX;
        }  
        
        Sum += Coefficient;
        TmpQy += this->NbrSiteY;
        Precision = Coefficient;
    }
    
    TmpQy = Qy - this->NbrSiteY;
    Precision = 1.0;
    
    while (Norm(Precision) > MACHINE_PRECISION)
    {
//         TmpQx = Qx + ((-TmpQy / 2) % this->NbrSiteX);
        TmpQx = Qx;
        Coefficient = 0.0;
        Precision = 1.0;
         while (Norm(Precision) > MACHINE_PRECISION)
        {
            Precision = this->ComputeDensityFormFactor(kx4, ky4, TmpQx, TmpQy, valleyIndex1) * this->ComputeDensityFormFactor(kx3, ky3, -TmpQx, -TmpQy, valleyIndex2) * this->GetVofQ(TmpQx, TmpQy);
            Coefficient += Precision;
//             cout << TmpQx << " " << TmpQy << " " << Precision << endl;
            TmpQx += this->NbrSiteX;
        }        
//         TmpQx = Qx + ((-TmpQy / 2) % this->NbrSiteX) - this->NbrSiteX;
        TmpQx = Qx - this->NbrSiteX;
        Precision = 1.0;
        while (Norm(Precision) > MACHINE_PRECISION)
        {
            Precision = this->ComputeDensityFormFactor(kx4, ky4, TmpQx, TmpQy, valleyIndex1) * this->ComputeDensityFormFactor(kx3, ky3, -TmpQx, -TmpQy, valleyIndex2) * this->GetVofQ(TmpQx, TmpQy);
            Coefficient += Precision;
//             cout << TmpQx << " " << TmpQy << " " << Precision << endl;
            TmpQx -= this->NbrSiteX;
        }  
        
        Sum += Coefficient;
        TmpQy -= this->NbrSiteY;
        Precision = Coefficient;
    }
    
  return Sum;          
    
}
