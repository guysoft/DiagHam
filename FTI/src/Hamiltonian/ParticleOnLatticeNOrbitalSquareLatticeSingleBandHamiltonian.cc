////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//    class of pyrochlore slab lattice model with interacting particles       //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 17/05/2012                      //
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
#include "Hamiltonian/ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <cmath>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;
using std::sin;
using std::cos;

// constructor
//
ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian::ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian()
{
  this->BandIndex = 0;
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrCellX = number of sites in the x direction
// nbrCellY = number of sites in the y direction
// uPotential = strength of the repulsive two body onsite interaction
// vPotential = strength of the repulsive two body neareast neighbor interaction
// bandIndex = index of the band to consider
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian::ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, int nbrCellY, 
															 double uPotential,  double vPotential,  
															 Abstract2DTightBindingModel* tightBindingModel, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrCellX;
  this->NbrSiteY = nbrCellY;
  this->LzMax = nbrCellX * nbrCellY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);

  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->BandIndex = 0;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      cout << "fast = ";
      PrintMemorySize(cout, TmpMemory)<< endl;
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian::~ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix[this->TightBindingModel->GetNbrStatePerBand()];
  if (this->FlatBand == false)
    {
      this->OneBodyInteractionFactors = new double [this->TightBindingModel->GetNbrStatePerBand()];
    }
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	if (this->FlatBand == false)
	  this->OneBodyInteractionFactors[Index] = 0.5 * this->TightBindingModel->GetEnergy(0, Index);
	OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
      }

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      this->NbrSectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrSectorIndicesPerSum = new int[this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	this->NbrSectorIndicesPerSum[i] = 0;      
      for (int k1a = 0; k1a < this->NbrSiteX; ++k1a)
	for (int k2a = 0; k2a < this->NbrSiteX; ++k2a)
	  for (int k1b = 0; k1b < this->NbrSiteY; ++k1b)
	    for (int k2b = 0; k2b < this->NbrSiteY; ++k2b) 
	      {
		int Index1 = (k1a * this->NbrSiteY) + k1b;
		int Index2 = (k2a * this->NbrSiteY) + k2b;
		if (Index1 < Index2)
		  ++this->NbrSectorIndicesPerSum[(((k1a + k2a) % this->NbrSiteX) *  this->NbrSiteY) + ((k1b + k2b) % this->NbrSiteY)];    
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
      for (int k1a = 0; k1a < this->NbrSiteX; ++k1a)
	for (int k2a = 0; k2a < this->NbrSiteX; ++k2a)
	  for (int k1b = 0; k1b < this->NbrSiteY; ++k1b)
	    for (int k2b = 0; k2b < this->NbrSiteY; ++k2b) 
	      {
		int Index1 = (k1a * this->NbrSiteY) + k1b;
		int Index2 = (k2a * this->NbrSiteY) + k2b;
		if (Index1 < Index2)
		  {
		    int TmpSum = (((k1a + k2a) % this->NbrSiteX) *  this->NbrSiteY) + ((k1b + k2b) % this->NbrSiteY);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }
      double FactorUAB = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUAC = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUBC = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));

      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      int k1a = Index1 / this->NbrSiteY;
	      int k1b = Index1 % this->NbrSiteY;
	      int k2a = Index2 / this->NbrSiteY;
	      int k2b = Index2 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int k3a = Index3 / this->NbrSiteY;
		  int k3b = Index3 % this->NbrSiteY;
		  int k4a = Index4 / this->NbrSiteY;
		  int k4b = Index4 % this->NbrSiteY;
		  
		  this->InteractionFactors[i][Index] = FactorUAB * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementAB(k2a, k2b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUAB * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementAB(k1a, k1b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUAB * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementAB(k2a, k2b, k3a, k3b);
 		  this->InteractionFactors[i][Index] += FactorUAB * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementAB(k1a, k1b, k3a, k3b);

 		  this->InteractionFactors[i][Index] += FactorUAC * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementAC(k2a, k2b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUAC * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementAC(k1a, k1b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUAC * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementAC(k2a, k2b, k3a, k3b);
 		  this->InteractionFactors[i][Index] += FactorUAC * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementAC(k1a, k1b, k3a, k3b);

 		  this->InteractionFactors[i][Index] += FactorUBC * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBC(k1a, k1b, k2a, k2b, k3a, k3b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUBC * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBC(k2a, k2b, k1a, k1b, k3a, k3b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUBC * (Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBC(k1a, k1b, k2a, k2b, k4a, k4b, k3a, k3b);
 		  this->InteractionFactors[i][Index] += FactorUBC * (Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1] * Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2]) * this->ComputeTwoBodyMatrixElementBC(k2a, k2b, k1a, k1b, k4a, k4b, k3a, k3b);

		  this->InteractionFactors[i][Index] *= -2.0;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  else
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
		if (Index1 <= Index2)
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
		if (Index1 <= Index2)
		  {
		    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }
      double FactorU = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      if ((this->FlatBand == false) || (this->VPotential != 0.0))
	FactorU *= this->UPotential;
      double FactorV = this->VPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      
      this->InteractionFactors = new Complex* [this->NbrSectorSums];
      for (int i = 0; i < this->NbrSectorSums; ++i)
	{
	  this->InteractionFactors[i] = new Complex[this->NbrSectorIndicesPerSum[i] * this->NbrSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrSectorIndicesPerSum[i]; ++j1)
	    {
	      int Index1 = this->SectorIndicesPerSum[i][j1 << 1];
	      int Index2 = this->SectorIndicesPerSum[i][(j1 << 1) + 1];
	      int kx1 = Index1 / this->NbrSiteY;
	      int ky1 = Index1 % this->NbrSiteY;
	      int kx2 = Index2 / this->NbrSiteY;
	      int ky2 = Index2 % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrSectorIndicesPerSum[i]; ++j2)
		{
		  int Index3 = this->SectorIndicesPerSum[i][j2 << 1];
		  int Index4 = this->SectorIndicesPerSum[i][(j2 << 1) + 1];
		  int kx3 = Index3 / this->NbrSiteY;
		  int ky3 = Index3 % this->NbrSiteY;
		  int kx4 = Index4 / this->NbrSiteY;
		  int ky4 = Index4 % this->NbrSiteY;
		  Complex Tmp = 0.0;
		  for (int Site1 = 0; Site1 < this->TightBindingModel->GetNbrBands(); ++Site1)
		    {
		  for (int Site2 = 0; Site2 < this->TightBindingModel->GetNbrBands(); ++Site2)
		    {
		      Tmp += (Conj(OneBodyBasis[Index1][BandIndex][Site1]) * OneBodyBasis[Index3][BandIndex][Site1] * Conj(OneBodyBasis[Index2][BandIndex][Site2]) * OneBodyBasis[Index4][BandIndex][Site2]);
		      Tmp += (Conj(OneBodyBasis[Index2][BandIndex][Site1]) * OneBodyBasis[Index3][BandIndex][Site1] * Conj(OneBodyBasis[Index1][BandIndex][Site2]) * OneBodyBasis[Index4][BandIndex][Site2]);
		      Tmp += (Conj(OneBodyBasis[Index1][BandIndex][Site1]) * OneBodyBasis[Index4][BandIndex][Site1] * Conj(OneBodyBasis[Index2][BandIndex][Site2]) * OneBodyBasis[Index3][BandIndex][Site2]);
		      Tmp += (Conj(OneBodyBasis[Index2][BandIndex][Site1]) * OneBodyBasis[Index4][BandIndex][Site1] * Conj(OneBodyBasis[Index1][BandIndex][Site2]) * OneBodyBasis[Index3][BandIndex][Site2]);
		    }
		    }
		  Tmp *= FactorU;
		  
		  		  
		  this->InteractionFactors[i][Index] =  Tmp;
		  
		  if(this->VPotential != 0.0)
		  {
		    Complex Tmp2 = 0.0;
		    int NbrLayers = (this->TightBindingModel->GetNbrBands() + 1)/4;
		    for (int Site = 0; Site < NbrLayers; ++Site)
		    {
		      int PosA = 4*Site;
		      int PosB = 4*Site+1;
		      int PosC = 4*Site+2;
		      Tmp2 = (Conj(OneBodyBasis[Index1][BandIndex][PosA]) * OneBodyBasis[Index3][BandIndex][PosA] * Conj(OneBodyBasis[Index2][BandIndex][PosB]) * OneBodyBasis[Index4][BandIndex][PosB]) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx4, ky4);
		      Tmp2 += (Conj(OneBodyBasis[Index2][BandIndex][PosA]) * OneBodyBasis[Index3][BandIndex][PosA] * Conj(OneBodyBasis[Index1][BandIndex][PosB]) * OneBodyBasis[Index4][BandIndex][PosB]) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx4, ky4);
		      Tmp2 += (Conj(OneBodyBasis[Index1][BandIndex][PosA]) * OneBodyBasis[Index4][BandIndex][PosA] * Conj(OneBodyBasis[Index2][BandIndex][PosB]) * OneBodyBasis[Index3][BandIndex][PosB]) * this->ComputeTwoBodyMatrixElementAB(kx2, ky2, kx3, ky3);
		      Tmp2 += (Conj(OneBodyBasis[Index2][BandIndex][PosA]) * OneBodyBasis[Index4][BandIndex][PosA] * Conj(OneBodyBasis[Index1][BandIndex][PosB]) * OneBodyBasis[Index3][BandIndex][PosB]) * this->ComputeTwoBodyMatrixElementAB(kx1, ky1, kx3, ky3);

 		  Tmp2 += (Conj(OneBodyBasis[Index1][BandIndex][PosA]) * OneBodyBasis[Index3][BandIndex][PosA] * Conj(OneBodyBasis[Index2][BandIndex][PosC]) * OneBodyBasis[Index4][BandIndex][PosC]) * this->ComputeTwoBodyMatrixElementAC(kx2, ky2, kx4, ky4);
 		  Tmp2 += (Conj(OneBodyBasis[Index2][BandIndex][PosA]) * OneBodyBasis[Index3][BandIndex][PosA] * Conj(OneBodyBasis[Index1][BandIndex][PosC]) * OneBodyBasis[Index4][BandIndex][PosC]) * this->ComputeTwoBodyMatrixElementAC(kx1, ky1, kx4, ky4);
 		  Tmp2 += (Conj(OneBodyBasis[Index1][BandIndex][PosA]) * OneBodyBasis[Index4][BandIndex][PosA] * Conj(OneBodyBasis[Index2][BandIndex][PosC]) * OneBodyBasis[Index3][BandIndex][PosC]) * this->ComputeTwoBodyMatrixElementAC(kx2, ky2, kx3, ky3);
 		  Tmp2 += (Conj(OneBodyBasis[Index2][BandIndex][PosA]) * OneBodyBasis[Index4][BandIndex][PosA] * Conj(OneBodyBasis[Index1][BandIndex][PosC]) * OneBodyBasis[Index3][BandIndex][PosC]) * this->ComputeTwoBodyMatrixElementAC(kx1, ky1, kx3, ky3);

 		  Tmp2 += (Conj(OneBodyBasis[Index1][BandIndex][PosB]) * OneBodyBasis[Index3][BandIndex][PosB] * Conj(OneBodyBasis[Index2][BandIndex][PosC]) * OneBodyBasis[Index4][BandIndex][PosC]) * this->ComputeTwoBodyMatrixElementBC(kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4);
 		  Tmp2 += (Conj(OneBodyBasis[Index2][BandIndex][PosB]) * OneBodyBasis[Index3][BandIndex][PosB] * Conj(OneBodyBasis[Index1][BandIndex][PosC]) * OneBodyBasis[Index4][BandIndex][PosC]) * this->ComputeTwoBodyMatrixElementBC(kx2, ky2, kx1, ky1, kx3, ky3, kx4, ky4);
 		  Tmp2 += (Conj(OneBodyBasis[Index1][BandIndex][PosB]) * OneBodyBasis[Index4][BandIndex][PosB] * Conj(OneBodyBasis[Index2][BandIndex][PosC]) * OneBodyBasis[Index3][BandIndex][PosC]) * this->ComputeTwoBodyMatrixElementBC(kx1, ky1, kx2, ky2, kx4, ky4, kx3, ky3);
 		  Tmp2 += (Conj(OneBodyBasis[Index2][BandIndex][PosB]) * OneBodyBasis[Index4][BandIndex][PosB] * Conj(OneBodyBasis[Index1][BandIndex][PosC]) * OneBodyBasis[Index3][BandIndex][PosC]) * this->ComputeTwoBodyMatrixElementBC(kx2, ky2, kx1, ky1, kx4, ky4, kx3, ky3);
		    }
		    Tmp2 *= FactorV;
		    
		  this->InteractionFactors[i][Index] += Tmp2;
		  }
		  
		  if (Index3 == Index4)
		    this->InteractionFactors[i][Index] *= 0.5;
		  if (Index1 == Index2)
		    this->InteractionFactors[i][Index] *= 0.5;
		  this->InteractionFactors[i][Index] *= 2.0;
		  
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;

  delete [] OneBodyBasis;
}

// conventions adopted for matrix elements:
// modified with respect to Tang, Mei, Wen:
// unit cell is triangle standing on base. bottom left corner site A and bottom right corner B, tip is site C.


// compute the matrix element for the two body interaction between two sites A and B 
//
// k1a = creation momentum along x for the B site
// k1b = creation momentum along y for the B site
// k2a = annihilation momentum along x for the B site
// k2b = annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementAB(int k1a, int k1b, int k2a, int k2b)
{
  Complex Tmp = 1 + Phase(this->KyFactor * ((double) (-k2a + k1a )));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A and C 
//
// k1a = creation momentum along x for the C site
// k1b = creation momentum along y for the C site
// k2a = annihilation momentum along x for the C site
// k2b = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementAC(int k1a, int k1b, int k2a, int k2b)
{
 Complex Tmp = 1 + Phase(this->KyFactor * ((double) (-k2b + k1b )));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B and C 
//
// k1a = creation momentum along x for the B site
// k1b = creation momentum along y for the B site
// k2a = creation momentum along x for the C site
// k2b = creation momentum along y for the C site
// k3a = annihilation momentum along x for the B site
// k3b = annihilation momentum along y for the B site
// k4a = annihilation momentum along x for the C site
// k4b = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementBC(int k1a, int k1b, int k2a, int k2b, int k3a, int k3b, int k4a, int k4b)
{
  Complex Tmp = 1 + Phase(this->KxFactor * ((double) (k4a - k2a )) - this->KyFactor * ((double) (k4b - k2b)));
  return Tmp;
}


// compute the matrix element for on-site two body interaction involving A sites
//
// return value = corresponding matrix element

Complex ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteAA()
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

Complex ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return 1.0;
}

// compute the matrix element for on-site two body interaction involving C sites
//
// kx1 = first creation momentum along x for the C site
// ky1 = first creation momentum along y for the C site
// kx2 = second creation momentum along x for the C site
// ky2 = second creation momentum along y for the C site
// kx3 = first annihilation momentum along x for the C site
// ky3 = first annihilation momentum along y for the C site
// kx4 = second annihilation momentum along x for the C site
// ky4 = second annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeNOrbitalSquareLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteCC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return 1.0;
}

