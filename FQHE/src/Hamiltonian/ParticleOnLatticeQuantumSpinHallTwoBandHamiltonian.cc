////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                   class of quatum spin Hall restricted to two band         //
//                                                                            //
//                        last modification : 27/02/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian.h"
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
// bandParameter = band parameter
// szSymmetryBreaking = amplitude of the Sz symmetry breaking term
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
												       int nbrSiteY, double bandParameter, double szSymmetryBreaking, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->HamiltonianShift = 0.0;
  this->BandParameter = bandParameter;
  this->SzSymmetryBreaking = szSymmetryBreaking;
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

ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::~ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
      {
	int Index = (kx1 * this->NbrSiteY) + ky1;
	double d1 = sin (2.0 * M_PI * ((double) kx1) / ((double) this->NbrSiteX));
	double d2 = sin (2.0 * M_PI * ((double) ky1) / ((double) this->NbrSiteY));
	double d3 = (this->BandParameter - cos (2.0 * M_PI * ((double) ky1) / ((double) this->NbrSiteY))
		     - cos (2.0 * M_PI * ((double) kx1) / ((double) this->NbrSiteX)));
	double TmpSzSymmetryBreaking = this->SzSymmetryBreaking * d1;
	HermitianMatrix TmpOneBobyHamiltonian(4, true);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 0, d3);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 1, Complex(d1, -d2));
	TmpOneBobyHamiltonian.SetMatrixElement(1, 1, -d3);
	TmpOneBobyHamiltonian.SetMatrixElement(2, 2, d3);
	TmpOneBobyHamiltonian.SetMatrixElement(2, 3, Complex(-d1, -d2));
	TmpOneBobyHamiltonian.SetMatrixElement(0, 3, TmpSzSymmetryBreaking);
	TmpOneBobyHamiltonian.SetMatrixElement(1, 2, TmpSzSymmetryBreaking);

	TmpOneBobyHamiltonian.SetMatrixElement(3, 3, -d3);
	ComplexMatrix TmpMatrix(4, 4, true);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	TmpMatrix[2][2] = 1.0;
	TmpMatrix[3][3] = 1.0;
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBobyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBobyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	OneBodyBasis[Index] = TmpMatrix;	
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1)  << " " << TmpDiag(2, 2)  << " " << TmpDiag(3, 3) << endl;
      }
 
  this->NbrInterSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
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
      this->NbrIntraSectorSums = this->NbrSiteX * this->NbrSiteY;
      this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	this->NbrIntraSectorIndicesPerSum[i] = 0;      
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
		  this->InteractionFactorsupup[i][Index] = -4.0  *  ((this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 0, 0) * (cos (2.0 * M_PI * ((double) kx2 - kx4) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky2 - ky4) / ((double) this->NbrSiteY))))
								     - (this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 0, 0) * (cos (2.0 * M_PI * ((double) kx1 - kx4) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky1 - ky4) / ((double) this->NbrSiteY))))
								     - (this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 0, 0) * (cos (2.0 * M_PI * ((double) kx2 - kx3) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky2 - ky3) / ((double) this->NbrSiteY))))
								     + (this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 0, 0) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 0, 0) * (cos (2.0 * M_PI * ((double) kx1 - kx3) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky1 - ky3) / ((double) this->NbrSiteY)))));
		  this->InteractionFactorsdowndown[i][Index] = -4.0  *  ((this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 1, 1) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 1, 1) * (cos (2.0 * M_PI * ((double) kx2 - kx4) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky2 - ky4) / ((double) this->NbrSiteY))))
								     - (this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 1, 1) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 1, 1) * (cos (2.0 * M_PI * ((double) kx1 - kx4) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky1 - ky4) / ((double) this->NbrSiteY))))
								     - (this->ComputeBasisContribution(OneBodyBasis, Index1, Index4, 1, 1) * this->ComputeBasisContribution(OneBodyBasis, Index2, Index3, 1, 1) * (cos (2.0 * M_PI * ((double) kx2 - kx3) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky2 - ky3) / ((double) this->NbrSiteY))))
								     + (this->ComputeBasisContribution(OneBodyBasis, Index2, Index4, 1, 1) * this->ComputeBasisContribution(OneBodyBasis, Index1, Index3, 1, 1) * (cos (2.0 * M_PI * ((double) kx1 - kx3) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky1 - ky3) / ((double) this->NbrSiteY)))));;
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
	      int kx1 = this->InterSectorIndicesPerSum[i][j1 << 1] / this->NbrSiteY;
	      int ky1 = this->InterSectorIndicesPerSum[i][j1 << 1] % this->NbrSiteY;
	      int kx2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1] / this->NbrSiteY;
	      int ky2 = this->InterSectorIndicesPerSum[i][(j1 << 1) + 1] % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrInterSectorIndicesPerSum[i]; ++j2)
		{
		  int kx3 = this->InterSectorIndicesPerSum[i][j2 << 1] / this->NbrSiteY;
		  int ky3 = this->InterSectorIndicesPerSum[i][j2 << 1] % this->NbrSiteY;
		  int kx4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1] / this->NbrSiteY;
		  int ky4 = this->InterSectorIndicesPerSum[i][(j2 << 1) + 1] % this->NbrSiteY;
		  this->InteractionFactorsupdown[i][Index] = -2.0 * this->UPotential *((cos (2.0 * M_PI * ((double) kx2 - kx4) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky2 - ky4) / ((double) this->NbrSiteY)))
										       + (cos (2.0 * M_PI * ((double) kx1 - kx3) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky1 - ky3) / ((double) this->NbrSiteY))));
		  ++TotalNbrInteractionFactors;
		  ++Index;
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the part of the interaction coefficient coming from the two band truncated basis
//
// basisMatrices = array of basis matrices
// momentumIndex1 = momentum index for the first particle
// momentumIndex2 = momentum index for the second particle
// bandIndex1 = band index of the first particle
// bandIndex2 = band index of the second particle
// return value = coefficient

Complex ParticleOnLatticeQuantumSpinHallTwoBandHamiltonian::ComputeBasisContribution(ComplexMatrix* basisMatrices, int momentumIndex1, int momentumIndex2, int bandIndex1, int bandIndex2) 
{
  ComplexMatrix& TmpMatrix1 = basisMatrices[momentumIndex1];
  ComplexMatrix& TmpMatrix2 = basisMatrices[momentumIndex2];
  return Conj(TmpMatrix1[bandIndex1][0]) * TmpMatrix2[bandIndex2][0] 
    + Conj(TmpMatrix1[bandIndex1][1]) * TmpMatrix2[bandIndex2][1]
    + Conj(TmpMatrix1[bandIndex1][2]) * TmpMatrix2[bandIndex2][2]
    + Conj(TmpMatrix1[bandIndex1][3]) * TmpMatrix2[bandIndex2][3];
}
