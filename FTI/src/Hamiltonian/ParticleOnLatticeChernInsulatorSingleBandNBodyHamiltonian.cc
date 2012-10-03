////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                   class of hamiltonian with particles on                   //
//               Chern insulator in the single band approximation             //
//                           and n-body interaction                           //
//                                                                            //
//                        last modification : 03/08/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian.h"
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

ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian::ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// bandParameter = band parameter
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian::ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, 
														     int nbrSiteY, double bandParameter, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->HamiltonianShift = 0.0;//4.0 * uPotential;
  this->NBodyValue = 3;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->BandParameter = bandParameter;
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactors = 0;
  this->FastMultiplicationFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  this->HermitianSymmetryFlag = true;
  this->TwoBodyFlag = false;
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

ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian::~ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian()
{
}
  

// evaluate all interaction factors
//   

void ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian::EvaluateInteractionFactors()
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
	HermitianMatrix TmpOneBobyHamiltonian(2, true);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 0, d3);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 1, Complex(d1, -d2));
	TmpOneBobyHamiltonian.SetMatrixElement(1, 1, -d3);
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
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << endl;
      }
  double** CosineTableX = new double*[this->NbrSiteX];
  for (int i = 0; i < this->NbrSiteX; ++i)
    {
      CosineTableX[i] = new double[this->NbrSiteX];
      for (int j = 0; j < this->NbrSiteX; ++j)
	{
	  CosineTableX[i][j] = cos (2.0 * M_PI * ((double) (i - j)) / ((double) this->NbrSiteX));
	}
    }
  double** CosineTableY = new double*[this->NbrSiteY];
  for (int i = 0; i < this->NbrSiteY; ++i)
    {
      CosineTableY[i] = new double[this->NbrSiteY];
      for (int j = 0; j < this->NbrSiteY; ++j)
	{
	  CosineTableY[i][j] = cos (2.0 * M_PI * ((double) (i - j)) / ((double) this->NbrSiteY));
	}
    }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
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
		if (Index1 < Index2)
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
		if (Index1 < Index2)
		  {
		    int TmpSum = (((kx1 + kx2) % this->NbrSiteX) *  this->NbrSiteY) + ((ky1 + ky2) % this->NbrSiteY);
		    this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->SectorIndicesPerSum[TmpSum][1 + (this->NbrSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    ++this->NbrSectorIndicesPerSum[TmpSum];    
		  }
	      }
      double Factor = 0.5 / ((double) this->NbrParticles);
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
		  this->InteractionFactors[i][Index] = ((Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1]) * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1])) * (CosineTableX[kx2][kx4] + CosineTableY[ky2][ky4]);
		  this->InteractionFactors[i][Index] -= ((Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1]) * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1])) * (CosineTableX[kx1][kx4] + CosineTableY[ky1][ky4]);
		  this->InteractionFactors[i][Index] -= ((Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]) * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1])) * (CosineTableX[kx2][kx3] + CosineTableY[ky2][ky3]);
		  this->InteractionFactors[i][Index] += ((Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] + Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]) * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] + Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1])) * (CosineTableX[kx1][kx3] + CosineTableY[ky1][ky3]);
		  this->InteractionFactors[i][Index] += 4.0 * Conj(OneBodyBasis[Index1][0][0]) * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][0] * OneBodyBasis[Index4][0][1];
		  this->InteractionFactors[i][Index] -= 4.0 * Conj(OneBodyBasis[Index2][0][0]) * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][0] * OneBodyBasis[Index4][0][1];
		  this->InteractionFactors[i][Index] -= 4.0 * Conj(OneBodyBasis[Index1][0][0]) * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][0] * OneBodyBasis[Index3][0][1];
		  this->InteractionFactors[i][Index] += 4.0 * Conj(OneBodyBasis[Index2][0][0]) * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][0] * OneBodyBasis[Index3][0][1];
		  this->InteractionFactors[i][Index] *= -4.0 * Factor;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  for (int i = 0; i < this->NbrSiteX; ++i)
    {
      delete[] CosineTableX[i];
    }
  delete[] CosineTableX;
  for (int i = 0; i < this->NbrSiteY; ++i)
    {
      delete[] CosineTableY[i];
    }
  delete[] CosineTableY;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = array of one body transformation matrices

void ParticleOnLatticeChernInsulatorSingleBandNBodyHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
}

