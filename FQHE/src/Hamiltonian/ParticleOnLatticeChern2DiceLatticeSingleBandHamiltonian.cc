////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//             class of dice lattice model with interacting particles         //
//                   in the single chern 2 band approximation                 // 
//                                                                            //
//                        last modification : 16/09/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian.h"
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
// t = nearest neighbor hopping amplitude
// epsilon = on site energy for site 3
// lambda = Rashba spin orbit coupling strength
// bfield1 = magnetic field strength on sites 1 and 2
// bfield3 = magnetic field strength on site 3
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrSiteX, 
												       int nbrSiteY, double uPotential, double t, double epsilon, double lambda, double bfield1, double bfield3, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);

  this->HamiltonianShift = 0.0;

  this->THopping = t;
  this->Epsilon = epsilon;
  this->Lambda = lambda;
  this->BField1 = bfield1;
  this->BField3 = bfield3;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;

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

ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::~ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
  if (this->FlatBand == false)
    this->OneBodyInteractionFactors = new double [this->NbrSiteX * this->NbrSiteY];
  this->ComputeOneBodyMatrices(OneBodyBasis);

  int BandIndex = 2;
 
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
      double FactorUOnSiteA = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUOnSiteB = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUOnSiteC = 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));

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

 		  this->InteractionFactors[i][Index] = FactorUOnSiteA * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementOnSiteA();
 		  this->InteractionFactors[i][Index] -= FactorUOnSiteA * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index3][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index4][BandIndex][1]) * this->ComputeTwoBodyMatrixElementOnSiteA();
 		  this->InteractionFactors[i][Index] -= FactorUOnSiteA * (Conj(OneBodyBasis[Index1][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index2][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementOnSiteA();
 		  this->InteractionFactors[i][Index] += FactorUOnSiteA * (Conj(OneBodyBasis[Index2][BandIndex][0]) * OneBodyBasis[Index4][BandIndex][0] * Conj(OneBodyBasis[Index1][BandIndex][1]) * OneBodyBasis[Index3][BandIndex][1]) * this->ComputeTwoBodyMatrixElementOnSiteA();

 		  this->InteractionFactors[i][Index] += FactorUOnSiteB * (Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2] * Conj(OneBodyBasis[Index2][BandIndex][3]) * OneBodyBasis[Index4][BandIndex][3]) * this->ComputeTwoBodyMatrixElementOnSiteB();
 		  this->InteractionFactors[i][Index] -= FactorUOnSiteB * (Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index3][BandIndex][2] * Conj(OneBodyBasis[Index1][BandIndex][3]) * OneBodyBasis[Index4][BandIndex][3]) * this->ComputeTwoBodyMatrixElementOnSiteB();
 		  this->InteractionFactors[i][Index] -= FactorUOnSiteB * (Conj(OneBodyBasis[Index1][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2] * Conj(OneBodyBasis[Index2][BandIndex][3]) * OneBodyBasis[Index3][BandIndex][3]) * this->ComputeTwoBodyMatrixElementOnSiteB();
 		  this->InteractionFactors[i][Index] += FactorUOnSiteB * (Conj(OneBodyBasis[Index2][BandIndex][2]) * OneBodyBasis[Index4][BandIndex][2] * Conj(OneBodyBasis[Index1][BandIndex][3]) * OneBodyBasis[Index3][BandIndex][3]) * this->ComputeTwoBodyMatrixElementOnSiteB();

 		  this->InteractionFactors[i][Index] += FactorUOnSiteC * (Conj(OneBodyBasis[Index1][BandIndex][4]) * OneBodyBasis[Index3][BandIndex][4] * Conj(OneBodyBasis[Index2][BandIndex][5]) * OneBodyBasis[Index4][BandIndex][5]) * this->ComputeTwoBodyMatrixElementOnSiteC();
 		  this->InteractionFactors[i][Index] -= FactorUOnSiteC * (Conj(OneBodyBasis[Index2][BandIndex][4]) * OneBodyBasis[Index3][BandIndex][4] * Conj(OneBodyBasis[Index1][BandIndex][5]) * OneBodyBasis[Index4][BandIndex][5]) * this->ComputeTwoBodyMatrixElementOnSiteC();
 		  this->InteractionFactors[i][Index] -= FactorUOnSiteC * (Conj(OneBodyBasis[Index1][BandIndex][4]) * OneBodyBasis[Index4][BandIndex][4] * Conj(OneBodyBasis[Index2][BandIndex][5]) * OneBodyBasis[Index3][BandIndex][5]) * this->ComputeTwoBodyMatrixElementOnSiteC();
 		  this->InteractionFactors[i][Index] += FactorUOnSiteC * (Conj(OneBodyBasis[Index2][BandIndex][4]) * OneBodyBasis[Index4][BandIndex][4] * Conj(OneBodyBasis[Index1][BandIndex][5]) * OneBodyBasis[Index3][BandIndex][5]) * this->ComputeTwoBodyMatrixElementOnSiteC();

		  this->InteractionFactors[i][Index] *= -2.0;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

// compute the matrix element for the on-site two body interaction for site A
//
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteA()
{
  return 1.0;
}

// compute the matrix element for the on-site two body interaction for site B 
//
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteB()
{
  return 1.0;
}

// compute the matrix element for the on-site two body interaction for site C 
//
// return value = corresponding matrix element

Complex ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteC()
{
  return 1.0;
}

// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = array of one body transformation matrices

void ParticleOnLatticeChern2DiceLatticeSingleBandHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
  double KX;
  double KY;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  KX = this->KxFactor * (((double) kx) + this->GammaX);
	  KY = this->KyFactor * (((double) ky) + this->GammaY);
	  int Index = (kx * this->NbrSiteY) + ky;

          double angle = 2.0 * M_PI / 3.0;
          Complex GammaK = this->THopping * (1.0 + Phase(-KX) + Phase(-KY));
          Complex GammaP = this->Lambda * (1.0 + Phase(-KX + angle) + Phase(-KY - angle));
          Complex GammaM = this->Lambda * (1.0 + Phase(-KX - angle) + Phase(-KY + angle));

	  HermitianMatrix TmpOneBodyHamiltonian(6, true);
 	  TmpOneBodyHamiltonian.SetMatrixElement(0, 4, Conj(GammaK));
 	  TmpOneBodyHamiltonian.SetMatrixElement(1, 5, Conj(GammaK));
 	  TmpOneBodyHamiltonian.SetMatrixElement(0, 5, I()*Conj(GammaP));
 	  TmpOneBodyHamiltonian.SetMatrixElement(1, 4, I()*Conj(GammaM));

 	  TmpOneBodyHamiltonian.SetMatrixElement(2, 4, GammaK);
 	  TmpOneBodyHamiltonian.SetMatrixElement(3, 5, GammaK);
 	  TmpOneBodyHamiltonian.SetMatrixElement(2, 5, I()*GammaM);
 	  TmpOneBodyHamiltonian.SetMatrixElement(3, 4, I()*GammaP);

	  TmpOneBodyHamiltonian.SetMatrixElement(0, 0, this->BField1);
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 1, -this->BField1);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 2, this->BField1);
	  TmpOneBodyHamiltonian.SetMatrixElement(3, 3, -this->BField1);
	  TmpOneBodyHamiltonian.SetMatrixElement(4, 4, this->Epsilon + this->BField3);
	  TmpOneBodyHamiltonian.SetMatrixElement(5, 5, this->Epsilon - this->BField3);

	  TmpOneBodyHamiltonian *= -1.0;

	  ComplexMatrix TmpMatrix(6, 6, true);
	  TmpMatrix.SetToIdentity();
	  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
	  oneBodyBasis[Index] = TmpMatrix;
	  if (this->FlatBand == false)
	    {
	      this->OneBodyInteractionFactors[Index] = 0.5 * TmpDiag(0, 0);
	  }
	  cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << " " << TmpDiag(2, 2) << " " 
	       << TmpDiag(3, 3) << " " << TmpDiag(4, 4) << " " << TmpDiag(5, 5) << endl;
	}
    }
}
