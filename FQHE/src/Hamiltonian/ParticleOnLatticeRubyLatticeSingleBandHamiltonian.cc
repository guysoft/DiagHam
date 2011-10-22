////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//           class of ruby lattice model with interacting particles           //
//                       in the single band approximation                     // 
//                                                                            //
//                        last modification : 17/10/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeRubyLatticeSingleBandHamiltonian.h"
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
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrCellX = number of sites in the x direction
// nbrCellY = number of sites in the y direction
// uPotential = strength of the repulsive two body neareast neighbor interaction
// tr = real part of the hopping amplitude between neareast neighbor sites with same parity
// ti = imaginary part of the hopping amplitude between neareast neighbor sites with same parity
// t1r = real part of the hopping amplitude next neareast neighbor sites with different parity
// t1i = real part of the hopping amplitude next neareast neighbor sites with different parity
// t4 = hopping amplitude along square diagonal
// mus = sublattice chemical potential on A1 sites
// gammaX = boundary condition twisting angle along e_a (measured in units of 2pi)
// gammaY = boundary condition twisting angle along e_b (measured in units of 2pi)
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ParticleOnLatticeRubyLatticeSingleBandHamiltonian(ParticleOnSphere* particles, int nbrParticles, int nbrCellX, 
												       int nbrCellY, double uPotential, double tr, double ti, double t1r, double t1i, double t4, double mus, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrCellX;
  this->NbrSiteY = nbrCellY;
  this->LzMax = nbrCellX * nbrCellY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);

  this->HamiltonianShift = 0.0;
  this->TrHopping = tr;
  this->TiHopping = ti;
  this->T1rHopping = t1r;
  this->T1iHopping = t1i;
  this->T4Hopping = t4;
  this->MuS = mus;
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
      cout << "fast = ";
      PrintMemorySize(cout, TmpMemory)<< endl;
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeRubyLatticeSingleBandHamiltonian::~ParticleOnLatticeRubyLatticeSingleBandHamiltonian()
{
}

// evaluate all interaction factors
//   

void ParticleOnLatticeRubyLatticeSingleBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->NbrSiteX * this->NbrSiteY];
  if (this->FlatBand == false)
    this->OneBodyInteractionFactors = new double [this->NbrSiteX * this->NbrSiteY];
  this->ComputeOneBodyMatrices(OneBodyBasis);
  
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
      double FactorUA1A2 = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUA1A3 = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUA1A5 = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUA1A6 = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUA2A3 = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUA2A4 = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUA2A6 = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUA3A4 = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUA3A5 = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUA4A5 = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUA4A6 = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));
      double FactorUA5A6 = this->UPotential * 0.5 / ((double) (this->NbrSiteX * this->NbrSiteY));


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
		  
		  this->InteractionFactors[i][Index] = FactorUA1A2 * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]) * this->ComputeTwoBodyMatrixElementA1A2();
 		  this->InteractionFactors[i][Index] -= FactorUA1A2 * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]) * this->ComputeTwoBodyMatrixElementA1A2();
 		  this->InteractionFactors[i][Index] -= FactorUA1A2 * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1]) * this->ComputeTwoBodyMatrixElementA1A2();
 		  this->InteractionFactors[i][Index] += FactorUA1A2 * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1]) * this->ComputeTwoBodyMatrixElementA1A2();

 		  this->InteractionFactors[i][Index] += FactorUA1A3 * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2]) * this->ComputeTwoBodyMatrixElementA1A3();
 		  this->InteractionFactors[i][Index] -= FactorUA1A3 * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2]) * this->ComputeTwoBodyMatrixElementA1A3();
 		  this->InteractionFactors[i][Index] -= FactorUA1A3 * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2]) * this->ComputeTwoBodyMatrixElementA1A3();
 		  this->InteractionFactors[i][Index] += FactorUA1A3 * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2]) * this->ComputeTwoBodyMatrixElementA1A3();

 		  this->InteractionFactors[i][Index] += FactorUA1A5 * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index2][0][4]) * OneBodyBasis[Index4][0][4]) * this->ComputeTwoBodyMatrixElementA1A5();
 		  this->InteractionFactors[i][Index] -= FactorUA1A5 * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index1][0][4]) * OneBodyBasis[Index4][0][4]) * this->ComputeTwoBodyMatrixElementA1A5();
 		  this->InteractionFactors[i][Index] -= FactorUA1A5 * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][4]) * OneBodyBasis[Index3][0][4]) * this->ComputeTwoBodyMatrixElementA1A5();
 		  this->InteractionFactors[i][Index] += FactorUA1A5 * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][4]) * OneBodyBasis[Index3][0][4]) * this->ComputeTwoBodyMatrixElementA1A5();

 		  this->InteractionFactors[i][Index] += FactorUA1A6 * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index2][0][5]) * OneBodyBasis[Index4][0][5]) * this->ComputeTwoBodyMatrixElementA1A6(k2a, k2b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUA1A6 * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index1][0][5]) * OneBodyBasis[Index4][0][5]) * this->ComputeTwoBodyMatrixElementA1A6(k1a, k1b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUA1A6 * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][5]) * OneBodyBasis[Index3][0][5]) * this->ComputeTwoBodyMatrixElementA1A6(k2a, k2b, k3a, k3b);
 		  this->InteractionFactors[i][Index] += FactorUA1A6 * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][5]) * OneBodyBasis[Index3][0][5]) * this->ComputeTwoBodyMatrixElementA1A6(k1a, k1b, k3a, k3b);

 		  this->InteractionFactors[i][Index] += FactorUA2A3 * (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2]) * this->ComputeTwoBodyMatrixElementA2A3(k2a, k2b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUA2A3 * (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2]) * this->ComputeTwoBodyMatrixElementA2A3(k1a, k1b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUA2A3 * (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2]) * this->ComputeTwoBodyMatrixElementA2A3(k2a, k2b, k3a, k3b);
 		  this->InteractionFactors[i][Index] += FactorUA2A3 * (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2]) * this->ComputeTwoBodyMatrixElementA2A3(k1a, k1b, k3a, k3b);

 		  this->InteractionFactors[i][Index] += FactorUA2A4 * (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index2][0][3]) * OneBodyBasis[Index4][0][3]) * this->ComputeTwoBodyMatrixElementA2A4();
 		  this->InteractionFactors[i][Index] -= FactorUA2A4 * (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index1][0][3]) * OneBodyBasis[Index4][0][3]) * this->ComputeTwoBodyMatrixElementA2A4();
 		  this->InteractionFactors[i][Index] -= FactorUA2A4 * (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][3]) * OneBodyBasis[Index3][0][3]) * this->ComputeTwoBodyMatrixElementA2A4();
 		  this->InteractionFactors[i][Index] += FactorUA2A4 * (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][3]) * OneBodyBasis[Index3][0][3]) * this->ComputeTwoBodyMatrixElementA2A4();

 		  this->InteractionFactors[i][Index] += FactorUA2A6 * (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index2][0][5]) * OneBodyBasis[Index4][0][5]) * this->ComputeTwoBodyMatrixElementA2A6();
 		  this->InteractionFactors[i][Index] -= FactorUA2A6 * (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index1][0][5]) * OneBodyBasis[Index4][0][5]) * this->ComputeTwoBodyMatrixElementA2A6();
 		  this->InteractionFactors[i][Index] -= FactorUA2A6 * (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][5]) * OneBodyBasis[Index3][0][5]) * this->ComputeTwoBodyMatrixElementA2A6();
 		  this->InteractionFactors[i][Index] += FactorUA2A6 * (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][5]) * OneBodyBasis[Index3][0][5]) * this->ComputeTwoBodyMatrixElementA2A6();

 		  this->InteractionFactors[i][Index] += FactorUA3A4 * (Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2] * Conj(OneBodyBasis[Index2][0][3]) * OneBodyBasis[Index4][0][3]) * this->ComputeTwoBodyMatrixElementA3A4(k2a, k2b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUA3A4 * (Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2] * Conj(OneBodyBasis[Index1][0][3]) * OneBodyBasis[Index4][0][3]) * this->ComputeTwoBodyMatrixElementA3A4(k1a, k1b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUA3A4 * (Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2] * Conj(OneBodyBasis[Index2][0][3]) * OneBodyBasis[Index3][0][3]) * this->ComputeTwoBodyMatrixElementA3A4(k2a, k2b, k3a, k3b);
 		  this->InteractionFactors[i][Index] += FactorUA3A4 * (Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2] * Conj(OneBodyBasis[Index1][0][3]) * OneBodyBasis[Index3][0][3]) * this->ComputeTwoBodyMatrixElementA3A4(k1a, k1b, k3a, k3b);

 		  this->InteractionFactors[i][Index] += FactorUA3A5 * (Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2] * Conj(OneBodyBasis[Index2][0][4]) * OneBodyBasis[Index4][0][4]) * this->ComputeTwoBodyMatrixElementA3A5();
 		  this->InteractionFactors[i][Index] -= FactorUA3A5 * (Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2] * Conj(OneBodyBasis[Index1][0][4]) * OneBodyBasis[Index4][0][4]) * this->ComputeTwoBodyMatrixElementA3A5();
 		  this->InteractionFactors[i][Index] -= FactorUA3A5 * (Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2] * Conj(OneBodyBasis[Index2][0][4]) * OneBodyBasis[Index3][0][4]) * this->ComputeTwoBodyMatrixElementA3A5();
 		  this->InteractionFactors[i][Index] += FactorUA3A5 * (Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2] * Conj(OneBodyBasis[Index1][0][4]) * OneBodyBasis[Index3][0][4]) * this->ComputeTwoBodyMatrixElementA3A5();

 		  this->InteractionFactors[i][Index] += FactorUA4A5 * (Conj(OneBodyBasis[Index1][0][3]) * OneBodyBasis[Index3][0][3] * Conj(OneBodyBasis[Index2][0][4]) * OneBodyBasis[Index4][0][4]) * this->ComputeTwoBodyMatrixElementA4A5();
 		  this->InteractionFactors[i][Index] -= FactorUA4A5 * (Conj(OneBodyBasis[Index2][0][3]) * OneBodyBasis[Index3][0][3] * Conj(OneBodyBasis[Index1][0][4]) * OneBodyBasis[Index4][0][4]) * this->ComputeTwoBodyMatrixElementA4A5();
 		  this->InteractionFactors[i][Index] -= FactorUA4A5 * (Conj(OneBodyBasis[Index1][0][3]) * OneBodyBasis[Index4][0][3] * Conj(OneBodyBasis[Index2][0][4]) * OneBodyBasis[Index3][0][4]) * this->ComputeTwoBodyMatrixElementA4A5();
 		  this->InteractionFactors[i][Index] += FactorUA4A5 * (Conj(OneBodyBasis[Index2][0][3]) * OneBodyBasis[Index4][0][3] * Conj(OneBodyBasis[Index1][0][4]) * OneBodyBasis[Index3][0][4]) * this->ComputeTwoBodyMatrixElementA4A5();

 		  this->InteractionFactors[i][Index] += FactorUA4A6 * (Conj(OneBodyBasis[Index1][0][3]) * OneBodyBasis[Index3][0][3] * Conj(OneBodyBasis[Index2][0][5]) * OneBodyBasis[Index4][0][5]) * this->ComputeTwoBodyMatrixElementA4A6();
 		  this->InteractionFactors[i][Index] -= FactorUA4A6 * (Conj(OneBodyBasis[Index2][0][3]) * OneBodyBasis[Index3][0][3] * Conj(OneBodyBasis[Index1][0][5]) * OneBodyBasis[Index4][0][5]) * this->ComputeTwoBodyMatrixElementA4A6();
 		  this->InteractionFactors[i][Index] -= FactorUA4A6 * (Conj(OneBodyBasis[Index1][0][3]) * OneBodyBasis[Index4][0][3] * Conj(OneBodyBasis[Index2][0][5]) * OneBodyBasis[Index3][0][5]) * this->ComputeTwoBodyMatrixElementA4A6();
 		  this->InteractionFactors[i][Index] += FactorUA4A6 * (Conj(OneBodyBasis[Index2][0][3]) * OneBodyBasis[Index4][0][3] * Conj(OneBodyBasis[Index1][0][5]) * OneBodyBasis[Index3][0][5]) * this->ComputeTwoBodyMatrixElementA4A6();

 		  this->InteractionFactors[i][Index] += FactorUA5A6 * (Conj(OneBodyBasis[Index1][0][4]) * OneBodyBasis[Index3][0][4] * Conj(OneBodyBasis[Index2][0][5]) * OneBodyBasis[Index4][0][5]) * this->ComputeTwoBodyMatrixElementA5A6(k2a, k2b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUA5A6 * (Conj(OneBodyBasis[Index2][0][4]) * OneBodyBasis[Index3][0][4] * Conj(OneBodyBasis[Index1][0][5]) * OneBodyBasis[Index4][0][5]) * this->ComputeTwoBodyMatrixElementA5A6(k1a, k1b, k4a, k4b);
 		  this->InteractionFactors[i][Index] -= FactorUA5A6 * (Conj(OneBodyBasis[Index1][0][4]) * OneBodyBasis[Index4][0][4] * Conj(OneBodyBasis[Index2][0][5]) * OneBodyBasis[Index3][0][5]) * this->ComputeTwoBodyMatrixElementA5A6(k2a, k2b, k3a, k3b);
 		  this->InteractionFactors[i][Index] += FactorUA5A6 * (Conj(OneBodyBasis[Index2][0][4]) * OneBodyBasis[Index4][0][4] * Conj(OneBodyBasis[Index1][0][5]) * OneBodyBasis[Index3][0][5]) * this->ComputeTwoBodyMatrixElementA5A6(k1a, k1b, k3a, k3b);

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
      if (this->FlatBand == false)
	FactorU *= this->UPotential;
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
 		  this->InteractionFactors[i][Index] = FactorU * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteA1A1();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteA1A1();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index3][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteA1A1();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][0]) * OneBodyBasis[Index4][0][0] * Conj(OneBodyBasis[Index1][0][0]) * OneBodyBasis[Index3][0][0]) * this->ComputeTwoBodyMatrixElementOnSiteA1A1();

 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteA2A2();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteA2A2();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteA2A2();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][1]) * OneBodyBasis[Index4][0][1] * Conj(OneBodyBasis[Index1][0][1]) * OneBodyBasis[Index3][0][1]) * this->ComputeTwoBodyMatrixElementOnSiteA2A2();

 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteA3A3();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteA3A3();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index4][0][2] * Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index3][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteA3A3();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][2]) * OneBodyBasis[Index4][0][2] * Conj(OneBodyBasis[Index1][0][2]) * OneBodyBasis[Index3][0][2]) * this->ComputeTwoBodyMatrixElementOnSiteA3A3();

 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][3]) * OneBodyBasis[Index3][0][3] * Conj(OneBodyBasis[Index2][0][3]) * OneBodyBasis[Index4][0][3]) * this->ComputeTwoBodyMatrixElementOnSiteA4A4();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][3]) * OneBodyBasis[Index3][0][3] * Conj(OneBodyBasis[Index1][0][3]) * OneBodyBasis[Index4][0][3]) * this->ComputeTwoBodyMatrixElementOnSiteA4A4();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][3]) * OneBodyBasis[Index4][0][3] * Conj(OneBodyBasis[Index2][0][3]) * OneBodyBasis[Index3][0][3]) * this->ComputeTwoBodyMatrixElementOnSiteA4A4();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][3]) * OneBodyBasis[Index4][0][3] * Conj(OneBodyBasis[Index1][0][3]) * OneBodyBasis[Index3][0][3]) * this->ComputeTwoBodyMatrixElementOnSiteA4A4();

 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][4]) * OneBodyBasis[Index3][0][4] * Conj(OneBodyBasis[Index2][0][4]) * OneBodyBasis[Index4][0][4]) * this->ComputeTwoBodyMatrixElementOnSiteA5A5();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][4]) * OneBodyBasis[Index3][0][4] * Conj(OneBodyBasis[Index1][0][4]) * OneBodyBasis[Index4][0][4]) * this->ComputeTwoBodyMatrixElementOnSiteA5A5();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][4]) * OneBodyBasis[Index4][0][4] * Conj(OneBodyBasis[Index2][0][4]) * OneBodyBasis[Index3][0][4]) * this->ComputeTwoBodyMatrixElementOnSiteA5A5();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][4]) * OneBodyBasis[Index4][0][4] * Conj(OneBodyBasis[Index1][0][4]) * OneBodyBasis[Index3][0][4]) * this->ComputeTwoBodyMatrixElementOnSiteA5A5();

 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][5]) * OneBodyBasis[Index3][0][5] * Conj(OneBodyBasis[Index2][0][5]) * OneBodyBasis[Index4][0][5]) * this->ComputeTwoBodyMatrixElementOnSiteA6A6();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][5]) * OneBodyBasis[Index3][0][5] * Conj(OneBodyBasis[Index1][0][5]) * OneBodyBasis[Index4][0][5]) * this->ComputeTwoBodyMatrixElementOnSiteA6A6();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index1][0][5]) * OneBodyBasis[Index4][0][5] * Conj(OneBodyBasis[Index2][0][5]) * OneBodyBasis[Index3][0][5]) * this->ComputeTwoBodyMatrixElementOnSiteA6A6();
 		  this->InteractionFactors[i][Index] += FactorU * (Conj(OneBodyBasis[Index2][0][5]) * OneBodyBasis[Index4][0][5] * Conj(OneBodyBasis[Index1][0][5]) * OneBodyBasis[Index3][0][5]) * this->ComputeTwoBodyMatrixElementOnSiteA6A6();


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
}

// compute the matrix element for the two body interaction between two sites A and B 
//
// k1a = creation momentum along x for the B site
// k1b = creation momentum along y for the B site
// k2a = annihilation momentum along x for the B site
// k2b = annihilation momentum along y for the B site
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementAB(int k1a, int k1b, int k2a, int k2b)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KxFactor * ((double) (k2a - k1a))));
  //Complex Tmp = Phase (0.5 * (this->KxFactor * ((double) (k2a - k1a))));
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites A1 and A2 
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementA1A2()
{
  return 1.0;
}

// compute the matrix element for the two body interaction between two sites A1 and A3 
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementA1A3()
{
  return 1.0;
}

// compute the matrix element for the two body interaction between two sites A1 and A5
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementA1A5()
{
  return 1.0;
}

// compute the matrix element for the two body interaction between two sites A1 and A6 
//
// k1x = creation momentum along x for the A6 site
// k1y = creation momentum along y for the A6 site
// k2x = annihilation momentum along x for the A6 site
// k2y = annihilation momentum along y for the A6 site
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementA1A6(int k1x, int k1y, int k2x, int k2y)
{
  return Phase(this->KxFactor * ((double) (k2x - k1x)));
}

// compute the matrix element for the two body interaction between two sites A2 and A3 
//
// k1x = creation momentum along x for the A3 site
// k1y = creation momentum along y for the A3 site
// k2x = annihilation momentum along x for the A3 site
// k2y = annihilation momentum along y for the A3 site
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementA2A3(int k1x, int k1y, int k2x, int k2y)
{
  return Phase((this->KxFactor * ((double) (k1x - k2x))) + (this->KyFactor * ((double) (k1y - k2y))));
}

// compute the matrix element for the two body interaction between two sites A2 and A4
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementA2A4()
{
  return 1.0;
}

// compute the matrix element for the two body interaction between two sites A2 and A6
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementA2A6()
{
  return 1.0;
}

// compute the matrix element for the two body interaction between two sites A3 and A4 
//
// k1x = creation momentum along x for the A4 site
// k1y = creation momentum along y for the A4 site
// k2x = annihilation momentum along x for the A4 site
// k2y = annihilation momentum along y for the A4 site
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementA3A4(int k1x, int k1y, int k2x, int k2y)
{
  return Phase(this->KxFactor * ((double) (k2x - k1x)));
}

// compute the matrix element for the two body interaction between two sites A3 and A5
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementA3A5()
{
  return 1.0;
}

// compute the matrix element for the two body interaction between two sites A4 and A5
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementA4A5()
{
  return 1.0;
}

// compute the matrix element for the two body interaction between two sites A4 and A6
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementA4A6()
{
  return 1.0;
}

// compute the matrix element for the two body interaction between two sites A5 and A6 
//
// k1x = creation momentum along x for the A6 site
// k1y = creation momentum along y for the A6 site
// k2x = annihilation momentum along x for the A6 site
// k2y = annihilation momentum along y for the A6 site
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementA5A6(int k1x, int k1y, int k2x, int k2y)
{
  return Phase((this->KxFactor * ((double) (k2x - k1x))) + (this->KyFactor * ((double) (k2y - k1y))));
}

// compute the matrix element for the two body interaction between two sites A and C 
//
// k1a = creation momentum along x for the C site
// k1b = creation momentum along y for the C site
// k2a = annihilation momentum along x for the C site
// k2b = annihilation momentum along y for the C site
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementAC(int k1a, int k1b, int k2a, int k2b)
{
  Complex Tmp = 2.0 * cos (0.5 * (this->KyFactor * ((double) (k2b - k1b))));
  //Complex Tmp = Phase (0.5 * (this->KyFactor * ((double) (k2b - k1b))));
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

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementBC(int k1a, int k1b, int k2a, int k2b, int k3a, int k3b, int k4a, int k4b)
{
  Complex Tmp = 2.0 * cos (0.5 * ((this->KxFactor * ((double) (k3a - k1a))) + (this->KyFactor * ((double) (k4b - k2b)))));
  //Complex Tmp = Phase(0.5 * ((this->KxFactor * ((double) (k3a - k1a))) + (this->KyFactor * ((double) (k4b - k2b)))));
  return Tmp;
}


// compute the matrix element for on-site two body interaction involving A1 sites
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteA1A1()
{
  return 1.0;
}

// compute the matrix element for on-site two body interaction involving A2 sites
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteA2A2()
{
  return 1.0;
}

// compute the matrix element for on-site two body interaction involving A3 sites
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteA3A3()
{
  return 1.0;
}

// compute the matrix element for on-site two body interaction involving A4 sites
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteA4A4()
{
  return 1.0;
}

// compute the matrix element for on-site two body interaction involving A5 sites
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteA5A5()
{
  return 1.0;
}

// compute the matrix element for on-site two body interaction involving A6 sites
//
// return value = corresponding matrix element

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteA6A6()
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

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteBB(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return Phase(0.5 * this->KxFactor * ((double) (kx4 + kx3 - kx2 -kx1)));
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

Complex ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeTwoBodyMatrixElementOnSiteCC(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return Phase(0.5 * this->KyFactor * ((double) (ky4 + ky3 - ky2 -ky1)));
}

// compute the one body transformation matrices and the optional one body band stucture contribution
//
// oneBodyBasis = array of one body transformation matrices

void ParticleOnLatticeRubyLatticeSingleBandHamiltonian::ComputeOneBodyMatrices(ComplexMatrix* oneBodyBasis)
{
  Complex CT (this->TrHopping, this->TiHopping);
  Complex CT1 (this->T1rHopping, this->T1iHopping);
  double KX;
  double KY;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  KX = this->KxFactor * (((double) kx) + this->GammaX);
	  KY = this->KyFactor * (((double) ky) + this->GammaY);
	  int Index = (kx * this->NbrSiteY) + ky;
	  
	  Complex PhaseX = Phase(KX);
	  Complex PhaseXY = Phase(KX + KY);
	  
	  HermitianMatrix TmpOneBodyHamiltonian(6, true);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 2, Conj(CT));
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 4, CT);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 4, Conj(CT));
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 3, Conj(CT));
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 5, CT);
	  TmpOneBodyHamiltonian.SetMatrixElement(3, 5, Conj(CT));

	  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, CT1);
	  TmpOneBodyHamiltonian.SetMatrixElement(3, 4, CT1);
	  TmpOneBodyHamiltonian.SetMatrixElement(0, 5, Conj(CT1) * Conj(PhaseX));
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 2, CT1 * PhaseXY);
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 3, CT1 * Conj(PhaseX));	
	  TmpOneBodyHamiltonian.SetMatrixElement(4, 5, CT1 * Conj(PhaseXY));

	  TmpOneBodyHamiltonian.SetMatrixElement(0, 3, this->T4Hopping * (1.0 + Conj(PhaseX)));
	  TmpOneBodyHamiltonian.SetMatrixElement(1, 4, this->T4Hopping * (1.0 + PhaseXY));
	  TmpOneBodyHamiltonian.SetMatrixElement(2, 5, this->T4Hopping * (Conj(PhaseX) + Conj(PhaseXY)));

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
