////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                to particles on a torus with magnetic translations          //
//                                                                            //
//                        last modification : 18/02/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeWithSpinChernInsulatorHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>


using std::cout;
using std::endl;
using std::ostream;



// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kineticFactor = multiplicative factor in front of the kinetic term
// uPotential = Hubbard potential strength
// bandParameter = band parameter
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeWithSpinChernInsulatorHamiltonian::ParticleOnLatticeWithSpinChernInsulatorHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles, int nbrSiteX, 
												       int nbrSiteY, double kineticFactor, double uPotential, double bandParameter, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->UPotential = 2.0 * uPotential / ((double) this->NbrParticles);
  this->HamiltonianShift = 0.0;//4.0 * uPotential;
  this->KineticFactor = kineticFactor;
  this->BandParameter = bandParameter;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
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
}

// destructor
//

ParticleOnLatticeWithSpinChernInsulatorHamiltonian::~ParticleOnLatticeWithSpinChernInsulatorHamiltonian()
{
}
  
// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeWithSpinChernInsulatorHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (ParticleOnSphereWithSpin*) hilbertSpace;
  this->EvaluateInteractionFactors();
}
 
// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeWithSpinChernInsulatorHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeWithSpinChernInsulatorHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnLatticeWithSpinChernInsulatorHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnLatticeWithSpinChernInsulatorHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSource, vDestination);
      delete TmpParticles;
    }
  return vDestination;
}
 
// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* ParticleOnLatticeWithSpinChernInsulatorHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
											       int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  if (this->FastMultiplicationFlag == false)
    {
      Complex* Coefficient2 = new Complex [nbrVectors];
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
      for (int i = firstComponent; i < LastComponent; ++i)
	this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2);
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSources, vDestinations, nbrVectors);
      delete[] Coefficient2;
      delete TmpParticles;
    }
  return vDestinations;
}

// evaluate all interaction factors
//   

void ParticleOnLatticeWithSpinChernInsulatorHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
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
      this->InteractionFactorsupup = new double* [this->NbrIntraSectorSums];
      this->InteractionFactorsdowndown = new double* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  this->InteractionFactorsupup[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  this->InteractionFactorsdowndown[i] = new double[this->NbrIntraSectorIndicesPerSum[i] * this->NbrIntraSectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrIntraSectorIndicesPerSum[i]; ++j1)
	    {
	      int kx1 = this->IntraSectorIndicesPerSum[i][j1 << 1] / this->NbrSiteY;
	      int ky1 = this->IntraSectorIndicesPerSum[i][j1 << 1] % this->NbrSiteY;
	      int kx2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1] / this->NbrSiteY;
	      int ky2 = this->IntraSectorIndicesPerSum[i][(j1 << 1) + 1] % this->NbrSiteY;
	      for (int j2 = 0; j2 < this->NbrIntraSectorIndicesPerSum[i]; ++j2)
		{
		  int kx3 = this->IntraSectorIndicesPerSum[i][j2 << 1] / this->NbrSiteY;
		  int ky3 = this->IntraSectorIndicesPerSum[i][j2 << 1] % this->NbrSiteY;
		  int kx4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1] / this->NbrSiteY;
		  int ky4 = this->IntraSectorIndicesPerSum[i][(j2 << 1) + 1] % this->NbrSiteY;
		  this->InteractionFactorsupup[i][Index] = -4.0  * this->UPotential * ((cos (2.0 * M_PI * ((double) kx2 - kx4) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky2 - ky4) / ((double) this->NbrSiteY)))
								   - (cos (2.0 * M_PI * ((double) kx1 - kx4) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky1 - ky4) / ((double) this->NbrSiteY)))
								   - (cos (2.0 * M_PI * ((double) kx2 - kx3) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky2 - ky3) / ((double) this->NbrSiteY)))
								   + (cos (2.0 * M_PI * ((double) kx1 - kx3) / ((double) this->NbrSiteX)) + cos (2.0 * M_PI * ((double) ky1 - ky3) / ((double) this->NbrSiteY))));
		  this->InteractionFactorsdowndown[i][Index] = this->InteractionFactorsupup[i][Index];
		  TotalNbrInteractionFactors += 2;
		  ++Index;
		}
	    }
	}
      this->InteractionFactorsupdown = new double* [this->NbrInterSectorSums];
      for (int i = 0; i < this->NbrInterSectorSums; ++i)
	{
	  this->InteractionFactorsupdown[i] = new double[this->NbrInterSectorIndicesPerSum[i] * this->NbrInterSectorIndicesPerSum[i]];
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
  this->OneBodyInteractionFactorsupup = new double [this->LzMax + 1];
  this->OneBodyInteractionFactorsdowndown = new double [this->LzMax + 1];
  this->OneBodyInteractionFactorsupdown = new Complex [this->LzMax + 1];
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
	cout << TmpDiag(0, 0) << " " << TmpDiag(1, 1) << endl;
       	this->OneBodyInteractionFactorsupup[Index] = this->KineticFactor * ((TmpDiag(0, 0) * SqrNorm(TmpMatrix[0][0])) 
									    + (TmpDiag(1, 1) * SqrNorm(TmpMatrix[1][0])));
       	this->OneBodyInteractionFactorsdowndown[Index] = this->KineticFactor * ((TmpDiag(0, 0) * SqrNorm(TmpMatrix[0][1])) 
										+ (TmpDiag(1, 1) * SqrNorm(TmpMatrix[1][1])));
       	this->OneBodyInteractionFactorsupdown[Index] = this->KineticFactor * ((TmpDiag(0, 0) * TmpMatrix[0][0] * Conj(TmpMatrix[0][1])) 
									      + (TmpDiag(1, 1) * TmpMatrix[1][0] * Conj(TmpMatrix[1][1])));
      }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

