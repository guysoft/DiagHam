////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                      Class author: Cecile Repellin                         //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to Laughlin qh on a sphere with      //
//        SU(2) spin with opposite magnetic field for each species            //
//                            and pairing                                     //
//                                                                            //
//                        last modification : 11/05/2016                      //
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


#include "Hamiltonian/ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;

//default constructor
ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// lzmax = maximum Lz value reached by a particle in the state
// architecture = architecture to use for precalculation
// onebodyPotentialUpUp =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin up, null pointer if none
// onebodyPotentialDownDown =  one-body potential (sorted from component on the lowest Lz state to component on the highest Lz state) for particles with spin down, null pointer if none
// onebodyPotentialPairing =  one-body pairing term (sorted from component on the lowest Lz state to component on the highest Lz state), on site, symmetric spin up / spin down
// chargingEnergy = factor in front of the charging energy (i.e 1/(2C))
// averageNumberParticles = average number of particles in the system
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing(QuasiholeOnSphereWithSpinAndPairing* particles, int lzmax, double* onebodyPotentialUpUp, double* onebodyPotentialDownDown, double* onebodyPotentialPairing, double chargingEnergy, double averageNumberParticles, AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
//   this->NbrParticles = 0;
//   this->FastMultiplicationFlag = false;
//   this->OneBodyTermFlag = false;
//   this->Architecture = architecture;
//   this->L2Hamiltonian = 0;
//   this->S2Hamiltonian = 0;
  
  this->ChargingEnergy = chargingEnergy;
//   this->EvaluateInteractionFactors();
  this->AverageNumberParticles = averageNumberParticles;
  this->HamiltonianShift =  0.0;
  long MinIndex;
  long MaxIndex;
//   this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
//   this->PrecalculationShift = (int) MinIndex;  
//   this->DiskStorageFlag = onDiskCacheFlag;
//   this->Memory = memory;
  if (onebodyPotentialUpUp != 0)
    {
      this->OneBodyInteractionFactorsupup = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsupup[i] = onebodyPotentialUpUp[i];
    }
  
  if (onebodyPotentialDownDown != 0)
    {
      this->OneBodyInteractionFactorsdowndown = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	this->OneBodyInteractionFactorsdowndown[i] = onebodyPotentialDownDown[i];
    }
  this->OneBodyInteractionFactorsPairing = 0;
  if (onebodyPotentialPairing != 0)
    {
      this->OneBodyInteractionFactorsPairing = new double [this->NbrLzValue];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->OneBodyInteractionFactorsPairing[i] = onebodyPotentialPairing[i];
	}
    }
//   if (precalculationFileName == 0)
//     {
//       if (memory > 0)
// 	{
// 	  long TmpMemory = this->FastMultiplicationMemory(memory);
// 	  if (TmpMemory < 1024)
// 	    cout  << "fast = " <<  TmpMemory << "b ";
// 	  else
// 	    if (TmpMemory < (1 << 20))
// 	      cout  << "fast = " << (TmpMemory >> 10) << "kb ";
// 	    else
// 	  if (TmpMemory < (1 << 30))
// 	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
// 	  else
// 	    {
// 	      cout  << "fast = " << (TmpMemory >> 30) << ".";
// 	      TmpMemory -= ((TmpMemory >> 30) << 30);
// 	      TmpMemory *= 100l;
// 	      TmpMemory >>= 30;
// 	      if (TmpMemory < 10l)
// 		cout << "0";
// 	      cout  << TmpMemory << " Gb ";
// 	    }
// 	  if (this->DiskStorageFlag == false)
// 	    {
// 	      this->EnableFastMultiplication();
// 	    }
// 	  else
// 	    {
// 	      char* TmpFileName = this->Architecture->GetTemporaryFileName();
// 	      this->EnableFastMultiplicationWithDiskStorage(TmpFileName);	      
// 	      delete[] TmpFileName;
// 	    }
// 	}
//     }
//   else
//     this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::~ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing() 
{
  if (this->OneBodyInteractionFactorsupup != 0)
    delete[] this->OneBodyInteractionFactorsupup;
  if (this->OneBodyInteractionFactorsdowndown != 0)
    delete[] this->OneBodyInteractionFactorsdowndown;
  if (this->OneBodyInteractionFactorsPairing != 0)
    delete [] this->OneBodyInteractionFactorsPairing;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (QuasiholeOnSphereWithSpinAndPairing*) hilbertSpace;
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift += shift;
}
  
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
						     int firstComponent, int nbrComponent)
{
  int index = firstComponent;
  int LastComponent = firstComponent + nbrComponent;
  int NbrElements;
  int* LeftIndices = 0;
  double* Coefficients = 0;
  double ChargeContribution;
  for (int i = 0; i < nbrComponent; ++i)
  {
    double TmpCoef = vSource[index];
    for (int lz = -this->LzMax; lz <= this->LzMax; lz += 2)
    {
      int shiftedLz = (lz + this->LzMax) / 2;
      if (this->OneBodyInteractionFactorsPairing != 0)
      {
	NbrElements = this->Particles->AuAd(index, lz, LeftIndices, Coefficients);
	for (int j = 0; j < NbrElements; ++j)
	{
	  vDestination[LeftIndices[j]] += (TmpCoef * Coefficients[j] * this->OneBodyInteractionFactorsPairing[shiftedLz]);
	}
	NbrElements = this->Particles->AduAdd(index, lz, LeftIndices, Coefficients);
	for (int j = 0; j < NbrElements; ++j)
	{
	  vDestination[LeftIndices[j]] += (TmpCoef * Coefficients[j] * this->OneBodyInteractionFactorsPairing[shiftedLz]);
	}
      }
      
      if (this->OneBodyInteractionFactorsupup != 0)
      {
	NbrElements = this->Particles->AduAu(index, lz, LeftIndices, Coefficients);
	for (int j = 0; j < NbrElements; ++j)
	  vDestination[LeftIndices[j]] += (TmpCoef * Coefficients[j] * this->OneBodyInteractionFactorsupup[shiftedLz]);
      }
      
      if (this->OneBodyInteractionFactorsdowndown != 0)
      {
	NbrElements = this->Particles->AddAd(index, lz, LeftIndices, Coefficients);
	for (int j = 0; j < NbrElements; ++j)
	  vDestination[LeftIndices[j]] += (TmpCoef * Coefficients[j] * this->OneBodyInteractionFactorsdowndown[shiftedLz]);
      }
    }
    
    ChargeContribution = this->Particles->GetTotalNumberOfParticles(index) - this->AverageNumberParticles;
    ChargeContribution *= ChargeContribution;
    ChargeContribution *= (vSource[index] * this->ChargingEnergy);
    vDestination[index] += ChargeContribution;
    
    ++index;
  }
  
  if (LeftIndices != 0)
    delete[] LeftIndices;
  if (Coefficients != 0)
    delete[] Coefficients;
  if (this->HamiltonianShift != 0.0)
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	vDestination[i] += this->HamiltonianShift * vSource[i];
    }
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int index = firstComponent;
  int LastComponent = firstComponent + nbrComponent;
  int NbrElements;
  int* LeftIndices = 0;
  double* Coefficients = 0;
  double ChargeContribution;
  double TmpOneBodyInteraction;
  for (int i = 0; i < nbrComponent; ++i)
  {
    for (int lz = -this->LzMax; lz <= this->LzMax; lz += 2)
    {
      int shiftedLz = (lz + this->LzMax) / 2;  
      if (this->OneBodyInteractionFactorsPairing != 0)
      {
	TmpOneBodyInteraction = this->OneBodyInteractionFactorsPairing[shiftedLz];    
	NbrElements = this->Particles->AuAd(index, lz, LeftIndices, Coefficients);
	for (int j = 0; j < NbrElements; ++j)
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][LeftIndices[j]] += (vSources[k][index] * Coefficients[j] * TmpOneBodyInteraction);
	NbrElements = this->Particles->AduAdd(index, lz, LeftIndices, Coefficients);
	for (int j = 0; j < NbrElements; ++j)
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][LeftIndices[j]] += (vSources[k][index] * Coefficients[j] * TmpOneBodyInteraction);
      }
      
      if (this->OneBodyInteractionFactorsupup != 0)
      {
	TmpOneBodyInteraction = this->OneBodyInteractionFactorsupup[shiftedLz];
	NbrElements = this->Particles->AduAu(index, lz, LeftIndices, Coefficients);
	for (int j = 0; j < NbrElements; ++j)
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][LeftIndices[j]] += (vSources[k][index] * Coefficients[j] * TmpOneBodyInteraction);
      }
      
      if (this->OneBodyInteractionFactorsdowndown != 0)
      {
	TmpOneBodyInteraction = this->OneBodyInteractionFactorsdowndown[shiftedLz];
	NbrElements = this->Particles->AddAd(index, lz, LeftIndices, Coefficients);
	for (int j = 0; j < NbrElements; ++j)
	  for (int k = 0; k < nbrVectors; ++k)
	    vDestinations[k][LeftIndices[j]] += (vSources[k][index] * Coefficients[j] * TmpOneBodyInteraction);
      }
    }
    ChargeContribution = this->Particles->GetTotalNumberOfParticles(index) - this->AverageNumberParticles;
    ChargeContribution *= ChargeContribution;
    ChargeContribution *= this->ChargingEnergy;
    for (int k = 0; k < nbrVectors; ++k)
      vDestinations[k][index] += (ChargeContribution * vSources[k][index]);
    
    ++index;
  }
  
  if (LeftIndices != 0)
    delete[] LeftIndices;
  if (Coefficients != 0)
    delete[] Coefficients;
  
  if (this->HamiltonianShift != 0.0)
    {
      for (int k= 0; k < nbrVectors; ++k)
	{
	  RealVector& TmpDestination = vDestinations[k];
	  RealVector& TmpSource = vSources[k];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    TmpDestination[i] += this->HamiltonianShift * TmpSource[i];
	}
    }
  return vDestinations;
}


// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space
AbstractHilbertSpace* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::GetHilbertSpace ()
{
  return this->Particles;
}