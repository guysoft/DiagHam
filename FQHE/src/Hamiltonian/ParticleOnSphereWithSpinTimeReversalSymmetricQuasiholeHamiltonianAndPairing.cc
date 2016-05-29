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
  
  this->ChargingEnergy = chargingEnergy;
  this->AverageNumberParticles = averageNumberParticles;
  this->HamiltonianShift =  0.0;
  long MinIndex;
  long MaxIndex;
  
  
  
  this->Architecture = architecture;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->Memory = memory;
  this->FastMultiplicationFlag = false;
  this->DiskStorageFlag = onDiskCacheFlag;

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
  
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      if (TmpMemory < 1024l)
	cout  << "fast = " <<  TmpMemory << "b ";
      else
	if (TmpMemory < (1l << 20))
	  cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	else
	  if (TmpMemory < (1l << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
      cout << endl;
      this->EnableFastMultiplication();
    }
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
  int LastComponent = firstComponent + nbrComponent;
  double TmpCoefficient;
  
  int MaximalNumberCouplingElements = ((QuasiholeOnSphereWithSpinAndPairing*) (this->Particles))->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  if (this->FastMultiplicationFlag == false)
    {
      cout << "recomputing hamiltonian elements at each step" << endl;
      QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
      int index = firstComponent;
      int NbrElements;
      double ChargeContribution;
      double TmpCoef;
      
      for (int i = 0; i < nbrComponent; ++i)
	{
	  TmpCoef = vSource[index];
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if (TmpCoefficient != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
		      NbrElements = TmpParticles->AduAdd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
		    }	
		}
	    }
	      
	  if (this->OneBodyInteractionFactorsupup != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
	      {
		NbrElements = TmpParticles->AduAu(index, lz, TmpLeftIndices, TmpInteractionElements);
		for (int j = 0; j < NbrElements; ++j)
		  vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
	      }
	    }
	  }
	      
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
	      {
		NbrElements = TmpParticles->AddAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		for (int j = 0; j < NbrElements; ++j)
		  vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
	      }
	    }
	  }
	  
	  if (this->ChargingEnergy != 0.0)
	  {
	    int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	    if (TmpTotalNbrParticles != this->AverageNumberParticles)
	    {
	      ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	      ChargeContribution *= ChargeContribution;
	      ChargeContribution *= (vSource[index] * this->ChargingEnergy);
	      vDestination[index] += ChargeContribution;
	    }
	  }
	  
	  ++index;
	}
      
      
      if (this->HamiltonianShift != 0.0)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    vDestination[i] += this->HamiltonianShift * vSource[i];
	}
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      TmpCoefficient = vSource[k];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * TmpCoefficient;
	      vDestination[k++] += this->HamiltonianShift * TmpCoefficient;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      cout << "Error: ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelAddMultiplyDiskStorage not implemented" << endl;
// 	      this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  
  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
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
  int LastComponent = firstComponent + nbrComponent;
  int MaximalNumberCouplingElements = ((QuasiholeOnSphereWithSpinAndPairing*) (this->Particles))->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  if (this->FastMultiplicationFlag == false)
    {
      QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
      int index = firstComponent;
      int NbrElements;
      double ChargeContribution;
      double TmpOneBodyInteraction;
      for (int i = 0; i < nbrComponent; ++i)
	{
	  for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsPairing != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsPairing[lz];    
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			  vDestinations[k][TmpLeftIndices[j]] += (vSources[k][index] * TmpInteractionElements[j] * TmpOneBodyInteraction);
		      NbrElements = TmpParticles->AduAdd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			  vDestinations[k][TmpLeftIndices[j]] += (vSources[k][index] * TmpInteractionElements[j] * TmpOneBodyInteraction);
		    }
		}
	      
	      if (this->OneBodyInteractionFactorsupup != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsupup[lz];
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AduAu(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			  vDestinations[k][TmpLeftIndices[j]] += (vSources[k][index] * TmpInteractionElements[j] * TmpOneBodyInteraction);
		    }
		}
	      
	      if (this->OneBodyInteractionFactorsdowndown != 0)
		{
		  TmpOneBodyInteraction = this->OneBodyInteractionFactorsdowndown[lz];
		  if (TmpOneBodyInteraction != 0.0)
		    {
		      NbrElements = TmpParticles->AddAd(index, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int k = 0; k < nbrVectors; ++k)
			  vDestinations[k][TmpLeftIndices[j]] += (vSources[k][index] * TmpInteractionElements[j] * TmpOneBodyInteraction);
		    }
		}
	    }
	    
	  if (this->ChargingEnergy != 0.0)
	  {
	    int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	    if (TmpTotalNbrParticles != this->AverageNumberParticles)
	    {
	      ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	      ChargeContribution *= ChargeContribution;
	      ChargeContribution *= this->ChargingEnergy;
	      for (int k = 0; k < nbrVectors; ++k)
		vDestinations[k][index] += (ChargeContribution * vSources[k][index]);
	    }
	  }
	  
	  ++index;
	}
      
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
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  double* Coefficient2 = new double [nbrVectors];
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int j;
	  int Pos;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  double Coefficient;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  Coefficient2[l] = vSources[l][k];
		  vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
		}
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	}
	else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      cout << "Error: ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelMultipleAddMultiplyDiskStorage not implemented" << endl;
// 	      this->LowLevelMultipleAddMultiplyDiskStorage(vSource, vDestination, nbrVectors, firstComponent, nbrComponent);
	    }
	}
    }
    
  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  return vDestinations;
}


// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::GetHilbertSpace ()
{
  return this->Particles;
}


// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that has to be precalcualted
// return value = number of non-zero matrix element

long ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::PartialFastMultiplicationMemory(int firstComponent, int nbrComponent)
{
  int NbrElements;
  double TmpCoefficient;
  long Memory = 0l;
  
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  int LastComponent = nbrComponent + firstComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
  {
    if (this->OneBodyInteractionFactorsPairing != 0)
    {
      for (int lz = 0; lz <= this->LzMax; ++lz)
	{
	  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
	  if (TmpCoefficient != 0.0)
	    {
	      NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
	      Memory += NbrElements;
	      this->NbrInteractionPerComponent[i - this->PrecalculationShift] += NbrElements;
	      
	      NbrElements = TmpParticles->AduAdd(i, lz, TmpLeftIndices, TmpInteractionElements);
	      Memory += NbrElements;
	      this->NbrInteractionPerComponent[i - this->PrecalculationShift] += NbrElements;
	    }	
	}
    }
    if (this->OneBodyInteractionFactorsupup != 0)
    {
      for (int lz = 0; lz <= this->LzMax; ++lz)
      {
	if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
	{
	  NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
	  Memory += NbrElements;
	  this->NbrInteractionPerComponent[i - this->PrecalculationShift] += NbrElements;
	}
      }
    }      
    if (this->OneBodyInteractionFactorsdowndown != 0)
      for (int lz = 0; lz <= this->LzMax; ++lz)
	if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
	{
	  NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
	  Memory += NbrElements;
	  this->NbrInteractionPerComponent[i - this->PrecalculationShift] += NbrElements;
	}
    
    int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
    if ((this->ChargingEnergy != 0.0) && (TmpTotalNbrParticles != this->AverageNumberParticles))
    {
      ++Memory;
      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
    }  
  }
  
  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  
  delete TmpParticles;
  return Memory;
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = index of the last component that has to be precalcualted

void ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{  
  double TmpCoefficient;
  double ChargeContribution;
  int NbrElements;
  int LastComponent = nbrComponent + firstComponent;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();

  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  
  long Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
    
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      int position = 0;
      if (this->OneBodyInteractionFactorsPairing != 0)
      {
	for (int lz = 0; lz <= this->LzMax; ++lz)
	{
	  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
	  if (TmpCoefficient != 0.0)
	    {
	      NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
	      for (int j = 0; j < NbrElements; ++j)
	      {
		this->InteractionPerComponentIndex[Pos][position] = TmpLeftIndices[j];
		this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficient * TmpInteractionElements[j];
		++position;
	      }
	      NbrElements = TmpParticles->AduAdd(i, lz, TmpLeftIndices, TmpInteractionElements);
	      for (int j = 0; j < NbrElements; ++j)
	      {
		this->InteractionPerComponentIndex[Pos][position] = TmpLeftIndices[j];
		this->InteractionPerComponentCoefficient[Pos][position] = TmpCoefficient * TmpInteractionElements[j];
		++position;
	      }
	    }
	}
      }
    if (this->OneBodyInteractionFactorsupup != 0)
    {
      for (int lz = 0; lz <= this->LzMax; ++lz)
      {
	if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
	{
	  NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
	  for (int j = 0; j < NbrElements; ++j)
	  {
	    this->InteractionPerComponentIndex[Pos][position] = TmpLeftIndices[j];
	    this->InteractionPerComponentCoefficient[Pos][position] = this->OneBodyInteractionFactorsupup[lz] * TmpInteractionElements[j];
	    ++position;
	  }
	}
      }
    }      
    if (this->OneBodyInteractionFactorsdowndown != 0)
      for (int lz = 0; lz <= this->LzMax; ++lz)
	if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
	{
	  NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
	  for (int j = 0; j < NbrElements; ++j)
	  {
	    this->InteractionPerComponentIndex[Pos][position] = TmpLeftIndices[j];
	    this->InteractionPerComponentCoefficient[Pos][position] = this->OneBodyInteractionFactorsdowndown[lz] * TmpInteractionElements[j];
	    ++position;
	  }
	}
    
    if (this->ChargingEnergy != 0.0)
    {
      int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
      if (TmpTotalNbrParticles != this->AverageNumberParticles)
      {
	ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	ChargeContribution *= ChargeContribution;
	ChargeContribution *= this->ChargingEnergy;
	this->InteractionPerComponentIndex[Pos][position] = i;
	this->InteractionPerComponentCoefficient[Pos][position] = ChargeContribution;
	++position;
      }
    }
    ++Pos;
   }
  
  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  
  delete TmpParticles;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  double TmpCoefficient;
  double TmpCoef;
  double ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  int l =  PosMod + firstComponent + this->PrecalculationShift;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
      Coefficient = vSource[l];
      for (j = 0; j < TmpNbrInteraction; ++j)
	vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
      vDestination[l] += this->HamiltonianShift * Coefficient;
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  int Index;  
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	{
	  TmpCoef = vSource[i];
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if (TmpCoefficient != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
		      NbrElements = TmpParticles->AduAdd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * TmpCoefficient);
		    }	
		}
	    }
	    
	  if (this->OneBodyInteractionFactorsupup != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	      {
		if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		  {
		    NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
		    for (int j = 0; j < NbrElements; ++j)
		      vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
		  }
	      }
	  }
	      
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
	      {
		NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		for (int j = 0; j < NbrElements; ++j)
		  vDestination[TmpLeftIndices[j]] += (TmpCoef * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
		}
	    }
	  }
	  
	  if (this->ChargingEnergy != 0.0)
	  {
	    int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	    if (TmpTotalNbrParticles != this->AverageNumberParticles)
	    {
	      ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	      ChargeContribution *= ChargeContribution;
	      ChargeContribution *= (TmpCoef * this->ChargingEnergy);
	      vDestination[i] += ChargeContribution;
	    }
	  }
      }
      }

  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  delete TmpParticles;
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* ParticleOnSphereWithSpinTimeReversalSymmetricQuasiholeHamiltonianAndPairing::LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double* Coefficient2 = new double [nbrVectors];
  
  double Coefficient;
  double TmpCoefficient;
  double TmpCoef;
  double ChargeContribution;
  int NbrElements;
  int TmpTotalNbrParticles;
  int Pos2;
  QuasiholeOnSphereWithSpinAndPairing* TmpParticles = (QuasiholeOnSphereWithSpinAndPairing*) this->Particles->Clone();
  int MaximalNumberCouplingElements = TmpParticles->GetMaximalNumberCouplingElements();
  int* TmpLeftIndices = new int [MaximalNumberCouplingElements];
  double* TmpInteractionElements = new double[MaximalNumberCouplingElements];
  
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  int l =  PosMod + firstComponent + this->PrecalculationShift;
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
      for (int k = 0; k < nbrVectors; ++k)
	{
	  Coefficient2[k] = vSources[k][l];
	  vDestinations[k][l] += this->HamiltonianShift * Coefficient2[k];
	}
      for (j = 0; j < TmpNbrInteraction; ++j)
	{
	  Pos2 = TmpIndexArray[j];
	  Coefficient = TmpCoefficientArray[j];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][Pos2] += Coefficient  * Coefficient2[k];
	    }
	}
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  int Index;  
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	{
	  if (this->OneBodyInteractionFactorsPairing != 0)
	    {
	      for (int lz = 0; lz <= this->LzMax; ++lz)
		{
		  TmpCoefficient = this->OneBodyInteractionFactorsPairing[lz];
		  if (TmpCoefficient != 0.0)
		    {
		      NbrElements = TmpParticles->AuAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][TmpLeftIndices[j]] += (vSources[l][i] * TmpInteractionElements[j] * TmpCoefficient);
		      NbrElements = TmpParticles->AduAdd(i, lz, TmpLeftIndices, TmpInteractionElements);
		      for (int j = 0; j < NbrElements; ++j)
			for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][TmpLeftIndices[j]] += (vSources[l][i] * TmpInteractionElements[j] * TmpCoefficient);
		    }	
		}
	    }
	    
	  if (this->OneBodyInteractionFactorsupup != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	      {
		if (this->OneBodyInteractionFactorsupup[lz] != 0.0)
		  {
		    NbrElements = TmpParticles->AduAu(i, lz, TmpLeftIndices, TmpInteractionElements);
		    for (int j = 0; j < NbrElements; ++j)
		      for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][TmpLeftIndices[j]] += (vSources[l][i] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsupup[lz]);
		  }
	      }
	  }
	      
	  if (this->OneBodyInteractionFactorsdowndown != 0)
	  {
	    for (int lz = 0; lz <= this->LzMax; ++lz)
	    {
	      if (this->OneBodyInteractionFactorsdowndown[lz] != 0.0)
	      {
		NbrElements = TmpParticles->AddAd(i, lz, TmpLeftIndices, TmpInteractionElements);
		for (int j = 0; j < NbrElements; ++j)
		  for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][TmpLeftIndices[j]] += (vSources[l][i] * TmpInteractionElements[j] * this->OneBodyInteractionFactorsdowndown[lz]);
		}
	    }
	  }
	  
	  if (this->ChargingEnergy != 0.0)
	  {
	    int TmpTotalNbrParticles = TmpParticles->GetTotalNumberOfParticles(i);
	    if (TmpTotalNbrParticles != this->AverageNumberParticles)
	    {
	      ChargeContribution = (double) (TmpTotalNbrParticles - this->AverageNumberParticles);
	      ChargeContribution *= ChargeContribution;
	      ChargeContribution *= (TmpCoef * this->ChargingEnergy);
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += ChargeContribution;
	    }
	  }
      }
      }

  delete[] Coefficient2;
  delete[] TmpLeftIndices;
  delete[] TmpInteractionElements;
  delete TmpParticles;
  return vDestinations;
}