////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                class of hamiltonian defined as a projector                 //
//                       for particles on a sphere defined                    //
//                                                                            //
//                        last modification : 14/08/2009                      //
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
#include "Hamiltonian/ParticleOnSphereProjectorHamiltonian.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"

  
#include <stdio.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// projectorState = state describing the projector
// projectorSpace = space associated to the projector state 
// projectorNbrParticles = number of particles for the projector state
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereProjectorHamiltonian::ParticleOnSphereProjectorHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
									   RealVector& projectorState, ParticleOnSphere* projectorSpace, int projectorNbrParticles, double l2Factor, 
									   AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
									   char* precalculationFileName )
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->OneBodyTermFlag = false;
  this->FullTwoBodyFlag = false;
  this->MaxNBody = projectorNbrParticles;
  this->ProjectorSpace = projectorSpace;
  this->ProjectorState = projectorState;

  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NBodyInteractionFactors = new double** [this->MaxNBody + 1];
  this->NbrSortedIndicesPerSum = new int* [this->MaxNBody + 1];
  this->SortedIndicesPerSum = new int** [this->MaxNBody + 1];
  this->MinSumIndices = new int [this->MaxNBody + 1];
  this->MaxSumIndices = new int [this->MaxNBody + 1];
  this->NBodySign = new double[this->MaxNBody + 1];

  this->NbrNIndices = new long[this->MaxNBody + 1];
  this->NIndices = new int*[this->MaxNBody + 1];
  this->NbrMIndices = new long*[this->MaxNBody + 1];
  this->MIndices = new int**[this->MaxNBody + 1];
  this->MNNBodyInteractionFactors = new double**[this->MaxNBody + 1];

  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->MinSumIndices[k] = 1;
      this->MaxSumIndices[k] = 0;      
      this->NBodyFlags[k] = false;
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic) && ((k & 1) == 0))
	{
	  this->NBodySign[k] = -1.0;
	}
      this->NbrNIndices[k] = 0;
      this->NIndices[k] = 0;
      this->NbrMIndices[k] = 0;
      this->MIndices[k] = 0;
      this->MNNBodyInteractionFactors[k] = 0;
    }

  this->NBodyFlags[this->MaxNBody] = true;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
  if (precalculationFileName == 0)
    {
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
	  if (this->DiskStorageFlag == false)
	    {
	      this->EnableFastMultiplication();
	    }
	  else
	    {
	      char* TmpFileName = this->Architecture->GetTemporaryFileName();
	      this->EnableFastMultiplicationWithDiskStorage(TmpFileName);	      
	      delete[] TmpFileName;
	    }
	}
      else
	{
	  this->FastMultiplicationFlag = false;
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);

  if (l2Factor != 0.0)
    {
      this->L2Operator = new ParticleOnSphereSquareTotalMomentumOperator(this->Particles, this->LzMax, l2Factor);
    }
  else
    {
      this->L2Operator = 0;
    }

}

// destructor
//

ParticleOnSphereProjectorHamiltonian::~ParticleOnSphereProjectorHamiltonian()
{
  for (int k = 1; k <= this->MaxNBody; ++k)
    if (this->NBodyFlags[k] == true)
      {
	for (int MinSum = this->MinSumIndices[k]; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	  {
	    delete[] this->SortedIndicesPerSum[k][MinSum];
	    if (this->MNNBodyInteractionFactors == 0)
	      delete[] this->NBodyInteractionFactors[k][MinSum];	      
	  }
	delete[] this->NbrSortedIndicesPerSum[k];
	delete[] this->SortedIndicesPerSum[k];
	if (this->MNNBodyInteractionFactors == 0)
	  delete[] this->NBodyInteractionFactors[k];
	else
	  {
	    for (int i = 0; i < this->NbrNIndices[k]; ++i)
	      {
		delete[] this->MNNBodyInteractionFactors[k][i];		
		delete[] this->MIndices[k][i];
	      }
	    delete[] this->NIndices[k];
	    delete[] this->NbrMIndices[k];
	    delete[] this->MIndices[k];
	    delete[] this->MNNBodyInteractionFactors[k];
	  }
      }

  delete[] this->NbrNIndices;
  delete[] this->NIndices;
  delete[] this->NbrMIndices;
  delete[] this->MIndices;
  delete[] this->MNNBodyInteractionFactors;

  delete[] this->NBodyFlags;
  delete[] this->NBodyInteractionFactors;
  delete[] this->SortedIndicesPerSum;
  delete[] this->NbrSortedIndicesPerSum;
  delete[] this->MinSumIndices;
  delete[] this->MaxSumIndices;
  delete[] this->NBodySign;
  if (this->L2Operator != 0)
    delete this->L2Operator;

  if (this->FastMultiplicationFlag == true)
    {
      if (this->DiskStorageFlag == false)
	{
	  long MinIndex;
	  long MaxIndex;
	  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
	  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
	  int ReducedDim = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if ((ReducedDim * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
	    ++ReducedDim;
	  for (int i = 0; i < ReducedDim; ++i)
	    {
	      delete[] this->InteractionPerComponentIndex[i];
	      delete[] this->InteractionPerComponentCoefficient[i];
	    }
	  delete[] this->InteractionPerComponentIndex;
	  delete[] this->InteractionPerComponentCoefficient;
	}
      else
	{
	  remove (this->DiskStorageFileName);
	  delete[] this->DiskStorageFileName;
	}
      delete[] this->NbrInteractionPerComponent;
    }
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* ParticleOnSphereProjectorHamiltonian::Clone ()
{
  return 0;
}


// evaluate all interaction factors
//   

void ParticleOnSphereProjectorHamiltonian::EvaluateInteractionFactors()
{
  int* TmpMonomial = new int[this->MaxNBody];
  int TmpSum = 0;
  double* IndexSymmetryFactor = new double [this->ProjectorSpace->GetHilbertSpaceDimension()];
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      FermionOnSphere* TmpSpace = (FermionOnSphere*) this->ProjectorSpace;
      TmpSpace->GetMonomial(0, TmpMonomial);
      this->MaxSumIndices[this->MaxNBody] = 0;
      for (int i = 0; i < this->MaxNBody; ++i)
	this->MaxSumIndices[this->MaxNBody] += TmpMonomial[i];
      this->MinSumIndices[this->MaxNBody] = this->MaxSumIndices[this->MaxNBody];
      this->SortedIndicesPerSum[this->MaxNBody] = new int*[TmpSum + 1];
      TmpSum = this->MaxSumIndices[this->MaxNBody];
      this->SortedIndicesPerSum[this->MaxNBody][TmpSum] = new int [this->ProjectorSpace->GetHilbertSpaceDimension() * this->MaxNBody];
      int* TmpSortedIndicesPerSum = this->SortedIndicesPerSum[this->MaxNBody][TmpSum];
      for (int i = 0; i < this->ProjectorSpace->GetHilbertSpaceDimension(); ++i)
	{
	  TmpSpace->GetMonomial(i, TmpMonomial);	  
	  for (int j = 0; j < this->MaxNBody; ++j)
	    {
	      (*TmpSortedIndicesPerSum) = TmpMonomial[j];
	      ++TmpSortedIndicesPerSum;
	    }
	  IndexSymmetryFactor[i] = 1.0;
	}
    }
  else  
    {
      double* SymmetryFactor = new double [this->MaxNBody + 1];
      SymmetryFactor[0] = 1.0;
      SymmetryFactor[1] = 1.0;
      for (int i = 2; i <= this->MaxNBody; ++i)
	SymmetryFactor[i] = SymmetryFactor[i - 1] / sqrt ((double) i);
      BosonOnSphereShort* TmpSpace = (BosonOnSphereShort*) this->ProjectorSpace;
      TmpSpace->GetMonomial(0, TmpMonomial);
      this->MaxSumIndices[this->MaxNBody] = 0;
      for (int i = 0; i < this->MaxNBody; ++i)
	TmpSum += TmpMonomial[i];
      this->MaxSumIndices[this->MaxNBody] = TmpSum;
      this->MinSumIndices[this->MaxNBody] = TmpSum;
      this->SortedIndicesPerSum[this->MaxNBody] = new int*[TmpSum + 1];
      this->SortedIndicesPerSum[this->MaxNBody][TmpSum] = new int [this->ProjectorSpace->GetHilbertSpaceDimension() * this->MaxNBody];
      int* TmpSortedIndicesPerSum = this->SortedIndicesPerSum[this->MaxNBody][TmpSum];
      for (int i = 0; i < this->ProjectorSpace->GetHilbertSpaceDimension(); ++i)
	{
	  TmpSpace->GetMonomial(i, TmpMonomial);	  
	  IndexSymmetryFactor[i] = 1.0;
	  int CurrentOccupation = 0;
	  int NbrOccupation = 0;
	  for (int j = 0; j < this->MaxNBody; ++j)
	    {
	      (*TmpSortedIndicesPerSum) = TmpMonomial[j];
	      ++TmpSortedIndicesPerSum;
	      if (CurrentOccupation != TmpMonomial[j])
		{
		  IndexSymmetryFactor[i] *= SymmetryFactor[NbrOccupation];
		  CurrentOccupation =  TmpMonomial[j];
		  NbrOccupation = 1;
		}
	      else
		++NbrOccupation;
	    }
	  IndexSymmetryFactor[i] *= SymmetryFactor[NbrOccupation];	  
	}
      delete[] SymmetryFactor;
    }

  this->NbrSortedIndicesPerSum[this->MaxNBody] = new int [TmpSum + 1];
  for (int MinSum = 0; MinSum < TmpSum; ++MinSum)
    {
      this->SortedIndicesPerSum[this->MaxNBody][MinSum] = 0;     
      this->NbrSortedIndicesPerSum[this->MaxNBody][MinSum] = 0;
    }
  this->NbrSortedIndicesPerSum[this->MaxNBody][TmpSum] = this->ProjectorSpace->GetHilbertSpaceDimension();

  this->NbrNIndices[this->MaxNBody] = this->ProjectorSpace->GetHilbertSpaceDimension();
  this->NIndices[this->MaxNBody] = new int[this->ProjectorSpace->GetHilbertSpaceDimension() * this->MaxNBody];
  this->NbrMIndices[this->MaxNBody] = new long[this->ProjectorSpace->GetHilbertSpaceDimension()];
  this->MIndices[this->MaxNBody] = new int*[this->ProjectorSpace->GetHilbertSpaceDimension()];
  this->MNNBodyInteractionFactors[this->MaxNBody] = new double* [this->ProjectorSpace->GetHilbertSpaceDimension()];

  int* TmpNIndices = this->NIndices[this->MaxNBody];
  int* TmpNIndices2 = this->SortedIndicesPerSum[this->MaxNBody][TmpSum];
  int TmpNbrNIndices = 0;	 
  for (int i = 0; i < this->ProjectorSpace->GetHilbertSpaceDimension(); ++i)
    {
      this->NbrMIndices[this->MaxNBody][TmpNbrNIndices] = this->ProjectorSpace->GetHilbertSpaceDimension();		    
      this->MIndices[this->MaxNBody][TmpNbrNIndices] = new int [this->ProjectorSpace->GetHilbertSpaceDimension() * this->MaxNBody];
      int* TmpMIndices = this->MIndices[this->MaxNBody][TmpNbrNIndices];
      int* TmpMIndices2 = this->SortedIndicesPerSum[this->MaxNBody][TmpSum];
      this->MNNBodyInteractionFactors[this->MaxNBody][TmpNbrNIndices] = new double [this->ProjectorSpace->GetHilbertSpaceDimension()];
      double* TmpInteraction = this->MNNBodyInteractionFactors[this->MaxNBody][TmpNbrNIndices];
      for (int j = 0; j < this->ProjectorSpace->GetHilbertSpaceDimension(); ++j)
	{
	  for (int l = 0; l < this->MaxNBody; ++l)
	    {
	      (*TmpMIndices) = (*TmpMIndices2);			
	      ++TmpMIndices;
	      ++TmpMIndices2;
	    }			
	  TmpInteraction[j] = this->ProjectorState[i] * this->ProjectorState[j] * IndexSymmetryFactor[i] * IndexSymmetryFactor[j];
	}
      for (int j = 0; j < this->MaxNBody; ++j)
	{
	  (*TmpNIndices) = (*TmpNIndices2);			
	  ++TmpNIndices;
	  ++TmpNIndices2;
	}
      ++TmpNbrNIndices;
    }


}

