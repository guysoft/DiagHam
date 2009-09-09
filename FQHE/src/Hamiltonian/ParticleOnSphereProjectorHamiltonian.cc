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

ParticleOnSphereProjectorHamiltonian::ParticleOnSphereProjectorHamiltonian((ParticleOnSphere* particles, int nbrParticles, int lzmax, 
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
  delete[] this->FourBodyPseudoPotential;

  delete[] this->NBodyFlags;
  delete[] this->NBodyInteractionFactors;
  delete[] this->SortedIndicesPerSum;
  delete[] this->NbrSortedIndicesPerSum;
  delete[] this->MinSumIndices;
  delete[] this->MaxSumIndices;
  delete[] this->NBodySign;
  if (this->L2Operator != 0)
    delete this->L2Operator;
  if (this->FullTwoBodyFlag == true)
    {
      delete[] this->InteractionFactors;
      delete[] this->M1Value;
      delete[] this->M2Value;
      delete[] this->M3Value;
      delete[] this->PseudoPotential;
    }
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
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      FermionOnSphere* TmpSpace = (FermionOnSphere*) this->ProjectorSpace;
      TmpSpace->GetMonomial(0, TmpMonomial);
      this->MaxSumIndices[this->MaxNBody] = 0;
      for (int i = 0; i < this->MaxNBody; ++i)
	this->MaxSumIndices[this->MaxNBody] += TmpMonomial[i];
      this->MinSumIndices[this->MaxNBody] = this->MaxSumIndices[this->MaxNBody];
      
      for (int MinSum = 0; MinSum < this->MaxSumIndices[this->MaxNBody]; ++MinSum)
	this->NbrSortedIndicesPerSum[this->MaxNBody][MinSum] = 0;
    }

  this->NbrSortedIndicesPerSum[this->MaxNBody][this->ProjectorSpaceTotalLz] = this->ProjectorSpace->GetHilbertSpaceDimension;
  this->NbrNIndices[this->MaxNBody] = this->ProjectorSpace->GetHilbertSpaceDimension;
  this->NIndices[this->MaxNBody] = new int[this->ProjectorSpace->GetHilbertSpaceDimension * this->MaxNBody];
  this->NbrMIndices[this->MaxNBody] = new long[this->ProjectorSpace->GetHilbertSpaceDimension];
  this->MIndices[this->MaxNBody] = new int*[this->ProjectorSpace->GetHilbertSpaceDimension];
  this->MNNBodyInteractionFactors[this->MaxNBody] = new double* [this->ProjectorSpace->GetHilbertSpaceDimension];
  TmpNbrNIndices = 0;	 
  int* TmpNIndices = this->NIndices[this->MaxNBody];
  for (int MinSum = 0; MinSum <= this->MaxSumIndices[this->MaxNBody]; ++MinSum)
    {
      int Lim = this->NbrSortedIndicesPerSum[this->MaxNBody][MinSum];
      if (Lim > 0)
	{
	  int* TmpNIndices2 = this->SortedIndicesPerSum[this->MaxNBody][MinSum];
	  int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2] + TmpNIndices2[3];
	  while (((this->MaxNBody * this->LzMax) - TmpMaxRealtiveMonentum)  < TmpSum)
	    --TmpMaxRealtiveMonentum;
	  TmpProjectorCoefficients[this->MaxNBody] = this->ComputeProjectorCoefficients(12, 1, TmpNIndices2, Lim);
	  for (int i = 5; i <= TmpMaxRealtiveMonentum; ++i)  
	    if (this->FourBodyPseudoPotential[i] != 0.0)
	      TmpProjectorCoefficients[i] = this->ComputeProjectorCoefficients(2 * i, 1, TmpNIndices2, Lim);
	  for (int i = 0; i < Lim; ++i)
	    {
	      this->NbrMIndices[this->MaxNBody][TmpNbrNIndices] = Lim;		    
	      this->MIndices[this->MaxNBody][TmpNbrNIndices] = new int [Lim * 4];
	      this->MNNBodyInteractionFactors[this->MaxNBody][TmpNbrNIndices] = new double [Lim];
	      int* TmpMIndices = this->MIndices[this->MaxNBody][TmpNbrNIndices];
	      int* TmpMIndices2 = this->SortedIndicesPerSum[this->MaxNBody][MinSum];
	      double* TmpInteraction = this->MNNBodyInteractionFactors[this->MaxNBody][TmpNbrNIndices];
	      for (int j = 0; j < Lim; ++j)
		{
		  for (int l = 0; l < this->MaxNBody; ++l)
		    {
		      (*TmpMIndices) = (*TmpMIndices2);			
		      ++TmpMIndices;
		      ++TmpMIndices2;
		    }			
		  double& TmpInteraction2 = TmpInteraction[j];
		  TmpInteraction2 = 0.0;
		  if (this->FourBodyPseudoPotential[this->MaxNBody] != 0.0)
		    TmpInteraction2 += this->FourBodyPseudoPotential[this->MaxNBody] * TmpProjectorCoefficients[this->MaxNBody][i] * TmpProjectorCoefficients[this->MaxNBody][j];
		  for (int k = 5; k <= TmpMaxRealtiveMonentum; ++k)  
		    if (this->FourBodyPseudoPotential[k] != 0.0)
		      TmpInteraction2 += this->FourBodyPseudoPotential[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
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
      delete[] TmpInteractionCoeffients;
    }
  else
    {
      this->MinSumIndices[this->MaxNBody] = 0;
      this->MaxSumIndices[this->MaxNBody] = this->LzMax * this->MaxNBody;
      double* TmpInteractionCoeffients = new double[this->MaxSumIndices[this->MaxNBody] + 1];
      double Coefficient;
      TmpInteractionCoeffients[0] = 1.0;
      TmpInteractionCoeffients[1] = 1.0;
      for (int i = 2; i <= this->MaxSumIndices[this->MaxNBody]; ++i)
	{
	  Coefficient = 1.0;
	  for (int j = 1; j < i; ++j)
	    {
	      double Coefficient2 = TmpInteractionCoeffients[j];
	      TmpInteractionCoeffients[j] += Coefficient;
	      Coefficient = Coefficient2;
	    }
	  TmpInteractionCoeffients[i] = 1.0;
	}
      Coefficient = this->MaxNBody.0 * M_PI / (((double) this->MaxSumIndices[this->MaxNBody]) + 1.0);
      double Radius = 2.0 / ((double) this->LzMax);
      for (int i = 2; i <= this->MaxNBody; ++i)
	{
	  Coefficient *= (double) (i * i);	  
	  Coefficient *= Radius;
	}
      for (int i = 0; i <= this->MaxSumIndices[this->MaxNBody]; ++i)
	TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);
      
      double** SortedIndicesPerSumSymmetryFactor;
      GetAllSymmetricIndices(this->NbrLzValue, this->MaxNBody, this->NbrSortedIndicesPerSum[this->MaxNBody], this->SortedIndicesPerSum[this->MaxNBody],
			     SortedIndicesPerSumSymmetryFactor);
      
      
      long TmpNbrNIndices = 0;
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[this->MaxNBody]; ++MinSum)
	TmpNbrNIndices += this->NbrSortedIndicesPerSum[this->MaxNBody][MinSum];
      this->NbrNIndices[this->MaxNBody] = TmpNbrNIndices;
      this->NIndices[this->MaxNBody] = new int[TmpNbrNIndices * this->MaxNBody];
      this->NbrMIndices[this->MaxNBody] = new long[TmpNbrNIndices];
      this->MIndices[this->MaxNBody] = new int*[TmpNbrNIndices];
      this->MNNBodyInteractionFactors[this->MaxNBody] = new double* [TmpNbrNIndices];
      TmpNbrNIndices = 0;	 
      int* TmpNIndices = this->NIndices[this->MaxNBody];
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[this->MaxNBody]; ++MinSum)
	{
	  int Lim = this->NbrSortedIndicesPerSum[this->MaxNBody][MinSum];
	  double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
	  int* TmpNIndices2 = this->SortedIndicesPerSum[this->MaxNBody][MinSum];
	  int TmpMaxRealtiveMonentum = 5;
	  if (this->MaxRelativeAngularMomentum <= TmpMaxRealtiveMonentum)
	    TmpMaxRealtiveMonentum = this->MaxRelativeAngularMomentum;
	  int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2] + TmpNIndices2[3];
	  while (((this->MaxNBody * this->LzMax) - TmpMaxRealtiveMonentum)  < TmpSum)
	    --TmpMaxRealtiveMonentum;
	  double** TmpProjectorCoefficients = new double* [TmpMaxRealtiveMonentum + 1];
	  if (this->FourBodyPseudoPotential[0] != 0.0)
	    TmpProjectorCoefficients[0] = this->ComputeProjectorCoefficients(0, 1, TmpNIndices2, Lim);
	  for (int i = 1; i <= TmpMaxRealtiveMonentum; ++i)  
	    if (this->FourBodyPseudoPotential[i] != 0.0)
	      TmpProjectorCoefficients[i] = this->ComputeProjectorCoefficients(2 * i, 1, TmpNIndices2, Lim);
	  for (int i = 0; i < Lim; ++i)
	    {
	      this->NbrMIndices[this->MaxNBody][TmpNbrNIndices] = Lim;		    
	      this->MIndices[this->MaxNBody][TmpNbrNIndices] = new int [Lim * this->MaxNBody];
	      this->MNNBodyInteractionFactors[this->MaxNBody][TmpNbrNIndices] = new double [Lim];
	      int* TmpMIndices = this->MIndices[this->MaxNBody][TmpNbrNIndices];
	      int* TmpMIndices2 = this->SortedIndicesPerSum[this->MaxNBody][MinSum];
	      double* TmpInteraction = this->MNNBodyInteractionFactors[this->MaxNBody][TmpNbrNIndices];
	      for (int j = 0; j < Lim; ++j)
		{
		  double Coefficient2 = TmpSymmetryFactors[j];
		  for (int l = 0; l < this->MaxNBody; ++l)
		    {
		      Coefficient2 *= TmpNormalizationCoeffients[(*TmpMIndices2)];		    
		      (*TmpMIndices) = (*TmpMIndices2);			
		      ++TmpMIndices;
		      ++TmpMIndices2;
		    }			
		  double& TmpInteraction2 = TmpInteraction[j];
		  TmpInteraction2 = 0.0;
		  if (this->FourBodyPseudoPotential[0] != 0.0)
		    TmpInteraction2 += this->FourBodyPseudoPotential[0] * TmpProjectorCoefficients[0][i] * TmpProjectorCoefficients[0][j];
		  for (int k = 2; k <= TmpMaxRealtiveMonentum; ++k)  
		    if (this->FourBodyPseudoPotential[k] != 0.0)
		      TmpInteraction2 += this->FourBodyPseudoPotential[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
		  TmpInteraction2 *= TmpSymmetryFactors[i] * TmpSymmetryFactors[j];
// 		  if ((MinSum == 1) && (i == 0) && (j == 0))
// 		    {
// 		      cout << MinSum << " " << i << " " << j << " " << TmpInteraction2 << endl;
// 		      for (int k = 2; k <= TmpMaxRealtiveMonentum; ++k)  
// 			if (this->FourBodyPseudoPotential[k] != 0.0)		      
// 			  cout << TmpProjectorCoefficients[k][i] << " " << TmpProjectorCoefficients[k][j] << endl;
// 		      cout << TmpSymmetryFactors[i] << " " << TmpSymmetryFactors[j] << endl;
// 		    }
		}
	      for (int j = 0; j < this->MaxNBody; ++j)
		{
		  (*TmpNIndices) = (*TmpNIndices2);			
		  ++TmpNIndices;
		  ++TmpNIndices2;
		}
	      ++TmpNbrNIndices;
	    }		
	  if (this->FourBodyPseudoPotential[0] != 0.0)
	    delete[] TmpProjectorCoefficients[0];
	  for (int i = 2; i <= TmpMaxRealtiveMonentum; ++i)  
	    if (this->FourBodyPseudoPotential[i] != 0.0)
	      delete[] TmpProjectorCoefficients[i];
	  delete[] TmpProjectorCoefficients;		
	}
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[this->MaxNBody]; ++MinSum)
	{
	  delete[] SortedIndicesPerSumSymmetryFactor[MinSum];
	}
      delete[] SortedIndicesPerSumSymmetryFactor;
      delete[] TmpInteractionCoeffients;
    }
  delete[] TmpNormalizationCoeffients;
}

