////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                          generic 4-body interaction                        //
//                                                                            //
//                        last modification : 12/08/2009                      //
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
#include "Hamiltonian/ParticleOnSphereGenericFiveBodyHamiltonian.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"

  
#include <stdio.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


// default constructor
//

ParticleOnSphereGenericFiveBodyHamiltonian::ParticleOnSphereGenericFiveBodyHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// fiveBodyPseudoPotential = array with the five-body pseudo-potentials sorted with respect to the relative angular momentum, 
//                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
// maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in FiveBodyPseudoPotential
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereGenericFiveBodyHamiltonian::ParticleOnSphereGenericFiveBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
											 double* fiveBodyPseudoPotential, int maxRelativeAngularMomentum, double l2Factor, 
											 AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
											 char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->OneBodyTermFlag = false;
  this->FullTwoBodyFlag = false;
  this->MaxNBody = 4;

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

  this->MaxRelativeAngularMomentum = maxRelativeAngularMomentum;
  this->NbrFiveBodyPseudoPotential = maxRelativeAngularMomentum + 1;
  this->FiveBodyPseudoPotential = new double[this->NbrFiveBodyPseudoPotential];
  for (int i = 0; i < this->NbrFiveBodyPseudoPotential; ++i)
    this->FiveBodyPseudoPotential[i] = fiveBodyPseudoPotential[i];

  this->NBodyFlags[4] = true;
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

// constructor from datas with a fully-defined two body interaction
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// fiveBodyPseudoPotential = array with the five-body pseudo-potentials sorted with respect to the relative angular momentum, 
//                            taking into account of additional degeneracy for relative momentum greater than 5 for bosons (8 for fermions)
// maxRelativeAngularMomentum =  maxixmum relative angular momentum that is used in FiveBodyPseudoPotential
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// pseudoPotential = array with the pseudo-potentials (ordered such that the first element corresponds to the delta interaction)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereGenericFiveBodyHamiltonian::ParticleOnSphereGenericFiveBodyHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
											 double* fiveBodyPseudoPotential, int maxRelativeAngularMomentum,
											 double l2Factor, double* pseudoPotential, 
											 AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
											 char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  if (pseudoPotential != 0)
    {
      this->PseudoPotential = new double [this->NbrLzValue];
      for (int i = 0; i < this->NbrLzValue; ++i)
	this->PseudoPotential[i] = pseudoPotential[this->LzMax - i];
      this->FullTwoBodyFlag = true;
    }
  else
    {
      this->FullTwoBodyFlag = false;
    }

  this->OneBodyTermFlag = false;
  this->FullTwoBodyFlag = true;
  this->MaxNBody = 4;
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

  this->PseudoPotential = new double [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->PseudoPotential[i] = pseudoPotential[this->LzMax - i];

  this->MaxRelativeAngularMomentum = maxRelativeAngularMomentum;
  this->NbrFiveBodyPseudoPotential = maxRelativeAngularMomentum;
  this->FiveBodyPseudoPotential = new double[this->NbrFiveBodyPseudoPotential + 1];
  for (int i = 0; i <= this->NbrFiveBodyPseudoPotential; ++i)
    this->FiveBodyPseudoPotential[i] = fiveBodyPseudoPotential[i];


  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->MinSumIndices[k] = 1;
      this->MaxSumIndices[k] = 0;      
      this->NBodyFlags[k] = false;
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic) && ((k & 1) == 0))
	this->NBodySign[k] = -1.0;
    }
  this->NBodyFlags[4] = true;
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
	    cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
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

ParticleOnSphereGenericFiveBodyHamiltonian::~ParticleOnSphereGenericFiveBodyHamiltonian()
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
  delete[] this->FiveBodyPseudoPotential;

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
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* ParticleOnSphereGenericFiveBodyHamiltonian::Clone ()
{
  return 0;
}


// evaluate all interaction factors
//   

void ParticleOnSphereGenericFiveBodyHamiltonian::EvaluateInteractionFactors()
{
  double* TmpNormalizationCoeffients = new double[this->NbrLzValue];
  double TmpFactor = ((double) this->NbrLzValue) / (4.0 * M_PI);
  double TmpBinomial = 1.0;
  TmpNormalizationCoeffients[0] = sqrt (TmpBinomial * TmpFactor);
  for (int i = 1; i < this->NbrLzValue; ++i)
    {
      TmpBinomial *= this->LzMax - ((double) i) + 1.0;
      TmpBinomial /= ((double) i);
      TmpNormalizationCoeffients[i] = sqrt (TmpBinomial * TmpFactor);
    }
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      double Coefficient;
      GetAllSkewSymmetricIndices(this->NbrLzValue, 4, this->NbrSortedIndicesPerSum[4], this->SortedIndicesPerSum[4]);
      this->MaxSumIndices[4] = (this->LzMax * 4) - 5;
      this->MinSumIndices[4] = 6;
      double* TmpInteractionCoeffients = new double[this->MaxSumIndices[4] + 1];
      TmpInteractionCoeffients[0] = 1.0;
      TmpInteractionCoeffients[1] = 1.0;
      for (int i = 2; i <= this->MaxSumIndices[4]; ++i)
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
      Coefficient = 4.0 * M_PI / (((double) this->MaxSumIndices[4]) + 1.0);
      double Radius = 2.0 / ((double) this->LzMax);
      for (int i = 2; i <= 4; ++i)
	{
	  Coefficient *= (double) (i * i);	  
	  Coefficient *= Radius;
	}
      for (int i = 0; i <= this->MaxSumIndices[4]; ++i)
	TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);
      
      long TmpNbrNIndices = 0;
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[4]; ++MinSum)
	TmpNbrNIndices += this->NbrSortedIndicesPerSum[4][MinSum];
      this->NbrNIndices[4] = TmpNbrNIndices;
      this->NIndices[4] = new int[TmpNbrNIndices * 4];
      this->NbrMIndices[4] = new long[TmpNbrNIndices];
      this->MIndices[4] = new int*[TmpNbrNIndices];
      this->MNNBodyInteractionFactors[4] = new double* [TmpNbrNIndices];
      TmpNbrNIndices = 0;	 
      int* TmpNIndices = this->NIndices[4];
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[4]; ++MinSum)
	{
	  int Lim = this->NbrSortedIndicesPerSum[4][MinSum];
	  if (Lim > 0)
	    {
	      int* TmpNIndices2 = this->SortedIndicesPerSum[4][MinSum];
	      int TmpMaxRealtiveMonentum = 8;
	      if (this->MaxRelativeAngularMomentum <= TmpMaxRealtiveMonentum)
		TmpMaxRealtiveMonentum = this->MaxRelativeAngularMomentum;
	      int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2] + TmpNIndices2[3];
	      while (((4 * this->LzMax) - TmpMaxRealtiveMonentum)  < TmpSum)
		--TmpMaxRealtiveMonentum;
	      double** TmpProjectorCoefficients = new double* [TmpMaxRealtiveMonentum + 1];
	      if (this->FiveBodyPseudoPotential[4] != 0.0)
		TmpProjectorCoefficients[4] = this->ComputeProjectorCoefficients(12, 1, TmpNIndices2, Lim);
	      for (int i = 5; i <= TmpMaxRealtiveMonentum; ++i)  
		if (this->FiveBodyPseudoPotential[i] != 0.0)
		  TmpProjectorCoefficients[i] = this->ComputeProjectorCoefficients(2 * i, 1, TmpNIndices2, Lim);
	      for (int i = 0; i < Lim; ++i)
		{
		  this->NbrMIndices[4][TmpNbrNIndices] = Lim;		    
		  this->MIndices[4][TmpNbrNIndices] = new int [Lim * 3];
		  this->MNNBodyInteractionFactors[4][TmpNbrNIndices] = new double [Lim];
		  int* TmpMIndices = this->MIndices[4][TmpNbrNIndices];
		  int* TmpMIndices2 = this->SortedIndicesPerSum[4][MinSum];
		  double* TmpInteraction = this->MNNBodyInteractionFactors[4][TmpNbrNIndices];
		  for (int j = 0; j < Lim; ++j)
		    {
		      for (int l = 0; l < 4; ++l)
			{
			  (*TmpMIndices) = (*TmpMIndices2);			
			  ++TmpMIndices;
			  ++TmpMIndices2;
			}			
		      double& TmpInteraction2 = TmpInteraction[j];
		      TmpInteraction2 = 0.0;
		      if (this->FiveBodyPseudoPotential[4] != 0.0)
			TmpInteraction2 += this->FiveBodyPseudoPotential[4] * TmpProjectorCoefficients[4][i] * TmpProjectorCoefficients[4][j];
		      for (int k = 5; k <= TmpMaxRealtiveMonentum; ++k)  
			if (this->FiveBodyPseudoPotential[k] != 0.0)
			  TmpInteraction2 += this->FiveBodyPseudoPotential[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
		    }
		  for (int j = 0; j < 4; ++j)
		    {
		      (*TmpNIndices) = (*TmpNIndices2);			
		      ++TmpNIndices;
		      ++TmpNIndices2;
		    }
		  ++TmpNbrNIndices;
		}
	      if (this->FiveBodyPseudoPotential[4] != 0.0)
		delete[] TmpProjectorCoefficients[4];
	      for (int i = 5; i <= TmpMaxRealtiveMonentum; ++i)  
		if (this->FiveBodyPseudoPotential[i] != 0.0)
		  delete[] TmpProjectorCoefficients[i];
	      delete[] TmpProjectorCoefficients;		
	    }
	}
      delete[] TmpInteractionCoeffients;
    }
  else
    {
      this->MinSumIndices[4] = 0;
      this->MaxSumIndices[4] = this->LzMax * 4;
      double* TmpInteractionCoeffients = new double[this->MaxSumIndices[4] + 1];
      double Coefficient;
      TmpInteractionCoeffients[0] = 1.0;
      TmpInteractionCoeffients[1] = 1.0;
      for (int i = 2; i <= this->MaxSumIndices[4]; ++i)
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
      Coefficient = 4.0 * M_PI / (((double) this->MaxSumIndices[4]) + 1.0);
      double Radius = 2.0 / ((double) this->LzMax);
      for (int i = 2; i <= 4; ++i)
	{
	  Coefficient *= (double) (i * i);	  
	  Coefficient *= Radius;
	}
      for (int i = 0; i <= this->MaxSumIndices[4]; ++i)
	TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);
      
      double** SortedIndicesPerSumSymmetryFactor;
      GetAllSymmetricIndices(this->NbrLzValue, 4, this->NbrSortedIndicesPerSum[4], this->SortedIndicesPerSum[4],
			     SortedIndicesPerSumSymmetryFactor);
      
      
      long TmpNbrNIndices = 0;
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[4]; ++MinSum)
	TmpNbrNIndices += this->NbrSortedIndicesPerSum[4][MinSum];
      this->NbrNIndices[4] = TmpNbrNIndices;
      this->NIndices[4] = new int[TmpNbrNIndices * 4];
      this->NbrMIndices[4] = new long[TmpNbrNIndices];
      this->MIndices[4] = new int*[TmpNbrNIndices];
      this->MNNBodyInteractionFactors[4] = new double* [TmpNbrNIndices];
      TmpNbrNIndices = 0;	 
      int* TmpNIndices = this->NIndices[4];
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[4]; ++MinSum)
	{
	  int Lim = this->NbrSortedIndicesPerSum[4][MinSum];
	  double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
	  int* TmpNIndices2 = this->SortedIndicesPerSum[4][MinSum];
	  int TmpMaxRealtiveMonentum = 5;
	  if (this->MaxRelativeAngularMomentum <= TmpMaxRealtiveMonentum)
	    TmpMaxRealtiveMonentum = this->MaxRelativeAngularMomentum;
	  int TmpSum = TmpNIndices2[0] + TmpNIndices2[1] + TmpNIndices2[2] + TmpNIndices2[3];
	  while (((4 * this->LzMax) - TmpMaxRealtiveMonentum)  < TmpSum)
	    --TmpMaxRealtiveMonentum;
	  double** TmpProjectorCoefficients = new double* [this->MaxRelativeAngularMomentum + 1];
	  if (this->FiveBodyPseudoPotential[0] != 0.0)
	    TmpProjectorCoefficients[0] = this->ComputeProjectorCoefficients(0, 1, TmpNIndices2, Lim);
	  for (int i = 1; i <= TmpMaxRealtiveMonentum; ++i)  
	    if (this->FiveBodyPseudoPotential[i] != 0.0)
	      TmpProjectorCoefficients[i] = this->ComputeProjectorCoefficients(2 * i, 1, TmpNIndices2, Lim);
	  if (this->MaxRelativeAngularMomentum >= 5)
	      TmpProjectorCoefficients[5] = this->ComputeProjectorCoefficients(8, 2, TmpNIndices2, Lim);
	  if (this->MaxRelativeAngularMomentum >= 6)
	      TmpProjectorCoefficients[6] = this->ComputeProjectorCoefficients(10, 1, TmpNIndices2, Lim);
	  if (this->MaxRelativeAngularMomentum >= 7)
	      TmpProjectorCoefficients[7] = this->ComputeProjectorCoefficients(10, 2, TmpNIndices2, Lim);
	    

	  for (int i = 0; i < Lim; ++i)
	    {
	      this->NbrMIndices[4][TmpNbrNIndices] = Lim;		    
	      this->MIndices[4][TmpNbrNIndices] = new int [Lim * 4];
	      this->MNNBodyInteractionFactors[4][TmpNbrNIndices] = new double [Lim];
	      int* TmpMIndices = this->MIndices[4][TmpNbrNIndices];
	      int* TmpMIndices2 = this->SortedIndicesPerSum[4][MinSum];
	      double* TmpInteraction = this->MNNBodyInteractionFactors[4][TmpNbrNIndices];
	      for (int j = 0; j < Lim; ++j)
		{
		  double Coefficient2 = TmpSymmetryFactors[j];
		  for (int l = 0; l < 4; ++l)
		    {
		      Coefficient2 *= TmpNormalizationCoeffients[(*TmpMIndices2)];		    
		      (*TmpMIndices) = (*TmpMIndices2);			
		      ++TmpMIndices;
		      ++TmpMIndices2;
		    }			
		  double& TmpInteraction2 = TmpInteraction[j];
		  TmpInteraction2 = 0.0;
		  if (this->FiveBodyPseudoPotential[0] != 0.0)
		    TmpInteraction2 += this->FiveBodyPseudoPotential[0] * TmpProjectorCoefficients[0][i] * TmpProjectorCoefficients[0][j];
		  for (int k = 2; k <= TmpMaxRealtiveMonentum; ++k)  
		    if (this->FiveBodyPseudoPotential[k] != 0.0)
		      TmpInteraction2 += this->FiveBodyPseudoPotential[k] * TmpProjectorCoefficients[k][i] * TmpProjectorCoefficients[k][j];
		  TmpInteraction2 *= TmpSymmetryFactors[i] * TmpSymmetryFactors[j];
		}
	      for (int j = 0; j < 4; ++j)
		{
		  (*TmpNIndices) = (*TmpNIndices2);			
		  ++TmpNIndices;
		  ++TmpNIndices2;
		}
	      ++TmpNbrNIndices;
	    }		
	  if (this->FiveBodyPseudoPotential[0] != 0.0)
	    delete[] TmpProjectorCoefficients[0];
	  for (int i = 2; i <= TmpMaxRealtiveMonentum; ++i)  
	    if (this->FiveBodyPseudoPotential[i] != 0.0)
	      delete[] TmpProjectorCoefficients[i];
	  delete[] TmpProjectorCoefficients;		
	}
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[4]; ++MinSum)
	{
	  delete[] SortedIndicesPerSumSymmetryFactor[MinSum];
	}
      delete[] SortedIndicesPerSumSymmetryFactor;
      delete[] TmpInteractionCoeffients;
    }
  delete[] TmpNormalizationCoeffients;
}

// compute all projector coefficient associated to a given relative angular momentum between 4 particles
//
// relativeMomentum = value of twice the relative angular momentum between the 4 particles
// degeneracyIndex = optional degeneracy index for relative angular momentum greater than 5 for bosons (8 for fermions)
// indices = array that contains all possible sets of indices (size of the array is 4 * nbrIndexSets)
// nbrIndexSets = number of sets

double* ParticleOnSphereGenericFiveBodyHamiltonian::ComputeProjectorCoefficients(int relativeMomentum, int degeneracyIndex, int* indices, int nbrIndexSets)
{
  double* TmpCoefficients = new double [nbrIndexSets];
  int JValue = (4 * this->LzMax) - relativeMomentum;

  int MaxJ = 3 * this->LzMax;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    MaxJ -= 2;
  ClebschGordanCoefficients Clebsh (this->LzMax, this->LzMax);
  ClebschGordanCoefficients* ClebshArray = 0;
  ClebschGordanCoefficients** ClebshArray2 = 0;
  if (degeneracyIndex == 2)
    {      
      ClebshArray2 = new ClebschGordanCoefficients*[2 * this->LzMax +1];
      for (int j = (2 * this->LzMax); j >= 0; --j)
	{
	  ClebshArray2[j] = new ClebschGordanCoefficients[2 * this->LzMax +1];
	  for (int k = (2 * this->LzMax); k >= 0; --k)
	    ClebshArray2[j][k] = ClebschGordanCoefficients(j, k);
	}
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	for (int i = 0; i < nbrIndexSets; ++i)
	  {
	    double Tmp = 0.0;
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[0], indices[1], indices[2], indices[3], JValue, ClebshArray2);
	    Tmp -= this->ComputeProjectorCoefficients4BodySecondChannel(indices[0], indices[1], indices[3], indices[2], JValue, ClebshArray2);
	    Tmp -= this->ComputeProjectorCoefficients4BodySecondChannel(indices[0], indices[2], indices[1], indices[3], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[0], indices[2], indices[3], indices[1], JValue, ClebshArray2);
	    Tmp -= this->ComputeProjectorCoefficients4BodySecondChannel(indices[0], indices[3], indices[2], indices[1], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[0], indices[3], indices[1], indices[2], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[1], indices[2], indices[0], indices[3], JValue, ClebshArray2);
	    Tmp -= this->ComputeProjectorCoefficients4BodySecondChannel(indices[1], indices[2], indices[3], indices[0], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[1], indices[3], indices[2], indices[0], JValue, ClebshArray2);
	    Tmp -= this->ComputeProjectorCoefficients4BodySecondChannel(indices[1], indices[3], indices[0], indices[2], JValue, ClebshArray2);	  
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[2], indices[3], indices[0], indices[1], JValue, ClebshArray2);
	    Tmp -= this->ComputeProjectorCoefficients4BodySecondChannel(indices[2], indices[3], indices[1], indices[0], JValue, ClebshArray2);
	    TmpCoefficients[i] = Tmp;
	    indices += 4;
	  }
      else
	for (int i = 0; i < nbrIndexSets; ++i)
	  {
	    double Tmp = 0.0;
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[0], indices[1], indices[2], indices[3], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[0], indices[1], indices[3], indices[2], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[0], indices[2], indices[1], indices[3], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[0], indices[2], indices[3], indices[1], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[0], indices[3], indices[2], indices[1], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[0], indices[3], indices[1], indices[2], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[1], indices[2], indices[0], indices[3], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[1], indices[2], indices[3], indices[0], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[1], indices[3], indices[2], indices[0], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[1], indices[3], indices[0], indices[2], JValue, ClebshArray2);	  
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[2], indices[3], indices[0], indices[1], JValue, ClebshArray2);
	    Tmp += this->ComputeProjectorCoefficients4BodySecondChannel(indices[2], indices[3], indices[1], indices[0], JValue, ClebshArray2);
	    TmpCoefficients[i] = Tmp;
	    indices += 4;
	  }
      for (int j = (2 * this->LzMax); j >= 0; --j)
	delete[] ClebshArray2[j];
      delete[] ClebshArray2;
    }
  else
    {
      ClebshArray = new ClebschGordanCoefficients[MaxJ + 1];
      for (int j = MaxJ; j >= 0; --j)
	ClebshArray[j] = ClebschGordanCoefficients(j, this->LzMax);
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	for (int i = 0; i < nbrIndexSets; ++i)
	  {
	    double Tmp = 0.0;
	    Tmp += this->ComputeProjectorCoefficients4Body(indices[0], indices[1], indices[2], indices[3], JValue, ClebshArray);
	    Tmp -= this->ComputeProjectorCoefficients4Body(indices[0], indices[1], indices[3], indices[2], JValue, ClebshArray);
	    Tmp -= this->ComputeProjectorCoefficients4Body(indices[0], indices[2], indices[1], indices[3], JValue, ClebshArray);
	    Tmp += this->ComputeProjectorCoefficients4Body(indices[0], indices[2], indices[3], indices[1], JValue, ClebshArray);
	    Tmp -= this->ComputeProjectorCoefficients4Body(indices[0], indices[3], indices[2], indices[1], JValue, ClebshArray);
	    Tmp += this->ComputeProjectorCoefficients4Body(indices[0], indices[3], indices[1], indices[2], JValue, ClebshArray);
	    Tmp += this->ComputeProjectorCoefficients4Body(indices[1], indices[2], indices[0], indices[3], JValue, ClebshArray);
	    Tmp -= this->ComputeProjectorCoefficients4Body(indices[1], indices[2], indices[3], indices[0], JValue, ClebshArray);
	    Tmp += this->ComputeProjectorCoefficients4Body(indices[1], indices[3], indices[2], indices[0], JValue, ClebshArray);
	    Tmp -= this->ComputeProjectorCoefficients4Body(indices[1], indices[3], indices[0], indices[2], JValue, ClebshArray);	  
	    Tmp += this->ComputeProjectorCoefficients4Body(indices[2], indices[3], indices[0], indices[1], JValue, ClebshArray);
	    Tmp -= this->ComputeProjectorCoefficients4Body(indices[2], indices[3], indices[1], indices[0], JValue, ClebshArray);
	    TmpCoefficients[i] = Tmp;
	    indices += 4;
	  }
      else
	for (int i = 0; i < nbrIndexSets; ++i)
	  {
	    double Tmp = 0.0;
	    Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[2], indices[3], indices[4], JValue, ClebshArray);
	    Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[2], indices[4], indices[3], JValue, ClebshArray);
	    Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[3], indices[2], indices[4], JValue, ClebshArray);
	    Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[3], indices[4], indices[2], JValue, ClebshArray);
	    Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[4], indices[3], indices[2], JValue, ClebshArray);
	    Tmp += this->ComputeProjectorCoefficients5Body(indices[0], indices[1], indices[4], indices[2], indices[3], JValue, ClebshArray);
	    TmpCoefficients[i] = Tmp;
	    indices += 4;
	  }
      delete[] ClebshArray;
    }
  return TmpCoefficients;
}

// compute a given projector coefficient for the 5-body interaction 
//
// m1 = first index
// m2 = second index
// m3 = third inde
// m4 = fourth index
// jValue = total angular momentum
// minJ = minimum angular momentum that can be reach by three particles
// return value = corresponding projector coefficient

double ParticleOnSphereGenericFiveBodyHamiltonian::ComputeProjectorCoefficients4Body(int m1, int m2, int m3, int m4, int jValue, 
										     ClebschGordanCoefficients* clebshArray)
{
  double Tmp = 0.0;
  double Tmp2 = 0.0;
  int TmpMinJ2;
  int Sum1 = ((m1 + m2) << 1)  - (2 * this->LzMax);
  int Sum2 = Sum1 + (m3 << 1)  - this->LzMax;
  int TmpMinJ = abs(Sum1);
  int MaxJ = 2 * this->LzMax;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    MaxJ -= 2;
  for (int j1 = MaxJ; j1 >= TmpMinJ; j1 -= 4)
    {
      Tmp2 = clebshArray[this->LzMax].GetCoefficient(((m1 << 1) - this->LzMax), ((m2 << 1)- this->LzMax), j1);
      TmpMinJ2 = abs(jValue - this->LzMax);
      if (TmpMinJ2 < abs(Sum2))
	TmpMinJ2 = abs(Sum2);
      for (int j2 = j1 + this->LzMax; j2 >= TmpMinJ2; j2 -= 2)
	{
	  if ((abs(Sum2 + ((m4 << 1) - this->LzMax)) <= jValue) && (abs(Sum1 + ((m3 << 1) - this->LzMax)) <= j2))
	    Tmp += Tmp2 * (clebshArray[j1].GetCoefficient(Sum1, ((m3 << 1) - this->LzMax), j2) * 
			   clebshArray[j2].GetCoefficient(Sum2, ((m4 << 1) - this->LzMax), jValue)); 
	}
    }
  return Tmp;
}

// compute a given projector coefficient for the 5-body interaction in the second channel
//
// m1 = first index
// m2 = second index
// m3 = third inde
// m4 = fourth index
// jValue = total angular momentum
// minJ = minimum angular momentum that can be reach by five particles
// return value = corresponding projector coefficient

double ParticleOnSphereGenericFiveBodyHamiltonian::ComputeProjectorCoefficients4BodySecondChannel(int m1, int m2, int m3, int m4, int jValue, 
												  ClebschGordanCoefficients** clebshArray)
{
  double Tmp = 0.0;
  double Tmp2 = 0.0;
  int Sum1 = ((m1 + m2) << 1)  - (2 * this->LzMax);
  int Sum2 = ((m3 + m4) << 1)  - (2 * this->LzMax);
  int MaxJ = 2 * this->LzMax;
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    MaxJ -= 2;
  int TmpMinJ1 = abs(Sum1);
  int TmpMinJ2 = abs(Sum2);
  for (int j1 = MaxJ; j1 >= TmpMinJ1; j1 -= 4)
    {
      Tmp2 = clebshArray[this->LzMax][this->LzMax].GetCoefficient(((m1 << 1) - this->LzMax), ((m2 << 1)- this->LzMax), j1);
      for (int j2 = MaxJ; j2 >= TmpMinJ2; j2 -= 4)
	{
 	  if ((abs(j1 - j2) <= jValue) && ((j1 + j2) >= jValue))
	    {
	      ClebschGordanCoefficients TmpClebsh (j1, j2);
	      Tmp += Tmp2 * (clebshArray[this->LzMax][this->LzMax].GetCoefficient(((m3 << 1) - this->LzMax), ((m4 << 1)- this->LzMax), j2) * 
			     clebshArray[j1][j2].GetCoefficient(Sum1, Sum2, jValue)); 
	    }
	}
    }


  return Tmp;
}
