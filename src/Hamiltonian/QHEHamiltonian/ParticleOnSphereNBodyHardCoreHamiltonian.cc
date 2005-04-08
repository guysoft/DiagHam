////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//                         n-body hard core interaction                       //
//                                                                            //
//                        last modification : 23/09/2004                      //
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
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereNBodyHardCoreHamiltonian.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>


using std::cout;
using std::endl;


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// architecture = architecture to use for precalculation
// nbrBody = number of particle that interact simultaneously through the hard core interaction
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereNBodyHardCoreHamiltonian::ParticleOnSphereNBodyHardCoreHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, int nbrBody,
										   AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag,
										   char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->NbrNbody = nbrBody;
  this->MaxNBody = this->NbrNbody;
  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NBodyInteractionFactors = new double** [this->MaxNBody + 1];
  this->NBodyInteractionWeightFactors = new double [this->MaxNBody + 1];
  this->NbrSortedIndicesPerSum = new int* [this->MaxNBody + 1];
  this->SortedIndicesPerSum = new int** [this->MaxNBody + 1];

  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->NBodyFlags[k] = false;
      this->NBodyInteractionWeightFactors[k] = 0.0;
    }
  this->NBodyFlags[this->NbrNbody] = true;
  this->NBodyInteractionWeightFactors[this->NbrNbody] = 1.0;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
  if (this->DiskStorageFlag == false)
    cout << "toto" << endl;
  else
    cout << "tata" << endl;
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
	      cout << "toto" << endl;
	      this->EnableFastMultiplication();
	    }
	  else
	    {
	      cout << "tata" << endl;
	      char* TmpFileName = this->Architecture->GetTemporaryFileName();
	      this->EnableFastMultiplicationWithDiskStorage(TmpFileName);	      
	      delete[] TmpFileName;
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// architecture = architecture to use for precalculation
// maxNbrBody = maximum number of particle that interact simultaneously through the hard core interaction
// nBodyFactors = weight of the different n-body interaction terms with respect to each other
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereNBodyHardCoreHamiltonian::ParticleOnSphereNBodyHardCoreHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
										   int maxNbrBody, double* nBodyFactors,
										   AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
										   char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->NbrNbody = maxNbrBody;
  this->MaxNBody = maxNbrBody;
  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NBodyInteractionFactors = new double** [this->MaxNBody + 1];
  this->NBodyInteractionWeightFactors = new double [this->MaxNBody + 1];
  this->NbrSortedIndicesPerSum = new int* [this->MaxNBody + 1];
  this->SortedIndicesPerSum = new int** [this->MaxNBody + 1];

  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      if (nBodyFactors[k] == 0.0)
	this->NBodyFlags[k] = false;
      else
	this->NBodyFlags[k] = true;
      this->NBodyInteractionWeightFactors[k] = nBodyFactors[k];
    }
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
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnSphereNBodyHardCoreHamiltonian::~ParticleOnSphereNBodyHardCoreHamiltonian()
{
  for (int k = 2; k <= this->MaxNBody; ++k)
    if (this->NBodyFlags[k] == true)
      {
	for (int MinSum = this->MinSumIndices; MinSum <= this->MaxSumIndices; ++MinSum)
	  {
	    delete[] this->SortedIndicesPerSum[k][MinSum];
	    delete[] this->NBodyInteractionFactors[k][MinSum];
	  }
	delete[] this->NbrSortedIndicesPerSum[k];
	delete[] this->SortedIndicesPerSum[k];
	delete[] this->NBodyInteractionFactors[k];
      }
  delete[] this->NBodyFlags;
  delete[] this->NBodyInteractionFactors;
  delete[] this->SortedIndicesPerSum;
  delete[] this->NbrSortedIndicesPerSum;
  delete[] this->NBodyInteractionWeightFactors;

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
      delete[] this->NbrInteractionPerComponent;
    }
}

// evaluate all interaction factors
//   

void ParticleOnSphereNBodyHardCoreHamiltonian::EvaluateInteractionFactors()
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
      // useless part (trivially equal to zero for fermions), to be used as an example for other interactions
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true) 
	  {
	    double SumCoefficient = 0.0;
	    double Coefficient;
	    GetAllSkewSymmetricIndices(this->NbrLzValue, k, this->NbrSortedIndicesPerSum[k], this->SortedIndicesPerSum[k]);
	    this->MaxSumIndices = (((this->NbrLzValue - 1) * this->NbrLzValue) - ((k - 1) * (k - 2)))/ 2;
	    this->MinSumIndices = (k * (k - 1)) / 2;
	    this->NBodyInteractionFactors[k] = new double* [MaxSumIndices + 1];
	    int Lim;
	    for (int MinSum = this->MinSumIndices; MinSum <= MaxSumIndices; ++MinSum)
	      {
		Lim = this->NbrSortedIndicesPerSum[k][MinSum];
		this->NBodyInteractionFactors[k][MinSum] = new double [Lim];
		double* TmpNBodyInteractionFactors = this->NBodyInteractionFactors[k][MinSum];		  
		int* TmpMIndices = this->SortedIndicesPerSum[k][MinSum];
		for (int i = 0; i < Lim; ++i)
		  {
		    Coefficient = SumCoefficient;
		    for (int l = 0; l < k; ++l)
		      Coefficient *= TmpNormalizationCoeffients[TmpMIndices[l]];		    
		    TmpNBodyInteractionFactors[i] = Coefficient * this->NBodyInteractionWeightFactors[k];
		    TmpMIndices += k;
		  }
	      }
	  }
    }
  else
    {
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true) 
	  {
	    this->MinSumIndices = 0;
	    this->MaxSumIndices = this->LzMax * k;
	    double* TmpInteractionCoeffients = new double[MaxSumIndices + 1];
	    double Coefficient;
	    TmpInteractionCoeffients[0] = 1.0;
	    TmpInteractionCoeffients[1] = 1.0;
	    for (int i = 2; i <= MaxSumIndices; ++i)
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
	    Coefficient = 4.0 * M_PI / (((double) MaxSumIndices) + 1.0);
	    double Radius = 2.0 / ((double) this->LzMax);
	    for (int i = 2; i <= k; ++i)
	      {
		Coefficient *= (double) (i * i);	  
		Coefficient *= Radius;
	      }
	    for (int i = 0; i <= MaxSumIndices; ++i)
	      TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);

	    double** SortedIndicesPerSumSymmetryFactor;
	    GetAllSymmetricIndices(this->NbrLzValue, k, this->NbrSortedIndicesPerSum[k], this->SortedIndicesPerSum[k],
				   SortedIndicesPerSumSymmetryFactor);
	    this->NBodyInteractionFactors[k] = new double* [MaxSumIndices + 1];
	    int Lim;
	    for (int MinSum = 0; MinSum <= MaxSumIndices; ++MinSum)
	      {
		Lim = this->NbrSortedIndicesPerSum[k][MinSum];
		double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
		this->NBodyInteractionFactors[k][MinSum] = new double [Lim];
		double* TmpNBodyInteractionFactors = this->NBodyInteractionFactors[k][MinSum];		  
		int* TmpMIndices = this->SortedIndicesPerSum[k][MinSum];
		for (int i = 0; i < Lim; ++i)
		  {
		    Coefficient = TmpSymmetryFactors[i] * TmpInteractionCoeffients[MinSum];
		    for (int l = 0; l < k; ++l)
		      Coefficient *= TmpNormalizationCoeffients[TmpMIndices[l]];		    
		    TmpNBodyInteractionFactors[i] = Coefficient * this->NBodyInteractionWeightFactors[k];
		    TmpMIndices += k;
		  }
	      }
	    for (int MinSum = 0; MinSum <= MaxSumIndices; ++MinSum)
	      {
		delete[] SortedIndicesPerSumSymmetryFactor[MinSum];
	      }
	    delete[] SortedIndicesPerSumSymmetryFactor;
	    delete[] TmpInteractionCoeffients;
	  }
    }
  delete[] TmpNormalizationCoeffients;
}

