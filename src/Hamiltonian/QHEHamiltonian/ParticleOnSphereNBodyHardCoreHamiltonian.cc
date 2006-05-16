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
#include "MathTools/ClebschGordanCoefficients.h"

#include <stdio.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


// default constructor
//

ParticleOnSphereNBodyHardCoreHamiltonian::ParticleOnSphereNBodyHardCoreHamiltonian()
{
}

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
  this->MinSumIndices = new int [this->MaxNBody + 1];
  this->MaxSumIndices = new int [this->MaxNBody + 1];
  this->NBodySign = new double[this->MaxNBody + 1];

  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->MinSumIndices[k] = 1;
      this->MaxSumIndices[k] = 0;      
      this->NBodyFlags[k] = false;
      this->NBodyInteractionWeightFactors[k] = 0.0;
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic) && ((k & 1) == 0))
	this->NBodySign[k] = -1.0;
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
  this->MinSumIndices = new int [this->MaxNBody + 1];
  this->MaxSumIndices = new int [this->MaxNBody + 1];
  this->NBodySign = new double[this->MaxNBody + 1];

  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->MinSumIndices[k] = 1;
      this->MaxSumIndices[k] = 0;      
      if (nBodyFactors[k] == 0.0)
	this->NBodyFlags[k] = false;
      else
	this->NBodyFlags[k] = true;
      this->NBodyInteractionWeightFactors[k] = nBodyFactors[k];
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic) && ((k & 1) == 0))
	this->NBodySign[k] = -1.0;
      if (this->NBodyInteractionWeightFactors[k] < 0.0)
	{
	  this->NBodyInteractionWeightFactors[k] *= -1.0;
	  this->NBodySign[k] *= -1.0;
	}
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
      else
	{
	  this->FastMultiplicationFlag = false;
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnSphereNBodyHardCoreHamiltonian::~ParticleOnSphereNBodyHardCoreHamiltonian()
{
  for (int k = 1; k <= this->MaxNBody; ++k)
    if (this->NBodyFlags[k] == true)
      {
	for (int MinSum = this->MinSumIndices[k]; MinSum <= this->MaxSumIndices[k]; ++MinSum)
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
  delete[] this->MinSumIndices;
  delete[] this->MaxSumIndices;
  delete[] this->NBodySign;
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
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true) 
	  {
	    double Coefficient;
	    GetAllSkewSymmetricIndices(this->NbrLzValue, k, this->NbrSortedIndicesPerSum[k], this->SortedIndicesPerSum[k]);
	    this->MaxSumIndices[k] = (((this->NbrLzValue - 1) * this->NbrLzValue) - ((k - 1) * (k - 2)))/ 2;
	    this->MinSumIndices[k] = (k * (k - 1)) / 2;
	    double* TmpInteractionCoeffients = new double[this->MaxSumIndices[k] + 1];
	    TmpInteractionCoeffients[0] = 1.0;
	    TmpInteractionCoeffients[1] = 1.0;
	    for (int i = 2; i <= this->MaxSumIndices[k]; ++i)
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
	    Coefficient = 4.0 * M_PI / (((double) this->MaxSumIndices[k]) + 1.0);
	    double Radius = 2.0 / ((double) this->LzMax);
	    for (int i = 2; i <= k; ++i)
	      {
		Coefficient *= (double) (i * i);	  
		Coefficient *= Radius;
	      }
	    for (int i = 0; i <= this->MaxSumIndices[k]; ++i)
	      TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);

	    this->NBodyInteractionFactors[k] = new double* [this->MaxSumIndices[k] + 1];
	    int Lim;
	    for (int MinSum = this->MinSumIndices[k]; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	      {
		Lim = this->NbrSortedIndicesPerSum[k][MinSum];
		this->NBodyInteractionFactors[k][MinSum] = new double [Lim];
		double* TmpNBodyInteractionFactors = this->NBodyInteractionFactors[k][MinSum];		  
		int* TmpMIndices = this->SortedIndicesPerSum[k][MinSum];
		double* TmpProjectorCoefficients = this->ComputeProjectorCoefficients(k, TmpMIndices, Lim);
		for (int i = 0; i < Lim; ++i)
		  {
		    Coefficient = TmpInteractionCoeffients[MinSum] * TmpProjectorCoefficients[i];
		    for (int l = 0; l < k; ++l)
		      {
			Coefficient *= TmpNormalizationCoeffients[TmpMIndices[l]];		    
		      }
		    TmpNBodyInteractionFactors[i] = TmpProjectorCoefficients[i] * this->NBodyInteractionWeightFactors[k];
		    TmpMIndices += k;
		  }
		delete[] TmpProjectorCoefficients;
	      }
	    delete[] TmpInteractionCoeffients;
	  }
    }
  else
    {
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true) 
	  {
	    this->MinSumIndices[k] = 0;
	    this->MaxSumIndices[k] = this->LzMax * k;
	    double* TmpInteractionCoeffients = new double[this->MaxSumIndices[k] + 1];
	    double Coefficient;
	    TmpInteractionCoeffients[0] = 1.0;
	    TmpInteractionCoeffients[1] = 1.0;
	    for (int i = 2; i <= this->MaxSumIndices[k]; ++i)
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
	    Coefficient = 4.0 * M_PI / (((double) this->MaxSumIndices[k]) + 1.0);
	    double Radius = 2.0 / ((double) this->LzMax);
	    for (int i = 2; i <= k; ++i)
	      {
		Coefficient *= (double) (i * i);	  
		Coefficient *= Radius;
	      }
	    for (int i = 0; i <= this->MaxSumIndices[k]; ++i)
	      TmpInteractionCoeffients[i] = sqrt(Coefficient / TmpInteractionCoeffients[i]);

	    double** SortedIndicesPerSumSymmetryFactor;
	    GetAllSymmetricIndices(this->NbrLzValue, k, this->NbrSortedIndicesPerSum[k], this->SortedIndicesPerSum[k],
				   SortedIndicesPerSumSymmetryFactor);
	    this->NBodyInteractionFactors[k] = new double* [this->MaxSumIndices[k] + 1];
 	    int Lim;
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices[k]; ++MinSum)
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
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	      {
		delete[] SortedIndicesPerSumSymmetryFactor[MinSum];
	      }
	    delete[] SortedIndicesPerSumSymmetryFactor;
	    delete[] TmpInteractionCoeffients;
	  }
    }
  delete[] TmpNormalizationCoeffients;
}

// compute all projector coefficient associated to a given a 
//
// nbrIndices = number of indices per set
// indices = array that contains all possible sets of indices (size of the array is nbrIndices * nbrIndexSets)
// nbrIndexSets = number of sets

double* ParticleOnSphereNBodyHardCoreHamiltonian::ComputeProjectorCoefficients(int nbrIndices, int* indices, int nbrIndexSets)
{
  double* TmpCoefficients = new double [nbrIndexSets];
  int JValue = (nbrIndices * (this->LzMax - (nbrIndices - 1)));
  switch (nbrIndices)
    {
    case 2:
      {
	ClebschGordanCoefficients Clebsh (this->LzMax, this->LzMax);
	for (int i = 0; i < nbrIndexSets; ++i)
	  {
	    TmpCoefficients[i] = Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[1] << 1) - this->LzMax), JValue);
	    indices += 2;
	  }
      }
      break;
    case 3:
      {
	int MaxJ = 2 * this->LzMax - 2;
	ClebschGordanCoefficients Clebsh (this->LzMax, this->LzMax);
	ClebschGordanCoefficients* ClebshArray = new ClebschGordanCoefficients[MaxJ + 1];
	int MinJ = JValue - this->LzMax;
	if (MinJ < 0)
	  MinJ = 0;
	for (int j = MaxJ; j >= MinJ; j -= 4)
	  ClebshArray[j] = ClebschGordanCoefficients(j, this->LzMax);
	for (int i = 0; i < nbrIndexSets; ++i)
	  {
	    double Tmp = 0.0;
	    int Sum = ((indices[0] + indices[1]) << 1)  - (2 * this->LzMax);
	    int TmpMinJ = MinJ;
	    if (TmpMinJ < abs(Sum))
	      TmpMinJ = abs(Sum);
	    for (int j = MaxJ; j >= TmpMinJ; j -= 4)
	      {
		Tmp += (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[1] << 1)- this->LzMax), j) * 
			ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue)); 
	      }
	    Sum = ((indices[1] + indices[2]) << 1)  - (2 * this->LzMax);
	    TmpMinJ = MinJ;
	    if (TmpMinJ < abs(Sum))
	      TmpMinJ = abs(Sum);
	    for (int j = MaxJ; j >= TmpMinJ; j -= 4)
	      {
		Tmp += (Clebsh.GetCoefficient(((indices[1] << 1) - this->LzMax), ((indices[2] << 1)- this->LzMax), j) * 
			ClebshArray[j].GetCoefficient(Sum, ((indices[0] << 1) - this->LzMax), JValue)); 
	      }
	    Sum = ((indices[2] + indices[0]) << 1)  - (2 * this->LzMax);
	    TmpMinJ = MinJ;
	    if (TmpMinJ < abs(Sum))
	      TmpMinJ = abs(Sum);
	    for (int j = MaxJ; j >= TmpMinJ; j -= 4)
	      {
		Tmp += (Clebsh.GetCoefficient(((indices[2] << 1) - this->LzMax), ((indices[0] << 1)- this->LzMax), j) * 
			ClebshArray[j].GetCoefficient(Sum, ((indices[1] << 1) - this->LzMax), JValue)); 
	      }
	    TmpCoefficients[i] = Tmp;
	    indices += 3;
	  }
	delete[] ClebshArray;
      }
      break;
    default:
      {
      }
    }
  return TmpCoefficients;
}
