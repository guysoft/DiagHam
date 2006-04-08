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
//   int Indices[2];
//   for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
//     {
//       this->Particles->PrintState(cout, i) << endl;
//       for (int m1 = 0; m1 <= this->LzMax; ++m1)
// 	for (int m2 = m1 + 1; m2 <= this->LzMax; ++m2)
// 	  {
// 	    Indices[0] = m1;
// 	    Indices[1] = m2;
// 	    double TmpCoef = this->Particles->ProdA(i, Indices, 2);
// 	    for (int m3 = 0; m3 <= this->LzMax; ++m3)
// 	      {
// 		int m4 = m1 + m2 - m3;
// 		if ((m4 > m3) && (m4 <= this->LzMax))
// 		  {
// 		    if (TmpCoef != 0.0)
// 		      {
// 			Indices[0] = m3;
// 			Indices[1] = m4;		  
// 			double Coef1; 
// 			int Index1 = this->Particles->ProdAd(Indices, 2, Coef1);
// 			Coef1 *= TmpCoef;
// 			double Coef2; 
// 			int Index2 = this->Particles->AdAdAA(i, m3, m4, m1, m2, Coef2);
// 			if ((Index1 != Index2) || (Coef1 != Coef2))
// 			  {
// 			    cout << "error at index " << i << " : " << m1 << " " << m2 << " " << m3 << " " << m4 << " gives " << Index1 << "(" << Index2 
// 				 << ") and " << Coef1 << "(" << Coef2<< ")" << endl;
// 			  }
// 			else
// 			  {
// 			    cout << "OK at index " << i << " : " << m1 << " " << m2 << " " << m3 << " " << m4 << endl;
// 			  }
// 		      }
// 		    else
// 		      {
// 			double Coef2; 
// 			int Index2 = this->Particles->AdAdAA(i, m3, m4, m1, m2, Coef2);
// 			if (Index2 != this->Particles->GetHilbertSpaceDimension())
// 			  {
// 			    cout << "(null) error at index " << i << " : " << m1 << " " << m2 << " " << m3 << " " << m4 << " gives " 
// 				 << this->Particles->GetHilbertSpaceDimension() <<  "(" << Index2 
// 				 << ") and 0 (" << Coef2<< ")" << endl;
// 			  }
// 			else
// 			  {
// 			    cout << "(null) OK at index " << i << " : " << m1 << " " << m2 << " " << m3 << " " << m4 << endl;
// 			  }
// 		      }
// 		  }
// 	      }
// 	  }
// //     }
//    cout << "MinSumIndices " << this->MinSumIndices << "  MaxSumIndices " << this->MaxSumIndices << endl;
//    for (int i = this->MinSumIndices; i <= this->MaxSumIndices; ++i)
//     {
//       cout << "sum = " << i << "    nbr indices = " << this->NbrSortedIndicesPerSum[3][i] << endl;
//       int* Indices2 = this->SortedIndicesPerSum[3][i];
//       for (int j = 0; j < this->NbrSortedIndicesPerSum[3][i]; ++j)
// 	{
// 	  cout << Indices2[0] << " " << Indices2[1] << " " << Indices2[2] << endl; 
// 	  Indices2 += 3;
// 	}
//     }
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
      // useless part (trivially equal to zero for fermions), to be used as an example for other interactions
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true) 
	  {
	    double Coefficient;
	    GetAllSkewSymmetricIndices(this->NbrLzValue, k, this->NbrSortedIndicesPerSum[k], this->SortedIndicesPerSum[k]);
	    this->MaxSumIndices = (((this->NbrLzValue - 1) * this->NbrLzValue) - ((k - 1) * (k - 2)))/ 2;
	    this->MinSumIndices = (k * (k - 1)) / 2;
	    double* TmpInteractionCoeffients = new double[MaxSumIndices + 1];
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

	    this->NBodyInteractionFactors[k] = new double* [MaxSumIndices + 1];
	    int Lim;
	    for (int MinSum = this->MinSumIndices; MinSum <= MaxSumIndices; ++MinSum)
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
		    //		    TmpNBodyInteractionFactors[i] = Coefficient * this->NBodyInteractionWeightFactors[k];
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

// compute all projector coefficient associated to a given a 
//
// nbrIndices = number of indices per set
// indices = array that contains all possible sets of indices (size of the array is nbrIndices * nbrIndexSets)
// nbrIndexSets = number of sets

double* ParticleOnSphereNBodyHardCoreHamiltonian::ComputeProjectorCoefficients(int nbrIndices, int* indices, int nbrIndexSets)
{
  double* TmpCoefficients = new double [nbrIndexSets];
  //  int JValue = (nbrIndices * (this->LzMax - (nbrIndices - 1)));
  int JValue = (nbrIndices * this->LzMax);
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
	int MaxJ = 2 * this->LzMax;
	ClebschGordanCoefficients Clebsh (this->LzMax, this->LzMax);
	ClebschGordanCoefficients* ClebshArray = new ClebschGordanCoefficients[MaxJ + 1];
	int MinJ = JValue - this->LzMax;
	if (MinJ < 0)
	  MinJ = 0;
	for (int j = MinJ; j <= MaxJ; j += 2)
	  ClebshArray[j] = ClebschGordanCoefficients(j, this->LzMax);
	for (int i = 0; i < nbrIndexSets; ++i)
	  {
	    double Tmp = 0.0;
	    int Sum = ((indices[0] + indices[1]) << 1)  - (2 * this->LzMax);
	    int j = MinJ;
	    if (j < abs(Sum))
	      j = abs(Sum);
	    for (; j <= MaxJ; j += 2)
	      {
// 		cout << ((indices[0] << 1) - this->LzMax) << " " << ((indices[1] << 1)- this->LzMax) << " " << j 
// 		     << " " << Sum << " " << ((indices[2] << 1) - this->LzMax) << " " << JValue << endl;
		Tmp += (Clebsh.GetCoefficient(((indices[0] << 1) - this->LzMax), ((indices[1] << 1)- this->LzMax), j) * 
			ClebshArray[j].GetCoefficient(Sum, ((indices[2] << 1) - this->LzMax), JValue)); 
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
