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
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereNBodyHardCoreHamiltonian::ParticleOnSphereNBodyHardCoreHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, int nbrBody,
										   AbstractArchitecture* architecture, long memory, 
										   char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->NbrNbody = nbrBody;
  this->MaxNBody = this->NbrNbody;
  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NBodyNbrInteractionFactors = new long [this->MaxNBody + 1];
  this->NBodyInteractionFactors = new double* [this->MaxNBody + 1];
  this->NBodyMValue = new int** [this->MaxNBody + 1];
  this->NBodyNValue = new int** [this->MaxNBody + 1];
  for (int k = 0; k <= this->MaxNBody; ++k)
    this->NBodyFlags[k] = false;
  this->NBodyFlags[this->NbrNbody] = true;

  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
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
	  if (memory > 0)
	    {
	      this->EnableFastMultiplication();
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
	for (int j = 0; j < this->NBodyNbrInteractionFactors[k]; ++j) 
	  {
	    delete[] this->NBodyMValue[k][j];
	    delete[] this->NBodyNValue[k][j];
	  }
	delete[] this->NBodyMValue[k];
	delete[] this->NBodyNValue[k];
	delete[] this->NBodyInteractionFactors[k];
      }
  delete[] this->NBodyFlags;
  delete[] this->NBodyNbrInteractionFactors;
  delete[] this->NBodyInteractionFactors;
  delete[] this->NBodyMValue;
  delete[] this->NBodyNValue;

  if (this->FastMultiplicationFlag == true)
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
	    int* NbrSortedIndicesPerSum;
	    int*** SortedIndicesPerSum;
	    GetAllSkewSymmetricIndices(this->NbrLzValue, k, NbrSortedIndicesPerSum, SortedIndicesPerSum);
	    long TmpNbrInteraction = 0;
	    int MaxSum = (((this->NbrLzValue - 1) * this->NbrLzValue) - ((k - 1) * (k - 2)))/ 2;
	    int MinSum = (k * (k - 1)) / 2;
	    for (int i = MinSum; i <= MaxSum; ++i)
	      TmpNbrInteraction += NbrSortedIndicesPerSum[i];
	    this->NBodyMValue[k] = new int* [TmpNbrInteraction * TmpNbrInteraction];
	    this->NBodyNValue[k] = new int* [TmpNbrInteraction * TmpNbrInteraction];
	    this->NBodyInteractionFactors[k] = new double [TmpNbrInteraction * TmpNbrInteraction];
	    double SumCoefficient = 0.0;
	    double Coefficient;
	    int Lim;
	    for (; MinSum <= MaxSum; ++MinSum)
	      {
		Lim = NbrSortedIndicesPerSum[MinSum];
		for (int i = 0; i < Lim; ++i)
		  {
		    int* TmpMIndices = SortedIndicesPerSum[MinSum][i];
		    Coefficient = SumCoefficient;
		    for (int l = 0; l < k; ++l)
		      Coefficient *= TmpNormalizationCoeffients[TmpMIndices[l]];
		    for (int j = 0; j < Lim; ++j)
		      {
			int* TmpNIndices = SortedIndicesPerSum[MinSum][j];
			this->NBodyMValue[k][TmpNbrInteraction] = new int[k];
			this->NBodyNValue[k][TmpNbrInteraction] = new int[k];
			double TmpCoef = Coefficient;
			for (int l = 0; l < k; ++l)
			  {
			    TmpCoef *= TmpNormalizationCoeffients[TmpNIndices[l]];
			    this->NBodyMValue[k][TmpNbrInteraction][l] = TmpMIndices[l];
			    this->NBodyNValue[k][TmpNbrInteraction][l] = TmpNIndices[l];
			  }
			this->NBodyInteractionFactors[k][TmpNbrInteraction] = TmpCoef;
			++TmpNbrInteraction;
		      }
		  }
	      }
	    MinSum = (k * (k - 1)) / 2;
	    for (; MinSum <= MaxSum; ++MinSum)
	      {
		for (int i = 0; i < NbrSortedIndicesPerSum[MinSum]; ++i)
		  delete[] SortedIndicesPerSum[MinSum][i];
		delete[] SortedIndicesPerSum[MinSum];
	      }
	    delete[] NbrSortedIndicesPerSum;
	    delete[] SortedIndicesPerSum;
	    this->NBodyNbrInteractionFactors[k] = TmpNbrInteraction;
	  }
    }
  else
    {
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true) 
	  {
	    int MaxSum = this->LzMax * k;
	    double* TmpInteractionCoeffients = new double[MaxSum + 1];
	    double Coefficient;
	    TmpInteractionCoeffients[0] = 1.0;
	    TmpInteractionCoeffients[1] = 1.0;
	    for (int i = 2; i <= MaxSum; ++i)
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
	    Coefficient = 4.0 * M_PI / (((double) MaxSum) + 1.0);
	    double Radius = 2.0 / ((double) this->LzMax);
	    for (int i = 2; i <= k; ++i)
	      {
		Coefficient *= (double) (i * i);	  
		Coefficient *= Radius;
	      }
	    for (int i = 0; i <= MaxSum; ++i)
	      TmpInteractionCoeffients[i] = Coefficient / TmpInteractionCoeffients[i];

	    int* NbrSortedIndicesPerSum;
	    int*** SortedIndicesPerSum;
	    double** SortedIndicesPerSumSymmetryFactor;
	    GetAllSymmetricIndices(this->NbrLzValue, k, NbrSortedIndicesPerSum, SortedIndicesPerSum,
				   SortedIndicesPerSumSymmetryFactor);
	    long TmpNbrInteraction = 0;
	    for (int i = 0; i <= MaxSum; ++i)
	      TmpNbrInteraction += (long) (NbrSortedIndicesPerSum[i] * NbrSortedIndicesPerSum[i]);
	    this->NBodyMValue[k] = new int* [TmpNbrInteraction];
	    this->NBodyNValue[k] = new int* [TmpNbrInteraction];
	    this->NBodyInteractionFactors[k] = new double [TmpNbrInteraction];
	    TmpNbrInteraction = 0;
	    double SumCoefficient = 1.0;
	    int Lim;
	    for (int MinSum = 0; MinSum <= MaxSum; ++MinSum)
	      {
		Lim = NbrSortedIndicesPerSum[MinSum];
		for (int i = 0; i < Lim; ++i)
		  {
		    int* TmpMIndices = SortedIndicesPerSum[MinSum][i];
		    Coefficient = SumCoefficient * SortedIndicesPerSumSymmetryFactor[MinSum][i] * TmpInteractionCoeffients[MinSum];
		    for (int l = 0; l < k; ++l)
		      Coefficient *= TmpNormalizationCoeffients[TmpMIndices[l]];		    
		    for (int j = 0; j < Lim; ++j)
		      {
			int* TmpNIndices = SortedIndicesPerSum[MinSum][j];
			this->NBodyMValue[k][TmpNbrInteraction] = new int[k];
			this->NBodyNValue[k][TmpNbrInteraction] = new int[k];
			double TmpCoef = Coefficient;
			for (int l = 0; l < k; ++l)
			  {
			    TmpCoef *= TmpNormalizationCoeffients[TmpNIndices[l]];
			    this->NBodyMValue[k][TmpNbrInteraction][l] = TmpMIndices[l];
			    this->NBodyNValue[k][TmpNbrInteraction][l] = TmpNIndices[l];
			  }
			this->NBodyInteractionFactors[k][TmpNbrInteraction] = TmpCoef * SortedIndicesPerSumSymmetryFactor[MinSum][j];
			++TmpNbrInteraction;
		      }
		  }
	      }
	    for (int MinSum = 0; MinSum <= MaxSum; ++MinSum)
	      {
		for (int i = 0; i < NbrSortedIndicesPerSum[MinSum]; ++i)
		  delete[] SortedIndicesPerSum[MinSum][i];
		delete[] SortedIndicesPerSum[MinSum];
	      }
	    delete[] NbrSortedIndicesPerSum;
	    delete[] SortedIndicesPerSum;
	    delete[] SortedIndicesPerSumSymmetryFactor;
	    delete[] TmpInteractionCoeffients;
	    this->NBodyNbrInteractionFactors[k] = TmpNbrInteraction;
	  }
    }
  delete[] TmpNormalizationCoeffients;
}

