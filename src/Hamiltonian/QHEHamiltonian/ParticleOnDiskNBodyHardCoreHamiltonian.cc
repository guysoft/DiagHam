////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated to particles on a disk with         //
//                         n-body hard core interaction                       //
//                                                                            //
//                        last modification : 04/09/2004                      //
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
#include "Hamiltonian/QHEHamiltonian/ParticleOnDiskNBodyHardCoreHamiltonian.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/FactorialCoefficient.h"

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

ParticleOnDiskNBodyHardCoreHamiltonian::ParticleOnDiskNBodyHardCoreHamiltonian(ParticleOnDisk* particles, int nbrParticles, int lzmax, int nbrBody,
									       AbstractArchitecture* architecture, long memory, 
									       char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = this->Particles->GetMaxLz();
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;

  this->NbrNbody = nbrBody;
  this->MaxNBody = this->NbrNbody;
  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NBodyInteractionFactors = new double** [this->MaxNBody + 1];
  this->NbrSortedIndicesPerSum = new int* [this->MaxNBody + 1];
  this->SortedIndicesPerSum = new int** [this->MaxNBody + 1];

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

ParticleOnDiskNBodyHardCoreHamiltonian::~ParticleOnDiskNBodyHardCoreHamiltonian()
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

void ParticleOnDiskNBodyHardCoreHamiltonian::EvaluateInteractionFactors()
{
  double* TmpNormalizationCoeffients = new double[this->NbrLzValue];
  double TmpFactor = ((double) this->NbrLzValue) / (4.0 * M_PI);
  double TmpBinomial = 1.0;
  TmpNormalizationCoeffients[0] = sqrt (TmpBinomial * TmpFactor);
  for (int i = 1; i < this->NbrLzValue; ++i)
    {
      TmpBinomial *= this->NbrLzValue - ((double) i) + 1.0;
      TmpBinomial /= ((double) i);
      TmpNormalizationCoeffients[i] = sqrt (TmpBinomial * TmpFactor);
    }
  if (this->Particles->GetParticleStatistic() == ParticleOnDisk::FermionicStatistic)
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
		    TmpNBodyInteractionFactors[i] = Coefficient;
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
	    this->MaxSumIndices = this->MaxMomentum * k;
	    double** SortedIndicesPerSumSymmetryFactor;
	    GetAllSymmetricIndices(this->NbrLzValue, k, this->NbrSortedIndicesPerSum[k], this->SortedIndicesPerSum[k],
				   SortedIndicesPerSumSymmetryFactor);
	    this->NBodyInteractionFactors[k] = new double* [MaxSumIndices + 1];
	    int Lim;
	    double Factor = 1.0;
	    for (int i = 1; i < k; ++i)
	      Factor *= 8.0 * M_PI;
	    Factor = 1.0 / Factor;
	    FactorialCoefficient Coef;
	    for (int MinSum = 0; MinSum <= MaxSumIndices; ++MinSum)
	      {
		Lim = this->NbrSortedIndicesPerSum[k][MinSum];
		double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
		this->NBodyInteractionFactors[k][MinSum] = new double [Lim];
		double* TmpNBodyInteractionFactors = this->NBodyInteractionFactors[k][MinSum];		  
		int* TmpMIndices = this->SortedIndicesPerSum[k][MinSum];
		for (int i = 0; i < Lim; ++i)
		  {
		    Coef.SetToOne();
		    Coef.FactorialMultiply(MinSum);
		    for (int l = 0; l < k; ++l)
		      Coef.FactorialDivide(TmpMIndices[l]);
		    Coef.PowerNDivide(k, MinSum);
		    TmpNBodyInteractionFactors[i] = sqrt(Coef.GetNumericalValue() * TmpSymmetryFactors[i] * Factor);
		    TmpMIndices += k;
		  }
	      }
	    for (int MinSum = 0; MinSum <= MaxSumIndices; ++MinSum)
	      {
		delete[] SortedIndicesPerSumSymmetryFactor[MinSum];
	      }
	    delete[] SortedIndicesPerSumSymmetryFactor;
	  }
    }
}

