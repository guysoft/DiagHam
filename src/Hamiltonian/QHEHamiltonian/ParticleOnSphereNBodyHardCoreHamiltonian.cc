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
/*
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

  int MaxSumIndices = this->MaxNBody * this->LzValue;
  double* TmpInteractionCoeffients = new double[MaxSumIndices + 1];
  TmpInteractionCoeffients[0] = sqrt (TmpBinomial * TmpFactor);
  for (int i = 1; i < MaxSumIndices; ++i)
    {
      TmpBinomial *= this->LzMax - ((double) i) + 1.0;
      TmpBinomial /= ((double) i);
      TmpInteractionCoeffients[i] = sqrt (TmpBinomial * TmpFactor);
    }

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      int MaxNbrMIndices = (this->MaxNBody * (this->MaxNBody - 1)) / 2; 
      double* TmpFactors = new double [MaxNbrMIndices * MaxNbrMIndices];
      int** TmpIndices = new int* [MaxNbrMIndices];
      for (int i = 0; i < MaxNbrMIndices; ++i)
	{
	  TmpIndices[i] = new int [this->NbrLzValue]
	}
      double* TmpInteractionFactors = new double[MaxNbrMIndices * MaxNbrMIndices];
      int* TmpMIndices;
      int* TmpNIndices;
      int SumIndices;
      int Pos;
      double Coefficient;
      long TmpNbrInteractionFactors;
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true) 
	  {
	    int NbrMIndices = (k * (k - 1)) / 2; 
	    int* TmpMIndices = TmpIndices[i];
	    Pos = 0;
	    double MaxCoefficient = 0.0;
	    for (int i = 0; i < NbrMIndices; ++i)
	      {
		TmpMIndices = TmpIndices[i];
		for (int j = 0; j < NbrMIndices; ++j)
		  {
		    TmpNIndices = TmpIndices[j];
		    SumIndices = 0;
		    for (int l = 0; l < k; ++l)
		      {
			SumIndices += TmpMIndices[l];
			SumIndices -= TmpNIndices[l];
		      }
		    if (SumIndices == 0)
		      {
			for (int l = 0; l < k; ++l)
			  SumIndices += TmpMIndices[l];
			double& TmpCoef = TmpInteractionFactors[Pos];
			TmpCoef = TmpInteractionCoeffients[SumIndices];			
			for (int l = 0; l < k; ++l)
			  {
			    TmpCoef *= TmpNormalizationCoeffients[TmpMIndices[l]];
			    TmpCoef *= TmpNormalizationCoeffients[TmpNIndices[l]];
			  }
			if ((MaxCoefficient == 0.0) || (MaxCoefficient < fabs(TmpCoef)))
			  {
			    MaxCoefficient = fabs(TmpCoef);
			  }
			++TmpNbrInteractionFactors;
		      }
		    else
		      {
			TmpInteractionFactors[Pos] = 0.0;
		      }
		    ++Pos;
		  }
	      }
	    this->NBodyMValue[k] = new int** [TmpNbrInteractionFactors];
	    this->NBodyNValue[k] = new int** [TmpNbrInteractionFactors];
	    this->NBodyInteractionFactors = new double [TmpNbrInteractionFactors];
	    MaxCoefficient *= MACHINE_PRECISION;
	    TmpNbrInteractionFactors = 0;
	    Pos = 0;
	    for (int i = 0; i < NbrMIndices; ++i)
	      {
		TmpMIndices = TmpIndices[i];
		for (int j = 0; j < NbrMIndices; ++j)
		  {
		    if (fabs(TmpInteractionFactors[Pos]) > MaxCoefficient)
		      {
			this->NBodyMValue[k][TmpNbrInteractionFactors]  = new int [k];
			this->NBodyNValue[k][TmpNbrInteractionFactors]  = new int [k];
			for (int l = 0; l < k; ++l)
			  {
			    this->NBodyMValue[k][TmpNbrInteractionFactors][l] = TmpMIndices[l];
			    this->NBodyNValue[k][TmpNbrInteractionFactors][l] = TmpIndices[j][l];
			  }
			this->NBodyInteractionFactors[TmpNbrInteractionFactors] = -TmpInteractionFactors[Pos];
			++TmpNbrInteractionFactors;
		      }
		    ++Pos;
		  }
	      }
	    this->NBodyNbrInteractionFactors[k] = TmpNbrInteractionFactors;
	  }
      delete[] TmpInteractionFactors;
    }
  else
    {
    }
  delete[] TmpNormalizationCoeffients;
  delete[] TmpInteractionCoeffients;
}
*/
