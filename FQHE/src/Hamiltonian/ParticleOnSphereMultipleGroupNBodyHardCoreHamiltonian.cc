////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//               forcing vanishing properties with multiple groups            //
//                        of n-body hard core interaction                     //
//                                                                            //
//                        last modification : 04/06/2011                      //
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
#include "Hamiltonian/ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"

#include <stdio.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


// default constructor
//

ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian::ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// nbrGroups = number of groups that have a given vanishing property
// nbrNBodys = number of particle that interact simultaneously through the hard core interaction for the each group
// l2Factor = multiplicative factor in front of an additional L^2 operator in the Hamiltonian (0 if none)
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian::ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, int nbrGroups, int* nbrNBodys, double l2Factor,
													     AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag,
													     char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;

  this->OneBodyTermFlag = false;
  this->FullTwoBodyFlag = false;

  this->NbrGroups = nbrGroups;
  this->NbrNBodys = new int [this->NbrGroups];
  this->NbrNbody = 0;
  for (int i = 0; i < this->NbrGroups; ++i)
    {
      this->NbrNBodys[i] = nbrNBodys[i];
      this->NbrNbody += this->NbrNBodys[i];
    }
  this->MaxNBody = this->NbrNbody;
  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NBodyInteractionFactors = new double** [this->MaxNBody + 1];
  this->NBodyInteractionWeightFactors = new double [this->MaxNBody + 1];
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
      this->NBodyInteractionWeightFactors[k] = 0.0;
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic) && ((k & 1) == 0))
	{
	  this->NBodySign[k] = -1.0;
	  if (k == 4)
	    this->NBodySign[k] = 1.0;
	}
      this->NbrNIndices[k] = 0;
      this->NIndices[k] = 0;
      this->NbrMIndices[k] = 0;
      this->MIndices[k] = 0;
      this->MNNBodyInteractionFactors[k] = 0;
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

  if (l2Factor != 0.0)
    {
      this->L2Operator = new ParticleOnSphereL2Hamiltonian(this->Particles, this->NbrParticles, this->LzMax, this->Particles->GetLzValue() , this->Architecture, l2Factor);
    }
  else
    {
      this->L2Operator = 0;
    }

}


// destructor
//

ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian::~ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian()
{
  delete[] this->NbrNBodys;
}

// evaluate all interaction factors
//   

void ParticleOnSphereMultipleGroupNBodyHardCoreHamiltonian::EvaluateInteractionFactors()
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
      cout << "fermionic case is not implemented" << endl;
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


	    long TmpNbrNIndices = 0;
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	      TmpNbrNIndices += this->NbrSortedIndicesPerSum[k][MinSum];
	    this->NbrNIndices[k] = TmpNbrNIndices;
	    this->NIndices[k] = new int[TmpNbrNIndices * k];
	    this->NbrMIndices[k] = new long[TmpNbrNIndices];
	    this->MIndices[k] = new int*[TmpNbrNIndices];
	    this->MNNBodyInteractionFactors[k] = new double* [TmpNbrNIndices];
	    TmpNbrNIndices = 0;	 
	    int* TmpNIndices = this->NIndices[k];
	    for (int MinSum = 0; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	      {
		int Lim = this->NbrSortedIndicesPerSum[k][MinSum];
		double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
		int* TmpNIndices2 = this->SortedIndicesPerSum[k][MinSum];
		for (int i = 0; i < Lim; ++i)
		  {
		    this->NbrMIndices[k][TmpNbrNIndices] = Lim;		    
		    this->MIndices[k][TmpNbrNIndices] = new int [Lim * k];
		    this->MNNBodyInteractionFactors[k][TmpNbrNIndices] = new double [Lim];
		    int* TmpMIndices = this->MIndices[k][TmpNbrNIndices];
		    int* TmpMIndices2 = this->SortedIndicesPerSum[k][MinSum];
		    double* TmpInteraction = this->MNNBodyInteractionFactors[k][TmpNbrNIndices];
		    Coefficient = TmpSymmetryFactors[i] * this->NBodyInteractionWeightFactors[k] * TmpInteractionCoeffients[MinSum] * TmpInteractionCoeffients[MinSum];
		    for (int l = 0; l < k; ++l)
		      Coefficient *= TmpNormalizationCoeffients[TmpNIndices2[l]];		    		      
		    for (int j = 0; j < Lim; ++j)
		      {
			double Coefficient2 = TmpSymmetryFactors[j];
			for (int l = 0; l < k; ++l)
			  {
			    Coefficient2 *= TmpNormalizationCoeffients[(*TmpMIndices2)];		    
			    (*TmpMIndices) = (*TmpMIndices2);			
			    ++TmpMIndices;
			    ++TmpMIndices2;
			  }			
			TmpInteraction[j] = Coefficient * Coefficient2;
		      }
		    for (int j = 0; j < k; ++j)
		      {
			(*TmpNIndices) = (*TmpNIndices2);			
			++TmpNIndices;
			++TmpNIndices2;
		      }
		    ++TmpNbrNIndices;
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
  if (this->FullTwoBodyFlag == true)
    {
      int Lim;
      int Min;
      int Pos = 0;
      ClebschGordanCoefficients Clebsch (this->LzMax, this->LzMax);
      int J = 2 * this->LzMax - 2;
      int m4;
      double ClebschCoef;
      double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
      
      int Sign = 1;
      if (this->LzMax & 1)
	Sign = 0;
      double MaxCoefficient = 0.0;
      
      if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
	{
	  for (int m1 = -this->LzMax; m1 <= this->LzMax; m1 += 2)
	    for (int m2 =  -this->LzMax; m2 < m1; m2 += 2)
	      {
		Lim = m1 + m2 + this->LzMax;
		if (Lim > this->LzMax)
		  Lim = this->LzMax;
		Min = m1 + m2 - this->LzMax;
		if (Min < -this->LzMax)
		  Min = -this->LzMax;
		for (int m3 = Min; m3 <= Lim; m3 += 2)
		  {
		    Clebsch.InitializeCoefficientIterator(m1, m2);
		    m4 = m1 + m2 - m3;
		    TmpCoefficient[Pos] = 0.0;
		    while (Clebsch.Iterate(J, ClebschCoef))
		      {
			if (((J >> 1) & 1) == Sign)
			  TmpCoefficient[Pos] += this->PseudoPotential[J >> 1] * ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		      }
		    if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		      MaxCoefficient = TmpCoefficient[Pos];
		    ++Pos;
		  }
	      }
	  this->NbrInteractionFactors = 0;
	  this->M1Value = new int [Pos];
	  this->M2Value = new int [Pos];
	  this->M3Value = new int [Pos];
	  this->InteractionFactors = new double [Pos];
	  cout << "nbr interaction = " << Pos << endl;
	  Pos = 0;
	  MaxCoefficient *= MACHINE_PRECISION;
	  double Factor = - 4.0;// / sqrt (0.5 * ((double) this->LzMax));
	  this->NbrM12Indices = (this->NbrLzValue * (this->NbrLzValue - 1)) / 2;
	  this->M1Value = new int [this->NbrM12Indices];
	  this->M2Value = new int [this->NbrM12Indices];
	  this->NbrM3Values = new int [this->NbrM12Indices];
	  this->M3Values = new int* [this->NbrM12Indices];
	  int TotalIndex = 0;
	  Pos = 0;
	  for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
	    for (int m2 = 0; m2 < m1; ++m2)
	      {
		Lim = m1 + m2;
		if (Lim > this->LzMax)
		  Lim = this->LzMax;
		Min = m1 + m2 - this->LzMax;
		if (Min < 0)
		  Min = 0;
		this->M1Value[TotalIndex] = m1;
		this->M2Value[TotalIndex] = m2;	    
		this->NbrM3Values[TotalIndex] = 0;
		for (int m3 = Min; m3 <= Lim; ++m3)
		  if ((2 * m3) > (m1 + m2))
		    ++this->NbrM3Values[TotalIndex];
		if (this->NbrM3Values[TotalIndex] > 0)
		  {
		    this->M3Values[TotalIndex] = new int [this->NbrM3Values[TotalIndex]];
		    int TmpIndex = 0;
		    for (int m3 = Min; m3 <= Lim; ++m3)
		      {
			if ((2 * m3) > (m1 + m2))
			  {
			    this->M3Values[TotalIndex][TmpIndex] = m3;
			    this->InteractionFactors[this->NbrInteractionFactors] = Factor * TmpCoefficient[Pos];
			    ++this->NbrInteractionFactors;
			    ++TmpIndex;
			  }
			++Pos;
		      }
		  }
		++TotalIndex;
	      }
	}
      else
	{
	  for (int m1 = -this->LzMax; m1 <= this->LzMax; m1 += 2)
	    for (int m2 =  -this->LzMax; m2 <= m1; m2 += 2)
	      {
		Lim = m1 + m2 + this->LzMax;
		if (Lim > this->LzMax)
		  Lim = this->LzMax;
		Min = m1 + m2 - this->LzMax;
		if (Min < -this->LzMax)
		  Min = -this->LzMax;
		for (int m3 = Min; m3 <= Lim; m3 += 2)
		  {
		    Clebsch.InitializeCoefficientIterator(m1, m2);
		    m4 = m1 + m2 - m3;
		    TmpCoefficient[Pos] = 0.0;
		    while (Clebsch.Iterate(J, ClebschCoef))
		      {
			if (((J >> 1) & 1) != Sign)
			  TmpCoefficient[Pos] += this->PseudoPotential[J >> 1] * ClebschCoef * Clebsch.GetCoefficient(m3, m4, J);
		      }
		    if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		      MaxCoefficient = TmpCoefficient[Pos];
		    ++Pos;
		  }
	      }
	  this->NbrInteractionFactors = 0;
	  this->M1Value = new int [Pos];
	  this->M2Value = new int [Pos];
	  this->M3Value = new int [Pos];
	  this->InteractionFactors = new double [Pos];
	  cout << "nbr interaction = " << Pos << endl;
	  Pos = 0;
	  MaxCoefficient *= MACHINE_PRECISION;
	  double Factor = 4.0 / sqrt (0.5 * ((double) this->LzMax));
	  for (int m1 = 0; m1 < this->NbrLzValue; ++m1)
	    {
	      for (int m2 = 0; m2 < m1; ++m2)
		{
		  Lim = m1 + m2;
		  if (Lim > this->LzMax)
		    Lim = this->LzMax;
		  Min = m1 + m2 - this->LzMax;
		  if (Min < 0)
		    Min = 0;
		  for (int m3 = Min; m3 <= Lim; ++m3)
		    {
		      if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
			{
			  if ((2 * m3) > (m1 + m2))
			    {
			      this->InteractionFactors[this->NbrInteractionFactors] = Factor * TmpCoefficient[Pos];
			      this->M1Value[this->NbrInteractionFactors] = m1;
			      this->M2Value[this->NbrInteractionFactors] = m2;
			      this->M3Value[this->NbrInteractionFactors] = m3;
			      ++this->NbrInteractionFactors;
			    }
			  else
			    if ((2 * m3) == (m1 + m2))
			      {
				this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * Factor * TmpCoefficient[Pos];
				this->M1Value[this->NbrInteractionFactors] = m1;
				this->M2Value[this->NbrInteractionFactors] = m2;
				this->M3Value[this->NbrInteractionFactors] = m3;
				++this->NbrInteractionFactors;
			      }
			}
		      ++Pos;
		    }
		}	
	      Lim = 2 * m1;
	      if (Lim > this->LzMax)
		Lim = this->LzMax;
	      Min = 2 * m1 - this->LzMax;
	      if (Min < 0)
		Min = 0;
	      for (int m3 = Min; m3 <= Lim; ++m3)
		{
		  if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      if (m3 > m1)
			{
			  this->InteractionFactors[this->NbrInteractionFactors] = 0.5 * Factor * TmpCoefficient[Pos];
			  this->M1Value[this->NbrInteractionFactors] = m1;
			  this->M2Value[this->NbrInteractionFactors] = m1;
			  this->M3Value[this->NbrInteractionFactors] = m3;
			  ++this->NbrInteractionFactors;
			}
		      else
			if (m3 == m1)
			  {
			    this->InteractionFactors[this->NbrInteractionFactors] = 0.25 * Factor * TmpCoefficient[Pos];
			    this->M1Value[this->NbrInteractionFactors] = m1;
			    this->M2Value[this->NbrInteractionFactors] = m1;
			    this->M3Value[this->NbrInteractionFactors] = m3;
			    ++this->NbrInteractionFactors;
			  }
		    }
		  ++Pos;
		}
	    }
	}
    }
}

