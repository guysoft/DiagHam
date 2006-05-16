////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a sphere with        //
//  n-body hard core interaction and two localized impurities at the poles    //
//                                                                            //
//                        last modification : 04/05/2006                      //
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
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "FunctionBasis/QHEFunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "Vector/RealVector.h"
#include "FunctionBasis/QHEFunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

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
// architecture = architecture to use for precalculation
// nbrBody = number of particle that interact simultaneously through the hard core interaction
// impurityPotential = potential strength associted to the impurities
// landauLevel = index of the Landau level where calculations have to be done (0 being the LLL)
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian::ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
														     int nbrBody, double impurityPotential, int landauLevel,
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

  this->ImpurityPotential = impurityPotential;
  this->LandauLevel = landauLevel;

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
  this->EvaluateInteractionFactorsWithTwoImpurities();
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
// impurityPotential = potential strength associted to the impurities
// landauLevel = index of the Landau level where calculations have to be done (0 being the LLL)
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian::ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian(ParticleOnSphere* particles, int nbrParticles, int lzmax, 
													       int maxNbrBody, double* nBodyFactors, double impurityPotential, int landauLevel,
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

  this->ImpurityPotential = impurityPotential;
  this->LandauLevel = landauLevel;

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
  this->EvaluateInteractionFactorsWithTwoImpurities();
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

ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian::~ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian()
{
}

// evaluate all interaction factors (including those arising from impurities)
//   

void ParticleOnSphereNBodyHardCoreWithTwoImpuritiesHamiltonian::EvaluateInteractionFactorsWithTwoImpurities()
{
  this->EvaluateInteractionFactors();
  if (this->ImpurityPotential < 0.0)
    this->NBodySign[1] *= -1.0;
  this->NBodyFlags[1] = true;
//   double TmpFactor = fabs(this->ImpurityPotential) * ((double) this->NbrLzValue) / (4.0 * M_PI);
  int Shift = this->LandauLevel + 1;
//   this->NbrSortedIndicesPerSum[1] = new int [2 * Shift];
//   this->SortedIndicesPerSum[1] = new int* [2 * Shift];  
//   this->NBodyInteractionFactors[1] = new double* [2 * Shift];
//   this->MinSumIndices[1] = 0;
//   this->MaxSumIndices[1] = (2 * this->LandauLevel) + 1;
  FactorialCoefficient Coef;
  ParticleOnSphereGenericLLFunctionBasis Basis (this->LzMax - (2 * this->LandauLevel),this-> LandauLevel);
  RealVector Value(2, true);
  Complex TmpValue;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      Value[0] = 0.0;
      Basis.GetFunctionValue(Value, TmpValue, i);
      cout << TmpValue << " ";
      Value[0] = M_PI;
      Basis.GetFunctionValue(Value, TmpValue, i);
      cout << TmpValue << " ";     
    }
  cout << endl;
  
  double TmpFactor = sqrt (fabs(this->ImpurityPotential));
  this->NbrSortedIndicesPerSum[1] = new int [2];
  this->SortedIndicesPerSum[1] = new int* [2];  
  this->NBodyInteractionFactors[1] = new double* [2];
  this->MinSumIndices[1] = 0;
  this->MaxSumIndices[1] = 1;
  this->NbrSortedIndicesPerSum[1][0] = 1;
  this->SortedIndicesPerSum[1][0] = new int [1];
  this->NBodyInteractionFactors[1][0] = new double[1];
  this->SortedIndicesPerSum[1][0][0] = this->LandauLevel;
  Value[0] = M_PI;
  Basis.GetFunctionValue(Value, TmpValue, this->LandauLevel);
  cout << TmpValue << " " << Norm(TmpValue) << endl;     
  this->NBodyInteractionFactors[1][0][0] = TmpFactor * Norm(TmpValue);
  this->NbrSortedIndicesPerSum[1][1] = 1;
  this->SortedIndicesPerSum[1][1] = new int [1];
  this->NBodyInteractionFactors[1][1] = new double[1];
  this->SortedIndicesPerSum[1][1][0] = this->LzMax - this->LandauLevel;
  Value[0] = 0.0;
  Basis.GetFunctionValue(Value, TmpValue, this->LzMax - this->LandauLevel);
  cout << TmpValue << " " << Norm(TmpValue) << endl;     
  this->NBodyInteractionFactors[1][1][0] = TmpFactor * Norm(TmpValue);
  cout << this->NBodyInteractionFactors[1][0][0] << " " << this->NBodyInteractionFactors[1][1][0] << endl;


//   for (int i = 0; i <= this->LandauLevel; ++i) 
//     {
//       this->NbrSortedIndicesPerSum[1][i] = 1;
//       this->SortedIndicesPerSum[1][i] = new int [1];
//       this->NBodyInteractionFactors[1][i] = new double[1];
//       this->SortedIndicesPerSum[1][i][0] = this->LandauLevel - i;
//       Coef.SetToOne();
//       Coef.PartialFactorialMultiply(this->LandauLevel - i + 1,this->LandauLevel);      
//       Coef.FactorialDivide(i);
//       Coef.PartialFactorialMultiply(this->LzMax - this->LandauLevel + 1, this->LzMax + i - this->LandauLevel);      
//       Coef.FactorialDivide(i);
//       this->NBodyInteractionFactors[1][i][0] = sqrt (TmpFactor * Coef.GetNumericalValue());
//       Value[0] = 0.0;
//       Basis.GetFunctionValue(Value, TmpValue, this->LandauLevel - i);
//       cout << sqrt (Coef.GetNumericalValue()) << " " << TmpValue << endl;


//       this->NbrSortedIndicesPerSum[1][i + Shift] = 1;
//       this->SortedIndicesPerSum[1][i + Shift] = new int [1];
//       this->NBodyInteractionFactors[1][i + Shift] = new double[1];
//       this->SortedIndicesPerSum[1][i + Shift][0] = this->LzMax - i;
//       Coef.SetToOne();
//       Coef.PartialFactorialMultiply(this->LandauLevel - i + 1,this->LandauLevel);      
//       Coef.FactorialDivide(i);
//       Coef.PartialFactorialMultiply(this->LzMax - this->LandauLevel + 1, this->LzMax - i);      
//       Coef.FactorialDivide(this->LandauLevel - i);
//       this->NBodyInteractionFactors[1][i + Shift][0] = sqrt (TmpFactor * Coef.GetNumericalValue());
//       Value[0] = M_PI;
//       Basis.GetFunctionValue(Value, TmpValue, this->LzMax - i);
//       cout << sqrt (Coef.GetNumericalValue()) << " " << TmpValue << endl;
// //      this->NBodyInteractionFactors[1][i + Shift][0] = this->NBodyInteractionFactors[1][i][0];
//     }





//   this->NbrSortedIndicesPerSum[1] = new int [this->LzMax + 1];
//   this->SortedIndicesPerSum[1] = new int* [this->LzMax + 1];  
//   this->NBodyInteractionFactors[1] = new double* [this->LzMax + 1];
//   this->MinSumIndices[1] = 0;
//   this->MaxSumIndices[1] = this->LzMax;


//   for (int i = 0; i <= this->LzMax; ++i) 
//     {
//       this->NbrSortedIndicesPerSum[1][i] = 1;
//       this->SortedIndicesPerSum[1][i] = new int [1];
//       this->NBodyInteractionFactors[1][i] = new double[1];
//       this->SortedIndicesPerSum[1][i][0] = i;
//       this->NBodyInteractionFactors[1][i][0] = sqrt((double) i);
//     }

}



