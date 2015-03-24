////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                        hollowcore n-body interaction                       //
//                                                                            //
//                        last modification : 17/01/2015                      //
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
#include "Hamiltonian/ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <algorithm>
#include <set>


// default constructor
//

ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian::ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// nbrNBody = value of the n (i.e. the n-body interaction)
// regenerateElementFlag = regenerate th interaction matrix elements instead of reading them from the harddrive
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian::ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio, int nbrNBody, bool regenerateElementFlag, AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->NBodyValue = nbrNBody;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->TwoBodyFlag = false;
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;
  this->OneBodyInteractionFactors = 0;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;
  this->EvaluateExponentialFactors();
  
  
  this->NbrEntryPrecalculatedInteractionCoefficients1 = this->MaxMomentum;
  for (int i = 1; i < this->NBodyValue; ++i)
    this->NbrEntryPrecalculatedInteractionCoefficients1 *= this->MaxMomentum;
  this->PrecalculatedInteractionCoefficients = new double* [this->NbrEntryPrecalculatedInteractionCoefficients1];
  this->NbrEntryPrecalculatedInteractionCoefficients2 = this->NBodyValue*(2*this->NBodyValue - 1);
  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
    {      
      this->PrecalculatedInteractionCoefficients[m1] = new double [this->NbrEntryPrecalculatedInteractionCoefficients2];
    }
  char* InteractionCoefficientFileName = new char [512];
  sprintf (InteractionCoefficientFileName, "%dbody_hollowcore_interactioncoefficient_2s_%d_ratio_%.10f.dat", this->NBodyValue, this->NbrLzValue, this->Ratio);
  if ((regenerateElementFlag == false) && (IsFile(InteractionCoefficientFileName)))
    {
      ifstream File;
      File.open(InteractionCoefficientFileName, ios::binary | ios::in);
      if (!File.is_open())
	{
	  cout << "cannot open " << InteractionCoefficientFileName << endl;
	}
      else
	{
	  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
	    ReadBlockLittleEndian(File, this->PrecalculatedInteractionCoefficients[m1], this->NbrEntryPrecalculatedInteractionCoefficients2);
	  File.close();
	}
    }
  else
    {
      ofstream File;
      File.open(InteractionCoefficientFileName, ios::binary | ios::out);
      if (!File.is_open())
	{
	  cout << "cannot create " << InteractionCoefficientFileName << endl;
	}
      else
	{ 
	  cout << "this->NbrEntryPrecalculatedInteractionCoefficients1 = "  << this->NbrEntryPrecalculatedInteractionCoefficients1 << endl;
	  
	  this->EvaluateInteractionCoefficientCreationUsingSymmetries();
// 	  int* TmpIndices;
// 	  int g;
// 	  int momentumTransfer;
// 	  for (int index = 0; index < this->NBodyValue * (2*this->NBodyValue - 1)* pow(this->NbrLzValue, this->NBodyValue); ++index)
// 	  {
// 	    TmpIndices = this->GetIndicesFromLinearizedIndex(index, g, momentumTransfer);
// 	    cout << index << " " << this->EvaluateLinearizedIndex (TmpIndices, g, momentumTransfer) << endl;
// // 	    for (int i = 0; i < this->NBodyValue ; ++i)
// // 	      cout << TmpIndices [i] << " " ;
// // 	    cout << g << " " << momentumTransfer << endl;
// 	  }
// 	  return;
// 	  this->EvaluateInteractionCoefficientCreation();
	  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
	  {
	    WriteBlockLittleEndian(File, this->PrecalculatedInteractionCoefficients[m1], this->NbrEntryPrecalculatedInteractionCoefficients2);
	  }

	  delete[] InteractionCoefficientFileName;
	  File.close();
	  
	  
	}
      
    }
  
  
  
  this->HamiltonianShift = 0.0;
  this->EvaluateInteractionFactors();
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  cout << "fast memory = ";
	  PrintMemorySize(cout,TmpMemory)<<endl;
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

ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian::~ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian()
{
}
  
// evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term (factor corresponding to the creation or the annihilation operators only) for each integer modulo the NBodyValue
//
// mIndices = array that contains the creation indices
// tmpSum = sum of all creation indices
// momentumTransfer = momentum transfer operated by the \prod_i a+_mi \prod_j a_nj operator, in units of the number of flux quanta
// return value = array of numerical coefficients  

double ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian::EvaluateIndividualInteractionCoefficientCreation(int* mIndices, int g, int momentumTransfer)
{
  double Coefficient;
  int AllDifferentIndicesFlag = 1;
  for (int i = 0; i < this->NBodyValue; ++i)
  {
    Coefficient = 0.0;
    for (int j = i + 1; j < this->NBodyValue; ++j)
      {
	AllDifferentIndicesFlag *= (mIndices[i] - mIndices[j]);
      }
  }
  if (AllDifferentIndicesFlag == 0)
    return Coefficient;
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  int tmpSum = 0;
  for (int i = 0; i < this->NBodyValue; ++i)
  {
    tmpSum += mIndices[i];
  }
  
  int* momFactor = new int [this->NBodyValue - 1];
  int* TmpIndices = new int [(this->NBodyValue - 1)];
  int* countIter = new int[(this->NBodyValue - 1)];
  for (int i = 0; i < this->NBodyValue - 1; ++i)
    momFactor[i] = this->NBodyValue * mIndices[i] -  tmpSum; 
  
  for (int j = 0; j < this->NBodyValue - 1; ++j)
    {
      int coef1 = (momFactor[j] / this->NbrLzValue) - (momFactor[j] / this->NbrLzValue) % this->NBodyValue;   
      TmpIndices[j] = (double) (((g + momentumTransfer) % this->NBodyValue) - coef1);	
      countIter[j] = 0;
    }
  Coefficient = this->EvaluateGaussianSum(0, TmpIndices, 0, countIter, momFactor);

  delete[] TmpIndices;
  delete[] momFactor;
  delete[] countIter;
  
  return Coefficient;
}

// evaluate the numerical coefficient in front of each \prod_i a+_mi \prod_j a_nj coupling term (factor corresponding to the creation or the annihilation operators only) for each integer modulo the NBodyValue
//
void ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficientCreation()
{
 //do not use symmetries to compute creation elements
 cout << "Generate creation elements without using the symmetries" << endl;
 for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
 {
    int* mIndices = new int [this->NBodyValue];
    int Tmp = m1;
    for (int i = 0; i < this->NBodyValue; ++i)
    {
      mIndices[i] = Tmp % this->MaxMomentum;
      Tmp /= this->MaxMomentum;
    }
    int TmpIndex = 0;
    for (int j = 0; j < 2*this->NBodyValue - 1; ++j)
    {
	int momentumTransfer = j - this->NBodyValue + 1;  
	int AllDifferentIndicesFlag = 1;
	for (int i = 0; i < this->NBodyValue; ++i)
	{
	  this->PrecalculatedInteractionCoefficients[m1][TmpIndex] = 0.0;
	  for (int k = i + 1; k < this->NBodyValue; ++k)
	    AllDifferentIndicesFlag *= (mIndices[i] - mIndices[k]);
	}
    if (AllDifferentIndicesFlag != 0)
    {
      for (int g = 0; g < this->NBodyValue; ++g)
      {
	int Index = this->EvaluateLinearizedIndex(mIndices, g, momentumTransfer);
	int sign;
	int CanonicalIndex = this->GetCanonicalIndex(Index, sign);
	int TmpG;
	int TmpMomentumTransfer;
	int* TmpIndices = this->GetIndicesFromLinearizedIndex (CanonicalIndex, TmpG, TmpMomentumTransfer);
// 	cout << "(" << Index << ") " ;
// 	for (int l = 0; l < this->NBodyValue; ++l)
// 	  cout << mIndices[l] << " " ;
// 	cout << g << " " << momentumTransfer << " | " ;
// 	cout << "(" << CanonicalIndex << ") " ;
// 	for (int l = 0; l < this->NBodyValue; ++l)
// 	  cout << TmpIndices[l] << " " ;
// 	cout << TmpG << " " << TmpMomentumTransfer << " : "  ;
	
	this->PrecalculatedInteractionCoefficients[m1][TmpIndex] = sign * this->EvaluateIndividualInteractionCoefficientCreation(TmpIndices, TmpG, TmpMomentumTransfer);
	double coefficient = this->EvaluateIndividualInteractionCoefficientCreation(mIndices, g, momentumTransfer);
	if (TmpIndex != (this->NBodyValue*j + g))
	  cout << TmpIndex << " " << (this->NBodyValue*j + g) << endl;
	++TmpIndex;
      }
    }
    else
	TmpIndex += this->NBodyValue;
    }
   }
}

// evaluate the numerical coefficient in front of each \prod_i a+_mi \prod_j a_nj coupling term (factor corresponding to the creation or the annihilation operators only) for each integer modulo the NBodyValue
//
void ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficientCreationUsingSymmetries()
{
  long NbrCreationElements = pow(this->NbrLzValue, this->NBodyValue) * this->NBodyValue * (2*this->NBodyValue - 1);
  int* sign = new int [NbrCreationElements];
  int nbrCanonicalElements = 0;
  int* mapIndexToCanonical = new int [NbrCreationElements];
  for (int i = 0; i < NbrCreationElements; ++i)
  {
     int TmpG;
     int TmpMomentumTransfer;
     int TmpG1;
     int TmpMomentumTransfer1;
//      int* TmpIndices1 = this->GetIndicesFromLinearizedIndex(i, TmpG, TmpMomentumTransfer);
     int canonicalIndex = this->GetCanonicalIndex(i, sign[i]);
     mapIndexToCanonical[i] = canonicalIndex;
     if (canonicalIndex == i)
       nbrCanonicalElements += 1;
// 	     int* TmpIndices2 = this->GetIndicesFromLinearizedIndex(canonicalIndex, TmpG1, TmpMomentumTransfer1);

// 	     cout << TmpIndices1[0] << " " << TmpIndices1[1] << " " << TmpG << " " << TmpMomentumTransfer << " | " << TmpIndices2[0] << " " << TmpIndices2[1] << " " << TmpG1 << " " << TmpMomentumTransfer1 << endl;
  }
  
  int* canonicalIndexList = new int [nbrCanonicalElements];
  int TmpIndex = 0;
  for (int i = 0; i < NbrCreationElements; ++i)
  {
//     cout << i << " " << mapIndexToCanonical[i] << endl;
     if (mapIndexToCanonical[i] == i)
     {
       canonicalIndexList[TmpIndex] = i;
       ++TmpIndex;       
     }
  }
  
  double* canonicalInteractionCoefficients = new double [nbrCanonicalElements];
  for (int i = 0; i < nbrCanonicalElements; ++i)
  {
    int TmpG;
    int TmpMomentumTransfer;
    int* TmpIndices1 = this->GetIndicesFromLinearizedIndex(canonicalIndexList[i], TmpG, TmpMomentumTransfer);
    canonicalInteractionCoefficients[i] = this->EvaluateIndividualInteractionCoefficientCreation(TmpIndices1, TmpG, TmpMomentumTransfer);
  }

  
  for (int m1 = 0; m1 < this->NbrEntryPrecalculatedInteractionCoefficients1; ++m1)
    {
      int* mIndices = new int [this->NBodyValue];
      int Tmp = m1;
      for (int i = 0; i < this->NBodyValue; ++i)
      {
	mIndices[i] = Tmp % this->MaxMomentum;
	Tmp /= this->MaxMomentum;
      }
      int TmpIndex = 0;
      for (int j = 0; j < 2*this->NBodyValue - 1; ++j)
      {
	int momentumTransfer = j - this->NBodyValue + 1;
	for (int g = 0; g < this->NBodyValue; ++g)
	{
	  int linearizedIndex = this->EvaluateLinearizedIndex(mIndices, g, momentumTransfer);
	  int canonicalIndex = mapIndexToCanonical[linearizedIndex];
	  int reducedCanonicalIndex = -1;
	  int l = 0;
	  while ((l < NbrCreationElements) && (reducedCanonicalIndex == -1))
	  {
	    if (canonicalIndexList[l] == canonicalIndex)
	      reducedCanonicalIndex = l;
	    ++l;
	  }
// 		  int TmpG;
// 		  int TmpMomentumTransfer;
// 		  int* TmpIndices1 = this->GetIndicesFromLinearizedIndex(linearizedCoefficient, TmpG, TmpMomentumTransfer);
// 		  cout << linearizedCoefficient << " " ;
// 		  for (int i =  0; i < this->NBodyValue; ++i)
// 		    cout << (TmpIndices1[i] - TmpIndices[m1][i]) << " " ;
// 		  cout << (TmpG - g) << " " << (TmpMomentumTransfer - momentumTransfer) << endl;
		  
// 	  cout << this->NbrEntryPrecalculatedInteractionCoefficients2 << " " << NbrCreationElements << " " << nbrCanonicalElements << " | " ;
// 	  cout << TmpIndex << " " << linearizedIndex << " " << reducedCanonicalIndex << endl;
	  this->PrecalculatedInteractionCoefficients[m1][TmpIndex] = sign[linearizedIndex] * canonicalInteractionCoefficients[reducedCanonicalIndex];
	  if (TmpIndex != (this->NBodyValue*j + g))
	    cout << TmpIndex << " " << (this->NBodyValue*j + g) << endl;
	  ++TmpIndex;
	}
      }
    }
}


// evaluate the N nested infinite sums of EvaluateInteractionCoefficientCreation
//
// nBodyValue = current index being incremented 
//TmpIndices = array containing the current value of all indices
// Sum = current value of the sum
// countIter = array of integers giving, for each TmpIndices[i], the number of times it has been incremented
// momFactor = array of indices that contains the information of the creation (or annihilation) indices
//return value = value of the coefficient
double ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian::EvaluateGaussianSum(int nBodyValue, int* TmpIndices, double Sum, int* countIter, int* momFactor)
{
  if (nBodyValue == this->NBodyValue - 1)
    return 0.0;
  
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  double normalizationCoefficient = pow(DoubleNbrLzValue,((double) (this->NBodyValue + 1)) / 4.0) * M_PI;
//   normalizationCoefficient = 1.0;
  double PIOnM = M_PI / DoubleNbrLzValue ;
  double Factor = 2.0*M_PI*this->Ratio / (DoubleNbrLzValue * ((double)(this->NBodyValue))* ((double)(this->NBodyValue)));
  normalizationCoefficient *= pow((sqrt(2.0 * (this->NBodyValue) * Factor)), ((double) (this->NBodyValue * (this->NBodyValue - 1) / 2)));
  int MinIter = 6;
  int TmpIndex = TmpIndices[0];
  double ExpFactor;
  double polynomialFactor;
  double sumRelativeMomenta;

  
  for (int i = nBodyValue; i < this->NBodyValue - 1; ++i)
    {
      TmpIndices[i] = TmpIndex;
      countIter[i] = 0;
    }
    
  ExpFactor = 0.0;
  polynomialFactor = 1.0;
  sumRelativeMomenta = 0.0;
  for (int j = 0; j < this->NBodyValue - 1; ++j)
    sumRelativeMomenta += ((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue;
  
  for (int j = 0; j < this->NBodyValue - 1; ++j)
    for (int k = j; k < this->NBodyValue - 1; ++k)
      ExpFactor += (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue)*(((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue);
  for (int j = 0; j < this->NBodyValue - 1; ++j)
  {
    for (int k = j + 1; k < this->NBodyValue - 1; ++k)
    {
      polynomialFactor *= ((((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue) - (((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue));
    }      
    polynomialFactor *= (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue + sumRelativeMomenta);
  }
  double Coefficient = normalizationCoefficient * polynomialFactor * exp(-Factor*ExpFactor);
  
  while ((abs(Coefficient) + abs(Sum) != abs(Sum)) || (countIter[nBodyValue] < MinIter))
  {
    countIter[nBodyValue] += 1;
    Sum += this->EvaluateGaussianSum(nBodyValue + 1, TmpIndices, 0.0, countIter, momFactor);
    if (nBodyValue == this->NBodyValue - 2)
      Sum += Coefficient;
    TmpIndices[nBodyValue] += this->NBodyValue;
    for (int i = nBodyValue + 1; i < this->NBodyValue - 1; ++i)
      TmpIndices[i] = TmpIndex;
    
    ExpFactor = 0.0;
    polynomialFactor = 1.0;
    sumRelativeMomenta = 0.0;
    for (int j = 0; j < this->NBodyValue - 1; ++j)
      sumRelativeMomenta += ((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue;
  
    for (int j = 0; j < this->NBodyValue - 1; ++j)
      for (int k = j; k < this->NBodyValue - 1; ++k)
	ExpFactor += (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue)*(((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue);
    for (int j = 0; j < this->NBodyValue - 1; ++j)
    {
      for (int k = j + 1; k < this->NBodyValue - 1; ++k)
      {
	polynomialFactor *= ((((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue) - (((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue));
      }      
      polynomialFactor *= (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue + sumRelativeMomenta);
    }
    Coefficient = normalizationCoefficient * polynomialFactor * exp(-Factor*ExpFactor);
    
  }
  
  TmpIndices[nBodyValue] = TmpIndex - this->NBodyValue;
  countIter[nBodyValue] = 0;
  for (int i = nBodyValue + 1; i < this->NBodyValue - 1; ++i)
  {
    TmpIndices[i] = TmpIndex;
    countIter[i] = 0;
  }
    
  ExpFactor = 0.0;
  polynomialFactor = 1.0;
  sumRelativeMomenta = 0.0;
  for (int j = 0; j < this->NBodyValue - 1; ++j)
    sumRelativeMomenta += ((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue;
  
  for (int j = 0; j < this->NBodyValue - 1; ++j)
    for (int k = j; k < this->NBodyValue - 1; ++k)
      ExpFactor += (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue)*(((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue);
  for (int j = 0; j < this->NBodyValue - 1; ++j)
  {
    for (int k = j + 1; k < this->NBodyValue - 1; ++k)
    {
      polynomialFactor *= ((((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue) - (((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue));
    }      
    polynomialFactor *= (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue + sumRelativeMomenta);
  }
  Coefficient = normalizationCoefficient * polynomialFactor * exp(-Factor*ExpFactor);
  
  
  while ((abs(Coefficient) + abs(Sum) != abs(Sum)) || (countIter[nBodyValue] < MinIter))
  {
    countIter[nBodyValue] += 1;
    Sum += this->EvaluateGaussianSum(nBodyValue + 1, TmpIndices, 0.0, countIter, momFactor);
    if (nBodyValue == this->NBodyValue - 2)
      Sum += Coefficient;
    TmpIndices[nBodyValue] -= this->NBodyValue;
    for (int i = nBodyValue + 1; i < this->NBodyValue - 1; ++i)
      TmpIndices[i] = TmpIndex;
    ExpFactor = 0.0;
    polynomialFactor = 1.0;
    sumRelativeMomenta = 0.0;
    for (int j = 0; j < this->NBodyValue - 1; ++j)
      sumRelativeMomenta += ((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue;
  
    for (int j = 0; j < this->NBodyValue - 1; ++j)
      for (int k = j; k < this->NBodyValue - 1; ++k)
	ExpFactor += (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue)*(((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue);
    for (int j = 0; j < this->NBodyValue - 1; ++j)
    {
      for (int k = j + 1; k < this->NBodyValue - 1; ++k)
      {
	polynomialFactor *= ((((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue) - (((double) (momFactor[k])) + ((double) (TmpIndices[k])) * DoubleNbrLzValue));
      }      
      polynomialFactor *= (((double) (momFactor[j])) + ((double) (TmpIndices[j])) * DoubleNbrLzValue + sumRelativeMomenta);
    }
    Coefficient = normalizationCoefficient * polynomialFactor * exp(-Factor*ExpFactor);
  }
  
  return Sum;
  
}
  

// Evaluate gaussian sum for the three body interaction (for test purposes only)
//
// momFactor = array of indices that contains the information of the creation (or annihilation) indices
// TmpIndices = array of indices that gives the initial indices that will be incremented in the sum
// return value = value of the sum

double ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian::EvaluateGaussianSum(int* momFactor, int* TmpIndices)
{
  double Sum = 0;
  double Precision = 1.0;
  int MinIter = 5;
  int countIter1 = 0;
  int countIter2 = 0;
  
  
  int TmpIndex1 = TmpIndices[0];
  int TmpIndex2 = TmpIndices[1];
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  double PIOnM = M_PI / DoubleNbrLzValue ;
  double Factor = 2.0*M_PI*this->Ratio / (DoubleNbrLzValue * ((double)(this->NBodyValue))* ((double)(this->NBodyValue)));
  double ExpFactor1;
  double ExpFactor2;
  
  ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
  ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
      
  double Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
  
  
  while ((Sum + Precision*Coefficient != Sum) || (countIter1 < MinIter))
    {
      countIter1 += 1;
      countIter2 = 0;
      while ((Sum + Precision*Coefficient != Sum) || (countIter2 < MinIter))
	{
	  countIter2 += 1;
	  Sum += Coefficient;
	  TmpIndices[0] += this->NBodyValue;
	  ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
	  ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
	  Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
	  
	}
      
      TmpIndices[0] = TmpIndex1 - this->NBodyValue;
      countIter2 = 0;
      ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
      ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
      Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
      
      while ((Sum + Precision*Coefficient != Sum) || (countIter2 < MinIter))
	{
	  countIter2 += 1;
	  Sum += Coefficient;
	  TmpIndices[0] -= this->NBodyValue;
	  ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
	  ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
	  Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
	  
	}
      
      TmpIndices[0] = TmpIndex1;
      TmpIndices[1] += this->NBodyValue;
      ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
      ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
      Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
    }
  
  countIter1 = 0;
  TmpIndices[0] = TmpIndex1;
  TmpIndices[1] = TmpIndex2 - this->NBodyValue;
  ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
  ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
  Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
  
  while ((Sum + Precision*Coefficient != Sum) || (countIter1 < MinIter))
    {
      countIter1 += 1;
      countIter2 = 0;
      while ((Sum + Precision*Coefficient != Sum) || (countIter2 < MinIter))
	{
	  countIter2 += 1;
	  Sum += Coefficient;
	  TmpIndices[0] += this->NBodyValue;
	  ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
	  ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
	  Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
	  
	}
      
      TmpIndices[0] = TmpIndex1 - this->NBodyValue;
      countIter2 = 0;
      ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
      ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
      Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
      
      while ((Sum + Precision*Coefficient != Sum) || (countIter2 < MinIter))
	{
	  countIter2 += 1;
	  Sum += Coefficient;
	  TmpIndices[0] -= this->NBodyValue;
	  ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
	  ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
	  Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
	  
	}
      
      TmpIndices[0] = TmpIndex1;
      TmpIndices[1] -= this->NBodyValue;
      ExpFactor1 = ((double) (momFactor[0])) + ((double) (TmpIndices[0])) * DoubleNbrLzValue;
      ExpFactor2 = ((double) (momFactor[1])) + ((double) (TmpIndices[1])) * DoubleNbrLzValue;
      Coefficient = exp(-Factor*(ExpFactor1*ExpFactor2 + ExpFactor1*ExpFactor1 + ExpFactor2*ExpFactor2));
    }
  
  return Sum;
}


	
// evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term
//
// creationCoefficients = array that contains the creation coefficients
// annihilationCoefficients = array that contains the annihilation coefficients
// nbrPermutations1 = number of permutations of the creation indexes
// nbrPermutations2 = number of permutations of the annihilation indexes
// return value = numerical coefficient  

double ParticleOnTorusNBodyHollowCoreWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(double* creationCoefficients, double* annihilationCoefficients)
{
  double Coefficient = 0.0;
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  FactorialCoefficient FactorialNBody;
  FactorialNBody.SetToOne();
  FactorialNBody.FactorialMultiply(this->NBodyValue);

 
  for (int j = 0; j < this->NBodyValue; ++j)
    Coefficient += creationCoefficients[j] * annihilationCoefficients[j];

  for (int i = 0; i < this->NBodyValue; ++i)
    Coefficient /= (M_PI * DoubleNbrLzValue);
  
  return (Coefficient/ (FactorialNBody.GetNumericalValue()));
}

  
