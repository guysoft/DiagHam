////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                          hardcore n-body interaction                       //
//                                                                            //
//                        last modification : 29/07/2014                      //
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


#ifndef PARTICLEONTORUSGENERICNBODYWITHMAGNETICTRANSLATIONSHAMILTONIAN_H
#define PARTICLEONTORUSGENERICNBODYWITHMAGNETICTRANSLATIONSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian.h"
#include "MathTools/FactorialCoefficient.h" 

#include <iostream>
#include <algorithm>


using std::ostream;



class ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian : public AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian
{

 protected:

  // array where the interaction coefficients are stored
  double** PrecalculatedInteractionCoefficients;
  // number of rows in PrecalculatedInteractionCoefficients
  int NbrEntryPrecalculatedInteractionCoefficients1;
  // number of columns in PrecalculatedInteractionCoefficients
  int NbrEntryPrecalculatedInteractionCoefficients2;

 public:

  // default constructor
  //
  ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // lzmax = maximum Lz value reached by a particle in the state
  // nbrNBody = value of the n (i.e. the n-body interaction)
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio, int nbrNBody,
								  AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
								  char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian();
  
  // clone hamiltonian without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractHamiltonian* Clone ();

	
  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);
  
 protected:
  
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors();
	
  // evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term
  //
  // creationCoefficients = array that contains the creation coefficients
  // annihilationCoefficients = array that contains the annihilation coefficients
  // return value = numerical coefficient  
  virtual double EvaluateInteractionCoefficient(double* creationCoefficients, double* annihilationCoefficients);
  
  // evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term in the case where they have been precalculated
  //
  //mIndices  = array that contains the creation indices
  //nIndices = array that contains the annihilation indices
  //return value = numerical coefficient
  virtual double EvaluateInteractionNIndexSymmetrizedCoefficient(int* mIndices, int* nIndices);

  // evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term (factor corresponding to the creation operators only) for each integer modulo the NBodyValue
  //
  // mIndices = array that contains the creation indices
  // momentumTransfer = momentum transfer operated by the \prod_i a+_mi \prod_j a_nj operator, in units of the number of flux quanta
  // return value = array of numerical coefficients 
  double* EvaluateInteractionCoefficientCreation(int* mIndices, int momentumTransfer);
  
  
  // evaluate the two nested Gaussian sum for a three body interaction (for test purposes)
  //
  // momFactor = array of indices that contains the information of the creation (or annihilation) indices
  // TmpIndices = array of indices that gives the initial indices that will be incremented in the sum
  // return value = value of the sum
  double EvaluateGaussianSum(int* momFactor, int* TmpIndices);
  
  
  // evaluate the N nested infinite sums of EvaluateInteractionCoefficientCreation
  //
  // nBodyValue = current index being incremented 
  //TmpIndices = array containing the current value of all indices
  // Sum = current value of the sum
  // countIter = array of integers giving, for each TmpIndices[i], the number of times it has been incremented
  // momFactor = array of indices that contains the information of the creation (or annihilation) indices
  //return value = value of the coefficient
  double EvaluateGaussianSum(int nBodyValue, int* TmpIndices, double Sum, int* countIter, int* momFactor);
  
   
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  virtual double EvaluateTwoBodyInteractionCoefficient(int m1, int m2, int m3, int m4);


};

inline double ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian::EvaluateInteractionNIndexSymmetrizedCoefficient(int* mIndices, int* nIndices)
{
  
  
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  FactorialCoefficient FactorialNBody;
  FactorialNBody.SetToOne();
  FactorialNBody.FactorialMultiply(this->NBodyValue);
  double FinalFactor = FactorialNBody.GetNumericalValue();
  for (int i = 0; i < this->NBodyValue; ++i)
    FinalFactor *= (M_PI * DoubleNbrLzValue);
 
  
 int* mIndices2 = new int[this->NBodyValue];
 int* nIndices2 = new int[this->NBodyValue];
 
 int momentumTransfer = 0;
 for (int i = 0; i < this->NBodyValue; ++i)
   momentumTransfer += (nIndices[i] - mIndices[i]);
 momentumTransfer /= this->NbrLzValue;
 int shiftedMomentumTransfer = momentumTransfer + this->NBodyValue - 1;
 int shift =  (this->NBodyValue - 1)*this->NBodyValue;
  
 double TmpInteraction; 
 int m1 = mIndices[0];
 int m2 = nIndices[0];
 
 int Factor = this->NbrLzValue;
 for (int k = 1; k < this->NBodyValue; ++k)
 {
  m1 += mIndices[k]*Factor;
  m2 += nIndices[k]*Factor;
  Factor *= this->NbrLzValue;
 }
 int m2Initial = m2;
 
 TmpInteraction = 0.0;
 for (int g = 0; g < this->NBodyValue; ++g)
 {
    TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][shift + g] * this->PrecalculatedInteractionCoefficients[m2][shiftedMomentumTransfer*this->NBodyValue + g];
 }
 
 for (int k = 0 ; k < this->NBodyValue; ++k)
   nIndices2[k] = nIndices[k];
 while (std::prev_permutation(nIndices2, nIndices2 + this->NBodyValue))
    {
      m2 = nIndices2[0];
      Factor = this->NbrLzValue;
      for (int k = 1 ; k < this->NBodyValue; ++k)
	{
	  m2 += nIndices2[k]*Factor;
	  Factor *= this->NbrLzValue;
	}
 
      for (int g = 0; g < this->NBodyValue; ++g)
	TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][shift + g] * this->PrecalculatedInteractionCoefficients[m2][shiftedMomentumTransfer*this->NBodyValue + g];
    }
 
 
 for (int k = 0 ; k < this->NBodyValue; ++k)
  mIndices2[k] = mIndices[k];
 while (std::prev_permutation(mIndices2, mIndices2 + this->NBodyValue))
  {
    m1 = mIndices2[0];
    Factor = this->NbrLzValue;
    for (int k = 1 ; k < this->NBodyValue; ++k)
      {
	m1 += mIndices2[k]*Factor;
	Factor *= this->NbrLzValue;
      }
    for (int g = 0; g < this->NBodyValue; ++g)
	TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][shift + g] * this->PrecalculatedInteractionCoefficients[m2][shiftedMomentumTransfer*this->NBodyValue + g];
    
    for (int k = 0 ; k < this->NBodyValue; ++k)
      nIndices2[k] = nIndices[k];
    while (std::prev_permutation(nIndices2, nIndices2 + this->NBodyValue))
    {
      
      m2 = nIndices2[0];
      Factor = this->NbrLzValue;
      for (int k = 1 ; k < this->NBodyValue; ++k)
	{
	  m2 += nIndices2[k]*Factor;
	  Factor *= this->NbrLzValue;
	}
 
      for (int g = 0; g < this->NBodyValue; ++g)
	TmpInteraction += this->PrecalculatedInteractionCoefficients[m1][shift + g] * this->PrecalculatedInteractionCoefficients[m2][shiftedMomentumTransfer*this->NBodyValue + g];
    }
  }    
  
  delete[] mIndices2;
  delete[] nIndices2;
  

  TmpInteraction /= FinalFactor;
  return TmpInteraction;
}


#endif
