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

#include <iostream>


using std::ostream;



class ParticleOnTorusNBodyHardCoreWithMagneticTranslationsHamiltonian : public AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian
{

 protected:

  // array where the interaction coefficients are stored
  double** PrecalculatedInteractionCoefficients;
  // number of entries in PrecalculatedInteractionCoefficients
  int NbrEntryPrecalculatedInteractionCoefficients;

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
  // nbrPermutations1 = number of permutations of the creation indexes
  // nbrPermutations2 = number of permutations of the annihilation indexes
  // return value = numerical coefficient  
  virtual double EvaluateInteractionCoefficient(double** creationCoefficients, double** annihilationCoefficients, int nbrPermutations1, int nbrPermutations2);
  
  // evaluate the numerical coefficient  in front of the \prod_i a+_mi \prod_j a_nj coupling term (factor corresponding to the creation operators only) for each integer modulo the NBodyValue
  //
  // mIndices = array that contains the creation indices
  // momentumTransfer = momentum transfer operated by the \prod_i a+_mi \prod_j a_nj operator, in units of the number of flux quanta
  // return value = array of numerical coefficients 
  double* EvaluateInteractionCoefficientCreation(int* mIndices, int momentumTransfer);
  
  
  // evaluate the two nested Gaussian sum for a three body interaction
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
  
  // Evaluate the number of permutations of a set of indices
  //
  // mIndices = array that contains the creation indices
  //return value = number of permutations
  int EvaluateNumberOfPermutations(int* mIndices);
  
  // evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
  //
  // m1 = first index
  // m2 = second index
  // m3 = third index
  // m4 = fourth index
  // return value = numerical coefficient
  virtual double EvaluateTwoBodyInteractionCoefficient(int m1, int m2, int m3, int m4);


};


#endif
