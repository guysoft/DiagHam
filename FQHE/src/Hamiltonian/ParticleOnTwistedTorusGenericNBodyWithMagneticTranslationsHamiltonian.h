////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a twisted torus      //
//                       with generic n-body interaction                      //
//                                                                            //
//                        last modification : 31/01/2015                      //
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


#ifndef PARTICLEONTWISTEDTORUSGENERICNBODYBODYWITHMAGNETICTRANSLATIONSHAMILTONIAN_H
#define PARTICLEONTWISTEDTORUSGENERICNBODYBODYWITHMAGNETICTRANSLATIONSHAMILTONIAN_H


#include "config.h"
#include "HilbertSpace/ParticleOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian.h"

#include <iostream>


using std::ostream;



class ParticleOnTwistedTorusGenericNBodyWithMagneticTranslationsHamiltonian : public AbstractQHEOnTorusWithMagneticTranslationsNBodyHamiltonian
{

 protected:

  // temporary arrays used during the interaction element evaluation
  double* QxValues;
  double* QyValues;
  double* Q2Values;
  double* CosineCoffients; 

  // angle (in radian) between the two fundamental cycles of the torus, along (L1 sin, L1 cos) and (0, L2)
  double Theta;
  // cosine of theta
  double CosTheta;
  // sine of theta
  double SinTheta;

 public:

  // default constructor
  //
  ParticleOnTwistedTorusGenericNBodyWithMagneticTranslationsHamiltonian();

  // constructor from default datas
  //
  // particles = Hilbert space associated to the system
  // nbrParticles = number of particles
  // maxMomentum = number of flux quanta
  // xMomentum = relative angular momentum along x 
  // ratio = torus aspect ratio (Lx/Ly)
  // theta =  angle (in pi units) between the two fundamental cycles of the torus, along (Lx sin theta, Lx cos theta) and (0, Ly)
  // nbrNBody = type of interaction i.e. the number of density operators that are involved in the interaction
  // architecture = architecture to use for precalculation
  // memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
  // onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
  // precalculationFileName = option file name where precalculation can be read instead of reevaluting them
  ParticleOnTwistedTorusGenericNBodyWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio, double theta,
									int nbrNBody, AbstractArchitecture* architecture, long memory = -1, bool onDiskCacheFlag = false, 
									char* precalculationFileName = 0);

  // destructor
  //
  ~ParticleOnTwistedTorusGenericNBodyWithMagneticTranslationsHamiltonian();
  
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
	
  // evaluate the numerical coefficient  in front of the Prod a^+_mi Prod a+_n coupling term
  //
  // mIndices = array containing the creation operator indices
  // nIndices = array containing the annihilation operator indices
  // return value = numerical coefficient  
  virtual Complex EvaluateInteractionCoefficient(int* mIndices, int* nIndices, long& nbrOperations);
  
  virtual Complex RecursiveEvaluateInteractionCoefficient(int xPosition, double currentSumQx, double currentSumQy, double currentSumQ2, double currentSumPhase, double& currentPrecision, long& nbrOperations);

  virtual double VFactor(double* q2Values);

};


#endif
