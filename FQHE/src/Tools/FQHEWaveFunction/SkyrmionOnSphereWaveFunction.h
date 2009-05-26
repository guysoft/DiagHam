////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of SkyrmionOnSphere state wave function                    //
//                                                                            //
//                        last modification : 20/04/2005                      //
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


#ifndef SKYRMIONONSPHEREWAVEFUNCTION_H
#define SKYRMIONONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "QHEWaveFunctionManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "Vector/RealVector.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"

#ifdef USE_HILBERT_SPACE

class AbstractQHEParticle;
class AbstractFunctionBasis;


class SkyrmionOnSphereWaveFunction: public Abstract1DComplexFunction
{

 protected:

  // Pointer to architecture
  AbstractArchitecture* Architecture;

  // Number of particles;
  int NbrParticles;

  // Nbr of flux in polarized part
  int PolarizedLzMax;

  // Angular momentum in polarized part
  int PolarizedLz;

  // data for bosonic hilbert-space
  int BosonLzMax;
  int BosonLz;
  int BosonSz;
  
  // vector that describes polarized state components in PolarizedSpace basis
  RealVector PolarizedState;

  // vector that describes bosonic spin texture components in BosonicSpace basis
  RealVector BosonicState;

  // Hilbert space associated to the polarized SkyrmionOnSphere state
  ParticleOnSphere* PolarizedSpace;

  // Hilbert space associated to the polarized SkyrmionOnSphere state
  ParticleOnSphereWithSpin *BosonicSpace;

  // one body real space basis to use 
  AbstractFunctionBasis* OneBodyBasis;  
  
  // flag indicating whether we are using an exact polarized wavefunction
  bool UseExact;

  // Analytic polarized wavefunction
  Abstract1DComplexFunction* AnalyticPolarizedWaveFunction;

 public:

  // constructor
  //
  // create object to be initialized from system options
  //
  SkyrmionOnSphereWaveFunction(AbstractArchitecture* architecture, OptionManager &manager, int nbrParticles,
			       int totalLzMax, int totalLz, int totalSz, QHEWaveFunctionManager *wfManager,
			       int basisType=ParticleOnSphereFunctionBasis::LeftHanded);


  // copy constructor
  //
  // function = reference on the wave function to copy
  SkyrmionOnSphereWaveFunction(const SkyrmionOnSphereWaveFunction& function);

  // destructor
  //
   ~SkyrmionOnSphereWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

  // test parity under reversal of all spins for the two parts of the function
  void TestSymmetries(RealVector& x);

  // add an option group containing all options related to the skyrmion wave functions
  //
  // manager = pointer to the option manager
  static void AddSkyrmionOptionGroup(OptionManager &manager, QHEWaveFunctionManager *wfManager);

};


#endif

#endif

