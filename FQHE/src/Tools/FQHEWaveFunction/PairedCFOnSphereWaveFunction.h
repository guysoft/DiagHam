////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                  //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                      with filled (pseudo) Landau levels                    //
//                                                                            //
//                        last modification : 18/05/2007                      //
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


#ifndef PAIREDCFONSPHEREWAVEFUNCTION_H
#define PAIREDCFONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "JainCFOnSphereOrbitals.h"


class PairedCFOnSphereWaveFunction: public Abstract1DComplexFunction
{

 protected:

  JainCFOnSphereOrbitals *Orbitals;

  double MooreReadCoefficient;
  double* CFCoefficients;
  int NbrLandauLevels;
  int NbrParticles;
  int AbsEffectiveFlux;

  GarbageFlag Flag;
  
  // factor used to multiply each element of the Slater matrix
  double ElementNorm;
  
  ComplexSkewSymmetricMatrix *Slater;

  // single-particle Jastrow factors
  Complex *Ji;

  // precalculated sums over Orbitals
  Complex **gAlpha;
  
 public:

  // constants

  enum
    {
      DefaultConventions = false,
      AssumeOldConventions = true
    };

  // default constructor
  //
  PairedCFOnSphereWaveFunction();

  // constructor
  //
  // nbrParticles = number of particles
  // nbrLandauLevel = number of Landau levels filled with composite fermions
  // nbrEffectiveFlux = number of flux quanta of the magnetic monopole field experienced by CF's
  // MooreReadCoefficient = prefactor of singular 1/z term in pair-wave function
  // CFCoefficients = prefactors of CF orbitals in shells 0, 1, 2, ... , nbrLandauLevel-1
  // correctPrefactors = flag that enables the correction of prefactors to adopt the conventions of previous code
  // jastrowPower = power to which the Jastrow factor has to be raised
  PairedCFOnSphereWaveFunction(int nbrParticles, int nbrLandauLevels, int nbrEffectiveFlux, double MooreReadCoefficient,
			       double * CFCoefficients, bool correctPrefactors=false, int jastrowPower=2);

  // copy constructor
  //
  // function = reference on the wave function to copy
  PairedCFOnSphereWaveFunction(const PairedCFOnSphereWaveFunction& function);

  // destructor
  //
  ~PairedCFOnSphereWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  virtual Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  virtual Complex operator ()(RealVector& x);

  // normalize the wave-function to one for the given particle positions
  // x = point where the function has to be evaluated
  void AdaptNorm(RealVector& x);



};

#endif //PAIREDCFONSPHEREWAVEFUNCTION
