////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                      with filled (pseudo) Landau levels                    //
//                                                                            //
//                        last modification : 16/09/2004                      //
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


#ifndef JAINCFFILLEDLEVELONSPHEREWAVEFUNCTION_H
#define JAINCFFILLEDLEVELONSPHEREWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"


class JainCFFilledLevelOnSphereWaveFunction: public Abstract1DComplexFunction
{

 protected:

  // number of particles
  int NbrParticles;

  // number of Landau levels filled with composite fermions
  int NbrLandauLevels;
  
  // power to which the Jastrow factor has to be raised
  int JastrowPower;

  // array containing prefactors of each projected monopole harmonic
  double** NormalizationPrefactors;

  // array containing constant factors that appears in the sum of projected monopole harmonic (except LLL)
  double** SumPrefactors;

 public:

  // constructor
  //
  // nbrParticles = number of particles
  // nbrLandauLevel = number of Landau levels filled with composite fermions
  // jastrowPower = power to which the Jastrow factor has to be raised
  JainCFFilledLevelOnSphereWaveFunction(int nbrParticles, int nbrLandauLevels, int jastrowPower);

  // copy constructor
  //
  // function = reference on the wave function to copy
  JainCFFilledLevelOnSphereWaveFunction(const JainCFFilledLevelOnSphereWaveFunction& function);

  // destructor
  //
   ~JainCFFilledLevelOnSphereWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

 private:

  // evaluate normalization factors of projected monopole harmonics
  //
  void EvaluateNormalizationPrefactors();

  // evaluate constant factors that appears in the sum of projected monopole harmonic (except LLL)
  //
  void EvaluateSumPrefactors();

 protected:

  // evaluate composite fermion monopole spherical harmonic 
  //
  // spinorUCoordinates = spinor u coordinates where the function has to be evalauted
  // spinorVCoordinates = spinor v coordinates where the function has to be evalauted
  // coordinate = index of the main coordinate (aka coordinate before project onto the lowest Landau level)
  // momentum = monopole spherical harmonic Lz momentum
  // landauLevel = index of the pseudo Landau level
  // return value = value of the monopole spherical harmonic at the givne point
  Complex EvaluateCFMonopoleHarmonic (Complex* spinorUCoordinates, Complex* spinorVCoordinates,
				      int coordinate, int momentum, int landauLevel);

};

#endif
