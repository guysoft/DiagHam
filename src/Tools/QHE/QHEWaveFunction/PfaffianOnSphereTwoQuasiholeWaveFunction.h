////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of Pfaffian wave function with two quasiholes on sphere         //
//                                                                            //
//                        last modification : 18/07/2006                      //
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


#ifndef PFAFFIANONSPHERETWOQUASIHOLEWAVEFUNCTION_H
#define PFAFFIANONSPHERETWOQUASIHOLEWAVEFUNCTION_H


#include "config.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"


class PfaffianOnSphereTwoQuasiholeWaveFunction: public Abstract1DComplexFunction
{

 protected:

  // number of particles
  int NbrParticles;

  // position of the first quasihole (spherical coordinates)
  double Theta1;
  double Phi1;
  // position of the second quasihole (spherical coordinates)
  double Theta2;
  double Phi2;


 public:

  // constructor
  //
  // nbrParticles = number of particles
  // theta1 = position of the first quasihole (spherical coordinates, theta angle)
  // phi1 = position of the first quasihole (spherical coordinates, phi angle)
  // theta2 = position of the second quasihole (spherical coordinates, theta angle)
  // phi2 = position of the second quasihole (spherical coordinates, phi angle)
  PfaffianOnSphereTwoQuasiholeWaveFunction(int nbrParticles, double theta1, double phi1, double theta2, double phi2);

  // copy constructor
  //
  // function = reference on the wave function to copy
  PfaffianOnSphereTwoQuasiholeWaveFunction(const PfaffianOnSphereTwoQuasiholeWaveFunction& function);

  // destructor
  //
   ~PfaffianOnSphereTwoQuasiholeWaveFunction();

  // clone function 
  //
  // return value = clone of the function 
  Abstract1DComplexFunction* Clone ();

  // evaluate function at a given point
  //
  // x = point where the function has to be evaluated
  // return value = function value at x  
  Complex operator ()(RealVector& x);

};

#endif
