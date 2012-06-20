////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of particle with SU3 spin on a torus                   //
//                taking into account magnetic translations                   //
//                                                                            //
//                        last modification : 18/06/2012                      //
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
#include "HilbertSpace/ParticleOnTorusWithSU3SpinAndMagneticTranslations.h"

#include <iostream>


using std::cout;
using std::endl;


// virtual destructor
//

ParticleOnTorusWithSU3SpinAndMagneticTranslations::~ParticleOnTorusWithSU3SpinAndMagneticTranslations ()
{
}

// apply a^+_m_1 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnTorusWithSU3SpinAndMagneticTranslations::Ad1A1 (int index, int m, int n, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_1 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnTorusWithSU3SpinAndMagneticTranslations::Ad1A2 (int index, int m, int n, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_1 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnTorusWithSU3SpinAndMagneticTranslations::Ad1A3 (int index, int m, int n, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_2 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnTorusWithSU3SpinAndMagneticTranslations::Ad2A1 (int index, int m, int n, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_2 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnTorusWithSU3SpinAndMagneticTranslations::Ad2A2 (int index, int m, int n, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_2 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnTorusWithSU3SpinAndMagneticTranslations::Ad2A3 (int index, int m, int n, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_3 a_n_1 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnTorusWithSU3SpinAndMagneticTranslations::Ad3A1 (int index, int m, int n, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_3 a_n_2 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnTorusWithSU3SpinAndMagneticTranslations::Ad3A2 (int index, int m, int n, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

// apply a^+_m_3 a_n_3 operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnTorusWithSU3SpinAndMagneticTranslations::Ad3A3 (int index, int m, int n, double& coefficient)
{
  return this->HilbertSpaceDimension;
}

