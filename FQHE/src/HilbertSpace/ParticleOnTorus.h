////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of particle on a torus                       //
//                                                                            //
//                        last modification : 18/07/2002                      //
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


#ifndef PARTICLEONTORUS_H
#define PARTICLEONTORUS_H


#include "config.h"
#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/ParticleOnSphere.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include <iostream>


using std::ostream;


class Matrix;


class ParticleOnTorus :  public ParticleOnSphere
{

 public:

  // virtual destructor
  //
  virtual ~ParticleOnTorus ();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic() = 0;

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient) = 0;

  // return matrix representation of the annihilation operator a_i
  //
  // i = operator index
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  virtual Matrix& A (int i, Matrix& M) = 0;

  // return matrix representation ofthw creation operator a^+_i
  //
  // i = operator index
  // M = matrix where representation has to be stored
  // return value = corresponding matrix
  virtual Matrix& Ad (int i, Matrix& M) = 0;

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Ky sector.
  // 
  // nbrBosonSector = number of particles that belong to the subsytem 
  // kySector = Ky sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kySector, RealVector& groundState);
  virtual HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kySector, ComplexVector& groundState);

  // apply a magnetic translation along x to a given state
  //
  // index = state index 
  // return value = translated state index
  virtual int ApplyXMagneticTranslation(int index);

};

#endif


