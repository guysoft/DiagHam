////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                           class of particle on disk                        //
//                                                                            //
//                        last modification : 30/01/2002                      //
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


#ifndef PARTICLEONDISK_H
#define PARTICLEONDISK_H


#include "config.h"
#include "Complex.h"
#include "HilbertSpace/QHEHilbertSpace/AbstractQHEParticle.h"


class ParticleOnDisk :  public AbstractQHEParticle
{

 public:

  // virtual destructor
  //
  virtual ~ParticleOnDisk ();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic() = 0;

  // get the maximum angular momentum that can be reached by a particle 
  //
  // return value = maximum momentum
  virtual int GetMaxLz() = 0;

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

  // apply a^+_m1 a^+_m2 a^+_m3 a_n1 a_n2 a_n3 operator to a given state (with m1+m2+m3=n1+n2+n3)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // m3 = third index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // n3 = third index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAdAAA (int index, int m1, int m2, int m3, int n1, int n2, int n3, double& coefficient);

  // apply a^+_m1 a^+_m2 a^+_m3 a^+_m4 a_n1 a_n2 a_n3 a_n4 operator to a given state (with m1+m2+m3+m4=n1+n2+n3+n4)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // m3 = third index for creation operator
  // m4 = fourth index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // n3 = third index for annihilation operator
  // n4 = fourth index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAdAdAAAA (int index, int m1, int m2, int m3, int m4, int n1, int n2, int n3, int n4, double& coefficient);

  // apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
  //
  // index = index of the state on which the operator has to be applied
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient) = 0;

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdA (int index, int m) = 0;

};

#endif


