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
//                        last modification : 30/01/2004                      //
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
#include "HilbertSpace/QHEHilbertSpace/ParticleOnDisk.h"


// virtual destructor
//

ParticleOnDisk::~ParticleOnDisk ()
{
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int ParticleOnDisk::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  int TmpM[2];
  int TmpN[2];
  TmpM[0] = m1;
  TmpM[1] = m2;
  TmpN[0] = n1;
  TmpN[1] = n2;
  return ProdAdProdA(index, TmpM, TmpN, 2, coefficient);
}

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

int ParticleOnDisk::AdAdAdAAA (int index, int m1, int m2, int m3, int n1, int n2, int n3, double& coefficient)
{
  int TmpM[3];
  int TmpN[3];
  TmpM[0] = m1;
  TmpM[1] = m2;
  TmpM[2] = m3;
  TmpN[0] = n1;
  TmpN[1] = n2;
  TmpN[2] = n3;
  return ProdAdProdA(index, TmpM, TmpN, 3, coefficient);
}

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

int ParticleOnDisk::AdAdAdAdAAAA (int index, int m1, int m2, int m3, int m4, int n1, int n2, int n3, int n4, double& coefficient)
{
  int TmpM[4];
  int TmpN[4];
  TmpM[0] = m1;
  TmpM[1] = m2;
  TmpM[2] = m3;
  TmpM[3] = m4;
  TmpN[0] = n1;
  TmpN[1] = n2;
  TmpN[2] = n3;
  TmpN[3] = n4;
  return ProdAdProdA(index, TmpM, TmpN, 4, coefficient);
}
