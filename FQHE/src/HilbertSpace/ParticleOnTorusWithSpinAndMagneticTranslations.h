////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//   class of particle with spin on a torus taking into account magnetic translations   //
//                                                                            //
//                        last modification : 27/11/2007                      //
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


#ifndef PARTICLEONTORUSWITHSPINANDMAGNETICTRANSLATIONS_H
#define PARTICLEONTORUSWITHSPINANDMAGNETICTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"

#include <iostream>


using std::ostream;


class Matrix;


class ParticleOnTorusWithSpinAndMagneticTranslations :  public AbstractHilbertSpace
{

 public:

  enum 
    {
      BosonicStatistic = 0x1,
      FermionicStatistic = 0x2,
    };
  
  // virtual destructor
  //
  virtual ~ParticleOnTorusWithSpinAndMagneticTranslations ();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic() = 0;

  // apply a^+_(d,m1) a^+_(d,m2) a_(d,n1) a_(d,n2) operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation) = 0;

  // apply a^+_(u,m1) a^+_(u,m2) a_(u,n1) a_(u,n2) operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation) = 0;

  // apply a^+_(u,m1) a^+_(d,m2) a_(d,n1) a_(u,n2) operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to be applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  virtual int AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation) = 0;
  
};


#endif //PARTICLEONTORUSWITHSPINANDMAGNETICTRANSLATIONS
