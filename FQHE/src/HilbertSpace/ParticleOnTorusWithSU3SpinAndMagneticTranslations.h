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


#ifndef PARTICLEONTORUSWITHSU3SPINANDMAGNETICTRANSLATIONS_H
#define PARTICLEONTORUSWITHSU3SPINANDMAGNETICTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractHilbertSpace.h"

#include <iostream>


using std::ostream;


class Matrix;


class ParticleOnTorusWithSU3SpinAndMagneticTranslations :  public AbstractHilbertSpace
{

 public:

  enum 
    {
      BosonicStatistic = 0x1,
      FermionicStatistic = 0x2,
    };
  
  // virtual destructor
  //
  virtual ~ParticleOnTorusWithSU3SpinAndMagneticTranslations ();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic() = 0;

  // apply a^+_m_1 a_m_1 operator to a given state (only state 1 Tz=+1/2, Y=+1/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_1 a_m_1
  virtual double Ad1A1 (int index, int m) = 0;

  // apply a^+_m_2 a_m_2 operator to a given state (only state 2 Tz=-1/2, Y=+1/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_2 a_m_2
  virtual double Ad2A2 (int index, int m) = 0;

  // apply a^+_m_3 a_m_3 operator to a given state (only state 3 Tz=0, Y=-2/3)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_3 a_m_3
  virtual double Ad3A3 (int index, int m) = 0;

  // apply a^+_m_1 a_n_1 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1A1 (int index, int m, int n, double& coefficient);

  // apply a^+_m_1 a_n_2 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1A2 (int index, int m, int n, double& coefficient);

  // apply a^+_m_1 a_n_3 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad1A3 (int index, int m, int n, double& coefficient);

  // apply a^+_m_2 a_n_1 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad2A1 (int index, int m, int n, double& coefficient);

  // apply a^+_m_2 a_n_2 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad2A2 (int index, int m, int n, double& coefficient);

  // apply a^+_m_2 a_n_3 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad2A3 (int index, int m, int n, double& coefficient);

  // apply a^+_m_3 a_n_1 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad3A1 (int index, int m, int n, double& coefficient);

  // apply a^+_m_3 a_n_2 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad3A2 (int index, int m, int n, double& coefficient);

  // apply a^+_m_3 a_n_3 operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int Ad3A3 (int index, int m, int n, double& coefficient);

  // apply a_n1_1 a_n2_1 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A1A1 (int index, int n1, int n2) = 0;

  // apply a_n1_1 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A1A2 (int index, int n1, int n2) = 0;

  // apply a_n1_1 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A1A3 (int index, int n1, int n2) = 0;

  // apply a_n1_2 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A2A2 (int index, int n1, int n2) = 0;

  // apply a_n1_2 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A2A3 (int index, int n1, int n2) = 0;

  // apply a_n1_3 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double A3A3 (int index, int n1, int n2) = 0;

  // apply a^+_m1_1 a^+_m2_1 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int Ad1Ad1 (int m1, int m2, double& coefficient, int& nbrTranslation) = 0;

  // apply a^+_m1_1 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int Ad1Ad2 (int m1, int m2, double& coefficient, int& nbrTranslation) = 0;

  // apply a^+_m1_1 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int Ad1Ad3 (int m1, int m2, double& coefficient, int& nbrTranslation) = 0;

  // apply a^+_m1_2 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int Ad2Ad2 (int m1, int m2, double& coefficient, int& nbrTranslation) = 0;

  // apply a^+_m1_2 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int Ad2Ad3 (int m1, int m2, double& coefficient, int& nbrTranslation) = 0;

  // apply a^+_m1_3 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = number of translations needed to find the canonical state
  // return value = index of the destination state 
  virtual int Ad3Ad3 (int m1, int m2, double& coefficient, int& nbrTranslation) = 0;

  
};


#endif
