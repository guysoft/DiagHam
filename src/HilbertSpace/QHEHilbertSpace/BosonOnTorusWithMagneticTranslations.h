////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on torus with magnetic translations            //
//                                                                            //
//                        last modification : 13/10/2003                      //
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


#ifndef BOSONONTORUSWITHMAGNETICTRANSLATIONS_H
#define BOSONONTORUSWITHMAGNETICTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/QHEHilbertSpace/ParticleOnTorusWithMagneticTranslations.h"

#include <iostream>


class BosonOnTorusState;


class BosonOnTorusWithMagneticTranslations :  public ParticleOnTorusWithMagneticTranslations
{

  // number of bosons
  int NbrBosons;
  // number of bosons plus 1
  int IncNbrBosons;
  // maximum momentum value reached by a boson
  int MaxMomentum;
  // number of different momentum value
  int NbrMomentum;
  // number of Lz values in a state
  int NbrLzValue;

  // momentum value in the x direction (modulo GCD of nbrBosons and maxMomentum)
  int XMomentum;
  // momentum value in the y direction (modulo GCD of nbrBosons and maxMomentum)
  int YMomentum;
  //  GCD of nbrBosons and maxMomentum
  int MomentumModulo;

  // array containing flag indicating if a state belonging to an orbit with a given number of member is compatible with x momentum constraint
  bool* CompatibilityXWithMomentum;
  // array containing flag indicating if a state belonging to an orbit with a given number of member is compatible with y momentum constraint
  bool** CompatibilityYWithMomentum;

  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;
  // number of state in each orbit
  int* NbrStateInOrbit;

  // array describing each state
  BosonOnTorusState* StateDescription;

  // BosonOnTorusState external parameter: reduced number of state (aka the number of unsigned long per state) minus 1
  int ReducedNbrState;
  // BosonOnTorusState external parameter: number of the state in the last unsigned long array describing the whole state
  int RemainderNbrState;

 public:

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // yMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)
  BosonOnTorusWithMagneticTranslations (int nbrBosons, int maxMomentum, int xMomentum, int yMomentum);

  // constructor from full datas
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
  // yMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)
  // hilbertSpaceDimension = Hilbert space dimension
  // stateDescription = array describing each state
  // stateMaxMomentum = array giving maximum Lz value reached for a boson in a given state
  BosonOnTorusWithMagneticTranslations (int nbrBosons, int maxMomentum, int xMomentum, int yMomentum, int hilbertSpaceDimension, 
					int** stateDescription, int* stateMaxMomentum);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnTorusWithMagneticTranslations(const BosonOnTorusWithMagneticTranslations& bosons);

  // destructor
  //
  ~BosonOnTorusWithMagneticTranslations ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnTorusWithMagneticTranslations& operator = (const BosonOnTorusWithMagneticTranslations& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  int GetParticleStatistic();

  // get momemtum value of a given state
  //
  // index = state index
  // return value = state momentum
  int GetMomentumValue(int index);

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
					 SubspaceSpaceConverter& converter);

  // apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of the destination state 
  int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

 private:

  // find state index
  //
  // stateDescription = array describing the state
  // lzmax = maximum Lz value reached by a boson in the state
  // return value = corresponding index
  int FindStateIndex(int* stateDescription, int lzmax);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson
  // return value = Hilbert space dimension
  int EvaluateHilbertSpaceDimension(int nbrBosons, int maxMomentum);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(int memory);

  // generate all states corresponding to the constraints
  // 
  // return value = hilbert space dimension
  int GenerateStates();

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson in the state
  // currentMaxMomentum = momentum maximum value for bosons that are still to be placed
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrBosons, int maxMomentum, int currentMaxMomentum, int pos);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // maxMomentum = momentum maximum value for a boson in the state
  // currentMaxMomentum = momentum maximum value for bosons that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentYMomentum = current value of the momentum in the y direction
  // return value = position from which new states have to be stored
  int RawGenerateStates(int nbrBosons, int maxMomentum, int currentMaxMomentum, int pos, int currentYMomentum);

};

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnTorusWithMagneticTranslations::GetParticleStatistic()
{
  return ParticleOnTorusWithMagneticTranslations::BosonicStatistic;
}

#endif


