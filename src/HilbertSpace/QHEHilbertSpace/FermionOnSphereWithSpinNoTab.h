////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//          Copyright (C) 2001-2005 Gunnar Moller and Nicolas Regnault        //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere with spin without            //
//                            sign precalculation table                       //
//                                                                            //
//                        last modification : 12/12/2005                      //
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


#ifndef FERMIONONSPHEREWITHSPINNOTAB_H
#define FERMIONONSPHEREWITHSPINNOTAB_H


#include "config.h"
#include "HilbertSpace/QHEHilbertSpace/ParticleOnSphereWithSpin.h"

#include <iostream>


class FermionOnSphereWithSpinNoTab :  public ParticleOnSphereWithSpin
{

  // number of fermions
  int NbrFermions;
  // number of fermions plus 1
  int IncNbrFermions;
  // momentum total value
  int TotalLz;
  // maximum Lz value reached by a fermion
  int LzMax;
  // number of Lz values in a state
  int NbrLzValue;
  // twice the total spin value
  int TotalSpin;

  // array describing each state
  unsigned long* StateDescription;
  // array giving maximum Lz value reached for a fermion in a given state
  int* StateLzMax;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given lzmax sector
  unsigned long LookUpTableMemorySize;
  // shift used in each lzmax sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  int** LookUpTable;


 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSpin = twce the total spin value
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpinNoTab (int nbrFermions, int totalLz, int lzMax, int totalSpin, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpinNoTab(const FermionOnSphereWithSpinNoTab& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpinNoTab ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpinNoTab& operator = (const FermionOnSphereWithSpinNoTab& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  int GetParticleStatistic();

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

  // apply a^+_m1_d a^+_m2_d a_n1_d a_n2_d operator to a given state (with m1+m2=n1+n2, only spin down)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // n1 = first index for annihilation operator (spin down)
  // n2 = second index for annihilation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_m1_u a^+_m2_u a_n1_u a_n2_u operator to a given state (with m1+m2=n1+n2, only spin up)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_m_d a_m_d operator to a given state (only spin down)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  double AddAd (int index, int m);

  // apply a^+_m_u a_m_u operator to a given state  (only spin up)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  double AduAu (int index, int m);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  // evaluate wave function in real space using a given basis and only for agiven range of components
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
				int firstComponent, int nbrComponent);                                
  
  // initialize evaluation of wave function in real space using a given basis and only for a given range of components and
  //
  // timeCoherence = true if time coherence has to be used
  void InitializeWaveFunctionEvaluation (bool timeCoherence = false);
  
  private:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  int FindStateIndex(unsigned long stateDescription, int lzmax);

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // totalSpin = twce the total spin value
  // return value = Hilbert space dimension
  int EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin);

  // evaluate Hilbert space dimension with shifted values for lzMax and totalLz
  //
  // nbrFermions = number of fermions
  // lzMax = two times momentum maximum value for a fermion plus one 
  // totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
  // return value = Hilbert space dimension  
  int ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // totalSz = spin total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int totalSz, int pos);

  // compute sign
  //
  // signs = 
  // return value = sign value (+1.0 or -1.0)
  double ComputeSign(unsigned long signs);

};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnSphereWithSpinNoTab::GetParticleStatistic()
{
  return ParticleOnSphereWithSpin::FermionicStatistic;
}

#endif


