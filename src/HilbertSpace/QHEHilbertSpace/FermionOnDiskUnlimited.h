////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of fermions on disk with no restriction on the         //
//               number of reachable states or the number of fermions         //
//                                                                            //
//                        last modification : 03/03/2004                      //
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


#ifndef FERMIONONDISKUNLIMITED_H
#define FERMIONONDISKUNLIMITED_H


#include "config.h"
#include "HilbertSpace/QHEHilbertSpace/ParticleOnDisk.h"
#include "HilbertSpace/QHEHilbertSpace/FermionOnSphereLongState.h"


class FermionOnDiskUnlimited:  public ParticleOnDisk
{

  // number of fermions
  int NbrFermions;
  // number of fermions plus 1
  int IncNbrFermions;
  // momentum total value
  int TotalLz;
  //  maximum momentum that can be reached by a fermion
  int LzMax;
  // number of one particle state with a different momentum
  int NbrLzValue;

  // array describing each state
  FermionOnSphereLongState* StateDescription;
  // array giving maximum Lz value reached for a fermion in a given state
  int* StateLzMax;
  // array giving the reduced number of state (aka the number of unsigned long per state) minus 1
  int* ReducedNbrState;

  // look_up table use to find states with a given lzmaxm and a given hash key
  int*** LookUpTable;
  // number of states with the same lzmax and a hash key
  int** NbrStateInLookUpTable;
  // mask corresponding to the number of bits to take into account in a hash key
  unsigned long HashKeyMask;

  // a table containing the parity of the number of one for each integer ranging from 0 to 2^MaximumSignLookUp - 1
  double* SignLookUpTable;
  // a table containing the mask on the bits to keep for each shift that is requested by sign evaluation
  unsigned long* SignLookUpTableMask;
  // number to evalute size of SignLookUpTable
  int MaximumSignLookUp;

  // a temporary state that is used during calculations
  FermionOnSphereLongState TemporaryState;
  // reduced number of state (aka the number of unsigned long per state) minus 1 for the temporary state
  int TemporaryStateReducedNbrState;

 public:

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = momentum total value
  FermionOnDiskUnlimited (int nbrFermions, int totalLz);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnDiskUnlimited(const FermionOnDiskUnlimited& fermions);

  // destructor
  //
  ~FermionOnDiskUnlimited ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnDiskUnlimited& operator = (const FermionOnDiskUnlimited& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  int GetParticleStatistic();

  // get the maximum angular momentum that can be reached by a particle 
  //
  // return value = maximum momentum
  int GetMaxLz();

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
  // return value = index of the destination state 
  int AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
  //
  // index = index of the state on which the operator has to be applied
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient);

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  double AdA (int index, int m);

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
  // stateDescription = reference on the state description
  // reducedNbrState =reduced number of state (aka the number of unsigned long per state) minus 1
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  int FindStateIndex(FermionOnSphereLongState& stateDescription, int reducedNbrState, int lzmax);

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  int EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz);

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  void GenerateLookUpTable(int memory);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int pos);

};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnDiskUnlimited::GetParticleStatistic()
{
  return ParticleOnDisk::FermionicStatistic;
}

// get the maximum angular momentum that can be reached by a particle 
//
// return value = maximum momentum

inline int FermionOnDiskUnlimited::GetMaxLz()
{
  return this->LzMax;
}

#endif


