////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of spin 1/2 chain with translastion invariance          //
//                                                                            //
//                        last modification : 29/01/2002                      //
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


#ifndef SPIN1_2CHAINWITHTRANSLATIONS_H
#define SPIN1_2CHAINWITHTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"
#include "Matrix/HermitianMatrix.h"

#include <iostream>


using std::ostream;


class Spin1_2ChainWithTranslations : public AbstractSpinChainWithTranslations
{

 protected: 

  // momentum 
  int Momentum;
  // array containing falg indicating if a state beloging to an orbit with a given number of member is compatible with momentum constraint
  bool* CompatibilityWithMomentum;

  // array containing rescaling factors when passing from one orbit to another
  double** RescalingFactors;
  // number of state in each orbit
  int* NbrStateInOrbit;

  // flag to indicate if the total sz component is fixed
  bool FixedSpinProjectionFlag;
  //  total sz component (if fixed)
  int Sz;
  
  // shift to apply to move a group spin using an elementary translation
  int StateShift;
  // shift to apply to move the spin from one end to the other one
  int ComplementaryStateShift;
  // mask to apply to get a group of spins
  unsigned long StateMask;

  // shift to apply to a state to obtain an index to the look-up table 
  int LookUpTableShift;
  // look-up table (LookUpTable[i] gives the index of the smallest state that greater than i <<  LookUpTableShift)
  long* LookUpTable;

  // array describing each state
  unsigned long* StateDescription;

 public:

  // default constructor
  //
  Spin1_2ChainWithTranslations ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // translationStep = indicates the step for an elementary translation
  // memorySize = memory size in bytes allowed for look-up table
  // memorySlice = maximum amount of memory that can be allocated to partially evalauted the states
  Spin1_2ChainWithTranslations (int chainLength, int momentum, int translationStep, int memorySize, int memorySlice);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // translationStep = indicates the step for an elementary translation
  // sz = twice the value of total Sz component
  // memorySize = memory size in bytes allowed for look-up table
  Spin1_2ChainWithTranslations (int chainLength, int momentum, int translationStep, int sz, int memorySize, int memorySlice);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin1_2ChainWithTranslations (const Spin1_2ChainWithTranslations& chain);

  // destructor
  //
  ~Spin1_2ChainWithTranslations ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin1_2ChainWithTranslations& operator = (const Spin1_2ChainWithTranslations& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // return Hilbert space dimension
  //
  // return value = Hilbert space dimension
  int GetHilbertSpaceDimension();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  AbstractQuantumNumber* GetQuantumNumber (int index);

  // return value of the value of the sum of the square of spin projection on (Oz) 
  //
  // index = index of the state to test
  // return value = twice spin projection on (Oz)
  double TotalSzSz (int index);

  // get the momentum of each state in the current Hilbert space
  //
  // return value = momentum value
  int GetMomentum();

  // return value of spin projection on (Oz) for a given state
  //
  // index = index of the state to test
  // return value = spin projection on (Oz)
  int TotalSz (int index);

  // return eigenvalue of Sz_i Sz_j associated to a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to consider
  // return value = corresponding eigenvalue
  double SziSzj (int i, int j, int state);

  // compute the parity (prod_i Sz_i) for a given state
  //
  // state = index of the state to be applied on Sz_i operator
  // return value = total Sz value
  virtual unsigned long Parity (int state);

  // return index of resulting state from application of P_ij operator on a given state
  //
  // i = first position
  // j = second position
  // state = index of the state to be applied on P_ij operator
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int Pij (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i S+_j operator on a given state
  //
  // i = position of S- operator
  // j = position of S+ operator
  // state = index of the state to be applied on S-_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SmiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S+_i operator on a given state (only valid if there is no constraint on total Sz)
  //
  // i = operator position
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int Spi (int i, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i operator on a given state (only valid if there is no constraint on total Sz)
  //
  // i = operator position
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int Smi (int i, int state, double& coefficient, int& nbrTranslation);
    
  // return index of resulting state from application of S+_i S+_j operator on a given state
  //
  // i = position of first S+ operator
  // j = position of second S+ operator
  // state = index of the state to be applied on S+_i S+_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SpiSpj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i S-_j operator on a given state
  //
  // i = position of first S- operator
  // j = position of second S- operator
  // state = index of the state to be applied on S-_i S-_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SmiSmj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S+_i S+_i operator on a given state
  //
  // i = position of first S+ operator
  // state = index of the state to be applied on S+_i S+_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SpiSpi (int i, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i S-_i operator on a given state
  //
  // i = position of the S- operator
  // state = index of the state to be applied on S-_i S-_i operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SmiSmi (int i, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S+_i Sz_j operator on a given state
  //
  // i = position of S+ operator
  // j = position of Sz operator
  // state = index of the state to be applied on S+_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SpiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // return index of resulting state from application of S-_i Sz_j operator on a given state
  //
  // i = position of S- operator
  // j = position of Sz operator
  // state = index of the state to be applied on S-_i Sz_j operator
  // coefficient = reference on double where numerical coefficient has to be stored
  // nbrTranslations = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
  // return value = index of resulting state
  int SmiSzj (int i, int j, int state, double& coefficient, int& nbrTranslation);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, SubspaceSpaceConverter& converter);

  // find state index
  //
  // state = state description
  // return value = corresponding index
  int FindStateIndex(unsigned long state);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  ostream& PrintState (ostream& Str, int state);

  // evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition.
  // 
  // nbrSpinUp = number of spin up that belong to the subsytem 
  // kSector = momentum of the subsystem
  // groundState = reference on the total system ground state
  // return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)
  HermitianMatrix EvaluatePartialDensityMatrixParticlePartition (int nbrSpinUpSector, int kSector, RealVector& groundState);

 protected:

  // return value of twice spin projection on (Oz) for a given state
  //
  // stateDescription = state to which the spin projection has to be evaluated
  // return value = twice spin projection on (Oz)
  int GetTotalSz (unsigned long stateDescription);

  // find the canonical form of a state
  //
  // state = state description
  // nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
  // return value = canonical form of the state
  unsigned long FindCanonicalForm(unsigned long state, int& nbrTranslation);

  // find the canonical form of a state and find how many translations are needed to obtain the same state
  //
  // stateDescription = state description
  // nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
  // nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state
  // return value = canonical form of the state
  unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslation, int& nbrTranslationToIdentity);

  // find how many translations are needed to obtain the same state
  //
  // stateDescription = unsigned integer describing the state
  // return value = number of translation needed to obtain the same state
  int FindNumberTranslation(unsigned long stateDescription);

  // constructor from pre-constructed datas
  //
  // hilbertSpaceDimension = Hilbert space dimension
  // chainDescription = array describing states
  // chainLength = number of spin 1
  // momemtum = total momentum of each state
  // sz = twice the value of total Sz component
  // fixedQuantumNumberFlag = true if hilbert space is restricted to a given quantum number
  // lookUpTableShift = shift to apply to a state to obtain an index to the look-up table 
  // complementaryStateShift = shift to apply to move the spin from one end to the other one
  Spin1_2ChainWithTranslations (int hilbertSpaceDimension, unsigned long* chainDescription, int chainLength, 
				int momentum, int sz, bool fixedQuantumNumberFlag, int lookUpTableShift, 
				int complementaryStateShift);

  // generate all states with constraint on total Sz
  //
  // position = current position in the Hilbert space basis
  // sitePosition = largest sit position that has to be filled 
  // currentNbrSpinUp = number of spin to put in each state
  // return value = new current position in the Hilbert space basis
  long GenerateStates(long position, int sitePosition, int currentNbrSpinUp);

  // create precalculation tables
  //
  void CreatePrecalculationTable();

  // create look-up table used to speed up index search
  //
  void CreateLookUpTable();

  // evaluate Hilbert space dimension
  //
  // nbrSpins = number of spins
  // sz = twice the z projection of the total momentum
  // return value = Hilbert space dimension
  long EvaluateHilbertSpaceDimension(int nbrSpins, int szMax);

};

// get the momentum of each state in the current Hilbert space
//
// return value = momentum value

inline int Spin1_2ChainWithTranslations::GetMomentum()
{
  return this->Momentum;
}

// find the canonical form of a state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// return value = canonical form of the state

inline unsigned long Spin1_2ChainWithTranslations::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState = stateDescription;
  int index = 1;  
  while (index < this->ChainLength)
    {
      stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->StateMask) << this->ComplementaryStateShift);
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = index;
	}
      ++index;
    }
  return CanonicalState;
}

// find the canonical form of a state and find how many translations are needed to obtain the same state
//
// stateDescription = state description
// nbrTranslation = reference on a integer where the number of translations needed to obtain the canonical form  will be stored
// nbrTranslationToIdentity = reference on the number of translation needed to obtain the same state
// return value = canonical form of the state

inline unsigned long Spin1_2ChainWithTranslations::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslation, int& nbrTranslationToIdentity)
{
  nbrTranslation = 0;
  nbrTranslationToIdentity = 1;
  unsigned long CanonicalState = stateDescription;
  unsigned long ReferenceState = stateDescription;
  stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->StateMask) << this->ComplementaryStateShift);
  while ((ReferenceState != stateDescription) && (nbrTranslationToIdentity < this->ChainLength))
    {
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = nbrTranslationToIdentity;
	}
      stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->StateMask) << this->ComplementaryStateShift);
      ++nbrTranslationToIdentity;
    }
  return CanonicalState;
}

// find how many translations are needed to obtain the same state
//
// stateDescription = unsigned integer describing the state
// return value = number of translation needed to obtain the same state

inline int Spin1_2ChainWithTranslations::FindNumberTranslation(unsigned long stateDescription)
{
  unsigned long TmpState = (stateDescription >> this->StateShift) | ((stateDescription & this->StateMask) << this->ComplementaryStateShift);
  int index = 1;  
  while (TmpState != stateDescription)
    {
      TmpState = (TmpState >> this->StateShift) | ((TmpState & this->StateMask) << this->ComplementaryStateShift);
      ++index;
    }
  return index;
}

#endif


