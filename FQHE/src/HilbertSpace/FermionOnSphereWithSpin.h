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


#ifndef FERMIONONSPHEREWITHSPIN_H
#define FERMIONONSPHEREWITHSPIN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"

#include <iostream>


class FermionOnSphere;


class FermionOnSphereWithSpin :  public ParticleOnSphereWithSpin
{

  friend class FermionOnSphereWithSpinSzProjection;
  friend class FermionOnSphereWithSpinLzSzSymmetry;
  friend class FermionOnSphereWithSpinSzSymmetry;
  friend class FermionOnSphereWithSpinLzSymmetry;
  friend class FermionOnSphereWithSU4Spin;
  friend class FermionOnSphereWithSpinAllSz;
  friend class FermionOnSphereWithSpinHaldaneBasis;
  friend class FermionOnSphereWithSpinAllSzLzSymmetry;


 protected:

  // number of fermions
  int NbrFermions;
  // number of fermions plus 1
  int IncNbrFermions;
  // momentum total value
  int TotalLz;
  // maximum Lz value reached by a fermion
  int LzMax;
  // number of fermions with spin up / down
  int NbrFermionsUp;
  int NbrFermionsDown;
  // number of Lz values in a stat
  int NbrLzValue;
  // twice the total spin value
  int TotalSpin;
  // highest bit in a given state description
  int HighestBit;

  // array describing each state
  unsigned long* StateDescription;
  // array giving maximum Lz value reached for a fermion in a given state
  int* StateHighestBit;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given lzmax sector
  unsigned long LookUpTableMemorySize;
  // shift used in each lzmax sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  int** LookUpTable;

  // a table containing ranging from 0 to 2^MaximumSignLookUp - 1
  double* SignLookUpTable;
  // a table containing the mask on the bits to keep for each shift that is requested by sign evaluation
  unsigned long* SignLookUpTableMask;
  // number to evalute size of SignLookUpTable
  int MaximumSignLookUp;

  // temporary state used when applying ProdA operator
  unsigned long ProdATemporaryState;
  // Lz maximum value associated to temporary state used when applying ProdA operator
  int ProdALzMax;

  // target space for operations leaving the Hilbert-space
  FermionOnSphereWithSpin *TargetSpace;


 public:

  // default constructor
  //
  FermionOnSphereWithSpin();

  // basic constructor
  // 
  // nbrFermions = number of fermions
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a fermion
  // totalSpin = twice the total spin value
  // memory = amount of memory granted for precalculations
  FermionOnSphereWithSpin (int nbrFermions, int totalLz, int lzMax, int totalSpin, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSphereWithSpin(const FermionOnSphereWithSpin& fermions);

  // destructor
  //
  ~FermionOnSphereWithSpin ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSphereWithSpin& operator = (const FermionOnSphereWithSpin& fermions);

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
  virtual List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  virtual AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  virtual AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter);

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  
  void SetTargetSpace(ParticleOnSphereWithSpin* targetSpace);

  
  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  
  int GetTargetHilbertSpaceDimension();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  bool WriteHilbertSpace (char* fileName);
  
  // apply creation operator to a word, using the conventions
  // for state-coding and quantum numbers of this space
  // state = word to be acted upon
  // m = Lz value of particle to be added
  // s = spin index of particle to be added (0=down, 1=up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  virtual unsigned long Ad (unsigned long state, int m, int s, double& coefficient);


  // apply a^+_m1_d a^+_m2_d a_n1_d a_n2_d operator to a given state (with m1+m2=n1+n2, only spin down)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // n1 = first index for annihilation operator (spin down)
  // n2 = second index for annihilation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_m1_u a^+_m2_u a_n1_u a_n2_u operator to a given state (with m1+m2=n1+n2, only spin up)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_d_m1 a^+_u_m2 a_d_n1 a_u_n2 operator to a given state (with m1+m2=n1+n2)
  //
  // index = index of the state on which the operator has to be applied
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient);

  // apply a^+_m_d a_m_d operator to a given state (only spin down)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AddAd (int index, int m);

  // apply a^+_m_u a_m_u operator to a given state  (only spin up)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AduAu (int index, int m);

  // apply a^+_m_u a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAu (int index, int m, int n, double& coefficient);

  // apply a^+_m_d a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAd (int index, int m, int n, double& coefficient);

  // apply a^+_m_u a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAd (int index, int m, int n, double& coefficient);

  // apply a^+_m_d a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAu (int index, int m, int n, double& coefficient);

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin up)
  // return value =  multiplicative factor 
  virtual double AuAu (int index, int n1, int n2);

  // apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin down)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double AdAd (int index, int n1, int n2);

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator (spin up)
  // n2 = second index for annihilation operator (spin down)
  // return value =  multiplicative factor 
  virtual double AuAd (int index, int n1, int n2);

  // apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin up)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdu (int m1, int m2, double& coefficient);

  // apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin down)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAdd (int m1, int m2, double& coefficient);

  // apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
  //
  // m1 = first index for creation operator (spin up)
  // m2 = second index for creation operator (spin down)
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdd (int m1, int m2, double& coefficient);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int* spinIndices, int nbrIndices);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // spinIndices = integer that gives the spin indices associated to each annihilation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int spinIndices, int nbrIndices);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // spinIndices = integer that gives the spin indices associated to each creation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int spinIndices, int nbrIndices, double& coefficient);

  // get the variance of the state
  // index = index of state to consider
  virtual int StateVariance (int index);

  // carefully test whether state is in Hilbert-space and find corresponding state index
  //
  // stateDescription = unsigned integer describing the state
  // highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
  // return value = corresponding index, or dimension of space, if not found
  virtual int CarefulFindStateIndex(unsigned long stateDescription, int highestBit);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // print a given State using the monomial notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintStateMonomial (ostream& Str, int state);

  // print a given State using the monomial notation, separating spin up from spin down
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintStateMonomialSeparatedSpin (ostream& Str, int state);

  // evaluate wave function in real space using a given basis and only for agiven range of components
  //
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = wave function evaluated at the given location
  virtual Complex EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
					int firstComponent, int nbrComponent);                                
  
  // initialize evaluation of wave function in real space using a given basis and only for a given range of components and
  //
  // timeCoherence = true if time coherence has to be used
  virtual void InitializeWaveFunctionEvaluation (bool timeCoherence = false);

  // create an SU(2) state from two U(1) state
  //
  // upState = vector describing the up spin part of the output state
  // upStateSpace = reference on the Hilbert space associated to the up spin part
  // downState = vector describing the down spin part of the output state
  // downStateSpace = reference on the Hilbert space associated to the down spin part  
  // return value = resluting SU(2) state
  virtual RealVector ForgeSU2FromU1(RealVector& upState, FermionOnSphere& upStateSpace, RealVector& downState, FermionOnSphere& downStateSpace);

  // create a U(1) state from an SU(2) state
  //
  // state = vector describing the SU(2) state
  // u1Space = reference on the Hilbert space associated to the U(1) state
  // return value = resulting U(1) state
  virtual RealVector ForgeU1FromSU2(RealVector& state, FermionOnSphere& u1Space);

  // convert a state such that its components are now expressed in the unnormalized basis
  //
  // state = reference to the state to convert
  // reference = set which component as to be normalized to 1
  // symmetryFactor = if true also remove the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertToUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);    

  // convert a state such that its components are now expressed in the normalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has been normalized to 1
  // symmetryFactor = if true also add the symmetry factors
  // return value = converted state
  virtual RealVector& ConvertFromUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // Evaluate the Density Matrix of the spin up fermions in a sector with a fixed lzUp 
  //
  // lzUp = twice total momentum of up fermions.
  // groundstate = reference on the total system groundstate
  // return value = density matrix of the subsystem of spins up fermions.
  virtual RealSymmetricMatrix EvaluatePartialDensityMatrixSpinSeparation (int lzUp, RealVector & groundstate);

  // compute the sign when moving all the up spin to the right from state index
  //
  // state = state whose spin splitting sign has to be computed
  // return = splitting sign
  virtual double GetSpinSeparationSignFromIndex(unsigned long index);

  // compute part of the Jack polynomial square normalization in a given range of indices
  //
  // state = reference on the unnormalized Jack polynomial
  // minIndex = first index to compute
  // nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
  // return value = quare normalization
  double JackSqrNormalization (RealVector& outputVector, long minIndex, long nbrComponents);

  // particle hole conjugate the spin down electrons, only (valid for N_phi=2N-1)
  // source: input state vector
  // target: output state vector
  virtual void ParticleHoleConjugateDownSpins(RealVector &source, RealVector &target);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescription, int lzmax);


  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // totalSpin = twce the total spin value
  // return value = Hilbert space dimension
  virtual int EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin);

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // totalLz = momentum total value
  // totalSpin = number of particles with spin up
  // return value = Hilbert space dimension      
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz, int totalSpin);

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion
  // currentLzMax = momentum maximum value for fermions that are still to be placed
  // totalLz = momentum total value
  // totalSz = spin total value
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  int OldGenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSz);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // lzMax = momentum maximum value for a fermion in the state
  // totalLz = momentum total value
  // totalSpin = number of particles with spin up
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int lzMax, int totalLz, int totalSpin, long pos);

  // compute sign
  //
  // signs = 
  // return value = sign value (+1.0 or -1.0)
  double ComputeSign(unsigned long signs);


  // compute the sign for permuting all electrons with spin up to the left of those with spin down
  // index = index of the state
  // indicesUp = location of the Fermions with spin up, counting as zero to LzMax
  // return value = sign value (+1.0 or -1.0)
  double GetStateSign(int index, int *IndicesDown);

  // compute the sign when moving all the up spin to the right
  //
  // state = state whose spin splitting sign has to be computed
  // return = splitting sign
  virtual double GetSpinSeparationSign(unsigned long state);

  // convert a state to its monomial representation
  //
  // initialState = initial bosonic state in its fermionic representation
  // finalStateUp = reference on the array where the monomial spin up representation has to be stored
  // finalStateDown = reference on the array where the monomial spin down representation has to be stored
  virtual void ConvertToMonomial(unsigned long initialState, int*& finalStateUp, int*& finalStateDown);

  // convert a  state from its monomial representation
  //
  // initialStateUp = array where the monomial spin up representation is stored
  // initialStateDown = array where the monomial spin down representation is stored
  // return value = state in its fermionic representation
  virtual unsigned long ConvertFromMonomial(int* initialStateUp, int* initialStateDown);


};

// get the particle statistic 
//
// return value = particle statistic

inline int FermionOnSphereWithSpin::GetParticleStatistic()
{
  return ParticleOnSphereWithSpin::FermionicStatistic;
}

// convert a state to its monomial representation
//
// initialState = initial bosonic state in its fermionic representation
// finalStateUp = reference on the array where the monomial spin up representation has to be stored
// finalStateDown = reference on the array where the monomial spin down representation has to be stored

inline void FermionOnSphereWithSpin::ConvertToMonomial(unsigned long initialState, int*& finalStateUp, int*& finalStateDown)
{
  int IndexUp = 0;
  int IndexDown = 0;
  unsigned long Tmp;
  for (int j = this->LzMax; j >= 0; --j)
    {
      Tmp = initialState >> (j * 2);
      if ((Tmp & 2ul) != 0ul)
	finalStateUp[IndexUp++] = j;
      if ((Tmp & 1ul) != 0ul)
	finalStateDown[IndexDown++] = j;      
    }
}

// convert a  state from its monomial representation
//
// initialStateUp = array where the monomial spin up representation is stored
// initialStateDown = array where the monomial spin down representation is stored
// return value = state in its fermionic representation

inline unsigned long FermionOnSphereWithSpin::ConvertFromMonomial(int* initialStateUp, int* initialStateDown)
{
  unsigned long TmpState = 0x0ul;  
  for (int j = 0; j < this->NbrFermionsUp; ++j)
    TmpState |= 0x2ul << (initialStateUp[j] << 1);
  for (int j = 0; j < this->NbrFermionsDown; ++j)
    TmpState |= 0x1ul << (initialStateDown[j] << 1);
  return TmpState;
}

// compute the sign when moving all the up spin to the right from state index
//
// state = state whose spin splitting sign has to be computed
// return = splitting sign

inline double FermionOnSphereWithSpin::GetSpinSeparationSignFromIndex(unsigned long index)
{
  return this->GetSpinSeparationSign(this->StateDescription[index]);
}

// compute the sign when moving all the up spin to the right
//
// state = state whose spin splitting sign has to be computed
// return = splitting sign

inline double FermionOnSphereWithSpin::GetSpinSeparationSign(unsigned long state)
{
  double Sign = 1.0;
  for (int j = 0; j <= this->LzMax; ++j)
    {
      unsigned long TmpState = state;
      if ((TmpState & (0x1ul << ((2 * j) + 1))) != 0x0ul)
	{
	  TmpState &= (0x1ul << ((2 * j) + 1)) - 0x1ul;
#ifdef __64_BITS__
	  TmpState &= 0x5555555555555555ul;	  
	  TmpState ^= TmpState >> 32;
#else
	  TmpState &= 0x55555555ul;
#endif	  
	  TmpState ^= TmpState >> 16;
	  TmpState ^= TmpState >> 8;
	  TmpState ^= TmpState >> 4;
	  TmpState ^= TmpState >> 2;
	  if ((TmpState & 0x1ul) != 0x0ul)
	    Sign *= -1.0;
	}
    }
  return Sign;
}

#endif


