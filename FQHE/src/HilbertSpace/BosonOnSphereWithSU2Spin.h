////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//      altenate version of the class for bosons on sphere with SU(2) spin    //
//                                                                            //
//                        last modification : 26/01/2012                      //
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


#ifndef BOSONONSPHEREWITHSU2SPIN_H
#define BOSONONSPHEREWITHSU2SPIN_H


#include "config.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"

#include <iostream>

using std::cout;
using std::endl;


class BosonOnSphereShort;
class FermionOnSphereWithSpin;


class BosonOnSphereWithSU2Spin :  public ParticleOnSphereWithSpin
{


  friend class BosonOnSphereWithSU4Spin;


 protected:

  // number of bosons
  int NbrBosons;
  // number of bosons plus 1
  int IncNbrBosons;
  // momentum total value
  int TotalLz;
  // maximum Lz value reached by a boson
  int LzMax;
  // number of Lz values in a state
  int NbrLzValue;
  // twice the total spin value
  int TotalSpin;
  // number of bosons with spin up
  int NbrBosonsUp;
  // number of bosons with spin down
  int NbrBosonsDown;

  // lzmax value for the fermionic states associated to the type up particles
  int NUpLzMax;
  // lzmax value for the fermionic states associated to the type down particles
  int NDownLzMax;
  // largest lzmax value for the fermionic states among N1LzMax, N2LzMax and N3LzMax
  int FermionicLzMax;

  // array that contains the state description, the first entry being StateDescriptionUp and the second entry being StateDescriptionDown
  unsigned long* StateDescriptionSigma[2];

  // temporay array describing the type up particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescriptionUp;
  // temporay array describing the type down particle occupation number of each state, only used during the Hilbert space generation
  unsigned long* StateDescriptionDown;

  // sorted array that contains each unique configuration for the type up particles
  unsigned long* UniqueStateDescriptionUp;
  // number of time each unique configuration for the type up particles appears in StateDescriptionUp
  int* UniqueStateDescriptionSubArraySizeUp;
  // number of unique configurations for the type up-plus particles
  long NbrUniqueStateDescriptionUp;
  // first time a type up appears in the Hilbert space
  int* FirstIndexUniqueStateDescriptionUp;

  // temporary state used when applying operators for type up particles
  unsigned long* TemporaryStateUp;
  // temporary state used when applying operators for type down particles
  unsigned long* TemporaryStateDown;
  // array that contains the temporary state, the first entry being TemporaryStateUp and the second entry being TemporaryStateDown
  unsigned long* TemporaryStateSigma[2];

  // temporary state used when applying ProdA operator for type up particles
  unsigned long* ProdATemporaryStateUp;
  // temporary state used when applying ProdA operator for type down particles
  unsigned long* ProdATemporaryStateDown;
  // array that contains the temporary state used when applying ProdA operator, the first entry being ProdATemporaryStateUp and the second entry being ProdATemporaryStateDown
  unsigned long* ProdATemporaryStateSigma[2];

  // pointer to the Hilbert space where the result of any operator lies
  BosonOnSphereWithSU2Spin* TargetSpace;

 public:

  // default constructor
  // 
  BosonOnSphereWithSU2Spin ();

  // basic constructor without any constraint on Sz
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  BosonOnSphereWithSU2Spin (int nbrBosons, int totalLz, int lzMax);

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // totalLz = twice the momentum total value
  // lzMax = twice the maximum Lz value reached by a boson
  // totalSpin = twice the total spin value
  // memory = amount of memory granted for precalculations
  BosonOnSphereWithSU2Spin (int nbrBosons, int totalLz, int lzMax, int totalSpin, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnSphereWithSU2Spin(const BosonOnSphereWithSU2Spin& bosons);

  // destructor
  //
  virtual ~BosonOnSphereWithSU2Spin ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnSphereWithSU2Spin& operator = (const BosonOnSphereWithSU2Spin& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();

  // get the number of orbitals
  //
  // return value = number of orbitals
  virtual int GetNbrOrbitals();

  // get the number of particles
  //
  // return value = number of particles
  virtual int GetNbrParticles();

  // set a different target space (for all basic operations)
  //
  // targetSpace = pointer to the target space
  virtual void SetTargetSpace(ParticleOnSphereWithSpin* targetSpace);

  // return Hilbert space dimension of the target space
  //
  // return value = Hilbert space dimension
  virtual int GetTargetHilbertSpaceDimension();

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

  // apply a^+_m_d a_m_d operator to a given state (only spin down isospin plus)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_dm a_m_dm
  virtual double AddAd (int index, int m);

  // apply a^+_m_u a_m_u operator to a given state  (only spin up isospin plus)
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m_um a_m_um
  virtual double AduAu (int index, int m);

  // apply a^+_m_u a_n_u operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAu (int index, int m, int n, double& coefficient);

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

  // apply a^+_m_d a_n_d operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation operator
  // n = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAd (int index, int m, int n, double& coefficient);

  // apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // sigma1 = SU(2) index for the first annihilation operator
  // sigma2 = SU(2) index for the second annihilation operator
  // return value =  multiplicative factor 
  virtual double AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2);

  // apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AuAu (int index, int n1, int n2);

  // apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AuAd (int index, int n1, int n2);

  // apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
  //
  // index = index of the state on which the operator has to be applied
  // n1 = first index for annihilation operator
  // n2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AdAd (int index, int n1, int n2);

  // apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is is 0 for up and 1 for down
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // sigma1 = SU(2) index for the first creation operator
  // sigma2 = SU(2) index for the second creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient);

  // apply a^+_m1_u a^+_m2_u operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdu (int m1, int m2, double& coefficient);

  // apply a^+_m1_u a^+_m2_d operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AduAdd (int m1, int m2, double& coefficient);

  // apply a^+_m1_d a^+_m2_d operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AddAdd (int m1, int m2, double& coefficient);

  // apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
  //
  // index = index of the state on which the operator has to be applied
  // n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // return value =  multiplicative factor 
  virtual double ProdA (int index, int* n, int* spinIndices, int nbrIndices);

  // apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
  //
  // m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
  // spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
  // nbrIndices = number of creation (or annihilation) operators
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient);

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
  virtual ostream& PrintStateMonomial (ostream& Str, long state);
 
  // convert a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
  //
  // initialState = state to transform  
  // targetState = vector where the transformed state has to be stored
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // firstComponent = index of the first component to compute in initialState
  // nbrComponents = number of consecutive components to compute
  virtual void TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, long firstComponent = 0l, long nbrComponents = 0l);

  // compute the transformation matrix from one SU(2) basis to another, transforming the one body basis in each momentum sector
  //
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // return value = transformation matrix
  virtual ComplexMatrix TransformationMatrixOneBodyBasis(ComplexMatrix* oneBodyBasis);

  // compute the projection matrix from the SU(2) Hilbert space to an U(1) Hilbert space
  // 
  // targetSpace = pointer to the U(1) Hilbert space
  // type = type of particles that has to be kept (0 for type up, 1 for type up-minus, 2 for type down, 3 for type down-minus)
  // return value = projection matrix
  virtual ComplexMatrix TransformationMatrixSU2ToU1(BosonOnSphereShort* targetSpace, int type = 0);

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
  //  virtual RealVector& ConvertFromUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // normalize Jack with respect to cylinder basis
  //
  // state = reference to the Jack state to normalize
  // aspect = aspect ratio of cylinder
  // return value = normalized state
  virtual RealVector& NormalizeJackToCylinder(RealVector& state, double aspect);

  // Compute the product of a spinful fermionic state with a Van der Monde determinant 
  //
  // fermionicState = reference on the spinful fermionic state
  // outputVector = reference on the vector where the result will be stored
  // fermionicSpace = pointer on the Hilbert Space associated to the spinful fermionic state
  // minIndex = first component to compute
  // nbrComponents = number of components to compute
  virtual void SlaterTimeSpinfulFermionicState(RealVector& fermionicState, RealVector& outputVector, FermionOnSphereWithSpin* fermionicSpace, 
					       int minIndex, int nbrComponents);

  protected:

  // find state index
  //
  // stateDescriptionUp = unsigned integer describing the fermionic state for type up particles
  // stateDescriptionDown = unsigned integer describing the fermionic state for type down particles
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // nbrNUp = number of particles with quantum number up
  // nbrNDown = number of particles with quantum number down
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int nbrNUp, int nbrNDown);

  // evaluate Hilbert space dimension without the Sz constraint
  //
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // return value = Hilbert space dimension
  virtual long ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz);

  // generate look-up table associated to current Hilbert space
  // 
  // memeory = memory size that can be allocated for the look-up table
  virtual void GenerateLookUpTable(unsigned long memory);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // lzMaxUp = momentum maximum value for a boson in the state with up
  // lzMaxDown = momentum maximum value for a boson in the state with down
  // totalLz = momentum total value
  // nbrNUp = number of particles with up
  // nbrNDown = number of particles with down
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int lzMaxUp, int lzMaxDown, int totalLz, 
			      int nbrNUp, int nbrNDown, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // lzMax = momentum maximum value for a boson
  // totalLz = momentum total value
  // currentFermionicPositionUp = current fermionic position within the state description for the spin up
  // currentFermionicPositionDown = current fermionic position within the state description for the spin down
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int lzMax, int totalLz, 
			      int currentFermionicPositionUp, int currentFermionicPositionDown, long pos);

  // convert a bosonic state into its fermionic counterpart
  //
  // initialState = reference on the array where initial bosonic state is stored
  // return value = corresponding fermionic state
  virtual unsigned long BosonToFermion(unsigned long*& initialState);

  // convert a bosonic state into its fermionic counterpart
  //
  // initialStateUp = reference on the array where initial bosonic state for the type up particles is stored
  // initialStateDown = reference on the array where initial bosonic state for the type down particles is stored
  // finalStateUp = reference on the corresponding fermionic state for the type up particles
  // finalStateDown = reference on the corresponding fermionic state for the type down particles
  virtual void BosonToFermion(unsigned long*& initialStateUp, unsigned long*& initialStateDown, 
			      unsigned long& finalStateUp, unsigned long& finalStateDown);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialState = initial fermionic state
  // initialStateLzMax= maximum lz value reached by the fermionic state
  // finalState = reference on the array where the bosonic state for the type up particles has to be stored
  virtual void FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState);

  // convert a fermionic state into its bosonic  counterpart
  //
  // initialStateUp = initial fermionic state for the type up particles
  // initialStateDown = initial fermionic state for the type down particles
  // finalStateUp = reference on the array where the bosonic state for the type up particles has to be stored
  // finalStateDown = reference on the array where the bosonic state for the type down particles has to be stored
  virtual void FermionToBoson(unsigned long initialStateUp, unsigned long initialStateDown, 
			      unsigned long*& finalStateUp, unsigned long*& finalStateDown);

  // convert a bosonic state to its monomial representation
  //
  // initialStateUp = initial spin up bosonic state in its fermionic representation 
  // initialStateDown = initial spin down bosonic state in its fermionic representation
  // finalStateUp = reference on the array where the monomial representation for the spin up has to be stored
  // finalStateDown = reference on the array where the monomial representation for the spin up has to be stored
  virtual void ConvertToMonomial(unsigned long initialStateUp, unsigned long initialStateDown, unsigned long*& finalStateUp, unsigned long*& finalStateDown);

  // convert a bosonic state from its monomial representation for single componentw
  //
  // initialStateUp = array where the monomial representation for the up spin is stored
  // initialStateDown = array where the monomial representation for the down spin is stored
  // finalStateUp = reference on the bosonic up state in its fermionic representation
  // finalStateDown = reference on the bosonic down state in its fermionic representation
  virtual unsigned long ConvertFromMonomial(unsigned long* initialStateUp, unsigned long* initialStateDown, 
					    unsigned long& finalStateUp, unsigned long& finalStateDown);

  // apply generic a^+_m1_i a^+_m2_j operator to the state produced using A*A* method (without destroying it)
  //
  // m1 = first index for creation operator
  // m2 = second index for creation operator
  // temporaryStatei= reference on the temporary array for the type i particles
  // temporaryStatej= reference on the temporary array for the type j particles
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state
  virtual int AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej, double& coefficient);

  // find state index
  //
  // stateDescriptionUp = array describing the bosonic state for type up particles
  // stateDescriptionDown = array describing the bosonic state for type down particles
  // return value = corresponding index
  virtual int FindStateIndex(unsigned long*& stateDescriptionUp, unsigned long*& stateDescriptionDown);

  // recursive part of the convertion from a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
  //
  // targetState = vector where the transformed state has to be stored
  // coefficient = current coefficient to assign
  // position = current particle consider in the n-body state
  // momentumIndices = array that gives the momentum partition of the initial n-body state
  // initialSU2Indices = array that gives the spin dressing the initial n-body state
  // currentSU2Indices = array that gives the spin dressing the current transformed n-body state
  // oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
  // occupationCoefficient = invert of the coefficient that comes from the initial state occupation number 
  // occupationCoefficientArray = array that provides 1/2 ln (N!)
  virtual void TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
				      int position, int* momentumIndices, int* initialSU2Indices, int* currentSU2Indices, ComplexMatrix* oneBodyBasis, 
				      double occupationCoefficient, double* occupationCoefficientArray);

  // Compute the product of a spinful Slater determinant with a Van der Monde determinant
  //
  // slaterUp = monomial representation of the Slater spin up part
  // slaterDown = monomial representation of the Slater spin up part
  // finalStatesUp = reference on the array where the spin up fermionic representation of the produced states will be stored
  // finalStatesUp = reference on the array where the spin down fermionic representation of the produced states will be stored
  // weigth = reference on the array where the coefficient of the states product will be stored
  // return value = number of produced states
  virtual unsigned long VanDerMondeTimesSlater (unsigned long* slaterUp, unsigned long* slaterDown, unsigned long*& finalStatesUp, 
						unsigned long*& finalStatesDown, long*& weigth);

};

// get the number of orbitals
//
// return value = number of orbitals

inline int BosonOnSphereWithSU2Spin::GetNbrOrbitals()
{
  return this->NbrLzValue;
}

// get the number of particles
//
// return value = number of particles

inline int BosonOnSphereWithSU2Spin::GetNbrParticles()
{
  return this->NbrBosons;
}

// get the particle statistic 
//
// return value = particle statistic

inline int BosonOnSphereWithSU2Spin::GetParticleStatistic()
{
  return AbstractQHEParticle::BosonicStatistic;
}

// convert a bosonic state into its fermionic counterpart
//
// initialState = reference on the array where initial bosonic state is stored
// return value = corresponding fermionic state

inline unsigned long BosonOnSphereWithSU2Spin::BosonToFermion(unsigned long*& initialState)
{
  unsigned long FinalState = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      FinalState |= ((1ul << initialState[i]) - 1ul) << Shift;
      Shift += initialState[i];
      ++Shift;
    }
  return FinalState;
}

// convert a bosonic state into its fermionic counterpart
//
// initialStateUp = reference on the array where initial bosonic state for the type up particles is stored
// initialStateDown = reference on the array where initial bosonic state for the type down particles is stored
// finalStateUp = reference on the corresponding fermionic state for the type up particles
// finalStateDown = reference on the corresponding fermionic state for the type down particles

inline void BosonOnSphereWithSU2Spin::BosonToFermion(unsigned long*& initialStateUp, unsigned long*& initialStateDown,
						     unsigned long& finalStateUp, unsigned long& finalStateDown)
{
  finalStateUp = 0x0ul;
  unsigned long Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalStateUp |= ((1ul << initialStateUp[i]) - 1ul) << Shift;
      Shift += initialStateUp[i];
      ++Shift;
    }
  finalStateDown = 0x0ul;
  Shift = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      finalStateDown |= ((1ul << initialStateDown[i]) - 1ul) << Shift;
      Shift += initialStateDown[i];
      ++Shift;
    }
}

// convert a fermionic state into its bosonic  counterpart
//
// initialState = initial fermionic state
// initialStateLzMax= maximum lz value reached by the fermionic state
// finalState = reference on the array where the bosonic state for the type up particles has to be stored

inline void BosonOnSphereWithSU2Spin::FermionToBoson(unsigned long initialState, int initialStateLzMax, unsigned long*& finalState)
{
  int FinalStateLzMax = 0;
  while ((initialStateLzMax >= 0) && ((initialState >> initialStateLzMax) == 0x0ul))
    --initialStateLzMax;
  while (initialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialState - 1ul) ^ (~initialState);
      TmpState &= ~(TmpState >> 1);
#ifdef  __64_BITS__
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0f0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ffff0000ul) != 0) << 4;      
      TmpPower |= ((TmpState & 0xffffffff00000000ul) != 0) << 5;      
#else
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ul) != 0) << 4;      
#endif
      finalState[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialState >>= TmpPower;
      ++FinalStateLzMax;
      initialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalState[FinalStateLzMax] = 0x0ul;
}


// convert a fermionic state into its bosonic  counterpart
//
// initialStateUp = initial fermionic state for the type up particles
// initialStateDown = initial fermionic state for the type down particles
// finalStateUp = reference on the array where the bosonic state for the type up particles has to be stored
// finalStateDown = reference on the array where the bosonic state for the type down particles has to be stored

inline void BosonOnSphereWithSU2Spin::FermionToBoson(unsigned long initialStateUp, unsigned long initialStateDown,
						     unsigned long*& finalStateUp, unsigned long*& finalStateDown)
{
  int FinalStateLzMax = 0;
  int InitialStateLzMax = this->NUpLzMax;
  while ((InitialStateLzMax >= 0) && ((initialStateUp >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStateUp - 1ul) ^ (~initialStateUp);
      TmpState &= ~(TmpState >> 1);
#ifdef  __64_BITS__
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0f0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ffff0000ul) != 0) << 4;      
      TmpPower |= ((TmpState & 0xffffffff00000000ul) != 0) << 5;      
#else
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ul) != 0) << 4;      
#endif
      finalStateUp[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateUp >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalStateUp[FinalStateLzMax] = 0x0ul;

  FinalStateLzMax = 0;
  InitialStateLzMax = this->NDownLzMax;
  while ((InitialStateLzMax >= 0) && ((initialStateDown >> InitialStateLzMax) == 0x0ul))
    --InitialStateLzMax;
  while (InitialStateLzMax >= 0)
    {
      unsigned long TmpState = (~initialStateDown - 1ul) ^ (~initialStateDown);
      TmpState &= ~(TmpState >> 1);
#ifdef  __64_BITS__
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0f0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ffff0000ul) != 0) << 4;      
      TmpPower |= ((TmpState & 0xffffffff00000000ul) != 0) << 5;      
#else
      unsigned int TmpPower = ((TmpState & 0xaaaaaaaaul) != 0);
      TmpPower |= ((TmpState & 0xccccccccul) != 0) << 1;
      TmpPower |= ((TmpState & 0xf0f0f0f0ul) != 0) << 2;
      TmpPower |= ((TmpState & 0xff00ff00ul) != 0) << 3;      
      TmpPower |= ((TmpState & 0xffff0000ul) != 0) << 4;      
#endif
      finalStateDown[FinalStateLzMax] = (unsigned long) TmpPower;
      ++TmpPower;
      initialStateDown >>= TmpPower;
      ++FinalStateLzMax;
      InitialStateLzMax -= TmpPower;
    }
  for (; FinalStateLzMax <= this->LzMax; ++FinalStateLzMax)
    finalStateDown[FinalStateLzMax] = 0x0ul;
}

// convert a bosonic state to its monomial representation
//
// initialStateUp = initial spin up bosonic state in its fermionic representation 
// initialStateDown = initial spin down bosonic state in its fermionic representation
// finalStateUp = reference on the array where the monomial representation for the spin up has to be stored
// finalStateDown = reference on the array where the monomial representation for the spin up has to be stored

inline void BosonOnSphereWithSU2Spin::ConvertToMonomial(unsigned long initialStateUp, unsigned long initialStateDown,
							unsigned long*& finalStateUp, unsigned long*& finalStateDown) 
{
  int InitialStateLzMax = this->NUpLzMax;
  int Index = 0;
  int TmpLz = this->NUpLzMax - this->NbrBosonsUp + 1;
  while (InitialStateLzMax >= 0)
    {
      while ((InitialStateLzMax >= 0) && (((initialStateUp >> InitialStateLzMax) & 0x1ul) != 0x0ul))
	{
	  finalStateUp[Index++] = TmpLz;
	  --InitialStateLzMax;
	}
      while ((InitialStateLzMax >= 0) && (((initialStateUp >> InitialStateLzMax) & 0x1ul) == 0x0ul))
	{
	  --TmpLz;
	  --InitialStateLzMax;
	}
    }
  InitialStateLzMax = this->NDownLzMax;
  TmpLz = this->NDownLzMax - this->NbrBosonsDown + 1;
  Index = 0;
  while (InitialStateLzMax >= 0)
    {
      while ((InitialStateLzMax >= 0) && (((initialStateDown >> InitialStateLzMax) & 0x1ul) != 0x0ul))
	{
	  finalStateDown[Index++] = TmpLz;
	  --InitialStateLzMax;
	}
      while ((InitialStateLzMax >= 0) && (((initialStateDown >> InitialStateLzMax) & 0x1ul) == 0x0ul))
	{
	  --TmpLz;
	  --InitialStateLzMax;
	}
    }
}

// convert a bosonic state from its monomial representation for single componentw
//
// initialStateUp = array where the monomial representation for the up spin is stored
// initialStateDown = array where the monomial representation for the down spin is stored
// finalStateUp = reference on the bosonic up state in its fermionic representation
// finalStateDown = reference on the bosonic down state in its fermionic representation

inline unsigned long BosonOnSphereWithSU2Spin::ConvertFromMonomial(unsigned long* initialStateUp, unsigned long* initialStateDown, 
								   unsigned long& finalStateUp, unsigned long& finalStateDown)
{
  finalStateUp = 0x0ul;
  for (int i = 0; i < this->NbrBosonsUp; ++i)
    finalStateUp |= 0x1ul << (initialStateUp[i] + ((unsigned long) (this->NbrBosonsUp - i)) - 1ul);
  finalStateDown = 0x0ul;
  for (int i = 0; i < this->NbrBosonsDown; ++i)
    finalStateDown |= 0x1ul << (initialStateDown[i] + ((unsigned long) (this->NbrBosonsDown - i)) - 1ul);
}

// apply generic a^+_m1_i a^+_m2_j operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// temporaryStatei= reference on the temporary array for the type i particles
// temporaryStatej= reference on the temporary array for the type j particles
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int BosonOnSphereWithSU2Spin::AdiAdj (int m1, int m2, unsigned long*& temporaryStatei, unsigned long*& temporaryStatej,
					     double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      this->TemporaryStateUp[i] = this->ProdATemporaryStateUp[i];
      this->TemporaryStateDown[i] = this->ProdATemporaryStateDown[i];
    }
  ++temporaryStatej[m2];
  coefficient = temporaryStatej[m2];
  ++temporaryStatei[m1];
  coefficient *= temporaryStatei[m1];
  coefficient = sqrt(coefficient);
  return this->FindStateIndex(this->TemporaryStateUp, this->TemporaryStateDown);
}

// apply a_n1_sigma1 a_n2_sigma2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call. Sigma is 0 for up and 1 for down
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// sigma1 = SU(2) index for the first annihilation operator
// sigma2 = SU(2) index for the second annihilation operator
// return value =  multiplicative factor 

inline double BosonOnSphereWithSU2Spin::AsigmaAsigma (int index, int n1, int n2, int sigma1, int sigma2)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateSigma[sigma1][n1] == 0) || (this->ProdATemporaryStateSigma[sigma2][n2] == 0) || 
      ((n1 == n2) && (sigma1 == sigma2) && (this->ProdATemporaryStateSigma[sigma1][n1] == 1)))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryStateSigma[sigma2][n2];
  --this->ProdATemporaryStateSigma[sigma2][n2];
  Coefficient *= this->ProdATemporaryStateSigma[sigma1][n1];
  --this->ProdATemporaryStateSigma[sigma1][n1];
  return sqrt(Coefficient);
}

// apply a^+_m1_sigma1 a^+_m2_sigma2 operator to the state produced using A*A* method (without destroying it). Sigma is is 0 for up and 1 for down
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// sigma1 = SU(3) index for the first creation operator
// sigma2 = SU(3) index for the second creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

inline int BosonOnSphereWithSU2Spin::AdsigmaAdsigma (int m1, int m2, int sigma1, int sigma2, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      this->TemporaryStateUp[i] = this->ProdATemporaryStateUp[i];
      this->TemporaryStateDown[i] = this->ProdATemporaryStateDown[i];
    }
  ++this->TemporaryStateSigma[sigma2][m2];
  coefficient = this->TemporaryStateSigma[sigma2][m2];
  ++this->TemporaryStateSigma[sigma1][m1];
  coefficient *= this->TemporaryStateSigma[sigma1][m1];
  coefficient = sqrt(coefficient);
  return this->FindStateIndex(this->TemporaryStateUp, this->TemporaryStateDown);
}

// find state index
//
// stateDescriptionUp = array describing the bosonic state for type up particles
// stateDescriptionDown = array describing the bosonic state for type down particles
// return value = corresponding index

inline int BosonOnSphereWithSU2Spin::FindStateIndex(unsigned long*& stateDescriptionUp, unsigned long*& stateDescriptionDown)
{
  unsigned long Tmp1;
  unsigned long Tmp2;
  this->BosonToFermion(stateDescriptionUp, stateDescriptionDown, Tmp1, Tmp2);
  return this->FindStateIndex(Tmp1, Tmp2);
}

#endif


