////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//   class of fermion on a torus taking into account magnetic translations    //
//                                                                            //
//                        last modification : 10/09/2003                      //
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
#include "HilbertSpace/QHEHilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "Matrix/Matrix.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/ArrayTools.h"

#include <math.h>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


// basic constructor
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
// yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)

FermionOnTorusWithMagneticTranslations::FermionOnTorusWithMagneticTranslations (int nbrFermions, int maxMomentum, 
										int xMomentum, int yMomentum)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->MaxMomentum = maxMomentum;  

  this->NbrMomentum = this->MaxMomentum + 1;
  this->MomentumModulo = FindGCD(this->NbrFermions, this->MaxMomentum);
  this->XMomentum = xMomentum % this->MomentumModulo;
  this->YMomentum = yMomentum % this->MaxMomentum;

  this->StateShift = this->MaxMomentum / this->MomentumModulo;
  this->MomentumIncrement = (this->NbrFermions * this->StateShift) % this->MomentumModulo;
  this->ComplementaryStateShift = this->MaxMomentum - this->StateShift;
  this->MomentumMask = ((unsigned long) 1);
  for (int i = 1; i < this->StateShift; ++i)
    {
      this->MomentumMask <<= 1;
      this->MomentumMask |= ((unsigned long) 1);
    }

//  cout << this->StateShift << " " << this->MomentumIncrement << " " << this->ComplementaryStateShift << " " << hex << this->MomentumMask << dec << endl;

  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum);
  cout << this->HilbertSpaceDimension << endl;
  this->HilbertSpaceDimension = this->GenerateStates();
  cout << this->HilbertSpaceDimension << endl;

  this->Flag.Initialize();
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrMomentum * sizeof(int);
//  UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
  UsedMemory +=  (1 << MaximumSignLookUp) * sizeof(double);
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnTorusWithMagneticTranslations::FermionOnTorusWithMagneticTranslations(const FermionOnTorusWithMagneticTranslations& fermions)
{
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;
  this->MomentumModulo = fermions.MomentumModulo;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->StateShift = fermions.StateShift;
  this->MomentumIncrement = fermions.MomentumIncrement;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->Flag = fermions.Flag;
}

// constructor from full datas
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
// yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)
// hilbertSpaceDimension = Hilbert space dimension
// stateDescription = array describing each state

FermionOnTorusWithMagneticTranslations::FermionOnTorusWithMagneticTranslations (int nbrFermions, int maxMomentum, int xMomentum, 
										int yMomentum, int hilbertSpaceDimension, 
										unsigned long* stateDescription)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  this->MomentumModulo = FindGCD(this->NbrFermions, this->MaxMomentum);
  this->StateShift = this->MaxMomentum / this->MomentumModulo;
  this->MomentumIncrement = (this->NbrFermions * this->StateShift) % this->MomentumModulo;
  this->ComplementaryStateShift = this->NbrMomentum - this->StateShift;
  this->MomentumMask = ((unsigned long) 1);
  for (int i = 1; i < this->StateShift; ++i)
    {
      this->MomentumMask <<= 1;
      this->MomentumMask |= ((unsigned long) 1);
    }
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = stateDescription;
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrMomentum * sizeof(int);
  UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
  UsedMemory +=  (1 << MaximumSignLookUp) * sizeof(double);
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
}


// destructor
//

FermionOnTorusWithMagneticTranslations::~FermionOnTorusWithMagneticTranslations ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;
/*      delete[] this->SignLookUpTable;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;*/

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;

    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnTorusWithMagneticTranslations& FermionOnTorusWithMagneticTranslations::operator = (const FermionOnTorusWithMagneticTranslations& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;
      delete[] this->SignLookUpTable;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;
  this->StateShift = fermions.StateShift;
  this->MomentumIncrement = fermions.MomentumIncrement;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->Flag = fermions.Flag;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnTorusWithMagneticTranslations::Clone()
{
  return new FermionOnTorusWithMagneticTranslations(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnTorusWithMagneticTranslations::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->XMomentum, this->MomentumModulo);
  TmpList += new PeriodicMomentumQuantumNumber (this->YMomentum, this->MomentumModulo);
  List<AbstractQuantumNumber*> TmpList2;
  TmpList2 += new VectorQuantumNumber (TmpList);
  return TmpList2;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnTorusWithMagneticTranslations::GetQuantumNumber (int index)
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new PeriodicMomentumQuantumNumber (this->XMomentum, this->MomentumModulo);
  TmpList += new PeriodicMomentumQuantumNumber (this->YMomentum, this->MomentumModulo);
  return new VectorQuantumNumber (TmpList);
}

// get momemtum value in the y direction of a given state
//
// index = state index
// return value = state momentum in the y direction

int FermionOnTorusWithMagneticTranslations::GetYMomentumValue(int index)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned long State = this->StateDescription[index];
  int Momentum = 0;
  for (int i = 0; i <= StateMaxMomentum; ++i)
    {
      Momentum += ((State >> i ) & ((unsigned long) 1)) * i;
    }
  return (Momentum % this->MomentumModulo);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnTorusWithMagneticTranslations::ExtractSubspace (AbstractQuantumNumber& q, 
									       SubspaceSpaceConverter& converter)
{
  if (q.GetQuantumNumberType() != (AbstractQuantumNumber::Vector | AbstractQuantumNumber::PeriodicMomentum))
    return 0;
  if (((VectorQuantumNumber&) q).GetQuantumNumbers().GetNbrElement() != 2)
    return 0;
  if ((this->XMomentum == ((PeriodicMomentumQuantumNumber*) (((VectorQuantumNumber&) q)[0]))->GetMomentum()) &&
      (this->YMomentum == ((PeriodicMomentumQuantumNumber*) (((VectorQuantumNumber&) q)[1]))->GetMomentum()))
    return this;
  else 
    return 0;
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2[MaxMomentum])
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnTorusWithMagneticTranslations::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  int StateMaxMomentum = this->StateMaxMomentum[index];
  unsigned long State = this->StateDescription[index];
  if ((n1 > StateMaxMomentum) || (n2 > StateMaxMomentum) || ((State & (((unsigned long) 1) << n1)) == 0) || 
      ((State & (((unsigned long) 1) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewMaxMomentum = StateMaxMomentum;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n2);
  if (NewMaxMomentum == n2)
    while ((TmpState >> NewMaxMomentum) == 0)
      --NewMaxMomentum;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(((unsigned long) 0x1) << n1);
  if (NewMaxMomentum == n1)
    while ((TmpState >> NewMaxMomentum) == 0)
      --NewMaxMomentum;
  if ((TmpState & (((unsigned long) 1) << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewMaxMomentum)
    {
      NewMaxMomentum = m2;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m2);
  if ((TmpState & (((unsigned long) 1) << m1))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewMaxMomentum)
    {
      NewMaxMomentum = m1;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (((unsigned long) 0x1) << m1);
  cout << "before: " << hex << TmpState << dec << endl;
  TmpState = this->FindCanonicalForm(TmpState, NewMaxMomentum, nbrTranslation, this->YMomentum);
  if (this->TestXMomentumConstraint(TmpState, NewMaxMomentum) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, NewMaxMomentum);
//  cout << m1 << " " << m2 << " " << n1 << " " << n2 << endl;
//  cout << hex << this->StateDescription[index] << " " << TmpState << " " << dec << " " << TmpIndex << " " << nbrTranslation << " " << coefficient;
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
  coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & ((unsigned long) 0x1))));
  cout << "sign: " << (1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & ((unsigned long) 0x1))))) << " " << TmpIndex  << " " << this->ReorderingSign[TmpIndex] << " " << nbrTranslation << endl;
//  this->CurrentNbrStateInOrbit = this->NbrStateInOrbit[TmpIndex];
//  this->CurrentNbrStateInOrbitRatio = this->NbrStateInOrbit[index] / this->NbrStateInOrbit[TmpIndex];
//  this->CurrentSignature = this->StateSignature[TmpIndex];
//  this->CurrentSignature = (this->CurrentSignature >> this->SignatureComplementaryTransaltions[nbrTranslation]) | 
//    ((this->CurrentSignature & this->SignatureMasks[nbrTranslation]) << nbrTranslation);
  unsigned long TmpMask = 0;
  for (int i = 0; i < nbrTranslation; ++i)
    TmpMask |= ((unsigned long) 0x1) << i;
//  cout << hex << this->CurrentSignature << " ";
//  this->CurrentSignature = (this->CurrentSignature >> nbrTranslation) | ((this->CurrentSignature & TmpMask) << (this->MomentumModulo - nbrTranslation));
//  cout << this->CurrentSignature << " " << this->StateSignature[index];  
//  this->CurrentSignature ^= this->StateSignature[index];  
//  cout << this->CurrentSignature << dec << endl;
  nbrTranslation *= this->StateShift;
//  cout  << " " << coefficient << endl;
  return TmpIndex;
}

// apply a^+_m a_m operator to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for creation operator
// return value =  resulting multiplicative factor 

double FermionOnTorusWithMagneticTranslations::AdA (int index, int m)
{
  if (this->StateMaxMomentum[index] < m)
    return 0.0;
  if ((this->StateDescription[index] & (((unsigned long) 0x1) << m)) == 0)
    return 0.0;
  else
    return 1.0;
}

// find canonical form of a state description
//
// stateDescription = unsigned integer describing the state
// maxMomentum = reference on the maximum momentum value that can be reached by a fermion in the stateDescription state (will be changed to the one of the canonical form)
// nbrTranslation = number of translation needed to obtain the canonical form
// yMomentum = state momentum value in the y direction
// return value = canonical form of a state description

unsigned long FermionOnTorusWithMagneticTranslations::FindCanonicalForm(unsigned long stateDescription, int& maxMomentum, int& nbrTranslation, int yMomentum)
{
  nbrTranslation = 0;
  unsigned long CanonicalState = stateDescription;
  int index = 1;  
  while (index < this->MomentumModulo)
    {
      stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = index;
	}
      ++index;
    }
  if (nbrTranslation != 0)
    {
      maxMomentum = this->MaxMomentum;
      stateDescription = ((unsigned long) 1) << this->MaxMomentum;
      while ((CanonicalState & stateDescription) ==0)      
	{
	  --maxMomentum;
	  stateDescription >>= 1;
	}
    }
  return CanonicalState;
}

// find how many translations on the x direction are needed to obtain the same state
//
// stateDescription = unsigned integer describing the state
// return value = number of translation needed to obtain the same state

int  FermionOnTorusWithMagneticTranslations::FindNumberXTranslation(unsigned long stateDescription)
{
  unsigned long TmpState = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
  int index = 1;  
  while (TmpState != stateDescription)
    {
      TmpState = (TmpState >> this->StateShift) | ((TmpState & this->MomentumMask) << this->ComplementaryStateShift);
      ++index;
    }
  return index;
}

// test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescription = unsigned integer describing the state
// maxMomentum = maximum momentum value that can be reached by a fermion in the stateDescription state
// return value = true if the state satisfy the x momentum constraint

bool FermionOnTorusWithMagneticTranslations::TestXMomentumConstraint(unsigned long stateDescription, int maxMomentum)
{
  if (this->NbrFermions & 1)
    {
      if (this->XMomentum == 0)
	return true;
      unsigned long TmpState = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
      int index = 1;  
      while (TmpState != stateDescription)
	{
	  TmpState = (TmpState >> this->StateShift) | ((TmpState & this->MomentumMask) << this->ComplementaryStateShift);
	  ++index;
	}
      if (((this->XMomentum * index) % this->MomentumModulo) == 0)
	return true;
      else
	return false;
    }
  else
    {
      unsigned long TmpState2 = stateDescription & this->MomentumMask;
      unsigned long TmpState = (stateDescription >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
      int TmpSignature = 0;
      int TmpNbrParticle = 0;
      int index = 1;  
      for (int l = 0; l < this->StateShift; ++l)
	{
	  if (TmpState2 & 1)
	    ++TmpNbrParticle;
	  TmpState2 >>= 1;
	}
      if (TmpNbrParticle & 1)
	++TmpSignature;
      while (TmpState != stateDescription)
	{
	  TmpNbrParticle = 0;
	  TmpState2 = TmpState & this->MomentumMask;
	  TmpState = (TmpState >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
	  for (int l = 0; l < this->StateShift; ++l)
	    {
	      if (TmpState2 & 1)
		++TmpNbrParticle;
	      TmpState2 >>= 1;
	    }
	  if (TmpNbrParticle & 1)
	    ++TmpSignature;
	  ++index;
	}
      cout << TmpSignature << " " << index << " " << hex << stateDescription << dec << endl;
      if ((((this->XMomentum * index) - ((this->MomentumModulo * TmpSignature) >> 1)) % this->MomentumModulo) == 0)
	return true;
      else
	return false;
    }
}

// find state index
//
// stateDescription = unsigned longeger describing the state
// maxMomentum = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnTorusWithMagneticTranslations::FindStateIndex(unsigned long stateDescription, int maxMomentum)
{
  int Pos = 0;
//  int Pos = this->LookUpTable[maxMomentum][stateDescription >> this->LookUpTableShift[maxMomentum]];
//  cout << maxMomentum << " " << hex << stateDescription << dec << " " << this->LookUpTableShift[maxMomentum] << Pos << endl;
  while ((Pos < this->HilbertSpaceDimension) && (this->StateDescription[Pos] != stateDescription))
    ++Pos;
  if (Pos == this->HilbertSpaceDimension)
    {
      cout << "error : " << hex <<stateDescription << dec << " " << maxMomentum << endl;
    }
/*  if (this->StateMaxMomentum[Pos] != maxMomentum)
    cout << "error !!! " << maxMomentum << " " << this->StateMaxMomentum[Pos] << " " << hex << stateDescription << dec << " " << Pos << endl;*/
  return Pos;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnTorusWithMagneticTranslations::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  for (int i = 0; i < this->MaxMomentum; ++i)
    Str << ((TmpState >> i) & ((unsigned long) 0x1)) << " ";
  Str << "  (" << hex << this->ReorderingSign[state] << dec << ")";
//  Str << " key = " << this->Keys[state] << " maxMomentum position = " << this->MaxMomentumPosition[Max * (this->NbrFermions + 1) + TmpState[Max]]
//  int TmpMaxMomentum = this->StateMaxMomentum[state];
//  int TmpTranslation = 0;
/*  Str << " position = " << FindStateIndex(TmpState, this->StateMaxMomentum[state]) << " momentum = " << this->GetYMomentumValue(state) << " canonical form = ";
  TmpState = this->FindCanonicalForm(this->StateDescription[state], TmpMaxMomentum, TmpTranslation, this->YMomentum);
  for (int i = 0; i < this->MaxMomentum; ++i)
    Str << ((TmpState >> i) & 0x1) << " ";
  Str << " position = " << FindStateIndex(TmpState, TmpMaxMomentum) << " momentum = " << TmpMaxMomentum   << " translation = " << TmpTranslation ; */
  return Str;
}

// generate all states corresponding to the constraints
// 
// return value = hilbert space dimension

int FermionOnTorusWithMagneticTranslations::GenerateStates()
{
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateMaxMomentum = new int [this->HilbertSpaceDimension];
//  return this->RawGenerateStates(this->NbrFermions, this->MaxMomentum - 1, this->MaxMomentum - 1, 0, 0);
  this->HilbertSpaceDimension = this->RawGenerateStates(this->NbrFermions, this->MaxMomentum - 1, this->MaxMomentum - 1, 0, 0);
  int* TmpNbrStateDescription = new int [this->MaxMomentum + 1];  
  for (int i = 0; i <= this->MaxMomentum; ++i)
    {
      TmpNbrStateDescription[i] = 0;
    }
  int NbrTranslation;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      cout << hex << this->StateDescription[i] << dec << endl;
      this->StateDescription[i] = this->FindCanonicalForm(this->StateDescription[i], this->StateMaxMomentum[i], NbrTranslation, this->YMomentum);
      ++TmpNbrStateDescription[this->StateMaxMomentum[i]];
    }
  unsigned long** TmpStateDescription = new unsigned long* [this->MaxMomentum + 1];  
  bool** TmpStateDescriptionFlag = new bool* [this->MaxMomentum + 1];  
  for (int i = 0; i <= this->MaxMomentum; ++i)
    {
      if (TmpNbrStateDescription[i] != 0)
	{
	  TmpStateDescription[i] = new unsigned long [TmpNbrStateDescription[i]];
	  TmpStateDescriptionFlag[i] = new bool [TmpNbrStateDescription[i]];
	}
      else
	{
	  TmpStateDescription[i] = 0;
	  TmpStateDescriptionFlag[i] = 0;
	}
      TmpNbrStateDescription[i] = 0;
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      TmpStateDescription[this->StateMaxMomentum[i]][TmpNbrStateDescription[this->StateMaxMomentum[i]]++] = this->StateDescription[i];      
    }

  int TmpHilbertSpaceDimension = 0;
  for (int i = 0; i <= this->MaxMomentum; ++i) 
    {
      int CurrentNbrState = TmpNbrStateDescription[i];
      if (CurrentNbrState > 0)
	{
	  unsigned long* TmpStateArray = TmpStateDescription[i];
	  bool* TmpStateArrayFlag = TmpStateDescriptionFlag[i];
	  SortArrayUpOrdering(TmpStateArray, CurrentNbrState);
	  int TmpNbr = CurrentNbrState;
	  if (this->TestXMomentumConstraint(TmpStateArray[0], i) == true)
	    {
	      TmpStateArrayFlag[0] = true;
	    }
	  else
	    {
	      TmpStateArrayFlag[0] = false;
	      --TmpNbr;	      
	    }
	  for (int j = 1; j < CurrentNbrState; ++j)
	    {
	      while ((j < CurrentNbrState) && (TmpStateArray[j - 1] == TmpStateArray[j]))
		{
		  TmpStateArrayFlag[j] = false;
		  --TmpNbr;
		  ++j;
		}
	      if (j < CurrentNbrState)
		{
		  if (this->TestXMomentumConstraint(TmpStateArray[j], i) == true)
		    {
		      TmpStateArrayFlag[j] = true;
		    }
		  else
		    {
		      TmpStateArrayFlag[j] = false;
		      --TmpNbr;	      
		    }
		}
	    }
	  TmpHilbertSpaceDimension += TmpNbr;
	}
    }
  delete[] this->StateMaxMomentum;
  delete[] this->StateDescription;
  this->StateDescription = new unsigned long [TmpHilbertSpaceDimension];
  this->StateSignature = new int [TmpHilbertSpaceDimension];
  this->StateMaxMomentum = new int [TmpHilbertSpaceDimension];
  this->NbrStateInOrbit = new int [TmpHilbertSpaceDimension];
  this->ReorderingSign = new unsigned long [TmpHilbertSpaceDimension];
  int Pos = 0;
  unsigned long TmpState;
  unsigned long TmpState2;
  int TmpSignature;
  unsigned long TmpReorderingSign;
  int TmpNbrParticle;
  for (int i = 0; i <= this->MaxMomentum; ++i) 
    {
      int CurrentNbrState = TmpNbrStateDescription[i];
      if (CurrentNbrState > 0)
	{
	  unsigned long* TmpStateArray = TmpStateDescription[i];
	  bool* TmpStateArrayFlag = TmpStateDescriptionFlag[i];
	  for (int j = 0; j < CurrentNbrState; ++j)
	    {
	      if (TmpStateArrayFlag[j] == true)
		{
		  this->StateDescription[Pos] = TmpStateArray[j];
		  this->StateMaxMomentum[Pos] = i;
		  this->NbrStateInOrbit[Pos] = this->FindNumberXTranslation(this->StateDescription[Pos]);
		  TmpState = this->StateDescription[Pos];
		  TmpSignature = 0;
		  TmpReorderingSign = (unsigned long) 0;
		  if ((this->NbrFermions & 1) == 0)
		    {
		      cout << Pos << ":" << endl;
		      for (int k = 1; k <= this->NbrStateInOrbit[Pos]; ++k)
			{
			  TmpState2 = TmpState & this->MomentumMask;
			  TmpState =  (TmpState >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
			  cout << hex << TmpState << " " << TmpState2 << dec << endl;
			  TmpNbrParticle = 0;
			  for (int l = 0; l < this->StateShift; ++l)
			    {
			      if (TmpState2 & 1)
				++TmpNbrParticle;
			      TmpState2 >>= 1;
			    }
			  if (TmpNbrParticle & 1)
			    {			 
			      TmpReorderingSign |= ((TmpReorderingSign << 1) ^ (((unsigned long) 0x1) << k));
			      ++TmpSignature;
			    }
			  else
			    {
			      TmpReorderingSign |= ((TmpReorderingSign << 1) & (((unsigned long) 0x1) << k));
			    }			  
			}
		    }
		  this->ReorderingSign[Pos] = TmpReorderingSign;
		  this->StateSignature[Pos] = TmpSignature;
		  ++Pos;
		}
	    }
	  delete[] TmpStateArray;
	  delete[] TmpStateArrayFlag;
	}
    }
  delete[] TmpStateDescription;
  delete[] TmpStateDescriptionFlag;
  delete[] TmpNbrStateDescription;

  return TmpHilbertSpaceDimension;
}

// generate all states corresponding to the constraints (without taking into the canonical form) 
// 
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// currentYMomentum = current value of the momentum in the y direction
// return value = position from which new states have to be stored

int FermionOnTorusWithMagneticTranslations::RawGenerateStates(int nbrFermions, int maxMomentum, int currentMaxMomentum, int pos, int currentYMomentum)
{
  if ((nbrFermions == 0) || (nbrFermions > (currentMaxMomentum + 1)))
    return pos;
  if (nbrFermions == 1)
    {
/*      int i = this->YMomentum - (currentYMomentum % this->MomentumModulo);
      if (i < 0)
	i += this->MomentumModulo;
      for (; i <= currentMaxMomentum; i += this->MomentumModulo)*/
      int i = this->YMomentum - (currentYMomentum % this->MaxMomentum);
      if (i < 0)
	i += this->MaxMomentum;
      for (; i <= currentMaxMomentum; i += this->MaxMomentum)
	{
	  this->StateDescription[pos] = ((unsigned long) 1) << i;
	  this->StateMaxMomentum[pos] = maxMomentum;
	  ++pos;
	}
      return pos;
    }
  int ReducedCurrentMaxMomentum = currentMaxMomentum - 1;
  int TmpPos = this->RawGenerateStates(nbrFermions - 1, maxMomentum, ReducedCurrentMaxMomentum, pos, currentYMomentum + currentMaxMomentum);
  unsigned long Mask = ((unsigned long) 1) << currentMaxMomentum;
  for (int i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (maxMomentum == currentMaxMomentum)
    return this->RawGenerateStates(nbrFermions, ReducedCurrentMaxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentYMomentum);
  else
    return this->RawGenerateStates(nbrFermions, maxMomentum, ReducedCurrentMaxMomentum, TmpPos, currentYMomentum);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnTorusWithMagneticTranslations::GenerateLookUpTable(int memory)
{
//  this->ReminderTable = new unsigned long [1 << (this->StateShift + 1)];
  // evaluate look-up table size
  memory /= (4 * this->NbrMomentum);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrMomentum)
    this->MaximumLookUpShift = this->NbrMomentum;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [this->NbrMomentum];
  this->LookUpTableShift = new int [this->NbrMomentum];
  for (int i = 0; i < this->NbrMomentum; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize];
  int CurrentMaxMomentum = this->StateMaxMomentum[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
  if (CurrentMaxMomentum < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentMaxMomentum] = 0;
  else
    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
  unsigned long CurrentLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  unsigned long TmpLookUpTableValue;
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentMaxMomentum != this->StateMaxMomentum[i])
	{
 	  CurrentMaxMomentum = this->StateMaxMomentum[i];
	  TmpLookUpTable = this->LookUpTable[CurrentMaxMomentum];
	  if (CurrentMaxMomentum < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentMaxMomentum] = 0;
	  else
	    this->LookUpTableShift[CurrentMaxMomentum] = CurrentMaxMomentum + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentMaxMomentum];
	  CurrentLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      CurrentLookUpTableValue = TmpLookUpTableValue;
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  
  // look-up tables for evaluating sign when applying creation/annihilation operators
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double [Size];
  int Count;
  int TmpNbr;
  for (int j = 0; j < Size; ++j)
    {
      Count = 0;
      TmpNbr = j;
      while (TmpNbr != 0)
	{
	  if (TmpNbr & 0x1)
	    ++Count;
	  TmpNbr >>= 1;
	}
      if (Count & 1)
	this->SignLookUpTable[j] = -1.0;
      else
	this->SignLookUpTable[j] = 1.0;
    }
#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#endif

  this->RescalingFactors = new double* [this->NbrMomentum];
  for (int i = 1; i <= this->MaxMomentum; ++i)
    {
      this->RescalingFactors[i] = new double [this->NbrMomentum];
      for (int j = 1; j <= this->MaxMomentum; ++j)
	{
	  this->RescalingFactors[i][j] = sqrt (((double) i) / ((double) j));
	}
    }

  this->TranslationSign = new unsigned long [this->HilbertSpaceDimension];
  unsigned long TmpTranslationSign;
  unsigned long TmpState;
  unsigned long TmpReminder;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      TmpState = this->StateDescription[i];
      TmpTranslationSign = 0;
      Count = 0;	  
      for (int j = 1; j < this->MomentumModulo; ++j)
	{
	  TmpReminder = (TmpState & this->MomentumMask) << this->ComplementaryStateShift;      
	  TmpState = (TmpState >> this->StateShift) | TmpReminder;
	  for (int k = 0; k < this->StateShift; ++k)
	    if ((TmpReminder & (((unsigned long) 1) << k)) != 0)
	      ++Count;
	  if ((Count & 1) == 0)
	    {
	      TmpTranslationSign |= ((unsigned long) 1) << j;
	    }
	}
      this->TranslationSign[i] = TmpTranslationSign;
    }
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// maxMomentum = momentum maximum value for a fermion
// return value = Hilbert space dimension

int FermionOnTorusWithMagneticTranslations::EvaluateHilbertSpaceDimension(int nbrFermions, int maxMomentum)
{
  FactorialCoefficient Dimension; 
  Dimension.PartialFactorialMultiply(maxMomentum - nbrFermions + 1, maxMomentum); 
  Dimension.FactorialDivide(nbrFermions);
  return Dimension.GetIntegerValue();
}

