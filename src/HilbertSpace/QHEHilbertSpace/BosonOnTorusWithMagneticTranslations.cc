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


#include "config.h"
#include "HilbertSpace/QHEHilbertSpace/BosonOnTorusWithMagneticTranslations.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "MathTools/IntegerAlgebraTools.h"

#include <math.h>


using std::cout;
using std::endl;


// basic constructor
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// xMomentum = momentum in the x direction (modulo GCD of nbrBosons and maxMomentum)
// yMomentum = momentum in the y direction (modulo GCD of nbrBosons and maxMomentum)

BosonOnTorusWithMagneticTranslations::BosonOnTorusWithMagneticTranslations (int nbrBosons, int maxMomentum, 
									    int xMomentum, int yMomentum)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;

  this->MaxMomentum = maxMomentum;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  this->MomentumModulo = FindGCD(this->NbrBosons, this->MaxMomentum);

  this->ReducedNbrState = GetReducedNbrState (this->NbrMomentum);
  this->RemainderNbrState = GetRemainderNbrState (this->NbrMomentum);
  this->TemporaryState.Resize(this->ReducedNbrState);

  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->MaxMomentum);
  cout << this->HilbertSpaceDimension << endl;
  this->HilbertSpaceDimension = this->GenerateStates();
  cout << this->HilbertSpaceDimension << endl;

  this->Flag.Initialize();
//  this->GenerateLookUpTable(0);

#ifdef __DEBUG__
/*  int UsedMemory = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    UsedMemory += (this->StateMaxMomentum[i] + 1) * sizeof(int) + sizeof(int*) + sizeof(int);
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;*/
#endif
}

// constructor from full datas (with no constraint on the total momentum)
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// hilbertSpaceDimension = Hilbert space dimension
// stateDescription = array describing each state
// stateMaxMomentum = array giving maximum Lz value reached for a fermion in a given state

BosonOnTorusWithMagneticTranslations::BosonOnTorusWithMagneticTranslations (int nbrBosons, int maxMomentum, int xMomentum, int yMomentum, int hilbertSpaceDimension, 
									    int** stateDescription, int* stateMaxMomentum)
{
/*  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->MaxMomentum = maxMomentum;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MomentumConstraintFlag = false;
  this->HilbertSpaceDimension = hilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = stateDescription;
  this->StateMaxMomentum = stateMaxMomentum;
  this->GenerateLookUpTable(1000000);
#ifdef __DEBUG__
  int UsedMemory = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    UsedMemory += (this->StateMaxMomentum[i] + 1) * sizeof(int) + sizeof(int*) + sizeof(int);
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
#endif
*/}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnTorusWithMagneticTranslations::BosonOnTorusWithMagneticTranslations(const BosonOnTorusWithMagneticTranslations& bosons)
{
/*  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrMomentum = bosons.NbrMomentum;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;
  this->Flag = bosons.Flag;*/
}

// destructor
//

BosonOnTorusWithMagneticTranslations::~BosonOnTorusWithMagneticTranslations ()
{
/*  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Keys;
      delete[] this->KeyMultiplicationTable;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;
      int Size = (this->MaxMomentum + 2) * this->IncNbrBosons;
      for (int i = 0; i < Size; ++i)
	{
	  if (this->KeyInvertSectorSize[i] > 0)
	    {
	      for (int j= 0; j < this->KeyInvertSectorSize[i]; ++j)
		delete[] this->KeyInvertIndices[i][j];
	      delete[] this->KeyInvertTable[i];
	      delete[] this->KeyInvertTableNbrIndices[i];
	      delete[] this->KeyInvertIndices[i];
	    }
	}
      
      delete[] this->KeyInvertSectorSize;
      delete[] this->KeyInvertTable;
      delete[] this->KeyInvertTableNbrIndices;
      delete[] this->KeyInvertIndices;
    }*/
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnTorusWithMagneticTranslations& BosonOnTorusWithMagneticTranslations::operator = (const BosonOnTorusWithMagneticTranslations& bosons)
{
/*  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;
    }
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrMomentum = bosons.NbrMomentum;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;
  this->Flag = bosons.Flag;*/
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnTorusWithMagneticTranslations::Clone()
{
  return new BosonOnTorusWithMagneticTranslations(*this);
}

// get momemtum value of a given state
//
// index = state index
// return value = state momentum

int BosonOnTorusWithMagneticTranslations::GetMomentumValue(int index)
{
/*  int StateMaxMomentum = this->StateMaxMomentum[index];
  int Momentum = 0;
  int* TmpStateDescription = this->StateDescription[index];
  for (int i = 0; i <= StateMaxMomentum; ++i)
    {
      Momentum += (TmpStateDescription[i] * i);
    }
  return (Momentum % this->MaxMomentum);*/
  return 0;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnTorusWithMagneticTranslations::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
/*  if (this->MomentumConstraintFlag == false)
    {
      for (int i = 0; i < this->MaxMomentum; ++i)
	TmpList += new PeriodicMomentumQuantumNumber (i, this->MaxMomentum);
    }
  else
    TmpList += new PeriodicMomentumQuantumNumber (this->MomentumConstraint, this->MaxMomentum);*/
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnTorusWithMagneticTranslations::GetQuantumNumber (int index)
{
/*  if (this->MomentumConstraintFlag == false)
    {
      return  new PeriodicMomentumQuantumNumber (this->GetMomentumValue(index), this->MaxMomentum);
    }
  else
    return new PeriodicMomentumQuantumNumber (this->MomentumConstraint, this->MaxMomentum);*/
  return 0;
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnTorusWithMagneticTranslations::ExtractSubspace (AbstractQuantumNumber& q, 
									     SubspaceSpaceConverter& converter)
{
  return 0;
/*  if (q.GetQuantumNumberType() != AbstractQuantumNumber::PeriodicMomentum)
    return 0;
  if (this->MomentumConstraintFlag == true)
    if (this->MomentumConstraint == ((PeriodicMomentumQuantumNumber&) q).GetMomentum())
      return this;
    else 
      return 0;*/
}

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

int BosonOnTorusWithMagneticTranslations::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
/*  int CurrentLzMax = this->StateMaxMomentum[index];
  int* State = this->StateDescription[index];
  if ((n1 > CurrentLzMax) || (n2 > CurrentLzMax) || (State[n1] == 0) || (State[n2] == 0) || ((n1 == n2) && (State[n1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = CurrentLzMax;
  if (NewLzMax < m1)
    NewLzMax = m1;
  if (NewLzMax < m2)
    NewLzMax = m2;
  int i = 0;
  int* TemporaryState = new int [NewLzMax + 1];
  for (; i <= CurrentLzMax; ++i)
    TemporaryState[i] = State[i];
  for (; i <= NewLzMax; ++i)
    TemporaryState[i] = 0;
  coefficient = TemporaryState[n2];
  --TemporaryState[n2];
  coefficient *= TemporaryState[n1];
  --TemporaryState[n1];
  ++TemporaryState[m2];
  coefficient *= TemporaryState[m2];
  ++TemporaryState[m1];
  coefficient *= TemporaryState[m1];
  coefficient = sqrt(coefficient);
  while (TemporaryState[NewLzMax] == 0)
    --NewLzMax;
  int DestIndex = this->FindStateIndex(TemporaryState, NewLzMax);
  delete[] TemporaryState;
  return DestIndex;*/
  return 0;
}

// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

int BosonOnTorusWithMagneticTranslations::FindStateIndex(int* stateDescription, int lzmax)
{
/*  int TmpKey = this->GenerateKey(stateDescription, lzmax);
  int Sector = lzmax * this->IncNbrBosons + stateDescription[lzmax];
  int TmpPos = 0;
  int TmpPos2 = this->KeyInvertSectorSize[Sector] - 1;
  int TmpPos3;
  int* TmpKeyInvertTable = this->KeyInvertTable[Sector];
  while (TmpPos2 != TmpPos)
    {
      TmpPos3 = (TmpPos2 + TmpPos) >> 1;
      if (TmpKey < TmpKeyInvertTable[TmpPos3])
	{
	  TmpPos2 = TmpPos3 - 1;
	}
      else
	if (TmpKey > TmpKeyInvertTable[TmpPos3])
	  {
	    TmpPos = TmpPos3 + 1;
	  }
	else
	  {
	    TmpPos2 = TmpPos3;
	    TmpPos = TmpPos3;		    
	  }
    }
  int i;
  int* TmpStateDescription;
  int Start;
  TmpPos3 = this->KeyInvertTableNbrIndices[Sector][TmpPos];
  int* TmpKeyInvertIndices = this->KeyInvertIndices[Sector][TmpPos];
  for (int j = 0; j < TmpPos3; ++j)
    {
      Start = TmpKeyInvertIndices[j];
      i = 0;
      TmpStateDescription = this->StateDescription[Start];
      while (i <= lzmax)
	{
	  if (stateDescription[i] != TmpStateDescription[i])
	    i = lzmax + 1;
	  ++i;
	}
      if (i == (lzmax + 1))
	return Start;
    }*/
  return this->HilbertSpaceDimension;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnTorusWithMagneticTranslations::PrintState (ostream& Str, int state)
{
  this->StateDescription[state].PrintState(Str, this->ReducedNbrState, this->RemainderNbrState);
/*  int Max = this->StateMaxMomentum[state];
  int i = 0;
  for (; i <= Max; ++i)
    Str << TmpState[i] << " ";
  for (; i < this->MaxMomentum; ++i)
    Str << "0 ";
  Str << " key = " << this->Keys[state] << " max momentum = " << this->MomentumMaxPosition[Max * (this->NbrBosons + 1) + TmpState[Max]]
      << " position = " << FindStateIndex(TmpState, Max) << " momentum = " << this->GetMomentumValue(state);*/
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a fermion in the state
// currentMaxMomentum = momentum maximum value for fermions that are still to be placed
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int BosonOnTorusWithMagneticTranslations::GenerateStates(int nbrBosons, int maxMomentum, int currentMaxMomentum, int pos)
{
  return 0;
}

// generate all states corresponding to the constraints
// 
// return value = hilbert space dimension

int BosonOnTorusWithMagneticTranslations::GenerateStates()
{
  this->CompatibilityWithXMomentum = new bool [this->IncNbrBosons];
  for (int i = 0; i <= this->NbrBosons; ++i)
    {
      if (((i * this->XMomentum) % this->MomentumModulo) == 0)
	this->CompatibilityWithXMomentum[i] = true;
      else
	this->CompatibilityWithXMomentum[i] = false;
    }
  int Lim = this->NbrBosons * this->MaxMomentum;
  this->CompatibilityWithYMomentum = new bool [Lim + 1];
  for (int i = 0; i <= Lim; ++i)
    {
      this->CompatibilityWithYMomentum[i]= false;	  
    }
  for (int i = this->YMomentum; i <= Lim; i += this->MomentumModulo)
    {
      this->CompatibilityWithYMomentum[i]= true;	  
    }

  int MaxHilbertSpaceDimension = this->HilbertSpaceDimension;
  this->TemporaryStateDescription = new BosonOnTorusState* [MaxHilbertSpaceDimension];
  BosonOnTorusState** InternalTemporaryStateDescription = new BosonOnTorusState* [MaxHilbertSpaceDimension];
  BosonOnTorusState* TmpState;
  this->HilbertSpaceDimension = -1;
  int NbrTranslation;
  if (this->CompatibilityWithYMomentum[0] == true)
    {
      TmpState = new BosonOnTorusState(this->ReducedNbrState);
      TmpState->SetOccupation (0, this->NbrBosons);
      TmpState->PutInCanonicalForm(this->TemporaryState, this->ReducedNbrState, this->RemainderNbrState, 
				   this->NbrMomentum, NbrTranslation, this->MomentumModulo);
      if (this->CompatibilityWithXMomentum[NbrTranslation] == true)
	{
	  InternalTemporaryStateDescription[0] = TmpState;
	  ++this->HilbertSpaceDimension;
	}
    }
  
  for (int i = 1; i <= this->MaxMomentum; ++i)
    {
      for (int j = 1; j < this->NbrBosons; ++j)
	{
	  int Pos = this->RawGenerateStates(this->NbrBosons - j, i - 1, 0, 0, j * i);
	  for (int k = 0; k < Pos; ++k)
	    {
	      this->TemporaryStateDescription[k]->SetOccupation (i, j);
	      this->TemporaryStateDescription[k]->PutInCanonicalForm(this->TemporaryState, this->ReducedNbrState, this->RemainderNbrState, 
								     this->NbrMomentum, NbrTranslation, this->MomentumModulo);
	      if ((this->CompatibilityWithXMomentum[NbrTranslation] == true) && 
		  ((this->HilbertSpaceDimension < 0) || (InternalTemporaryStateDescription[this->HilbertSpaceDimension]->Lesser(*(this->TemporaryStateDescription[k]), 
																this->ReducedNbrState))))
		{
		  ++this->HilbertSpaceDimension;
		  InternalTemporaryStateDescription[this->HilbertSpaceDimension] = this->TemporaryStateDescription[k];
		  this->TemporaryStateDescription[k] = new BosonOnTorusState (this->ReducedNbrState);
		}
	      else
		this->TemporaryStateDescription[k]->EmptyState(this->ReducedNbrState);
	    }
	}
      if (this->CompatibilityWithYMomentum[i * this->NbrBosons] == true)
	{
	  TmpState = new BosonOnTorusState(this->ReducedNbrState);
	  TmpState->SetOccupation (i, this->NbrBosons);
	  TmpState->PutInCanonicalForm(this->TemporaryState, this->ReducedNbrState, this->RemainderNbrState, 
				       this->NbrMomentum, NbrTranslation, this->MomentumModulo);
	  if ((this->CompatibilityWithXMomentum[NbrTranslation] == true) && 
	      ((this->HilbertSpaceDimension < 0) || (InternalTemporaryStateDescription[this->HilbertSpaceDimension]->Lesser(*TmpState, this->ReducedNbrState))))
	    {
	      ++this->HilbertSpaceDimension;
	      InternalTemporaryStateDescription[this->HilbertSpaceDimension] = TmpState;
	    }
	}
    }  

  ++this->HilbertSpaceDimension;
  for (int i = 0; i < MaxHilbertSpaceDimension; ++i)
    delete this->TemporaryStateDescription[i];
  delete[] this->TemporaryStateDescription;
  this->TemporaryStateDescription = 0;
  if (this->HilbertSpaceDimension > 0)
    {
      this->StateDescription = new BosonOnTorusState [this->HilbertSpaceDimension];
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  this->StateDescription[i].TransfertState(*(InternalTemporaryStateDescription[i]));
	  this->StateMaxMomentum[i] = this->StateDescription[i].GetHighestIndex(this->ReducedNbrState);
	}
    }
  for (int i = 0; i < MaxHilbertSpaceDimension; ++i)
    delete InternalTemporaryStateDescription[i];
  delete[] InternalTemporaryStateDescription;
  return this->HilbertSpaceDimension;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson in the state
// currentMaxMomentum = momentum maximum value for bosons that are still to be placed
// pos = position in StateDescription array where to store states
// currentYMomentum = current value of the momentum in the y direction
// return value = position from which new states have to be stored

int BosonOnTorusWithMagneticTranslations::RawGenerateStates(int nbrBosons, int maxMomentum, int currentMaxMomentum, int pos, int currentYMomentum)
{
  if (nbrBosons == 0)
    {
      if (this->CompatibilityWithYMomentum[currentYMomentum] == true)
	{
	  return pos + 1;
	}
      else
	{
 	  return pos;
	}
    }
  if (currentMaxMomentum == maxMomentum)
    {
      if (this->CompatibilityWithYMomentum[currentYMomentum + (maxMomentum * nbrBosons)] == true)
	{
	  this->TemporaryStateDescription[pos]->SetOccupation (maxMomentum, nbrBosons);
	  return pos + 1;
	}
      else
	return pos;
    }

  int TmpNbrBosons = nbrBosons;
  int IncCurrentMaxMomentum = currentMaxMomentum + 1;
  int TmpPos = pos;
  while (TmpNbrBosons > 0)
    {
      TmpPos = this->RawGenerateStates(nbrBosons - TmpNbrBosons, maxMomentum, IncCurrentMaxMomentum, pos, currentYMomentum + (TmpNbrBosons * currentMaxMomentum));
      for (int i = pos; i < TmpPos; i++)
	this->TemporaryStateDescription[i]->SetOccupation (currentMaxMomentum, TmpNbrBosons);      
      --TmpNbrBosons;
      pos = TmpPos;
    }
  if (currentMaxMomentum == 0)
    {
      return pos;
    }
  else
    return this->RawGenerateStates(nbrBosons, maxMomentum, IncCurrentMaxMomentum, pos, currentYMomentum);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnTorusWithMagneticTranslations::GenerateLookUpTable(int memory)
{
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// return value = Hilbert space dimension

int BosonOnTorusWithMagneticTranslations::EvaluateHilbertSpaceDimension(int nbrBosons, int maxMomentum)
{
  FactorialCoefficient Dimension; 
  Dimension.PartialFactorialMultiply(maxMomentum, maxMomentum + nbrBosons - 1); 
  Dimension.FactorialDivide(nbrBosons);
  return Dimension.GetIntegerValue();
}

