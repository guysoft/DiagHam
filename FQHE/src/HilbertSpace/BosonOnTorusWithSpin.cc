////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of bosons on a torus with spin                   //
//                                                                            //
//                        last modification : 03/04/2012                      //
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
#include "HilbertSpace/BosonOnTorusWithSpin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "HilbertSpace/BosonOnTorusShort.h" 
#include "GeneralTools/ArrayTools.h"
#include "Vector/ComplexVector.h"

#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;


// constructor with a constraint on the total Ky momentum
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

BosonOnTorusWithSpin::BosonOnTorusWithSpin (int nbrBosons, int maxMomentum, int kyMomentum, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->SzFlag = false;
  this->TotalLz = maxMomentum;
  this->TotalSpin = 0;
  this->NbrBosonsUp = 0;
  this->NbrBosonsDown = 0;
  this->KyMomentum = kyMomentum;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = maxMomentum;
  this->Flag.Initialize();
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];

  this->NUpLzMax = this->LzMax + this->NbrBosons;
  this->NDownLzMax = this->LzMax + this->NbrBosons;
  this->FermionicLzMax = this->NUpLzMax;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
      this->Flag.Initialize();
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, 0, this->LzMax + this->NbrBosons + 1, this->LzMax + this->NbrBosons + 1, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      cout  << "Dimension = " << this->LargeHilbertSpaceDimension << endl;
      for (long i = 0; i < TmpLargeHilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState = this->StateDescriptionUp[i];
	  unsigned long Tmp = 0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & 0x1ul;
	  this->StateDescriptionUp[i] >>= this->NbrBosons - Tmp; 
	  TmpState = this->StateDescriptionDown[i];
	  Tmp = 0l;
	  for (int j = 0; j <= this->FermionicLzMax; ++j)
	    Tmp += (TmpState >> j) & 0x1ul;
	  this->StateDescriptionDown[i] >>= this->NbrBosons - Tmp; 
	}
      SortDoubleElementArrayDownOrdering<unsigned long>(this->StateDescriptionUp, this->StateDescriptionDown, TmpLargeHilbertSpaceDimension);
      this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (4 * sizeof(unsigned long));
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
}

// constructor with a constraint on total spin momentum and total momentum
// 
// nbrBosons = number of bosons
// maxMomentum = momentum maximum value for a boson
// totalSpin = twice the total spin along z
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

BosonOnTorusWithSpin::BosonOnTorusWithSpin (int nbrBosons, int maxMomentum, int totalSpin, int kyMomentum, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->SzFlag = true;
  this->TotalLz = maxMomentum;
  this->TotalSpin = totalSpin;
  this->NbrBosonsUp = (this->NbrBosons + this->TotalSpin) / 2;
  this->NbrBosonsDown = (this->NbrBosons - this->TotalSpin) / 2;
  this->KyMomentum = kyMomentum;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = maxMomentum;
  this->Flag.Initialize();
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];

  this->NUpLzMax = this->LzMax + this->NbrBosons;
  this->NDownLzMax = this->LzMax + this->NbrBosons;
  this->FermionicLzMax = this->NUpLzMax;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, 0, this->NbrBosonsUp);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
      this->Flag.Initialize();
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, 0, this->LzMax + this->NbrBosonsUp + 1, this->LzMax + this->NbrBosonsDown + 1, this->NbrBosonsUp, 0l);
      cout  << "Dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
      SortDoubleElementArrayDownOrdering<unsigned long>(this->StateDescriptionUp, this->StateDescriptionDown, TmpLargeHilbertSpaceDimension);
      this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
      long UsedMemory = 0;
      UsedMemory += (long) this->HilbertSpaceDimension * (4 * sizeof(unsigned long));
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
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnTorusWithSpin::BosonOnTorusWithSpin(const BosonOnTorusWithSpin& bosons)
{
  this->KyMomentum = bosons.KyMomentum;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->SzFlag = bosons.SzFlag;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpLzMax = bosons.NUpLzMax;
  this->NDownLzMax = bosons.NDownLzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;
}

// destructor
//

BosonOnTorusWithSpin::~BosonOnTorusWithSpin ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnTorusWithSpin& BosonOnTorusWithSpin::operator = (const BosonOnTorusWithSpin& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescriptionUp;
      delete[] this->StateDescriptionDown;
      delete[] this->UniqueStateDescriptionUp;
      delete[] this->UniqueStateDescriptionSubArraySizeUp;
      delete[] this->FirstIndexUniqueStateDescriptionUp;
    }
  delete[] this->TemporaryStateUp;
  delete[] this->TemporaryStateDown;
  delete[] this->ProdATemporaryStateUp;
  delete[] this->ProdATemporaryStateDown;
  this->KyMomentum = bosons.KyMomentum;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->SzFlag = bosons.SzFlag;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpLzMax = bosons.NUpLzMax;
  this->NDownLzMax = bosons.NDownLzMax;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnTorusWithSpin::Clone()
{
  return new BosonOnTorusWithSpin(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnTorusWithSpin::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescriptionUp[state], this->StateDescriptionDown[state],
		       this->TemporaryStateUp, this->TemporaryStateDown); 

  unsigned long Tmp;
  Str << " | ";
  for (int i = 0; i <= this->LzMax; ++i)
    {
      Str << "(" << this->TemporaryStateUp[i] << "," << this->TemporaryStateDown[i] << ") | ";
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// currentFermionicPositionUp = current fermionic position within the state description for the spin up
// currentFermionicPositionDown = current fermionic position within the state description for the spin down
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnTorusWithSpin::GenerateStates(int nbrBosons, int currentKy, int currentTotalKy, int currentFermionicPositionUp, int currentFermionicPositionDown, long pos)
{
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->KyMomentum)
	{
 	  this->StateDescriptionUp[pos] = 0x0ul;
 	  this->StateDescriptionDown[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKy < 0)
    return pos;
  for (int i = nbrBosons; i >= 0; --i)
    {
      unsigned long MaskUp = ((0x1ul << i) - 0x1ul) << (currentFermionicPositionUp - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
	  long TmpPos = this->GenerateStates(nbrBosons - i - j, currentKy - 1, currentTotalKy + ((i + j) * currentKy), currentFermionicPositionUp - i - 1, currentFermionicPositionDown - j - 1, pos);
	  unsigned long MaskDown = ((0x1ul << j) - 0x1ul) << (currentFermionicPositionDown - j - 1);
	  for (; pos < TmpPos; ++pos)
	    {
	      this->StateDescriptionUp[pos] |= MaskUp;
	      this->StateDescriptionDown[pos] |= MaskDown;
	    }
	}
    }
  return pos;
};

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// currentFermionicPositionUp = current fermionic position within the state description for the spin up
// currentFermionicPositionDown = current fermionic position within the state description for the spin down
// nbrSpinUp = number of particles with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnTorusWithSpin::GenerateStates(int nbrBosons, int currentKy, int currentTotalKy, 
					  int currentFermionicPositionUp, int currentFermionicPositionDown, int nbrSpinUp, long pos)
{
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrBosons))
    return 0l;
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->KyMomentum)
	{
 	  this->StateDescriptionUp[pos] = 0x0ul;
 	  this->StateDescriptionDown[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKy < 0)
    return pos;
  for (int i = nbrSpinUp; i >= 0; --i)
    {
      unsigned long MaskUp = ((0x1ul << i) - 0x1ul) << (currentFermionicPositionUp - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
	  long TmpPos = this->GenerateStates(nbrBosons - i - j, currentKy - 1, currentTotalKy + ((i + j) * currentKy), currentFermionicPositionUp - i - 1, currentFermionicPositionDown - j - 1, nbrSpinUp - i, pos);
	  unsigned long MaskDown = ((0x1ul << j) - 0x1ul) << (currentFermionicPositionDown - j - 1);
	  for (; pos < TmpPos; ++pos)
	    {
	      this->StateDescriptionUp[pos] |= MaskUp;
	      this->StateDescriptionDown[pos] |= MaskDown;
	    }
	}
    }
  return pos;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long BosonOnTorusWithSpin::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy)
{
  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->KyMomentum)
	return 1l;
      else	
	return 0l;
    }
  if (currentKy < 0)
    return 0l;
  long Count = 0;
  if (nbrBosons == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if (((j + currentTotalKy) % this->NbrLzValue) == this->KyMomentum)
	    Count += 2l;
	}
      return Count;
    }
  for (int i = nbrBosons; i >= 0; --i)
    Count += (((long) i) + 1l) * this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKy - 1, currentTotalKy + (i * currentKy));
  return Count;
}

// evaluate Hilbert space dimension for a given total spin momentum
//
// nbrBosons = number of bosons
// currentKy = current momentum along y for a single particle
// currentTotalKy = current total momentum along y
// nbrSpinUp = number of particles with spin up
// return value = Hilbert space dimension

long BosonOnTorusWithSpin::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKy, int currentTotalKy, int nbrSpinUp)
{
  if ((nbrSpinUp < 0) || (nbrSpinUp > nbrBosons))
    return 0l;

  if (nbrBosons == 0)
    {
      if ((currentTotalKy % this->NbrLzValue) == this->KyMomentum)
	return 1l;
      else	
	return 0l;
    }
  if (currentKy < 0)
    return 0l;
  long Count = 0;
  if (nbrBosons == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if (((j + currentTotalKy) % this->NbrLzValue) == this->KyMomentum)
	    Count++;
	}
      return Count;
    }
  for (int i = nbrBosons; i >= 0; --i)
    for (int j = i; j >= 0; --j)
      Count += this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKy - 1, currentTotalKy + (i * currentKy), nbrSpinUp -j);
  return Count;
}


// project out any configurations that have particles on levels other than lll
//
// inputVector = vector to apply the projection to
// outputVector = projected vector
// finalSpace = reference to output vector space

void BosonOnTorusWithSpin::ProjectionInTheLowestLevel(RealVector &inputVector, RealVector & outputVector, BosonOnTorusShort *finalSpace)
{
  unsigned long Etat;
  int Idx; 
  for(int i = 0 ; i < finalSpace->GetHilbertSpaceDimension() ; i++)
    {
      Etat = finalSpace->StateDescription[i]; 
      Idx = this->FindStateIndex(0,Etat);
      if ( Idx < this->HilbertSpaceDimension ) 
	{
	  outputVector[i] = inputVector[Idx];
	}
    }
  cout <<"Norm after projection" <<outputVector.Norm()<<endl;
}


void BosonOnTorusWithSpin::ProjectionInTheLowestLevel(ComplexVector &inputVector, ComplexVector & outputVector, BosonOnTorusShort *finalSpace)
{
  unsigned long Etat;
  int Idx; 
  for(int i = 0 ; i < finalSpace->GetHilbertSpaceDimension() ; i++)
    {
      Etat = finalSpace->StateDescription[i]; 
      Idx = this->FindStateIndex(Etat,0);
      if ( Idx < this->HilbertSpaceDimension ) 
	{
	  cout <<i<< " "<<Idx<<endl;
	  outputVector[i] = inputVector[Idx];
	}
    }
  cout <<"Norm after projection" <<outputVector.Norm()<<endl;
}
