////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//     class for bosons on sphere with SU(2) spin and Sz<->-Sz symmetry       //
//                                                                            //
//                        last modification : 23/09/2016                      //
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
#include "HilbertSpace/BosonOnSphereWithSU2SpinSzSymmetry.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/Permutations.h"

#include <math.h>
#include <cstdlib>
#include <algorithm>
#include <map>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
// 

BosonOnSphereWithSU2SpinSzSymmetry::BosonOnSphereWithSU2SpinSzSymmetry ()
{
  this->SzParitySign = 1.0;
}

// basic constructor without any constraint on Sz
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity

BosonOnSphereWithSU2SpinSzSymmetry::BosonOnSphereWithSU2SpinSzSymmetry (int nbrBosons, int totalLz, int lzMax, bool minusSzParity)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  if (minusSzParity == false)
    this->SzParitySign = 1.0;
  else
    this->SzParitySign = -1.0;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;

  this->NbrBosonsUp = 0;
  this->NbrBosonsDown = 0;
  this->NUpLzMax = this->LzMax + this->NbrBosons - 1;
  this->NDownLzMax = this->LzMax + this->NbrBosons - 1;
  this->FermionicLzMax = this->NUpLzMax;
  if (this->NDownLzMax > this->FermionicLzMax)
    this->FermionicLzMax = this->NDownLzMax;
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1);
  
  this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, 
						       this->FermionicLzMax, this->FermionicLzMax, 0l);

  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2SpinSzSymmetry!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->GenerateSatetsWithDiscreateSymmetry();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->TargetSpace = this;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  

  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(10000000);
      
#ifdef __DEBUG__
      int UsedMemory = 0;
      UsedMemory += this->HilbertSpaceDimension * (4 * sizeof(unsigned long));
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

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson
// totalSpin = twice the total spin value (not taken into account)
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity
// memory = amount of memory granted for precalculations

BosonOnSphereWithSU2SpinSzSymmetry::BosonOnSphereWithSU2SpinSzSymmetry (int nbrBosons, int totalLz, int lzMax, int totalSpin, bool minusSzParity, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  if (minusSzParity == false)
    this->SzParitySign = 1.0;
  else
    this->SzParitySign = -1.0;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;

  this->NbrBosonsUp = this->NbrBosons + this->TotalSpin;
  this->NbrBosonsDown = this->NbrBosons - this->TotalSpin;
  if ((this->NbrBosonsUp < 0) || ((this->NbrBosonsUp & 0x1) != 0) ||
      (this->NbrBosonsDown < 0) || ((this->NbrBosonsDown & 0x1) != 0))
    this->LargeHilbertSpaceDimension = 0l;
  else
    {
      this->NbrBosonsUp >>= 1;
      this->NbrBosonsDown >>= 1;
      this->NUpLzMax = this->LzMax + this->NbrBosonsUp - 1;
      this->NDownLzMax = this->LzMax + this->NbrBosonsDown - 1;
      this->FermionicLzMax = this->NUpLzMax;
      if (this->NDownLzMax > this->FermionicLzMax)
	this->FermionicLzMax = this->NDownLzMax;
      this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, this->NbrBosonsUp, this->NbrBosonsDown);
    }
  this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, 
						       this->NbrBosonsUp, this->NbrBosonsDown, 0l);
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2SpinSzSymmetry!" << endl;
      exit(1);
    }

  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->GenerateSatetsWithDiscreateSymmetry();
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  

  this->TargetSpace = this;

  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->GenerateLookUpTable(memory);
      
#ifdef __DEBUG__
      int UsedMemory = 0;
      UsedMemory += this->HilbertSpaceDimension * (4 * sizeof(unsigned long));
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

BosonOnSphereWithSU2SpinSzSymmetry::BosonOnSphereWithSU2SpinSzSymmetry(const BosonOnSphereWithSU2SpinSzSymmetry& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpLzMax = bosons.NUpLzMax;
  this->NDownLzMax = bosons.NDownLzMax;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->TotalSpin = bosons.TotalSpin;
  this->SzParitySign = bosons.SzParitySign;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

BosonOnSphereWithSU2SpinSzSymmetry::~BosonOnSphereWithSU2SpinSzSymmetry ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereWithSU2SpinSzSymmetry& BosonOnSphereWithSU2SpinSzSymmetry::operator = (const BosonOnSphereWithSU2SpinSzSymmetry& bosons)
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
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalSpin = bosons.TotalSpin;
  this->SzParitySign = bosons.SzParitySign;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->NUpLzMax = bosons.NUpLzMax;
  this->NDownLzMax = bosons.NDownLzMax;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->FermionicLzMax = bosons.FermionicLzMax;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->TemporaryStateSigma[0] = this->TemporaryStateUp;
  this->TemporaryStateSigma[1] = this->TemporaryStateDown;
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateSigma[0] = this->ProdATemporaryStateUp;
  this->ProdATemporaryStateSigma[1] = this->ProdATemporaryStateDown;
  this->StateDescriptionUp = bosons.StateDescriptionUp;
  this->StateDescriptionDown = bosons.StateDescriptionDown;
  this->UniqueStateDescriptionUp = bosons.UniqueStateDescriptionUp;
  this->UniqueStateDescriptionSubArraySizeUp = bosons.UniqueStateDescriptionSubArraySizeUp;
  this->StateDescriptionSigma[0] = this->StateDescriptionUp;
  this->StateDescriptionSigma[1] = this->StateDescriptionDown;
  this->NbrUniqueStateDescriptionUp = bosons.NbrUniqueStateDescriptionUp;
  this->FirstIndexUniqueStateDescriptionUp = bosons.FirstIndexUniqueStateDescriptionUp;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereWithSU2SpinSzSymmetry::Clone()
{
  return new BosonOnSphereWithSU2SpinSzSymmetry(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnSphereWithSU2SpinSzSymmetry::SetTargetSpace(ParticleOnSphereWithSpin* targetSpace)
{
  this->TargetSpace = (BosonOnSphereWithSU2SpinSzSymmetry*) targetSpace;
}

// generate the Hilbert space with the discrete symmetry constraint
//

void BosonOnSphereWithSU2SpinSzSymmetry::GenerateSatetsWithDiscreateSymmetry()

{
  long TmpHilbertSpaceDimension = 0l;
  if (this->SzParitySign > 0.0)
    {
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if (this->StateDescriptionUp[i] >= this->StateDescriptionDown[i])
	    {
	      ++TmpHilbertSpaceDimension;
	    }
	}
    }
  else
    {
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  if (this->StateDescriptionUp[i] > this->StateDescriptionDown[i])
	    {
	      ++TmpHilbertSpaceDimension;
	    }
	}
    }
  if (TmpHilbertSpaceDimension > 0l)
    {
      unsigned long* TmpStateDescriptionUp = new unsigned long [TmpHilbertSpaceDimension];
      unsigned long* TmpStateDescriptionDown = new unsigned long [TmpHilbertSpaceDimension];
      TmpHilbertSpaceDimension = 0;
      if (this->SzParitySign > 0.0)
	{
	  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      if (this->StateDescriptionUp[i] >= this->StateDescriptionDown[i])
		{
		  TmpStateDescriptionUp[TmpHilbertSpaceDimension] = this->StateDescriptionUp[i];
		  TmpStateDescriptionDown[TmpHilbertSpaceDimension] = this->StateDescriptionDown[i];
		  ++TmpHilbertSpaceDimension;
		}
	    }
	}
      else
	{
	  for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      if (this->StateDescriptionUp[i] > this->StateDescriptionDown[i])
		{
		  TmpStateDescriptionUp[TmpHilbertSpaceDimension] = this->StateDescriptionUp[i];
		  TmpStateDescriptionDown[TmpHilbertSpaceDimension] = this->StateDescriptionDown[i];
		  ++TmpHilbertSpaceDimension;
		}
	    }	 
	}	 
      delete[] this->StateDescriptionUp;
      delete[] this->StateDescriptionDown;
      this->StateDescriptionUp = TmpStateDescriptionUp;
      this->StateDescriptionDown = TmpStateDescriptionDown;
    }
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
}

// find state index
//
// stateDescriptionUp = unsigned integer describing the fermionic state for type up particles
// stateDescriptionDown = unsigned integer describing the fermionic state for type down particles
// return value = corresponding index

int BosonOnSphereWithSU2SpinSzSymmetry::FindStateIndex(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown)
{
  int PosMin = 0;
  int PosMax = this->NbrUniqueStateDescriptionUp - 1;
  int PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->UniqueStateDescriptionUp[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescriptionUp))
    {
       if (CurrentState > stateDescriptionUp)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->UniqueStateDescriptionUp[PosMid];
    }
  if (CurrentState != stateDescriptionUp)
    {
      PosMid = PosMax;
      if (this->UniqueStateDescriptionUp[PosMid] != stateDescriptionUp)
	return this->HilbertSpaceDimension;
     }

  PosMin = this->FirstIndexUniqueStateDescriptionUp[PosMid];
  PosMax = PosMin + this->UniqueStateDescriptionSubArraySizeUp[PosMid] - 1;
  PosMid = (PosMin + PosMax) >> 1;
  CurrentState = this->StateDescriptionDown[PosMid];
  while ((PosMax > PosMin) && (CurrentState != stateDescriptionDown))
    {
       if (CurrentState > stateDescriptionDown)
	 {
	   PosMin = PosMid + 1;
	 }
       else
 	{
 	  PosMax = PosMid - 1;
	} 
       PosMid = (PosMin + PosMax) >> 1;
       CurrentState = this->StateDescriptionDown[PosMid];
    }
  if (CurrentState == stateDescriptionDown)
    return PosMid;
  if ( this->StateDescriptionDown[PosMax] != stateDescriptionDown)
    return this->HilbertSpaceDimension;
  else
    return PosMax;
}


// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2SpinSzSymmetry::AduAu (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  if (this->TemporaryStateUp[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  if (this->StateDescriptionUp[index] == this->StateDescriptionDown[index])
    this->ProdATemporaryOrbitFactor = 1.0;
  else
    this->ProdATemporaryOrbitFactor = M_SQRT2;
  coefficient = (double) this->TemporaryStateUp[n];
  --this->TemporaryStateUp[n];
  ++this->TemporaryStateUp[m];
  coefficient *= (double) this->TemporaryStateUp[m];
  coefficient = sqrt(coefficient);  
  unsigned long TmpUp = this->BosonToFermion(this->TemporaryStateUp);
  return this->SymmetrizeAdAdResult(TmpUp, this->StateDescriptionDown[index], coefficient);  
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2SpinSzSymmetry::AduAd (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  if (this->TemporaryStateDown[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  if (this->StateDescriptionUp[index] == this->StateDescriptionDown[index])
    this->ProdATemporaryOrbitFactor = 1.0;
  else
    this->ProdATemporaryOrbitFactor = M_SQRT2;
  coefficient = (double) this->TemporaryStateDown[n];
  --this->TemporaryStateDown[n];
  ++this->TemporaryStateUp[m];
  coefficient *= (double) this->TemporaryStateUp[m];
  coefficient = sqrt(coefficient); 
  unsigned long TmpUp = this->BosonToFermion(this->TemporaryStateUp);
  unsigned long TmpDown = this->BosonToFermion(this->TemporaryStateDown);
  return this->SymmetrizeAdAdResult(TmpUp, TmpDown, coefficient);  
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2SpinSzSymmetry::AddAu (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  if (this->TemporaryStateUp[n] == 0)
    { 
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;      
    }
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  if (this->StateDescriptionUp[index] == this->StateDescriptionDown[index])
    this->ProdATemporaryOrbitFactor = 1.0;
  else
    this->ProdATemporaryOrbitFactor = M_SQRT2;
  coefficient = (double) this->TemporaryStateUp[n];
  --this->TemporaryStateUp[n];
  ++this->TemporaryStateDown[m];
  coefficient *= (double) this->TemporaryStateDown[m];
  coefficient = sqrt(coefficient);  
  unsigned long TmpUp = this->BosonToFermion(this->TemporaryStateUp);
  unsigned long TmpDown = this->BosonToFermion(this->TemporaryStateDown);
  return this->SymmetrizeAdAdResult(TmpUp, TmpDown, coefficient);  
}

// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2SpinSzSymmetry::AddAd (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  if (this->TemporaryStateDown[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  if (this->StateDescriptionUp[index] == this->StateDescriptionDown[index])
    this->ProdATemporaryOrbitFactor = 1.0;
  else
    this->ProdATemporaryOrbitFactor = M_SQRT2;
  coefficient = (double) this->TemporaryStateDown[n];
  --this->TemporaryStateDown[n];
  ++this->TemporaryStateDown[m];
  coefficient *= (double) this->TemporaryStateDown[m];
  coefficient = sqrt(coefficient);  
  unsigned long TmpDown = this->BosonToFermion(this->TemporaryStateDown);
  return this->SymmetrizeAdAdResult(this->StateDescriptionUp[index], TmpDown, coefficient);  
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU2SpinSzSymmetry::AuAu (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  if ((this->ProdATemporaryStateUp[n1] == 0) || (this->ProdATemporaryStateUp[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryStateUp[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  if (this->StateDescriptionUp[index] == this->StateDescriptionDown[index])
    this->ProdATemporaryOrbitFactor = 1.0;
  else
    this->ProdATemporaryOrbitFactor = M_SQRT2;
  double Coefficient = this->ProdATemporaryStateUp[n2];
  --this->ProdATemporaryStateUp[n2];
  Coefficient *= this->ProdATemporaryStateUp[n1];
  --this->ProdATemporaryStateUp[n1];
  return sqrt(Coefficient);
}

// apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU2SpinSzSymmetry::AuAd (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateUp[n1] == 0) || (this->ProdATemporaryStateDown[n2] == 0))
    {
      return 0.0;
    }
  if (this->StateDescriptionUp[index] == this->StateDescriptionDown[index])
    this->ProdATemporaryOrbitFactor = 1.0;
  else
    this->ProdATemporaryOrbitFactor = M_SQRT2;
  double Coefficient = this->ProdATemporaryStateDown[n2];
  --this->ProdATemporaryStateDown[n2];
  Coefficient *= this->ProdATemporaryStateUp[n1];
  --this->ProdATemporaryStateUp[n1];
  return sqrt(Coefficient);
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU2SpinSzSymmetry::AdAd (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateDown[n1] == 0) || (this->ProdATemporaryStateDown[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryStateDown[n1] == 1)))    
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  if (this->StateDescriptionUp[index] == this->StateDescriptionDown[index])
    this->ProdATemporaryOrbitFactor = 1.0;
  else
    this->ProdATemporaryOrbitFactor = M_SQRT2;
  double Coefficient = this->ProdATemporaryStateDown[n2];
  --this->ProdATemporaryStateDown[n2];
  Coefficient *= this->ProdATemporaryStateDown[n1];
  --this->ProdATemporaryStateDown[n1];
  return sqrt(Coefficient);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphereWithSU2SpinSzSymmetry::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  double Coefficient = 1.0;
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  if (this->StateDescriptionUp[index] == this->StateDescriptionDown[index])
    this->ProdATemporaryOrbitFactor = 1.0;
  else
    this->ProdATemporaryOrbitFactor = M_SQRT2;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (spinIndices[i] == 0)
	{
	  unsigned long& Tmp = this->ProdATemporaryStateDown[n[i]];
	  if (Tmp == 0x0ul)
	    {
	      return 0.0;
	    }
	  Coefficient *= Tmp;
	  --Tmp;
	}
      else
	{
	  unsigned long& Tmp = this->ProdATemporaryStateUp[n[i]];
	  if (Tmp == 0x0ul)
	    {
	      return 0.0;
	    }
	  Coefficient *= Tmp;
	  --Tmp;
	}
    }
  return sqrt(Coefficient);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2SpinSzSymmetry::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      this->TemporaryStateUp[i] = this->ProdATemporaryStateUp[i];
      this->TemporaryStateDown[i] = this->ProdATemporaryStateDown[i];
    }
  coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (spinIndices[i] == 0)
	{
	  unsigned long& Tmp = this->TemporaryStateDown[m[i]];
	  ++Tmp;
	  coefficient *= Tmp;
	}
      else
	{
	  unsigned long& Tmp = this->TemporaryStateUp[m[i]];
	  ++Tmp;
	  coefficient *= Tmp;
	}
    }
  coefficient = sqrt(coefficient);
  return this->SymmetrizeAdAdResult(this->TemporaryStateUp, this->TemporaryStateDown, coefficient);  
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated
// szSector = Sz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix BosonOnSphereWithSU2SpinSzSymmetry::EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int lzSector, int szSector, RealVector& groundState, bool removeBinomialCoefficient)
{
  int nbrOrbitalA = this->LzMax + 1;
  int nbrOrbitalB = this->LzMax + 1;  
  if (nbrParticleSector == 0)
    {
      if ((lzSector == 0) && (szSector == 0))
	{
	  RealMatrix TmpEntanglementMatrix(1, this->HilbertSpaceDimension, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    TmpEntanglementMatrix.SetMatrixElement(0, i, groundState[i]);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  if (nbrParticleSector == this->NbrBosons)
    {
      if ((lzSector == this->TotalLz) && (szSector == this->TotalSpin))
	{
	  RealMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    TmpEntanglementMatrix.SetMatrixElement(i, 0, groundState[i]);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    } 
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;
  int ComplementaryLzSector = this->TotalLz - lzSector;
  if ((abs(ComplementarySzSector) > ComplementaryNbrParticles) || (abs(ComplementaryLzSector) > (this->LzMax * ComplementaryNbrParticles)))
    {
      RealMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;  
    }

  double* LogFactorials = new double[this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2 ; i <= this->NbrBosons; ++i)
    LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
  double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrParticles] - LogFactorials[nbrParticleSector];
  if (removeBinomialCoefficient == true)
    TmpLogBinomial = 0.0;


  BosonOnSphereWithSU2Spin SubsytemSpace(nbrParticleSector, lzSector, this->LzMax, szSector);
  BosonOnSphereWithSU2Spin ComplementarySubsytemSpace(ComplementaryNbrParticles, ComplementaryLzSector, this->LzMax, ComplementarySzSector);  
  RealMatrix TmpEntanglementMatrix(SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySubsytemSpace.GetHilbertSpaceDimension(), true);
  long TmpNbrNonZeroElements = 0l;
  unsigned long** TmpSubsytemSpaceOccupationNumbersUp = new unsigned long* [SubsytemSpace.GetHilbertSpaceDimension()];
  unsigned long** TmpSubsytemSpaceOccupationNumbersDown = new unsigned long* [SubsytemSpace.GetHilbertSpaceDimension()];
  unsigned long** TmpSubsytemSpaceMonomialUp = new unsigned long* [SubsytemSpace.GetHilbertSpaceDimension()];
  unsigned long** TmpSubsytemSpaceMonomialDown = new unsigned long* [SubsytemSpace.GetHilbertSpaceDimension()];
  double* TmpSubsytemLogFactorials = new double [SubsytemSpace.GetHilbertSpaceDimension()];
  unsigned long* TmpMonomialUp1 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpMonomialDown1 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpMonomialUp3 = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomialDown3 = new unsigned long [this->NbrBosons];

  for (int i = 0; i < SubsytemSpace.GetHilbertSpaceDimension(); ++i)
    {
      TmpSubsytemSpaceOccupationNumbersUp[i] = new unsigned long [this->NbrLzValue];
      TmpSubsytemSpaceOccupationNumbersDown[i] = new unsigned long [this->NbrLzValue];
      TmpSubsytemSpaceMonomialUp[i] = new unsigned long [nbrParticleSector];
      TmpSubsytemSpaceMonomialDown[i] = new unsigned long [nbrParticleSector];
      SubsytemSpace.FermionToBoson(SubsytemSpace.StateDescriptionUp[i], SubsytemSpace.StateDescriptionDown[i], 
				      TmpSubsytemSpaceOccupationNumbersUp[i], TmpSubsytemSpaceOccupationNumbersDown[i]);
      SubsytemSpace.ConvertToMonomial(SubsytemSpace.StateDescriptionUp[i], SubsytemSpace.StateDescriptionDown[i], 
					 TmpSubsytemSpaceMonomialUp[i], TmpSubsytemSpaceMonomialDown[i]);
      unsigned long* TmpOccupationNumberUp = TmpSubsytemSpaceOccupationNumbersUp[i];
      unsigned long* TmpOccupationNumberDown = TmpSubsytemSpaceOccupationNumbersDown[i];
      double TmpFactor = 0.0;
      for (int k = 0; k <= SubsytemSpace.LzMax; ++k)
	{
	  TmpFactor += LogFactorials[TmpOccupationNumberUp[k]];
	  TmpFactor += LogFactorials[TmpOccupationNumberDown[k]];
	}
      TmpSubsytemLogFactorials[i] = TmpFactor;      
    }
  for (int MinIndex = 0; MinIndex < ComplementarySubsytemSpace.GetHilbertSpaceDimension(); ++MinIndex)    
    {
      ComplementarySubsytemSpace.ConvertToMonomial(ComplementarySubsytemSpace.StateDescriptionUp[MinIndex], ComplementarySubsytemSpace.StateDescriptionDown[MinIndex], 
						   TmpMonomialUp1, TmpMonomialDown1);
       for (int k = 0; k < ComplementaryNbrParticles; ++k)
	 {
	   TmpMonomialUp1[k] += this->LzMax + 1 - nbrOrbitalA;
	   TmpMonomialDown1[k] += this->LzMax + 1 - nbrOrbitalA;
	 }   
      ComplementarySubsytemSpace.FermionToBoson(ComplementarySubsytemSpace.StateDescriptionUp[MinIndex], ComplementarySubsytemSpace.StateDescriptionDown[MinIndex],  
						ComplementarySubsytemSpace.TemporaryStateUp, ComplementarySubsytemSpace.TemporaryStateDown);
      double ComplementarySubsytemSpaceFactorial = 0.0;
      for (int k = 0; k <= ComplementarySubsytemSpace.LzMax; ++k)
	{
	  ComplementarySubsytemSpaceFactorial += LogFactorials[ComplementarySubsytemSpace.TemporaryStateUp[k]];
	  ComplementarySubsytemSpaceFactorial += LogFactorials[ComplementarySubsytemSpace.TemporaryStateDown[k]];
	}
      for (int j = 0; j < SubsytemSpace.GetHilbertSpaceDimension(); ++j)
	{
	  unsigned long* TmpMonomialUp2 = TmpSubsytemSpaceMonomialUp[j];
	  unsigned long* TmpMonomialDown2 = TmpSubsytemSpaceMonomialDown[j];
	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementarySubsytemSpace.NbrBosonsUp) && (TmpIndex3 < SubsytemSpace.NbrBosonsUp)) 
	    {
	      while ((TmpIndex2 < ComplementarySubsytemSpace.NbrBosonsUp) && (TmpMonomialUp2[TmpIndex3] <= TmpMonomialUp1[TmpIndex2]))
		{
		  TmpMonomialUp3[TmpIndex4] = TmpMonomialUp1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementarySubsytemSpace.NbrBosonsUp)
		{
		  while ((TmpIndex3 < SubsytemSpace.NbrBosonsUp) && (TmpMonomialUp1[TmpIndex2] <= TmpMonomialUp2[TmpIndex3]))
		    {
		      TmpMonomialUp3[TmpIndex4] = TmpMonomialUp2[TmpIndex3];
		      ++TmpIndex3;
		      ++TmpIndex4;		  
		    }
		}
	    }
	  while (TmpIndex2 < ComplementarySubsytemSpace.NbrBosonsUp)
	    {
	      TmpMonomialUp3[TmpIndex4] = TmpMonomialUp1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  while (TmpIndex3 < SubsytemSpace.NbrBosonsUp)
	    {
	      TmpMonomialUp3[TmpIndex4] = TmpMonomialUp2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }
	  TmpIndex2 = 0;
	  TmpIndex3 = 0;
	  TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementarySubsytemSpace.NbrBosonsDown) && (TmpIndex3 < SubsytemSpace.NbrBosonsDown)) 
	    {
	      while ((TmpIndex2 < ComplementarySubsytemSpace.NbrBosonsDown) && (TmpMonomialDown2[TmpIndex3] <= TmpMonomialDown1[TmpIndex2]))
		{
		  TmpMonomialDown3[TmpIndex4] = TmpMonomialDown1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementarySubsytemSpace.NbrBosonsDown)
		{
		  while ((TmpIndex3 < SubsytemSpace.NbrBosonsDown) && (TmpMonomialDown1[TmpIndex2] <= TmpMonomialDown2[TmpIndex3]))
		    {
		      TmpMonomialDown3[TmpIndex4] = TmpMonomialDown2[TmpIndex3];
		      ++TmpIndex3;
		      ++TmpIndex4;		  
		    }
		}
	    }
	  while (TmpIndex2 < ComplementarySubsytemSpace.NbrBosonsDown)
	    {
	      TmpMonomialDown3[TmpIndex4] = TmpMonomialDown1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  while (TmpIndex3 < SubsytemSpace.NbrBosonsDown)
	    {
	      TmpMonomialDown3[TmpIndex4] = TmpMonomialDown2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }

	  unsigned long TmpStateUp;
	  unsigned long TmpStateDown;
	  this->ConvertFromMonomial(TmpMonomialUp3, TmpMonomialDown3, TmpStateUp, TmpStateDown);
	  double TmpCoefficient = 1.0;
	  int TmpPos =  this->SymmetrizeAdAdResult(TmpStateUp, TmpStateDown, TmpCoefficient);  
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpStateUp, TmpStateDown, this->TemporaryStateUp, this->TemporaryStateDown);
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->LzMax; ++k)
		{
		  TmpFactorial += LogFactorials[this->TemporaryStateUp[k]];
		  TmpFactorial += LogFactorials[this->TemporaryStateDown[k]];
		}
	      TmpFactorial -= ComplementarySubsytemSpaceFactorial + TmpSubsytemLogFactorials[j] + TmpLogBinomial;
	      TmpFactorial *= 0.5; 	      
	      ++TmpNbrNonZeroElements;
	      double Tmp = TmpCoefficient * exp(TmpFactorial) * groundState[TmpPos];
	      TmpEntanglementMatrix.SetMatrixElement(j, MinIndex, Tmp);
	    }
	}
    }

  for (int i = 0; i < SubsytemSpace.GetHilbertSpaceDimension(); ++i)
    {
      delete[] TmpSubsytemSpaceOccupationNumbersUp[i];
      delete[] TmpSubsytemSpaceOccupationNumbersDown[i];
      delete[] TmpSubsytemSpaceMonomialUp[i];
      delete[] TmpSubsytemSpaceMonomialDown[i];
    }
  delete[] TmpSubsytemSpaceOccupationNumbersUp;
  delete[] TmpSubsytemSpaceOccupationNumbersDown;
  delete[] TmpSubsytemSpaceMonomialUp;
  delete[] TmpSubsytemSpaceMonomialDown;
  delete[] LogFactorials;
  if (TmpNbrNonZeroElements > 0l)
    {
      return TmpEntanglementMatrix;
    }
  else
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
}
   
// evaluate a entanglement matrix of a subsystem of the whole system described by a given ground state, using a generic real space partition. 
// The entanglement matrix is only evaluated in a given Lz sector and computed from precalculated particle entanglement matrix
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// szSector = Sz sector in which the density matrix has to be evaluated 
// nbrOrbitalA = number of orbitals that have to be kept for the A part
// weightOrbitalAUp = weight of each orbital in the A part with spin up (starting from the leftmost orbital)
// weightOrbitalADown = weight of each orbital in the A part with spin down (starting from the leftmost orbital)
// nbrOrbitalB = number of orbitals that have to be kept for the B part
// weightOrbitalBUp = weight of each orbital in the B part with spin up (starting from the leftmost orbital)
// weightOrbitalBDown = weight of each orbital in the B part with spin down (starting from the leftmost orbital)
// entanglementMatrix = reference on the entanglement matrix (will be overwritten)
// return value = reference on the entanglement matrix

RealMatrix& BosonOnSphereWithSU2SpinSzSymmetry::EvaluateEntanglementMatrixGenericRealSpacePartitionFromParticleEntanglementMatrix (int nbrParticleSector, int lzSector, int szSector,
																   int nbrOrbitalA, double* weightOrbitalAUp, double* weightOrbitalADown, 
																   int nbrOrbitalB, double* weightOrbitalBUp, double* weightOrbitalBDown, RealMatrix& entanglementMatrix)
{
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementarySzSector = this->TotalSpin - szSector;
  int TotalLzDisk = ConvertLzFromSphereToDisk(this->TotalLz, this->NbrBosons, this->LzMax);
  int LzADisk = ConvertLzFromSphereToDisk(lzSector, nbrParticleSector, nbrOrbitalA - 1);
  if ((LzADisk < 0) || (LzADisk > ((nbrOrbitalA - 1) * nbrParticleSector)))
    {
      return entanglementMatrix;	  
    }
  int LzBDisk = (TotalLzDisk - LzADisk) - ComplementaryNbrParticles * (this->LzMax + 1 - nbrOrbitalB);
  if ((LzBDisk < 0) || (LzBDisk > ((nbrOrbitalB - 1) * ComplementaryNbrParticles)))
    {
      return entanglementMatrix;	  
    }
  int ComplementaryLzSector = ConvertLzFromDiskToSphere(LzBDisk, ComplementaryNbrParticles, nbrOrbitalB - 1);
  BosonOnSphereWithSU2Spin SubsytemSpace(nbrParticleSector, lzSector, nbrOrbitalA - 1, szSector);
  BosonOnSphereWithSU2Spin ComplementarySubsytemSpace(ComplementaryNbrParticles, ComplementaryLzSector, nbrOrbitalB - 1, ComplementarySzSector);  

  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.GetHilbertSpaceDimension() << endl;
  unsigned long* TmpMonomialUp1 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpMonomialDown1 = new unsigned long [ComplementaryNbrParticles];
  unsigned long* TmpMonomialUp3 = new unsigned long [nbrParticleSector];
  unsigned long* TmpMonomialDown3 = new unsigned long [nbrParticleSector];
  
  for (int i = 0; i < SubsytemSpace.GetHilbertSpaceDimension(); ++i)
    {
      SubsytemSpace.ConvertToMonomial(SubsytemSpace.StateDescriptionUp[i], SubsytemSpace.StateDescriptionDown[i], TmpMonomialUp3, TmpMonomialDown3);
      double Tmp = 1.0;
      for (int j = 0; j < SubsytemSpace.NbrBosonsUp; j++)
	{
	  Tmp *= weightOrbitalAUp[TmpMonomialUp3[j]];
	}
      for (int j = 0; j < SubsytemSpace.NbrBosonsDown; j++)
	{
	  Tmp *= weightOrbitalADown[TmpMonomialDown3[j]];
	}
      for (int j = 0; j < ComplementarySubsytemSpace.GetHilbertSpaceDimension(); ++j)          
	entanglementMatrix(i, j) *= Tmp;      
    }
  
  for (int MinIndex = 0; MinIndex < ComplementarySubsytemSpace.GetHilbertSpaceDimension(); ++MinIndex)    
    {
      ComplementarySubsytemSpace.ConvertToMonomial(ComplementarySubsytemSpace.StateDescriptionUp[MinIndex], ComplementarySubsytemSpace.StateDescriptionDown[MinIndex], TmpMonomialUp1, TmpMonomialDown1);
      double FormFactor = 1.0;
      for (int i = 0; i < ComplementarySubsytemSpace.NbrBosonsUp; i++)
	FormFactor *= weightOrbitalBUp[TmpMonomialUp1[i]];
      for (int i = 0; i < ComplementarySubsytemSpace.NbrBosonsDown; i++)
	FormFactor *= weightOrbitalBDown[TmpMonomialDown1[i]];
      for (int j = 0; j < SubsytemSpace.GetHilbertSpaceDimension(); ++j)
	entanglementMatrix(j, MinIndex) *= FormFactor; 
    }
  
  delete[] TmpMonomialUp1;
  delete[] TmpMonomialDown1;
  delete[] TmpMonomialUp3;
  delete[] TmpMonomialDown3;
  
  return entanglementMatrix;
}
