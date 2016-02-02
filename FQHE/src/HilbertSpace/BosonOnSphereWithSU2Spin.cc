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


#include "config.h"
#include "HilbertSpace/BosonOnSphereWithSU2Spin.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"

#include <math.h>
#include <cstdlib>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
// 

BosonOnSphereWithSU2Spin::BosonOnSphereWithSU2Spin ()
{
}

// basic constructor without any constraint on Sz
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson

BosonOnSphereWithSU2Spin::BosonOnSphereWithSU2Spin (int nbrBosons, int totalLz, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
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
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2Spin!" << endl;
      exit(1);
    } 
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->TargetSpace = this;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  



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

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson
// totalSpin = twice the total spin value
// memory = amount of memory granted for precalculations

BosonOnSphereWithSU2Spin::BosonOnSphereWithSU2Spin (int nbrBosons, int totalLz, int lzMax, int totalSpin, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = totalSpin;
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
//   for (int i = 0; i < TmpHilbertSpaceDimension; ++i)
//     cout << i << " : " << hex << this->StateDescriptionUp[i] << " " << this->StateDescriptionDown[i] << dec << endl;
  if (TmpHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU2Spin!" << endl;
      exit(1);
    }
  this->LargeHilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  

  this->TargetSpace = this;

  this->GenerateLookUpTable(memory);
  //   for (int i = 0; i < this->HilbertSpaceDimension; ++i)	
  //      {
  //        cout << i << " : ";
  //        this->PrintState(cout, i);
  //        cout << this->FindStateIndex(this->StateDescriptionUp[i], this->StateDescriptionDown[i]);
  //        cout << endl;
//        unsigned long Tmp1;
//        unsigned long Tmp2;
//        this->FermionToBoson(this->StateDescriptionUp[i], this->NUpLzMax, this->TemporaryStateUp);
//        this->FermionToBoson(this->StateDescriptionDown[i], this->NDownLzMax, this->TemporaryStateDown);
//        this->BosonToFermion(this->TemporaryStateUp, this->TemporaryStateDown, Tmp1, Tmp2);
//        cout << this->StateDescriptionUp[i] << " " << Tmp1 << " & " << this->StateDescriptionDown[i] << " " << Tmp2 << endl;
//        for (int j = 0; j < this->NbrLzValue; ++j)
// 	 {
// 	   cout << " (" << this->TemporaryStateUp[j] << "," << this->TemporaryStateDown[j] << ")";
// 	 }
//        cout << endl;
//     }
  
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

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSphereWithSU2Spin::BosonOnSphereWithSU2Spin(const BosonOnSphereWithSU2Spin& bosons)
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

BosonOnSphereWithSU2Spin::~BosonOnSphereWithSU2Spin ()
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
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereWithSU2Spin& BosonOnSphereWithSU2Spin::operator = (const BosonOnSphereWithSU2Spin& bosons)
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

AbstractHilbertSpace* BosonOnSphereWithSU2Spin::Clone()
{
  return new BosonOnSphereWithSU2Spin(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereWithSU2Spin::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereWithSU2Spin::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereWithSU2Spin::ExtractSubspace (AbstractQuantumNumber& q, 
								 SubspaceSpaceConverter& converter)
{
  return 0;
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnSphereWithSU2Spin::SetTargetSpace(ParticleOnSphereWithSpin* targetSpace)
{
  this->TargetSpace = (BosonOnSphereWithSU2Spin*) targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int BosonOnSphereWithSU2Spin::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// apply a^+_m_u a_m_u operator to a given state  (only spin up isospin plus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_um a_m_um

double  BosonOnSphereWithSU2Spin::AduAu (int index, int m)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  return (double) (this->TemporaryStateUp[m]);  
}

// apply a^+_m_d a_m_d operator to a given state (only spin down isospin plus)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_dm a_m_dm

double BosonOnSphereWithSU2Spin::AddAd (int index, int m)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  return (double) (this->TemporaryStateDown[m]);  
}

// find state index
//
// stateDescriptionUp = unsigned integer describing the fermionic state for type up particles
// stateDescriptionDown = unsigned integer describing the fermionic state for type down particles
// return value = corresponding index

int BosonOnSphereWithSU2Spin::FindStateIndex(unsigned long stateDescriptionUp, unsigned long stateDescriptionDown)
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
    PosMid = PosMax;

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
  if (CurrentState != stateDescriptionDown)
    return PosMax;
  else
    return PosMid;
}



// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereWithSU2Spin::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescriptionUp[state], this->StateDescriptionDown[state],
		       this->TemporaryStateUp, this->TemporaryStateDown); 

  unsigned long Tmp;
  Str << " | ";
  for (int i = this->LzMax; i >=0 ; --i)
    {
      Str << "(" << this->TemporaryStateUp[i] << "," << this->TemporaryStateDown[i] << ") | ";
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMaxUp = momentum maximum value for a boson in the state with up
// lzMaxDown = momentum maximum value for a boson in the state with down
// totalLz = momentum total value
// nbrNUp = number of particles with quantum up
// nbrNDown = number of particles with quantum down
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereWithSU2Spin::GenerateStates(int nbrBosons, int lzMaxUp, int lzMaxDown, int totalLz, 
					      int nbrNUp, int nbrNDown, long pos)
{
  if ((nbrBosons < 0) || (totalLz < 0) || (nbrNUp < 0) || (nbrNDown < 0))
    return pos;
  if ((nbrBosons == 0) && (totalLz == 0))
    {
      this->StateDescriptionUp[pos] = 0x0ul;
      this->StateDescriptionDown[pos] = 0x0ul;
      return (pos + 1l);
    }
  if ((lzMaxUp < 0) || (lzMaxDown < 0))
    return pos;

  long TmpPos;
  unsigned long MaskUp;
  unsigned long MaskDown;

  if (nbrNUp == 0)
    {
      for (int l = nbrNDown; l > 0; --l)
	{
	  TmpPos = this->GenerateStates(nbrBosons - l, 0, lzMaxDown - 1, totalLz - (lzMaxDown * l), 
					0, nbrNDown - l, pos); 
	  MaskDown = ((0x1ul << l) - 1ul) << (lzMaxDown + nbrNDown - l);
	  for (; pos < TmpPos; ++pos)
	    {
	      this->StateDescriptionDown[pos] |= MaskDown;
	    }
	}
      pos = this->GenerateStates(nbrBosons, 0, lzMaxDown - 1, totalLz, 0, nbrNDown, pos);
      return pos;
    }
  
  TmpPos = this->GenerateStates(nbrBosons - nbrNUp, 0, lzMaxDown, totalLz - (lzMaxUp * nbrNUp), 0, nbrNDown, pos); 
  MaskUp = ((0x1ul << nbrNUp) - 1ul) << lzMaxUp;
  for (; pos < TmpPos; ++pos)
    this->StateDescriptionUp[pos] |= MaskUp;
  for (int i = nbrNUp - 1; i > 0; --i)
    {
      TmpPos = this->GenerateStates(nbrBosons - i, lzMaxUp - 1, lzMaxDown, totalLz - (lzMaxUp * i), nbrNUp - i, nbrNDown, pos); 
      MaskUp = ((0x1ul << i) - 1ul) << (lzMaxUp + nbrNUp - i);
      for (; pos < TmpPos; ++pos)
	{
	  this->StateDescriptionUp[pos] |= MaskUp;
	}
    }
  pos = this->GenerateStates(nbrBosons, lzMaxUp - 1, lzMaxDown, totalLz, nbrNUp, nbrNDown, pos);
  return pos;
};


// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// currentFermionicPositionUp = current fermionic position within the state description for the spin up
// currentFermionicPositionDown = current fermionic position within the state description for the spin down
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereWithSU2Spin::GenerateStates(int nbrBosons, int lzMax, int totalLz, 
					      int currentFermionicPositionUp, int currentFermionicPositionDown, long pos)
{
  if ((nbrBosons < 0) || (totalLz < 0))
    return pos;
  if ((nbrBosons == 0) && (totalLz == 0))
    {
      this->StateDescriptionUp[pos] = 0x0ul;
      this->StateDescriptionDown[pos] = 0x0ul;
      return (pos + 1l);
    }
  if (lzMax < 0)
    return pos;

  for (int i = nbrBosons; i >= 0; --i)
    {
      unsigned long MaskUp = ((0x1ul << i) - 0x1ul)  << (currentFermionicPositionUp - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
	  long TmpPos = this->GenerateStates(nbrBosons - i - j, lzMax - 1, totalLz - (lzMax * (i + j)), 
					     currentFermionicPositionUp - i - 1, currentFermionicPositionDown - j - 1, pos); 
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


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnSphereWithSU2Spin::GenerateLookUpTable(unsigned long memory)
{  
  long TmpUniquePartition = 1l;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->StateDescriptionUp[i - 1] == this->StateDescriptionUp[i]))
	{
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	++TmpUniquePartition;
    }

  this->NbrUniqueStateDescriptionUp = TmpUniquePartition;
  this->UniqueStateDescriptionUp = new unsigned long [this->NbrUniqueStateDescriptionUp];
  this->UniqueStateDescriptionSubArraySizeUp = new int [this->NbrUniqueStateDescriptionUp];
  this->FirstIndexUniqueStateDescriptionUp = new int [this->NbrUniqueStateDescriptionUp];
  TmpUniquePartition = 0l;
  this->UniqueStateDescriptionUp[0l] = this->StateDescriptionUp[0l];
  this->UniqueStateDescriptionSubArraySizeUp[0] = 1;
  this->FirstIndexUniqueStateDescriptionUp[0] = 0;
  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
    {
      while ((i < this->LargeHilbertSpaceDimension) && (this->StateDescriptionUp[i - 1] == this->StateDescriptionUp[i]))
	{
	  ++this->UniqueStateDescriptionSubArraySizeUp[TmpUniquePartition];
	  ++i;
	}
      if (i < this->LargeHilbertSpaceDimension)
	{
	  ++TmpUniquePartition;
	  this->UniqueStateDescriptionUp[TmpUniquePartition] = this->StateDescriptionUp[i];
	  this->UniqueStateDescriptionSubArraySizeUp[TmpUniquePartition] = 1; 
	  this->FirstIndexUniqueStateDescriptionUp[TmpUniquePartition] = i;
	}
    }
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// nbrNUp = number of particles with quantum number up
// nbrNDown = number of particles with quantum number down
// return value = Hilbert space dimension

long BosonOnSphereWithSU2Spin::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, 
								    int nbrNUp, int nbrNDown)
{
  if ((nbrBosons < 0) || (totalLz < 0) || (nbrNUp < 0) || (nbrNDown < 0))
    return 0l;
  if ((nbrBosons == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0)
    return 0l;
  if (nbrBosons == 1)
    {
      if (lzMax >= totalLz)
	return 1l;
      else
	return 0l;
    }
  long Tmp = 0l;
  for (int i = nbrNUp; i >= 0; --i)
    for (int j = nbrNDown; j >= 0; --j)
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - (i + j), lzMax - 1, totalLz - (lzMax * (i + j)), 
							nbrNUp - i, nbrNDown - j);
  return  Tmp;
}

// evaluate Hilbert space dimension without the Sz constraint
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension

long BosonOnSphereWithSU2Spin::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  if ((nbrBosons < 0) || (totalLz < 0))
    return 0l;
  if ((nbrBosons == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0)
    return 0l;
  if (nbrBosons == 1)
    {
      if (lzMax >= totalLz)
	return 2l;
      else
	return 0l;
    }
  long Tmp = 0l;
  for (int i = nbrBosons; i >= 0; --i)
    Tmp += (i + 1) * this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - i, lzMax - 1, totalLz - (lzMax * i));
  return Tmp;
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2Spin::AduAu (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  if (this->TemporaryStateUp[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUp[n];
  --this->TemporaryStateUp[n];
  ++this->TemporaryStateUp[m];
  coefficient *= (double) this->TemporaryStateUp[m];
  coefficient = sqrt(coefficient);  
  return this->TargetSpace->FindStateIndex(this->BosonToFermion(this->TemporaryStateUp), this->StateDescriptionDown[index]);  
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2Spin::AduAd (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  
  if (this->TemporaryStateDown[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  /*cout<<"Index = "<<index<<" m = "<<m << " n " <<n<<endl;
  cout<<" Up = ";
  for (int i = 0; i <=  this->LzMax; i++ )
    {
      cout <<this->TemporaryStateUp[i]<<" ";
    }
  cout <<endl;
  cout <<"Down = ";
  for (int i = 0; i <=  this->LzMax; i++ )
    {
      cout <<this->TemporaryStateDown[i]<<" ";
    }
  cout <<endl;*/
  coefficient = (double) this->TemporaryStateDown[n];
  --this->TemporaryStateDown[n];
  ++this->TemporaryStateUp[m];
  coefficient *= (double) this->TemporaryStateUp[m];
  coefficient = sqrt(coefficient); 
  return this->TargetSpace->FindStateIndex(this->BosonToFermion(this->TemporaryStateUp), this->BosonToFermion(this->TemporaryStateDown));
  /*unsigned long Tmp = this->FindStateIndex(this->BosonToFermion(this->TemporaryStateUp), this->BosonToFermion(this->TemporaryStateDown));
  cout <<"Tmp = "<<Tmp<<endl;
  cout<<" Up = ";
  for (int i = 0; i <=  this->LzMax; i++ )
    {
      cout <<this->TemporaryStateUp[i]<<" ";
    }
  cout <<endl;
  cout <<"Down = ";
  for (int i = 0; i <=  this->LzMax; i++ )
    {
      cout <<this->TemporaryStateDown[i]<<" ";
    } 
  cout <<endl;
  cout <<"end"<<endl;
  
  return Tmp;  */
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2Spin::AddAu (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->TemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  if (this->TemporaryStateUp[n] == 0)
    { 
      coefficient = 0.0;
      return this->TargetSpace->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateUp[n];
  --this->TemporaryStateUp[n];
  ++this->TemporaryStateDown[m];
  coefficient *= (double) this->TemporaryStateDown[m];
  coefficient = sqrt(coefficient);  
  return this->TargetSpace->FindStateIndex(this->TargetSpace->BosonToFermion(this->TemporaryStateUp), this->TargetSpace->BosonToFermion(this->TemporaryStateDown));  
}

// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereWithSU2Spin::AddAd (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->TemporaryStateDown);
  if (this->TemporaryStateDown[n] == 0)
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryStateDown[n];
  --this->TemporaryStateDown[n];
  ++this->TemporaryStateDown[m];
  coefficient *= (double) this->TemporaryStateDown[m];
  coefficient = sqrt(coefficient);  
  return this->TargetSpace->FindStateIndex(this->StateDescriptionUp[index], this->BosonToFermion(this->TemporaryStateDown));  
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU2Spin::AuAu (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  if ((this->ProdATemporaryStateUp[n1] == 0) || (this->ProdATemporaryStateUp[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryStateUp[n1] == 1)))
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
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

double BosonOnSphereWithSU2Spin::AuAd (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateUp[n1] == 0) || (this->ProdATemporaryStateDown[n2] == 0))
    {
      return 0.0;
    }
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

double BosonOnSphereWithSU2Spin::AdAd (int index, int n1, int n2)
{
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
  if ((this->ProdATemporaryStateDown[n1] == 0) || (this->ProdATemporaryStateDown[n2] == 0)
      || ((n1 == n2) && (this->ProdATemporaryStateDown[n1] == 1)))    
    {
      return 0.0;
    }
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  double Coefficient = this->ProdATemporaryStateDown[n2];
  --this->ProdATemporaryStateDown[n2];
  Coefficient *= this->ProdATemporaryStateDown[n1];
  --this->ProdATemporaryStateDown[n1];
  return sqrt(Coefficient);
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU2Spin::AduAdu (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUp, this->TemporaryStateUp, coefficient);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU2Spin::AduAdd (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateUp, this->TemporaryStateDown, coefficient);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU2Spin::AddAdd (int m1, int m2, double& coefficient)
{
  return this->AdiAdj(m1, m2, this->TemporaryStateDown, this->TemporaryStateDown, coefficient);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphereWithSU2Spin::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  double Coefficient = 1.0;
  this->FermionToBoson(this->StateDescriptionUp[index], this->NUpLzMax, this->ProdATemporaryStateUp);
  this->FermionToBoson(this->StateDescriptionDown[index], this->NDownLzMax, this->ProdATemporaryStateDown);
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

int BosonOnSphereWithSU2Spin::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
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
  return this->FindStateIndex(this->TemporaryStateUp, this->TemporaryStateDown);
}

// convert a state from one SU(2) basis to another, transforming the one body basis in each momentum sector
//
// initialState = state to transform  
// targetState = vector where the transformed state has to be stored
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// firstComponent = index of the first component to compute in initialState
// nbrComponents = number of consecutive components to compute

void BosonOnSphereWithSU2Spin::TransformOneBodyBasis(ComplexVector& initialState, ComplexVector& targetState, ComplexMatrix* oneBodyBasis, long firstComponent, long nbrComponents)
{
  int* TmpMomentumIndices = new int [this->NbrBosons];
  int* TmpSU2Indices = new int [this->NbrBosons];
  int* TmpSU2Indices2 = new int [this->NbrBosons];
  double* OccupationCoefficientArray = new double [this->NbrBosons + 1];
  OccupationCoefficientArray[0] = 0.0;
  for (int i = 1; i <= this->NbrBosons; ++i)
    OccupationCoefficientArray[i] = OccupationCoefficientArray[i - 1] + 0.5 * log((double) i);
  targetState.ClearVector();
  long LastComponent = firstComponent + nbrComponents;
  if (nbrComponents == 0)
    LastComponent = this->LargeHilbertSpaceDimension;
  for (long i = firstComponent; i < LastComponent; ++i)
    {
      this->FermionToBoson(this->StateDescriptionUp[i], this->StateDescriptionDown[i],
			   this->TemporaryStateUp, this->TemporaryStateDown); 
      double OccupationCoefficient = 0.0;
      int TmpIndex = 0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  for (int l = 0; l < this->TemporaryStateDown[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 1;
	      ++TmpIndex;	      
	    }
	  for (int l = 0; l < this->TemporaryStateUp[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 0;
	      ++TmpIndex;	      
	    }
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateUp[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateDown[j]];
	}
      this->TransformOneBodyBasisRecursive(targetState, initialState[i], 0, TmpMomentumIndices, TmpSU2Indices, TmpSU2Indices2, oneBodyBasis, OccupationCoefficient, OccupationCoefficientArray);
    }
  delete[] OccupationCoefficientArray;
  delete[] TmpMomentumIndices;
  delete[] TmpSU2Indices;
  delete[] TmpSU2Indices2;
}

// compute the transformation matrix from one SU(2) basis to another, transforming the one body basis in each momentum sector
//
// oneBodyBasis = array that gives the unitary matrices associated to each one body transformation, one per momentum sector
// return value = transformation matrix

ComplexMatrix BosonOnSphereWithSU2Spin::TransformationMatrixOneBodyBasis(ComplexMatrix* oneBodyBasis)
{
  int* TmpMomentumIndices = new int [this->NbrBosons];
  int* TmpSU2Indices = new int [this->NbrBosons];
  int* TmpSU2Indices2 = new int [this->NbrBosons];
  ComplexMatrix TmpMatrix(this->HilbertSpaceDimension, this->HilbertSpaceDimension, true);
  double* OccupationCoefficientArray = new double [this->NbrBosons + 1];
  OccupationCoefficientArray[0] = 0.0;
  for (int i = 1; i <= this->NbrBosons; ++i)
    OccupationCoefficientArray[i] = OccupationCoefficientArray[i - 1] + 0.5 * log((double) i);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->FermionToBoson(this->StateDescriptionUp[i], this->StateDescriptionDown[i], 
			   this->TemporaryStateUp, this->TemporaryStateDown); 
      int TmpIndex = 0;
      double OccupationCoefficient = 0.0;
      for (int j = this->LzMax; j >= 0; --j)
	{
	  for (int l = 0; l < this->TemporaryStateDown[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 1;
	      ++TmpIndex;	      
	    }
	  for (int l = 0; l < this->TemporaryStateUp[j]; ++l)
	    {
	      TmpMomentumIndices[TmpIndex] = j;
	      TmpSU2Indices[TmpIndex] = 0;
	      ++TmpIndex;	      
	    }
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateUp[j]];
	  OccupationCoefficient -= OccupationCoefficientArray[this->TemporaryStateDown[j]];
	}
      this->TransformOneBodyBasisRecursive(TmpMatrix[i], 1.0, 0, TmpMomentumIndices, TmpSU2Indices, TmpSU2Indices2, oneBodyBasis,
					   OccupationCoefficient, OccupationCoefficientArray);
    }
  delete[] TmpMomentumIndices;
  delete[] TmpSU2Indices;
  delete[] TmpSU2Indices2;
  delete[] OccupationCoefficientArray;
  return TmpMatrix;
}

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

void BosonOnSphereWithSU2Spin::TransformOneBodyBasisRecursive(ComplexVector& targetState, Complex coefficient,
							      int position, int* momentumIndices, int* initialSU2Indices, int* currentSU2Indices, ComplexMatrix* oneBodyBasis,
							      double occupationCoefficient, double* occupationCoefficientArray) 
{
  if (position == this->NbrBosons)
    {
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  this->TemporaryStateUp[i] = 0ul;
	  this->TemporaryStateDown[i] = 0ul;
	}
      for (int i = 0; i < this->NbrBosons; ++i)
	{
	  switch (currentSU2Indices[i])
	    {
	    case 0:
	      this->TemporaryStateUp[momentumIndices[i]]++;
	      break;
	    case 1:
	      this->TemporaryStateDown[momentumIndices[i]]++;
	      break;
	    }
	}
      int Index = this->FindStateIndex(this->TemporaryStateUp, this->TemporaryStateDown);
      if (Index < this->HilbertSpaceDimension)
	{
	  
	  for (int i = 0; i <= this->LzMax; ++i)
	    {
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateUp[i]];
	      occupationCoefficient += occupationCoefficientArray[this->TemporaryStateDown[i]];
	    }
	  targetState[Index] += coefficient * exp (occupationCoefficient);
	}
      return;      
    }
  else
    {
      currentSU2Indices[position] = 0;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][1 - initialSU2Indices[position]][1]), position + 1, momentumIndices, initialSU2Indices, currentSU2Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
      currentSU2Indices[position] = 1;
      this->TransformOneBodyBasisRecursive(targetState, coefficient * (oneBodyBasis[momentumIndices[position]][1 - initialSU2Indices[position]][0]), position + 1, momentumIndices, initialSU2Indices, currentSU2Indices, oneBodyBasis, occupationCoefficient, occupationCoefficientArray);
    }
}

// compute the projection matrix from the SU(2) Hilbert space to an U(1) Hilbert space
// 
// targetSpace = pointer to the U(1) Hilbert space
// type = type of particles that has to be kept (0 for type up, 1 for type down)
// return value = projection matrix

ComplexMatrix BosonOnSphereWithSU2Spin::TransformationMatrixSU2ToU1(BosonOnSphereShort* targetSpace, int type)
{
  ComplexMatrix TmpMatrix (targetSpace->HilbertSpaceDimension, this->HilbertSpaceDimension, true);
  unsigned long* TmpStateDescription;
  unsigned long* TmpStateDescriptionOther;
  if (type == 0)
    {
      TmpStateDescription = this->StateDescriptionUp;
      TmpStateDescriptionOther = this->StateDescriptionDown;
    }
  else
    {
      TmpStateDescription = this->StateDescriptionDown;
      TmpStateDescriptionOther = this->StateDescriptionUp;
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (TmpStateDescriptionOther[i] == 0x0ul)
	{
	  unsigned long TmpState = TmpStateDescription[i];
	  int TmpLzMax = this->FermionicLzMax;
	  while ((TmpState >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  int Index = targetSpace->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (Index < targetSpace->HilbertSpaceDimension)
	    {
	      TmpMatrix[i][Index] = 1.0;
	    }
	}
    }
  return TmpMatrix;
}
