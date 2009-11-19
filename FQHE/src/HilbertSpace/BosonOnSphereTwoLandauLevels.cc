////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of bosons on sphere including two                 //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 09/09/2009                      //
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
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereTwoLandauLevels.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"

#include <cmath>
#include <bitset>
#include <cstdlib>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::bitset;

#define WANT_LAPACK

#ifdef __LAPACK__
#ifdef WANT_LAPACK
#define  __USE_LAPACK_HERE__
#endif
#endif


// default constructor
//

BosonOnSphereTwoLandauLevels::BosonOnSphereTwoLandauLevels()
{
}

// basic constructor with no contraint on the number of particles per spin component
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMaxUp = twice the maximum Lz value reached by a boson with a spin up
// lzMaxDown = twice the maximum Lz value reached by a boson with a spin down
// memory = amount of memory granted for precalculations

BosonOnSphereTwoLandauLevels::BosonOnSphereTwoLandauLevels (int nbrBosons, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  this->NbrBosonsUp = 0;
  this->NbrBosonsDown = 0;
  this->LzMaxUp = lzMaxUp;
  this->LzMaxDown = lzMaxDown;
  if (this->LzMaxUp >= this->LzMaxDown)
    {
      this->LzMax = this->LzMaxUp;
      this->LzShiftUp = 0;
      this->LzShiftDown = (this->LzMaxUp - this->LzMaxDown) >> 1;
    }
  else
    {
      this->LzMax = this->LzMaxDown;
      this->LzShiftDown = 0;
      this->LzShiftUp = (this->LzMaxDown - this->LzMaxUp) >> 1;
    }
  this->UpStateShift = nbrBosons + lzMaxDown;
  this->LzTotalShift = this->LzMaxDown + this->LzMaxUp;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  int TmpDimension = this->GenerateStates(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, 0);
  if (TmpDimension != this->HilbertSpaceDimension)
    {
      cout << "Mismatch in State-count and State Generation in BosonOnSphereTwoLandauLevels! " << this->HilbertSpaceDimension << " " << TmpDimension  << endl;
  for (int i = 0; i < TmpDimension; ++i)
    this->PrintState(cout, i) << endl;
       exit(1);
    }

  //  this->HilbertSpaceDimension = this->GenerateStates(this->NbrBosonsUp, this->NbrBosonsDown, this->LzMaxUp, this->LzMaxDown, );
  this->GenerateLookUpTable(memory);
  
#ifdef __DEBUG__
  long UsedMemory = 0;
  UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;
  UsedMemory = this->NbrLzValue * sizeof(int);
  UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
  cout << "memory requested for lookup table = ";
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

BosonOnSphereTwoLandauLevels::BosonOnSphereTwoLandauLevels(const BosonOnSphereTwoLandauLevels& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->LzMaxUp = bosons.LzMaxUp;
  this->LzMaxDown = bosons.LzMaxDown;
  this->LzShiftUp = bosons.LzShiftUp;
  this->LzShiftDown = bosons.LzShiftDown;
  this->LzTotalShift = bosons.LzTotalShift;
  this->UpStateShift = bosons.UpStateShift;
  this->NbrLzValue = bosons.NbrLzValue;
  this->TotalSpin = bosons.TotalSpin;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;  
  this->SignLookUpTable = bosons.SignLookUpTable;
  this->SignLookUpTableMask = bosons.SignLookUpTableMask;
  this->MaximumSignLookUp = bosons.MaximumSignLookUp;
}

// destructor
//

BosonOnSphereTwoLandauLevels::~BosonOnSphereTwoLandauLevels ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereTwoLandauLevels& BosonOnSphereTwoLandauLevels::operator = (const BosonOnSphereTwoLandauLevels& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->LzMaxUp = bosons.LzMaxUp;
  this->LzMaxDown = bosons.LzMaxDown;
  this->LzShiftUp = bosons.LzShiftUp;
  this->LzShiftDown = bosons.LzShiftDown;
  this->LzTotalShift = bosons.LzTotalShift;
  this->NbrLzValue = bosons.NbrLzValue;
  this->TotalSpin = bosons.TotalSpin;
  this->NbrBosonsUp = bosons.NbrBosonsUp;
  this->NbrBosonsDown = bosons.NbrBosonsDown;
  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereTwoLandauLevels::Clone()
{
  return new BosonOnSphereTwoLandauLevels(*this);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereTwoLandauLevels::ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter)
{
  return 0;
}


// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereTwoLandauLevels::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereTwoLandauLevels::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereTwoLandauLevels::PrintState (ostream& Str, int state)
{
  return Str;
}

// evaluate Hilbert space dimension without constraint on the number of particles per level
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension

long BosonOnSphereTwoLandauLevels::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  if ((nbrBosons < 0) || (totalLz < 0))
    return 0l;
  if ((nbrBosons == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0) 
    return 0l;
    
  if ((lzMax == 0) && (totalLz != 0))
    return 0l;

  long Tmp = 0l;
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      for (int i = nbrBosons; i >= 1; --i)
	Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - i, lzMax - 1, totalLz - (i * lzMax));
      if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
	for (int i = nbrBosons & 0x0fffffffe; i >= 2; i -= 2)
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - i, lzMax - 1, totalLz - (i * lzMax));
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))    
    for (int i = nbrBosons; i >= 1; --i) 
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - i, lzMax - 1, totalLz - (i * lzMax));
  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax - 1, totalLz);
  return Tmp;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereTwoLandauLevels::GenerateStates(int nbrBosons, int lzMax, int totalLz, long pos)
{
  if (nbrBosons < 0)
    return pos;
  if ((nbrBosons == 0) && (totalLz == 0))
    {
      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
    
  if (lzMax < 0)
    return pos;
  
  if ((lzMax == 0) && (totalLz != 0))
    return pos;

  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      long TmpPos = this->GenerateStates(nbrBosons - 1, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = (0x1ul << lzMax) << this->UpStateShift;
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if (((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp)) &&
      ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown)))
    {
      long TmpPos = this->GenerateStates(nbrBosons - 2, lzMax - 1, totalLz - (2 * lzMax), pos);
      unsigned long Mask = ((0x1ul << lzMax) << this->UpStateShift) | (0x1ul << lzMax);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
    {
      long TmpPos = this->GenerateStates(nbrBosons - 1, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x1ul << lzMax;
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }

  return this->GenerateStates(nbrBosons, lzMax - 1, totalLz, pos);
}

// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

RealVector BosonOnSphereTwoLandauLevels::ForgeSU2FromU1(RealVector& upState, BosonOnSphere& upStateSpace, RealVector& downState, BosonOnSphere& downStateSpace)
{
  RealVector FinalState(this->HilbertSpaceDimension, true);
//   for (int j = 0; j < upStateSpace.HilbertSpaceDimension; ++j)
//     {
//       unsigned long TmpUpState = upStateSpace.StateDescription[j] << this->LzShiftUp;
//       int TmpPos = upStateSpace.LzMax + this->LzShiftUp;
//       while (TmpPos > 0)
// 	{
// 	  unsigned long Tmp = TmpUpState & (0x1ul << TmpPos);
// 	  TmpUpState |= Tmp << TmpPos;
// 	  TmpUpState ^= Tmp;
// 	  --TmpPos;
// 	}
//       TmpUpState <<= 1;
//       double TmpComponent = upState[j];
//       int Max = 63;
//       while ((TmpUpState & (0x1ul << Max)) == 0x0ul)
// 	--Max;
//       int Min = 0;
//       while ((TmpUpState & (0x1ul << Min)) == 0x0ul)
// 	++Min;
//       unsigned long TmpUpStateMask = (0x1ul << Max) - 1;
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	if ((this->StateDescription[i] & TmpUpState) == TmpUpState)
// 	  {	    
// 	    unsigned long TmpUpState3 = this->StateDescription[i] & TmpUpStateMask;
// 	    unsigned long TmpUpState2 = TmpUpState3;
// #ifdef  __64_BITS__
// 	    TmpUpState3 &= 0x5555555555555555ul;
// 	    TmpUpState2 &= 0xaaaaaaaaaaaaaaaaul;
// #else
// 	    TmpUpState3 &= 0x55555555ul;
// 	    TmpUpState2 &= 0xaaaaaaaaul;
// #endif	    
// 	    unsigned long Sign = 0x0;
// 	    int Pos = this->LzMax << 1;
// 	    while ((Pos > 0) && ((TmpUpState3 & (0x1ul << Pos)) == 0x0ul))
// 	      Pos -= 2;
// 	    while (Pos > 0)
// 	      {
// 		unsigned long TmpUpState4 = TmpUpState2 & ((0x1ul << Pos) - 1ul);
// #ifdef  __64_BITS__
// 		TmpUpState4 ^= TmpUpState4 >> 32;
// #endif	
// 		TmpUpState4 ^= TmpUpState4 >> 16;
// 		TmpUpState4 ^= TmpUpState4 >> 8;
// 		TmpUpState4 ^= TmpUpState4 >> 4;
// 		TmpUpState4 ^= TmpUpState4 >> 2;
// 		TmpUpState4 ^= TmpUpState4 >> 1;
// 		Sign ^= TmpUpState4;
// 		Pos -= 2;
// 		while ((Pos > 0) && ((TmpUpState3 & (0x1ul << Pos)) == 0x0ul))
// 		  Pos -= 2;
// 	      }
// 	    if ((Sign & 0x1ul) == 0x0ul)
// 	      FinalState[i] = TmpComponent;
// 	    else
// 	      FinalState[i] = -TmpComponent;
// 	  }
//     }

//   for (int j = 0; j < downStateSpace.HilbertSpaceDimension; ++j)
//     {
//       unsigned long TmpDownState = downStateSpace.StateDescription[j] << this->LzShiftDown;
//       int TmpPos = downStateSpace.LzMax + this->LzShiftDown;
//       while (TmpPos > 0)
// 	{
// 	  unsigned long Tmp = TmpDownState & (0x1ul << TmpPos);
// 	  TmpDownState |= Tmp << TmpPos;
// 	  TmpDownState ^= Tmp;
// 	  --TmpPos;
// 	}
//       double TmpComponent = downState[j];
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	if ((this->StateDescription[i] & TmpDownState) == TmpDownState)
// 	  {
// 	    FinalState[i] *= TmpComponent;
// 	  }
//     }

  return FinalState;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnSphereTwoLandauLevels::GenerateLookUpTable(unsigned long memory)
{
}

