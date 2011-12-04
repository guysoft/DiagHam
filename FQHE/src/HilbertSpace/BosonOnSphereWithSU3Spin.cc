////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                    class of bosons on sphere with SU(3) spin               //
//                                                                            //
//                        last modification : 03/12/2011                      //
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
#include "HilbertSpace/BosonOnSphereWithSU3Spin.h"
#include "HilbertSpace/BosonOnSphere.h"
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

BosonOnSphereWithSU3Spin::BosonOnSphereWithSU3Spin ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a boson
// totalTz = twice the total Tz value
// totalY = three time the total Y value
// memory = amount of memory granted for precalculations

BosonOnSphereWithSU3Spin::BosonOnSphereWithSU3Spin (int nbrBosons, int totalLz, int lzMax, int totalTz, int totalY,
						    unsigned long memory)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->TotalY = totalY;
  this->TotalTz = totalTz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
  int N1 = (2 * nbrBosons) + totalY + (3 * totalTz);
  int N2 = (2 * nbrBosons) + totalY - (3 * totalTz);
  int N3 = nbrBosons - totalY;
  if ((N1 < 0) || (N2 < 0) || (N3 < 0) || ((N1 % 6) != 0) || ((N2 % 6) != 0) || ((N3 % 3) != 0))
    this->HilbertSpaceDimension = 0;
  else
    {
      N1 /= 6;
      N2 /= 6;
      N3 /= 3;
      this->HilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, N1, N2, N3);
    }
  this->StateDescription1 = new unsigned long [this->HilbertSpaceDimension];
  this->StateDescription2 = new unsigned long [this->HilbertSpaceDimension];
  this->StateDescription3 = new unsigned long [this->HilbertSpaceDimension];
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->LzMax, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, 
						       N1, N2, N3, 0l);
  if (TmpHilbertSpaceDimension != this->HilbertSpaceDimension)
    {
      cout << TmpHilbertSpaceDimension << " " << this->HilbertSpaceDimension << endl;
      cout << "Mismatch in State-count and State Generation in BosonOnSphereWithSU3Spin!" << endl;
      exit(1);
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    cout << i << " : " << hex << this->StateDescription1[i] << " " << this->StateDescription2[i] << " " << this->StateDescription3[i] << dec << endl;
  this->HilbertSpaceDimension = TmpHilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  
  //  this->GenerateLookUpTable(memory);
//   for (int i = 0; i < this->HilbertSpaceDimension; ++i)	
//     this->PrintState(cout, i) << endl;
// #ifdef __DEBUG__
//   int UsedMemory = 0;
//   UsedMemory += this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
//   cout << "memory requested for Hilbert space = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;
//   UsedMemory = this->NbrLzValue * sizeof(int);
//   UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
//   cout << "memory requested for lookup table = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;

// #endif
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSphereWithSU3Spin::BosonOnSphereWithSU3Spin(const BosonOnSphereWithSU3Spin& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->TotalTz = bosons.TotalTz;
  this->TotalY = bosons.TotalY;
  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;  
  this->SignLookUpTable = bosons.SignLookUpTable;
  this->SignLookUpTableMask = bosons.SignLookUpTableMask;
  this->MaximumSignLookUp = bosons.MaximumSignLookUp;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

BosonOnSphereWithSU3Spin::~BosonOnSphereWithSU3Spin ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
//       if (this->StateHighestBit != 0)
// 	delete[] this->StateHighestBit;
//       delete[] this->LookUpTableShift;
//       for (int i = 0; i < 2*this->NbrLzValue; ++i)
// 	delete[] this->LookUpTable[i];
//       delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereWithSU3Spin& BosonOnSphereWithSU3Spin::operator = (const BosonOnSphereWithSU3Spin& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      //      delete[] this->StateDescription;
      //      delete[] this->StateHighestBit;
    }
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->TotalTz = bosons.TotalTz;
  this->TotalY = bosons.TotalY;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->StateDescription = bosons.StateDescription;
  this->StateHighestBit = bosons.StateHighestBit;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;  
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereWithSU3Spin::Clone()
{
  return new BosonOnSphereWithSU3Spin(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereWithSU3Spin::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereWithSU3Spin::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereWithSU3Spin::ExtractSubspace (AbstractQuantumNumber& q, 
								 SubspaceSpaceConverter& converter)
{
  return 0;
}

// apply a^+_m_1 a_m_1 operator to a given state (only state 1 Tz=+1/2, Y=+1/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_1 a_m_1

double  BosonOnSphereWithSU3Spin::Ad1A1 (int index, int m)
{
  if ((this->StateDescription[index] & (0x1l << (m * 3))) != 0)
    return 1.0;
  else
    return 0.0;
}

// apply a^+_m_2 a_m_2 operator to a given state (only state 2 Tz=-1/2, Y=+1/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_2 a_m_2

double BosonOnSphereWithSU3Spin::Ad2A2 (int index, int m)
{
  if ((this->StateDescription[index] & (0x2l << (m * 3))) != 0)
    return 1.0;
  else
    return 0.0;
}
 
// apply a^+_m_3 a_m_3 operator to a given state (only state 3 Tz=0, Y=-2/3)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m_3 a_m_3

double BosonOnSphereWithSU3Spin::Ad3A3 (int index, int m)
{
  if ((this->StateDescription[index] & (0x4l << (m * 3))) != 0)
    return 1.0;
  else
    return 0.0;
}

// find state index
//
// stateDescription1 = unsigned integer describing the state in the 1 sector
// stateDescription2 = unsigned integer describing the state in the 2 sector
// stateDescription3 = unsigned integer describing the state in the 3 sector
// lzmax = maximum Lz value reached by a boson in the state
// return value = corresponding index

int BosonOnSphereWithSU3Spin::FindStateIndex(unsigned long stateDescription1, unsigned long stateDescription2, unsigned long stateDescription3, int lzmax)
{
  unsigned long CurrentState = stateDescription1 >> this->LookUpTableShift[lzmax];
  int PosMin = this->LookUpTable[lzmax][CurrentState];
  int PosMax = this->LookUpTable[lzmax][CurrentState+ 1];
  int PosMid = (PosMin + PosMax) >> 1;
  CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription1))
    {
      if (CurrentState > stateDescription1)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription1)
    return PosMid;
  else
    return PosMin;
}



// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereWithSU3Spin::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  Str << " | ";
  for (int i = this->NbrLzValue-1; i >=0 ; --i)
    {
      Tmp = ((TmpState >> (i * 3)) & ((unsigned long) 0x7ul));
      switch (Tmp)
	{
	case 0x7ul :
	  Str << "123 ";
	  break;
	case 0x6ul :
	  Str << " 23 ";
	  break;
	case 0x5ul :
	  Str << "1 3 ";
	  break;
	case 0x4ul :
	  Str << "  3 ";
	  break;
	case 0x3ul :
	  Str << "12  ";
	  break;
	case 0x2ul :
	  Str << " 2  ";
	  break;
	case 0x1ul :
	  Str << "1   ";
	  break;
	case 0x0ul :
	  Str << " 0  ";
	  break;
	}
      Str << "| ";
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson in the state
// totalLz = momentum total value
// nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
// nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
// nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereWithSU3Spin::GenerateStates(int nbrBosons, int lzMax, int totalLz, int nbrN1, int nbrN2, int nbrN3, long pos)
{
  if ((nbrBosons < 0) || (totalLz < 0) || (nbrN1 < 0) || (nbrN2 < 0) || (nbrN3 < 0))
    return pos;
  if ((nbrBosons == 0) && (totalLz == 0))
    {
      this->StateDescription1[pos] = 0x0ul;
      this->StateDescription2[pos] = 0x0ul;
      this->StateDescription3[pos] = 0x0ul;
      return (pos + 1l);
    }
  if (lzMax < 0)
    return pos;

  if (nbrBosons == 1) 
    {
      if (lzMax >= totalLz)
	{
	  this->StateDescription1[pos] = ((unsigned long) nbrN1) << totalLz;
	  this->StateDescription2[pos] = ((unsigned long) nbrN2) << totalLz;
	  this->StateDescription3[pos] = ((unsigned long) nbrN3) << totalLz;
	  return (pos + 1l);
	}
      else
	return pos;
    }

  long TmpPos;
  unsigned long Mask1;
  unsigned long Mask2;
  unsigned long Mask3;
  for (int i = nbrN1; i >= 0; --i)
    for (int j = nbrN2; j >= 0; --j)
      for (int k = nbrN3; k >= 0; --k)
	{
	  TmpPos = this->GenerateStates(nbrBosons - (i + j + k), lzMax - 1, totalLz - (lzMax * (i + j + k)), 
					nbrN1 - i, nbrN2 - j, nbrN3 - k, pos); 
	  Mask1 = ((0x1ul << i) - 1ul) << (lzMax + nbrN1 - i);
	  Mask2 = ((0x1ul << i) - 1ul) << (lzMax + nbrN2 - j);
	  Mask3 = ((0x1ul << i) - 1ul) << (lzMax + nbrN3 - k);
	  for (; pos < TmpPos; ++pos)
	    {
	      this->StateDescription1[pos] |= Mask1;
	      this->StateDescription2[pos] |= Mask1;
	      this->StateDescription3[pos] |= Mask1;
	    }
	}
  return pos;
};


// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void BosonOnSphereWithSU3Spin::GenerateLookUpTable(unsigned long memory)
{
  // get every highest bit poisition
  unsigned long TmpPosition = this->StateDescription[0];
  int CurrentHighestBit = (this->LzMax + 1) * 3 - 1;
  while ((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul)
    --CurrentHighestBit;  
  int MaxHighestBit = CurrentHighestBit;
  this->StateHighestBit[0] = CurrentHighestBit;
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      TmpPosition = this->StateDescription[i];
      while ((TmpPosition & (0x1ul << CurrentHighestBit)) == 0x0ul)
	--CurrentHighestBit;  
      this->StateHighestBit[i] = CurrentHighestBit;
   }

  // evaluate look-up table size
  memory /= (sizeof(int*) * MaxHighestBit);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > MaxHighestBit)
    this->MaximumLookUpShift = MaxHighestBit;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new int* [MaxHighestBit + 1];
  this->LookUpTableShift = new int [MaxHighestBit + 1];
  for (int i = 0; i <= MaxHighestBit; ++i)
    this->LookUpTable[i] = new int [this->LookUpTableMemorySize + 1];

  CurrentHighestBit = this->StateHighestBit[0];
  int* TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
  if (CurrentHighestBit < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentHighestBit] = 0;
  else
    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentHighestBit];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = this->StateDescription[0] >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentHighestBit != this->StateHighestBit[i])
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
 	  CurrentHighestBit = this->StateHighestBit[i];
	  TmpLookUpTable = this->LookUpTable[CurrentHighestBit];
	  if (CurrentHighestBit < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentHighestBit] = 0;
	  else
	    this->LookUpTableShift[CurrentHighestBit] = CurrentHighestBit + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentHighestBit];
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  CurrentLookUpTableValue = this->LookUpTableMemorySize;
	  while (CurrentLookUpTableValue > TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = this->StateDescription[i] >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  while (CurrentLookUpTableValue > 0)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->HilbertSpaceDimension - 1;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[0] = this->HilbertSpaceDimension - 1;

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
    this->SignLookUpTableMask[i] = 0xfffful;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = 0xfffful;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = 0xfffful >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = 0x0ul;
#endif
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// nbrN1 = number of particles with quantum number Tz=+1/2 and Y=+1/3
// nbrN2 = number of particles with quantum number Tz=-1/2 and Y=+1/3
// nbrN3 = number of particles with quantum number Tz=0 and Y=-2/3
// return value = Hilbert space dimension

long BosonOnSphereWithSU3Spin::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int nbrN1, int nbrN2, int nbrN3)
{
  if ((nbrBosons < 0) || (totalLz < 0) || (lzMax < 0) || (nbrN1 < 0) || (nbrN2 < 0) || (nbrN3 < 0) || (lzMax < 0))
    return 0l;
  if ((nbrBosons == 0) && (totalLz == 0))
    return 1l;
  if (nbrBosons == 1)
    {
      if (lzMax >= totalLz)
	return 1l;
      else
	return 0l;
    }
  long Tmp = 0l;
  for (int i = nbrN1; i >= 0; --i)
    for (int j = nbrN2; j >= 0; --j)
      for (int k = nbrN3; k >= 0; --k)
	Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - (i + j + k), lzMax - 1, totalLz - (lzMax * (i + j + k)), 
							  nbrN1 - i, nbrN2 - j, nbrN3 - k);
  return  Tmp;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// totalTz = twice the total Tz value
// totalY = three time the total Y value
// return value = Hilbert space dimension

int BosonOnSphereWithSU3Spin::EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz, int totalTz, int totalY)
{
  int N1 = (2 * nbrBosons) + totalY + (3 * totalTz);
  int N2 = (2 * nbrBosons) + totalY - (3 * totalTz);
  int N3 = nbrBosons - totalY;
  if ((N1 >= 0) && (N2 >= 0) && (N3 >= 0) && ((N1 % 6) == 0) && ((N2 % 6) == 0) && ((N3 % 3) == 0))
    {
      N1 /= 6;
      N2 /= 6;
      N3 /= 3;
      return this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax, (totalLz + (nbrBosons * lzMax)) >> 1, N1, N2, N3);
    }
  else
    return 0;
}

// apply a_n1_1 a_n2_1 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU3Spin::A1A1 (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 *= 3;
  n2 *= 3;
 if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_1 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU3Spin::A1A2 (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 *= 3;
  n2 *= 3;
  ++n2;
 if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_1 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU3Spin::A1A3 (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 *= 3;
  n2 *= 3;
  n2 += 2;
 if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_2 a_n2_2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU3Spin::A2A2 (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 *= 3;
  ++n1;
  n2 *= 3;
  ++n2;
 if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_2 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU3Spin::A2A3 (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 *= 3;
  ++n1;
  n2 *= 3;
  n2 += 2;
 if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a_n1_3 a_n2_3 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next Ad*Ad* call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor

double BosonOnSphereWithSU3Spin::A3A3 (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 *= 3;
  n1 += 2;
  n2 *= 3;
  n2 += 2;
 if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0) || (n1 == n2))
    return 0.0;
  this->ProdALzMax = this->StateHighestBit[index];
  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(0x1ul << n1);
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply a^+_m1_1 a^+_m2_1 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU3Spin::Ad1Ad1 (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 *= 3;
  m2 *= 3;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, TmpState, TmpState, NewLzMax);
}

// apply a^+_m1_1 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU3Spin::Ad1Ad2 (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 *= 3;
  m2 *= 3;
  ++m2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, TmpState, TmpState, NewLzMax);
}

// apply a^+_m1_1 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU3Spin::Ad1Ad3 (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 *= 3;
  m2 *= 3;
  m2 += 2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, TmpState, TmpState, NewLzMax);
}

// apply a^+_m1_2 a^+_m2_2 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU3Spin::Ad2Ad2 (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 *= 3;
  ++m1;
  m2 *= 3;
  ++m2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, TmpState, TmpState, NewLzMax);
}

// apply a^+_m1_2 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU3Spin::Ad2Ad3 (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 *= 3;
  ++m1;
  m2 *= 3;
  m2 += 2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, TmpState, TmpState, NewLzMax);
}

// apply a^+_m1_3 a^+_m2_3 operator to the state produced using A*A* method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

int BosonOnSphereWithSU3Spin::Ad3Ad3 (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 *= 3;
  m1 += 2;
  m2 *= 3;
  m2 += 2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  coefficient = 1.0;
  int NewLzMax = this->ProdALzMax;
  if (m2 > NewLzMax)
    NewLzMax = m2;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
    }
  TmpState |= (0x1ul << m2);
  if (m1 > NewLzMax)
    NewLzMax = m1;
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
    }
  TmpState |= (0x1ul << m1);
  return this->FindStateIndex(TmpState, TmpState, TmpState, NewLzMax);
}

