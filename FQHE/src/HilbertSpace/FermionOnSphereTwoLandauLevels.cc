////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on sphere including two                //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 19/05/2009                      //
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
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/ArrayTools.h"

#include <cmath>
#include <bitset>
#include <cstdlib>
#include <algorithm>

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

FermionOnSphereTwoLandauLevels::FermionOnSphereTwoLandauLevels()
{
}

// basic constructor with no contraint on the number of particles per spin component
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMaxUp = twice the maximum Lz value reached by a fermion with a spin up
// lzMaxDown = twice the maximum Lz value reached by a fermion with a spin down
// memory = amount of memory granted for precalculations

FermionOnSphereTwoLandauLevels::FermionOnSphereTwoLandauLevels (int nbrFermions, int totalLz, int lzMaxUp, int lzMaxDown, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
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
  this->LzTotalShift = this->LzMaxDown + this->LzMaxUp;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateFullHilbertSpaceDimension(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  int TmpDimension = this->GenerateFullStates(this->NbrFermions, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 0);
  if (TmpDimension != this->HilbertSpaceDimension)
    {
      cout << "Mismatch in State-count and State Generation in FermionOnSphereTwoLandauLevels! " << this->HilbertSpaceDimension << " " << TmpDimension  << endl;
  for (int i = 0; i < TmpDimension; ++i)
    this->PrintState(cout, i) << endl;
       exit(1);
    }

  //  this->HilbertSpaceDimension = this->GenerateStates(this->NbrFermionsUp, this->NbrFermionsDown, this->LzMaxUp, this->LzMaxDown, );
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
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereTwoLandauLevels::FermionOnSphereTwoLandauLevels(const FermionOnSphereTwoLandauLevels& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->LzMaxUp = fermions.LzMaxUp;
  this->LzMaxDown = fermions.LzMaxDown;
  this->LzShiftUp = fermions.LzShiftUp;
  this->LzShiftDown = fermions.LzShiftDown;
  this->LzTotalShift = fermions.LzTotalShift;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
}

// destructor
//

FermionOnSphereTwoLandauLevels::~FermionOnSphereTwoLandauLevels ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereTwoLandauLevels& FermionOnSphereTwoLandauLevels::operator = (const FermionOnSphereTwoLandauLevels& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->LzMaxUp = fermions.LzMaxUp;
  this->LzMaxDown = fermions.LzMaxDown;
  this->LzShiftUp = fermions.LzShiftUp;
  this->LzShiftDown = fermions.LzShiftDown;
  this->LzTotalShift = fermions.LzTotalShift;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereTwoLandauLevels::Clone()
{
  return new FermionOnSphereTwoLandauLevels(*this);
}

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereTwoLandauLevels::AduAu (int index, int m)
{
  if ((this->StateDescription[index] & (0x2l << ((m + this->LzShiftUp) << 1))) != 0)
    return 1.0;
  else
    return 0.0;
}

// apply a^+_d_m a_d_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_d_m a_d_m

double FermionOnSphereTwoLandauLevels::AddAd (int index, int m)
{
  if ((this->StateDescription[index] & (0x1l << ((m + this->LzShiftDown) << 1))) != 0)
    return 1.0;
  else
    return 0.0;
}



// apply a^+_m_u a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereTwoLandauLevels::AduAu (int index, int m, int n, double& coefficient)
{  
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m += this->LzShiftUp;
  n += this->LzShiftUp;
  m = (m<<1) + 1;
  n = (n<<1) + 1;  
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(((unsigned long) (0x1ul)) << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if ((State & (((unsigned long) (0x1ul)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (((unsigned long) (0x1ul)) << m);
  return this->FindStateIndex(State, NewLargestBit);
}

// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereTwoLandauLevels::AddAd (int index, int m, int n, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m += this->LzShiftDown;
  n += this->LzShiftDown;
  m <<= 1;
  n <<= 1;
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(((unsigned long) (0x1ul)) << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if ((State & (((unsigned long) (0x1ul)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (((unsigned long) (0x1ul)) << m);
  return this->FindStateIndex(State, NewLargestBit);
}


// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereTwoLandauLevels::AduAd (int index, int m, int n, double& coefficient)
  {
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m += this->LzShiftUp;
  n += this->LzShiftDown;
  m = (m<<1) + 1;
  n <<= 1;
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(((unsigned long) (0x1ul)) << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if ((State & (((unsigned long) (0x1ul)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (((unsigned long) (0x1ul)) << m);
  return this->FindStateIndex(State, NewLargestBit);
}



// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnSphereTwoLandauLevels::AddAu (int index, int m, int n, double& coefficient)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m += this->LzShiftDown;
  n += this->LzShiftUp;
  m <<= 1;
  n = (n<<1) + 1;  
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(((unsigned long) (0x1ul)) << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if ((State & (((unsigned long) (0x1ul)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLargestBit)
    {
      NewLargestBit = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(State >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(State >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (((unsigned long) (0x1ul)) << m);
  return this->FindStateIndex(State, NewLargestBit);
}


// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 

double FermionOnSphereTwoLandauLevels::AuAu (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 += this->LzShiftUp;
  n2 += this->LzShiftUp;
  n1 <<= 1;
  ++n1;
  n2 <<= 1;
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

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnSphereTwoLandauLevels::AdAd (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 += this->LzShiftDown;
  n2 += this->LzShiftDown;
  n1 <<= 1;
  n2 <<= 1;
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

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnSphereTwoLandauLevels::AuAd (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription[index];
  n1 += this->LzShiftUp;
  n2 += this->LzShiftDown;
  n1 <<= 1;
  ++n1;
  n2 <<= 1;
  if (((this->ProdATemporaryState & (0x1ul << n1)) == 0) || ((this->ProdATemporaryState & (0x1ul << n2)) == 0))
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
  while ( ((this->ProdATemporaryState >> this->ProdALzMax) == 0) && (this->ProdALzMax>0))
    --this->ProdALzMax;
  return Coefficient;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereTwoLandauLevels::AduAdu (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 += this->LzShiftUp;
  m2 += this->LzShiftUp;
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  ++m2;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
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
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereTwoLandauLevels::AddAdd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 += this->LzShiftDown;
  m2 += this->LzShiftDown;
  m1 <<= 1;
  m2 <<= 1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    return this->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
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
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereTwoLandauLevels::AduAdd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  m1 += this->LzShiftUp;
  m2 += this->LzShiftDown;
  m1 <<= 1;
  ++m1;
  m2 <<= 1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0))
    return this->HilbertSpaceDimension;
  int NewLzMax = this->ProdALzMax;
  coefficient = 1.0;
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
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereTwoLandauLevels::ProdA (int index, int* n, int* spinIndices, int nbrIndices)
{
  this->ProdALzMax = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (spinIndices[i] == 0)
	Index = ((n[i] + this->LzShiftDown) << 1) + spinIndices[i];
      else
 	Index = ((n[i] + this->LzShiftUp) << 1) + spinIndices[i];
     if ((this->ProdATemporaryState & (0x1l << Index)) == 0)
	{
	  return 0.0;
	}
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> Index) & this->SignLookUpTableMask[Index]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
#ifdef  __64_BITS__
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      this->ProdATemporaryState &= ~(0x1l << Index);
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;

  return Coefficient;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// spinIndices = integer that gives the spin indices associated to each annihilation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereTwoLandauLevels::ProdA (int index, int* n, int spinIndices, int nbrIndices)
{
  this->ProdALzMax = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (((spinIndices >> i) & 0x1) == 0)
	Index = ((n[i] + this->LzShiftDown) << 1) + ((spinIndices >> i) & 0x1);
      else
	Index = ((n[i] + this->LzShiftUp) << 1) + ((spinIndices >> i) & 0x1);
      if ((this->ProdATemporaryState & (0x1l << Index)) == 0)
	{
	  return 0.0;
	}
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> Index) & this->SignLookUpTableMask[Index]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
#ifdef  __64_BITS__
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      this->ProdATemporaryState &= ~(0x1l << Index);
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;

  return Coefficient;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = array of spin indixes associated to each annihilation operators first index corresponding to the leftmost operator, 0 stands for spin down and 1 stands for spin up)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereTwoLandauLevels::ProdAd (int* m, int* spinIndices, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  unsigned long TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (spinIndices[i] == 0)
	Index = ((m[i] + this->LzShiftDown) << 1) + spinIndices[i];
      else
 	Index = ((m[i] + this->LzShiftUp) << 1) + spinIndices[i];
      if ((TmpState & (0x1l << Index)) != 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
	}
      TmpState |= (0x1l << Index);
    }
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// spinIndices = integer that gives the spin indices associated to each creation operators, first index corresponding to the rightmost bit (i.e. 2^0), 0 stands for spin down and 1 stands for spin up
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereTwoLandauLevels::ProdAd (int* m, int spinIndices, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  unsigned long TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (((spinIndices >> i) & 0x1) == 0)
	Index = ((m[i] + this->LzShiftDown) << 1) + ((spinIndices >> i) & 0x1);
      else
	Index = ((m[i] + this->LzShiftUp) << 1) + ((spinIndices >> i) & 0x1);
      if ((TmpState & (0x1l << Index)) != 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      if (Index > NewLzMax)
	{
	  NewLzMax = Index;
	}
      else
	{
	  coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 16))  & this->SignLookUpTableMask[Index + 16]];
#ifdef  __64_BITS__
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
	  coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
	}
      TmpState |= (0x1l << Index);
    }
  return this->FindStateIndex(TmpState, NewLzMax);
}

// generate all states corresponding to the constraints
// 
// nbrFermionsUp = number of fermions with spin up
// nbrFermionsDown = number of fermions with spin down
// lzMaxUp = momentum maximum value for a fermion
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereTwoLandauLevels::GenerateStates(int nbrFermionsUp, int nbrFermionsDown, int lzMax, int totalLz, long pos)
{
  cout << "warning : untested code" << endl;
  if ((nbrFermionsUp < 0) || (nbrFermionsDown < 0) || (totalLz < 0))
    return pos;
  if ((nbrFermionsUp == 0) && (totalLz == 0) && (nbrFermionsDown == 0))
      {
	this->StateDescription[pos] = 0x0ul;
	return (pos + 1l);
      }
    
  if (lzMax < 0)
    return pos;
  
  if ((nbrFermionsUp == 1) && (nbrFermionsDown == 0)) 
    {
      if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
	{
	  this->StateDescription[pos] = 0x2ul << (totalLz << 1);
	  return (pos + 1l);
	}
      else
	return pos;
    }
  if ((nbrFermionsDown == 1) && (nbrFermionsUp == 0)) 
    {
      if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
	{
	  this->StateDescription[pos] = 0x1ul << (totalLz << 1);
	  return (pos + 1l);
	}
      else
	return pos;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return pos;


  if (((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp)) &&
      ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown)))
    {
      long TmpPos = this->GenerateStates(nbrFermionsUp - 1, nbrFermionsDown - 1, lzMax - 1, totalLz - (2 * lzMax), pos);
      unsigned long Mask = 0x3ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      long TmpPos = this->GenerateStates(nbrFermionsUp - 1, nbrFermionsDown, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x2ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
    {
      long TmpPos = this->GenerateStates(nbrFermionsUp, nbrFermionsDown  - 1, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x1ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  return this->GenerateStates(nbrFermionsUp, nbrFermionsDown, lzMax - 1, totalLz, pos);
};

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereTwoLandauLevels::GenerateFullStates(int nbrFermions, int lzMax, int totalLz, long pos)
{
  if (nbrFermions < 0)
    return pos;
  if ((nbrFermions == 0) && (totalLz == 0))
    {
      this->StateDescription[pos] = 0x0ul;
      return (pos + 1l);
    }
    
  if (lzMax < 0)
    return pos;
  
  if (nbrFermions == 1) 
     {
       if (lzMax >= totalLz)
	 {
	   if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
	     {
	       this->StateDescription[pos] = 0x2ul << (totalLz << 1);
	       ++pos;
	     }
	   if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
	     {
	       this->StateDescription[pos] = 0x1ul << (totalLz << 1);
	       ++pos;
	     }
	 }
       return pos;
     }

  if ((lzMax == 0) && (totalLz != 0))
    return pos;

  if (((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp)) &&
      ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown)))
    {
      long TmpPos = this->GenerateFullStates(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax), pos);
      unsigned long Mask = 0x3ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      long TmpPos = this->GenerateFullStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x2ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
    {
      long TmpPos = this->GenerateFullStates(nbrFermions - 1, lzMax - 1, totalLz - lzMax,  pos);
      unsigned long Mask = 0x1ul << ((lzMax) << 1);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }

  return this->GenerateFullStates(nbrFermions, lzMax - 1, totalLz, pos);
}

// evaluate Hilbert space dimension
//
// nbrFermionsUp = number of fermions with spin up
// nbrFermionsDown = number of fermions with spin down
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereTwoLandauLevels::ShiftedEvaluateHilbertSpaceDimension(int nbrFermionsUp, int nbrFermionsDown, int lzMax, int totalLz)
{
  cout << "warning : untested code" << endl;
  if ((nbrFermionsUp < 0) || (nbrFermionsDown < 0) || (totalLz < 0))
    return 0l;
  if ((nbrFermionsUp == 0) && (nbrFermionsDown == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0) 
    return 0l;
    
  if ((nbrFermionsUp == 1) && (nbrFermionsDown == 0)) 
    {
      if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
	return 1l;
      else
	return 0l;
    }

  if ((nbrFermionsUp == 0) && (nbrFermionsDown == 1)) 
    {
      if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
	return 1l;
      else
	return 0l;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return 0l;

  long Tmp = 0l;
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
	Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp - 1, nbrFermionsDown - 1, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp - 1, nbrFermionsDown, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
    Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp, nbrFermionsDown - 1, lzMax - 1, totalLz - lzMax);
  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermionsUp, nbrFermionsDown, lzMax- 1, totalLz);
  return Tmp;
}

// evaluate Hilbert space dimension without constraint on the number of particles per level
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereTwoLandauLevels::ShiftedEvaluateFullHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions < 0) || (totalLz < 0))
    return 0l;
  if ((nbrFermions == 0) && (totalLz == 0))
    return 1l;
  if (lzMax < 0) 
    return 0l;
    
  if (nbrFermions == 1) 
    {
      long Tmp = 0l;
      if (lzMax >= totalLz)
	{
	  if (((this->LzMaxUp + this->LzShiftUp) >= totalLz) && (totalLz >= this->LzShiftUp))
	    ++Tmp;
	  if (((this->LzMaxDown + this->LzShiftDown) >= totalLz) && (totalLz >= this->LzShiftDown))
	    ++Tmp;
	}
      return Tmp;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return 0l;

  long Tmp = 0l;
  if ((lzMax <= (this->LzMaxUp + this->LzShiftUp)) && (lzMax >= this->LzShiftUp))
    {
      if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))
	Tmp += this->ShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 2, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += this->ShiftedEvaluateFullHilbertSpaceDimension(nbrFermions -1, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (this->LzMaxDown + this->LzShiftDown)) && (lzMax >= this->LzShiftDown))    
    Tmp += this->ShiftedEvaluateFullHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax);
  Tmp += this->ShiftedEvaluateFullHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz);
  return Tmp;
}


// create an SU(2) state from two U(1) state
//
// upState = vector describing the up spin part of the output state
// upStateSpace = reference on the Hilbert space associated to the up spin part
// downState = vector describing the down spin part of the output state
// downStateSpace = reference on the Hilbert space associated to the down spin part  
// return value = resluting SU(2) state

RealVector FermionOnSphereTwoLandauLevels::ForgeSU2FromU1(RealVector& upState, FermionOnSphere& upStateSpace, RealVector& downState, FermionOnSphere& downStateSpace)
{
  RealVector FinalState(this->HilbertSpaceDimension, true);
  for (int j = 0; j < upStateSpace.HilbertSpaceDimension; ++j)
    {
      unsigned long TmpUpState = upStateSpace.StateDescription[j] << this->LzShiftUp;
      int TmpPos = upStateSpace.LzMax + this->LzShiftUp;
      while (TmpPos > 0)
	{
	  unsigned long Tmp = TmpUpState & (0x1ul << TmpPos);
	  TmpUpState |= Tmp << TmpPos;
	  TmpUpState ^= Tmp;
	  --TmpPos;
	}
      TmpUpState <<= 1;
      double TmpComponent = upState[j];
      int Max = 63;
      while ((TmpUpState & (0x1ul << Max)) == 0x0ul)
	--Max;
      int Min = 0;
      while ((TmpUpState & (0x1ul << Min)) == 0x0ul)
	++Min;
      unsigned long TmpUpStateMask = (0x1ul << Max) - 1;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	if ((this->StateDescription[i] & TmpUpState) == TmpUpState)
	  {	    
	    unsigned long TmpUpState3 = this->StateDescription[i] & TmpUpStateMask;
	    unsigned long TmpUpState2 = TmpUpState3;
#ifdef  __64_BITS__
	    TmpUpState3 &= 0x5555555555555555ul;
	    TmpUpState2 &= 0xaaaaaaaaaaaaaaaaul;
#else
	    TmpUpState3 &= 0x55555555ul;
	    TmpUpState2 &= 0xaaaaaaaaul;
#endif	    
	    unsigned long Sign = 0x0;
	    int Pos = this->LzMax << 1;
	    while ((Pos > 0) && ((TmpUpState3 & (0x1ul << Pos)) == 0x0ul))
	      Pos -= 2;
	    while (Pos > 0)
	      {
		unsigned long TmpUpState4 = TmpUpState2 & ((0x1ul << Pos) - 1ul);
#ifdef  __64_BITS__
		TmpUpState4 ^= TmpUpState4 >> 32;
#endif	
		TmpUpState4 ^= TmpUpState4 >> 16;
		TmpUpState4 ^= TmpUpState4 >> 8;
		TmpUpState4 ^= TmpUpState4 >> 4;
		TmpUpState4 ^= TmpUpState4 >> 2;
		TmpUpState4 ^= TmpUpState4 >> 1;
		Sign ^= TmpUpState4;
		Pos -= 2;
		while ((Pos > 0) && ((TmpUpState3 & (0x1ul << Pos)) == 0x0ul))
		  Pos -= 2;
	      }
	    if ((Sign & 0x1ul) == 0x0ul)
	      FinalState[i] = TmpComponent;
	    else
	      FinalState[i] = -TmpComponent;
	  }
    }

  for (int j = 0; j < downStateSpace.HilbertSpaceDimension; ++j)
    {
      unsigned long TmpDownState = downStateSpace.StateDescription[j] << this->LzShiftDown;
      int TmpPos = downStateSpace.LzMax + this->LzShiftDown;
      while (TmpPos > 0)
	{
	  unsigned long Tmp = TmpDownState & (0x1ul << TmpPos);
	  TmpDownState |= Tmp << TmpPos;
	  TmpDownState ^= Tmp;
	  --TmpPos;
	}
      double TmpComponent = downState[j];
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	if ((this->StateDescription[i] & TmpDownState) == TmpDownState)
	  {
	    FinalState[i] *= TmpComponent;
	  }
    }

  return FinalState;
}

// compute the product of a bosonic state and a fermionic state belonging in two Landau levels
//
// bosonState = real vector where the bosonic state is stored
// fermionState = real vector where the fermionic state is stored
// outputVector = real vector where the result has to be stored
// finalStates = array where the obtained states are stored in their fermionic representation
// weigth = array where the coefficients for each obtained state are stored
// bosonSpace = pointer to the bosonic Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereTwoLandauLevels::BosonicStateTimeFermionicState(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, unsigned long * finalStates, 
								    double * weigth, BosonOnSphereShort * bosonSpace, FermionOnSphereTwoLandauLevels * finalSpace, int firstComponent, int nbrComponent)
{
  unsigned long* Monomial = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  int MaxComponent = firstComponent + nbrComponent;
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j]!=0)
	{
	  this->ConvertToMonomial(j, Slater);		
	  for (int i = firstComponent; i < MaxComponent; i++)
	    {
	      if(bosonState[i] != 0)
		{
		  bosonSpace->GetMonomial(i, Monomial);
		  for (int Index=0; Index < this->NbrFermions;Index++)
		    Monomial[Index]*=2;
		  unsigned int Limit=this->MonomialsTimesSlater(Slater,Monomial, finalStates, weigth,finalSpace);
		  for (unsigned int Index = 0; Index<Limit; Index++)
		    {
		      int TmpLzMax = 2 * finalSpace->LzMaxUp + 1;
		      while ((finalStates[Index] >> TmpLzMax) == 0x0ul)
			--TmpLzMax;
		      outputVector[finalSpace->FindStateIndex(finalStates[Index],TmpLzMax)] += bosonState[i] * fermionState[j] * weigth[Index];
		    }
		}
	    }
	}
    }
}

// compute the product of a monomial and a Slater determinant belonging in two Landau levels
// 
// slater = array where the Slater determinant is stored in its monomial representation
// monomial = array where the monomial is stored in its monomial representation
// finalStates = array where the obtained states are stored in their fermionic representation
// weigth = array where the coefficients for each obtained state are stored
// finalSpace = pointer to the final HilbertSpace
// return value = number of different obtained states

unsigned int FermionOnSphereTwoLandauLevels::MonomialsTimesSlater(unsigned long* slater,unsigned long* monomial, unsigned long*& finalStates, double*& weigth,FermionOnSphereTwoLandauLevels* finalSpace)
{
  unsigned int NbrNonZero = 0;
  unsigned long * State =new unsigned long[this->NbrFermions];
  unsigned long * TmpSlater = new unsigned long [this->NbrFermions];
  int NbrPermutation = 0;
  bool Bool = true;
  double Coef = 1.0;
  int k = 1;
  for (int Index = 0; Index < this->NbrFermions; Index++)
    State[Index] = slater[Index] + monomial[Index];
  finalSpace->GeneratesDifferentState(finalStates, weigth, slater, State, this, 0, NbrNonZero, Coef);
  while (std::prev_permutation(monomial, monomial + this->NbrFermions))
    {		
      NbrPermutation = 0;
      Bool = true;
      for (int Index = 0; Index < this->NbrFermions;Index++)
	{
	  State[Index] = slater[Index] + monomial[Index];
	  TmpSlater[Index] = slater[Index];
	}
      SortArrayDownOrdering(State, TmpSlater, this->NbrFermions, NbrPermutation);
      k = 1;
      while((k<this->NbrFermions)&&(Bool))
	{
	  if(State[k-1] == State[k])
	    {
	      Bool=false;
	    }
	  k++;
	}
      if(Bool)
	{			
	  ((NbrPermutation & 1) == 0)?Coef=1.0:Coef=-1.0;
	  finalSpace->GeneratesDifferentState(finalStates, weigth,TmpSlater, State, this, 0, NbrNonZero, Coef);
	}
    }
  delete [] State;
  delete [] TmpSlater;
  return NbrNonZero;
}

// generate the different states that appear in the product of a monomial and a Slater determinant in the two Landau levels
//
// finalStates = array where the obtained states has to be stored in their fermionic representation
// weigth = array where the coefficients for each obtained state has to be stored
// slater = array where the Slater determinant is stored in its monomial representation
// state = array where the obtained state is stored in its monomial representation
// slaterSpace = pointer to the Hilbert Space which the Slater determinant belongs to
// index = index of the particle being examinate
// nbrNonZero = number of different obtained states
// coef = coefficient of the state being generate

void FermionOnSphereTwoLandauLevels::GeneratesDifferentState(unsigned long * finalStates,double * weigth,unsigned long * slater,unsigned long * state,FermionOnSphereTwoLandauLevels * slaterSpace,int index,unsigned int & nbrNonZero,double coef)
{
  while((((state[index]&0x01ul)==0x0ul)||(state[index]==((this->LzMaxUp<<1)+0x1ul)))&&(index<(this->NbrFermions-1)))
    index++;
  if((index==this->NbrFermions-1)||((state[index]>>1)==0))
    {
      unsigned long TmpState = this->ConvertFromMonomial(state);
      nbrNonZero+=SearchInArrayAndSetWeight(TmpState,finalStates,weigth,nbrNonZero,coef);
      if(((state[index]&0x01ul)==1)&&((state[index]>>1)!=0))
	{
	  coef*= -((double)((double)(slater[index]>>1)/(double)slaterSpace->LzMaxUp)-((double)(state[index]>>1)/(double)this->LzMaxUp));
	  state[index]--;
	  unsigned long TmpState = this->ConvertFromMonomial(state);
	  nbrNonZero += SearchInArrayAndSetWeight(TmpState,finalStates,weigth,nbrNonZero,coef);
	  state[index]++;
	}
      return;
    }
  this->GeneratesDifferentState(finalStates,weigth,slater,state,slaterSpace,index+1,nbrNonZero,coef);
  
  if(state[index]!=(state[index+1]+1))
    {
      coef*= -((double)((double)(slater[index]>>1)/(double)slaterSpace->LzMaxUp)-((double)(state[index]>>1)/(double)this->LzMaxUp));
      state[index]--;
      this->GeneratesDifferentState(finalStates,weigth,slater,state,slaterSpace,index+1,nbrNonZero,coef);
      state[index]++;
    }
}


// compute the projection of the product of a bosonic state and a fermionic state
//
// bosonState = real vector where the bosonic state is stored
// fermionState = real vector where the fermionic state is stored
// outputVector = real vector where the result has to be stored
// finalStates = array where the obtained states are stored in their fermionic representation
// weigth = array where the coefficients for each obtained state are stored
// bosonSpace = pointer to the bosonic Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereTwoLandauLevels::BosonicStateTimeFermionicState(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, unsigned long* finalStates, 
								    double* weigth, BosonOnSphereShort* bosonSpace,FermionOnSphere* finalSpace, int firstComponent,int nbrComponent)
{
  unsigned long* Monomial = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  int NbrMax = firstComponent+nbrComponent;
  int NbrVariable = 0;
  unsigned long* Variable = new unsigned long[this->NbrFermions];
  for (int j = 0; j < this->HilbertSpaceDimension; j++)
    {
      if(fermionState[j] != 0)
	{
	  this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
	  for (int i = firstComponent; i < NbrMax; i++)
	    {
	      if(bosonState[i] != 0)
		{
		  bosonSpace->GetMonomial(i, Monomial);
		  unsigned int Limit = this->MonomialsTimesSlaterProjection(Slater, Monomial, Variable, NbrVariable, finalStates, weigth, finalSpace);
		  for (unsigned int Index=  0; Index < Limit; Index++)
		    {
		      int TmpLzMax = finalSpace->LzMax;
		      while (((finalStates[Index] >> TmpLzMax) & 0x1ul) == 0x0ul)
			--TmpLzMax;
		      outputVector[finalSpace->FindStateIndex(finalStates[Index],TmpLzMax)] += bosonState[i] * fermionState[j] * weigth[Index];
		    }
		}
	    }
	}
    }
}

// compute the projection of the product of a bosonic state and a fermionic state using the lz->-lz symmetry
//
// bosonState = real vector where the bosonic state is stored
// fermionState = real vector where the fermionic state is stored
// outputVector = real vector where the result has to be stored
// finalStates = array where the obtained states are stored in their fermionic representation
// weigth = array where the coefficients for each obtained state are stored
// bosonSpace = pointer to the bosonic Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereTwoLandauLevels::BosonicStateTimeFermionicStateSymmetric(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, unsigned long* finalStates, double* weigth, 
									     BosonOnSphereShort* bosonSpace,FermionOnSphere* finalSpace, int firstComponent,int nbrComponent)
{
  unsigned long* Monomial = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  int NbrMax = firstComponent+nbrComponent;
  int NbrVariable=0;
  unsigned long* Variable = new unsigned long[this->NbrFermions];
  for (int j = 0 ; j < this->HilbertSpaceDimension ; j++)
    {
      if(fermionState[j]!=0)
	{
	  this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
	  for (int i = firstComponent ; i < NbrMax ; i++)
	    {
	      if(bosonState[i]!=0)
		{
		  unsigned long TmpState=bosonSpace->FermionBasis->GetSymmetricState(bosonSpace->FermionBasis->StateDescription[i]);
		  int BTmpLzMax = bosonSpace->LzMax+this->NbrFermions-1;
		  while (((TmpState >> BTmpLzMax) & 0x1ul) == 0x0ul)
		    --BTmpLzMax;
		  if(bosonSpace->FermionBasis->FindStateIndex(TmpState,BTmpLzMax)>i)
		    {
		      bosonSpace->GetMonomial(i,Monomial);
		      unsigned int Limit = this->MonomialsTimesSlaterProjection(Slater,Monomial,Variable,NbrVariable,finalStates,weigth,finalSpace);
		      for (unsigned int Index=0; Index < Limit;Index++)
			{
			  int TmpLzMax = finalSpace->LzMax;
			  while (((finalStates[Index] >> TmpLzMax) & 0x1ul) == 0x0ul)
			    --TmpLzMax;
			  outputVector[finalSpace->FindStateIndex(finalStates[Index],TmpLzMax)]+=bosonState[i]*fermionState[j]*weigth[Index];
			  unsigned long TmpState = finalSpace->GetSymmetricState(finalStates[Index]);
			  TmpLzMax = finalSpace->LzMax;
			  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
			    --TmpLzMax;
			  outputVector[finalSpace->FindStateIndex(TmpState,TmpLzMax)]+=bosonState[i]*fermionState[j]*weigth[Index];
			}
		    }
		  if(bosonSpace->FermionBasis->FindStateIndex(TmpState,BTmpLzMax)==i)
		    {
		      bosonSpace->GetMonomial(i, Monomial);
		      unsigned int Limit=this->MonomialsTimesSlaterProjection(Slater,Monomial,Variable,NbrVariable,finalStates,weigth,finalSpace);
		      for (unsigned int Index=0; Index<Limit;Index++)
			{
			  int TmpLzMax = finalSpace->LzMax;
			  while (((finalStates[Index] >> TmpLzMax) & 0x1ul) == 0x0ul)
			    --TmpLzMax;
			  outputVector[finalSpace->FindStateIndex(finalStates[Index],TmpLzMax)]+=bosonState[i]*fermionState[j]*weigth[Index];
			}
		    }
		}
	    }
	}
    }
}

// compute the product and the projection of a Slater determinant and a monomial 
// 
// slater = array where the slater is stored in its monomial representation
// monomial = array where the monomial is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// finalStates = array where the obtained states are stored in their fermionic representation
// weigth = array where the coefficients for each obtained state are stored
// finalSpace = pointer to the final HilbertSpace
// return value = number of different obtained states

unsigned int FermionOnSphereTwoLandauLevels::MonomialsTimesSlaterProjection(unsigned long* slater,unsigned long* monomial,unsigned long * variable,int nbrVariable, unsigned long*& finalStates, 
									    double*& weigth,FermionOnSphere* finalSpace)
{
  unsigned int NbrNonZero = 0;
  unsigned long* State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0;
  bool Bool = true;
  double Coef = 1.0;
  long TmpLzMaxUp = this->LzMaxUp;
  long TmpFinalLzMaxUp = 2l + finalSpace->LzMax;
  double InverseFactor = 1.0 / (((double) TmpLzMaxUp) * ((double) TmpFinalLzMaxUp));
  for (int i = 0; i < this->NbrFermions ; i++)
    State[i] = slater[i] + monomial[i];
  for(int k = 0 ; (k < nbrVariable) && (Coef != 0.0); k++)
    {
      long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
      if (Numerator == 0l)
	Coef = 0.0;
      else
	Coef *= ((double) Numerator) * InverseFactor;
    }
  
  unsigned long Mask;
  unsigned long Sign = 0ul;
  if(Coef != 0.0)
    {
      for(int i = 0; (i < this->NbrFermions)&&(Bool); i++)
	{
	  Mask=(1ul << (State[i]-1));
	  if ( (TmpState & Mask) != 0ul)
	    Bool = false;
	  unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef _64_BITS__
	  TmpState2 ^= TmpState2 >> 32;
#endif
	  TmpState2 ^= TmpState2 >> 16;
	  TmpState2 ^= TmpState2 >> 8;
	  TmpState2 ^= TmpState2 >> 4;
	  TmpState2 ^= TmpState2 >> 2;
	  TmpState2 ^= TmpState2 >> 1;
	  Sign ^= TmpState2;
	  TmpState |= Mask;
	}
      if(Bool)
	{
	  ((Sign & 0x1ul) == 0ul)?:Coef*=-1.0;
	  NbrNonZero+=SearchInArrayAndSetWeight(TmpState,finalStates,weigth,NbrNonZero,Coef);
	}
    }
  while (std::prev_permutation(monomial,monomial+this->NbrFermions))
    {
      Coef=1.0;
      for(int i = 0; i < this->NbrFermions; i++)
	{
	  State[i] = slater[i] + monomial[i];
	}
      for(int k = 0; (k < nbrVariable) && (Coef != 0.0); k++)
	{
	  long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
	  if (Numerator == 0l)
	    Coef = 0.0;
	  else
	    Coef *= ((double) Numerator) * InverseFactor;
	}
      if (Coef != 0.0)
	{
	  Bool = true;
	  TmpState = 0ul;
	  Sign = 0ul;
	  for (int i=0; (i < this->NbrFermions)&&(Bool);i++)
	    {
	      Mask = (1ul << (State[i] - 1));
	      if((TmpState&Mask) != 0)
		Bool=false;
	      unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef  __64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif	
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  if(Bool)
	    {
	      ((Sign & 0x1ul) == 0ul)?:Coef*=-1.0;
	      NbrNonZero += SearchInArrayAndSetWeight(TmpState,finalStates,weigth,NbrNonZero,Coef);
	    }
	}
    }
  delete [] State;
  return NbrNonZero;
}

// compute the projection of the product of a bosonic state and a fermionic state
//
// lllFermionState = real vector where the lowest Landau level fermionic state is stored
// fermionState = real vector where the two Landau level fermionic state is stored
// outputVector = real vector where the result has to be stored
// finalStates = array where the obtained states are stored in their fermionic representation
// weigth = array where the coefficients for each obtained state are stored
// lllFermionSpace = pointer to the lowest Landau level Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereTwoLandauLevels::LLLFermionicStateTimeFermionicState(RealVector& lllFermionState, RealVector& fermionState, RealVector& outputVector, unsigned long* finalStates, double* weigth, 
									 FermionOnSphere* lllFermionSpace,BosonOnSphereShort* finalSpace, int firstComponent,int nbrComponent)
{
        unsigned long* LLLSlater = new unsigned long[this->NbrFermions];
        unsigned long* Slater = new unsigned long[this->NbrFermions];
        int NbrMax = firstComponent+nbrComponent;
        int NbrVariable = 0;
        FactorialCoefficient Coefficient;
        unsigned long* Variable = new unsigned long[this->NbrFermions];
        for (int j = 0; j < this->HilbertSpaceDimension; j++)
        {
                if(fermionState[j] != 0)
                {
                        NbrVariable=0;
                        this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
                        for (int i = firstComponent; i < NbrMax; i++)
                        {
                                if(lllFermionState[i] != 0)
                                {
                                        lllFermionSpace->GetMonomial(i, LLLSlater);
                                        unsigned int Limit=this->SlaterTimesSlaterProjection(Slater,LLLSlater,Variable,NbrVariable,finalStates,weigth,finalSpace);
                                        for (unsigned int Index=0; Index<Limit;Index++)
                                        {
                                                int FTmpLzMax = finalSpace->LzMax+this->NbrFermions-1;
                                                while (((finalStates[Index] >> FTmpLzMax) & 0x1ul) == 0x0ul)
                                                        --FTmpLzMax;
                                                finalSpace->FermionToBoson(finalStates[Index],FTmpLzMax,finalSpace->TemporaryState,finalSpace->TemporaryStateLzMax);
                                                Coefficient.SetToOne();
                                                for(int p=0;p<finalSpace->TemporaryStateLzMax+1;p++)
                                                {
                                                        Coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
                                                }
                                                outputVector[finalSpace->FermionBasis->FindStateIndex(finalStates[Index],FTmpLzMax)] += lllFermionState[i]*fermionState[j]*weigth[Index]*Coefficient.GetIntegerValue();
                                        }
                                }
                        }
                }
        }
}





// compute the projection of the product of a bosonic state and a fermionic state using lz->-lz symmetry
//
// lllFermionState = real vector where the lowest Landau level fermionic state is stored
// fermionState = real vector where the two Landau level fermionic state is stored
// outputVector = real vector where the result has to be stored
// finalStates = array where the obtained states are stored in their fermionic representation
// weigth = array where the coefficients for each obtained state are stored
// lllFermionSpace = pointer to the lowest Landau level Hilbert space
// finalSpace = pointer to the final Hilbert space
// firstComponent = first component to be computed
// nbrComponent = number of components to be computed

void FermionOnSphereTwoLandauLevels::LLLFermionicStateTimeFermionicStateSymmetric(RealVector& lllFermionState, RealVector& fermionState, RealVector& outputVector, unsigned long* finalStates, double* weigth,
										  FermionOnSphere* lllFermionSpace, BosonOnSphereShort* finalSpace, int firstComponent, int nbrComponent)
{
  unsigned long* LLLSlater = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  int NbrMax = firstComponent+nbrComponent;
  int NbrVariable = 0;
  unsigned long* Variable = new unsigned long[this->NbrFermions];
  FactorialCoefficient coefficient;
  for (int j = 0 ; j < this->HilbertSpaceDimension ; j++)
    {
      if(fermionState[j] != 0)
	{
	  NbrVariable = 0;
	  this->ConvertToMonomialVariable(this->StateDescription[j], Slater, NbrVariable, Variable);
	  for (int i = firstComponent ; i < NbrMax ; i++)
	    {
	      if(lllFermionState[i]!=0)
		{
		  unsigned long TmpState=lllFermionSpace->GetSymmetricState(lllFermionSpace->StateDescription[i]);
		  int TmpLzMax = lllFermionSpace->LzMax;
		  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		    --TmpLzMax;
		  if(lllFermionSpace->FindStateIndex(TmpState,TmpLzMax) > i)
		    {
		      lllFermionSpace->GetMonomial(i,LLLSlater);
		      unsigned int Limit = this->SlaterTimesSlaterProjection(Slater, LLLSlater, Variable, NbrVariable, finalStates, weigth, finalSpace);
		      for (unsigned int Index=0; Index < Limit;Index++)
			{
			  int FTmpLzMax = finalSpace->LzMax+this->NbrFermions-1;
			  while (((finalStates[Index] >> FTmpLzMax) & 0x1ul) == 0x0ul)
			    --FTmpLzMax;
			  finalSpace->FermionToBoson(finalStates[Index], FTmpLzMax, finalSpace->TemporaryState, finalSpace->TemporaryStateLzMax);
			  coefficient.SetToOne();
			  for(int p = 0;p < (finalSpace->TemporaryStateLzMax + 1); p++)
			    {
			      coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
			    }
			  outputVector[finalSpace->FermionBasis->FindStateIndex(finalStates[Index], FTmpLzMax)] += lllFermionState[i] * fermionState[j] * weigth[Index] * coefficient.GetIntegerValue();
			  unsigned long TmpState = finalSpace->FermionBasis->GetSymmetricState(finalStates[Index]);
			  FTmpLzMax = finalSpace->LzMax + this->NbrFermions - 1;
			  while (((TmpState >> FTmpLzMax) & 0x1ul) == 0x0ul)
			    --FTmpLzMax;
			  finalSpace->FermionToBoson(finalStates[Index], FTmpLzMax, finalSpace->TemporaryState, finalSpace->TemporaryStateLzMax);
			  coefficient.SetToOne();
			  for(int p = 0;p < (finalSpace->TemporaryStateLzMax + 1); p++)
			    {
			      coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
			    }
			  outputVector[finalSpace->FermionBasis->FindStateIndex(TmpState,FTmpLzMax)] += lllFermionState[i] * fermionState[j] * weigth[Index] * coefficient.GetIntegerValue();
			}
		    }
		  if(lllFermionSpace->FindStateIndex(TmpState,TmpLzMax) == i)
		    {
		      lllFermionSpace->GetMonomial(i, LLLSlater);
		      unsigned int Limit = this->SlaterTimesSlaterProjection(Slater, LLLSlater, Variable, NbrVariable, finalStates, weigth, finalSpace);
		      for (unsigned int Index=  0; Index < Limit; Index++)
			{
			  int FTmpLzMax = finalSpace->LzMax+this->NbrFermions-1;
			  while (((finalStates[Index] >> FTmpLzMax) & 0x1ul) == 0x0ul)
			    --FTmpLzMax;
			  finalSpace->FermionToBoson(finalStates[Index], FTmpLzMax, finalSpace->TemporaryState, finalSpace->TemporaryStateLzMax);
			  coefficient.SetToOne();
			  for(int p = 0; p < (finalSpace->TemporaryStateLzMax + 1); p++)
			    {
			      coefficient.FactorialMultiply(finalSpace->TemporaryState[p]);
			    }			  
			  outputVector[finalSpace->FermionBasis->FindStateIndex(finalStates[Index],FTmpLzMax)] += lllFermionState[i] * fermionState[j] * weigth[Index] * coefficient.GetIntegerValue();
			}
		    }
		}
	    }
	}
    }
}


// compute the product and the projection of a Slater determinant in the LLL and a Slater determinant in two Landau levels
//
// slater = array where the slater determinant in the two landau levels is stored in its monomial representation
// lllslater = array where the slater determinant in the LLL is stored in its monomial representation
// variable = reference on the array where the indice of fermions in the second Landau level is stored
// nbrVariable = number of fermions in the second Landau level
// finalStates = array where the obtained states are stored in their fermionic representation
// weigth = array where the coefficients for each obtained state are stored
// finalSpace = pointer to the final HilbertSpace
// return value = number of different obtained states

unsigned int FermionOnSphereTwoLandauLevels::SlaterTimesSlaterProjection(unsigned long* slater,unsigned long* lllslater,unsigned long * variable,int nbrVariable, unsigned long*& finalStates,
									 double*& weigth, BosonOnSphereShort* finalSpace)
{
  unsigned int NbrStates = 0;
  unsigned long * State = new unsigned long[this->NbrFermions];
  unsigned long TmpState = 0ul;
  double Coef = 1.0;
  long TmpLzMaxUp = this->LzMaxUp;
  long TmpFinalLzMaxUp = 2l + finalSpace->LzMax;
  double InverseFactor = 1.0 / (((double) TmpLzMaxUp) * ((double) TmpFinalLzMaxUp));
  for (int i = 0; i < this->NbrFermions ; i++)
    State[i] = slater[i] + lllslater[i];
  for(int k = 0 ; (k < nbrVariable) && (Coef != 0.0); k++)
    {
      long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
      if (Numerator == 0l)
	Coef = 0.0;
      else
	Coef *= ((double) Numerator) * InverseFactor;
    }
  
  unsigned long Mask=0ul;
  unsigned long Sign = 0ul;
  if(Coef != 0.0)
    {
      for (int i=0; (i < this->NbrFermions);i++)
	{
	  State[i]--;
	}
      NbrStates += SearchInArrayAndSetWeight(finalSpace->ConvertFromMonomial(State), finalStates, weigth, NbrStates, Coef);
    }
  while (std::prev_permutation(lllslater, lllslater + this->NbrFermions))
    {
      Coef = 1.0;
      for(int i = 0 ; i < this->NbrFermions; i++)
	State[i] = slater[i] + lllslater[i];
      for(int k = 0 ; (k < nbrVariable) && (Coef != 0.0); k++)
	{
	  long Numerator = -(slater[variable[k]] * TmpFinalLzMaxUp) + (State[variable[k]] * TmpLzMaxUp);
	  if (Numerator == 0l)
	    Coef = 0.0;
	  else
	    Coef *= ((double) Numerator) * InverseFactor;
	}
      if( Coef != 0.0 )
	{
	  TmpState = 0ul;
	  Sign = 0ul;
	  for (int i=0; (i < this->NbrFermions);i++)
	    {
	      State[i]--;
	      Mask = (1ul << lllslater[i]);
	      unsigned long TmpState2 = TmpState & (Mask - 1ul);
#ifdef  __64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState |= Mask;
	    }
	  SortArrayDownOrdering(State,this->NbrFermions);
	  ((Sign & 0x1ul) == 0ul)?:Coef*=-1.0;
	  NbrStates += SearchInArrayAndSetWeight(finalSpace->ConvertFromMonomial(State), finalStates, weigth, NbrStates, Coef);
	}
    }
  delete [] State;

  return NbrStates;
// original line was
//  return NbrNonZero;

}
