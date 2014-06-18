////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//                                  in real space                             //
//                        Class author: Cecile Repellin                       //
//                                                                            //
//                        last modification : 17/06/2014                      //
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
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"

#include <math.h>
#include <cstdlib>
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

FermionOnLatticeWithSpinRealSpace::FermionOnLatticeWithSpinRealSpace ()
{
  this->NbrFermions = 0;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSite = 0;
  this->LzMax = this->NbrSite;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->StateHighestBit = 0;  
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinRealSpace::FermionOnLatticeWithSpinRealSpace (int nbrFermions, int nbrSite, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions);
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSite - 1, 0l);
      if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space, " << TmpLargeHilbertSpaceDimension << " generated states, should be " << this->LargeHilbertSpaceDimension << endl;
	}
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
}

// basic constructor when Sz is preserved
// 
// nbrFermions = number of fermions
//nbrSite = number of sites
// nbrSpinUp = number of particles with spin up
// memory = amount of memory granted for precalculations

// FermionOnSquareLatticeWithSpinMomentumSpace::FermionOnSquareLatticeWithSpinMomentumSpace (int nbrFermions, int nbrSpinUp, int nbrSite, unsigned long memory)
// {
//   
// }

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnLatticeWithSpinRealSpace::FermionOnLatticeWithSpinRealSpace(const FermionOnLatticeWithSpinRealSpace& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSite = fermions.NbrSite;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;
  this->HighestBit = fermions.HighestBit;
  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;  
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

FermionOnLatticeWithSpinRealSpace::~FermionOnLatticeWithSpinRealSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeWithSpinRealSpace& FermionOnLatticeWithSpinRealSpace::operator = (const FermionOnLatticeWithSpinRealSpace& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrSite = fermions.NbrSite;
  this->NbrLzValue = fermions.NbrLzValue;
  this->SzFlag = fermions.SzFlag;
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

AbstractHilbertSpace* FermionOnLatticeWithSpinRealSpace::Clone()
{
  return new FermionOnLatticeWithSpinRealSpace(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnLatticeWithSpinRealSpace::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      Tmp = (TmpState >> (i << 1)) & 03ul;
      switch (Tmp)
	{
	case 0x0ul:
	  Str << "0 ";
	  break;
	case 0x1ul:
	  Str << "d ";
	  break;
	case 0x2ul:
	  Str << "u ";
	  break;
	case 0x3ul:
	  Str << "X ";
	  break;
	}
    }
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinRealSpace::GenerateStates(int nbrFermions, int currentSite, long pos)
{
  if (nbrFermions == 0)
    {
      this->StateDescription[pos] = 0x0ul;	  
      return (pos + 1l);
    }
  if ((currentSite < 0) || (nbrFermions < 0))
    return pos;
  if (nbrFermions == 1)
    {
      for (int j = currentSite; j >= 0; --j)
	{
	  this->StateDescription[pos] = 0x2ul << (j << 1);
	  ++pos;
	  this->StateDescription[pos] = 0x1ul << (j << 1);
	  ++pos;
	}
      return pos;
    }
  long TmpPos = this->GenerateStates(nbrFermions - 2, currentSite - 1, pos);
  unsigned long Mask = 0x3ul << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, pos);
  Mask = 0x2ul << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
   TmpPos = this->GenerateStates(nbrFermions - 1, currentSite - 1, pos);
   Mask = 0x1ul << ((currentSite) << 1);
   for (; pos < TmpPos; ++pos)
     this->StateDescription[pos] |= Mask;
  return this->GenerateStates(nbrFermions, currentSite - 1, pos);
};

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// nbrSpinUp = number of fermions with spin up
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

// long FermionOnLatticeWithSpinRealSpace::GenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, long pos)
// {
//   
// }

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// return value = Hilbert space dimension
long FermionOnLatticeWithSpinRealSpace::EvaluateHilbertSpaceDimension(int nbrFermions)
{
  BinomialCoefficients binomials(2*this->NbrSite);
  long dimension = binomials(2*this->NbrSite, this->NbrFermions);
  return dimension;
}


// evaluate Hilbert space dimension with a fixed number of fermions with spin up
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrSpinUp = number of fermions with spin up
// return value = Hilbert space dimension
/*
long FermionOnLatticeWithSpinRealSpace::EvaluateHilbertSpaceDimension(int nbrFermions,int nbrSpinUp)
{
  
}*/

