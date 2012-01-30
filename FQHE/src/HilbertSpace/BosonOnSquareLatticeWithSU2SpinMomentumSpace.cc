////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                class of bosons on a square lattice with SU(2) spin         //
//                                in momentum space                           //
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
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"
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
#include "GeneralTools/ArrayTools.h"

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


// basic constructor
// 
// nbrBosons = number of bosons
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

BosonOnSquareLatticeWithSU2SpinMomentumSpace::BosonOnSquareLatticeWithSU2SpinMomentumSpace (int nbrBosons, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->TotalSpin = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->Flag.Initialize();
  this->TemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->TemporaryStateDown = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateUp = new unsigned long[this->NbrLzValue];
  this->ProdATemporaryStateDown = new unsigned long[this->NbrLzValue];

  this->NUpLzMax = this->LzMax + this->NbrBosons - 1;
  this->NDownLzMax = this->LzMax + this->NbrBosons - 1;
  this->FermionicLzMax = this->NUpLzMax;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->StateDescriptionUp = new unsigned long [this->LargeHilbertSpaceDimension];
      this->StateDescriptionDown = new unsigned long [this->LargeHilbertSpaceDimension];
      this->Flag.Initialize();
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->LzMax + this->NbrBosons, this->LzMax + this->NbrBosons, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
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

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSquareLatticeWithSU2SpinMomentumSpace::BosonOnSquareLatticeWithSU2SpinMomentumSpace(const BosonOnSquareLatticeWithSU2SpinMomentumSpace& bosons)
{
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
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

BosonOnSquareLatticeWithSU2SpinMomentumSpace::~BosonOnSquareLatticeWithSU2SpinMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSquareLatticeWithSU2SpinMomentumSpace& BosonOnSquareLatticeWithSU2SpinMomentumSpace::operator = (const BosonOnSquareLatticeWithSU2SpinMomentumSpace& bosons)
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
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
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

AbstractHilbertSpace* BosonOnSquareLatticeWithSU2SpinMomentumSpace::Clone()
{
  return new BosonOnSquareLatticeWithSU2SpinMomentumSpace(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSquareLatticeWithSU2SpinMomentumSpace::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->StateDescriptionUp[state], this->StateDescriptionDown[state], TemporaryStateUp, TemporaryStateDown); 

  unsigned long Tmp;
  Str << "[";
  for (int i = 0; i <= this->LzMax; ++i)
    {
      int TmpKx = i / this->NbrSiteY;
      int TmpKy = i % this->NbrSiteY;
      if (this->TemporaryStateUp[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryStateUp[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << ",up)";
	}
      if (this->TemporaryStateDown[i] > 0)
	{
	  for (int j = 0; j < this->TemporaryStateDown[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << ",down)";
	}
    }
  Str << "]";
  //  Str << " = " << this->StateDescriptionUp[state] << "," << this->StateDescriptionDown[state];
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentFermionicPositionUp = current fermionic position within the state description for the spin up
// currentFermionicPositionDown = current fermionic position within the state description for the spin down
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSquareLatticeWithSU2SpinMomentumSpace::GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int currentFermionicPositionUp, int currentFermionicPositionDown, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
 	  this->StateDescriptionUp[pos] = 0x0ul;
 	  this->StateDescriptionDown[pos] = 0x0ul;
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  for (int i = nbrBosons; i >= 0; --i)
    {
      unsigned long MaskUp = ((0x1ul << i) - 0x1ul) << (currentFermionicPositionUp - i - 1);
      for (int j = nbrBosons - i; j >= 0; --j)
	{
	  long TmpPos = this->GenerateStates(nbrBosons - i - j, currentKx, currentKy - 1, currentTotalKx + ((i + j) * currentKx), currentTotalKy + ((i + j) * currentKy), currentFermionicPositionUp - i - 1, currentFermionicPositionDown - j - 1, pos);
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


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long BosonOnSquareLatticeWithSU2SpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrBosons == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    Count += 2l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += 2l;
	    }
	}
      return Count;
    }
  for (int i = nbrBosons; i >= 0; --i)
    Count += (((long) i) + 1l) * this->EvaluateHilbertSpaceDimension(nbrBosons - i, currentKx, currentKy - 1, currentTotalKx + (i * currentKx), currentTotalKy + (i * currentKy));
  return Count;
}
