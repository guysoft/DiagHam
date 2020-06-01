////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on a square lattice with SU(8) spin         //
//                                in momentum space                           //
//                                                                            //
//                        last modification : 11/05/2020                      //
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
#include "HilbertSpace/FermionOnSquareLatticeWithSU8SpinMomentumSpace.h"
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


// basic constructor
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU8SpinMomentumSpace::FermionOnSquareLatticeWithSU8SpinMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->SU4Flag = false;
  this->NbrFermions1 = 0;
  this->NbrFermions2 = 0;
  this->NbrFermions3 = 0;
  this->NbrFermions4 = 0;
  this->NbrFermions5 = 0;
  this->NbrFermions6 = 0;
  this->NbrFermions7 = 0;
  this->NbrFermions8 = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	this->PrintState(cout, i) << " " << hex << this->StateDescription[i] << dec << endl;
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

// constructor when preserving only a SU(4) degree of freedom
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// nbrParticles12 = number of particles with sigma=1 or sigma=2
// nbrParticles34 = number of particles with sigma=3 or sigma=4
// nbrParticles56 = number of particles with sigma=5 or sigma=6
// nbrParticles78 = number of particles with sigma=7 or sigma=8
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU8SpinMomentumSpace::FermionOnSquareLatticeWithSU8SpinMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum,
												int nbrParticles12, int nbrParticles34, int nbrParticles56, int nbrParticles78,
												unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->SU4Flag = true;
  this->NbrFermions1 = 0;
  this->NbrFermions2 = 0;
  this->NbrFermions3 = 0;
  this->NbrFermions4 = 0;
  this->NbrFermions5 = 0;
  this->NbrFermions6 = 0;
  this->NbrFermions7 = 0;
  this->NbrFermions8 = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0,
									 nbrParticles12, nbrParticles34, nbrParticles56, nbrParticles78);
  cout << this->LargeHilbertSpaceDimension << endl;
  
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0,
								nbrParticles12, nbrParticles34, nbrParticles56, nbrParticles78, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
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

// constructor when conserving all 8 Cartan related quantum numbers
// 
// nbrFermions = number of fermions
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// nbrParticleSigma = array that provides the number of particles with a given internal degree of freedom
// memory = amount of memory granted for precalculations

FermionOnSquareLatticeWithSU8SpinMomentumSpace::FermionOnSquareLatticeWithSU8SpinMomentumSpace (int nbrFermions, int nbrSiteX, int nbrSiteY,
												int kxMomentum, int kyMomentum,
												int* nbrParticleSigma, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->SU4Flag = true;
  this->NbrFermions1 = nbrParticleSigma[0];
  this->NbrFermions2 = nbrParticleSigma[1];
  this->NbrFermions3 = nbrParticleSigma[2];
  this->NbrFermions4 = nbrParticleSigma[3];
  this->NbrFermions5 = nbrParticleSigma[4];
  this->NbrFermions6 = nbrParticleSigma[5];
  this->NbrFermions7 = nbrParticleSigma[6];
  this->NbrFermions8 = nbrParticleSigma[7];
  this->TotalLz = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0,
									 this->NbrFermions1, this->NbrFermions2, this->NbrFermions3, this->NbrFermions4,
									 this->NbrFermions5, this->NbrFermions6, this->NbrFermions7, this->NbrFermions8);
  
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
      this->StateHighestBit = new int [this->HilbertSpaceDimension];  
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0,
								this->NbrFermions1, this->NbrFermions2,
								this->NbrFermions3, this->NbrFermions4,
								this->NbrFermions5, this->NbrFermions6,
								this->NbrFermions7, this->NbrFermions8, 0l);
      if (this->LargeHilbertSpaceDimension != TmpLargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space " << this->LargeHilbertSpaceDimension << " " << TmpLargeHilbertSpaceDimension << endl;
	}
//       for (int i = 0; i < this->HilbertSpaceDimension; ++i)
// 	this->PrintState(cout, i) << endl;
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


// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSquareLatticeWithSU8SpinMomentumSpace::FermionOnSquareLatticeWithSU8SpinMomentumSpace(const FermionOnSquareLatticeWithSU8SpinMomentumSpace& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->SzFlag = fermions.SzFlag;
  this->SU4Flag = fermions.SU4Flag;
  this->NbrFermions1 = fermions.NbrFermions1;
  this->NbrFermions2 = fermions.NbrFermions2;
  this->NbrFermions3 = fermions.NbrFermions3;
  this->NbrFermions4 = fermions.NbrFermions4;
  this->NbrFermions5 = fermions.NbrFermions5;
  this->NbrFermions6 = fermions.NbrFermions6;
  this->NbrFermions7 = fermions.NbrFermions7;
  this->NbrFermions8 = fermions.NbrFermions8;
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
}

// destructor
//

FermionOnSquareLatticeWithSU8SpinMomentumSpace::~FermionOnSquareLatticeWithSU8SpinMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSquareLatticeWithSU8SpinMomentumSpace& FermionOnSquareLatticeWithSU8SpinMomentumSpace::operator = (const FermionOnSquareLatticeWithSU8SpinMomentumSpace& fermions)
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
  this->NbrSiteX = fermions.NbrSiteX;
  this->NbrSiteY = fermions.NbrSiteY;
  this->KxMomentum = fermions.KxMomentum;
  this->KyMomentum = fermions.KyMomentum;
  this->NbrLzValue = fermions.NbrLzValue;
  this->SzFlag = fermions.SzFlag;
  this->SU4Flag = fermions.SU4Flag;
  this->NbrFermions1 = fermions.NbrFermions1;
  this->NbrFermions2 = fermions.NbrFermions2;
  this->NbrFermions3 = fermions.NbrFermions3;
  this->NbrFermions4 = fermions.NbrFermions4;
  this->NbrFermions5 = fermions.NbrFermions5;
  this->NbrFermions6 = fermions.NbrFermions6;
  this->NbrFermions7 = fermions.NbrFermions7;
  this->NbrFermions8 = fermions.NbrFermions8;
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

AbstractHilbertSpace* FermionOnSquareLatticeWithSU8SpinMomentumSpace::Clone()
{
  return new FermionOnSquareLatticeWithSU8SpinMomentumSpace(*this);
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSquareLatticeWithSU8SpinMomentumSpace::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  unsigned long Tmp;
  Str << "[";
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      Tmp = (TmpState >> (i << 3));
      int TmpKx = i / this->NbrSiteY;
      int TmpKy = i % this->NbrSiteY;
      for (int j = 0; j < 8; ++j)
	{
	  if ((Tmp & (0x1ul << j)) != 0ul)
	    Str << "(" << TmpKx << "," << TmpKy << "," << (j+1) << ")";
	}
    }
  Str << "]";
  return Str;
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU8SpinMomentumSpace::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions < 0)
    return pos;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    {
	      this->StateDescription[pos] = 0x80ul << (((currentKx * this->NbrSiteY) + j) << 3);
	      ++pos;
	      this->StateDescription[pos] = 0x40ul << (((currentKx * this->NbrSiteY) + j) << 3);
	      ++pos;
	      this->StateDescription[pos] = 0x20ul << (((currentKx * this->NbrSiteY) + j) << 3);
	      ++pos;
	      this->StateDescription[pos] = 0x10ul << (((currentKx * this->NbrSiteY) + j) << 3);
	      ++pos;
	      this->StateDescription[pos] = 0x8ul << (((currentKx * this->NbrSiteY) + j) << 3);
	      ++pos;
	      this->StateDescription[pos] = 0x4ul << (((currentKx * this->NbrSiteY) + j) << 3);
	      ++pos;
	      this->StateDescription[pos] = 0x2ul << (((currentKx * this->NbrSiteY) + j) << 3);
	      ++pos;
	      this->StateDescription[pos] = 0x1ul << (((currentKx * this->NbrSiteY) + j) << 3);
	      ++pos;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  this->StateDescription[pos] = 0x80ul << (((i * this->NbrSiteY) + j) << 3);
		  ++pos;
		  this->StateDescription[pos] = 0x40ul << (((i * this->NbrSiteY) + j) << 3);
		  ++pos;
		  this->StateDescription[pos] = 0x20ul << (((i * this->NbrSiteY) + j) << 3);
		  ++pos;
		  this->StateDescription[pos] = 0x10ul << (((i * this->NbrSiteY) + j) << 3);
		  ++pos;
		  this->StateDescription[pos] = 0x8ul << (((i * this->NbrSiteY) + j) << 3);
		  ++pos;
		  this->StateDescription[pos] = 0x4ul << (((i * this->NbrSiteY) + j) << 3);
		  ++pos;
		  this->StateDescription[pos] = 0x2ul << (((i * this->NbrSiteY) + j) << 3);
		  ++pos;
		  this->StateDescription[pos] = 0x1ul << (((i * this->NbrSiteY) + j) << 3);
		  ++pos;
		}
	    }
	}
      return pos;
    }


  long TmpPos;
  unsigned long Mask;
  int TmpNbrParticles;
  
  for (int i = 255; i >= 0; --i)
    {
      TmpNbrParticles = ((i & 1) + ((i & 2) >> 1) + ((i & 4) >> 2) + ((i & 8) >> 3)
			 + ((i & 16) >> 4) + ((i & 32) >> 5) + ((i & 64) >> 6) + ((i & 128) >> 7));
      TmpPos = this->GenerateStates(nbrFermions - TmpNbrParticles, currentKx, currentKy - 1,
				    currentTotalKx + (TmpNbrParticles * currentKx), currentTotalKy + (TmpNbrParticles * currentKy), pos);
      Mask = ((unsigned long) i) << (((currentKx * this->NbrSiteY) + currentKy) << 3);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
   }
  return pos;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrParticles1 = number of particles with sigma=1
// nbrParticles2 = number of particles with sigma=2
// nbrParticles3 = number of particles with sigma=3
// nbrParticles4 = number of particles with sigma=4
// nbrParticles5 = number of particles with sigma=5
// nbrParticles6 = number of particles with sigma=6
// nbrParticles7 = number of particles with sigma=7
// nbrParticles8 = number of particles with sigma=8
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU8SpinMomentumSpace::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy,
								    int nbrParticles1, int nbrParticles2, int nbrParticles3, int nbrParticles4,
								    int nbrParticles5, int nbrParticles6, int nbrParticles7, int nbrParticles8, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0)
      || (nbrParticles1 < 0) || (nbrParticles2 < 0) || (nbrParticles3 < 0) || (nbrParticles4 < 0)
      || (nbrParticles5 < 0) || (nbrParticles6 < 0) || (nbrParticles7 < 0) || (nbrParticles8 < 0)
      || (nbrParticles1 > nbrFermions) || (nbrParticles2 > nbrFermions) || (nbrParticles3 > nbrFermions) || (nbrParticles4 > nbrFermions)
      || (nbrParticles5 > nbrFermions) || (nbrParticles6 > nbrFermions) || (nbrParticles7 > nbrFermions) || (nbrParticles8 > nbrFermions))
    {
      return 0l;
    }
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
        {
            this->StateDescription[pos] = 0x0ul;	  
            return (pos + 1l);
        }
      else
	{
          return pos;
	}
    }
  if (currentKx < 0)
    return pos;
  if (nbrFermions == 1)
    {
      unsigned long TmpMask = ((unsigned long) (nbrParticles1 | (nbrParticles2 << 1) | (nbrParticles3 << 2) | (nbrParticles4 << 3)
						| (nbrParticles5 << 4) | (nbrParticles6 << 5) | (nbrParticles7 << 6)
						| (nbrParticles8 << 7)));
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) &&
	      (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    {
	      this->StateDescription[pos] = TmpMask << (((currentKx * this->NbrSiteY) + j) << 3);
	      ++pos;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) &&
		  (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  this->StateDescription[pos] = TmpMask << (((i * this->NbrSiteY) + j) << 3);
		  ++pos;
		}
	    }
	}
      return pos;
    }

  long TmpPos;
  unsigned long Mask;
  int TmpNbrParticles;  
  for (int i = 255; i >= 0; --i)
    {
      TmpNbrParticles = ((i & 1) + ((i & 2) >> 1) + ((i & 4) >> 2) + ((i & 8) >> 3)
			 + ((i & 16) >> 4) + ((i & 32) >> 5) + ((i & 64) >> 6) + ((i & 128) >> 7));
      TmpPos = this->GenerateStates(nbrFermions - TmpNbrParticles, currentKx, currentKy - 1,
				    currentTotalKx + (TmpNbrParticles * currentKx),
				    currentTotalKy + (TmpNbrParticles * currentKy),
				    nbrParticles1 - (i & 1), nbrParticles2 - ((i & 2) >> 1),
				    nbrParticles3 - ((i & 4) >> 2), nbrParticles4 - ((i & 8) >> 3),
				    nbrParticles5 - ((i & 16) >> 4), nbrParticles6 - ((i & 32) >> 5),
				    nbrParticles7 - ((i & 64) >> 6), nbrParticles8 - ((i & 128) >> 7), pos);
      Mask = ((unsigned long) i) << (((currentKx * this->NbrSiteY) + currentKy) << 3);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
   }
  return pos;
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrParticles12 = number of particles with sigma=1 or sigma=2
// nbrParticles34 = number of particles with sigma=3 or sigma=4
// nbrParticles56 = number of particles with sigma=5 or sigma=6
// nbrParticles78 = number of particles with sigma=7 or sigma=8
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSquareLatticeWithSU8SpinMomentumSpace::GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy,
								    int nbrParticles12, int nbrParticles34, int nbrParticles56, int nbrParticles78, long pos)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0)
      || (nbrParticles12 < 0) || (nbrParticles34 < 0) || (nbrParticles56 < 0) || (nbrParticles78 < 0)
      || (nbrParticles12 > nbrFermions) || (nbrParticles34 > nbrFermions) || (nbrParticles56 > nbrFermions) || (nbrParticles78 > nbrFermions))
    {
      return 0l;
    }
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
        {
            this->StateDescription[pos] = 0x0ul;	  
            return (pos + 1l);
        }
      else
          return pos;
    }
  if (currentKx < 0)
    return pos;
  if (nbrFermions == 1)
    {
      unsigned long TmpMask = ((unsigned long) (nbrParticles12 | (nbrParticles34 << 2) | (nbrParticles56 << 4) | (nbrParticles78 << 6)));
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) &&
	      (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    {
	      this->StateDescription[pos] = TmpMask << ((((currentKx * this->NbrSiteY) + j) << 3) + 1);
	      ++pos;
	      this->StateDescription[pos] = TmpMask << (((currentKx * this->NbrSiteY) + j) << 3);
	      ++pos;
	    }
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) &&
		  (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		{
		  this->StateDescription[pos] = TmpMask << ((((i * this->NbrSiteY) + j) << 3) + 1);
		  ++pos;
		  this->StateDescription[pos] = TmpMask << (((i * this->NbrSiteY) + j) << 3);
		  ++pos;
		}
	    }
	}
      return pos;
    }

  long TmpPos;
  unsigned long Mask;
  int TmpNbrParticles;  
  for (int i = 255; i >= 0; --i)
    {
      TmpNbrParticles = ((i & 1) + ((i & 2) >> 1) + ((i & 4) >> 2) + ((i & 8) >> 3)
			 + ((i & 16) >> 4) + ((i & 32) >> 5) + ((i & 64) >> 6) + ((i & 128) >> 7));
      TmpPos = this->GenerateStates(nbrFermions - TmpNbrParticles, currentKx, currentKy - 1,
				    currentTotalKx + (TmpNbrParticles * currentKx),
				    currentTotalKy + (TmpNbrParticles * currentKy),
				    nbrParticles12 - ((i & 1) + ((i & 2) >> 1)),
				    nbrParticles34 - (((i & 4) >> 2) + ((i & 8) >> 3)),
				    nbrParticles56 - (((i & 16) >> 4) + ((i & 32) >> 5)),
				    nbrParticles78 - (((i & 64) >> 6) + ((i & 128) >> 7)), pos);
      Mask = ((unsigned long) i) << (((currentKx * this->NbrSiteY) + currentKy) << 3);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
   }
  return pos;
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU8SpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrFermions < 0)
    return 0l;
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  return 1l;
	}
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    Count += 8l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Count += 8l;
	    }
	}
      return Count;
    }
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions - 8, currentKx, currentKy - 1, currentTotalKx + (8 * currentKx), currentTotalKy + (8 * currentKy));
  Count += (8 * this->EvaluateHilbertSpaceDimension(nbrFermions - 7, currentKx, currentKy - 1, currentTotalKx + (7 * currentKx), currentTotalKy + (7 * currentKy)));
  Count += (28 * this->EvaluateHilbertSpaceDimension(nbrFermions - 6, currentKx, currentKy - 1, currentTotalKx + (6 * currentKx), currentTotalKy + (6 * currentKy)));
  Count += (56 * this->EvaluateHilbertSpaceDimension(nbrFermions - 5, currentKx, currentKy - 1, currentTotalKx + (5 * currentKx), currentTotalKy + (5 * currentKy)));
  Count += (70 * this->EvaluateHilbertSpaceDimension(nbrFermions - 4, currentKx, currentKy - 1, currentTotalKx + (4 * currentKx), currentTotalKy + (4 * currentKy)));
  Count += (56 * this->EvaluateHilbertSpaceDimension(nbrFermions - 3, currentKx, currentKy - 1, currentTotalKx + (3 * currentKx), currentTotalKy + (3 * currentKy)));
  Count += (28 * this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentKx, currentKy - 1, currentTotalKx + (2 * currentKx), currentTotalKy + (2 * currentKy)));
  Count += (8 * this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentKx, currentKy - 1, currentTotalKx + currentKx, currentTotalKy + currentKy));
  Count += this->EvaluateHilbertSpaceDimension(nbrFermions, currentKx, currentKy - 1, currentTotalKx, currentTotalKy);
  return Count;
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrParticles1 = number of particles with sigma=1
// nbrParticles2 = number of particles with sigma=2
// nbrParticles3 = number of particles with sigma=3
// nbrParticles4 = number of particles with sigma=4
// nbrParticles5 = number of particles with sigma=5
// nbrParticles6 = number of particles with sigma=6
// nbrParticles7 = number of particles with sigma=7
// nbrParticles8 = number of particles with sigma=8
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU8SpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy,
										   int nbrParticles1, int nbrParticles2, int nbrParticles3, int nbrParticles4,
										   int nbrParticles5, int nbrParticles6, int nbrParticles7, int nbrParticles8)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0)
      || (nbrParticles1 < 0) || (nbrParticles2 < 0) || (nbrParticles3 < 0) || (nbrParticles4 < 0)
      || (nbrParticles5 < 0) || (nbrParticles6 < 0) || (nbrParticles7 < 0) || (nbrParticles8 < 0)
      || (nbrParticles1 > nbrFermions) || (nbrParticles2 > nbrFermions) || (nbrParticles3 > nbrFermions) || (nbrParticles4 > nbrFermions)
      || (nbrParticles5 > nbrFermions) || (nbrParticles6 > nbrFermions) || (nbrParticles7 > nbrFermions) || (nbrParticles8 > nbrFermions))
    {
      return 0l;
    }
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  return 1l;
	}
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  
  unsigned long Tmp = 0l;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    Tmp += 1l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Tmp += 1l;
	    }
	}
      return Tmp;
    }
  
  for (int i = 255; i >= 0; --i)
    {
      int TmpNbrParticles = ((i & 1) + ((i & 2) >> 1) + ((i & 4) >> 2) + ((i & 8) >> 3)
			     + ((i & 16) >> 4) + ((i & 32) >> 5) + ((i & 64) >> 6) + ((i & 128) >> 7));
      Tmp += this->EvaluateHilbertSpaceDimension(nbrFermions - TmpNbrParticles, currentKx, currentKy - 1,
						 currentTotalKx + (TmpNbrParticles * currentKx), currentTotalKy + (TmpNbrParticles * currentKy),
						 nbrParticles1 - (i & 1), nbrParticles2 - ((i & 2) >> 1),
						 nbrParticles3 - ((i & 4) >> 2), nbrParticles4 - ((i & 8) >> 3),
						 nbrParticles5 - ((i & 16) >> 4), nbrParticles6 - ((i & 32) >> 5),
						 nbrParticles7 - ((i & 64) >> 6), nbrParticles8 - ((i & 128) >> 7));
    }
  return Tmp;
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// nbrParticles12 = number of particles with sigma=1 or sigma=2
// nbrParticles34 = number of particles with sigma=3 or sigma=4
// nbrParticles56 = number of particles with sigma=5 or sigma=6
// nbrParticles78 = number of particles with sigma=7 or sigma=8
// return value = Hilbert space dimension

long FermionOnSquareLatticeWithSU8SpinMomentumSpace::EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy,
										   int nbrParticles12, int nbrParticles34, int nbrParticles56, int nbrParticles78)
{
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if ((nbrFermions < 0)
      || (nbrParticles12 < 0) || (nbrParticles34 < 0) || (nbrParticles56 < 0) || (nbrParticles78 < 0)
      || (nbrParticles12 > nbrFermions) || (nbrParticles34 > nbrFermions) || (nbrParticles56 > nbrFermions) || (nbrParticles78 > nbrFermions))
    {
      return 0l;
    }
  if (nbrFermions == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  return 1l;
	}
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  
  unsigned long Tmp = 0l;
  if (nbrFermions == 1)
    {
      for (int j = currentKy; j >= 0; --j)
	{
	  if ((((currentKx + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
	    Tmp += 2l;
	}
      for (int i = currentKx - 1; i >= 0; --i)
	{
	  for (int j = this->NbrSiteY - 1; j >= 0; --j)
	    {
	      if ((((i + currentTotalKx) % this->NbrSiteX) == this->KxMomentum) && (((j + currentTotalKy) % this->NbrSiteY) == this->KyMomentum))
		Tmp += 2l;
	    }
	}
      return Tmp;
    }
  
  for (int i = 255; i >= 0; --i)
    {
      int TmpNbrParticles = ((i & 1) + ((i & 2) >> 1) + ((i & 4) >> 2) + ((i & 8) >> 3)
			     + ((i & 16) >> 4) + ((i & 32) >> 5) + ((i & 64) >> 6) + ((i & 128) >> 7));
      Tmp += this->EvaluateHilbertSpaceDimension(nbrFermions - TmpNbrParticles, currentKx, currentKy - 1,
						 currentTotalKx + (TmpNbrParticles * currentKx), currentTotalKy + (TmpNbrParticles * currentKy),
						 nbrParticles12 - ((i & 1) + ((i & 2) >> 1)),
						 nbrParticles34 - (((i & 4) >> 2) + ((i & 8) >> 3)),
						 nbrParticles56 - (((i & 16) >> 4) + ((i & 32) >> 5)),
						 nbrParticles78 - (((i & 64) >> 6) + ((i & 128) >> 7)));
    }
  return Tmp;
}
