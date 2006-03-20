////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on sphere without total Lz restriction         //
//                                                                            //
//                        last modification : 16/03/2006                      //
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
#include "HilbertSpace/QHEHilbertSpace/FullBosonOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SpinQuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"

#include <math.h>


using std::cout;
using std::endl;


// basic constructor
// 
// nbrBosons = number of bosons
// lzMax = maximum Lz value reached by a boson

FullBosonOnSphere::FullBosonOnSphere (int nbrBosons, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->LzMax);
  this->TemporaryState = new int [this->NbrLzValue];
  this->ProdATemporaryState = new int [this->NbrLzValue];
  this->Flag.Initialize();
  this->StateDescription = new int* [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];
  this->StateTotalLz = new int [this->HilbertSpaceDimension];
  this->GenerateStates(this->NbrBosons, this->LzMax, this->LzMax, this->ShiftedTotalLz, 0);
  this->KeyMultiplicationTable = new int [this->LzMax + 1];
  this->GenerateLookUpTable(0);
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;

#ifdef __DEBUG__
  int UsedMemory = 0;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    UsedMemory += (this->StateLzMax[i] + 1) * sizeof(int) + sizeof(int*);
  UsedMemory += (this->TotalLz + 1) * sizeof(int);
  UsedMemory += this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += ((this->TotalLz + 1) * this->IncNbrBosons) * sizeof(int);
  UsedMemory += this->HilbertSpaceDimension * sizeof(int);
/*  cout << "memory requested for Hilbert space = ";
  if (UsedMemory >= 1024)
    if (UsedMemory >= 1048576)
      cout << (UsedMemory >> 20) << "Mo" << endl;
    else
      cout << (UsedMemory >> 10) << "ko" <<  endl;
  else
    cout << UsedMemory << endl;*/
#endif
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

FullBosonOnSphere::FullBosonOnSphere(const FullBosonOnSphere& bosons)
{
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->StateTotalLz = bosons.StateTotalLz;
  this->Flag = bosons.Flag;
  this->Keys = bosons.Keys;
  this->KeyMultiplicationTable = bosons.KeyMultiplicationTable;
  this->LzMaxPosition = bosons.LzMaxPosition;
  this->KeyInvertSectorSize = bosons.KeyInvertSectorSize;
  this->KeyInvertTable = bosons.KeyInvertTable;
  this->KeyInvertTableNbrIndices = bosons.KeyInvertTableNbrIndices;
  this->KeyInvertIndices = bosons.KeyInvertIndices;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new int [this->NbrLzValue];
  this->ProdATemporaryState = new int [this->NbrLzValue];
}

// destructor
//

FullBosonOnSphere::~FullBosonOnSphere ()
{
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Keys;
      delete[] this->KeyMultiplicationTable;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->StateTotalLz;
      int Size = (this->LzMax + 2) * this->IncNbrBosons;
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
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
      delete[] LzMaxPosition;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FullBosonOnSphere& FullBosonOnSphere::operator = (const FullBosonOnSphere& bosons)
{
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->Keys;
      delete[] this->KeyMultiplicationTable;
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	delete[] this->StateDescription[i];
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
      delete[] this->StateTotalLz;
      int Size = (this->LzMax + 2) * this->IncNbrBosons;
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
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
    }
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->LzMax = bosons.LzMax;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->StateTotalLz = bosons.StateTotalLz;
  this->Flag = bosons.Flag;
  this->Keys = bosons.Keys;
  this->KeyMultiplicationTable = bosons.KeyMultiplicationTable;
  this->LzMaxPosition = bosons.LzMaxPosition;
  this->KeyInvertSectorSize = bosons.KeyInvertSectorSize;
  this->KeyInvertTable = bosons.KeyInvertTable;
  this->KeyInvertTableNbrIndices = bosons.KeyInvertTableNbrIndices;
  this->KeyInvertIndices = bosons.KeyInvertIndices;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new int [this->NbrLzValue];
  this->ProdATemporaryState = new int [this->NbrLzValue];

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FullBosonOnSphere::Clone()
{
  return new FullBosonOnSphere(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FullBosonOnSphere::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FullBosonOnSphere::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FullBosonOnSphere::ExtractSubspace (AbstractQuantumNumber& q, 
							  SubspaceSpaceConverter& converter)
{
  return 0;
}

// find state index
//
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// totalLz = state total Lz value
// return value = corresponding index

int FullBosonOnSphere::FindStateIndex(int* stateDescription, int lzmax, int totalLz)
{
  int TmpKey = this->GenerateKey(stateDescription, lzmax);
  int Sector = lzmax * this->IncNbrBosons + stateDescription[lzmax];
  int TmpPos = 0;
  int TmpPos2 = this->KeyInvertSectorSize[Sector] - 1;
  int TmpPos3;
  int* TmpKeyInvertTable = this->KeyInvertTable[totalLz][Sector];
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
  int* TmpKeyInvertIndices = this->KeyInvertIndices[totalLz][Sector][TmpPos];

  TmpPos2 =this->KeyInvertTableNbrIndices[totalLz][Sector][TmpPos] - 1;
  TmpPos = 0;
  while (TmpPos2 != TmpPos)
    {
      TmpPos3 = (TmpPos2 + TmpPos) >> 1;
      Start = TmpKeyInvertIndices[TmpPos3];
      TmpStateDescription = this->StateDescription[Start];
      i = lzmax;
      while ((i >= 0) && (stateDescription[i] ==  TmpStateDescription[i]))
        --i;
      if (i == -1)
        {
          return Start;
        }
      if (stateDescription[i] < TmpStateDescription[i])
        {
          TmpPos = TmpPos3 + 1;
        }
      else
        {
           TmpPos2 = TmpPos3 - 1;
        }
    }
  return TmpKeyInvertIndices[TmpPos];
}

// generate all states
// 
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson in the state
// return value = Hilbert space size

int FullBosonOnSphere::GenerateStates(int nbrBosons, int lzMax)
{
  int TmpMaxLz =  nbrBosons * lzmax;
  int Pos = 0;
  int i = 0;
  if (((nbrBosons & 1) == 1) && ((lzmax & 1) == 1))
    i = 1; 
  for (; i <= TmpMaxLz; i += 2)
    {
      int TmpPos = this->GenerateStates(nbrBosons, lzMax, lzMax, i, Pos);
      for (int j = Pos; j < TmpPos; ++j)
	this->StateTotalLz[j] = Pos;
      Pos = TmpPos;
    }
  return TmpPos;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FullBosonOnSphere::GenerateLookUpTable(int memory)
{
  this->Keys = new int [this->HilbertSpaceDimension];
  for (int i = 0; i <= this->LzMax; ++i)
    this->KeyMultiplicationTable[i] = i * i * i * i;//this->IncNbrBosons;
  int Size = (this->LzMax + 2) * this->IncNbrBosons;
  this->LzMaxPosition = new int [Size];
  this->KeyInvertSectorSize = new int [Size];
  this->KeyInvertTable = new int* [Size];
  this->KeyInvertTableNbrIndices = new int* [Size];
  this->KeyInvertIndices = new int** [Size];
  for (int i = 0; i < Size; ++i)
    this->KeyInvertSectorSize[i] =0;
  int CurrentLzMax = this->LzMax;
  int CurrentNbrLzMax = 1;
  this->LzMaxPosition[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 0; 
  this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 0;   
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->Keys[i] = this->GenerateKey(this->StateDescription[i], this->StateLzMax[i]);
      if (CurrentLzMax != this->StateLzMax[i])
	{
	  CurrentLzMax = this->StateLzMax[i];
	  CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	  this->LzMaxPosition[CurrentLzMax * this->IncNbrBosons + CurrentNbrLzMax] = i; 
	  this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	}
      else
	if (this->StateDescription[i][CurrentLzMax] != CurrentNbrLzMax)
	  {
	    CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	    this->LzMaxPosition[CurrentLzMax * this->IncNbrBosons + CurrentNbrLzMax] = i;
	    this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	  }
	else
	  {
	    ++this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax]; 
	  }
    }
  
  for (int i = 0; i < Size; ++i)
    {
      if (this->KeyInvertSectorSize[i] > 0)
	{
	  this->KeyInvertTable[i] = new int [this->KeyInvertSectorSize[i]];
	  this->KeyInvertTableNbrIndices[i] = new int [this->KeyInvertSectorSize[i]];
	}
      else
	{
	  this->KeyInvertTable[i] = 0; 
	  this->KeyInvertTableNbrIndices[i] = 0;
	}
    }

  CurrentLzMax = this->StateLzMax[0];
  CurrentNbrLzMax = this->StateDescription[0][CurrentLzMax];
  int CurrentKeyInvertSectorSize = this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
  int* TmpKeyInvertTable = this->KeyInvertTable[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  int* TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  TmpKeyInvertTable[0] = this->Keys[0];
  TmpKeyInvertTableNbrIndices[0] = 1;
  for (int i = 1; i < this->HilbertSpaceDimension; ++i)
    {
      if (CurrentLzMax != this->StateLzMax[i])
	{
	  CurrentLzMax = this->StateLzMax[i];
	  CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	  CurrentKeyInvertSectorSize = this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	  this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	  TmpKeyInvertTable = this->KeyInvertTable[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	  TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	  TmpKeyInvertTable[0] = this->Keys[i];
	  TmpKeyInvertTableNbrIndices[0] = 1;
	}
      else
	if (this->StateDescription[i][CurrentLzMax] != CurrentNbrLzMax)
	  {
	    CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	    CurrentKeyInvertSectorSize = this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	    this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax] = 1;
	    TmpKeyInvertTable = this->KeyInvertTable[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	    TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	    TmpKeyInvertTable[0] = this->Keys[i];
	    TmpKeyInvertTableNbrIndices[0] = 1;
	  }
	else
	  {
	    int Lim = this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	    int j = 0;
	    bool Flag = false;
	    int TmpKey = this->Keys[i];
	    for (; ((j < Lim) && (Flag == false)); ++j)
	      if (TmpKeyInvertTable[j] == TmpKey)
		{
		  Flag = true;
		  --j;
		}
	    if (Flag == true)
	      {
		++TmpKeyInvertTableNbrIndices[j];
	      }
	    else
	      {
		TmpKeyInvertTable[this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax]] = TmpKey;
		TmpKeyInvertTableNbrIndices[this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax]] = 1;
		++this->KeyInvertSectorSize[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
	      }
	  }
    }
 
  // sort key invert table and resize arrays
  int TmpPos;
  int TmpSize;
  for (int i = 0; i < Size; ++i)
    {
      if (this->KeyInvertSectorSize[i] > 0)
	{
	  int Lim = this->KeyInvertSectorSize[i];
	  TmpKeyInvertTable = this->KeyInvertTable[i];
	  TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[i];
	  int* TmpKeyInvertTable2 = new int [Lim];
	  int* TmpKeyInvertTableNbrIndices2 = new int [Lim];
	  int Tmp;
	  for (int j = 0; j < Lim; ++j)
	    {
	      Tmp = TmpKeyInvertTable[j];
	      TmpPos = j;
	      for (int k = j + 1; k < Lim; ++k)
		if (Tmp > TmpKeyInvertTable[k])
		  {
		    Tmp = TmpKeyInvertTable[k];
		    TmpPos = k;
		  }
	      TmpKeyInvertTable2[j] = Tmp;
	      TmpKeyInvertTableNbrIndices2[j] = TmpKeyInvertTableNbrIndices[TmpPos];
	      TmpKeyInvertTableNbrIndices[TmpPos] = TmpKeyInvertTableNbrIndices[j];
	      TmpKeyInvertTable[TmpPos] = TmpKeyInvertTable[j];
	    }
	  delete[] TmpKeyInvertTable;
	  delete[] TmpKeyInvertTableNbrIndices;
	  this->KeyInvertTable[i] = TmpKeyInvertTable2;
	  this->KeyInvertTableNbrIndices[i] = TmpKeyInvertTableNbrIndices2;
	  this->KeyInvertIndices[i] = new int* [Lim];
	  for (int j = 0; j < Lim; ++j)
	    {
	      this->KeyInvertIndices[i][j] = new int [this->KeyInvertTableNbrIndices[i][j]];
	      this->KeyInvertTableNbrIndices[i][j] = 0;
	    }
	}
    }

  // find all hilbert space indices that have the same key in each sector
  CurrentLzMax = this->LzMax;
  CurrentNbrLzMax = 1;
  TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  int** TmpKeyInvertIndices = this->KeyInvertIndices[CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax];
  int TmpPos2;
  int TmpPos3;
  int TmpPos4 = CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax;
  TmpSize = this->KeyInvertSectorSize[TmpPos4];
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      int TmpKey = this->Keys[i];
      if (CurrentLzMax != this->StateLzMax[i])
	{
	  CurrentLzMax = this->StateLzMax[i];
	  CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	  TmpPos4 = CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax;
	  TmpKeyInvertIndices = this->KeyInvertIndices[TmpPos4];
	  TmpKeyInvertTable = this->KeyInvertTable[TmpPos4];
	  TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[TmpPos4];
	  TmpSize = this->KeyInvertSectorSize[TmpPos4];
	  TmpPos = 0;
	  TmpPos2 = TmpSize - 1;
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
	  TmpKeyInvertIndices[TmpPos][0] = i;
	  TmpKeyInvertTableNbrIndices[TmpPos] = 1;
	}
      else
	if (this->StateDescription[i][CurrentLzMax] != CurrentNbrLzMax)
	  {
	    CurrentNbrLzMax = this->StateDescription[i][CurrentLzMax];
	    TmpPos4 = CurrentLzMax * (this->IncNbrBosons) + CurrentNbrLzMax;
	    TmpKeyInvertIndices = this->KeyInvertIndices[TmpPos4];
	    TmpKeyInvertTable = this->KeyInvertTable[TmpPos4];
	    TmpKeyInvertTableNbrIndices = this->KeyInvertTableNbrIndices[TmpPos4];
	    TmpSize = this->KeyInvertSectorSize[TmpPos4];
	    TmpPos = 0;
	    TmpPos2 = TmpSize - 1;
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
	    TmpKeyInvertIndices[TmpPos][0] = i;
	    TmpKeyInvertTableNbrIndices[TmpPos] = 1;
	  }
	else
	  {
	    TmpPos = 0;
	    TmpPos2 = TmpSize - 1;
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
	    TmpKeyInvertIndices[TmpPos][TmpKeyInvertTableNbrIndices[TmpPos]] = i;
	    ++TmpKeyInvertTableNbrIndices[TmpPos];
	  }
    }

}

// generate look-up table associated to current Hilbert space
// 
// stateDescription = array describing the state
// lzmax = maximum Lz value reached by a boson in the state
// return value = key associated to the state

int FullBosonOnSphere::GenerateKey(int* stateDescription, int lzmax)
{
  int Key = 0;
  for (int i = 0; i <= lzmax; ++i)
    {
      Key += this->KeyMultiplicationTable[i] * stateDescription[i];
    }
  return Key;
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a boson
// totalLz = momentum total value
// return value = Hilbert space dimension

int FullBosonOnSphere::EvaluateHilbertSpaceDimension(int nbrBosons, int lzMax)
{
  if (((nbrBosons & 1) != 0) && ((lzMax & 1) != 0))
    return ( >> 1);
  else
    return (( + this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax, (lzMax * nbrBosons) >> 1)) >> 1);
}



