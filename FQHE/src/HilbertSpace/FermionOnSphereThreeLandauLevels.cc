////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                   Copyright (C) 2001-2005 Nicolas Regnault                 //
//                                                                            //
//                                                                            //
//                  class of fermions on sphere including three               //
//                                  Landau levels                             //
//                                                                            //
//                        last modification : 20/04/2010                      //
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
#include "HilbertSpace/FermionOnSphereThreeLandauLevels.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/ArrayTools.h"

#include <math.h>
#include <cstdlib>
#include <algorithm>

using std::cout;
using std::endl;
using std::hex;
using std::dec;


// default constructor
// 

FermionOnSphereThreeLandauLevels::FermionOnSphereThreeLandauLevels ()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// memory = amount of memory granted for precalculations

FermionOnSphereThreeLandauLevels::FermionOnSphereThreeLandauLevels (int nbrFermions, int totalLz, int lzMax, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalY = 0;
  this->TotalTz = 0;
  this->LzMax = lzMax + 4;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->Flag.Initialize();
  this->HilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax - 4, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1);
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];
  long TmpHilbertSpaceDimension = this->GenerateStates(this->NbrFermions, this->LzMax - 4, this->LzMax, (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1, 0l);
  if (TmpHilbertSpaceDimension != this->HilbertSpaceDimension)
    {
      cout << "Mismatch in State-count (" << this->HilbertSpaceDimension << ") and State Generation (" << TmpHilbertSpaceDimension << ") in FermionOnSphereThreeLandauLevels!" << endl;
      exit(1);
    }
  this->HilbertSpaceDimension = TmpHilbertSpaceDimension;
  cout << "Hilbert space dimension = " << this->HilbertSpaceDimension << endl;  
  this->GenerateLookUpTable(memory);
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
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereThreeLandauLevels::FermionOnSphereThreeLandauLevels(const FermionOnSphereThreeLandauLevels& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalTz = fermions.TotalTz;
  this->TotalY = fermions.TotalY;
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

FermionOnSphereThreeLandauLevels::~FermionOnSphereThreeLandauLevels ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      if (this->StateHighestBit != 0)
	delete[] this->StateHighestBit;
      delete[] this->LookUpTableShift;
      for (int i = 0; i < 2*this->NbrLzValue; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereThreeLandauLevels& FermionOnSphereThreeLandauLevels::operator = (const FermionOnSphereThreeLandauLevels& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->TotalTz = fermions.TotalTz;
  this->TotalY = fermions.TotalY;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
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

AbstractHilbertSpace* FermionOnSphereThreeLandauLevels::Clone()
{
  return new FermionOnSphereThreeLandauLevels(*this);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// nbrFluxQuanta = number of flux quanta
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereThreeLandauLevels::GenerateStates(int nbrFermions, int nbrFluxQuanta, int lzMax, int totalLz, long pos)
{
  if ((nbrFermions < 0) || (totalLz < 0))
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
	  if (((nbrFluxQuanta + 4) >= totalLz) && (totalLz >= 0))
	    {
	      this->StateDescription[pos] = 0x4ul << (totalLz * 3);
	      ++pos;
	    }
	  if (((nbrFluxQuanta + 3) >= totalLz) && (totalLz >= 1))
	    {
	      this->StateDescription[pos] = 0x2ul << (totalLz * 3);
	      ++pos;
	    }
	  if (((nbrFluxQuanta + 2) >= totalLz) && (totalLz >= 2))
	    {
	      this->StateDescription[pos] = 0x1ul << (totalLz * 3);
	      ++pos;
	    }
	}
      return pos;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return pos;

  long TmpPos = pos;
  unsigned long Mask = 0x0ul;
  if ((lzMax <= (nbrFluxQuanta + 4)) && (lzMax >= 0))
    {
      if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 1))
	{
	  if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
	    {
	      TmpPos = this->GenerateStates(nbrFermions - 3, nbrFluxQuanta, lzMax - 1, totalLz - (3 * lzMax), pos);
	      Mask = 0x7ul << (lzMax * 3);
	      for (; pos < TmpPos; ++pos)
		this->StateDescription[pos] |= Mask;
	    }
	  TmpPos = this->GenerateStates(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax), pos);
	  Mask = 0x6ul << (lzMax * 3);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= Mask;
	}
      if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
	{
	  TmpPos = this->GenerateStates(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax), pos);
	  Mask = 0x5ul << (lzMax * 3);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= Mask;
	}
      TmpPos = this->GenerateStates(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax, pos);
      Mask = 0x4ul << (lzMax * 3);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 1))
    {
      if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
	{
	  TmpPos = this->GenerateStates(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax), pos);
	  Mask = 0x3ul << (lzMax * 3);
	  for (; pos < TmpPos; ++pos)
	    this->StateDescription[pos] |= Mask;
	}
      TmpPos = this->GenerateStates(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax, pos);
      Mask = 0x2ul << (lzMax * 3);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))    
    {
      TmpPos = this->GenerateStates(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax, pos);
      Mask = 0x1ul << (lzMax * 3);
      for (; pos < TmpPos; ++pos)
	this->StateDescription[pos] |= Mask;
    }
  return this->GenerateStates(nbrFermions, nbrFluxQuanta, lzMax - 1, totalLz, pos);
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// nbrFluxQuanta = number of flux quanta
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereThreeLandauLevels::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int nbrFluxQuanta, int lzMax, int totalLz)
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
	  if (((nbrFluxQuanta + 4) >= totalLz) && (totalLz >= 0))
	    ++Tmp;
	  if (((nbrFluxQuanta + 3) >= totalLz) && (totalLz >= 1))
	    ++Tmp;
	  if (((nbrFluxQuanta + 2) >= totalLz) && (totalLz >= 2))
	    ++Tmp;
	}
      return Tmp;
    }

  if ((lzMax == 0) && (totalLz != 0))
    return 0l;

  long Tmp = 0l;
  if ((lzMax <= (nbrFluxQuanta + 4)) && (lzMax >= 0))
    {
      if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 1))
	{
	  if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
	    Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 3, nbrFluxQuanta, lzMax - 1, totalLz - (3 * lzMax));
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
	}
      if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
	Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (nbrFluxQuanta + 3)) && (lzMax >= 1))
    {
      if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))
	Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 2, nbrFluxQuanta, lzMax - 1, totalLz - (2 * lzMax));
      Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
    }
  if ((lzMax <= (nbrFluxQuanta + 2)) && (lzMax >= 2))    
    Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, nbrFluxQuanta, lzMax - 1, totalLz - lzMax);
  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, nbrFluxQuanta, lzMax - 1, totalLz);
  return Tmp;
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

void FermionOnSphereThreeLandauLevels::BosonicStateTimeFermionicState(RealVector& bosonState, RealVector& fermionState, RealVector& outputVector, unsigned long* finalStates, double* weigth, 
								      BosonOnSphereShort* bosonSpace, FermionOnSphere* finalSpace, int firstComponent,int nbrComponent)
{
  unsigned long* Monomial = new unsigned long[this->NbrFermions];
  unsigned long* Slater = new unsigned long[this->NbrFermions];
  int NbrMax = firstComponent+nbrComponent;
  int NbrVariable=0;
  unsigned long* Variable = new unsigned long[this->NbrFermions];
  for (int j=0; j<this->HilbertSpaceDimension;j++)
    {
      if(fermionState[j]!=0)
	{
	  this->ConvertToMonomialVariable(this->StateDescription[j], Slater,NbrVariable,Variable);
	  for (int i=firstComponent;i<NbrMax;i++)
	    {
	      if(bosonState[i]!=0)
		{
		  bosonSpace->GetMonomial(i,Monomial);
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

unsigned int FermionOnSphereThreeLandauLevels::MonomialsTimesSlaterProjection(unsigned long* slater, unsigned long* monomial, unsigned long* variable, int nbrVariable, 
									      unsigned long*& finalStates, double*& weigth, FermionOnSphere* finalSpace)
{
  unsigned int NbrNonZero=0;
  unsigned long * State = new unsigned long[this->NbrFermions];
  unsigned long TmpState=0;
  bool Bool=true;
  double C = 1.0;
  long PowerIn;
  long PowerOut;
  long Numerator;
  long AlphaIn = this->LzMax * (this->LzMax-1);
  long AlphaOut = (finalSpace->LzMax+4) * (finalSpace->LzMax+3);
  for (int i = 0; i < this->NbrFermions ; i++)
    State[i]=slater[i]+monomial[i];
  for(int k = 0 ; (k < nbrVariable) && (C != 0.0); k++)
    {
      PowerIn = (long) slater[variable[k]>>1];
      PowerOut = (long) State[variable[k]>>1];
      if((variable[k] & 0x1ul) == 0ul)
	{
	  Numerator=(PowerIn-0x1ul)*(0x2ul+finalSpace->LzMax)-(PowerOut-0x1ul)*(this->LzMax-0x2ul);
	  if(Numerator == 0x0l)
	    C = 0.0;
	  else
	    C *= ((double)Numerator / (double)((this->LzMax-0x2ul)*(0x2ul+finalSpace->LzMax)));
	}
      else
	{
	  Numerator= ((AlphaOut*PowerIn*(PowerIn-0x1ul)-AlphaIn*PowerOut*(PowerOut-0x1ul))*(0x2ul+finalSpace->LzMax)-
		       ((this->LzMax-1)*(2*PowerIn-this->LzMax)*AlphaOut-(finalSpace->LzMax+3)*(2*PowerOut-(finalSpace->LzMax+0x4ul))*AlphaIn)*(PowerOut-0x1ul));
	  if(Numerator == 0x0l)
	    C = 0.0;
	  else
	    C*=((double)Numerator/(double)(AlphaOut*(0x2ul+finalSpace->LzMax)));
	}
    }
  unsigned long Mask;
  unsigned long Sign = 0ul;
  if (C != 0.0)
    {
      for(int i = 0 ; (i < this->NbrFermions)&&(Bool); i++)
	{
	  Mask = (1ul << (State[i]-0x2ul));
	  if ( (TmpState&Mask) != 0ul)
	    Bool = false;
	  unsigned long TmpState2 = TmpState&(Mask-0x1ul);
	  #ifdef _64_BITS__
	  TmpState2 ^= TmpState2 >> 32;
	  #endif
	  TmpState2 ^= TmpState2 >> 16;
	  TmpState2 ^= TmpState2 >> 8;
	  TmpState2 ^= TmpState2 >> 4;
	  TmpState2 ^= TmpState2 >> 2;
	  TmpState2 ^= TmpState2 >> 1;
	  Sign ^= TmpState2;
	  TmpState|=Mask;
	}
      if(Bool)
	{
	  if ((Sign & 0x1ul) == 1ul) C*=-1.0;
	  NbrNonZero += SearchInArrayAndSetWeight(TmpState, finalStates, weigth, NbrNonZero, C);
	}
    }
  while (std::prev_permutation(monomial,monomial+this->NbrFermions))
    {
      C = 1.0;
      for(int i = 0 ; i<this->NbrFermions;i++)
	{
	  State[i] = slater[i] + monomial[i];
	}
      for(int k = 0 ; (k < nbrVariable) && (C != 0.0); k++)
	{
	  PowerIn = (long) slater[variable[k]>>1];
	  PowerOut = (long) State[variable[k]>>1];
	  if((variable[k] & 0x1ul) == 0ul)
	    {
	      Numerator = (PowerIn-0x1ul)*(0x2ul+finalSpace->LzMax) - (PowerOut-0x1ul)*(this->LzMax-0x2ul);
	      if(Numerator == 0x0l)
		C=0.0;
	      else
		C *= ((double)Numerator/(double)((this->LzMax-0x2ul)*(0x2ul+finalSpace->LzMax)));
	    }
	  else
	    {
	      Numerator= ((AlphaOut*PowerIn*(PowerIn-0x1ul)-AlphaIn*PowerOut*(PowerOut-0x1ul))*(0x2ul+finalSpace->LzMax) - 
			   ((this->LzMax-1)*(2*PowerIn-this->LzMax)*AlphaOut-(finalSpace->LzMax+3)*(2*PowerOut-(finalSpace->LzMax+0x4ul))*AlphaIn)*(PowerOut-0x1ul));
	      if(Numerator==0x0l)
		C = 0.0;
	      else
		C *= ((double)Numerator / (double)(AlphaOut*(2+finalSpace->LzMax)));
	    }
	}
      if (C != 0.0)
	{
	  Bool = true;
	  TmpState = 0;
	  Sign = 0ul;
	  for (int i=0; (i < this->NbrFermions)&&(Bool);i++)
	    {
	      Mask = (1ul << (State[i]-0x2ul));
	      if((TmpState&Mask) != 0)
		Bool=false;
	      unsigned long TmpState2 = TmpState&(Mask-0x1ul);
#ifdef  __64_BITS__
	      TmpState2 ^= TmpState2 >> 32;
#endif
	      TmpState2 ^= TmpState2 >> 16;
	      TmpState2 ^= TmpState2 >> 8;
	      TmpState2 ^= TmpState2 >> 4;
	      TmpState2 ^= TmpState2 >> 2;
	      TmpState2 ^= TmpState2 >> 1;
	      Sign ^= TmpState2;
	      TmpState|=Mask;
	    }
	  if(Bool)
	    {
	      if ((Sign & 0x1ul) == 1ul) C*=-1.0;
	      NbrNonZero += SearchInArrayAndSetWeight(TmpState,finalStates,weigth,NbrNonZero,C);
	    }
	}
    }
  delete [] State;
  return NbrNonZero;
}
