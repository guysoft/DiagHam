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
  if (this->LzMaxUp >= this->LzMaxDown) //we assume this is the case and take it that the lower Landau level is labelled down.
    {
      this->LzMax = this->LzMaxUp;
      //not sure of following shifts. Everything should be shifted by lzMaxUp 
      this->LzShiftUp = 0; 
      this->LzShiftDown = (this->LzMaxUp - this->LzMaxDown) >> 1;
    }
  else 
    {
      cout << "lzMaxUp should be greater than lzMaxdown!" << endl;
      this->LzMax = this->LzMaxDown;
      this->LzShiftDown = 0;
      this->LzShiftUp = (this->LzMaxDown - this->LzMaxUp) >> 1;
    }
      
  this->UpStateShift = nbrBosons + lzMaxDown;
  this->LzTotalShift = this->LzMaxDown + this->LzMaxUp;
  this->NbrLzValue = this->LzMaxUp + this->LzMaxDown + 2;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->ShiftedEvaluateHilbertSpaceDimension(this->NbrBosons, 0, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  this->Flag.Initialize();
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateLzMax = new int [this->HilbertSpaceDimension];  //this will store two values, the index of the split beween Landau levels and the maximum lowest bit used. 
  this->TemporaryState = new unsigned long [this->NbrLzValue]; //to keep bosonic description in
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue]; //to keep bosonic description in
  int TmpDimension = this->GenerateStates(this->NbrBosons, 0, (this->TotalLz + (this->NbrBosons * this->LzMax)) >> 1, 0);
  if (TmpDimension != this->HilbertSpaceDimension)
    {
      cout << "Mismatch in State-count and State Generation in BosonOnSphereTwoLandauLevels! " << this->HilbertSpaceDimension << " " << TmpDimension  << endl;
  for (int i = 0; i < TmpDimension; ++i)
    this->PrintState(cout, i) << endl;
       exit(1);
    }

  //  this->HilbertSpaceDimension = this->GenerateStates(this->NbrBosonsUp, this->NbrBosonsDown, this->LzMaxUp, this->LzMaxDown, );
  //this->GenerateLookUpTable(memory);
  
#ifdef __DEBUG__
  this->LookUpTableMemorySize = 0;
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
  this->StateLzMax = bosons.StateLzMax;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;  
  this->SignLookUpTable = bosons.SignLookUpTable;
  this->SignLookUpTableMask = bosons.SignLookUpTableMask;
  this->MaximumSignLookUp = bosons.MaximumSignLookUp;
  this->ProdATemporaryState = bosons.ProdATemporaryState;
  this->TemporaryState = bosons.TemporaryState;
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
      delete[] this->StateLzMax;
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
  this->StateLzMax = bosons.StateLzMax;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;  
  this->SignLookUpTable = bosons.SignLookUpTable;
  this->SignLookUpTableMask = bosons.SignLookUpTableMask;
  this->MaximumSignLookUp = bosons.MaximumSignLookUp;
  this->ProdATemporaryState = bosons.ProdATemporaryState;
  this->TemporaryState = bosons.TemporaryState;
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
  unsigned long Value = this->StateDescription[state];
  Str << this->StateDescription[state] << ":\t";
#ifdef  __64_BITS__
  const int SHIFT = 63;
#else
  const int SHIFT = 31;
#endif
  const unsigned long MASK = 0x1ul << SHIFT;

   for ( int i = 1; i <= SHIFT + 1; i++ ) 
   {
     if ( (Value & MASK) > 0 ) 
       {
	 Str << '1';
       }
     else
       {
         Str << '0';
       }
     Value <<= 1;

      if ( i % 8 == 0 )
         Str << ' ';
   }
  Str << "  " << this->StateLzMax[state] << endl;
  
  Value = this->StateDescription[state];
  /*int occupations[this->LzMaxUp*2];
  for ( int i = 0 ; i < this->LzMaxUp*2 ; i++ )
    {
      occupations[i] = 0;
    }
  
  int bit = 63;
  int pos = 0;
  unsigned long mask;
  bool last = false;
  while ( bit >= 0 && pos < this->LzMaxUp*2 && bit >= (63 - this->StateLzMax[state]) ) 
    {
      mask = 0x1ul << bit; 
      if ( (mask & value) > 0 ) 
        {
	  occupations[pos]++;
	  last = true;
	}	
      else 
        {
	  if ( last ) last = false;
	  pos++;
	}
      bit--;
    }
  
  for ( int i = 0 ; i <= this->LzMaxUp ; i++ ) {
      Str << occupations[this->LzMaxUp-i] << " " ;
  }
  Str << endl;
  Str << "  ";
  for ( int i = 0 ; i <= this->LzMaxDown ; i++ ) {
      Str << occupations[this->LzMaxUp + 1 + LzMaxDown - i] << " " ;
  }*/
  
  this->FermionToBoson(Value, this->StateLzMax[state], this->TemporaryState, this->TemporaryStateLzMax);
  for ( int i = 0 ; i <= this->LzMaxUp ; i++ ) {
  	  if ( (this->LzMaxUp-i) <=  this->TemporaryStateLzMax ) {
	      Str << this->TemporaryState[this->LzMaxUp-i] << " " ;
	  } else {
		  Str << "0" << " " ;
	  }
  }
  Str << endl;
  Str << "  ";
  for ( int i = 0 ; i <= this->LzMaxDown ; i++ ) {
      if ( (this->LzMaxUp + 1 + LzMaxDown - i) <=  this->TemporaryStateLzMax ) 
      {
      	Str << this->TemporaryState[this->LzMaxUp + 1 + LzMaxDown - i] << " " ;
      } else {
      	Str << "0" << " " ;
      }
  }  
  
  Str << endl;
  Value = this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax);
  
   for ( int i = 1; i <= SHIFT + 1; i++ ) 
   {
     if ( (Value & MASK) > 0 ) 
       {
	 Str << '1';
       }
     else
       {
         Str << '0';
       }
     Value <<= 1;

      if ( i % 8 == 0 )
         Str << ' ';
   }
  Str << "  " << this->StateLzMax[state] << ", index: " << this->FindStateIndex(this->StateDescription[state]) << endl;
  
  cout << "Norm: " <<  this->GetConfigNorm(state) << endl;
  
  this->PrintStateMonomial(Str, state) << endl;
  
  unsigned long* MonomialRep;
  MonomialRep = new unsigned long[this->NbrBosons];
  this->GetMonomial((long)state, MonomialRep);
  Str << endl <<  "Converted back from monomial: " << this->ConvertFromMonomial(MonomialRep) << endl;
  delete []MonomialRep;

  return Str;
}


// print a given State
//
// Str = reference on current output stream 
// state = binary representation of state to print
// return value = reference on current output stream 

ostream& BosonOnSphereTwoLandauLevels::PrintStateBinary (ostream& Str, unsigned long state)
{
  unsigned long Value = state;
#ifdef  __64_BITS__
  const int SHIFT = 63;
#else
  const int SHIFT = 31;
#endif
  const unsigned long MASK = 0x1ul << SHIFT;

   for ( int i = 1; i <= SHIFT + 1; i++ ) 
   {
     if ( (Value & MASK) > 0 ) 
       {
	 Str << '1';
       }
     else
       {
         Str << '0';
       }
     Value <<= 1;

      if ( i % 8 == 0 )
         Str << ' ';
   }

  return Str;
}

// print monomial representation of a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereTwoLandauLevels::PrintStateMonomial (ostream& Str, int state)
{
  unsigned long* MonomialRep;
  MonomialRep = new unsigned long[this->NbrBosons];
  this->GetMonomial((long)state, MonomialRep);
  for ( int i = 0; i < this->NbrBosons; i++ )
    {
	if ( i > 0 ) Str << ", ";
	Str << MonomialRep[i] ;
    }
  delete []MonomialRep;
  return Str;
}

// evaluate Hilbert space dimension without constraint on the number of particles per level
//
// nbrBosons = number of bosons
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension
long BosonOnSphereTwoLandauLevels::ShiftedEvaluateHilbertSpaceDimension(int nbrBosons, int lzMax, int totalLz)
{
  if ((nbrBosons == 0) && (totalLz == 0)) //all bosons gone and correct total totalLz been realised. Fill in correct entries on way back up.
    {
      return 1l;
    }
    
  if (nbrBosons < 0 || totalLz < 0  )//|| (this->MaxLzLeft(nbrBosons,pos) < totalLz) ) //if state not working out. 
    return 0l;
    
  if (lzMax < 0) //if the position is negative. This should never happen.
    return 0l;
  
  if ((lzMax == (this->LzMaxUp + this->LzMaxDown + 2)) && (totalLz != 0)) //if at the final position and still hav not satisfied the totallz.
    return 0l;
  
  int currentLz = 0; // this is the lz value of the current position.
  if ( lzMax <= this->LzMaxUp )  //if its on the Up (second) Landau level.
    {
      currentLz = this->LzMaxUp - lzMax; 
    }
  else //if its on the down (lowest) Landau level
    {
      currentLz = this->LzMaxDown + 1 - (lzMax - (this->LzMaxUp+1));
    }
    
  long Tmp = 0;
  if ( lzMax < (this->LzMaxUp + this->LzMaxDown + 2) ) //if the position has not reached the end.
    {
      //Can place between 0 and nbrBosons on this spot and then move on.
      for ( int i = nbrBosons ; i >= 1 ; i-- ) 
        {
	  //place i bosons
	  Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons - i, lzMax+1, totalLz - (currentLz*i));
	}
	Tmp += this->ShiftedEvaluateHilbertSpaceDimension(nbrBosons, lzMax+1, totalLz);
    }
    return Tmp;
}


// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = this is the possition offset from the max lz on second Landau level. Goes from 0 to (lzMaxUp + lzMaxdown + 1)
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long BosonOnSphereTwoLandauLevels::GenerateStates(int nbrBosons, int lzMax, int totalLz, long pos)
{
  /* 
    This recursive function generates all the bosonic states with nbrBosons over two Landau levels 
    with totalLz. The lowest Landau level has (this->lzMaxDown + 1) places with shifted Lz going from 1 to (this->lzMaxDown + 1).
    The second Landau level has (this->lzMaxUp + 1) places with shifted Lz going from 0 to (this->lzMaxUp)
  
    Each state is stored as an unsigned long integer and should be generated and stored in decreasing numerical order. 
    The bits of the unsigned long integer are split in half and the most significant bits used for the second Landau level while the 
    least significant bits are used for the lowest Landau level. 
    
    The configuration on each Landau level within each subset is stored starting from the most significant bit corresponding to the 
    largest Lz in that Landau level. 
  */
  
  if ((nbrBosons == 0) && (totalLz == 0)) //all bosons gone and correct total totalLz been realised. Fill in correct entries on way back up.
    {
      this->StateDescription[pos] = 0x0ul;
      this->StateLzMax[pos] = ((lzMax + this->NbrBosons - 1 - 1 ) );
      return (pos + 1l);
    }
    
  if (nbrBosons < 0 || totalLz < 0  )//|| (this->MaxLzLeft(nbrBosons,pos) < totalLz) ) //if state not working out. 
    return pos;
    
  if (lzMax < 0) //if the position is negative. This should never happen.
    return pos;
  
  if ((lzMax == (this->LzMaxUp + this->LzMaxDown + 2)) && (totalLz != 0)) //if at the final position and still hav not satisfied the totallz.
    return pos;
  
  int currentLz = 0; // this is the lz value of the current position.
  if ( lzMax <= this->LzMaxUp )  //if its on the Up (second) Landau level.
    {
      currentLz = this->LzMaxUp - lzMax; 
    }
  else //if its on the down (lowest) Landau level
    {
      currentLz = this->LzMaxDown + 1 - (lzMax - (this->LzMaxUp+1));
    }
    
  long TmpPos = pos;
  if ( lzMax < (this->LzMaxUp + this->LzMaxDown + 2) ) //if the position has not reached the end.
    {
      //Can place between 0 and nbrBosons on this spot and then move on.
      for ( int i = nbrBosons ; i >= 1 ; i-- ) 
        {
	  //place i bosons
	  TmpPos = this->GenerateStates(nbrBosons - i, lzMax+1, totalLz - (currentLz*i),  pos);
	  unsigned long Mask = 0l;
	  unsigned long bosons = (0x1ul << i) - 1;
#ifdef  __64_BITS__
	  Mask = (bosons << (63 - (i - 1 ) - (lzMax + this->NbrBosons - nbrBosons) )); //this places bosons starting from the MSB.
#else
	  Mask = (bosons << (31 - (i - 1 ) - (lzMax + this->NbrBosons - nbrBosons) )); //this places bosons starting from the MSB.
#endif	  
	  for (; pos < TmpPos; ++pos)
	    {
	      this->StateDescription[pos] |= Mask;
	    }
	  pos = TmpPos;
	}
	return this->GenerateStates(nbrBosons, lzMax+1, totalLz,  TmpPos);
    }
    return TmpPos;
}

// works out the maximum possible totallz that is left
//
// nbrBosons = the number of bosons that are left
// pos = the index of the position we are on in filling with bosons where 0 is the largest Lz on the second Landau level
// return value = maximum possible totallz that can be result

long BosonOnSphereTwoLandauLevels::MaxLzLeft(int nbrBosons, int pos) {
  long MaxLzLeft = 0; 
  if ( pos > this->LzMaxUp ) //if its on the first Landau level 
    {
      MaxLzLeft = nbrBosons * (this->LzMaxDown + 1 - (pos - this->LzMaxUp));
    }
  else //else if tis on the second Landau level
    {
      MaxLzLeft = nbrBosons * (this->LzMaxUp - pos);
    }
  return MaxLzLeft;
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
// void BosonOnSphereTwoLandauLevels::GenerateLookUpTable(unsigned long memory)
// {
// }

// find state index
//
// stateDescription = unsigned integer describing the state
// return value = corresponding index

int BosonOnSphereTwoLandauLevels::FindStateIndex(unsigned long stateDescription) 
{
	int start, end, mid;
	
	start = 0;					//index of start of range
	end = this->HilbertSpaceDimension;		//index of end of range + 1 
	
	while ( (end - start) > 0 ) 
	  {
	  	mid = (start + end) >> 1 ; 		//work out the mid-point
	  	if ( stateDescription > this->StateDescription[mid] )	
	  	  {
	  	  	end = mid;
	  	  }
	  	else if ( stateDescription < this->StateDescription[mid] )
	  	  {
	  	  	start = mid + 1;	  	 
	  	  }
	  	else
	  	  {
	  	  	return mid;
	  	  }
	  }	
	return this->HilbertSpaceDimension;
}


// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = Lz of first operator on LLL
// n2 = Lz of second operator on LLL
// return value =  multiplicative factor 

double BosonOnSphereTwoLandauLevels::AdAd (int index, int n1, int n2)
{
  int n1_index, n2_index;
  n1_index = this->GetIndexFromLzD(n1);
  n2_index = this->GetIndexFromLzD(n2);
  this->FermionToBoson(this->StateDescription[index],this->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
  if ((n1_index > this->ProdATemporaryStateLzMax) || (n2_index > this->ProdATemporaryStateLzMax) || 
      (this->ProdATemporaryState[n1_index] == 0) || (this->ProdATemporaryState[n2_index] == 0) || ((n1_index == n2_index) && (this->ProdATemporaryState[n1_index] == 1)))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryState[n2_index];
  --this->ProdATemporaryState[n2_index];
  Coefficient *= this->ProdATemporaryState[n1_index];
  --this->ProdATemporaryState[n1_index];
  for (int i = this->ProdATemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0ul;
  return sqrt(Coefficient);
}


// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = Lz of first operator on second LL
// n2 = Lz of second operator on second LL
// return value =  multiplicative factor 

double BosonOnSphereTwoLandauLevels::AuAu (int index, int n1, int n2)
{
  int n1_index, n2_index;
  n1_index = this->GetIndexFromLzU(n1);
  n2_index = this->GetIndexFromLzU(n2);
  this->FermionToBoson(this->StateDescription[index],this->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
  if ((n1_index > this->ProdATemporaryStateLzMax) || (n2_index > this->ProdATemporaryStateLzMax) || 
      (this->ProdATemporaryState[n1_index] == 0) || (this->ProdATemporaryState[n2_index] == 0) || ((n1_index == n2_index) && (this->ProdATemporaryState[n1_index] == 1)))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryState[n2_index];
  --this->ProdATemporaryState[n2_index];
  Coefficient *= this->ProdATemporaryState[n1_index];
  --this->ProdATemporaryState[n1_index];
  for (int i = this->ProdATemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0ul;
  return sqrt(Coefficient);
}


// apply a_n1_u a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = Lz of first operator on second LL
// n2 = Lz of second operator on LLL
// return value =  multiplicative factor 

double BosonOnSphereTwoLandauLevels::AuAd (int index, int n1, int n2)
{
  int n1_index, n2_index;
  n1_index = this->GetIndexFromLzU(n1);
  n2_index = this->GetIndexFromLzD(n2);
  this->FermionToBoson(this->StateDescription[index],this->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
  if ((n1_index > this->ProdATemporaryStateLzMax) || (n2_index > this->ProdATemporaryStateLzMax) || 
      (this->ProdATemporaryState[n1_index] == 0) || (this->ProdATemporaryState[n2_index] == 0) || ((n1_index == n2_index) && (this->ProdATemporaryState[n1_index] == 1)))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryState[n2_index];
  --this->ProdATemporaryState[n2_index];
  Coefficient *= this->ProdATemporaryState[n1_index];
  --this->ProdATemporaryState[n1_index];
  for (int i = this->ProdATemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0ul;
  return sqrt(Coefficient);
}


// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereTwoLandauLevels::AduAdu (int m1, int m2, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  int m1_index = this->GetIndexFromLzU(m1);
  int m2_index = this->GetIndexFromLzU(m2);
  ++this->TemporaryState[m2_index];
  coefficient = this->TemporaryState[m2_index];
  ++this->TemporaryState[m1_index];
  coefficient *= this->TemporaryState[m1_index];
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->NbrLzValue - 1;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax));
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnSphereTwoLandauLevels::AddAdd (int m1, int m2, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  int m1_index = this->GetIndexFromLzD(m1);
  int m2_index = this->GetIndexFromLzD(m2);
  ++this->TemporaryState[m2_index];
  coefficient = this->TemporaryState[m2_index];
  ++this->TemporaryState[m1_index];
  coefficient *= this->TemporaryState[m1_index];
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->NbrLzValue - 1;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax));
}

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it) //
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int BosonOnSphereTwoLandauLevels::AduAdd (int m1, int m2, double& coefficient)
{
   for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  int m1_index = this->GetIndexFromLzU(m1);
  int m2_index = this->GetIndexFromLzD(m2);
  ++this->TemporaryState[m2_index];
  coefficient = this->TemporaryState[m2_index];
  ++this->TemporaryState[m1_index];
  coefficient *= this->TemporaryState[m1_index];
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->NbrLzValue - 1;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax));
}

// calculate the pseudo potentials denoted V^S_J in the literature
//
// S double the max angular momentum
// J the relative angular momentum
// return value = The pseudo potential V^S_J

double BosonOnSphereTwoLandauLevels::CalculatePseudoPotential(int S, int J) 
{
   return  ( (double)this->NChooseC(4*S - 2*J,2*S - J) / (double)this->NChooseC(4*S + 2, 2*S + 1) ) * ( (double)this->NChooseC(4*S + 2*J + 2, 2*S + J + 1) / (double)this->NChooseC(4*S + 2, 2*S + 1) );
}


// calculate the number of ways of choosing c elements from n
// 
// n the number of options
// c the number of choices
// return value = number of ways of choosing c elements from n

unsigned long BosonOnSphereTwoLandauLevels::NChooseC(int n , int c ) 
{
    if ( c == n || c == 0 ) return 1ul;
    if ( c == 1ul ) return n;
    return this->NChooseC(n,c-1) * (n-(c-1))/c;
}

// apply a^+_m_d a_m_d operator to a given state (only spin down)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereTwoLandauLevels::AddAd (int index, int m)
{
    int m_index = this->GetIndexFromLzD(m);
    this->FermionToBoson(this->StateDescription[index],this->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
    if ((m_index > this->TemporaryStateLzMax) || (this->TemporaryState[m_index] == 0) )
      {
	return 0.0;
      }
    else 
      {
	return (double)this->TemporaryState[m_index];
      }
}

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereTwoLandauLevels::AduAu (int index, int m)
{
    int m_index = this->GetIndexFromLzU(m);
    this->FermionToBoson(this->StateDescription[index],this->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
    if ((m_index > this->TemporaryStateLzMax) || (this->TemporaryState[m_index] == 0) )
      {
	return 0.0;
      }
    else 
      {
	return (double)this->TemporaryState[m_index];
      }
}

// apply a^+_m_u a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereTwoLandauLevels::AduAu (int index, int m, int n, double& coefficient)
{
    int m_index = this->GetIndexFromLzU(m);
    int n_index = this->GetIndexFromLzU(n);
    this->FermionToBoson(this->StateDescription[index],this->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
    if ((n_index > this->TemporaryStateLzMax) || (this->TemporaryState[n_index] == 0) )
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
    else 
      {
	coefficient = this->TemporaryState[n_index];
	this->TemporaryState[n_index]--;
	if ( m_index > this->TemporaryStateLzMax ) 
	{
	  for ( int index = this->TemporaryStateLzMax + 1 ; index <= m_index ; index++ ) {
	      this->TemporaryState[index] = 0;
	  }
	  this->TemporaryStateLzMax = m_index;
	}
	this->TemporaryState[m_index]++;
	coefficient *= this->TemporaryState[m_index];
	coefficient = sqrt(coefficient);
	return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax));
      }
}

// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereTwoLandauLevels::AddAd (int index, int m, int n, double& coefficient)
{
    int m_index = this->GetIndexFromLzD(m);
    int n_index = this->GetIndexFromLzD(n);
    this->FermionToBoson(this->StateDescription[index],this->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
    if ((n_index > this->TemporaryStateLzMax) || (this->TemporaryState[n_index] == 0) )
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
    else 
      {
	coefficient = this->TemporaryState[n_index];
	this->TemporaryState[n_index]--;
	if ( m_index > this->TemporaryStateLzMax ) 
	{
	  for (int index = this->TemporaryStateLzMax + 1 ; index <= m_index ; index++ ) {
	      this->TemporaryState[index] = 0;
	  }
	  this->TemporaryStateLzMax = m_index;
	}
	this->TemporaryState[m_index]++;
	coefficient *= this->TemporaryState[m_index];
	coefficient = sqrt(coefficient);
	return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax));
      }
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereTwoLandauLevels::AduAd (int index, int m, int n, double& coefficient)
{
    int m_index = this->GetIndexFromLzU(m);
    int n_index = this->GetIndexFromLzD(n);
    this->FermionToBoson(this->StateDescription[index],this->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
    if ((n_index > this->TemporaryStateLzMax) || (this->TemporaryState[n_index] == 0) )
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
    else 
      {
	coefficient = this->TemporaryState[n_index];
	this->TemporaryState[n_index]--;
	if ( m_index > this->TemporaryStateLzMax ) 
	{
	  for (int index = this->TemporaryStateLzMax + 1 ; index <= m_index ; index++ ) {
	      this->TemporaryState[index] = 0;
	  }
	  this->TemporaryStateLzMax = m_index;
	}
	this->TemporaryState[m_index]++;
	coefficient *= this->TemporaryState[m_index];
	coefficient = sqrt(coefficient);
	return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax));
      }
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereTwoLandauLevels::AduAd (int index, int m, double& coefficient)
{
    int m_index = this->GetIndexFromLzU(m);
    int n_index = this->GetIndexFromLzD(m);
    this->FermionToBoson(this->StateDescription[index],this->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
    if ((n_index > this->TemporaryStateLzMax) || (this->TemporaryState[n_index] == 0) )
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
    else 
      {
	coefficient = this->TemporaryState[n_index];
	this->TemporaryState[n_index]--;
	if ( m_index > this->TemporaryStateLzMax ) 
	{
	  for (int index = this->TemporaryStateLzMax + 1 ; index <= m_index ; index++ ) {
	      this->TemporaryState[index] = 0;
	  }
	  this->TemporaryStateLzMax = m_index;
	}
	this->TemporaryState[m_index]++;
	coefficient *= this->TemporaryState[m_index];
	coefficient = sqrt(coefficient);
	return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax));
      }
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereTwoLandauLevels::AddAu (int index, int m, int n, double& coefficient)
{
    int m_index = this->GetIndexFromLzD(m);
    int n_index = this->GetIndexFromLzU(n);
    this->FermionToBoson(this->StateDescription[index],this->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
    if ((n_index > this->TemporaryStateLzMax) || (this->TemporaryState[n_index] == 0) )
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
    else 
      {
	coefficient = this->TemporaryState[n_index];
	this->TemporaryState[n_index]--;
	if ( m_index > this->TemporaryStateLzMax ) 
	{
	  for (int index = this->TemporaryStateLzMax + 1 ; index <= m_index ; index++ ) {
	      this->TemporaryState[index] = 0;
	  }
	  this->TemporaryStateLzMax = m_index;
	}
	this->TemporaryState[m_index]++;
	coefficient *= this->TemporaryState[m_index];
	coefficient = sqrt(coefficient);
	return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax));
      }
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereTwoLandauLevels::AddAu (int index, int m, double& coefficient)
{
    int m_index = this->GetIndexFromLzD(m);
    int n_index = this->GetIndexFromLzU(m);
    this->FermionToBoson(this->StateDescription[index],this->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
    if ((n_index > this->TemporaryStateLzMax) || (this->TemporaryState[n_index] == 0) )
      {
	coefficient = 0.0;
	return this->HilbertSpaceDimension;
      }
    else 
      {
	coefficient = this->TemporaryState[n_index];
	this->TemporaryState[n_index]--;
	if ( m_index > this->TemporaryStateLzMax ) 
	{
	  for (int index = this->TemporaryStateLzMax + 1 ; index <= m_index ; index++ ) {
	      this->TemporaryState[index] = 0;
	  }
	this->TemporaryStateLzMax = m_index;
      }
      this->TemporaryState[m_index]++;
      coefficient *= this->TemporaryState[m_index];
      coefficient = sqrt(coefficient);
      return this->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax)); 
    }
}

// project out any configurations that have particles on levels other than lll
//
// inputVector = vector to apply the projection to
// outputVector = projected vector
// finalSpace = reference to output vector space
  
void  BosonOnSphereTwoLandauLevels::ProjectionInTheLowestLevel(RealVector &inputVector, RealVector & outputVector, BosonOnSphereShort *finalSpace)
{
  unsigned long Etat, Tmp;
  int Idx; 
  for(int i = 0 ; i < finalSpace->GetHilbertSpaceDimension() ; i++)
    {
      //need to translate the format on a single LL into the 2LL picture. 
#ifdef  __64_BITS__  
      int s = 63; // extra shift needed at end
#else      
      int s = 31; // extra shift needed at end
#endif      
      Etat = finalSpace->FermionBasis->StateDescription[i] << (s - (this->LzMaxDown + this->NbrBosons) -(this->LzMaxUp) ); 
      Idx = this->FindStateIndex(Etat);
      if ( Idx < this->HilbertSpaceDimension ) 
	{
	  outputVector[i] = inputVector[Idx];
	}
    }
}


// compute the number of particles in each Landau level
//
// state = ID of the state to handle
// LLOccupationConfiguration = array where the decomposition will be store

void  BosonOnSphereTwoLandauLevels::LandauLevelOccupationNumber(int state, int* lLOccupationConfiguration)
{
  unsigned long State = this->StateDescription[state];
#ifdef  __64_BITS__  
  int idx = 63;  
  while ( idx >= ( 63 - lLOccupationConfiguration[1] - this->LzMaxUp ) ) 
#else    
  int idx = 31;  
  while ( idx >= ( 31 - lLOccupationConfiguration[1] - this->LzMaxUp ) ) 
#endif
    {
      if ( (State >> idx ) & 0x1ul ) 
	{
	  lLOccupationConfiguration[1]++;
	}
      idx--;
    }
  lLOccupationConfiguration[0] = this->NbrBosons - lLOccupationConfiguration[1];
}