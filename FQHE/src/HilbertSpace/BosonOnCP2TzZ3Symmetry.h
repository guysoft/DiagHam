////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Cecile Repellin                 //
//                                                                            //
//                                                                            //
//                          class of bosons on CP2                            //
//                     including Tz <-> -Tz and Z3 symmetry                   //
//                                                                            //
//                        last modification : 24/01/2013                      //
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


#ifndef BOSONONCP2TZZ3SYMMETRY_H
#define BOSONONCP2TZZ3SYMMETRY_H

#include "config.h"
#include "HilbertSpace/BosonOnCP2TzSymmetry.h"

#include <iostream>


#define M_SQRT3   1.73205080756888
#define M_SQRT1_3 0.577350269189626


class BosonOnCP2TzZ3Symmetry : public BosonOnCP2TzSymmetry
{

  public:

  // default constructor
  // 
  BosonOnCP2TzZ3Symmetry ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuanta = number of flux quanta (p)
  // totalJz = total value of jz
  // totalKz = total value of kz
  // minusTzParity = select the  Tz <-> -Tz symmetric sector with negative parity
  // memory = amount of memory granted for precalculations
  BosonOnCP2TzZ3Symmetry (int nbrBosons, int nbrFluxQuanta, int totalTz, int totalY, bool minusTzParity, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  BosonOnCP2TzZ3Symmetry(const BosonOnCP2TzZ3Symmetry& bosons);

  // destructor
  //
  ~BosonOnCP2TzZ3Symmetry ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnCP2TzZ3Symmetry& operator = (const BosonOnCP2TzZ3Symmetry& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // save Hilbert space description to disk
  //
  // fileName = name of the file where the Hilbert space description has to be saved
  // return value = true if no error occured
  virtual bool WriteHilbertSpace (char* fileName);
  
  
 protected:
   
  // get canonical expression of a given state and its symmetry
  //
  // initialState = ID of the state (in the unsymmetrized basis table) that has to be converted to its canonical expression
  // tzSymmetry = reference to an integer that describes the symmetry of the state: 1 if the state is symmetric, 0 if not 
  // return value = corresponding canonical state
  virtual unsigned long GetCanonicalState (int initialState, int& symmetrySignature);
  
  // get canonical expression of a given state and its symmetry
  //
  // tzSymmetry = reference to an integer that describes the symmetry of the state: 1 if the state is symmetric, 0 if not 
  // canonicalFlag = reference on an integer that says if the state was canonical (0) or not (1)
  // return value = corresponding canonical state
  virtual unsigned long GetCanonicalStateFromTemporaryBosonicPartition (int& symmetrySignature, int& tzCanonicalFlag);
  
  // factorized code that is used to symmetrize result of the AdxAdy operations
  //
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state
  virtual int SymmetrizeAdAdResult(double& coefficient);

};

inline unsigned long BosonOnCP2TzZ3Symmetry::GetCanonicalState (int initialState, int& symmetrySignature)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[initialState], this->FermionBasis->StateLzMax[initialState], this->TemporaryState, this->TemporaryStateLzMax);
  unsigned long canonicalState = this->FermionBasis->StateDescription[initialState];
  int TzStateBosonicLzMax = 0;
  int Rot1StateBosonicLzMax = 0;
  int TzRot1StateBosonicLzMax = 0;
  int Rot2StateBosonicLzMax = 0;
  int TzRot2StateBosonicLzMax = 0;
  unsigned long* TzStateBosonic;
  unsigned long* Rot1StateBosonic;
  unsigned long* TzRot1StateBosonic;
  unsigned long* Rot2StateBosonic;
  unsigned long* TzRot2StateBosonic;
  TzStateBosonic = new unsigned long[this->NbrLzValue];
  Rot1StateBosonic = new unsigned long[this->NbrLzValue];
  TzRot1StateBosonic = new unsigned long[this->NbrLzValue];
  Rot2StateBosonic = new unsigned long[this->NbrLzValue];
  TzRot2StateBosonic = new unsigned long[this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
  {
   TzStateBosonic[i] = 0x0ul; 
   Rot1StateBosonic[i] = 0x0ul; 
   TzRot1StateBosonic[i] = 0x0ul; 
   Rot2StateBosonic[i] = 0x0ul; 
   TzRot2StateBosonic[i] = 0x0ul; 
  }
  for (int index = 0; index <= this->TemporaryStateLzMax; ++index)
  {
   if (this->TemporaryState[index] > 0)
    {
      int indexTzState = this->GetLinearizedIndex(-this->quantumNumberTz[index], this->quantumNumberY[index], 1);
      int indexRot1State = this->GetLinearizedIndex(this->NbrFluxQuanta - 2*this->quantumNumberR[index] - this->quantumNumberS[index], this->NbrFluxQuanta - 3*this->quantumNumberS[index], 1);
      int indexTzRot1State = this->GetLinearizedIndex(-(this->NbrFluxQuanta - 2*this->quantumNumberR[index] - this->quantumNumberS[index]), this->NbrFluxQuanta - 3*this->quantumNumberS[index], 1);
      int indexRot2State = this->GetLinearizedIndex(2*this->quantumNumberS[index] + this->quantumNumberR[index] - this->NbrFluxQuanta, this->NbrFluxQuanta - 3*this->quantumNumberR[index], 1);
      int indexTzRot2State = this->GetLinearizedIndex(-(2*this->quantumNumberS[index] + this->quantumNumberR[index] - this->NbrFluxQuanta), this->NbrFluxQuanta - 3*this->quantumNumberR[index], 1);
      
      TzStateBosonic[indexTzState] = this->TemporaryState[index];
      Rot1StateBosonic[indexRot1State] = this->TemporaryState[index];
      TzRot1StateBosonic[indexTzRot1State] = this->TemporaryState[index];
      Rot2StateBosonic[indexRot2State] = this->TemporaryState[index];
      TzRot2StateBosonic[indexTzRot2State] = this->TemporaryState[index];
      
      if ( indexTzState > TzStateBosonicLzMax)
	TzStateBosonicLzMax = indexTzState;
      if ( indexRot1State > Rot1StateBosonicLzMax)
	Rot1StateBosonicLzMax = indexRot1State;
      if ( indexTzRot1State > TzRot1StateBosonicLzMax)
	TzRot1StateBosonicLzMax = indexTzRot1State;
      if ( indexRot2State > Rot2StateBosonicLzMax)
	Rot2StateBosonicLzMax = indexRot2State;
      if ( indexTzRot2State > TzRot2StateBosonicLzMax)
	TzRot2StateBosonicLzMax = indexTzRot2State;
    }
  }
  
  if (this->BosonToFermion(TzStateBosonic, TzStateBosonicLzMax) > canonicalState)
    canonicalState = this->BosonToFermion(TzStateBosonic, TzStateBosonicLzMax);
  if (this->BosonToFermion(Rot1StateBosonic, Rot1StateBosonicLzMax) > canonicalState)
    canonicalState = this->BosonToFermion(Rot1StateBosonic, Rot1StateBosonicLzMax);
  if (this->BosonToFermion(TzRot1StateBosonic, TzRot1StateBosonicLzMax) > canonicalState)
    canonicalState = this->BosonToFermion(TzRot1StateBosonic, TzRot1StateBosonicLzMax);
  if (this->BosonToFermion(Rot2StateBosonic, Rot2StateBosonicLzMax) > canonicalState)
    canonicalState = this->BosonToFermion(Rot2StateBosonic, Rot2StateBosonicLzMax);
  if (this->BosonToFermion(TzRot2StateBosonic, TzRot2StateBosonicLzMax) > canonicalState)
    canonicalState = this->BosonToFermion(TzRot2StateBosonic, TzRot2StateBosonicLzMax);
    
  symmetrySignature = 0;
  if (this->BosonToFermion(Rot1StateBosonic, Rot1StateBosonicLzMax) == this->FermionBasis->StateDescription[initialState])
    symmetrySignature = 2;
  if (this->BosonToFermion(TzStateBosonic, TzStateBosonicLzMax) == this->FermionBasis->StateDescription[initialState])
  {
    symmetrySignature += 1;
    return canonicalState;
  }
  if (this->BosonToFermion(TzRot1StateBosonic, TzRot1StateBosonicLzMax) == this->BosonToFermion(Rot1StateBosonic, Rot1StateBosonicLzMax))
  {
    symmetrySignature += 1;
    return canonicalState;
  }
  if (this->BosonToFermion(TzRot2StateBosonic, TzRot2StateBosonicLzMax) == this->BosonToFermion(Rot2StateBosonic, Rot2StateBosonicLzMax))
  {
    symmetrySignature += 1;
    return canonicalState;
  }
//   cout << "init " << symmetrySignature << " " << initialState << endl;
  return canonicalState;
}

inline unsigned long BosonOnCP2TzZ3Symmetry::GetCanonicalStateFromTemporaryBosonicPartition (int& symmetrySignature, int& tzCanonicalFlag)
{
  unsigned long initialState = this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax);
  unsigned long canonicalState = this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax);
  int TzStateBosonicLzMax = 0;
  int Rot1StateBosonicLzMax = 0;
  int TzRot1StateBosonicLzMax = 0;
  int Rot2StateBosonicLzMax = 0;
  int TzRot2StateBosonicLzMax = 0;
  unsigned long* TzStateBosonic;
  unsigned long* Rot1StateBosonic;
  unsigned long* TzRot1StateBosonic;
  unsigned long* Rot2StateBosonic;
  unsigned long* TzRot2StateBosonic;
  TzStateBosonic = new unsigned long[this->NbrLzValue];
  Rot1StateBosonic = new unsigned long[this->NbrLzValue];
  TzRot1StateBosonic = new unsigned long[this->NbrLzValue];
  Rot2StateBosonic = new unsigned long[this->NbrLzValue];
  TzRot2StateBosonic = new unsigned long[this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
  {
   TzStateBosonic[i] = 0x0ul; 
   Rot1StateBosonic[i] = 0x0ul; 
   TzRot1StateBosonic[i] = 0x0ul; 
   Rot2StateBosonic[i] = 0x0ul; 
   TzRot2StateBosonic[i] = 0x0ul; 
  }
  for (int index = 0; index <= this->TemporaryStateLzMax; ++index)
  {
   if (this->TemporaryState[index] > 0)
    {
      int indexTzState = this->GetLinearizedIndex(-this->quantumNumberTz[index], this->quantumNumberY[index], 1);
      int indexRot1State = this->GetLinearizedIndex(this->NbrFluxQuanta - 2*this->quantumNumberR[index] - this->quantumNumberS[index], this->NbrFluxQuanta - 3*this->quantumNumberS[index], 1);
      int indexTzRot1State = this->GetLinearizedIndex(-(this->NbrFluxQuanta - 2*this->quantumNumberR[index] - this->quantumNumberS[index]), this->NbrFluxQuanta - 3*this->quantumNumberS[index], 1);
      int indexRot2State = this->GetLinearizedIndex(2*this->quantumNumberS[index] + this->quantumNumberR[index] - this->NbrFluxQuanta, this->NbrFluxQuanta - 3*this->quantumNumberR[index], 1);
      int indexTzRot2State = this->GetLinearizedIndex(-(2*this->quantumNumberS[index] + this->quantumNumberR[index] - this->NbrFluxQuanta), this->NbrFluxQuanta - 3*this->quantumNumberR[index], 1);
      
//       cout << index << " " << indexTzState << " " << indexRot1State << " " << indexTzRot1State << " " << indexRot2State << " " << indexTzRot2State << endl;
      
      TzStateBosonic[indexTzState] = this->TemporaryState[index];
      Rot1StateBosonic[indexRot1State] = this->TemporaryState[index];
      TzRot1StateBosonic[indexTzRot1State] = this->TemporaryState[index];
      Rot2StateBosonic[indexRot2State] = this->TemporaryState[index];
      TzRot2StateBosonic[indexTzRot2State] = this->TemporaryState[index];
      
      if ( indexTzState > TzStateBosonicLzMax)
	TzStateBosonicLzMax = indexTzState;
      if ( indexRot1State > Rot1StateBosonicLzMax)
	Rot1StateBosonicLzMax = indexRot1State;
      if ( indexTzRot1State > TzRot1StateBosonicLzMax)
	TzRot1StateBosonicLzMax = indexTzRot1State;
      if ( indexRot2State > Rot2StateBosonicLzMax)
	Rot2StateBosonicLzMax = indexRot2State;
      if ( indexTzRot2State > TzRot2StateBosonicLzMax)
	TzRot2StateBosonicLzMax = indexTzRot2State;
    }
  }
  tzCanonicalFlag = 0;
  symmetrySignature = 0;
//   int TmpTzSymmetry = 0;
  
  if (this->BosonToFermion(TzStateBosonic, TzStateBosonicLzMax) > canonicalState)
  {
    canonicalState = this->BosonToFermion(TzStateBosonic, TzStateBosonicLzMax);
    tzCanonicalFlag = 1;
  }
  if (this->BosonToFermion(Rot1StateBosonic, Rot1StateBosonicLzMax) > canonicalState)
  {
    canonicalState = this->BosonToFermion(Rot1StateBosonic, Rot1StateBosonicLzMax);
    tzCanonicalFlag = 0;
  }
  if (this->BosonToFermion(TzRot1StateBosonic, TzRot1StateBosonicLzMax) > canonicalState)
  {
    canonicalState = this->BosonToFermion(TzRot1StateBosonic, TzRot1StateBosonicLzMax);
    tzCanonicalFlag = 1;
  }
  if (this->BosonToFermion(Rot2StateBosonic, Rot2StateBosonicLzMax) > canonicalState)
  {
    canonicalState = this->BosonToFermion(Rot2StateBosonic, Rot2StateBosonicLzMax);
    tzCanonicalFlag = 0;
  }
  if (this->BosonToFermion(TzRot2StateBosonic, TzRot2StateBosonicLzMax) > canonicalState)
  {
    canonicalState = this->BosonToFermion(TzRot2StateBosonic, TzRot2StateBosonicLzMax);
    tzCanonicalFlag = 1;
  }
 //   cout << symmetrySignature << endl;
 if (this->BosonToFermion(Rot1StateBosonic, Rot1StateBosonicLzMax) == initialState)
   symmetrySignature = 2;
 if (this->BosonToFermion(TzStateBosonic, TzStateBosonicLzMax) == initialState)
 {
   symmetrySignature += 1;
   return canonicalState;
 }
 if (this->BosonToFermion(TzRot1StateBosonic, TzRot1StateBosonicLzMax) == this->BosonToFermion(Rot1StateBosonic, Rot1StateBosonicLzMax))
 {
   symmetrySignature += 1;
   return canonicalState;
 }
 if (this->BosonToFermion(TzRot2StateBosonic, TzRot2StateBosonicLzMax) == this->BosonToFermion(Rot2StateBosonic, Rot2StateBosonicLzMax))
 {
   symmetrySignature += 1;
   return canonicalState;
 }
 return canonicalState;
}


// factorized code that is used to symmetrize result of the AdxAdy operations
//
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state

inline int BosonOnCP2TzZ3Symmetry::SymmetrizeAdAdResult(double& coefficient)
{
  unsigned long TmpState;
  int symmetrySignature;
  int tzCanonicalFlag;
  TmpState = this->GetCanonicalStateFromTemporaryBosonicPartition(symmetrySignature, tzCanonicalFlag);
  int TmpStateLzMax = this->NbrLzValue + this->NbrBosons - 1;
  while (((TmpState >> TmpStateLzMax) & 1) == 0)
    --TmpStateLzMax;
  if ((symmetrySignature & 2) != 0 && ((symmetrySignature & 1) != 0))
    {
      if (this->TzParitySign < 0.0)
	return this->HilbertSpaceDimension;
      if ((this->ProdASignature & 1) == 0)
	coefficient *= M_SQRT2;
      if ((this->ProdASignature & 2) == 0)
	coefficient *= M_SQRT3;
//       cout << symmetrySignature << " " << this->FindStateIndex(TmpState, TmpStateLzMax) << endl;
      return this->FindStateIndex(TmpState, TmpStateLzMax);
    }
  if ((symmetrySignature & 2) != 0)
    {
      if (this->TzParitySign < 0.0)
	return this->HilbertSpaceDimension;
      if ((this->ProdASignature & 1) != 0)
	coefficient *= M_SQRT1_2;
      if ((this->ProdASignature & 2) == 0)
	coefficient *= M_SQRT3;
//       cout << symmetrySignature << " " << this->FindStateIndex(TmpState, TmpStateLzMax) << endl;
      return this->FindStateIndex(TmpState, TmpStateLzMax);
    }
  
  if ((symmetrySignature & 1) != 0)
    {
      if (this->TzParitySign < 0.0)
	return this->HilbertSpaceDimension;
      if ((this->ProdASignature & 1) == 0)
	coefficient *= M_SQRT2;
      if ((this->ProdASignature & 2) != 0)
	coefficient *= M_SQRT1_3;
//       cout << symmetrySignature << " " << this->FindStateIndex(TmpState, TmpStateLzMax) << endl;
      return this->FindStateIndex(TmpState, TmpStateLzMax);
    }
  if (tzCanonicalFlag == 1)
    coefficient *= this->TzParitySign;
  if ((this->ProdASignature & 1) != 0)
    coefficient *= M_SQRT1_2;
  if ((this->ProdASignature & 2) != 0)
    coefficient *= M_SQRT1_3;
//   cout << symmetrySignature << " " << this->FindStateIndex(TmpState, TmpStateLzMax) << endl;
  return this->FindStateIndex(TmpState, TmpStateLzMax);
}
#endif


