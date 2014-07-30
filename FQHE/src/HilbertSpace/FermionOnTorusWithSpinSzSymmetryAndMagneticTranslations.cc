////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of fermion with spin on a torus taking               //
//       into account magnetic translations and the Sz<->-Sz symmetry         //
//                                                                            //
//                        last modification : 30/07/2014                      //
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
#include "HilbertSpace/FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorusWithSpin.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"
#include "QuantumNumber/VectorQuantumNumber.h"
#include "MathTools/FactorialCoefficient.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"
#include "Matrix/Matrix.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "GeneralTools/ArrayTools.h"

#include <cmath>
#include <cstdlib>


using std::cout;
using std::endl;
using std::dec;
using std::hex;
using std::abs;


// default constructor
// 

FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations ()
{
  this->SzParitySign = 1.0;
}  


// basic constructor
// 
// nbrFermions= number of fermions
// totalSpin = twice the total spin value
// maxMomentum = momentum maximum value for a fermion
// xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
// yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity

FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations (int nbrFermions, int totalSpin, int maxMomentum, 
														  int xMomentum, int yMomentum, bool minusSzParity)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (this->NbrFermions+this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions-this->TotalSpin)/2;  
  this->MaxMomentum = maxMomentum;  

  this->NbrMomentum = this->MaxMomentum + 1;
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = FindGCD(this->NbrFermions, this->MaxMomentum);
  cout << "MomentumModulo=" << MomentumModulo<<endl;
  this->XMomentum = xMomentum % this->MomentumModulo;
  this->YMomentum = yMomentum % this->MaxMomentum;
  this->SzParitySign = 1.0;
  if (minusSzParity == true)
    this->SzParitySign = -1.0;

  this->StateShift = 2 * (this->MaxMomentum / this->MomentumModulo);
  this->MomentumIncrement = (this->NbrFermions * this->StateShift / 2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = 0x1ul;
  for (int i = 1; i < this->StateShift; ++i)
    {
      this->MomentumMask <<= 1;
      this->MomentumMask |= 0x1ul;
    }

  this->MaximumSignLookUp = 16;
  this->GenerateSignLookUpTable();
  this->HilbertSpaceDimension = ShiftedEvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum - 1, 0, this->NbrFermionsUp);
  cout << "intermediate Hilbert space dimension = " << HilbertSpaceDimension << endl;
  this->HilbertSpaceDimension = this->GenerateStates();
  cout << "Hilbert space dimension = "<< HilbertSpaceDimension << endl;
  this->Flag.Initialize();

  if (this->HilbertSpaceDimension !=0)
    this->GenerateLookUpTable(1000000);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrMomentum * sizeof(int);
  //  UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
  UsedMemory +=  (1 << MaximumSignLookUp) * sizeof(double);
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

// basic constructor without constraint on Sz
// 
// nbrFermions= number of fermions
// maxMomentum = momentum maximum value for a fermion
// xMomentum = momentum in the x direction (modulo GCD of nbrFermions and maxMomentum)
// yMomentum = momentum in the y direction (modulo GCD of nbrFermions and maxMomentum)  
// minusSzParity = select the  Sz <-> -Sz symmetric sector with negative parity

FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations (int nbrFermions, int maxMomentum, int xMomentum, int yMomentum,
														  bool minusSzParity)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalSpin = 0;
  this->NbrFermionsUp = this->NbrFermions;
  this->NbrFermionsDown = 0;  
  this->MaxMomentum = maxMomentum;  

  this->NbrMomentum = this->MaxMomentum + 1;
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = FindGCD(this->NbrFermions, this->MaxMomentum);
  cout << "MomentumModulo=" << MomentumModulo<<endl;
  this->XMomentum = xMomentum % this->MomentumModulo;
  this->YMomentum = yMomentum % this->MaxMomentum;
  this->SzParitySign = 1.0;
  if (minusSzParity == true)
    this->SzParitySign = -1.0;

  this->StateShift = 2*(this->MaxMomentum / this->MomentumModulo);
  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2*this->MaxMomentum - this->StateShift;
  this->MomentumMask = 0x1ul;
  for (int i = 1; i < this->StateShift; ++i)
    {
      this->MomentumMask <<= 1;
      this->MomentumMask |= 0x1ul;
    }

  this->MaximumSignLookUp = 16;
  this->GenerateSignLookUpTable();
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->MaxMomentum - 1, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  cout << "intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  this->HilbertSpaceDimension = this->GenerateStates(true);
  cout << "Hilbert space dimension = "<< HilbertSpaceDimension << endl;
  this->Flag.Initialize();
  if (this->HilbertSpaceDimension != 0)
    this->GenerateLookUpTable(1000000);
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
#ifdef __DEBUG__
  int UsedMemory = 0;
  UsedMemory += 2 * this->HilbertSpaceDimension * sizeof(int);
  UsedMemory += this->NbrMomentum * sizeof(int);
  //  UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
  UsedMemory +=  (1 << MaximumSignLookUp) * sizeof(double);
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

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations(const FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations& fermions)
{
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;
  this->TotalSpin = fermions.TotalSpin;
  this->IncNbrFermions = fermions.IncNbrFermions;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->NbrFermionStates = fermions.NbrFermionStates;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;
  this->SzParitySign = fermions.SzParitySign;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;
  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::~FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations ()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
      delete[] this->NbrStateInOrbit;      
      delete[] this->ReorderingSign;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;
      delete[] this->SignLookUpTableMask;      
  
      if (this->HilbertSpaceDimension != 0)
	{
	  delete[] this->LookUpTableShift;
	  for (int i = 0; i < this->NbrFermionStates; ++i)
	    delete[] this->LookUpTable[i];
	  delete[] this->LookUpTable;	  
	  
	  for (int i = 1; i <= this->MaxMomentum ; ++i)
	    delete[] this->RescalingFactors[i];
	  delete[] this->RescalingFactors;	  
	}
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations& FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::operator = (const FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;

      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalSpin = fermions.TotalSpin;

  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->NbrFermionStates = fermions.NbrFermionStates;
  this->MomentumModulo = fermions.MomentumModulo;
  this->XMomentum = fermions.XMomentum;
  this->YMomentum = fermions.YMomentum;
  this->SzParitySign = fermions.SzParitySign;

  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::Clone()
{
  return new FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations(*this);
}


// apply a^+_(d,m1) a^+_(d,m2) a_(d,n1) a_(d,n2) operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AddAddAdAd (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  n1<<=1;
  n2<<=1;  
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((0x1ul) << n1)) == 0) || 
      ((State & ((0x1ul) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(0x1ul << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(0x1ul << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  m1<<=1;
  m2<<=1;
  if ((TmpState & (0x1ul << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewHighestBit)
    {
      NewHighestBit = m2;
    }
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
  if ((TmpState & (0x1ul << m1))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewHighestBit)
    {
      NewHighestBit = m1;
    }
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
  TmpState = this->FindCanonicalForm(TmpState, NewHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, NewHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
  coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
  nbrTranslation *= this->StateShift/2;
  return TmpIndex;
}

// apply a^+_(u,m1) a^+_(u,m2) a_(u,n1) a_(u,n2) operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AduAduAuAu (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  n1<<=1;
  n1+=1;
  n2<<=1;
  n2+=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((0x1ul) << n1)) == 0) || 
      ((State & ((0x1ul) << n2)) == 0) || (n1 == n2) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(0x1ul << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(0x1ul << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;  
  m1<<=1;
  m1+=1;
  m2<<=1;
  m2+=1;
  if ((TmpState & (0x1ul << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewHighestBit)
    {
      NewHighestBit = m2;
    }
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
  if ((TmpState & (0x1ul << m1))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewHighestBit)
    {
      NewHighestBit = m1;
    }
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
  TmpState = this->FindCanonicalForm(TmpState, NewHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, NewHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
  coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
  nbrTranslation *= this->StateShift/2;
  return TmpIndex;
}

// apply a^+_(u,m1) a^+_(d,m2) a_(d,n1) a_(u,n2) operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to be applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 
int FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AddAduAdAu (int index, int m1, int m2, int n1, int n2, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  n1<<=1;
  n2<<=1;
  n2+=1;
  if ((n1 > StateHighestBit) || (n2 > StateHighestBit) || ((State & ((0x1ul) << n1)) == 0) || ((State & ((0x1ul) << n2)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewHighestBit = StateHighestBit;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(0x1ul << n2);
  if (NewHighestBit == n2)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(0x1ul << n1);
  if (NewHighestBit == n1)
    while ((NewHighestBit)&&(TmpState >> NewHighestBit) == 0)
      --NewHighestBit;
  m1<<=1;
  m2<<=1;
  m2+=1;
  if ((TmpState & (0x1ul << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m2 > NewHighestBit)
    {
      NewHighestBit = m2;
    }
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
  if ((TmpState & (0x1ul << m1))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m1 > NewHighestBit)
    {
      NewHighestBit = m1;
    }
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
  TmpState = this->FindCanonicalForm(TmpState, NewHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, NewHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, NewHighestBit);
  coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
  coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
  nbrTranslation *= this->StateShift/2;
  return TmpIndex;
}


// apply a^+_m_d a_m_d operator to a given state (only spin down)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m
double FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AddAd (int index, int m)
{
  if (this->StateHighestBit[index] < 2*m)
    return 0.0;
  if ((this->StateDescription[index] & (0x1ul << (m<<1))) == 0)
    return 0.0;
  else
    return 1.0;
}
  

// apply a^+_m_u a_m_u operator to a given state  (only spin up)
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m
double FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AduAu (int index, int m)
{
  if (this->StateHighestBit[index] < 2*m+1)
    return 0.0;
  if ((this->StateDescription[index] & (0x1ul << (2*m+1))) == 0)
    return 0.0;
  else
    return 1.0;
}

// apply a^+_m_u a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AduAu (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m = (m << 1) + 1;
  n = (n << 1) + 1;
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(0x1ul << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if ((State & (0x1ul << m))!= 0)
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
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (0x1ul << m);
  State = this->FindCanonicalForm(State, NewLargestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(State, NewLargestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(State, NewLargestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
      nbrTranslation *= this->StateShift / 2;
    }
  return TmpIndex;
}

// apply a^+_m_d a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AddAd (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m <<= 1;
  n <<= 1;
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(0x1ul << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if ((State & (0x1ul << m))!= 0)
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
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (0x1ul << m);
  State = this->FindCanonicalForm(State, NewLargestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(State, NewLargestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(State, NewLargestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
      nbrTranslation *= this->StateShift / 2;
    }
  return TmpIndex;
}

// apply a^+_m_u a_n_d operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AduAd (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m = (m << 1) + 1;
  n <<= 1;
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0) )
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(0x1ul << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if ((State & (0x1ul << m))!= 0)
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
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (0x1ul << m);
  State = this->FindCanonicalForm(State, NewLargestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(State, NewLargestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(State, NewLargestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
      nbrTranslation *= this->StateShift / 2;
    }
  return TmpIndex;
}

// apply a^+_m_d a_n_u operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation/annihilation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslation = reference on the number of translations to applied to the resulting state to obtain the return orbit describing state
// return value = index of the destination state 

int FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AddAu (int index, int m, int n, double& coefficient, int& nbrTranslation)
{
  int StateHighestBit = this->StateHighestBit[index];
  unsigned long State = this->StateDescription[index];
  m <<= 1;
  n = (n << 1) + 1;  
  if ((n > StateHighestBit) || ((State & (0x1ul << n)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLargestBit = StateHighestBit;
  coefficient = this->SignLookUpTable[(State >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(State >> (n + 16)) & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(State >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(State >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  State &= ~(0x1ul << n);
  if (NewLargestBit == n)
    while ((State >> NewLargestBit) == 0)
      --NewLargestBit;

  if ((State & (0x1ul << m))!= 0)
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
      coefficient *= this->SignLookUpTable[(State >> (m + 16)) & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(State >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(State >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  State |= (0x1ul << m);
  State = this->FindCanonicalForm(State, NewLargestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(State, NewLargestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(State, NewLargestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[index]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
      nbrTranslation *= this->StateShift / 2;
    }
  return TmpIndex;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdu call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin up)
// return value =  multiplicative factor 
double FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AuAu (int index, int n1, int n2)
{
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n1+=1;
  n2<<=1;
  n2+=1;
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((0x1ul) << n1)) == 0) || 
      ((ProdATemporaryState & ((0x1ul) << n2)) == 0) || (n1 == n2))
    return 0.0;  
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  ProdATemporaryState &= ~(0x1ul << n2);
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  ProdATemporaryState &= ~(0x1ul << n1);
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  return coefficient;
}

// apply a_n1_d a_n2_d operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AddAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin down)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AdAd (int index, int n1, int n2)
{  
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n2<<=1;
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((0x1ul) << n1)) == 0) || 
      ((ProdATemporaryState & ((0x1ul) << n2)) == 0) || (n1 == n2))
    return 0.0;
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  ProdATemporaryState &= ~(0x1ul << n2);
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  ProdATemporaryState &= ~(0x1ul << n1);
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  return coefficient;
}

// apply a_n1_u a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 

double FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AuAd (int index, int n1, int n2)
{
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n1+=1;
  n2<<=1;
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((0x1ul) << n1)) == 0) || 
      ((ProdATemporaryState & ((0x1ul) << n2)) == 0) || (n1 == n2))
    return 0.0;
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  ProdATemporaryState &= ~(0x1ul << n2);
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  ProdATemporaryState &= ~(0x1ul << n1);
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  return coefficient;
}

// apply a_n1_d a_n2_u operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be kept in cache until next AduAdd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator (spin up)
// n2 = second index for annihilation operator (spin down)
// return value =  multiplicative factor 
double FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AdAu (int index, int n1, int n2)
{
  this->ProdAHighestBit = this->StateHighestBit[index];
  this->ProdATemporaryState = this->StateDescription[index];
  this->ProdAIndex = index;
  n1<<=1;
  n2<<=1;
  n2+=1;
  if ((n1 > ProdAHighestBit) || (n2 > ProdAHighestBit) || ((ProdATemporaryState & ((0x1ul) << n1)) == 0) || 
      ((ProdATemporaryState & ((0x1ul) << n2)) == 0) || (n1 == n2))
    return 0.0;
  double coefficient = this->SignLookUpTable[(ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  ProdATemporaryState &= ~(0x1ul << n2);
  if (ProdAHighestBit == n2)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  ProdATemporaryState &= ~(0x1ul << n1);
  if (ProdAHighestBit == n1)
    while ((ProdAHighestBit)&&(ProdATemporaryState >> ProdAHighestBit) == 0)
      --ProdAHighestBit;
  return coefficient;
}

// apply a^+_m1_u a^+_m2_u operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin up)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AduAdu (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->ProdATemporaryState;
  int TmpAHighestBit = ProdAHighestBit;
  m1<<=1;
  m1+=1;
  m2<<=1;
  m2+=1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  if (m2 > TmpAHighestBit)
    TmpAHighestBit = m2;
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
  
  if (m1 > TmpAHighestBit)
    TmpAHighestBit = m1;
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
  TmpState = this->FindCanonicalForm(TmpState, TmpAHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, TmpAHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, TmpAHighestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[ProdAIndex]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
      nbrTranslation *= this->StateShift/2;
    }
  return TmpIndex;
}

// apply a^+_m1_d a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin down)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AddAdd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->ProdATemporaryState;
  int TmpAHighestBit = ProdAHighestBit;
  m1<<=1;
  m2<<=1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  if (m2 > TmpAHighestBit)
    TmpAHighestBit = m2;
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
  
  if (m1 > TmpAHighestBit)
    TmpAHighestBit = m1;
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
  TmpState = this->FindCanonicalForm(TmpState, TmpAHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, TmpAHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, TmpAHighestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[ProdAIndex]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
      nbrTranslation *= this->StateShift/2;
    }
  return TmpIndex;
}

  

// apply a^+_m1_u a^+_m2_d operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator (spin up)
// m2 = second index for creation operator (spin down)
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
int FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::AduAdd (int m1, int m2, double& coefficient, int& nbrTranslation)
{
  unsigned long TmpState = this->ProdATemporaryState;
  int TmpAHighestBit = ProdAHighestBit;
  m1<<=1;
  m1+=1;
  m2<<=1;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  if (m2 > TmpAHighestBit)
    TmpAHighestBit = m2;
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
  
  if (m1 > TmpAHighestBit)
    TmpAHighestBit = m1;
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
  TmpState = this->FindCanonicalForm(TmpState, TmpAHighestBit, nbrTranslation);
  if (this->TestXMomentumConstraint(TmpState, TmpAHighestBit) == false)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int TmpIndex = this->FindStateIndex(TmpState, TmpAHighestBit);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= this->RescalingFactors[this->NbrStateInOrbit[ProdAIndex]][this->NbrStateInOrbit[TmpIndex]];
      coefficient *= 1.0 - (2.0 * ((double) ((this->ReorderingSign[TmpIndex] >> nbrTranslation) & 0x1ul)));
      nbrTranslation *= this->StateShift/2;
    }
  return TmpIndex;
}

// convert a state defined in the (Kx,Ky) basis into a state in the Ky basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::ConvertFromKxKyBasis(ComplexVector& state, ParticleOnSphere* space)
{
  FermionOnTorusWithSpin* TmpSpace = (FermionOnTorusWithSpin*) space;
  ComplexVector TmpVector (TmpSpace->HilbertSpaceDimension, true);
  Complex* FourrierCoefficients = new Complex [this->MomentumModulo];
  for (int i = 0; i < this->MomentumModulo; ++i)
    FourrierCoefficients[i] = Phase (-2.0 * M_PI * ((double) (i * this->XMomentum)) / ((double) this->MomentumModulo));
  for (int i = 0; i < TmpSpace->HilbertSpaceDimension; ++i)
    {
      unsigned long TmpState = TmpSpace->StateDescription[i];
      int NbrTranslation = 0;
      int TmpMaxMomentum = TmpSpace->StateHighestBit[i];
      TmpState = this->FindCanonicalFormAndTestXMomentumConstraint(TmpState, TmpMaxMomentum, NbrTranslation);
      if (NbrTranslation >= 0)
	{
	  int Pos = this->FindStateIndex(TmpState, TmpMaxMomentum);
	  if (Pos < this->HilbertSpaceDimension)
	    {
	      TmpVector[i] =  (state[Pos] * (1.0 - (2.0 * ((double) ((this->ReorderingSign[Pos] >> NbrTranslation) & 0x1ul)))) * 
			       FourrierCoefficients[NbrTranslation] / sqrt((double) this->NbrStateInOrbit[Pos]));
	    }
	}
    }
  delete[] FourrierCoefficients;
  return TmpVector;
}


// find canonical form of a state description
//
// stateDescription = unsigned integer describing the state
// stateHighestBit = reference on the highest non-zero bit in the canonical representation
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = canonical form of a state description

unsigned long FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::FindCanonicalForm(unsigned long stateDescription, int& stateHighestBit, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;
  int index = 1;
  stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);  
  while (stateDescriptionReference != stateDescription)
    {
      if (stateDescription < CanonicalState)
	{
	  CanonicalState = stateDescription;
	  nbrTranslation = index;
	}
      stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
      ++index;
    }
  stateHighestBit = 2*this->MaxMomentum-1;
  stateDescription = 0x1ul << stateHighestBit;
  while ((CanonicalState & stateDescription) == 0)
    {
      --stateHighestBit;
      stateDescription >>= 1;
    }
  if (nbrTranslation != 0)
    nbrTranslation = index - nbrTranslation;
  return CanonicalState;
}

// find how many translations on the x direction are needed to obtain the same state
//
// stateDescription = unsigned integer describing the state
// return value = number of translation needed to obtain the same state

int  FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::FindNumberXTranslation(unsigned long stateDescription)
{
  unsigned long TmpState = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
  int index = 1;  
  while (TmpState != stateDescription)
    {
      TmpState = (TmpState >> this->StateShift) | ((TmpState & this->MomentumMask) << this->ComplementaryStateShift);
      ++index;
    }
  return index;
}

// test if a state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescription = unsigned integer describing the state
// maxMomentum = maximum momentum value that can be reached by a fermion in the stateDescription state
// return value = true if the state satisfy the x momentum constraint

bool FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::TestXMomentumConstraint(unsigned long stateDescription, int maxMomentum)
{
  if (this->NbrFermions & 1)
    {
      if (this->XMomentum == 0)
	return true;
      unsigned long TmpState = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
      int index = 1;  
      while (TmpState != stateDescription)
	{
	  TmpState = (TmpState >> this->StateShift) | ((TmpState & this->MomentumMask) << this->ComplementaryStateShift);
	  ++index;
	}
      if (((this->XMomentum * index) % this->MomentumModulo) == 0)
	return true;
      else
	return false;
    }
  else
    {
      unsigned long TmpState2 = stateDescription & this->MomentumMask;
      unsigned long TmpState = (stateDescription >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
      int TmpSignature = 0;
      int index = 1;  
#ifndef  __64_BITS__
      TmpSignature += (this->NbrParticleLookUpTable[TmpState2 & 0xfffful] + this->NbrParticleLookUpTable[(TmpState2 >> 16) & 0xfffful]) & 1;
#else
      TmpSignature += (this->NbrParticleLookUpTable[TmpState2 & 0xfffful] + this->NbrParticleLookUpTable[(TmpState2 >> 16) & 0xfffful]
		       + this->NbrParticleLookUpTable[(TmpState2 >> 32) & 0xfffful] + this->NbrParticleLookUpTable[(TmpState2 >> 48) & 0xfffful]) & 1;
#endif

      while (TmpState != stateDescription)
	{
	  TmpState2 = TmpState & this->MomentumMask;
	  TmpState = (TmpState >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
#ifndef  __64_BITS__
	  TmpSignature += (this->NbrParticleLookUpTable[TmpState2 & 0xfffful] + this->NbrParticleLookUpTable[(TmpState2 >> 16) & 0xfffful]) & 1;
#else
	  TmpSignature += (this->NbrParticleLookUpTable[TmpState2 & 0xfffful] + this->NbrParticleLookUpTable[(TmpState2 >> 16) & 0xfffful]
			   + this->NbrParticleLookUpTable[(TmpState2 >> 32) & 0xfffful] + this->NbrParticleLookUpTable[(TmpState2 >> 48) & 0xfffful]) & 1;
#endif
	  ++index;
	}
      if ((((this->XMomentum * index) - ((this->MomentumModulo * TmpSignature) >> 1)) % this->MomentumModulo) == 0)
	return true;
      else
	return false;
    }
}

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to the x momentum constraint
//
// stateDescription = unsigned integer describing the state
// stateHighestBit = reference on the highest non-zero bit in the canonical representation
// nbrTranslation = number of translation needed to obtain the canonical form
// return value = canonical form of a state description and -1 in nbrTranslation if the state does not fit the x momentum constraint

unsigned long FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::FindCanonicalFormAndTestXMomentumConstraint(unsigned long stateDescription, int& stateHighestBit, int& nbrTranslation)
{
  nbrTranslation = 0;
  unsigned long CanonicalState = stateDescription;
  unsigned long stateDescriptionReference = stateDescription;
  int index = 1;  
  if (this->NbrFermions & 1)
    {
      stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
      while (stateDescriptionReference != stateDescription)
	{
	  if (stateDescription < CanonicalState)
	    {
	      CanonicalState = stateDescription;
	      nbrTranslation = index;
	    }
	  stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
	  ++index;
	}
      if (((this->XMomentum * index) % this->MomentumModulo) != 0)
	{
	  nbrTranslation = -1;
	  return CanonicalState;
	}
    }
  else
    {
      unsigned long TmpState = stateDescription & this->MomentumMask;
      stateDescription = (stateDescription >> this->StateShift) | (TmpState << this->ComplementaryStateShift);
      int TmpSignature = 0;
#ifndef  __64_BITS__
      TmpSignature += (this->NbrParticleLookUpTable[TmpState & 0xfffful] + this->NbrParticleLookUpTable[(TmpState >> 16) & 0xfffful]) & 1;
#else
      TmpSignature += (this->NbrParticleLookUpTable[TmpState & 0xfffful] + this->NbrParticleLookUpTable[(TmpState >> 16) & 0xfffful]
		       + this->NbrParticleLookUpTable[(TmpState >> 32) & 0xfffful] + this->NbrParticleLookUpTable[(TmpState >> 48) & 0xfffful]) & 1;
#endif

      while (stateDescription != stateDescriptionReference)
	{
	  if (stateDescription < CanonicalState)
	    {
	      CanonicalState = stateDescription;
	      nbrTranslation = index;
	    }
	  TmpState = stateDescription & this->MomentumMask;
	  stateDescription = (stateDescription >> this->StateShift) | (TmpState << this->ComplementaryStateShift);
#ifndef  __64_BITS__
	  TmpSignature += (this->NbrParticleLookUpTable[TmpState & 0xfffful] + this->NbrParticleLookUpTable[(TmpState >> 16) & 0xfffful]) & 1;
#else
	  TmpSignature += (this->NbrParticleLookUpTable[TmpState & 0xfffful] + this->NbrParticleLookUpTable[(TmpState >> 16) & 0xfffful]
			   + this->NbrParticleLookUpTable[(TmpState >> 32) & 0xfffful] + this->NbrParticleLookUpTable[(TmpState >> 48) & 0xfffful]) & 1;
#endif
	  ++index;
	}
      if ((((this->XMomentum * index) - ((this->MomentumModulo * TmpSignature) >> 1)) % this->MomentumModulo) != 0)
	{
	  nbrTranslation = -1;
	  return CanonicalState;
	}
    }
  if (nbrTranslation != 0)
    {
      stateHighestBit = 2*this->MaxMomentum-1;
      stateDescription = 0x1ul << stateHighestBit;
      while ((CanonicalState & stateDescription) ==0)      
	{
	  --stateHighestBit;
	  stateDescription >>= 1;
	}
      nbrTranslation = index - nbrTranslation;
    }
  return CanonicalState;
}

// generate all states corresponding to the constraints
// 
// fullSzFlag = if true, does not apply the Sz contraint
// fullKyFlag = if true, does not apply the Ky contraint
// return value = hilbert space dimension

int FermionOnTorusWithSpinSzSymmetryAndMagneticTranslations::GenerateStates(bool fullSzFlag, bool fullKyFlag)
{
  int maxHighestBit = 2*(this->MaxMomentum + 1);
  this->StateDescription = new unsigned long [this->HilbertSpaceDimension];
  this->StateHighestBit = new int [this->HilbertSpaceDimension];  
  if (fullKyFlag == false)
    {
      if (fullSzFlag == false)
	this->HilbertSpaceDimension = this->RawGenerateStates(this->NbrFermions, this->MaxMomentum - 1, 0, this->NbrFermionsUp, 0l);
      else
	this->HilbertSpaceDimension = this->RawGenerateStates(this->NbrFermions, 2 * this->MaxMomentum - 1, 2 * this->MaxMomentum - 1, 0l, 0);    
    }
  else
    {
      this->HilbertSpaceDimension = this->RawGenerateStates(this->NbrFermions, this->MaxMomentum - 1, 0l);
      int TmpMaxHighestBit = 2*(this->MaxMomentum + 1);
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	{
	  while (((this->StateDescription[i] >> TmpMaxHighestBit) & 0x1ul) == 0x0ul)
	    --TmpMaxHighestBit;
	  this->StateHighestBit[i] = TmpMaxHighestBit;
	}
    }
  int* TmpNbrStateDescription = new int [maxHighestBit];  
  for (int i = 0; i < maxHighestBit; ++i)
    {
      TmpNbrStateDescription[i] = 0;
    }
  int NbrTranslation;
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->StateDescription[i] = this->FindCanonicalForm(this->StateDescription[i], this->StateHighestBit[i], NbrTranslation);
      ++TmpNbrStateDescription[this->StateHighestBit[i]];
    }  
  unsigned long** TmpStateDescription = new unsigned long* [maxHighestBit];  
  bool** TmpStateDescriptionFlag = new bool* [maxHighestBit];  
  for (int i = 0; i < maxHighestBit; ++i)
    {
      if (TmpNbrStateDescription[i] != 0)
	{
	  TmpStateDescription[i] = new unsigned long [TmpNbrStateDescription[i]];
	  TmpStateDescriptionFlag[i] = new bool [TmpNbrStateDescription[i]];
	}
      else
	{
	  TmpStateDescription[i] = 0;
	  TmpStateDescriptionFlag[i] = 0;
	}
      TmpNbrStateDescription[i] = 0;
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    TmpStateDescription[this->StateHighestBit[i]][TmpNbrStateDescription[this->StateHighestBit[i]]++] = this->StateDescription[i];      
  int TmpHilbertSpaceDimension = 0;
  for (int i = 0; i < maxHighestBit; ++i) 
    {
      int CurrentNbrState = TmpNbrStateDescription[i];
      if (CurrentNbrState > 0)
	{
	  unsigned long* TmpStateArray = TmpStateDescription[i];
	  bool* TmpStateArrayFlag = TmpStateDescriptionFlag[i];
	  SortArrayUpOrdering(TmpStateArray, CurrentNbrState);
	  int TmpNbr = CurrentNbrState;
	  if (this->TestXMomentumConstraint(TmpStateArray[0], i) == true)
	    {
	      TmpStateArrayFlag[0] = true;
	    }
	  else
	    {
	      TmpStateArrayFlag[0] = false;
	      --TmpNbr;
	    }
	  for (int j = 1; j < CurrentNbrState; ++j)
	    {
	      while ((j < CurrentNbrState) && (TmpStateArray[j - 1] == TmpStateArray[j]))
		{
		  TmpStateArrayFlag[j] = false;
		  --TmpNbr;
		  ++j;
		}
	      if (j < CurrentNbrState)
		{
		  if (this->TestXMomentumConstraint(TmpStateArray[j], i) == true)
		    {
		      TmpStateArrayFlag[j] = true;
		    }
		  else
		    {
		      TmpStateArrayFlag[j] = false;
		      --TmpNbr;	      
		    }
		}
	    }
	  TmpHilbertSpaceDimension += TmpNbr;
	}
    }
  delete[] this->StateHighestBit;
  delete[] this->StateDescription;
  if (TmpHilbertSpaceDimension > 0)
    {
      this->StateDescription = new unsigned long [TmpHilbertSpaceDimension];
      this->StateHighestBit = new int [TmpHilbertSpaceDimension];
      this->NbrStateInOrbit = new int [TmpHilbertSpaceDimension];
      this->ReorderingSign = new unsigned long [TmpHilbertSpaceDimension];
      int Pos = 0;
      unsigned long TmpState;
      unsigned long TmpState2;
      int TmpSignature;
      unsigned long TmpReorderingSignUp;
      unsigned long TmpReorderingSignDown;
      unsigned long TmpReorderingSign;
      int TmpNbrParticleUp;
      int TmpNbrParticleDown;
      int TmpNbrParticle;
      for (int i = 0; i <maxHighestBit; ++i) 
	{
	  int CurrentNbrState = TmpNbrStateDescription[i];
	  if (CurrentNbrState > 0)
	    {
	      unsigned long* TmpStateArray = TmpStateDescription[i];
	      bool* TmpStateArrayFlag = TmpStateDescriptionFlag[i];
	      for (int j = 0; j < CurrentNbrState; ++j)
		{
		  // old code:
		  if (TmpStateArrayFlag[j] == true)
		    {
		      this->StateDescription[Pos] = TmpStateArray[j];
		      this->StateHighestBit[Pos] = i;
		      this->NbrStateInOrbit[Pos] = this->FindNumberXTranslation(this->StateDescription[Pos]);		  
		      TmpState = this->StateDescription[Pos];
		      TmpSignature = 0;
		      TmpReorderingSign = 0x0ul;
		      if ((this->NbrFermions & 1) == 0)
			{
			  for (int k = 1; k <= this->NbrStateInOrbit[Pos]; ++k)
			    {
			      TmpState2 = TmpState & this->MomentumMask;
			      TmpState =  (TmpState >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
			      TmpNbrParticle = 0;
			      for (int l = 0; l < this->StateShift; ++l)
				{
				  if (TmpState2 & 1)
				    ++TmpNbrParticle;
				  TmpState2 >>= 1;
				}
			      if (TmpNbrParticle & 1)
				{			 
				  TmpReorderingSign |= (((TmpReorderingSign << 1) & (0x1ul << k))) ^ (0x1ul << k);
				  ++TmpSignature;
				}
			      else
				{
				  TmpReorderingSign |= ((TmpReorderingSign << 1) & (0x1ul << k));
				}			  
			    }
			}
		      // do not assign this result, just keep for comparison
		      // this->ReorderingSign[Pos] = TmpReorderingSign;
		      // ++Pos;
		    }
		  // end old code
		  // new code:
		  if (TmpStateArrayFlag[j] == true)
		    {
		      this->StateDescription[Pos] = TmpStateArray[j];
		      this->StateHighestBit[Pos] = i;
		      this->NbrStateInOrbit[Pos] = this->FindNumberXTranslation(this->StateDescription[Pos]);		  
		      TmpState = this->StateDescription[Pos];
		      TmpSignature = 0;
		      TmpReorderingSignUp = 0x0ul;
		      TmpReorderingSignDown = 0x0ul;
		      if (((this->NbrFermionsUp & 1) == 0)||((this->NbrFermionsDown & 1) == 0))
			{
			  for (int k = 1; k <= this->NbrStateInOrbit[Pos]; ++k)
			    {
			      TmpState2 = TmpState & this->MomentumMask;
			      TmpState =  (TmpState >> this->StateShift) | (TmpState2 << this->ComplementaryStateShift);
			      TmpNbrParticleUp = 0;
			      TmpNbrParticleDown = 0;
			      for (int l = 0; l < this->StateShift; l+=2)
				{
				  if (TmpState2 & 0x1ul)
				    ++TmpNbrParticleDown;
				  if (TmpState2 & 0x2ul)
				    ++TmpNbrParticleUp;
				  TmpState2 >>= 2;
				}
			      if ((TmpNbrParticleUp & 1)&&((this->NbrFermionsUp & 1) == 0))
				{			 
				  TmpReorderingSignUp |= (((TmpReorderingSignUp << 1) & (0x1ul << k))) ^ (0x1ul << k);
				  ++TmpSignature;
				}
			      else
				{
				  TmpReorderingSignUp |= ((TmpReorderingSignUp << 1) & (0x1ul << k));
				}
			      if ((TmpNbrParticleDown & 1)&&((this->NbrFermionsDown & 1) == 0))
				{			 
				  TmpReorderingSignDown |= (((TmpReorderingSignDown << 1) & (0x1ul << k))) ^ (0x1ul << k);
				  ++TmpSignature;
				}
			      else
				{
				  TmpReorderingSignDown |= ((TmpReorderingSignDown << 1) & (0x1ul << k));
				}			  
			      
			    }
			}
		      this->ReorderingSign[Pos] = TmpReorderingSign;// TmpReorderingSignUp^TmpReorderingSignDown;
		      // cout << "old: "<<TmpReorderingSign<<", new: "<<this->ReorderingSign[Pos]<<endl;
		      ++Pos;
		    }
		  
		}
	      delete[] TmpStateArray;
	      delete[] TmpStateArrayFlag;
	    }
	}
      delete[] TmpStateDescription;
      delete[] TmpStateDescriptionFlag;
      delete[] TmpNbrStateDescription;
    }

  return TmpHilbertSpaceDimension;
}

