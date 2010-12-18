////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of bosons on sphere for system size such that            //
//         LzMax + NbrBosons - 1 < 63 or 31 (64 bits or 32bits systems)       //
//                                                                            //
//                        last modification : 10/10/2007                      //
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
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h" 
#include "GeneralTools/StringTools.h"
#include "GeneralTools/ArrayTools.h"

#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"

#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


// default constructor
//

BosonOnSphereShort::BosonOnSphereShort ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a boson

BosonOnSphereShort::BosonOnSphereShort (int nbrBosons, int totalLz, int lzMax)
{
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  this->NbrLzValue = this->LzMax + 1;
  if (nbrBosons > 0)
    this->FermionBasis = new FermionOnSphere(nbrBosons, totalLz, lzMax + nbrBosons - 1);
  else
    this->FermionBasis = new FermionOnSphere(nbrBosons, totalLz, lzMax);
  this->HilbertSpaceDimension = this->FermionBasis->HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = this->FermionBasis->LargeHilbertSpaceDimension;

  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->Flag.Initialize();
  int TmpLzMax = this->LzMax;
  this->TargetSpace = this;

  if (this->ShiftedTotalLz < TmpLzMax)
    {
      TmpLzMax = this->ShiftedTotalLz;	  
    }
  this->KeptCoordinates = new int;
  (*(this->KeptCoordinates)) = -1;
  this->Minors = 0;

}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSphereShort::BosonOnSphereShort(const BosonOnSphereShort& bosons)
{
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  this->Minors = bosons.Minors;
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];  
}

// destructor
//

BosonOnSphereShort::~BosonOnSphereShort ()
{
  if (this->FermionBasis != 0)
    delete this->FermionBasis;
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
    }
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSphereShort& BosonOnSphereShort::operator = (const BosonOnSphereShort& bosons)
{
  delete[] this->TemporaryState;
  delete[] this->ProdATemporaryState;
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      if (this->Minors != 0)
	{
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    if (this->Minors[i] != 0)
	      delete[] this->Minors[i];
	  delete[] this->Minors;
	}
      delete this->KeptCoordinates;
    }
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->ShiftedTotalLz = bosons. ShiftedTotalLz;
  this->LzMax = bosons.LzMax;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  this->KeptCoordinates = bosons.KeptCoordinates;
  this->Minors = bosons.Minors;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];

  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSphereShort::Clone()
{
  return new BosonOnSphereShort(*this);
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void BosonOnSphereShort::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  this->TargetSpace = (BosonOnSphereShort*) targetSpace;
  this->FermionBasis->SetTargetSpace(this->TargetSpace->FermionBasis);
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int BosonOnSphereShort::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> BosonOnSphereShort::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* BosonOnSphereShort::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* BosonOnSphereShort::ExtractSubspace (AbstractQuantumNumber& q, 
						      SubspaceSpaceConverter& converter)
{
  return 0;
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShort::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  if ((n1 > this->TemporaryStateLzMax) || (n2 > this->TemporaryStateLzMax) || (this->TemporaryState[n1] == 0) || (this->TemporaryState[n2] == 0) || ((n1 == n2) && (this->TemporaryState[n1] == 1)))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  for (int i = this->TemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = 0ul;
  coefficient = this->TemporaryState[n2];
  --this->TemporaryState[n2];
  coefficient *= this->TemporaryState[n1];
  --this->TemporaryState[n1];
  ++this->TemporaryState[m2];
  coefficient *= this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->TargetSpace->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShort::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  --nbrIndices;
  for (int i = 0; i <= nbrIndices; ++i)
    {
      if ((n[i] > this->TemporaryStateLzMax) || (this->TemporaryState[n[i]] == 0))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
    }
  for (int i = this->TemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = 0ul;
  coefficient = 1.0;
  for (int i = nbrIndices; i >= 0; --i)
    {
      if (this->TemporaryState[n[i]] == 0)
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension; 	    
	}
      coefficient *= (double) this->TemporaryState[n[i]];
      --this->TemporaryState[n[i]];
    }
  for (int i = nbrIndices; i >= 0; --i)
    {
      ++this->TemporaryState[m[i]];
      coefficient *= (double) this->TemporaryState[m[i]];
    }
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->TargetSpace->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double BosonOnSphereShort::AA (int index, int n1, int n2)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
  if ((n1 > this->ProdATemporaryStateLzMax) || (n2 > this->ProdATemporaryStateLzMax) || 
      (this->ProdATemporaryState[n1] == 0) || (this->ProdATemporaryState[n2] == 0) || ((n1 == n2) && (this->ProdATemporaryState[n1] == 1)))
    {
      return 0.0;
    }
  double Coefficient = this->ProdATemporaryState[n2];
  --this->ProdATemporaryState[n2];
  Coefficient *= this->ProdATemporaryState[n1];
  --this->ProdATemporaryState[n1];
  for (int i = this->ProdATemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0ul;
  return sqrt(Coefficient);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double BosonOnSphereShort::ProdA (int index, int* n, int nbrIndices)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->ProdATemporaryState, this->ProdATemporaryStateLzMax);
  int TmpCoefficient = 1;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      if (n[i] > this->ProdATemporaryStateLzMax)
	return 0.0;
      unsigned long& Tmp = this->ProdATemporaryState[n[i]];
      if (Tmp == 0)
	return 0.0;
      TmpCoefficient *= Tmp;
      --Tmp;
    }
  for (int i = this->ProdATemporaryStateLzMax + 1; i < this->NbrLzValue; ++i)
    this->ProdATemporaryState[i] = 0;
  return sqrt((double) TmpCoefficient);
}

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShort::AdAd (int m1, int m2, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  ++this->TemporaryState[m2];
  coefficient = this->TemporaryState[m2];
  ++this->TemporaryState[m1];
  coefficient *= this->TemporaryState[m1];
  coefficient = sqrt(coefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->TargetSpace->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShort::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->TemporaryState[i] = this->ProdATemporaryState[i];
  int TmpCoefficient = 1;
  for (int i = 0; i < nbrIndices; ++i)
    TmpCoefficient *= ++this->TemporaryState[m[i]];
  coefficient = sqrt((double) TmpCoefficient);
  this->TemporaryStateLzMax = this->LzMax;
  while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
    --this->TemporaryStateLzMax;
  return this->TargetSpace->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereShort::AdA (int index, int m)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  if (this->TemporaryStateLzMax < m)  
    return 0.0;
  return (double) (this->TemporaryState[m]);  
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double BosonOnSphereShort::AdA (long index, int m)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  if (this->TemporaryStateLzMax < m)  
    return 0.0;
  return (double) (this->TemporaryState[m]);  
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int BosonOnSphereShort::AdA (int index, int m, int n, double& coefficient)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  if ((this->TemporaryStateLzMax < n)  || (this->TemporaryState[n] == 0))
    { 
      coefficient = 0.0;
      return this->HilbertSpaceDimension;      
    }
  coefficient = (double) this->TemporaryState[n];
  --this->TemporaryState[n];
  if ((this->TemporaryStateLzMax == n) && (this->TemporaryState[n] == 0))
    {
      while (this->TemporaryState[this->TemporaryStateLzMax] == 0)
	--this->TemporaryStateLzMax;
    }
  if (this->TemporaryStateLzMax < m) 
    {
      for (int i = this->TemporaryStateLzMax + 1; i <= m; ++i)
	this->TemporaryState[i] = 0;
      this->TemporaryStateLzMax = m;
    }
  ++this->TemporaryState[m];
  coefficient *= (double) this->TemporaryState[m];
  coefficient = sqrt(coefficient);  
  int TmpIndex = this->TargetSpace->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
  if (TmpIndex == this->TargetSpace->HilbertSpaceDimension)
    coefficient = 0.0;
  return TmpIndex;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereShort::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[state], this->FermionBasis->StateLzMax[state], this->TemporaryState, this->TemporaryStateLzMax);
  int i = 0;
  for (; i <= this->TemporaryStateLzMax; ++i)
    Str << this->TemporaryState[i] << " ";
  for (; i <= this->LzMax; ++i)
    Str << "0 ";
  return Str;
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSphereShort::PrintStateMonomial (ostream& Str, int state)
{
  //  this->FermionToBoson(this->FermionBasis->StateDescription[state], this->FermionBasis->StateLzMax[state], this->TemporaryState, this->TemporaryStateLzMax);
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  //  this->ConvertToMonomial(this->TemporaryState, this->TemporaryStateLzMax, TmpMonomial);
  this->ConvertToMonomial(this->FermionBasis->StateDescription[state], this->FermionBasis->StateLzMax[state], TmpMonomial);
  Str << "[";
  if (TmpMonomial[0] != 0)
    Str << TmpMonomial[0];
  for (int i = 1; (i < this->NbrBosons) && (TmpMonomial[i] > 0); ++i)
    Str << "," << TmpMonomial[i];
  Str << "]";
  delete[] TmpMonomial;
  return Str;
}


// find state index from a string
//
// stateDescription = string describing the state
// return value = corresponding index, -1 if an error occured

int BosonOnSphereShort::FindStateIndex(char* stateDescription)
{
   char** TmpDescription;
  if (SplitLine(stateDescription, TmpDescription, ' ') != (this->LzMax + 1))
    return -1;
  int TmpNbrParticles = 0;
  int TmpTotalLz = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      int Tmp = atoi(TmpDescription[i]);
      this->TemporaryState[i] = Tmp;
      TmpTotalLz += (i * Tmp);
      TmpNbrParticles += Tmp;
      delete[] TmpDescription[i];
    }
  delete[] TmpDescription;
  if ((TmpNbrParticles != this->NbrBosons) || (TmpTotalLz != ((this->TotalLz + this->NbrBosons * this->LzMax) >> 1)))
    return -1;
  int NewLzMax = this->LzMax;
  while (this->TemporaryState[NewLzMax] == 0)
    --NewLzMax;
 return this->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, NewLzMax), NewLzMax + this->NbrBosons - 1);
}

// evaluate wave function in real space using a given basis and only for a given range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex BosonOnSphereShort::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis, 
						  int firstComponent, int nbrComponent)
{
  Complex Value;
  Complex Tmp;
  ComplexMatrix Perm(this->NbrBosons, this->NbrBosons);
  ComplexMatrix Functions(this->LzMax + 1, this->NbrBosons);
  RealVector TmpCoordinates(2);
  int* Indices = new int [this->NbrBosons];
  int Pos;
  int Lz;
  for (int j = 0; j < this->NbrBosons; ++j)
    {
      TmpCoordinates[0] = position[j << 1];
      TmpCoordinates[1] = position[1 + (j << 1)];
      for (int i = 0; i <= this->LzMax; ++i)
	{
	  basis.GetFunctionValue(TmpCoordinates, Tmp, i);
	  Functions[j].Re(i) = Tmp.Re;
	  Functions[j].Im(i) = Tmp.Im;
	}
    }
  double* Factors = new double [this->NbrBosons + 1];
  Factors[0] = 1.0;
  Factors[1] = 1.0;
  for (int i = 2; i <= this->NbrBosons; ++i)
    Factors[i] = Factors[i - 1] / sqrt((double) i);
  double TmpFactor;
  int* ChangeBitSign;
  int* ChangeBit;
  int TmpStateDescription;
  Perm.EvaluateFastPermanentPrecalculationArray(ChangeBit, ChangeBitSign);
  int LastComponent = firstComponent + nbrComponent;
  for (int k = firstComponent; k < LastComponent; ++k)
    {
      Pos = 0;
      Lz = 0;
      TmpFactor = state[k] * Factors[this->NbrBosons];
      this->FermionToBoson(this->FermionBasis->StateDescription[k], this->FermionBasis->StateLzMax[k], this->TemporaryState, this->TemporaryStateLzMax);
      while (Pos < this->NbrBosons)
	{
	  TmpStateDescription = this->TemporaryState[Lz];
	  if (TmpStateDescription != 0)
	    {
	      TmpFactor *= Factors[TmpStateDescription];
	      for (int j = 0; j < TmpStateDescription; ++j)
		{
		  Indices[Pos] = Lz;
		  ++Pos;
		}
	    }
	  ++Lz;
	}
      for (int i = 0; i < this->NbrBosons; ++i)
	{
	  ComplexVector& TmpColum2 = Functions[i];	  
	  for (int j = 0; j < this->NbrBosons; ++j)
	    {
	      Perm[i].Re(j) = TmpColum2.Re(Indices[j]);
	      Perm[i].Im(j) = TmpColum2.Im(Indices[j]);
	    }
	}
      Value += Perm.FastPermanent(ChangeBit, ChangeBitSign) * TmpFactor;
    }
  delete[] ChangeBitSign;
  delete[] ChangeBit;
  delete[] Factors;
  delete[] Indices;
  return Value;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  BosonOnSphereShort::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, RealVector& groundState)
{  
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrBosonSector == 0))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  if (subsytemSize > this->LzMax)
    {
      if ((lzSector == this->TotalLz) && (nbrBosonSector == this->NbrBosons))
	{
	  RealSymmetricMatrix TmpDensityMatrix(this->HilbertSpaceDimension);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    for (int j = i; j < this->HilbertSpaceDimension; ++j)
	      TmpDensityMatrix.SetMatrixElement(i, j, groundState[i] * groundState[j]);
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  if ((ShiftedLzComplementarySector < (NbrBosonsComplementarySector * subsytemSize)) || (ShiftedLzComplementarySector > (NbrBosonsComplementarySector * (this->LzMax))))
    {
      RealSymmetricMatrix TmpDensityMatrix;
      return TmpDensityMatrix;	  
    }

  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0.0;
 	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  unsigned long  TmpState2 = 0x0;
	  for (int i = 0; i < nbrBosonSector; ++i)
	    TmpState2 |= 0x1ul << i;
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | TmpState2;
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		TmpValue += groundState[TmpPos] * groundState[TmpPos];	
	    }

	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}      
    }
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  double TmpValue = 0;
 	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		TmpValue += groundState[TmpPos] * groundState[TmpPos];	
	    }	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int MinIndex = 0;
  // int MaxIndex = this->HilbertSpaceDimension - 1;
  if (nbrBosonSector == 1)
    {
      double TmpValue = 0.0;
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | (0x1ul << ShiftedLzSector);
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    TmpValue += groundState[TmpPos] * groundState[TmpPos];	
	}
      RealSymmetricMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
      return TmpDensityMatrix;
    }
  if (NbrBosonsComplementarySector == 0)
    {
      if (ShiftedLzComplementarySector != 0)
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
      BosonOnSphere TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
      MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      double TmpValue;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  TmpValue = groundState[MinIndex + i];
	  for (int j = i; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    TmpDensityMatrix.SetMatrixElement(i, j, TmpValue * groundState[MinIndex + j]);
	}
      return TmpDensityMatrix;
    }


  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;

  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpComplementaryState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]]);
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// reconstruct a state that contains only a certain subset of Schmidt eigenvalues of the given ground state
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// eigenvalueCut = only keep Schmidt levels that are larger than e^{-eigenvalueCut}
// groundState = reference on the total system ground state
// rebuiltSchmitGroundState = reference on the final state
// diagonalizedDensityMatrix = set of density matrix (Schmidt) eigenvalues
// transformationMatrix = the "truncation" matrix that connects the coefficients  (in the N-body basis) of the ground state and the final (truncated) state
// return value = reconstructed ground state vector

RealVector& BosonOnSphereShort::EvaluatePartialSchmidtDecomposition (int subsytemSize, int nbrBosonSector, int lzSector, double eigenvalueCut, RealVector& groundState, RealVector& rebuiltSchmidtGroundState, RealDiagonalMatrix& diagonalizedDensityMatrix, RealMatrix& transformationMatrix)
{
  if (subsytemSize <= 0)
    {
      cout << "Requested zero orbitals in the subsystem... exiting" << endl;
      return rebuiltSchmidtGroundState;
    }
  if (subsytemSize > this->LzMax)
    {
      cout << "Requested the entire system... exiting" << endl;
      for (long i = 0l; i < groundState.GetLargeVectorDimension(); ++i)
	rebuiltSchmidtGroundState[i] = groundState[i];
      return rebuiltSchmidtGroundState;
    }

  //Evaluate how many eigenvalues to keep 
  int NbrKeptEigenvalues = 0;  
  for (int i = 0; i < diagonalizedDensityMatrix.GetNbrRow(); ++i)
    if (diagonalizedDensityMatrix[i] >=  eigenvalueCut)
      ++NbrKeptEigenvalues;
  cout << "Keeping "<<NbrKeptEigenvalues<<" / "<<diagonalizedDensityMatrix.GetNbrRow()<<" eigenvalues"<<endl;
  if (NbrKeptEigenvalues == 0)
    return rebuiltSchmidtGroundState;
 
 
  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  if ((ShiftedLzComplementarySector < (NbrBosonsComplementarySector * subsytemSize)) || (ShiftedLzComplementarySector > (NbrBosonsComplementarySector * (this->LzMax))))
    {
      return rebuiltSchmidtGroundState;
    }

  //Only one orbital is kept in the subsystem, lzSector has to be 0
  //and reduced density matrix has only one element
  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
 	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  unsigned long  TmpState2 = 0x0;
	  for (int i = 0; i < nbrBosonSector; ++i)
	    TmpState2 |= 0x1ul << i;
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | TmpState2;
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		  rebuiltSchmidtGroundState[TmpPos] += groundState[TmpPos];
	    }
         return rebuiltSchmidtGroundState;
	}
      else
	{
           return rebuiltSchmidtGroundState;
	}      
    }

  //Subsystem contains 0 bosons, lzSector can therefore only be 0
  //and the reduced density matrix contains only one element
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		  rebuiltSchmidtGroundState[TmpPos] += groundState[TmpPos];	
	    }	  
         return rebuiltSchmidtGroundState;
	}
      else
	{
         return rebuiltSchmidtGroundState;
	}
    }

  //Case when there is only 1 boson in the subsystem A
  //Reduced density matrix only has one element
  if (nbrBosonSector == 1)
    {
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | (0x1ul << ShiftedLzSector);
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
   	      rebuiltSchmidtGroundState[TmpPos] += groundState[TmpPos];	
	}
     return rebuiltSchmidtGroundState;
    }

  //Subsystem A contains all the bosons
  if (NbrBosonsComplementarySector == 0)
    {
      if (ShiftedLzComplementarySector != 0)
	{
         return rebuiltSchmidtGroundState;
	}
      BosonOnSphere TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
      cout << " N_A=N, subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;

      //Build the "truncation" matrix
      RealMatrix TmpMatrix (TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension);
      for (int i = 0; i < TmpMatrix.GetNbrRow(); ++i)
        for (int j = 0; j < TmpMatrix.GetNbrRow(); ++j)
          {
            double Tmp = 0.0;
	    for (int k = 0; k < NbrKeptEigenvalues; ++k)
	      Tmp += transformationMatrix(i, k) * transformationMatrix(j, k);
	    TmpMatrix(i, j) = Tmp;
          }

      int MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  double Tmp = 0.0;
	  for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
              Tmp += groundState[MinIndex + j] * TmpMatrix(i, j);
   	  rebuiltSchmidtGroundState[MinIndex + i] = Tmp;	
	}
      return rebuiltSchmidtGroundState;
    }


  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  cout << "N_A<N, subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;

  RealMatrix TmpMatrix (TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpDestinationHilbertSpace.HilbertSpaceDimension);
  for (int i = 0; i < TmpMatrix.GetNbrRow(); ++i)
    for (int j = 0; j < TmpMatrix.GetNbrRow(); ++j)
      {
	double Tmp = 0.0;
	for (int k = 0; k < NbrKeptEigenvalues; ++k)
	  Tmp += transformationMatrix(i, k) * transformationMatrix(j, k);
	TmpMatrix(i, j) = Tmp;
      }


  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpComplementaryState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
              double Tmp = 0.0;
	      for (int k = 0; k < Pos; ++k)
                    Tmp += groundState[TmpStatePosition[k]] * TmpMatrix(Pos2, TmpStatePosition2[k]); 
              rebuiltSchmidtGroundState[TmpStatePosition[j]] = Tmp;
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  return rebuiltSchmidtGroundState;  
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particles
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax to -Lzmax+subsytemSize-1)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem

HermitianMatrix BosonOnSphereShort::EvaluatePartialDensityMatrix (int subsytemSize, int nbrBosonSector, int lzSector, ComplexVector& groundState)
{
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrBosonSector == 0))
	{
	  HermitianMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  if (subsytemSize > this->LzMax)
    {
      if ((lzSector == this->TotalLz) && (nbrBosonSector == this->NbrBosons))
	{
	  HermitianMatrix TmpDensityMatrix(this->HilbertSpaceDimension);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    for (int j = i; j < this->HilbertSpaceDimension; ++j)
	      TmpDensityMatrix.SetMatrixElement(i, j, Conj(groundState[i]) * groundState[j]);
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedLzSector = (lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1;
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  if ((ShiftedLzComplementarySector < (NbrBosonsComplementarySector * subsytemSize)) || (ShiftedLzComplementarySector > (NbrBosonsComplementarySector * (this->LzMax))))
    {
      HermitianMatrix TmpDensityMatrix;
      return TmpDensityMatrix;	  
    }

  if (subsytemSize == 1)
    {
      if (lzSector == 0)
	{
	  Complex TmpValue = 0.0;
 	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  unsigned long  TmpState2 = 0x0;
	  for (int i = 0; i < nbrBosonSector; ++i)
	    TmpState2 |= 0x1ul << i;
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | TmpState2;
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		TmpValue += Conj(groundState[TmpPos]) * groundState[TmpPos];	
	    }

	  HermitianMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}      
    }
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  Complex TmpValue = 0;
 	  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
	  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	    {
	      unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
	      int TmpLzMax = this->FermionBasis->LzMax;
	      while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
		--TmpLzMax;
	      int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	      if (TmpPos != this->HilbertSpaceDimension)
		TmpValue += Conj(groundState[TmpPos]) * groundState[TmpPos];	
	    }	  HermitianMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int MinIndex = 0;
  // int MaxIndex = this->HilbertSpaceDimension - 1;
  if (nbrBosonSector == 1)
    {
      Complex TmpValue = 0.0;
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  unsigned long TmpState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector) | (0x1ul << ShiftedLzSector);
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    TmpValue += Conj(groundState[TmpPos]) * groundState[TmpPos];	
	}
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);	    
      return TmpDensityMatrix;
    }
  if (NbrBosonsComplementarySector == 0)
    {
      if (ShiftedLzComplementarySector != 0)
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
      BosonOnSphere TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
      cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
      HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
      MinIndex = this->HilbertSpaceDimension - TmpDestinationHilbertSpace.HilbertSpaceDimension;
      Complex TmpValue;
      for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
	{
	  TmpValue = Conj(groundState[MinIndex + i]);
	  for (int j = i; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	    TmpDensityMatrix.SetMatrixElement(i, j, TmpValue * groundState[MinIndex + j]);
	}
      return TmpDensityMatrix;
    }


  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  HermitianMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;

  BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - nbrBosonSector, 2 * ShiftedLzComplementarySector - ((this->NbrBosons - nbrBosonSector) * (this->LzMax + subsytemSize)), this->LzMax - subsytemSize);

  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      unsigned long TmpComplementaryState = TmpHilbertSpace.FermionBasis->StateDescription[MinIndex] << (subsytemSize + nbrBosonSector);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpState = TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] | TmpComplementaryState;
	  int TmpLzMax = this->FermionBasis->LzMax;
	  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]]);
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]]);
		  }
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}
  
// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given Lz sector and fixed number of particle. The geometrical cut is a stripe.
// 
// subsytemSize = number of states that belong to the subsytem (ranging from -Lzmax+shitedCut to -Lzmax+shitedCut+subsytemSize-1)
// shiftedCut = first orbital belonging to the subsystem (with angular momentum -Lzmax+shitedCut)
// nbrBosonSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// lzSector = Lz sector in which the density matrix has to be evaluated 
// return value = density matrix of the subsytem  (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnSphereShort::EvaluateShiftedPartialDensityMatrix (int subsytemSize, int nbrShiftedOrbitals, int nbrBosonSector, int lzSector, RealVector& groundState)
{
  if (nbrShiftedOrbitals == 0)
    return this->EvaluatePartialDensityMatrix(subsytemSize, nbrBosonSector, lzSector, groundState);
  if (subsytemSize <= 0)
    {
      if ((lzSector == 0) && (nbrBosonSector == 0))
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }
  if (subsytemSize > this->LzMax)
    {
      RealSymmetricMatrix TmpDensityMatrix;
      return TmpDensityMatrix;	  
    }

  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  int ShiftedLzSector = ((lzSector + nbrBosonSector * (subsytemSize - 1)) >> 1) + (nbrBosonSector * nbrShiftedOrbitals);
  int ShiftedLzComplementarySector = ShiftedTotalLz - ShiftedLzSector;
  int NbrBosonsComplementarySector = this->NbrBosons - nbrBosonSector;
  if ((ShiftedLzComplementarySector < 0) || (ShiftedLzComplementarySector > (NbrBosonsComplementarySector * (this->LzMax))))
    {
      RealSymmetricMatrix TmpDensityMatrix;
      return TmpDensityMatrix;	  
    }

  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, subsytemSize - 1);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  long TmpNbrNonZeroElements = 0;

  for (int NbrBosonsLeft = 0; NbrBosonsLeft <= NbrBosonsComplementarySector; ++NbrBosonsLeft)
    {
      int MaxShiftedLzLeft = NbrBosonsLeft * (nbrShiftedOrbitals - 1);
      for (int ShiftedLzLeft = 0; ShiftedLzLeft <= MaxShiftedLzLeft; ++ShiftedLzLeft)
	{
	  int NbrBosonsRight = NbrBosonsComplementarySector - NbrBosonsLeft;
	  int ShiftedLzRight = ShiftedLzComplementarySector - ShiftedLzLeft;
	  if ((ShiftedLzRight >= (NbrBosonsRight * (nbrShiftedOrbitals + subsytemSize))) && (ShiftedLzRight <= (NbrBosonsRight * this->LzMax)))
	    {
	      BosonOnSphereShort TmpHilbertSpaceLeft(NbrBosonsLeft, 2 * ShiftedLzLeft - (NbrBosonsLeft * (nbrShiftedOrbitals - 1)), 
						     nbrShiftedOrbitals - 1);
	      BosonOnSphereShort TmpHilbertSpaceRight(NbrBosonsRight, 2 * (ShiftedLzRight - (NbrBosonsRight * (nbrShiftedOrbitals + subsytemSize))) - (NbrBosonsRight * (this->LzMax - nbrShiftedOrbitals - subsytemSize)), this->LzMax - nbrShiftedOrbitals - subsytemSize);	      
	      for (int MinIndexLeft = 0; MinIndexLeft < TmpHilbertSpaceLeft.HilbertSpaceDimension; ++MinIndexLeft)    
		{
		  unsigned long TmpComplementaryStateLeft = TmpHilbertSpaceLeft.FermionBasis->StateDescription[MinIndexLeft];
		  for (int MinIndexRight = 0; MinIndexRight < TmpHilbertSpaceRight.HilbertSpaceDimension; ++MinIndexRight)    
		    {
		      int Pos = 0;
		      unsigned long TmpComplementaryState = (TmpHilbertSpaceRight.FermionBasis->StateDescription[MinIndexRight] << (subsytemSize + nbrBosonSector + nbrShiftedOrbitals + NbrBosonsLeft)) | TmpComplementaryStateLeft;
		      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
			{
			  unsigned long TmpState = ((TmpDestinationHilbertSpace.FermionBasis->StateDescription[j] << (nbrShiftedOrbitals + NbrBosonsLeft))
						    | TmpComplementaryState);
			  int TmpLzMax = this->FermionBasis->LzMax;
			  while (((TmpState >> TmpLzMax) & 0x1ul) == 0x0ul)
			    --TmpLzMax;
			  int TmpPos = this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
			  if (TmpPos != this->HilbertSpaceDimension)
			    {
			      TmpStatePosition[Pos] = TmpPos;
			      TmpStatePosition2[Pos] = j;
			      ++Pos;
			    }
			}
		      if (Pos != 0)
			{
			  ++TmpNbrNonZeroElements;
			  for (int j = 0; j < Pos; ++j)
			    {
			      int Pos2 = TmpStatePosition2[j];
			      double TmpValue = groundState[TmpStatePosition[j]];
			      for (int k = 0; k < Pos; ++k)
				if (TmpStatePosition2[k] >= Pos2)
				  TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]]);
			    }
			}
		    }
		}
	    }
	}
    }
  delete[] TmpStatePosition2;
  delete[] TmpStatePosition;
  if (TmpNbrNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// Compute the row and column dimension of the orbital entanglement matrix of 2 cuts counting only those rows/columns that are not completely zero
// Also returns from the set of indices in the reduced density matrix corresponding to rows with atleast one non-zero entry
// Columns contain the hilbert space of B and C which are traced out
// SizeB = number of orbitals in part B, i.e. in the cap around Lzmax/2.
// SizeA = number of orbitals in the bulk of the sphere 
// NbrBosonsA = number of particles that belong to A
// groundState = reference on the total system ground state
// LzA = Lz sector of A in which the density matrix has to be evaluated as measured on a sphere with only A
// return value = pointer with the 1st element being the row dimension
//                2nd element is the column dimension of the oem 
//                3rd element onward gives the positions of the rows in the oem/reduced density matrix which are not filled with 0's
//                (returns 0 if there is a probem/there is no hilbert space)

long* BosonOnSphereShort::Compute2CutEntanglementMatrixDimensions (int SizeB, int SizeA, int NbrBosonsA, int LzA, RealVector& groundState)
{
  //Lz of full ground state on a disk
  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  //Lz of A on a disk
  int ShiftedLzA = (LzA + NbrBosonsA * (SizeA-1) ) >> 1;
  
  
  long DimA;
  int SizeC = this->LzMax +1 -SizeA -SizeB;
  
  BosonOnSphereShort TmpHilbertSpaceA(NbrBosonsA, LzA, SizeA - 1);
  
  if(NbrBosonsA>0)
    {
      DimA = TmpHilbertSpaceA.HilbertSpaceDimension;
    }
  else 
    if(NbrBosonsA==0)
      DimA=1;
    else
      {
	cout<<"The number of bosons in A is negative!";
	return 0; 
      }
  int* TmpIsSqueezedStateA = new int [DimA];
  //Initialize it
  
  long* Dims = new long [2+DimA];
  long PosninDims = 2;
  Dims[0] =0;
  Dims[1] =0;
  for(long i=0; i<DimA; i++)
    {
      TmpIsSqueezedStateA[i] = 0;
      Dims[2+i] =0;
    }
  for(int NbrBosonsB=0; NbrBosonsB <= (this->NbrBosons - NbrBosonsA); NbrBosonsB++)
    {
      int ShiftedLzBMax= (SizeB-1)*NbrBosonsB;
      int NbrBosonsC = this->NbrBosons - NbrBosonsA - NbrBosonsB;
      int ShiftedLzCMin = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - ShiftedLzBMax - (SizeA+SizeB)*NbrBosonsC;
      if(ShiftedLzCMin > (SizeC-1)*NbrBosonsC)
	continue;
      int ShiftedLzCMax = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - (SizeA+SizeB)*NbrBosonsC;
      if(ShiftedLzCMax<0)
	continue;
      int ShiftedLzBTemp = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - (this->LzMax)*NbrBosonsC;
      int ShiftedLzBMin = (ShiftedLzBTemp > 0) ? ShiftedLzBTemp : 0;
      int ShiftedLzBPossMax = (ShiftedLzCMax < ShiftedLzBMax) ? ShiftedLzCMax : ShiftedLzBMax;
      for(int ShiftedLzB=ShiftedLzBMin; ShiftedLzB <= ShiftedLzBPossMax; ShiftedLzB++)
	{
	  //2*Lz of B when only B is on a sphere
	  long DimB, DimC;
	  int SphereLzB = (ShiftedLzB<<1) - ( NbrBosonsB* (SizeB-1) );
	  
	  //Lz of C on a disk  
	  int ShiftedLzC = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - ShiftedLzB - (SizeA+SizeB)*NbrBosonsC;
	  
	  if( (ShiftedLzC<0) || (ShiftedLzC > (SizeC-1)*NbrBosonsC) )
	    {
	      //cout<<"No part C on a disk possible for values of NbrBosonsB="<<NbrBosonsB<<" and ShiftedLzB="<<ShiftedLzB<<endl;
	      continue;
	    }
	  
	  BosonOnSphereShort *TmpHilbertSpaceB=0;
	  if(NbrBosonsB>0)
	    {
	      TmpHilbertSpaceB= new BosonOnSphereShort (NbrBosonsB, SphereLzB, SizeB - 1 );
	      DimB = TmpHilbertSpaceB->GetHilbertSpaceDimension();
	    }
	  else 
	    if(NbrBosonsB==0)
	      DimB=1;
	    else
	      {
		cout<<"The number of bosons in B is negative!";
		return 0; 
	      }
	  
	  
	  //2*Lz of C when only C is on a sphere  
	  int SphereLzC = (ShiftedLzC<<1) - ( NbrBosonsC* (SizeC-1) );
	  BosonOnSphereShort *TmpHilbertSpaceC = 0;  
	  
	  if(NbrBosonsC>0)
	    {
	      TmpHilbertSpaceC = new BosonOnSphereShort (NbrBosonsC, SphereLzC, SizeC - 1 );
	      DimC = TmpHilbertSpaceC->GetHilbertSpaceDimension();
	    }
	  else 
	    if(NbrBosonsC==0)
	      DimC=1;
	    else
	      {
		cout<<"The number of bosons in C is negative!";
		return 0; 
	      }
	  for ( long i=0; i<DimB ; ++i)
	    {
	      for ( long j=0; j<DimC ; ++j)
		{
		  //Shifting the C Hilbert space element left
		  unsigned long TmpStateC, TmpStateB;
		  if(NbrBosonsC>0)
		    TmpStateC= TmpHilbertSpaceC->FermionBasis->StateDescription[j] << (SizeA + SizeB + NbrBosonsA + NbrBosonsB);
		  else
		    TmpStateC = 0;
		  
		  if(NbrBosonsB>0)	
		    TmpStateB = TmpHilbertSpaceB->FermionBasis->StateDescription[i];
		  else
		    TmpStateB = 0;
		  int Pos = 0;
		  for ( long k=0; k<DimA ; ++k)
		    {
		      unsigned long TmpStateA;
		      if(NbrBosonsA>0)	
			TmpStateA = TmpHilbertSpaceA.FermionBasis->StateDescription[k] << (SizeB + NbrBosonsB);
		      else
			TmpStateA = 0;
		      
		      unsigned long TmpFullState = (TmpStateB | TmpStateA) | TmpStateC;
		      
		      int TmpLzMax = this->FermionBasis->LzMax;
		      while (((TmpFullState >> TmpLzMax) & 0x1ul) == 0x0ul)
			--TmpLzMax;
		      int TmpPos = this->FermionBasis->FindStateIndex(TmpFullState, TmpLzMax);
		      
		      if ((TmpPos != this->HilbertSpaceDimension) && ( groundState[TmpPos] != 0) )
			{
			  if(Pos ==0)
			    {
			      Dims[1]++;
			      Pos++;
			    }
			  if(TmpIsSqueezedStateA[k] == 0)
			    {  
			      TmpIsSqueezedStateA[k] = 1;
			      Dims[PosninDims] = k;
			      PosninDims++;
			      Dims[0]++;
			    }
			}
		    }
		  
		}
            }
	  
	  //Finishing the braces of the big loop running over possible NbrParticlesB and ShiftedLzB		
	  if(NbrBosonsC>0)
	    delete TmpHilbertSpaceC;
          if(NbrBosonsB>0)
	    delete TmpHilbertSpaceB;
	  
	}
    }
  
  delete [] TmpIsSqueezedStateA; 
  return Dims;            
}
  			
  			
// evaluate a density matrix with 2 cuts of the whole system described by the RealVector groundState. The reduced density matrix is evaluated for a given Lz sector and number of particles
//
// SizeB = number of orbitals in part B, i.e. in the cap around Lzmax/2.
// SizeA = number of orbitals in the bulk of the sphere 
// NbrBosonsA = number of particles that belong to A
// groundState = reference on the total system ground state
// LzA = Lz sector of A in which the density matrix has to be evaluated as measured on a sphere with only A
// return value = density matrix of the subsytem (return a zero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnSphereShort::Evaluate2CutReducedDensityMatrix (int SizeB, int SizeA, int NbrBosonsA, int LzA, RealVector& groundState)
{
  //Lz of full ground state on a disk
  int ShiftedTotalLz = (this->TotalLz + this->NbrBosons * this->LzMax) >> 1;
  //Lz of A on a disk
  int ShiftedLzA = (LzA + NbrBosonsA * (SizeA-1) ) >> 1;
  //Make a Hilbert Space for A of dimension DimA only if there are some bosons in it
  long DimA;
  int SizeC = this->LzMax +1 -SizeA -SizeB;
  
  BosonOnSphereShort TmpHilbertSpaceA(NbrBosonsA, LzA, SizeA - 1);
  
  if(NbrBosonsA > 0)
    {
      DimA = TmpHilbertSpaceA.HilbertSpaceDimension;
    }
  else 
    if(NbrBosonsA==0)
      DimA=1;
    else
      {
	cout<<"The number of bosons in A is negative!";
	RealSymmetricMatrix TmpDensityMatrixZero;
	return TmpDensityMatrixZero; 
      }
  RealSymmetricMatrix TmpDensityMatrix(DimA, true);
  
  int* TmpFullStatePosition = new int [DimA];
  int* TmpStateAPosition = new int [DimA];
  int TmpNonZeroElements = 0;
  
  for(int NbrBosonsB=0; NbrBosonsB <= (this->NbrBosons - NbrBosonsA); NbrBosonsB++)
    {
      int ShiftedLzBMax= (SizeB-1)*NbrBosonsB;
      int NbrBosonsC = this->NbrBosons - NbrBosonsA - NbrBosonsB;
      int ShiftedLzCMin = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - ShiftedLzBMax - (SizeA+SizeB)*NbrBosonsC;
      if(ShiftedLzCMin > (SizeC-1)*NbrBosonsC)
	continue;
      int ShiftedLzCMax = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - (SizeA+SizeB)*NbrBosonsC;
      if(ShiftedLzCMax<0)
	continue;
      
      for(int ShiftedLzB=0; ShiftedLzB <= ShiftedLzBMax; ShiftedLzB++)
	{
	  //2*Lz of B when only B is on a sphere
	  long DimB, DimC;
	  int SphereLzB = (ShiftedLzB<<1) - ( NbrBosonsB* (SizeB-1) );
	  
	  //Lz of C on a disk  
	  int ShiftedLzC = ShiftedTotalLz - (ShiftedLzA + SizeB*NbrBosonsA) - ShiftedLzB - (SizeA+SizeB)*NbrBosonsC;
	  
	  if( (ShiftedLzC<0) || (ShiftedLzC > (SizeC-1)*NbrBosonsC) )
	    {
	      cout<<"No part C on a disk possible for values of NbrBosonsB="<<NbrBosonsB<<" and ShiftedLzB="<<ShiftedLzB<<endl;
	      continue;
	    }
	  BosonOnSphereShort *TmpHilbertSpaceB=0;
	  if(NbrBosonsB>0) 
	    {
	      TmpHilbertSpaceB= new BosonOnSphereShort (NbrBosonsB, SphereLzB, SizeB - 1 );
	      DimB = TmpHilbertSpaceB->GetHilbertSpaceDimension();
	    }
	  else 
	    if(NbrBosonsB==0)
	      DimB=1;
	    else
	      {
		cout<<"The number of bosons in B is negative!";
		RealSymmetricMatrix TmpDensityMatrixZero;
		return TmpDensityMatrixZero; 
	      }
	  
	  //2*Lz of C when only C is on a sphere  
	  int SphereLzC = (ShiftedLzC<<1) - ( NbrBosonsC* (SizeC-1) );
	  BosonOnSphereShort *TmpHilbertSpaceC = 0;  
	  
	  if(NbrBosonsC>0)
	    {
	      TmpHilbertSpaceC = new BosonOnSphereShort (NbrBosonsC, SphereLzC, SizeC - 1 );
	      DimC = TmpHilbertSpaceC->GetHilbertSpaceDimension();
	    }
	  else 
	    if(NbrBosonsC==0)
	      DimC=1;
	    else
	      {
		cout<<"The number of bosons in C is negative!";
		RealSymmetricMatrix TmpDensityMatrixZero;
		return TmpDensityMatrixZero; 
	      }
	  for (long i = 0; i < DimB ; ++i)
	    {
	      for ( long j = 0; j < DimC ; ++j)
		{
		  //Shifting the C Hilbert space element left
		  unsigned long TmpStateC, TmpStateB;
		  if(NbrBosonsC>0)
		    TmpStateC= TmpHilbertSpaceC->FermionBasis->StateDescription[j] << (SizeA + SizeB + NbrBosonsA + NbrBosonsB);
		  else
		    TmpStateC = 0;
		  
		  if(NbrBosonsB>0)	
		    TmpStateB = TmpHilbertSpaceB->FermionBasis->StateDescription[i];
		  else
		    TmpStateB = 0;
		  int Pos=0;
		  
		  for ( long k=0; k<DimA ; ++k)
		    {
		      unsigned long TmpStateA;
		      if(NbrBosonsA>0)	
			TmpStateA = TmpHilbertSpaceA.FermionBasis->StateDescription[k] << (SizeB + NbrBosonsB);
		      else
			TmpStateA = 0;
		      
		      unsigned long TmpFullState = (TmpStateB | TmpStateA) | TmpStateC;
		      
		      int TmpLzMax = this->FermionBasis->LzMax;
		      while (((TmpFullState >> TmpLzMax) & 0x1ul) == 0x0ul)
			--TmpLzMax;
		      int TmpPos = this->FermionBasis->FindStateIndex(TmpFullState, TmpLzMax);
		      
		      if (TmpPos != this->HilbertSpaceDimension)
			{
			  //this->PrintState(cout,TmpPos)<<" "<<groundState[TmpPos]<<" " << k<<endl;
			  TmpFullStatePosition[Pos] = TmpPos;
			  //cout<<TmpPos<<endl;
			  TmpStateAPosition[Pos] = k;
			  ++Pos;
			}
		    }
		  
		  if(Pos!=0)
		    {
		      //Then for the given B and C, we have some contribution to the reduced density matrix
		      ++TmpNonZeroElements;
		      for(long k=0; k<Pos; k++)
			{
			  long RowA=TmpStateAPosition[k];
			  double coeff_k=groundState[TmpFullStatePosition[k]];
			  for(long m=0; m<Pos; m++) 
			    if (TmpStateAPosition[m] >= RowA)
			      TmpDensityMatrix.AddToMatrixElement(RowA,TmpStateAPosition[m],coeff_k*groundState[TmpFullStatePosition[m]]);
			  //cout<<"trying to set something"<<coeff_k;
			}
		    }
		  
		}
	    }
	  //Finishing the braces of the big loop running over possible NbrParticlesB and ShiftedLzB		
  		  if(NbrBosonsC>0)
		    delete TmpHilbertSpaceC;
		  if(NbrBosonsB>0)
		    delete TmpHilbertSpaceB;
		  
	}
    }
  
  delete[] TmpStateAPosition;
  delete[] TmpFullStatePosition;
  if (TmpNonZeroElements > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix  BosonOnSphereShort::EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector& groundState, AbstractArchitecture* architecture)
{  
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  if (nbrBosonSector == this->NbrBosons)
    {
      if (lzSector == this->TotalLz)
	{
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, 1.0);
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int ComplementaryNbrBosonSector = this->NbrBosons - nbrBosonSector;

  if (nbrBosonSector == 1)
    {
      double TmpValue = 0.0;
      unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
      unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - 1, this->TotalLz - lzSector, this->LzMax);
      unsigned long ShiftedLzVSector = (lzSector + this->LzMax) >> 1;
      FactorialCoefficient Factorial;
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
	  TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateLzMax);
	  int TmpIndex2 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial1[TmpIndex2] >= ShiftedLzVSector))
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  TmpMonomial3[TmpIndex4] = ShiftedLzVSector;
	  ++TmpIndex4;
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      Factorial.SetToOne();
	      for (int k = 0; k <= TmpHilbertSpace.TemporaryStateLzMax; ++k)
		if (TmpHilbertSpace.TemporaryState[k] > 1)
		  Factorial.FactorialDivide(TmpHilbertSpace.TemporaryState[k]);
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		if (this->TemporaryState[k] > 1)
		  Factorial.FactorialMultiply(this->TemporaryState[k]);
	      Factorial.BinomialDivide(this->NbrBosons, nbrBosonSector);
	      TmpValue += groundState[TmpPos] * groundState[TmpPos] * (Factorial.GetNumericalValue());	
	    }
	}
      
      RealSymmetricMatrix TmpDensityMatrix(1);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
      return TmpDensityMatrix;
    }

  if (abs(this->TotalLz - lzSector) > (ComplementaryNbrBosonSector * LzMax))
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }

  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, this->LzMax);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  BosonOnSphereShort TmpHilbertSpace(ComplementaryNbrBosonSector, this->TotalLz - lzSector, this->LzMax);

  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &TmpDestinationHilbertSpace, &TmpHilbertSpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnSphereShort::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
 									    RealVector& groundState,  RealSymmetricMatrix* densityMatrix)
{
   BosonOnSphereShort* TmpHilbertSpace =  (BosonOnSphereShort*) complementaryHilbertSpace;
   BosonOnSphereShort* TmpDestinationHilbertSpace =  (BosonOnSphereShort*) destinationHilbertSpace;
   int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
   int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
   unsigned long* TmpMonomial2 = new unsigned long [NbrBosonSector];
   unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
   unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
   int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   int MaxIndex = minIndex + nbrIndex;
   long TmpNbrNonZeroElements = 0l;

   double* LogFactorials = new double[this->NbrBosons + 1];
   LogFactorials[0] = 0.0;
   LogFactorials[1] = 0.0;
   for (int i = 2 ; i <= this->NbrBosons; ++i)
     LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
   double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
   double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrBosonSector] - LogFactorials[NbrBosonSector];

   for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->FermionBasis->StateDescription[i], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[i], TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateLzMax);
      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->TemporaryStateLzMax; ++k)
	TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState[k]];
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }
   for (; minIndex < MaxIndex; ++minIndex)    
     {
      int Pos = 0;
      TmpHilbertSpace->ConvertToMonomial(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpMonomial1);
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpHilbertSpace->TemporaryState, TmpHilbertSpace->TemporaryStateLzMax);
      double TmpHilbertSpaceFactorial = 0.0;
      for (int k = 0; k <= TmpHilbertSpace->TemporaryStateLzMax; ++k)
	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState[k]];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  TmpDestinationHilbertSpace->ConvertToMonomial(TmpDestinationHilbertSpace->FermionBasis->StateDescription[j], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[j], TmpMonomial2);
	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpIndex3 < NbrBosonSector)) 
	    {
	      while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementaryNbrBosonSector)
		{
		  while ((TmpIndex3 < NbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
		    {
		      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
		      ++TmpIndex3;
		      ++TmpIndex4;		  
		    }
		}
	    }
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  while (TmpIndex3 < NbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }

	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		TmpFactorial += LogFactorials[this->TemporaryState[k]];
	      TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
	      TmpFactorial *= 0.5; 
	      
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      TmpStateCoefficient[Pos] = exp(TmpFactorial);
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
		    ++TmpNbrNonZeroElements;
		  }
	    }
	}
     }
   delete[] TmpMonomial2;
   delete[] TmpMonomial1;
   delete[] TmpMonomial3;
   delete[] TmpStatePosition;
   delete[] TmpStatePosition2;
   delete[] TmpStateCoefficient;
   delete[] TmpDestinationLogFactorials;
   return TmpNbrNonZeroElements;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using real space partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
// thetaTop =  inclination angle defining one edge of the cut in degrees
// thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in degrees
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnSphereShort::EvaluatePartialDensityMatrixRealSpacePartition (int nbrBosonSector, int lzSector, double thetaTop, double thetaBottom, double phiRange, RealVector& groundState, AbstractArchitecture* architecture)
{
  if ((abs(lzSector) > (nbrBosonSector * this->LzMax)) || (thetaBottom <= thetaTop) || (phiRange <= 0.0))
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }

  thetaTop *= M_PI / 180.0;
  thetaBottom *= M_PI / 180.0;
  phiRange /= 360.0;
  
  double* IncompleteBetaThetaTop = 0;
  double* IncompleteBetaThetaBottom = 0;

  this->EvaluatePartialDensityMatrixRealSpacePartitionCoefficient(this->LzMax, thetaTop, thetaBottom, IncompleteBetaThetaTop, IncompleteBetaThetaBottom);

  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  unsigned long* TmpMonomial1 = new unsigned long [this->NbrBosons];
	  double TmpValue = 0.0;
	  for (int MinIndex = 0; MinIndex < this->HilbertSpaceDimension; ++MinIndex)    
	    {
	      this->ConvertToMonomial(this->FermionBasis->StateDescription[MinIndex], this->FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
	      double FormFactor = 0.0;
	      for (int i=0; i < this->NbrBosons; i++)
		FormFactor += log(1.0 - IncompleteBetaThetaBottom[TmpMonomial1[i]] + IncompleteBetaThetaTop[TmpMonomial1[i]] + (1.0 - phiRange) * (IncompleteBetaThetaBottom[TmpMonomial1[i]] - IncompleteBetaThetaTop[TmpMonomial1[i]]) );
	      FormFactor = exp(FormFactor);
	      TmpValue += groundState[MinIndex] * groundState[MinIndex] * FormFactor;	
	    }
	  RealSymmetricMatrix TmpDensityMatrix(1);
	  TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue); 
	  
	  delete[] IncompleteBetaThetaTop;
	  delete[] IncompleteBetaThetaBottom;
	  delete[] TmpMonomial1;
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  if (nbrBosonSector == this->NbrBosons)
    {
      if (lzSector == this->TotalLz)
	{
	  RealSymmetricMatrix TmpDensityMatrix(this->HilbertSpaceDimension, true);
	  unsigned long* TmpMonomial1 = new unsigned long [this->NbrBosons];
	  double* TmpStateCoefficient = new double [this->HilbertSpaceDimension];
	  for( int i = 0; i < this->HilbertSpaceDimension; i++)
	    {
	      TmpStateCoefficient[i] = 0.5 * this->NbrBosons * log(phiRange);
	      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial1);
	      for( int j=0; j<this->NbrBosons; j++)
		{
		  TmpStateCoefficient [i] += 0.5*log( IncompleteBetaThetaBottom[TmpMonomial1[j]] - IncompleteBetaThetaTop[TmpMonomial1[j]]);
		}
	      TmpStateCoefficient[i] = exp(TmpStateCoefficient[i]);
	    }
	  
	  for(int pos1 = 0; pos1 < this->HilbertSpaceDimension; pos1++)
	    for(int pos2 = pos1; pos2 < this->HilbertSpaceDimension; pos2++)
	      {
		TmpDensityMatrix.SetMatrixElement(pos1, pos2, groundState[pos1]*groundState[pos2]*TmpStateCoefficient[pos1]*TmpStateCoefficient[pos2]);
	      }
	  delete[] TmpMonomial1;
	  delete[] TmpStateCoefficient;
	  delete[] IncompleteBetaThetaTop;
	  delete[] IncompleteBetaThetaBottom;
	  return TmpDensityMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;	  
	}
    }

  int ComplementaryNbrBosonSector = this->NbrBosons - nbrBosonSector;

  if (nbrBosonSector == 1)
    {
      double TmpValue = 0.0;
      unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
      unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
      BosonOnSphereShort TmpHilbertSpace(this->NbrBosons - 1, this->TotalLz - lzSector, this->LzMax);
      unsigned long ShiftedLzVSector = (lzSector + this->LzMax) >> 1;
      double TmpStateCoefficient = phiRange * (IncompleteBetaThetaBottom[ShiftedLzVSector] - IncompleteBetaThetaTop[ShiftedLzVSector]);
      FactorialCoefficient Factorial;
      for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
	{
	  TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
	  TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateLzMax);
	  double FormFactor = 0.0;
	  for (int i = 0; i < ComplementaryNbrBosonSector; i++)
	    FormFactor += log(1.0 - IncompleteBetaThetaBottom[TmpMonomial1[i]] + IncompleteBetaThetaTop[TmpMonomial1[i]] + (1.0 - phiRange) * (IncompleteBetaThetaBottom[TmpMonomial1[i]] - IncompleteBetaThetaTop[TmpMonomial1[i]]));
	  FormFactor = exp(FormFactor);
	  int TmpIndex2 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial1[TmpIndex2] >= ShiftedLzVSector))
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  TmpMonomial3[TmpIndex4] = ShiftedLzVSector;
	  ++TmpIndex4;
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      Factorial.SetToOne();
	      for (int k = 0; k <= TmpHilbertSpace.TemporaryStateLzMax; ++k)
		if (TmpHilbertSpace.TemporaryState[k] > 1)
		  Factorial.FactorialDivide(TmpHilbertSpace.TemporaryState[k]);
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		if (this->TemporaryState[k] > 1)
		  Factorial.FactorialMultiply(this->TemporaryState[k]);
	      TmpValue += groundState[TmpPos] * groundState[TmpPos] * (Factorial.GetNumericalValue()) * FormFactor * TmpStateCoefficient;	
	    }
	}
      
      RealSymmetricMatrix TmpDensityMatrix(1);
      TmpDensityMatrix.SetMatrixElement(0, 0, TmpValue);
      delete[] IncompleteBetaThetaTop;
      delete[] IncompleteBetaThetaBottom;
      delete [] TmpMonomial1;
      delete [] TmpMonomial3;
      return TmpDensityMatrix;
    }

  if (abs(this->TotalLz - lzSector) > (ComplementaryNbrBosonSector * LzMax))
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }

  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, this->LzMax);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  RealSymmetricMatrix TmpDensityMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, true);
  BosonOnSphereShort TmpHilbertSpace(ComplementaryNbrBosonSector, this->TotalLz - lzSector, this->LzMax);

   FQHESphereParticleEntanglementSpectrumOperation Operation(this, &TmpDestinationHilbertSpace, &TmpHilbertSpace, groundState, TmpDensityMatrix, IncompleteBetaThetaBottom, IncompleteBetaThetaTop, phiRange);
   Operation.ApplyOperation(architecture);
   if (Operation.GetNbrNonZeroMatrixElements() > 0)	
     return TmpDensityMatrix;
   else
     {
       RealSymmetricMatrix TmpDensityMatrixZero;
       return TmpDensityMatrixZero;
     }


//   int* TmpStatePosition = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   double* TmpStateCoefficientOccFactorial = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   long TmpNbrNonZeroElements = 0;
//   unsigned long* TmpMonomial2 = 0;
//   unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
//   unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];

//   double* LogFactorials = new double[this->NbrBosons + 1];
//   LogFactorials[0] = 0.0;
//   LogFactorials[1] = 0.0;
//   for (int i = 2 ; i <= this->NbrBosons; ++i)
//     LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 

//   FactorialCoefficient Factorial;
//   unsigned long** TmpDestinationHilbertSpaceOccupationNumbers = new unsigned long* [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   unsigned long** TmpDestinationHilbertSpaceMonomial = new unsigned long* [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   int* TmpDestinationHilbertSpaceLzMax = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
//   for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
//     {
//       TmpDestinationHilbertSpaceOccupationNumbers[i] = new unsigned long [this->NbrLzValue];
//       TmpDestinationHilbertSpaceMonomial[i] = new unsigned long [nbrBosonSector];
//       TmpDestinationHilbertSpace.FermionToBoson(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpDestinationHilbertSpaceOccupationNumbers[i], TmpDestinationHilbertSpaceLzMax[i]);
//       TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpDestinationHilbertSpaceMonomial[i]);
//       unsigned long* TmpOccupationNumber = TmpDestinationHilbertSpaceOccupationNumbers[i];
//       int TmpLzMax = TmpDestinationHilbertSpaceLzMax[i];
//       double TmpFactor = 0.0;
//       for (int k = 0; k <= TmpLzMax; ++k)
// 	TmpFactor += LogFactorials[TmpOccupationNumber[k]];
//       TmpDestinationLogFactorials[i] =  TmpFactor;
//     }

//   //Compute the coefficients multiplying rhoA in TmpStateCoefficient
//   for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
//     {
//       TmpStateCoefficient [i] = 0.5 * nbrBosonSector * log(phiRange);
      
//       for( int j=0; j<nbrBosonSector; j++)
// 	{
// 	  TmpStateCoefficient [i] += 0.5 * log( IncompleteBetaThetaBottom[TmpDestinationHilbertSpaceMonomial[i][j]] - IncompleteBetaThetaTop[TmpDestinationHilbertSpaceMonomial[i][j]]);
// 	}
//       TmpStateCoefficient[i] = exp(TmpStateCoefficient[i]);
//     }

//   for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
//     {
//       int Pos = 0;
//       TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
//       double FormFactor = 0.0;
//       for (int i=0; i < ComplementaryNbrBosonSector; i++)
// 	FormFactor += log(1.0 - IncompleteBetaThetaBottom[TmpMonomial1[i]] + IncompleteBetaThetaTop[TmpMonomial1[i]] + (1.0 - phiRange) * (IncompleteBetaThetaBottom[TmpMonomial1[i]] - IncompleteBetaThetaTop[TmpMonomial1[i]]) );
//       FormFactor = exp(FormFactor);
//       TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateLzMax);
//       double TmpHilbertSpaceFactorial = 0.0;
//       for (int k = 0; k <= TmpHilbertSpace.TemporaryStateLzMax; ++k)
// 	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace.TemporaryState[k]];
//       for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
// 	{
// 	  TmpMonomial2 = TmpDestinationHilbertSpaceMonomial[j];
// 	  int TmpIndex2 = 0;
// 	  int TmpIndex3 = 0;
// 	  int TmpIndex4 = 0;
// 	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpIndex3 < nbrBosonSector)) 
// 	    {
// 	      while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
// 		{
// 		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
// 		  ++TmpIndex2;
// 		  ++TmpIndex4;		  
// 		}
// 	      if (TmpIndex2 < ComplementaryNbrBosonSector)
// 		{
// 		  while ((TmpIndex3 < nbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
// 		    {
// 		      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
// 		      ++TmpIndex3;
// 		      ++TmpIndex4;		  
// 		    }
// 		}
// 	    }
// 	  while (TmpIndex2 < ComplementaryNbrBosonSector)
// 	    {
// 	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
// 	      ++TmpIndex2;
// 	      ++TmpIndex4;		  
// 	    }
// 	  while (TmpIndex3 < nbrBosonSector)
// 	    {
// 	      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
// 	      ++TmpIndex3;
// 	      ++TmpIndex4;		  
// 	    }

// 	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
// 	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
// 	  if (TmpPos != this->HilbertSpaceDimension)
// 	    {	      
// 	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
// 	      double TmpFactorial = 0.0;	      
// 	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
// 		TmpFactorial += LogFactorials[this->TemporaryState[k]];
// 	      TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j];
// 	      TmpFactorial *= 0.5; 
	      
// 	      TmpStatePosition[Pos] = TmpPos;
// 	      TmpStatePosition2[Pos] = j;
// 	      TmpStateCoefficientOccFactorial[Pos] = exp(TmpFactorial);
// 	      ++Pos;
// 	    }
// 	}
//       if (Pos != 0)
// 	{
// 	  ++TmpNbrNonZeroElements;
// 	  for (int j = 0; j < Pos; ++j)
// 	    {
// 	      int Pos2 = TmpStatePosition2[j];
// 	      double TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient[Pos2] * TmpStateCoefficientOccFactorial[j];
// 	      for (int k = 0; k < Pos; ++k)
// 		if (TmpStatePosition2[k] >= Pos2)
// 		  {
// 		    TmpDensityMatrix.AddToMatrixElement(Pos2, TmpStatePosition2[k],  FormFactor * TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[TmpStatePosition2[k]] * TmpStateCoefficientOccFactorial[k]);
// 		  }
// 	    }
// 	}
//     }
//   for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
//     {
//       delete[] TmpDestinationHilbertSpaceOccupationNumbers[i];
//       delete[] TmpDestinationHilbertSpaceMonomial[i];
//     }
//   delete[] TmpDestinationHilbertSpaceOccupationNumbers;
//   delete[] TmpDestinationHilbertSpaceLzMax;
//   delete[] TmpDestinationHilbertSpaceMonomial;
//   delete[] TmpStatePosition2;
//   delete[] TmpStatePosition;
//   delete[] TmpStateCoefficient;
//   delete[] TmpMonomial1;
//   delete[] TmpMonomial3;
//   if (TmpNbrNonZeroElements > 0)	
//     return TmpDensityMatrix;
//   else
//     {
//       RealSymmetricMatrix TmpDensityMatrixZero;
//       return TmpDensityMatrixZero;
//     }
}

// core part of the evaluation density matrix real space partition calculation
// 
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// incompleteBetaThetaTop = pointer to the array where the top part coefficients are stored
// incompleteBetaThetaBotton = pointer on the pointer to the array where the bottom part coefficients are stored
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
// return value = number of components that have been added to the density matrix

long BosonOnSphereShort::EvaluatePartialDensityMatrixRealSpacePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
									     RealVector& groundState,  RealSymmetricMatrix* densityMatrix, double* incompleteBetaThetaBottom, double* incompleteBetaThetaTop, double phiRange)
{
  BosonOnSphereShort* TmpHilbertSpace =  (BosonOnSphereShort*) complementaryHilbertSpace;
  BosonOnSphereShort* TmpDestinationHilbertSpace =  (BosonOnSphereShort*) destinationHilbertSpace;
  int ComplementaryNbrBosonSector = TmpHilbertSpace->NbrBosons;
  int NbrBosonSector = TmpDestinationHilbertSpace->NbrBosons;
  unsigned long* TmpMonomial2 = new unsigned long [NbrBosonSector];
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
  unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double* TmpStateCoefficient = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  double* TmpStateCoefficientOccFactorial = new double [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;

  double* LogFactorials = new double[this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2 ; i <= this->NbrBosons; ++i)
    LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 


  for (int i = 0; i < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace->FermionToBoson(TmpDestinationHilbertSpace->FermionBasis->StateDescription[i], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[i], TmpDestinationHilbertSpace->TemporaryState, TmpDestinationHilbertSpace->TemporaryStateLzMax);
      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpDestinationHilbertSpace->TemporaryStateLzMax; ++k)
	TmpFactor += LogFactorials[TmpDestinationHilbertSpace->TemporaryState[k]];
      TmpDestinationLogFactorials[i] =  TmpFactor;

      TmpStateCoefficient[i] = 0.5 * NbrBosonSector * log(phiRange);
      
      TmpDestinationHilbertSpace->ConvertToMonomial(TmpDestinationHilbertSpace->FermionBasis->StateDescription[i], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[i], TmpMonomial2);
     for(int j = 0; j < NbrBosonSector; j++)
	{
	  TmpStateCoefficient[i] += 0.5 * log( incompleteBetaThetaBottom[TmpMonomial2[j]] - incompleteBetaThetaTop[TmpMonomial2[j]]);
	}
      TmpStateCoefficient[i] = exp(TmpStateCoefficient[i]);
    }

  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      TmpHilbertSpace->ConvertToMonomial(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpMonomial1);
      double FormFactor = 0.0;
      for (int i=0; i < ComplementaryNbrBosonSector; i++)
	FormFactor += log(1.0 - incompleteBetaThetaBottom[TmpMonomial1[i]] + incompleteBetaThetaTop[TmpMonomial1[i]] + (1.0 - phiRange) * (incompleteBetaThetaBottom[TmpMonomial1[i]] - incompleteBetaThetaTop[TmpMonomial1[i]]) );
      FormFactor = exp(FormFactor);
      TmpHilbertSpace->FermionToBoson(TmpHilbertSpace->FermionBasis->StateDescription[minIndex], TmpHilbertSpace->FermionBasis->StateLzMax[minIndex], TmpHilbertSpace->TemporaryState, TmpHilbertSpace->TemporaryStateLzMax);
      double TmpHilbertSpaceFactorial = 0.0;
      for (int k = 0; k <= TmpHilbertSpace->TemporaryStateLzMax; ++k)
	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace->TemporaryState[k]];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  TmpDestinationHilbertSpace->ConvertToMonomial(TmpDestinationHilbertSpace->FermionBasis->StateDescription[j], TmpDestinationHilbertSpace->FermionBasis->StateLzMax[j], TmpMonomial2);
	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpIndex3 < NbrBosonSector)) 
	    {
	      while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementaryNbrBosonSector)
		{
		  while ((TmpIndex3 < NbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
		    {
		      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
		      ++TmpIndex3;
		      ++TmpIndex4;		  
		    }
		}
	    }
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  while (TmpIndex3 < NbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }

	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {	      
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		TmpFactorial += LogFactorials[this->TemporaryState[k]];
	      TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j];
	      TmpFactorial *= 0.5; 
	      
	      TmpStatePosition[Pos] = TmpPos;
	      TmpStatePosition2[Pos] = j;
	      TmpStateCoefficientOccFactorial[Pos] = exp(TmpFactorial);
	      ++Pos;
	    }
	}
      if (Pos != 0)
	{
	  ++TmpNbrNonZeroElements;
	  for (int j = 0; j < Pos; ++j)
	    {
	      int Pos2 = TmpStatePosition2[j];
	      double TmpValue = groundState[TmpStatePosition[j]] * TmpStateCoefficient[Pos2] * TmpStateCoefficientOccFactorial[j];
	      for (int k = 0; k < Pos; ++k)
		if (TmpStatePosition2[k] >= Pos2)
		  {
		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k],  FormFactor * TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[TmpStatePosition2[k]] * TmpStateCoefficientOccFactorial[k]);
		  }
	    }
	}
    }
  delete[] TmpMonomial2;
  delete[] TmpMonomial1;
  delete[] TmpMonomial3;
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  delete[] TmpDestinationLogFactorials;
  delete[] TmpStateCoefficientOccFactorial;
  return TmpNbrNonZeroElements;
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

RealMatrix BosonOnSphereShort::EvaluatePartialEntanglementMatrixParticlePartition (int nbrBosonSector, int lzSector, RealVector& groundState, bool removeBinomialCoefficient)
{
  if (nbrBosonSector == 0)
    {
      if (lzSector == 0)
	{
	  RealMatrix TmpEntanglementMatrix(1, this->HilbertSpaceDimension, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    TmpEntanglementMatrix.SetMatrixElement(0, i, groundState[i]);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealSymmetricMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }

  if (nbrBosonSector == this->NbrBosons)
    {
      if (lzSector == this->TotalLz)
	{
	  RealMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    TmpEntanglementMatrix.SetMatrixElement(i, 0, groundState[i]);
	  return TmpEntanglementMatrix;
	}
      else
	{
	  RealMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;	  
	}
    }

  int ComplementaryNbrBosonSector = this->NbrBosons - nbrBosonSector;

  if (abs(this->TotalLz - lzSector) > (ComplementaryNbrBosonSector * LzMax))
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }

  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, this->LzMax);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  long TmpNbrNonZeroElements = 0;
  unsigned long* TmpMonomial2 = 0;
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
  unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];

  double* LogFactorials = new double[this->NbrBosons + 1];
  LogFactorials[0] = 0.0;
  LogFactorials[1] = 0.0;
  for (int i = 2 ; i <= this->NbrBosons; ++i)
    LogFactorials[i] = LogFactorials[i - 1] + log((double) i); 
  double TmpLogBinomial = LogFactorials[this->NbrBosons] - LogFactorials[ComplementaryNbrBosonSector] - LogFactorials[nbrBosonSector];
  if (removeBinomialCoefficient == true)
    TmpLogBinomial = 0.0;

  BosonOnSphereShort TmpHilbertSpace(ComplementaryNbrBosonSector, this->TotalLz - lzSector, this->LzMax);
  FactorialCoefficient Factorial;
  unsigned long** TmpDestinationHilbertSpaceOccupationNumbers = new unsigned long* [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  unsigned long** TmpDestinationHilbertSpaceMonomial = new unsigned long* [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  int* TmpDestinationHilbertSpaceLzMax = new int [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  double* TmpDestinationLogFactorials = new double [TmpDestinationHilbertSpace.HilbertSpaceDimension];
  RealMatrix TmpEntanglementMatrix(TmpDestinationHilbertSpace.HilbertSpaceDimension, TmpHilbertSpace.HilbertSpaceDimension, true);
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpaceOccupationNumbers[i] = new unsigned long [this->NbrLzValue];
      TmpDestinationHilbertSpaceMonomial[i] = new unsigned long [nbrBosonSector];
      TmpDestinationHilbertSpace.FermionToBoson(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpDestinationHilbertSpaceOccupationNumbers[i], TmpDestinationHilbertSpaceLzMax[i]);
      TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpDestinationHilbertSpaceMonomial[i]);
      unsigned long* TmpOccupationNumber = TmpDestinationHilbertSpaceOccupationNumbers[i];
      int TmpLzMax = TmpDestinationHilbertSpaceLzMax[i];
      double TmpFactor = 0.0;
      for (int k = 0; k <= TmpLzMax; ++k)
	TmpFactor += LogFactorials[TmpOccupationNumber[k]];
      TmpDestinationLogFactorials[i] =  TmpFactor;
    }
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      int Pos = 0;
      TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
      TmpHilbertSpace.FermionToBoson(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpHilbertSpace.TemporaryState, TmpHilbertSpace.TemporaryStateLzMax);
      double TmpHilbertSpaceFactorial = 0.0;
      for (int k = 0; k <= TmpHilbertSpace.TemporaryStateLzMax; ++k)
	TmpHilbertSpaceFactorial += LogFactorials[TmpHilbertSpace.TemporaryState[k]];
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	{
	  TmpMonomial2 = TmpDestinationHilbertSpaceMonomial[j];
	  int TmpIndex2 = 0;
	  int TmpIndex3 = 0;
	  int TmpIndex4 = 0;
	  while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpIndex3 < nbrBosonSector)) 
	    {
	      while ((TmpIndex2 < ComplementaryNbrBosonSector) && (TmpMonomial2[TmpIndex3] <= TmpMonomial1[TmpIndex2]))
		{
		  TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
		  ++TmpIndex2;
		  ++TmpIndex4;		  
		}
	      if (TmpIndex2 < ComplementaryNbrBosonSector)
		{
		  while ((TmpIndex3 < nbrBosonSector) && (TmpMonomial1[TmpIndex2] <= TmpMonomial2[TmpIndex3]))
		    {
		      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
		      ++TmpIndex3;
		      ++TmpIndex4;		  
		    }
		}
	    }
	  while (TmpIndex2 < ComplementaryNbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial1[TmpIndex2];
	      ++TmpIndex2;
	      ++TmpIndex4;		  
	    }
	  while (TmpIndex3 < nbrBosonSector)
	    {
	      TmpMonomial3[TmpIndex4] = TmpMonomial2[TmpIndex3];
	      ++TmpIndex3;
	      ++TmpIndex4;		  
	    }

	  unsigned long TmpState = this->ConvertFromMonomial(TmpMonomial3);
	  int TmpPos = this->FermionBasis->FindStateIndex(TmpState,  TmpMonomial3[0] + this->NbrBosons - 1);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      this->FermionToBoson(TmpState, TmpMonomial3[0] + this->NbrBosons - 1, this->TemporaryState, this->TemporaryStateLzMax);
	      double TmpFactorial = 0.0;	      
	      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
		TmpFactorial += LogFactorials[this->TemporaryState[k]];
	      TmpFactorial -= TmpHilbertSpaceFactorial + TmpDestinationLogFactorials[j] + TmpLogBinomial;
	      TmpFactorial *= 0.5; 	      
	      ++TmpNbrNonZeroElements;
	      double Tmp = exp(TmpFactorial) * groundState[TmpPos];
	      TmpEntanglementMatrix.SetMatrixElement(j, MinIndex, Tmp);
	    }
	}
    }
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      delete[] TmpDestinationHilbertSpaceOccupationNumbers[i];
      delete[] TmpDestinationHilbertSpaceMonomial[i];
    }
  delete[] TmpDestinationHilbertSpaceOccupationNumbers;
  delete[] TmpDestinationHilbertSpaceLzMax;
  delete[] TmpDestinationHilbertSpaceMonomial;
  delete[] TmpMonomial1;
  delete[] TmpMonomial3;
  if (TmpNbrNonZeroElements > 0)	
    return TmpEntanglementMatrix;
  else
    {
      RealMatrix TmpEntanglementMatrixZero;
      return TmpEntanglementMatrixZero;
    }
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using real space partition. The density matrix is only evaluated in a given Lz sector.
// and computed from precalculated entanglement matrix
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees
// thetaTop =  inclination angle defining one edge of the cut in degrees
// thetaBottom = inclination angle defining the bottom edge of the cut. thetaBottom>thetaTop in degrees
// entanglementMatrix = reference on the entanglement matrix
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix BosonOnSphereShort::EvaluatePartialDensityMatrixRealSpacePartitionFromEntanglementMatrix (int nbrBosonSector, int lzSector, double thetaTop, double thetaBottom, double phiRange, RealMatrix& entanglementMatrix)
{
  if ((thetaBottom <= thetaTop) || (phiRange <= 0.0))
    {
      RealSymmetricMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }

  thetaTop *= M_PI / 180.0;
  thetaBottom *= M_PI / 180.0;
  phiRange /= 360.0;
  
  double* IncompleteBetaThetaTop = 0;
  double* IncompleteBetaThetaBottom = 0;

  this->EvaluatePartialDensityMatrixRealSpacePartitionCoefficient(this->LzMax, thetaTop, thetaBottom, IncompleteBetaThetaTop, IncompleteBetaThetaBottom);

  int ComplementaryNbrBosonSector = this->NbrBosons - nbrBosonSector;
  BosonOnSphereShort TmpDestinationHilbertSpace(nbrBosonSector, lzSector, this->LzMax);
  cout << "subsystem Hilbert space dimension = " << TmpDestinationHilbertSpace.HilbertSpaceDimension << endl;
  unsigned long* TmpMonomial1 = new unsigned long [ComplementaryNbrBosonSector];
  unsigned long* TmpMonomial3 = new unsigned long [this->NbrBosons];

  BosonOnSphereShort TmpHilbertSpace(ComplementaryNbrBosonSector, this->TotalLz - lzSector, this->LzMax);
  for (int i = 0; i < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++i)
    {
      TmpDestinationHilbertSpace.ConvertToMonomial(TmpDestinationHilbertSpace.FermionBasis->StateDescription[i], TmpDestinationHilbertSpace.FermionBasis->StateLzMax[i], TmpMonomial3);
      double Tmp = 0.0;
      Tmp = 0.5 * nbrBosonSector * log(phiRange);      
      for( int j = 0; j < nbrBosonSector; j++)
	{
	  Tmp += log( IncompleteBetaThetaBottom[TmpMonomial3[j]] - IncompleteBetaThetaTop[TmpMonomial3[j]]);
	}
      Tmp = exp(0.5 * Tmp);
      for (int j = 0; j < TmpHilbertSpace.HilbertSpaceDimension; ++j)          
	entanglementMatrix(i, j) *= Tmp;      
    }
  
  for (int MinIndex = 0; MinIndex < TmpHilbertSpace.HilbertSpaceDimension; ++MinIndex)    
    {
      TmpHilbertSpace.ConvertToMonomial(TmpHilbertSpace.FermionBasis->StateDescription[MinIndex], TmpHilbertSpace.FermionBasis->StateLzMax[MinIndex], TmpMonomial1);
      double FormFactor = 0.0;
      for (int i=0; i < ComplementaryNbrBosonSector; i++)
	FormFactor += log(1.0 - IncompleteBetaThetaBottom[TmpMonomial1[i]] + IncompleteBetaThetaTop[TmpMonomial1[i]] + (1.0 - phiRange) * (IncompleteBetaThetaBottom[TmpMonomial1[i]] - IncompleteBetaThetaTop[TmpMonomial1[i]]) );
      FormFactor = exp(0.5 * FormFactor);
      for (int j = 0; j < TmpDestinationHilbertSpace.HilbertSpaceDimension; ++j)
	entanglementMatrix(j, MinIndex) *= FormFactor; 
    }
     
  delete[] TmpMonomial1;
  delete[] TmpMonomial3;

  RealSymmetricMatrix TmpDensityMatrix(entanglementMatrix, true);
  return TmpDensityMatrix;
}

// convert a state such that its components are now expressed in the unnormalized basis
//
// state = reference to the state to convert
// reference = set which component as to be normalized to 1
// symmetryFactor = if true also remove the symmetry factors
// return value = converted state

RealVector& BosonOnSphereShort::ConvertToUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = 1.0;
  if (reference >= 0l)
    {
      Factor /= state[reference];
    }
  else
    {
      reference = 0l;
    }
  this->ConvertToMonomial(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], TmpMonomialReference);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      SqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      InvSqrtCoefficients[k] = 1.0 / SqrtCoefficients[k];
    }
  FactorialCoefficient ReferenceFactorial;
  FactorialCoefficient Factorial;
  this->FermionToBoson(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], 
		       this->TemporaryState, this->TemporaryStateLzMax);
  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
    if (this->TemporaryState[k] > 1)
      ReferenceFactorial.FactorialDivide(this->TemporaryState[k]);
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;
      while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	{
	  while ((Index1 < this->NbrBosons) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
	    {
	      ++Index1;
	      ++Index2;
	    }
	  while ((Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }	  
	}
      while (Index1 < this->NbrBosons)
	{
	  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	  ++Index1;
	}
      while (Index2 < this->NbrBosons)
	{
	  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	  ++Index2;
	}
      if (symmetryFactor == true)
	{
	  Factorial = ReferenceFactorial;
	  this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialMultiply(this->TemporaryState[k]);
	  Coefficient *= sqrt(Factorial.GetNumericalValue());
	}
      state[i] *= Coefficient;
    }
  return state;
}

// convert a state such that its components are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// symmetryFactor = if true also add the symmetry factors
// return value = converted state

RealVector& BosonOnSphereShort::ConvertFromUnnormalizedMonomial(RealVector& state, long reference, bool symmetryFactor)
{
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = 1.0;
  if (reference >= 0l)
    Factor = 1.0;
  else
    reference = 0l;
  this->ConvertToMonomial(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], TmpMonomialReference);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      SqrtCoefficients[k] = 1.0 / InvSqrtCoefficients[k];
    }
  FactorialCoefficient ReferenceFactorial;
  FactorialCoefficient Factorial;
  this->FermionToBoson(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], 
		       this->TemporaryState, this->TemporaryStateLzMax);
  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
    if (this->TemporaryState[k] > 1)
      ReferenceFactorial.FactorialMultiply(this->TemporaryState[k]);
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;
      if (symmetryFactor == true)
	{
	  while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	    {
	      while ((Index1 < this->NbrBosons) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
		{
		  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
		  ++Index1;
		}
	      while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
		{
		  ++Index1;
		  ++Index2;
		}
	      while ((Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
		{
		  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
		  ++Index2;
		}	  
	    }
	  while (Index1 < this->NbrBosons)
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while (Index2 < this->NbrBosons)
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }
	  Factorial = ReferenceFactorial;
	  this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialDivide(this->TemporaryState[k]);
	  Coefficient *= sqrt(Factorial.GetNumericalValue());
	}
      else
	{
	  Factorial = ReferenceFactorial;
	  this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			       this->TemporaryState, this->TemporaryStateLzMax);
	  for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	    if (this->TemporaryState[k] > 1)
	      Factorial.FactorialDivide(this->TemporaryState[k]);
	  Coefficient *= sqrt(Factorial.GetNumericalValue());
	}
      state[i] *= Coefficient;
    }
  state /= state.Norm();
  return state;
}

// convert a state such that its components, given in the conformal limit,  are now expressed in the normalized basis
//
// state = reference to the state to convert
// reference = set which component has been normalized to 1
// return value = converted state

RealVector& BosonOnSphereShort::ConvertFromConformalLimit(RealVector& state, long reference)
{
  unsigned long* TmpMonomialReference = new unsigned long [this->NbrBosons];
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  double Factor = 1.0;
  if (reference >= 0l)
    Factor = 1.0;
  else
    reference = 0l;
  this->ConvertToMonomial(this->FermionBasis->StateDescription[reference], this->FermionBasis->StateLzMax[reference], TmpMonomialReference);
  double* SqrtCoefficients = new double [this->LzMax + 1];
  double* InvSqrtCoefficients = new double [this->LzMax + 1];
  BinomialCoefficients Binomials(this->LzMax);
  for (int k = 0; k <= this->LzMax; ++k)
    {
      InvSqrtCoefficients[k] = sqrt(Binomials.GetNumericalCoefficient(this->LzMax, k));
      SqrtCoefficients[k] = 1.0 / InvSqrtCoefficients[k];
    }
  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
    {
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
      int Index1 = 0;
      int Index2 = 0;
      double Coefficient = Factor;

      while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons))
	{
	  while ((Index1 < this->NbrBosons) && (TmpMonomialReference[Index1] > TmpMonomial[Index2]))
	    {
	      Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	      ++Index1;
	    }
	  while ((Index1 < this->NbrBosons) && (Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] == TmpMonomial[Index2]))
	    {
	      ++Index1;
	      ++Index2;
	    }
	  while ((Index2 < this->NbrBosons) && (TmpMonomialReference[Index1] < TmpMonomial[Index2]))
	    {
	      Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	      ++Index2;
	    }	  
	}
      while (Index1 < this->NbrBosons)
	{
	  Coefficient *= InvSqrtCoefficients[TmpMonomialReference[Index1]];
	  ++Index1;
	}
      while (Index2 < this->NbrBosons)
	{
	  Coefficient *= SqrtCoefficients[TmpMonomial[Index2]];
	  ++Index2;
	}
      state[i] *= Coefficient;
    }
  state /= state.Norm();
  return state;
}

// fuse two states which belong to different Hilbert spaces 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// leftVector = reference on the vector whose Hilbert space will be fuse to the left
// rightVector = reference on the vector whose Hilbert space will be fuse to the right
// padding = number of unoccupied one body states that have to be inserted between the fused left and right spaces
// leftSpace = point to the Hilbert space that will be fuse to the left
// rightSpace = point to the Hilbert space that will be fuse to the right
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// coefficient = optional multiplicative factor to apply to the fused state 
// return value = reference on the fused state

RealVector& BosonOnSphereShort::FuseStates (RealVector& outputVector, RealVector& leftVector, RealVector& rightVector, int padding, 
					   ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace,
					    bool symmetrizedFlag, double coefficient)
{
  BosonOnSphereShort* LeftSpace = (BosonOnSphereShort*) leftSpace;
  BosonOnSphereShort* RightSpace = (BosonOnSphereShort*) rightSpace;
   int StateShift = RightSpace->FermionBasis->LzMax + padding + 2;
  for (long i = 0; i <  LeftSpace->LargeHilbertSpaceDimension; ++i)
    {
      unsigned long TmpState1 = LeftSpace->FermionBasis->StateDescription[i] << StateShift;
      double Coefficient = coefficient * leftVector[i];
      int TmpLzMax = this->FermionBasis->LzMax;
      while ((TmpState1 >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      if (symmetrizedFlag == false)
	{
	  for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState2 = RightSpace->FermionBasis->StateDescription[j];
	      TmpState2 |= TmpState1;
	      double Coefficient2 = Coefficient;
	      Coefficient2 *= rightVector[j];	  
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
	      outputVector[TmpIndex] = Coefficient2;
	    }
	}
      else
	{
	  for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState2 = RightSpace->FermionBasis->StateDescription[j];
	      TmpState2 |= TmpState1;
	      double Coefficient2 = Coefficient;
	      Coefficient2 *= rightVector[j];	  
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
	      outputVector[TmpIndex] = Coefficient2;
	      unsigned long TmpState3 = this->FermionBasis->GetSymmetricState(TmpState2);
	      if (TmpState3 != TmpState2)
		{
		  int TmpLzMax2 = this->FermionBasis->LzMax;
		  while ((TmpState3 >> TmpLzMax2) == 0x0ul)
		    --TmpLzMax2;
		  TmpIndex = this->FermionBasis->FindStateIndex(TmpState3, TmpLzMax2);
		  outputVector[TmpIndex] = Coefficient2;      
		}
	    }
	}
    }
  return outputVector;
}

// fuse two states which belong to different Hilbert spaces 
//
// outputVector = reference on the vector which will contain the fused states (without zeroing components which do not occur in the fusion)
// leftVector = reference on the vector whose Hilbert space will be fuse to the left
// rightVector = reference on the vector whose Hilbert space will be fuse to the right
// padding = number of unoccupied one body states that have to be inserted between the fused left and right spaces
// leftSpace = point to the Hilbert space that will be fuse to the left
// rightSpace = point to the Hilbert space that will be fuse to the right
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// coefficient = optional multiplicative factor to apply to the fused state 
// return value = reference on the fused state

LongRationalVector& BosonOnSphereShort::FuseStates (LongRationalVector& outputVector, LongRationalVector& leftVector, LongRationalVector& rightVector, int padding, 
						    ParticleOnSphere* leftSpace, ParticleOnSphere* rightSpace,
						    bool symmetrizedFlag, LongRational& coefficient)
{
  BosonOnSphereShort* LeftSpace = (BosonOnSphereShort*) leftSpace;
  BosonOnSphereShort* RightSpace = (BosonOnSphereShort*) rightSpace;
   int StateShift = RightSpace->FermionBasis->LzMax + padding + 2;
  for (long i = 0; i <  LeftSpace->LargeHilbertSpaceDimension; ++i)
    {
      unsigned long TmpState1 = LeftSpace->FermionBasis->StateDescription[i] << StateShift;
      LongRational Coefficient = coefficient * leftVector[i];
      int TmpLzMax = this->FermionBasis->LzMax;
      while ((TmpState1 >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      if (symmetrizedFlag == false)
	{
	  for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState2 = RightSpace->FermionBasis->StateDescription[j];
	      TmpState2 |= TmpState1;
	      LongRational Coefficient2 = Coefficient;
	      Coefficient2 *= rightVector[j];	  
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
	      outputVector[TmpIndex] = Coefficient2;
	    }
	}
      else
	{
	  for (long j = 0; j < RightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      unsigned long TmpState2 = RightSpace->FermionBasis->StateDescription[j];
	      TmpState2 |= TmpState1;
	      LongRational Coefficient2 = Coefficient;
	      Coefficient2 *= rightVector[j];	  
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax);
	      outputVector[TmpIndex] = Coefficient2;
	      unsigned long TmpState3 = this->FermionBasis->GetSymmetricState(TmpState2);
	      if (TmpState3 != TmpState2)
		{
		  int TmpLzMax2 = this->FermionBasis->LzMax;
		  while ((TmpState3 >> TmpLzMax2) == 0x0ul)
		    --TmpLzMax2;
		  TmpIndex = this->FermionBasis->FindStateIndex(TmpState3, TmpLzMax2);
		  outputVector[TmpIndex] = Coefficient2;      
		}
	    }
	}
    }
  return outputVector;
}

// use product rule to produce part of the components of a system from a smaller one
//
// outputVector = reference on the vector which will contain the product rule state  (without zeroing components which do not occur in the fusion)
// inputVector = reference on the vector associated to the smaller system
// inputSpace = pointer to the Hilbert space of the smaller system
// commonPattern = array describing the shared leftmost pattern between the n-body states in both the smaller and larger system sizes
// commonPatterSize = number of elements in the commonPattern array
// addedPattern = array describing the pattern that has to be inserted to go from the smaller system to the larger one
// addedPatterSize = number of elements in the addedPattern array
// coefficient = multiplicqtive fqctor to go fron the component of the smaller system to the larger one
// symmetrizedFlag = assume that the target state has to be invariant under the Lz<->-Lz symmetry
// return value = reference on the product rule state

RealVector& BosonOnSphereShort::ProductRules (RealVector& outputVector, RealVector& inputVector, ParticleOnSphere* inputSpace, 
					      int* commonPattern, int commonPatterSize, int* addedPattern, int addedPatterSize,
					      double coefficient, bool symmetrizedFlag)
{
  BosonOnSphereShort* InputSpace = (BosonOnSphereShort*) inputSpace;
  int NbrParticlesCommonPattern = 0;
  unsigned long InputPattern = 0x0ul;
  int TmpIndex = 0;
  for (int i = commonPatterSize - 1; i >= 0; --i)
    {
      for (int j = 0; j < commonPattern[i]; ++j)
	{
	  InputPattern |= 0x1ul << TmpIndex;
	  ++TmpIndex;
	}
      ++TmpIndex;
      NbrParticlesCommonPattern += commonPattern[i];
    }
  int NbrParticlesAddedPattern = 0;
  unsigned long OutputPattern = 0x0ul;
  TmpIndex = 0;
  for (int i = addedPatterSize - 1; i >= 0; --i)
    {
      for (int j = 0; j < addedPattern[i]; ++j)
	{
	  OutputPattern |= 0x1ul << TmpIndex;
	  ++TmpIndex;
	}
      ++TmpIndex;
      NbrParticlesAddedPattern += addedPattern[i];
    }
  unsigned long InputMask = ((0x1ul << (commonPatterSize + NbrParticlesCommonPattern)) - 1ul) << (InputSpace->FermionBasis->LzMax - (commonPatterSize + NbrParticlesCommonPattern) + 1);
  unsigned long InputMask2 = ~InputMask;  
  OutputPattern |= InputPattern << (addedPatterSize + NbrParticlesAddedPattern);
  OutputPattern <<= (InputSpace->FermionBasis->LzMax - (commonPatterSize + NbrParticlesCommonPattern) + 2);
  InputPattern <<= (InputSpace->FermionBasis->LzMax - (commonPatterSize + NbrParticlesCommonPattern) + 2); 
  cout << hex << InputMask << " " << OutputPattern << " " << InputPattern << dec << endl;
  int OutputLzMax = this->FermionBasis->LzMax;
  while (((OutputPattern >> OutputLzMax) & 0x1ul) == 0x0ul)
    --OutputLzMax;
  long Count = 0l;
  if (symmetrizedFlag == false)
    {
      for (int i = 0; i <  InputSpace->HilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState1 = InputSpace->FermionBasis->StateDescription[i];
	  //	  cout << i << " " << hex << TmpState1 << dec << endl;
	  if ((TmpState1 & InputMask) == InputPattern)
	    {
	      TmpState1 &= InputMask2;
	      TmpState1 |= OutputPattern;
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState1, OutputLzMax);
	      double& TmpCoef = outputVector[TmpIndex];
	      if (TmpCoef == 0.0)
		++Count;
	      TmpCoef = coefficient * inputVector[i];	  
	    }
	}
    }
  else
    {
      for (int i = 0; i <  InputSpace->HilbertSpaceDimension; ++i)
	{
	  unsigned long TmpState1 = InputSpace->FermionBasis->StateDescription[i];
	  if ((TmpState1 & InputMask) == InputPattern)
	    {
	      TmpState1 &= InputMask2;
	      TmpState1 |= OutputPattern;
	      int TmpIndex = this->FermionBasis->FindStateIndex(TmpState1, OutputLzMax);
	      double& TmpCoef = outputVector[TmpIndex];
	      if (TmpCoef == 0.0)
		++Count;
	      double TmpCoef3 = coefficient * inputVector[i];
	      TmpCoef = TmpCoef3;	  
	      unsigned long TmpState2 = this->FermionBasis->GetSymmetricState(TmpState1);
	      if (TmpState2 != TmpState1)
		{
		  int TmpLzMax2 = this->FermionBasis->LzMax;
		  while ((TmpState2 >> TmpLzMax2) == 0x0ul)
		    --TmpLzMax2;
		  TmpIndex = this->FermionBasis->FindStateIndex(TmpState2, TmpLzMax2);
		  double& TmpCoef2 = outputVector[TmpIndex];
                  if (TmpCoef2 == 0.0)
                    ++Count;
                  TmpCoef2 = TmpCoef3;
		}
	    }

	}
    }
  cout << "nbr of newly added components : " << Count << endl;
  return outputVector;
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state = reference on the unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

double BosonOnSphereShort::JackSqrNormalization (RealVector& outputVector, long minIndex, long nbrComponents)
{
  double SqrNorm = 0.0;
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = ((unsigned long) this->LzMax) >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  for (long i = minIndex; i < MaxIndex; ++i)
    {
      Factorial.SetToOne();
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
      this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			   this->TemporaryState, this->TemporaryStateLzMax);
      for (int k = 0; k < this->NbrBosons; ++k)
	{
	  if (HalfLzMax < TmpMonomial[k])
	    Factorial.PartialFactorialMultiply(HalfLzMax + 1, TmpMonomial[k]);
	  else
	    if (HalfLzMax > TmpMonomial[k])
	      Factorial.PartialFactorialDivide(TmpMonomial[k] + 1, HalfLzMax);
	}	      
      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	if (this->TemporaryState[k] > 1)
	  Factorial.FactorialDivide(this->TemporaryState[k]);
      SqrNorm +=(outputVector[i] * outputVector[i]) * Factorial.GetNumericalValue();
      if ((i & 0x3fffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  cout << endl;
  delete[] TmpMonomial;
  return SqrNorm;
}

// compute part of the Jack polynomial square normalization in a given range of indices
//
// state = reference on the unnormalized Jack polynomial
// minIndex = first index to compute 
// nbrComponents = number of indices to compute (0 if they all have to be computed from minIndex)
// return value = quare normalization 

LongRational BosonOnSphereShort::JackSqrNormalization (LongRationalVector& outputVector, long minIndex, long nbrComponents)
{
  LongRational SqrNorm = 0l;
  unsigned long* TmpMonomial = new unsigned long [this->NbrBosons];
  FactorialCoefficient Factorial;
  unsigned long HalfLzMax = ((unsigned long) this->LzMax) >> 1;
  long MaxIndex = minIndex + nbrComponents;
  if (MaxIndex == minIndex)
    MaxIndex = this->LargeHilbertSpaceDimension;
  for (long i = minIndex; i < MaxIndex; ++i)
    {
      Factorial.SetToOne();
      this->ConvertToMonomial(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], TmpMonomial);
      this->FermionToBoson(this->FermionBasis->StateDescription[i], this->FermionBasis->StateLzMax[i], 
			   this->TemporaryState, this->TemporaryStateLzMax);
      for (int k = 0; k < this->NbrBosons; ++k)
	{
	  if (HalfLzMax < TmpMonomial[k])
	    Factorial.PartialFactorialMultiply(HalfLzMax + 1, TmpMonomial[k]);
	  else
	    if (HalfLzMax > TmpMonomial[k])
	      Factorial.PartialFactorialDivide(TmpMonomial[k] + 1, HalfLzMax);
	}	      
      for (int k = 0; k <= this->TemporaryStateLzMax; ++k)
	if (this->TemporaryState[k] > 1)
	  Factorial.FactorialDivide(this->TemporaryState[k]);
      SqrNorm += (outputVector[i] * outputVector[i]) * Factorial.GetLongRationalValue();
      if ((i & 0x3fffl) == 0l)
	{
	  cout << i << " / " << this->LargeHilbertSpaceDimension << " (" << ((i * 100l) / this->LargeHilbertSpaceDimension) << "%)           \r";
	  cout.flush();
	}
    }
  cout << endl;
  delete[] TmpMonomial;
  return SqrNorm;
}

// get Lz component of a component
//
// j = index of the component in Hilbert space
// return value = twice the  component
int BosonOnSphereShort::GetLzValue(int j)
{
  return this->TotalLz;
}

// compute all Kostka coefficients for a given Schur polynomial 
//
// index = index of the partition that describe the Schur polynomial 
// kostkaVector = vector where kostka numbers will be stored
// return value = number of kostka coefficients

long BosonOnSphereShort::KostkaForOneSchur(long index, RealVector& kostkaVector)
{	
  unsigned long * InitialState = new unsigned long[this->NbrLzValue];
  int finalStateLzMax;
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], InitialState, finalStateLzMax);
  int* ReferenceState = new int[this->NbrLzValue];
  for(int i = 0; i <= finalStateLzMax; i++)
    {
      ReferenceState[i] = (int) InitialState[i];
    }
  for(int i = finalStateLzMax + 1; i < this->NbrLzValue; i++)
    {
      ReferenceState[i] = 0;
    }
  BosonOnSphereHaldaneBasisShort* SqueezedSpace = new BosonOnSphereHaldaneBasisShort(this->NbrBosons, this->TotalLz, this->LzMax, 
										     ReferenceState);
  RealVector JackVector(SqueezedSpace->GetHilbertSpaceDimension(), true);
  SqueezedSpace->GenerateJackPolynomial(JackVector, 1.0);
  int Dimension = SqueezedSpace->GetHilbertSpaceDimension();
  
  for (int i = 0; i < Dimension; i++)
    {
      int TmpLzMax = this->LzMax + this->NbrBosons - 1;
      unsigned long TmpState = SqueezedSpace->FermionBasis->StateDescription[i];
      while ((TmpState >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      kostkaVector[this->FermionBasis->FindStateIndex(TmpState, TmpLzMax)] = JackVector[i];
    }
  delete SqueezedSpace;
  delete[] ReferenceState;
  return Dimension;
}

// divide a set of fermionic vectors by a Jastrow factor to get their bosonic counterpart
//
// sourceVector = array of fermionic statesc be divided by a jastrow factor
// outputVector = array of where bosonic states will be stored
// firstComponent = index of the first component to transform 
// nbrComponent = number of components to compute
// nbrStates  = number of states to handle

void BosonOnSphereShort::DivideByJastrow(RealVector* sourceVector, RealVector* OutputVector, long firstComponent,
					 long nbrComponent, int nbrStates)
{
  long IndiceMax = firstComponent + nbrComponent;
  RealVector KostkaVector (this->GetHilbertSpaceDimension(), true);
  for(long i = firstComponent; i < IndiceMax; i++)
    {
      this->KostkaForOneSchur(i,KostkaVector);
      for(long j = i; j < this->GetHilbertSpaceDimension(); j++)
	{
	  for(int k = 0; k < nbrStates; k++)	
	    OutputVector[k][j] += sourceVector[k][i] * KostkaVector[j];
	}
      KostkaVector.ClearVector();
    }
}

void BosonOnSphereShort::FuseParticlesInState(RealVector& firstState, RealVector& outputVector, BosonOnSphereShort* finalSpace, 
					      long minIndex, long nbrComponents)
{
  long MaxIndex= minIndex + nbrComponents;  
  if (nbrComponents == 0l)
    MaxIndex = finalSpace->GetLargeHilbertSpaceDimension();
  long* WeigthVector = new long[finalSpace->GetLargeHilbertSpaceDimension()];
  unsigned long* IndicesVector = new unsigned long[finalSpace->GetLargeHilbertSpaceDimension()];
  for(long Index = minIndex; Index < MaxIndex; Index++)
    {
      int Size = this->FuseParticlesInMonomial(Index, finalSpace, WeigthVector, IndicesVector);
      for(int i = 0; i < Size; i++)
	{
	  int TmpLzMax = finalSpace->LzMax + finalSpace->NbrBosons - 1;
	  while ((IndicesVector[i] >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  outputVector[finalSpace->FermionBasis->FindStateIndex(IndicesVector[i], TmpLzMax)] += firstState[Index] * WeigthVector[i];
	}
    }
}

// fuse particles two by two in a given monomial
//
// index = monomial index
// finalSpace = space where the fused state lies
// weigthVector = weigths of each fused component
// indicesVector = indices of each fused component
// return value = number of generated components when fusing

int BosonOnSphereShort::FuseParticlesInMonomial(long index, BosonOnSphereShort* finalSpace, long* weigthVector, unsigned long* indicesVector)
{
  FactorialCoefficient Coefficient;
  Coefficient.SetToOne();
  Coefficient.Power2Multiply(finalSpace->NbrBosons);
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  unsigned long* TmpMonomial = new unsigned long[this->NbrBosons];
  unsigned long* FinalMonomial = new unsigned long[finalSpace->NbrBosons];
  this->GetMonomial(index,TmpMonomial);
  int Size = 0;
  Size = this->GeneratePairs(TmpMonomial, weigthVector, indicesVector, FinalMonomial, this->NbrBosons, 0, finalSpace);
  for(int l = 0; l <= TemporaryStateLzMax; l++)
    {
      Coefficient.FactorialDivide(this->TemporaryState[l]);		
    }
  FactorialCoefficient Coefficient1;
  for(int i = 0; i < Size; i++)
    {
      Coefficient1 = Coefficient;
      int TmpLzMax = finalSpace->LzMax + finalSpace->NbrBosons - 1;
      while ((indicesVector[i] >> TmpLzMax) == 0x0ul)
	--TmpLzMax;
      finalSpace->FermionToBoson(indicesVector[i], TmpLzMax, finalSpace->TemporaryState, finalSpace->TemporaryStateLzMax);
      for(int p = 0; p <= finalSpace->TemporaryStateLzMax; p++)
	{
	  Coefficient1.FactorialMultiply(finalSpace->TemporaryState[p]);
	}
      if(weigthVector[i] < 0)
	{
	  Coefficient1 *= (-weigthVector[i]);
	  weigthVector[i] = (-Coefficient1.GetIntegerValue());
	}
      else
	{
	  Coefficient1 *= weigthVector[i];
	  weigthVector[i] = Coefficient1.GetIntegerValue();
	}
    }
  return Size;
}
	
int BosonOnSphereShort::GeneratePairs(unsigned long* Monomial, long* weigthVector, unsigned long* indicesVector, unsigned long* FinalMonomial, int remainder, int count, BosonOnSphereShort* finalSpace)
{
  if(remainder == 0)
    {
      unsigned long* TMonomial = new unsigned long[finalSpace->NbrBosons];
      for(int pll = 0; pll < finalSpace->NbrBosons; pll++)
	TMonomial[pll] = FinalMonomial[pll];
      SortArrayDownOrdering(TMonomial,finalSpace->NbrBosons);
      unsigned long Result = finalSpace->ConvertFromMonomial(TMonomial);
      for(int k = 0; k < count; k++)
	{
	  if(Result == indicesVector[k])
	    {
	      weigthVector[k]++;
	      return 0;
	    }
	}
      indicesVector[count] = Result;
      weigthVector[count] = 1l;
      return 1;
    }
  for(int i = 1; i < remainder; i++)
    {
      FinalMonomial[finalSpace->NbrBosons - (remainder / 2)] = Monomial[0] + Monomial[i];
      unsigned long* TmpMonomial = new unsigned long[remainder - 2];
      for(int j = 0; j < (i - 1); j++)
	{
	  TmpMonomial[j] = Monomial[j+1];
	}
      for(int j = i; j < (remainder - 1); j++)
	{
	  TmpMonomial[j - 1] = Monomial[j + 1];
	}
      if(remainder == 2)
	count += this->GeneratePairs(TmpMonomial, weigthVector, indicesVector, FinalMonomial, remainder - 2, count, finalSpace);
      else
	count = this->GeneratePairs(TmpMonomial, weigthVector, indicesVector, FinalMonomial, remainder - 2, count, finalSpace);
    }
  return count;
}


// symmetrized a product of two uncoupled states 
//
// outputVector = reference on the vector which will contain the symmetrozed state
// leftVector = reference on the vector associated to the first color
// rightVector = reference on the vector associated to the second color
// leftSpace = pointer to the Hilbert space of the first color
// rightSpace = pointer to the Hilbert space of the second color
// unnormalizedBasisFlag = assume evrything has to be done in the unnormalized basis
// return value = symmetrized state

RealVector BosonOnSphereShort::SymmetrizeU1U1State (RealVector& leftVector, RealVector& rightVector, 
						    BosonOnSphereShort* leftSpace, BosonOnSphereShort* rightSpace, bool unnormalizedBasisFlag)
{
  RealVector SymmetrizedVector (this->LargeHilbertSpaceDimension);
  FactorialCoefficient Factorial1;
  FactorialCoefficient Factorial2;
  if (unnormalizedBasisFlag == true)
    {
      for (long i = 0l; i < leftSpace->LargeHilbertSpaceDimension; ++i)
	{
	  this->FermionToBoson(leftSpace->FermionBasis->StateDescription[i], leftSpace->FermionBasis->StateLzMax[i], 
			       leftSpace->TemporaryState, leftSpace->TemporaryStateLzMax);
	  for (int k = leftSpace->TemporaryStateLzMax + 1;  k <= leftSpace->LzMax; ++k)
	    leftSpace->TemporaryState[k] = 0;
	  double TmpCoefficient = leftVector[i];
	  Factorial1.SetToOne();
	  Factorial1.Power2Divide(leftSpace->NbrBosons);
	  for (int k = 0; k <= leftSpace->TemporaryStateLzMax; ++k)
	    if (leftSpace->TemporaryState[k] > 1)
	      Factorial1.FactorialDivide(leftSpace->TemporaryState[k]);
	  
	  for (long j = 0l; j < rightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      this->FermionToBoson(rightSpace->FermionBasis->StateDescription[j], rightSpace->FermionBasis->StateLzMax[j], 
				   rightSpace->TemporaryState, rightSpace->TemporaryStateLzMax);
	      int k = 0;
	      for (; k <= rightSpace->TemporaryStateLzMax; ++k)
		this->TemporaryState[k] = leftSpace->TemporaryState[k] + rightSpace->TemporaryState[k];
	      this->TemporaryStateLzMax = rightSpace->TemporaryStateLzMax;
	      if (leftSpace->TemporaryStateLzMax > rightSpace->TemporaryStateLzMax)
		{
		  for (; k <= leftSpace->TemporaryStateLzMax; ++k)
		    this->TemporaryState[k] = leftSpace->TemporaryState[k];
		  this->TemporaryStateLzMax = leftSpace->TemporaryStateLzMax;
		}
	      int TmpPos = this->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
	      if (TmpPos < this->HilbertSpaceDimension)
		{
		  Factorial2 = Factorial1;
		  for (k = 0; k <= rightSpace->TemporaryStateLzMax; ++k)
		    if (rightSpace->TemporaryState[k] > 1)
		      Factorial2.FactorialDivide(rightSpace->TemporaryState[k]);
		  for (k = 0; k <= this->TemporaryStateLzMax; ++k)
		    if (this->TemporaryState[k] > 1)
		      Factorial2.FactorialMultiply(this->TemporaryState[k]);	      
		  SymmetrizedVector[TmpPos] += Factorial2.GetNumericalValue() * TmpCoefficient * rightVector[j];
		}
	    }
	}
    }
  else
    {
      for (long i = 0l; i < leftSpace->LargeHilbertSpaceDimension; ++i)
	{
	  this->FermionToBoson(leftSpace->FermionBasis->StateDescription[i], leftSpace->FermionBasis->StateLzMax[i], 
			       leftSpace->TemporaryState, leftSpace->TemporaryStateLzMax);
	  for (int k = leftSpace->TemporaryStateLzMax + 1;  k <= leftSpace->LzMax; ++k)
	    leftSpace->TemporaryState[k] = 0;
	  double TmpCoefficient = leftVector[i];
	  Factorial1.SetToOne();
	  Factorial1.Power2Divide(leftSpace->NbrBosons);
	  for (int k = 0; k <= leftSpace->TemporaryStateLzMax; ++k)
	    if (leftSpace->TemporaryState[k] > 1)
	      Factorial1.FactorialDivide(leftSpace->TemporaryState[k]);
	  
	  for (long j = 0l; j < rightSpace->LargeHilbertSpaceDimension; ++j)
	    {
	      this->FermionToBoson(rightSpace->FermionBasis->StateDescription[j], rightSpace->FermionBasis->StateLzMax[j], 
				   rightSpace->TemporaryState, rightSpace->TemporaryStateLzMax);
	      int k = 0;
	      for (; k <= rightSpace->TemporaryStateLzMax; ++k)
		this->TemporaryState[k] = leftSpace->TemporaryState[k] + rightSpace->TemporaryState[k];
	      this->TemporaryStateLzMax = rightSpace->TemporaryStateLzMax;
	      if (leftSpace->TemporaryStateLzMax > rightSpace->TemporaryStateLzMax)
		{
		  for (; k <= leftSpace->TemporaryStateLzMax; ++k)
		    this->TemporaryState[k] = leftSpace->TemporaryState[k];
		  this->TemporaryStateLzMax = leftSpace->TemporaryStateLzMax;
		}
	      int TmpPos = this->FermionBasis->FindStateIndex(this->BosonToFermion(this->TemporaryState, this->TemporaryStateLzMax), this->TemporaryStateLzMax + this->NbrBosons - 1);
	      if (TmpPos < this->HilbertSpaceDimension)
		{
		  Factorial2 = Factorial1;
		  for (k = 0; k <= rightSpace->TemporaryStateLzMax; ++k)
		    if (rightSpace->TemporaryState[k] > 1)
		      Factorial2.FactorialDivide(rightSpace->TemporaryState[k]);
		  for (k = 0; k <= this->TemporaryStateLzMax; ++k)
		    if (this->TemporaryState[k] > 1)
		      Factorial2.FactorialMultiply(this->TemporaryState[k]);	      
		  SymmetrizedVector[TmpPos] += sqrt(Factorial2.GetNumericalValue()) * TmpCoefficient * rightVector[j];
		}
	    }
	}
      SymmetrizedVector /= SymmetrizedVector.Norm();
    }
  return SymmetrizedVector;
}

// request whether state with given index satisfies a general Pauli exclusion principle
// index = state index
// pauliK = number of particles allowed in consecutive orbitals
// pauliR = number of consecutive orbitals

bool BosonOnSphereShort::HasPauliExclusions(int index, int pauliK, int pauliR)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[index], this->FermionBasis->StateLzMax[index], this->TemporaryState, this->TemporaryStateLzMax);
  int Max = this->TemporaryStateLzMax + 2 - pauliR;
  
  for (int m = 0; m < Max; ++m)
    {
      int Count = TemporaryState[m];
      for (int Offset = 1; Offset < pauliR; ++Offset)
	Count += TemporaryState[m + Offset];
      if (Count > pauliK)
	return false;
    }
  return true;
}
