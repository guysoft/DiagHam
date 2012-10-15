////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of fermions on sphere that allow to use MPS with operator       //
//                                                                            //
//                        last modification : 09/10/2012                      //
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
#include "HilbertSpace/FermionOnSphereMPSWrapper.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Vector/RealVector.h"
#include "Matrix/RealMatrix.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h" 

#include <math.h>
#include <stdlib.h>
#include <fstream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constuctor
//

FermionOnSphereMPSWrapper::FermionOnSphereMPSWrapper()
{
}

// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// referenceState = array that describes the root configuration
// rowIndex = row index of the MPS element that has to be evaluated (-1 if the trace has to be considered instead of a single matrix element)
// columnIndex = column index of the MPS element that has to be evaluated
// bMatrices = array that gives the B matrices 
// memory = amount of memory granted for precalculations

FermionOnSphereMPSWrapper::FermionOnSphereMPSWrapper (int nbrFermions, int& totalLz, int lzMax, int* referenceState,  
						      int rowIndex, int columnIndex, SparseComplexMatrix* bMatrices, unsigned long memory)
{
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->TotalLz = 0;
  this->StateDescription = 0x0ul;
  int TmpIndex = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      this->StateDescription |= ((unsigned long) (referenceState[i] & 1)) << i;
      if (referenceState[i] == 1)
	{
	  this->TotalLz += i;
	  this->StateLzMax = i;
	}
    }
  this->TotalLz = ((this->TotalLz << 1) - (this->LzMax * this->NbrFermions));
  totalLz = this->TotalLz;
  this->LargeHilbertSpaceDimension = 1l;
  this->HilbertSpaceDimension = 1;
  this->Flag.Initialize();
  this->StateDescription = 0l;
  this->MaximumSignLookUp = 16;
  this->GenerateLookUpTable(memory);
  this->MPSRowIndex = rowIndex;
  this->MPSColumnIndex = columnIndex;

  int NbrBMatrices = 2;
  SparseComplexMatrix* SparseConjugateBMatrices = new SparseComplexMatrix[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      SparseConjugateBMatrices[i] = bMatrices[i].HermitianTranspose();
    }
  SparseComplexMatrix** SparseTensorProductBMatrices = new SparseComplexMatrix*[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      SparseTensorProductBMatrices[i] = new SparseComplexMatrix[NbrBMatrices];
      for (int j = 0; j < NbrBMatrices; ++j)
	{
	  SparseTensorProductBMatrices[i][j] = TensorProduct(bMatrices[i], SparseConjugateBMatrices[j]);
	}
    }
  delete[]SparseConjugateBMatrices;


  this->NormalizedB1B1 = new SparseComplexMatrix [this->LzMax + 1];
  this->NormalizedB0B1 = new SparseComplexMatrix [this->LzMax + 1];
  this->NormalizedB1B0 = new SparseComplexMatrix [this->LzMax + 1];
  this->NormalizedB0B0B1B1 = new SparseComplexMatrix [this->LzMax + 1];
  
  long TmpMemory = (((long) SparseTensorProductBMatrices[1][1].GetNbrRow()) * 
		    ((long) SparseTensorProductBMatrices[1][1].GetNbrColumn())) / 100l;
  cout << "Requested memory for sparse matrix multiplications = " << ((TmpMemory * (2l * sizeof(double) + sizeof(int))) >> 20) << "Mb" << endl;
  this->TmpMatrixElements = new Complex [TmpMemory];
  this->TmpColumnIndices = new int [TmpMemory];
  this->TmpElements = new Complex [SparseTensorProductBMatrices[1][1].GetNbrRow()];

  SparseComplexMatrix TmpMatrixNorm;
  double TmpBinomial = 1.0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      this->NormalizedB0B0B1B1[i] = SparseComplexMatrixLinearCombination(1.0, SparseTensorProductBMatrices[0][0], TmpBinomial, SparseTensorProductBMatrices[1][1]);
      this->NormalizedB1B1[i].Copy(SparseTensorProductBMatrices[1][1]);
      this->NormalizedB1B1[i] *= TmpBinomial;
      this->NormalizedB0B1[i].Copy(SparseTensorProductBMatrices[0][1]);
      this->NormalizedB0B1[i] *= sqrt(TmpBinomial);
      this->NormalizedB1B0[i].Copy(SparseTensorProductBMatrices[1][0]);
      this->NormalizedB1B0[i] *= sqrt(TmpBinomial);
      TmpBinomial *= (double) (i + 1);
      if (i < this->LzMax)
	TmpBinomial /= (double) (this->LzMax - i); 
      if (i > 0)     
	TmpMatrixNorm.Multiply(this->NormalizedB0B0B1B1[i], this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements);
      else
	TmpMatrixNorm.Copy(this->NormalizedB0B0B1B1[i]);
    }
  Complex Tmp;
  TmpMatrixNorm.GetMatrixElement(this->MPSRowIndex, this->MPSColumnIndex, Tmp);
  this->StateNormalization = Tmp.Re;

  for (int i = 0; i < NbrBMatrices; ++i)
    delete[] SparseTensorProductBMatrices[i];
  delete[] SparseTensorProductBMatrices;

}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereMPSWrapper::FermionOnSphereMPSWrapper(const FermionOnSphereMPSWrapper& fermions)
{
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->ShiftedTotalLz = fermions.ShiftedTotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->Flag = fermions.Flag;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->MPSRowIndex = fermions.MPSRowIndex;
  this->MPSColumnIndex = fermions.MPSColumnIndex;
  this->NormalizedB0B0B1B1 = fermions.NormalizedB0B0B1B1;
  this->NormalizedB1B1 = fermions.NormalizedB1B1;
  this->NormalizedB0B1 = fermions.NormalizedB0B1;
  this->NormalizedB1B0 = fermions.NormalizedB1B0;
  this->StateNormalization = fermions.StateNormalization;
  long TmpMemory = (((long) this->NormalizedB0B0B1B1[0].GetNbrRow()) * 
		    ((long)  this->NormalizedB0B0B1B1[0].GetNbrColumn())) / 100l;
  this->TmpMatrixElements = new Complex [TmpMemory];
  this->TmpColumnIndices = new int [TmpMemory];
  this->TmpElements = new Complex [this->NormalizedB0B0B1B1[0].GetNbrRow()];
}

// destructor
//

FermionOnSphereMPSWrapper::~FermionOnSphereMPSWrapper ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->NormalizedB0B0B1B1;
      delete[] this->NormalizedB1B1;
      delete[] this->NormalizedB0B1;
      delete[] this->NormalizedB1B0;
    }
  delete[] this->TmpMatrixElements;
  delete[] this->TmpColumnIndices;
  delete[] this->TmpElements;
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereMPSWrapper& FermionOnSphereMPSWrapper::operator = (const FermionOnSphereMPSWrapper& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      delete[] this->NormalizedB0B0B1B1;
      delete[] this->NormalizedB1B1;
      delete[] this->NormalizedB0B1;
      delete[] this->NormalizedB1B0;
    }
  delete[] this->TmpMatrixElements;
  delete[] this->TmpColumnIndices;
  delete[] this->TmpElements;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateLzMax = fermions.StateLzMax;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->MPSRowIndex = fermions.MPSRowIndex;
  this->MPSColumnIndex = fermions.MPSColumnIndex;
  this->NormalizedB0B0B1B1 = fermions.NormalizedB0B0B1B1;
  this->NormalizedB1B1 = fermions.NormalizedB1B1;
  this->NormalizedB0B1 = fermions.NormalizedB0B1;
  this->NormalizedB1B0 = fermions.NormalizedB1B0;
  this->StateNormalization = fermions.StateNormalization;
  long TmpMemory = (((long) this->NormalizedB0B0B1B1[0].GetNbrRow()) * 
		    ((long)  this->NormalizedB0B0B1B1[0].GetNbrColumn())) / 100l;
  this->TmpMatrixElements = new Complex [TmpMemory];
  this->TmpColumnIndices = new int [TmpMemory];
  this->TmpElements = new Complex [this->NormalizedB0B0B1B1[0].GetNbrRow()];
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereMPSWrapper::Clone()
{
  return new FermionOnSphereMPSWrapper(*this);
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnSphereMPSWrapper::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnSphereMPSWrapper::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnSphereMPSWrapper::ExtractSubspace (AbstractQuantumNumber& q, 
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

int FermionOnSphereMPSWrapper::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  if ((m1 == m2) || (n1 == n2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  unsigned long TmpState = this->StateDescription;
  coefficient = this->SignLookUpTable[(TmpState >> n2) & this->SignLookUpTableMask[n2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  TmpState &= ~(0x1ul << n2);
  coefficient *= this->SignLookUpTable[(TmpState >> n1) & this->SignLookUpTableMask[n1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  TmpState &= ~(0x1ul << n1);
  coefficient *= this->SignLookUpTable[(TmpState >> m2) & this->SignLookUpTableMask[m2]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 16))  & this->SignLookUpTableMask[m2 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 32)) & this->SignLookUpTableMask[m2 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m2 + 48)) & this->SignLookUpTableMask[m2 + 48]];
#endif
  TmpState |= (0x1ul << m2);
  coefficient *= this->SignLookUpTable[(TmpState >> m1) & this->SignLookUpTableMask[m1]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 16))  & this->SignLookUpTableMask[m1 + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 32)) & this->SignLookUpTableMask[m1 + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (m1 + 48)) & this->SignLookUpTableMask[m1 + 48]];
#endif
  TmpState |= (0x1ul << m1);
  SparseComplexMatrix TmpMatrix;
  if (m1 == 0)
    {
      if ((m1 == n1) || (m1 == n2))
	TmpMatrix.Copy(this->NormalizedB1B1[0]);
      else
	TmpMatrix.Copy(this->NormalizedB0B1[0]);
    }
  else
    {
      if (m2 == 0)
	{
	  if ((m2 == n1) || (m2 == n2))
	    TmpMatrix.Copy(this->NormalizedB1B1[0]);
	  else
	    TmpMatrix.Copy(this->NormalizedB0B1[0]);
	}
      else
	{
	  if ((0 == n1) || (0 == n2))
	    {
	      TmpMatrix.Copy(this->NormalizedB1B0[0]);
	    }
	  else
	    {
	      TmpMatrix.Copy(this->NormalizedB0B0B1B1[0]);
	    }
	}
    }
  for (int j = 1; j <= this->LzMax; ++j)
    {
      if (j == m1)
	{
	  if ((m1 == n1) || (m1 == n2))
	    TmpMatrix.Multiply(this->NormalizedB1B1[j], this->TmpMatrixElements, this->TmpColumnIndices,  this->TmpElements);
	  else
	    TmpMatrix.Multiply(this->NormalizedB0B1[j], this->TmpMatrixElements, this->TmpColumnIndices,  this->TmpElements);
	}
      else
	{
	  if (j == m2)
	    {
	      if ((m2 == n1) || (m2 == n2))
		TmpMatrix.Multiply(this->NormalizedB1B1[j], this->TmpMatrixElements, this->TmpColumnIndices,  this->TmpElements);
	      else
		TmpMatrix.Multiply(this->NormalizedB0B1[j], this->TmpMatrixElements, this->TmpColumnIndices,  this->TmpElements);
	    }
	  else
	    {
	      if ((j == n1) || (j == n2))
		{
		  TmpMatrix.Multiply(this->NormalizedB1B0[j], this->TmpMatrixElements, this->TmpColumnIndices,  this->TmpElements);
		}
	      else
		{
		  TmpMatrix.Multiply(this->NormalizedB0B0B1B1[j], this->TmpMatrixElements, this->TmpColumnIndices,  this->TmpElements);
		}
	    }
	}
    }
  Complex Tmp = 0.0;
  TmpMatrix.GetMatrixElement(this->MPSRowIndex, this->MPSColumnIndex, Tmp);
  coefficient *= Tmp.Re / this->StateNormalization;
  return 0;
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereMPSWrapper::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  --nbrIndices;
  for (int i = 0; i < nbrIndices; ++i)
    {
      if ((n[i] > this->StateLzMax) || ((this->StateDescription & (((unsigned long) (0x1)) << n[i])) == 0))
	{
	  coefficient = 0.0;
	  return this->HilbertSpaceDimension;
	}
      for (int j = i + 1; j <= nbrIndices; ++j)
	if ((n[i] == n[j]) || (m[i] == m[j]))
	  {
	    coefficient = 0.0;
	    return this->HilbertSpaceDimension; 	    
	  }
    }
  if (n[nbrIndices] > this->StateLzMax)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }

  int NewLzMax = this->StateLzMax;
  unsigned long TmpState = this->StateDescription;

  int Index;
  coefficient = 1.0;
  for (int i = nbrIndices; i >= 0; --i)
    {
      Index = n[i];
      coefficient *= this->SignLookUpTable[(TmpState >> Index) & this->SignLookUpTableMask[Index]];
      coefficient *= this->SignLookUpTable[(TmpState >> (Index+ 16))  & this->SignLookUpTableMask[Index+ 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 32)) & this->SignLookUpTableMask[Index + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (Index + 48)) & this->SignLookUpTableMask[Index + 48]];
#endif
      TmpState &= ~(((unsigned long) (0x1)) << Index);
      if (NewLzMax == Index)
	while ((TmpState >> NewLzMax) == 0)
	  --NewLzMax;
    }
  for (int i = nbrIndices; i >= 0; --i)
    {
      Index = m[i];
      if ((TmpState & (((unsigned long) (0x1)) << Index))!= 0)
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
      TmpState |= (((unsigned long) (0x1)) << Index);
    }
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereMPSWrapper::ProdA (int index, int* n, int nbrIndices)
{
  this->ProdALzMax = this->StateLzMax;
  this->ProdATemporaryState = this->StateDescription;
  int Index;
  double Coefficient = 1.0;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = n[i];
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
  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdALzMax = 0;
      return Coefficient;      
    }
  while (((this->ProdATemporaryState >> this->ProdALzMax) == 0) && (this->ProdALzMax > 0))
    --this->ProdALzMax;

  return Coefficient;
}

// apply a_n1 a_n2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
//
// index = index of the state on which the operator has to be applied
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// return value =  multiplicative factor 

double FermionOnSphereMPSWrapper::AA (int index, int n1, int n2)
{
  this->ProdATemporaryState = this->StateDescription;

  if (((ProdATemporaryState & (((unsigned long) (0x1)) << n1)) == 0) 
      || ((ProdATemporaryState & (((unsigned long) (0x1)) << n2)) == 0) || (n1 == n2))
    return 0.0;

  this->ProdALzMax = this->StateLzMax;

  double Coefficient = this->SignLookUpTable[(this->ProdATemporaryState >> n2) & this->SignLookUpTableMask[n2]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 16))  & this->SignLookUpTableMask[n2 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 32)) & this->SignLookUpTableMask[n2 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n2 + 48)) & this->SignLookUpTableMask[n2 + 48]];
#endif
  this->ProdATemporaryState &= ~(((unsigned long) (0x1)) << n2);
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> n1) & this->SignLookUpTableMask[n1]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 16))  & this->SignLookUpTableMask[n1 + 16]];
#ifdef  __64_BITS__
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 32)) & this->SignLookUpTableMask[n1 + 32]];
  Coefficient *= this->SignLookUpTable[(this->ProdATemporaryState >> (n1 + 48)) & this->SignLookUpTableMask[n1 + 48]];
#endif
  this->ProdATemporaryState &= ~(((unsigned long) (0x1)) << n1);

  if (this->ProdATemporaryState == 0x0ul)
    {
      this->ProdALzMax = 0;
      return Coefficient;      
    }
  while ((this->ProdATemporaryState >> this->ProdALzMax) == 0)
    --this->ProdALzMax;
  return Coefficient;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereMPSWrapper::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  coefficient = 1.0;
  unsigned long TmpState = this->ProdATemporaryState;
  int NewLzMax = this->ProdALzMax;
  int Index;
  for (int i = nbrIndices - 1; i >= 0; --i)
    {
      Index = m[i];
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

// apply a^+_m1 a^+_m2 operator to the state produced using AA method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereMPSWrapper::AdAd (int m1, int m2, double& coefficient)
{
  unsigned long TmpState = this->ProdATemporaryState;
  if ((TmpState & (((unsigned long) (0x1)) << m2))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
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
  TmpState |= (((unsigned long) (0x1)) << m2);
  if ((TmpState & (((unsigned long) (0x1)) << m1))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
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
  TmpState |= (((unsigned long) (0x1)) << m1);
  return this->FindStateIndex(TmpState, NewLzMax);
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereMPSWrapper::AdA (int index, int m)
{
  SparseComplexMatrix TmpMatrix;
  if (m == 0)
    {
      TmpMatrix.Copy(this->NormalizedB1B1[0]);
      for (int j = 1; j <= this->LzMax; ++j)
	TmpMatrix.Multiply(this->NormalizedB0B0B1B1[j], this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements);
    }
  else
    {
      TmpMatrix.Copy(this->NormalizedB0B0B1B1[0]);
      for (int j = 1; j < m; ++j)
	TmpMatrix.Multiply(this->NormalizedB0B0B1B1[j], this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements);
      TmpMatrix.Multiply(NormalizedB1B1[m], TmpMatrixElements, TmpColumnIndices, TmpElements);
      for (int j = m + 1; j <= this->LzMax; ++j)
	TmpMatrix.Multiply(this->NormalizedB0B0B1B1[j], this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements);
   }
  Complex Tmp = 0.0;
  TmpMatrix.GetMatrixElement(this->MPSRowIndex, this->MPSColumnIndex, Tmp);
  return Tmp.Re / this->StateNormalization;
}




// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

// attention: check sign returned by this function!
int FermionOnSphereMPSWrapper::AdA (int index, int m, int n, double& coefficient)
{
  int StateLzMax = this->StateLzMax;
  unsigned long State = this->StateDescription;
  if ((n > StateLzMax) || ((State & (((unsigned long) (0x1)) << n)) == 0))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  int NewLzMax = StateLzMax;
  unsigned long TmpState = State;
  coefficient = this->SignLookUpTable[(TmpState >> n) & this->SignLookUpTableMask[n]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 16))  & this->SignLookUpTableMask[n + 16]];
#ifdef  __64_BITS__
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 32)) & this->SignLookUpTableMask[n + 32]];
  coefficient *= this->SignLookUpTable[(TmpState >> (n + 48)) & this->SignLookUpTableMask[n + 48]];
#endif
  TmpState &= ~(((unsigned long) (0x1)) << n);
  if ((TmpState != 0x0ul))
    {
      while ((TmpState >> NewLzMax) == 0)
	--NewLzMax;
    }
  else
    NewLzMax = 0;
  if ((TmpState & (((unsigned long) (0x1)) << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  if (m > NewLzMax)
    {
      NewLzMax = m;
    }
  else
    {
      coefficient *= this->SignLookUpTable[(TmpState >> m) & this->SignLookUpTableMask[m]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 16))  & this->SignLookUpTableMask[m + 16]];
#ifdef  __64_BITS__
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 32)) & this->SignLookUpTableMask[m + 32]];
      coefficient *= this->SignLookUpTable[(TmpState >> (m + 48)) & this->SignLookUpTableMask[m + 48]];
#endif
    }
  TmpState |= (((unsigned long) (0x1)) << m);
  return this->FindStateIndex(TmpState, NewLzMax);
}


// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

int FermionOnSphereMPSWrapper::FindStateIndex(unsigned long stateDescription, int lzmax)
{
  if (stateDescription == this->StateDescription)
    return 0;
  else
    return 1;
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereMPSWrapper::PrintState (ostream& Str, int state)
{
  Str << "MPS : ";
  unsigned long TmpState = this->StateDescription;
  for (int i = 0; i < this->NbrLzValue; ++i)
    Str << ((TmpState >> i) & ((unsigned long) 0x1)) << " ";
  return Str;
}

// print a given State using the monomial notation
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereMPSWrapper::PrintStateMonomial (ostream& Str, long state)
{
  Str << "MPS : ";
  unsigned long TmpState = this->StateDescription;
  Str << "[";
  int i = this->LzMax;
  while (((TmpState >> i) & 0x1ul) == 0x0ul)
    --i;
  Str << i;
  --i;
  for (; i >=0; --i)
    if (((TmpState >> i) & 0x1ul) != 0x0ul)
      Str << "," << i;
  Str << "]";
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion in the state
// currentLzMax = momentum maximum value for fermions that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

int FermionOnSphereMPSWrapper::GenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, int pos)
{
  return 1;
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnSphereMPSWrapper::GenerateLookUpTable(unsigned long memory)
{
  this->GenerateSignLookUpTable();
}

// generate look-up table for sign calculation
// 

void FermionOnSphereMPSWrapper::GenerateSignLookUpTable()
{
  // look-up tables for evaluating sign when applying creation/annihilation operators
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double [Size];
  int Count;
  int TmpNbr;
  for (int j = 0; j < Size; ++j)
    {
      Count = 0;
      TmpNbr = j;
      while (TmpNbr != 0)
	{
	  if (TmpNbr & 0x1)
	    ++Count;
	  TmpNbr >>= 1;
	}
      if (Count & 1)
	this->SignLookUpTable[j] = -1.0;
      else
	this->SignLookUpTable[j] = 1.0;
    }
#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#endif
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereMPSWrapper::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  return 1l;
}

