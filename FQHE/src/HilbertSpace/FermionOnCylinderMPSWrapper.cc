////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of fermions on cylinder that allow to use MPS with operator     //
//                                                                            //
//                        last modification : 15/10/2012                      //
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
#include "HilbertSpace/FermionOnCylinderMPSWrapper.h"
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

FermionOnCylinderMPSWrapper::FermionOnCylinderMPSWrapper()
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

FermionOnCylinderMPSWrapper::FermionOnCylinderMPSWrapper (int nbrFermions, int& totalLz, int lzMax, int* referenceState,
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
  this->GenerateLookUpTable(10000000);
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


  this-> NbrPrecalculatedMatrixProducts = this->LzMax + 1;

  this->NormalizedB1B1 = new SparseComplexMatrix [1];
  this->NormalizedB0B1 = new SparseComplexMatrix [1];
  this->NormalizedB1B0 = new SparseComplexMatrix [1];
  this->NormalizedB0B0B1B1 = new SparseComplexMatrix [NbrPrecalculatedMatrixProducts];
  
  long TmpMemory = (((long) SparseTensorProductBMatrices[1][1].GetNbrRow()) * 
		    ((long) SparseTensorProductBMatrices[1][1].GetNbrColumn())) / 100l;
  cout << "Requested memory for sparse matrix multiplications = " << ((TmpMemory * (2l * sizeof(double) + sizeof(int))) >> 20) << "Mb" << endl;
  this->TmpMatrixElements = new Complex [TmpMemory];
  this->TmpColumnIndices = new int [TmpMemory];
  this->TmpElements = new Complex [SparseTensorProductBMatrices[1][1].GetNbrRow()];

  this->NormalizedB0B0B1B1[0] = SparseComplexMatrixLinearCombination(1.0, SparseTensorProductBMatrices[0][0], 1.0, SparseTensorProductBMatrices[1][1]);
  this->NormalizedB1B1[0].Copy(SparseTensorProductBMatrices[1][1]);
  this->NormalizedB0B1[0].Copy(SparseTensorProductBMatrices[0][1]);
  this->NormalizedB1B0[0].Copy(SparseTensorProductBMatrices[1][0]);

  unsigned long PrecalculationMemory = this->NormalizedB0B0B1B1[0].GetAllocatedMemory();
  this->NbrPrecalculatedMatrixProducts = 0;
  for (int i = 1; (i <= this->LzMax) && (PrecalculationMemory < memory); ++i)
    {
      this->NormalizedB0B0B1B1[i].Copy(this->NormalizedB0B0B1B1[i - 1]);
      this->NormalizedB0B0B1B1[i].Multiply(this->NormalizedB0B0B1B1[0], this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements);
      PrecalculationMemory += this->NormalizedB0B0B1B1[0].GetAllocatedMemory();
      ++this->NbrPrecalculatedMatrixProducts;
    }
  cout << "Requested memory for precalculations = " << (PrecalculationMemory >> 20) << "Mb" << endl;
  SparseComplexMatrix TmpMatrixNorm;
  TmpMatrixNorm.Copy(this->NormalizedB0B0B1B1[this->NbrPrecalculatedMatrixProducts - 1]);
  int NbrSteps = (this->LzMax + 1) / this->NbrPrecalculatedMatrixProducts;
  for (int i = 1; i < NbrSteps; ++i)
    TmpMatrixNorm.Multiply(this->NormalizedB0B0B1B1[this->NbrPrecalculatedMatrixProducts - 1], this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements);
  NbrSteps = (this->LzMax + 1) % this->NbrPrecalculatedMatrixProducts;
  if (NbrSteps > 0)
    TmpMatrixNorm.Multiply(this->NormalizedB0B0B1B1[NbrSteps - 1], this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements);
    
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

FermionOnCylinderMPSWrapper::FermionOnCylinderMPSWrapper(const FermionOnCylinderMPSWrapper& fermions)
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
  this->NbrPrecalculatedMatrixProducts = fermions.NbrPrecalculatedMatrixProducts;
  long TmpMemory = (((long) this->NormalizedB0B0B1B1[0].GetNbrRow()) * 
		    ((long)  this->NormalizedB0B0B1B1[0].GetNbrColumn())) / 100l;
  this->TmpMatrixElements = new Complex [TmpMemory];
  this->TmpColumnIndices = new int [TmpMemory];
  this->TmpElements = new Complex [this->NormalizedB0B0B1B1[0].GetNbrRow()];
}

// destructor
//

FermionOnCylinderMPSWrapper::~FermionOnCylinderMPSWrapper ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnCylinderMPSWrapper& FermionOnCylinderMPSWrapper::operator = (const FermionOnCylinderMPSWrapper& fermions)
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
  this->NbrPrecalculatedMatrixProducts = fermions.NbrPrecalculatedMatrixProducts;
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

AbstractHilbertSpace* FermionOnCylinderMPSWrapper::Clone()
{
  return new FermionOnCylinderMPSWrapper(*this);
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

int FermionOnCylinderMPSWrapper::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
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
	    TmpMatrix.Multiply(this->NormalizedB1B1[0], this->TmpMatrixElements, this->TmpColumnIndices,  this->TmpElements);
	  else
	    TmpMatrix.Multiply(this->NormalizedB0B1[0], this->TmpMatrixElements, this->TmpColumnIndices,  this->TmpElements);
	}
      else
	{
	  if (j == m2)
	    {
	      if ((m2 == n1) || (m2 == n2))
		TmpMatrix.Multiply(this->NormalizedB1B1[0], this->TmpMatrixElements, this->TmpColumnIndices,  this->TmpElements);
	      else
		TmpMatrix.Multiply(this->NormalizedB0B1[0], this->TmpMatrixElements, this->TmpColumnIndices,  this->TmpElements);
	    }
	  else
	    {
	      if ((j == n1) || (j == n2))
		{
		  TmpMatrix.Multiply(this->NormalizedB1B0[0], this->TmpMatrixElements, this->TmpColumnIndices,  this->TmpElements);
		}
	      else
		{
		  TmpMatrix.Multiply(this->NormalizedB0B0B1B1[0], this->TmpMatrixElements, this->TmpColumnIndices,  this->TmpElements);
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

int FermionOnCylinderMPSWrapper::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
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

double FermionOnCylinderMPSWrapper::ProdA (int index, int* n, int nbrIndices)
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

double FermionOnCylinderMPSWrapper::AA (int index, int n1, int n2)
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

int FermionOnCylinderMPSWrapper::ProdAd (int* m, int nbrIndices, double& coefficient)
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

int FermionOnCylinderMPSWrapper::AdAd (int m1, int m2, double& coefficient)
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

double FermionOnCylinderMPSWrapper::AdA (int index, int m)
{
  SparseComplexMatrix TmpMatrix;
  if (m == 0)
    {
      TmpMatrix.Copy(this->NormalizedB1B1[0]);
      int NbrSteps = this->LzMax / this->NbrPrecalculatedMatrixProducts;
      for (int j = 0; j < NbrSteps; ++j)
	TmpMatrix.Multiply(this->NormalizedB0B0B1B1[this->NbrPrecalculatedMatrixProducts - 1], this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements);
      NbrSteps = this->LzMax % this->NbrPrecalculatedMatrixProducts;
      if (NbrSteps > 0)
	TmpMatrix.Multiply(this->NormalizedB0B0B1B1[NbrSteps - 1], this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements);
    }
  else
    {
      int NbrSteps = m / this->NbrPrecalculatedMatrixProducts;
      if (NbrSteps == 0)
	TmpMatrix.Copy(this->NormalizedB0B0B1B1[m - 1]);
      else
	{
	  TmpMatrix.Copy(this->NormalizedB0B0B1B1[this->NbrPrecalculatedMatrixProducts - 1]);
	  for (int j = 1; j < NbrSteps; ++j)
	    {
	      TmpMatrix.Multiply(this->NormalizedB0B0B1B1[this->NbrPrecalculatedMatrixProducts - 1], this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements);
	    }
	  NbrSteps = m % this->NbrPrecalculatedMatrixProducts;
	  if (NbrSteps > 0)
	    {
	      TmpMatrix.Multiply(this->NormalizedB0B0B1B1[NbrSteps - 1], this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements);
	    }
	}
      TmpMatrix.Multiply(NormalizedB1B1[0], TmpMatrixElements, TmpColumnIndices, TmpElements);
      if (m < this->LzMax)
	{
	  NbrSteps = (this->LzMax - m) / this->NbrPrecalculatedMatrixProducts;
	  for (int j = 0; j < NbrSteps; ++j)
	    TmpMatrix.Multiply(this->NormalizedB0B0B1B1[this->NbrPrecalculatedMatrixProducts - 1], this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements);
	  NbrSteps = (this->LzMax - m) % this->NbrPrecalculatedMatrixProducts;
	  if (NbrSteps > 0)
	    TmpMatrix.Multiply(this->NormalizedB0B0B1B1[NbrSteps - 1], this->TmpMatrixElements, this->TmpColumnIndices, this->TmpElements);
	}
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
int FermionOnCylinderMPSWrapper::AdA (int index, int m, int n, double& coefficient)
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


