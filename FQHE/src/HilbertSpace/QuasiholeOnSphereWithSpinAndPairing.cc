////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2005 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//              class of fermions on sphere with spin and pairing             //
//      (i.e. Sz conservation but no conservation of the particle number)     //
//                                                                            //
//                        last modification : 12/04/2016                      //
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
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/QuasiholeOnSphereWithSpinAndPairing.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <cmath>
#include <bitset>
#include <cstdlib>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::bitset;


// default constructor
//

QuasiholeOnSphereWithSpinAndPairing::QuasiholeOnSphereWithSpinAndPairing()
{
}

// basic constructor
// 
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// totalSpin = twice the total spin value
// memory = amount of memory granted for precalculations

QuasiholeOnSphereWithSpinAndPairing::QuasiholeOnSphereWithSpinAndPairing (int kValue, int rValue, int totalLz, int lzMax, int totalSpin, unsigned long memory)
{
  this->NbrFermions = 0;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (this->NbrFermions + this->TotalSpin)/2;
  this->NbrFermionsDown = (this->NbrFermions - this->TotalSpin)/2;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  
  this->KValue = kValue;
  this->RValue = rValue;
  this->FermionFactor = 1;
  this->NbrFermionsUpMax = (this->KValue * (this->LzMax + 1 + this->RValue)) / (this->RValue + (this->KValue * this->FermionFactor)); 
  int MaxTotalLz = this->GetMaximalLzSingleLayer(this->NbrFermionsUpMax);
  this->NbrQuasiholeEntriesSingleLayer = this->GetLinearIndexSingleLayer(this->NbrFermionsUpMax, MaxTotalLz) + 1;
  this->NbrQuasiholesPerNPerLzSingleLayer = new int [this->NbrQuasiholeEntriesSingleLayer];
  this->SingleLayerIndices = new int [this->NbrQuasiholeEntriesSingleLayer];
  
  int TmpIndex = 0;
  for (int TmpNbrParticles = 0; TmpNbrParticles <= this->NbrFermionsUpMax; ++TmpNbrParticles)
  {
    MaxTotalLz = this->GetMaximalLzSingleLayer(TmpNbrParticles);
    for (int TmpLz = -MaxTotalLz; TmpLz <= MaxTotalLz; TmpLz += 2)
    {
      this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex] = 0;
      this->SingleLayerIndices[TmpIndex] = 0;
      ++TmpIndex;  
    }
  }

  char* TmpFileName = new char [512];
  sprintf (TmpFileName, "fermions_qh_states_k_%d_r_%d_nphi_%d.dat", this->KValue, this->RValue, this->LzMax);
  MultiColumnASCIIFile DimensionQuasiholeSubspaces;
  if (DimensionQuasiholeSubspaces.Parse(TmpFileName) == false)
    {
      DimensionQuasiholeSubspaces.DumpErrors(cout) << endl;
    } 
  
  int NbrNonZeroElements = DimensionQuasiholeSubspaces.GetNbrLines();
  int* NbrFermionsSingleLayer = DimensionQuasiholeSubspaces.GetAsIntegerArray (1);
  int* LzSingleLayer = DimensionQuasiholeSubspaces.GetAsIntegerArray (2);
  int* DimensionsSingleLayer = DimensionQuasiholeSubspaces.GetAsIntegerArray (3);
  
  
  for (int i = 0; i < NbrNonZeroElements; ++i)
  {
    TmpIndex = this->GetLinearIndexSingleLayer(NbrFermionsSingleLayer[i], LzSingleLayer[i]);
    this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex] = DimensionsSingleLayer[i];
  }
  delete[] NbrFermionsSingleLayer;
  delete[] LzSingleLayer;
  delete[] DimensionsSingleLayer;
  
  TmpIndex = 0;
  int SingleLayerDimension = 0 ;
  for (; TmpIndex < this->NbrQuasiholeEntriesSingleLayer; ++ TmpIndex)
  {
    this->SingleLayerIndices[TmpIndex] = SingleLayerDimension;
    SingleLayerDimension += this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex];
  }
  
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension();
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;

  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  this->NbrFermionUpFullSpace = new int [this->LargeHilbertSpaceDimension];
  this->LzValueUpFullSpace = new int [this->LargeHilbertSpaceDimension];
  
  this->LargeHilbertSpaceDimension = this->GenerateStates();
   
  RealSymmetricMatrix SingleLayerAnnihilationMatrix(SingleLayerDimension, true);
  RealSymmetricMatrix* SingleLayerAdAMatrix = new RealSymmetricMatrix[this->LzMax + 1];
  for (int i = 0; i < this->LzMax + 1; ++i)
    SingleLayerAdAMatrix[i] = RealSymmetricMatrix (SingleLayerDimension, true);
  RealMatrix TmpSingleLayerInteraction;
  int TmpIndexRight;
  int TmpIndexLeft;
  int TmpIndexRight1;
  int TmpIndexLeft1;
  double TmpInteraction;
  for (int TmpRightNbrParticles = 0; TmpRightNbrParticles <= this->NbrFermionsUpMax; ++TmpRightNbrParticles)
  {
    int MaxRightTotalLz = this->GetMaximalLzSingleLayer(TmpRightNbrParticles);
    for (int TmpRightLz = -MaxRightTotalLz; TmpRightLz <= MaxRightTotalLz; TmpRightLz += 2)
    {
      TmpIndexRight1 = this->GetLinearIndexSingleLayer(TmpRightNbrParticles, TmpRightLz);
      for (int OperatorLzValue = -this->LzMax; OperatorLzValue <= this->LzMax; OperatorLzValue += 2)
	{
	  int TmpLeftNbrParticles = TmpRightNbrParticles - 1;
	  int MaxLeftTotalLz = this->GetMaximalLzSingleLayer(TmpLeftNbrParticles);
	  int TmpLeftLz = TmpRightLz - OperatorLzValue;
	  if ((TmpLeftNbrParticles >= 0) && (TmpLeftLz >= -MaxLeftTotalLz) && (TmpLeftLz <= MaxLeftTotalLz))
	    {
	      TmpIndexLeft1 = this->GetLinearIndexSingleLayer(TmpLeftNbrParticles, TmpLeftLz);
	      sprintf (TmpFileName, "fermions_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_c_%d.mat", this->KValue, this->RValue, TmpRightNbrParticles, this->LzMax, TmpRightLz, OperatorLzValue);
	      if (TmpSingleLayerInteraction.ReadMatrix(TmpFileName) == false)
	      {
	      } 
	      else
	      {
		if (abs(this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexLeft1] - TmpSingleLayerInteraction.GetNbrRow()) > 0)
		  cout << this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexLeft1] << " " << TmpSingleLayerInteraction.GetNbrRow() << endl;
		if (abs(this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexRight1] - TmpSingleLayerInteraction.GetNbrColumn()) > 0)
		  cout << this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexRight1] << " " << TmpSingleLayerInteraction.GetNbrColumn() << endl;
		
		for (int i = 0; i < TmpSingleLayerInteraction.GetNbrRow(); ++i)
		{
		  for (int j = 0; j < TmpSingleLayerInteraction.GetNbrColumn(); ++j)
		  {
		    TmpIndexLeft = this->SingleLayerIndices[TmpIndexLeft1] + i;
		    TmpIndexRight = this->SingleLayerIndices[TmpIndexRight1] + j;
		    TmpSingleLayerInteraction.GetMatrixElement(i, j, TmpInteraction);
		    SingleLayerAnnihilationMatrix.AddToMatrixElement(TmpIndexLeft, TmpIndexRight, TmpInteraction);
		  }
		}
	      }
	    }
	    
	  
	  int ShiftedOperatorLzValue = (OperatorLzValue + this->LzMax) / 2;	  
	  sprintf (TmpFileName, "fermions_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d.mat", this->KValue, this->RValue, TmpRightNbrParticles, this->LzMax, TmpRightLz, OperatorLzValue);
	  
	  if (TmpSingleLayerInteraction.ReadMatrix(TmpFileName) == false)
	  {
	  } 
	  else
	  {
// 	    cout << this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexRight1] << " " << TmpSingleLayerInteraction.GetNbrRow() << endl;
// 	    cout << this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexRight1] << " " << TmpSingleLayerInteraction.GetNbrColumn() << endl;
		
	    for (int i = 0; i < TmpSingleLayerInteraction.GetNbrRow(); ++i)
	    {
	      for (int j = 0; j < TmpSingleLayerInteraction.GetNbrColumn(); ++j)
	      {
		TmpIndexLeft = this->SingleLayerIndices[TmpIndexRight1] + i;
		TmpIndexRight = this->SingleLayerIndices[TmpIndexRight1] + j;
		TmpSingleLayerInteraction.GetMatrixElement(i, j, TmpInteraction);
// 		cout << TmpIndexLeft << " " << TmpIndexRight << " " << ShiftedOperatorLzValue << " " << TmpInteraction << endl;
		SingleLayerAdAMatrix[ShiftedOperatorLzValue].AddToMatrixElement(TmpIndexLeft, TmpIndexRight, TmpInteraction);
	      }
	    }
	   }
	  
	  
	}
  
    }
  }

  
  delete[] TmpFileName;
  
  this->AnnihilationElementsOneLayer = SparseRealMatrix ((SingleLayerAnnihilationMatrix), 1.0e-13);
  this->AdAElementsOneLayer = new SparseRealMatrix[this->LzMax + 1];
  for (int i = 0; i <= this->LzMax; ++i)
    this->AdAElementsOneLayer[i] = SparseRealMatrix (SingleLayerAdAMatrix[i], 1.0e-13);
  delete[] SingleLayerAdAMatrix;
  
  
  
  cout.precision(14);
  
//   for (int TmpRightNbrParticles = 0; TmpRightNbrParticles <= this->NbrFermionsUpMax; ++TmpRightNbrParticles)
//   {
//     int MaxRightTotalLz = this->GetMaximalLzSingleLayer(TmpRightNbrParticles);
//     for (int TmpRightLz = -MaxRightTotalLz; TmpRightLz <= MaxRightTotalLz; TmpRightLz += 2)
//     {
//       TmpIndexRight1 = this->GetLinearIndexSingleLayer(TmpRightNbrParticles, TmpRightLz);
//       cout << TmpRightNbrParticles << " " << TmpRightLz << " " << this->SingleLayerIndices[TmpIndexRight1] << endl;
//       
//     }
//   }
//   cout << SingleLayerDimension << endl;
//   double Tmp;
//   for (int i = 0; i < SingleLayerDimension; ++i)
//   {
//     this->AdAElementsOneLayer[0].GetMatrixElement(i, i, Tmp);
//     cout << Tmp << " " ;
//     cout << endl;
//   }
  
  
  this->StateDescription = 0;
  this->StateHighestBit = 0;
  
  this->Flag.Initialize();
  this->TargetSpace = this;
  

  
//   int* LeftIndices;
//   double* InteractionElements;
//   int NbrElements;
//   for (int i = 0; i < this->HilbertSpaceDimension; ++i)
//   {
//     for (int OperatorLzValue = -this->LzMax; OperatorLzValue <= this->LzMax; OperatorLzValue += 2)
//     {
//       NbrElements = this->AduAu(i, OperatorLzValue, LeftIndices, InteractionElements);
//       cout  << " cd" << OperatorLzValue << "c" << OperatorLzValue << " |" << i << "> : (" << NbrElements << ") " ; 
//       for (int j = 0; j < NbrElements; ++j)
//       {
// // 	cout << j << " " ;
// 	cout << LeftIndices[j] << "/" << InteractionElements[j] << " ";  
//       }
//       cout << endl;
//     }
//   }
  
  
  

//   this->GenerateLookUpTable(memory);


#ifdef __DEBUG__
//   long UsedMemory = 0;
//   UsedMemory += this->LargeHilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
//   cout << "memory requested for Hilbert space = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;
//   UsedMemory = this->NbrLzValue * sizeof(int);
//   if (this->NbrFermions > 0)
//     UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
//   cout << "memory requested for lookup table = ";
//   if (UsedMemory >= 1024)
//     if (UsedMemory >= 1048576)
//       cout << (UsedMemory >> 20) << "Mo" << endl;
//     else
//       cout << (UsedMemory >> 10) << "ko" <<  endl;
//   else
//     cout << UsedMemory << endl;

#endif
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

QuasiholeOnSphereWithSpinAndPairing::QuasiholeOnSphereWithSpinAndPairing(const QuasiholeOnSphereWithSpinAndPairing& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
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
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

QuasiholeOnSphereWithSpinAndPairing::~QuasiholeOnSphereWithSpinAndPairing ()
{
  delete[] SingleLayerIndices;
  delete[] NbrFermionUpFullSpace;
  delete[] LzValueUpFullSpace;
  delete[] AdAElementsOneLayer;
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

QuasiholeOnSphereWithSpinAndPairing& QuasiholeOnSphereWithSpinAndPairing::operator = (const QuasiholeOnSphereWithSpinAndPairing& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;
    }
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
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

AbstractHilbertSpace* QuasiholeOnSphereWithSpinAndPairing::Clone()
{
  return new QuasiholeOnSphereWithSpinAndPairing(*this);
}


  
// generate all states corresponding to the constraints
// 
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// totalSpin = twice the total spin
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long QuasiholeOnSphereWithSpinAndPairing::GenerateStates()
{
  long TmpDimension = 0;
  int TmpIndex = 0;
  int TmpComplementaryIndex;
  int TmpDimensionSector;
  
  int TmpComplementaryNbrParticles;
  int TmpComplementaryLz;
  
  int TmpIndex1 = 0;
  for (int TmpNbrParticles = 0; TmpNbrParticles <= this->NbrFermionsUpMax; ++TmpNbrParticles)
  {
    int MaxTotalLz = this->GetMaximalLzSingleLayer(TmpNbrParticles);
    TmpComplementaryNbrParticles = TmpNbrParticles - this->TotalSpin;
    if ((TmpComplementaryNbrParticles >= 0) && (TmpComplementaryNbrParticles <= this->NbrFermionsUpMax))
    {
      int MaxComplementaryLz = this->GetMaximalLzSingleLayer(TmpComplementaryNbrParticles);
      for (int TmpLz = -MaxTotalLz; TmpLz <= MaxTotalLz; TmpLz += 2)
      {
	TmpComplementaryLz = TmpLz - this->TotalLz;
	if (abs(TmpComplementaryLz) <= MaxComplementaryLz)
	{
	  TmpIndex = this->GetLinearIndexSingleLayer(TmpNbrParticles, TmpLz);
	  TmpComplementaryIndex = this->GetLinearIndexSingleLayer(TmpComplementaryNbrParticles, TmpComplementaryLz);
	  TmpDimensionSector = this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex] * this->NbrQuasiholesPerNPerLzSingleLayer[ TmpComplementaryIndex];
	  TmpDimension += TmpDimensionSector; 
	  for (int i = 0; i < TmpDimensionSector; ++i)
	  {
	    this->NbrFermionUpFullSpace[TmpIndex1] = TmpNbrParticles;
	    this->LzValueUpFullSpace[TmpIndex1] = TmpLz;
	    ++TmpIndex1;
	  
	  }
	}
      }
    }
  }
  return TmpDimension;      
}


// evaluate Hilbert space dimension
//
// posMax = highest position for next particle to be placed
// totalLz = momentum total value
// totalSpin = twice the total spin
// return value = Hilbert space dimension

long QuasiholeOnSphereWithSpinAndPairing::EvaluateHilbertSpaceDimension()
{
  long TmpDimension = 0;
  int TmpIndex = 0;
  int TmpComplementaryIndex;
  
  int TmpComplementaryNbrParticles;
  int TmpComplementaryLz;
  for (int TmpNbrParticles = 0; TmpNbrParticles <= this->NbrFermionsUpMax; ++TmpNbrParticles)
  {
    int MaxTotalLz = this->GetMaximalLzSingleLayer(TmpNbrParticles);
    TmpComplementaryNbrParticles = TmpNbrParticles - this->TotalSpin;
    if ((TmpComplementaryNbrParticles >= 0) && (TmpComplementaryNbrParticles <= this->NbrFermionsUpMax))
    {
      int MaxComplementaryLz = this->GetMaximalLzSingleLayer(TmpComplementaryNbrParticles);
      for (int TmpLz = -MaxTotalLz; TmpLz <= MaxTotalLz; TmpLz += 2)
      {
	TmpComplementaryLz = TmpLz - this->TotalLz;
	if (abs(TmpComplementaryLz) <= MaxComplementaryLz)
	{
	  TmpComplementaryIndex = this->GetLinearIndexSingleLayer(TmpComplementaryNbrParticles, TmpComplementaryLz);
	  TmpIndex = this->GetLinearIndexSingleLayer(TmpNbrParticles, TmpLz);
// 	  this->FullSystemIndices[TmpIndex] = TmpDimension;
	  TmpDimension += this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex] * this->NbrQuasiholesPerNPerLzSingleLayer[ TmpComplementaryIndex];  
	}
      }
    }
  }
  return TmpDimension; 
}


// find state index from the value of each number in one layer
//
// nbrParticlesUp = number of particles with up spin
// lzValueUp = value of the angular momentum for up spins
// alpha = index of state with up spin
// beta = index of state with down spin
// return value = corresponding index, -1 if an error occured
int QuasiholeOnSphereWithSpinAndPairing::FindStateIndex(int nbrParticlesUp, int lzValueUp, int alpha, int beta)
{
  int TmpStateIndex;
  int TmpMinIndex = 0;
  int TmpMaxIndex = this->LargeHilbertSpaceDimension;
  int TmpValue;
  while (TmpMinIndex < TmpMaxIndex)
  {
    TmpStateIndex = (TmpMaxIndex + TmpMinIndex) / 2;
    TmpValue = this->NbrFermionUpFullSpace[TmpStateIndex];
    if (TmpValue == nbrParticlesUp)
      break;
    if (TmpValue > nbrParticlesUp)
      TmpMaxIndex = TmpStateIndex;
    if (TmpValue < nbrParticlesUp)
      TmpMinIndex = TmpStateIndex;   
  }
  TmpMinIndex = TmpStateIndex;
  TmpMaxIndex = TmpStateIndex;
  while ((TmpMinIndex > 0) && (this->NbrFermionUpFullSpace[TmpMinIndex - 1] == nbrParticlesUp))
    --TmpMinIndex;
  while ((TmpMaxIndex < this->LargeHilbertSpaceDimension - 1) && (this->NbrFermionUpFullSpace[TmpMaxIndex + 1] == nbrParticlesUp))
    ++TmpMaxIndex;
  
  int TmpMinIndex1 = TmpMinIndex;
  ++TmpMaxIndex;
  while (TmpMinIndex < TmpMaxIndex)
  {
    TmpStateIndex = (TmpMaxIndex + TmpMinIndex) / 2;
    TmpValue = this->LzValueUpFullSpace[TmpStateIndex];
    if (TmpValue == lzValueUp)
      break;
    if (TmpValue > lzValueUp)
      TmpMaxIndex = TmpStateIndex;
    if (TmpValue < lzValueUp)
      TmpMinIndex = TmpStateIndex;   
  }
  
  while ((TmpStateIndex > 0) && (TmpStateIndex > TmpMinIndex1) && (this->LzValueUpFullSpace[TmpStateIndex - 1] == lzValueUp))
    --TmpStateIndex;
  TmpStateIndex += (alpha * (this->NbrQuasiholesPerNPerLzSingleLayer[this->GetLinearIndexSingleLayer(nbrParticlesUp, lzValueUp)]) + beta);
  return TmpStateIndex;
}

// find the values of beta in each layer for a given state
//
// index = state index
// alpha = reference on the value of alpha (up spin)
// beta = reference on the value of beta (down spsin)
void QuasiholeOnSphereWithSpinAndPairing::FindBetaIndices(int index, int& alpha, int& beta)
{
  int TmpMinIndex = index;
  int nbrParticlesUp = this->NbrFermionUpFullSpace[index];
  int lzValueUp = this->LzValueUpFullSpace[index];
  
  while ((TmpMinIndex > 0) && (this->NbrFermionUpFullSpace[TmpMinIndex - 1] == nbrParticlesUp) && (this->LzValueUpFullSpace[TmpMinIndex - 1] == lzValueUp))
    --TmpMinIndex;
//   cout << TmpMaxIndex << " " << this->LargeHilbertSpaceDimension << " " << this->NbrFermionUpFullSpace[TmpMaxIndex + 1] << " " << this->LzValueUpFullSpace[TmpMaxIndex + 1] <<  " " << nbrParticlesUp << " " << lzValueUp << endl;
  
 int TmpNbr = this->NbrQuasiholesPerNPerLzSingleLayer[this->GetLinearIndexSingleLayer(nbrParticlesUp, lzValueUp)];
 if (TmpNbr != 0)
 {
    alpha = (index - TmpMinIndex) / TmpNbr;
    beta = (index - TmpMinIndex) % TmpNbr;
  }
  else
  {
    alpha = 0;
    beta = 0;
  }
//   cout << index << " " << TmpMinIndex << " " << TmpMaxIndex << " " << alpha << " " << beta << endl;
}

// apply a_u_m a_d_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for destruction operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state
int QuasiholeOnSphereWithSpinAndPairing::AuAd (int index, int m, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements;
  
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index] - 1;
  int LzUp = this->LzValueUpFullSpace[index] - m;  
  if ((NbrParticlesUp < 0) || (abs(LzUp) > this->GetMaximalLzSingleLayer(NbrParticlesUp)))
    return 0;
  
  int TmpIndexUpLeft = this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUp);
  int TmpIndexUpRight = this->GetLinearIndexSingleLayer(NbrParticlesUp + 1, LzUp + m);
  
  int NbrParticlesDown = NbrParticlesUp - this->TotalSpin;
  int LzDown = LzUp - this->TotalLz;
  if ((NbrParticlesDown < 0) || (abs(LzDown) > this->GetMaximalLzSingleLayer(NbrParticlesDown)))
    return 0;
  int TmpIndexDownLeft = this->GetLinearIndexSingleLayer(NbrParticlesDown, LzDown);
  int TmpIndexDownRight = this->GetLinearIndexSingleLayer(NbrParticlesDown + 1, LzDown + m);
  
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  
  NumberCouplingElements = this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexUpLeft] * this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft];
  
  if (NumberCouplingElements == 0)
    return 0;
  
  leftIndices = new int [NumberCouplingElements];
  interactionElements = new double [NumberCouplingElements];
  
  for (int i = 0; i < this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexUpLeft]; ++i)
  {
    for (int j = 0; j < this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft]; ++j)
    {
      int TmpIndex = i * this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft] + j;
      leftIndices[TmpIndex] = this->FindStateIndex(NbrParticlesUp, LzUp, i, j);
//       cout << NumberCouplingElements << " " << i << " " << j << " " << TmpIndex << " " << this->FindStateIndex(NbrParticlesUp, LzUp, i, j) << " " endl;
      this->AnnihilationElementsOneLayer.GetMatrixElement((this->SingleLayerIndices[TmpIndexUpLeft] + i), (this->SingleLayerIndices[TmpIndexUpRight] + BetaUp), interactionElements[TmpIndex]);
      
      
      this->AnnihilationElementsOneLayer.GetMatrixElement((this->SingleLayerIndices[TmpIndexDownLeft] + j), (this->SingleLayerIndices[TmpIndexDownRight] + BetaDown), Tmp);
      
      interactionElements[TmpIndex] = interactionElements[TmpIndex] * Tmp;
    }
  }
  
  return NumberCouplingElements;    
}

// apply a^\dagger_d_m a_d_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for destruction operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state
int QuasiholeOnSphereWithSpinAndPairing::AddAd (int index, int m, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements = 0;
  int ShiftedOperatorLzValue = (m + this->LzMax)/  2;
  
  int NbrParticlesDown = this->NbrFermionUpFullSpace[index] - this->TotalSpin;
  int LzDown = this->LzValueUpFullSpace[index] - this->TotalLz;
  int TmpIndexDown = this->GetLinearIndexSingleLayer(NbrParticlesDown, LzDown);  
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  
  NumberCouplingElements = this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDown];
  
  leftIndices = new int [NumberCouplingElements];
  interactionElements = new double [NumberCouplingElements];
  
 
  for (int i = 0; i < this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDown]; ++i)
  {
    leftIndices[i] = this->FindStateIndex(NbrParticlesDown, LzDown, i, BetaDown);
    this->AdAElementsOneLayer[ShiftedOperatorLzValue].GetMatrixElement((this->SingleLayerIndices[TmpIndexDown] + i), (this->SingleLayerIndices[TmpIndexDown] + BetaDown), interactionElements[i]);
    
  }
    
  return NumberCouplingElements;    
}


// apply a^\dagger_u_m a_u_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for destruction operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state
int QuasiholeOnSphereWithSpinAndPairing::AduAu (int index, int m, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements = 0;
  int ShiftedOperatorLzValue = (m + this->LzMax)/  2;
  
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index];
  int LzUp = this->LzValueUpFullSpace[index];
  int TmpIndexUp = this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUp);  
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  
  NumberCouplingElements = this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexUp];
  
  leftIndices = new int [NumberCouplingElements];
  interactionElements = new double [NumberCouplingElements];
  
 
  for (int i = 0; i < this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexUp]; ++i)
  {
    leftIndices[i] = this->FindStateIndex(NbrParticlesUp, LzUp, BetaUp, i);
    this->AdAElementsOneLayer[ShiftedOperatorLzValue].GetMatrixElement((this->SingleLayerIndices[TmpIndexUp] + i), (this->SingleLayerIndices[TmpIndexUp] + BetaUp), interactionElements[i]);
    
  }
    
  return NumberCouplingElements;    
}
