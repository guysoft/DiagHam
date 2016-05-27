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
// directory = optional path to data files
// totalNbrParticles = if positive, fix the total number of particles
// memory = amount of memory granted for precalculations

QuasiholeOnSphereWithSpinAndPairing::QuasiholeOnSphereWithSpinAndPairing (int kValue, int rValue, int totalLz, int lzMax, int totalSpin, 
									  const char* directory, int totalNbrParticles, unsigned long memory)
{
  this->NbrFermions = 0;
  if (totalNbrParticles >= 0)
    {
      this->NbrFermions = totalNbrParticles;
    }
  this->TotalLz = totalLz;
  this->TotalSpin = totalSpin;
  this->NbrFermionsUp = (this->NbrFermions + this->TotalSpin) / 2;
  this->NbrFermionsDown = (this->NbrFermions - this->TotalSpin) / 2;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  
  this->KValue = kValue;
  this->RValue = rValue;
  this->FermionFactor = 1;
  this->NbrFermionsUpMin = 0;
  this->NbrFermionsUpMax = (this->KValue * (this->LzMax + 1 + this->RValue)) / (this->RValue + (this->KValue * this->FermionFactor)); 
  if (totalNbrParticles >= 0)
    {
      this->NbrFermionsUpMin = this->NbrFermionsUp;
      this->NbrFermionsUpMax = this->NbrFermionsUp;      
    }
  int MaxTotalLz = this->GetMaximalLzSingleLayer(this->NbrFermionsUpMax);
  this->NbrQuasiholeEntriesSingleLayer = this->GetLinearIndexSingleLayer(this->NbrFermionsUpMax, MaxTotalLz) + 1;
  this->NbrQuasiholesPerNPerLzSingleLayer = new int [this->NbrQuasiholeEntriesSingleLayer];
  this->SingleLayerIndices = new int [this->NbrQuasiholeEntriesSingleLayer];
  this->FirstIndexWithNbrParticlesUpLzValueUp = new int[this->NbrQuasiholeEntriesSingleLayer];
  
  int TmpIndex = 0;
  for (int TmpNbrParticles = this->NbrFermionsUpMin; TmpNbrParticles <= this->NbrFermionsUpMax; ++TmpNbrParticles)
    {
      MaxTotalLz = this->GetMaximalLzSingleLayer(TmpNbrParticles);
      for (int TmpLz = -MaxTotalLz; TmpLz <= MaxTotalLz; TmpLz += 2)
	{
	  this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex] = 0;
	  this->SingleLayerIndices[TmpIndex] = 0;
	  ++TmpIndex;  
	}
    }
  
  char* TmpFileName;
  if (directory != 0)
    { 
      TmpFileName = new char [512 + strlen(directory)];
      sprintf (TmpFileName, "%s/fermions_qh_states_k_%d_r_%d_nphi_%d.dat", directory, this->KValue, this->RValue, this->LzMax);
    }
  else
    {
      TmpFileName = new char [512 + strlen(directory)];
      sprintf (TmpFileName, "fermions_qh_states_k_%d_r_%d_nphi_%d.dat", this->KValue, this->RValue, this->LzMax);
    }
  MultiColumnASCIIFile DimensionQuasiholeSubspaces;
  if (DimensionQuasiholeSubspaces.Parse(TmpFileName) == false)
    {
      cout << "Error: " << TmpFileName << " cannot be found" << endl;
      this->HilbertSpaceDimension = 0;
      return;
    } 
  
  int NbrNonZeroElements = DimensionQuasiholeSubspaces.GetNbrLines();
  int* NbrFermionsSingleLayer = DimensionQuasiholeSubspaces.GetAsIntegerArray (1);
  int* LzSingleLayer = DimensionQuasiholeSubspaces.GetAsIntegerArray (2);
  int* DimensionsSingleLayer = DimensionQuasiholeSubspaces.GetAsIntegerArray (3);
  
  
  for (int i = 0; i < NbrNonZeroElements; ++i)
    {
//        if ((NbrFermionsSingleLayer[i] >= this->NbrFermionsUpMin) && 
//  	  (NbrFermionsSingleLayer[i] <= this->NbrFermionsUpMax))
	{
	  TmpIndex = this->GetLinearIndexSingleLayer(NbrFermionsSingleLayer[i], LzSingleLayer[i]);
	  this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex] = DimensionsSingleLayer[i];
	}
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

//   for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
//     {
//       this->PrintState(cout, i) << endl;
//     }

  if (this->LargeHilbertSpaceDimension > 0l)
    {   
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
      for (int TmpRightNbrParticles = this->NbrFermionsUpMin; TmpRightNbrParticles <= this->NbrFermionsUpMax; ++TmpRightNbrParticles)
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
		      if (directory != 0)
			{
			  sprintf (TmpFileName, "%s/fermions_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_c_%d.mat", directory, this->KValue, this->RValue, 
				   TmpRightNbrParticles, this->LzMax, TmpRightLz, OperatorLzValue);
			}
		      else
			{ 
			  sprintf (TmpFileName, "fermions_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_c_%d.mat", this->KValue, this->RValue, TmpRightNbrParticles, this->LzMax, TmpRightLz, OperatorLzValue);
			}
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
			      for (int j = i; j < TmpSingleLayerInteraction.GetNbrColumn(); ++j)
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
		  if (directory != 0)
		    {
		      sprintf (TmpFileName, "%s/fermions_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d.mat", directory, this->KValue, this->RValue, 
			       TmpRightNbrParticles, this->LzMax, TmpRightLz, OperatorLzValue);
		    }
		  else
		    {
		      sprintf (TmpFileName, "fermions_qh_k_%d_r_%d_n_%d_nphi_%d_lz_%d_cdc_%d.mat", this->KValue, this->RValue, TmpRightNbrParticles, this->LzMax, TmpRightLz, OperatorLzValue);
		    }
		  if (TmpSingleLayerInteraction.ReadMatrix(TmpFileName) == false)
		    {
		    } 
		  else
		    {
		      // 	    cout << this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexRight1] << " " << TmpSingleLayerInteraction.GetNbrRow() << endl;
		      // 	    cout << this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexRight1] << " " << TmpSingleLayerInteraction.GetNbrColumn() << endl;
		      
		      for (int i = 0; i < TmpSingleLayerInteraction.GetNbrRow(); ++i)
			{
			  for (int j = i; j < TmpSingleLayerInteraction.GetNbrColumn(); ++j)
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
      this->AnnihilationElementsOneLayer = SparseRealMatrix ((SingleLayerAnnihilationMatrix), 1.0e-13);
      this->AdAElementsOneLayer = new SparseRealMatrix[this->LzMax + 1];
      for (int i = 0; i <= this->LzMax; ++i)
	this->AdAElementsOneLayer[i] = SparseRealMatrix (SingleLayerAdAMatrix[i], 1.0e-13);
      delete[] SingleLayerAdAMatrix;
    }

  delete[] TmpFileName;
  
  this->Flag.Initialize();

 
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
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->KValue = fermions.KValue;
  this->RValue = fermions.RValue;
  this->FermionFactor = fermions.FermionFactor;
  this->NbrQuasiholeEntriesSingleLayer = fermions.NbrQuasiholeEntriesSingleLayer;
  this->NbrQuasiholesPerNPerLzSingleLayer = fermions.NbrQuasiholesPerNPerLzSingleLayer;
  this->SingleLayerIndices = fermions.SingleLayerIndices;
  this->NbrFermionUpFullSpace = fermions.NbrFermionUpFullSpace;
  this->LzValueUpFullSpace = fermions.LzValueUpFullSpace;
  this->AnnihilationElementsOneLayer = fermions.AnnihilationElementsOneLayer;
  this->AdAElementsOneLayer = fermions.AdAElementsOneLayer;
  this->NbrFermionsUpMin = fermions.NbrFermionsUpMin;
  this->NbrFermionsUpMax = fermions.NbrFermionsUpMax; 
  this->FirstIndexWithNbrParticlesUpLzValueUp = fermions.FirstIndexWithNbrParticlesUpLzValueUp;
  this->MaximalNumberCouplingElements = fermions.MaximalNumberCouplingElements;
}

// destructor
//

QuasiholeOnSphereWithSpinAndPairing::~QuasiholeOnSphereWithSpinAndPairing ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] SingleLayerIndices;
      delete[] NbrFermionUpFullSpace;
      delete[] LzValueUpFullSpace;
      delete[] AdAElementsOneLayer;
      delete[] FirstIndexWithNbrParticlesUpLzValueUp;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

QuasiholeOnSphereWithSpinAndPairing& QuasiholeOnSphereWithSpinAndPairing::operator = (const QuasiholeOnSphereWithSpinAndPairing& fermions)
{
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->Flag = fermions.Flag;
  this->NbrFermions = fermions.NbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->TotalSpin = fermions.TotalSpin;
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->KValue = fermions.KValue;
  this->RValue = fermions.RValue;
  this->FermionFactor = fermions.FermionFactor;
  this->NbrQuasiholeEntriesSingleLayer = fermions.NbrQuasiholeEntriesSingleLayer;
  this->NbrQuasiholesPerNPerLzSingleLayer = fermions.NbrQuasiholesPerNPerLzSingleLayer;
  this->SingleLayerIndices = fermions.SingleLayerIndices;
  this->NbrFermionUpFullSpace = fermions.NbrFermionUpFullSpace;
  this->LzValueUpFullSpace = fermions.LzValueUpFullSpace;
  this->AnnihilationElementsOneLayer = fermions.AnnihilationElementsOneLayer;
  this->AdAElementsOneLayer = fermions.AdAElementsOneLayer;
  this->NbrFermionsUpMin = fermions.NbrFermionsUpMin;
  this->NbrFermionsUpMax = fermions.NbrFermionsUpMax; 
  this->FirstIndexWithNbrParticlesUpLzValueUp = fermions.FirstIndexWithNbrParticlesUpLzValueUp;
  this->MaximalNumberCouplingElements = fermions.MaximalNumberCouplingElements;
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
// return value = Hilbert space dimension      

long QuasiholeOnSphereWithSpinAndPairing::GenerateStates()
{
  long TmpDimension = 0;
  int TmpIndex = 0;
  int TmpComplementaryIndex;
  int TmpDimensionSector;
  this->MaximalNumberCouplingElements = 0;
  
  int TmpComplementaryNbrParticles;
  int TmpComplementaryLz;
  
  int TmpIndex1 = 0;
  for (int TmpNbrParticles = this->NbrFermionsUpMin; TmpNbrParticles <= this->NbrFermionsUpMax; ++TmpNbrParticles)
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
		  if (TmpDimensionSector > this->MaximalNumberCouplingElements)
		    this->MaximalNumberCouplingElements = TmpDimensionSector;
		  this->FirstIndexWithNbrParticlesUpLzValueUp[TmpIndex] = TmpDimension;
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
  this->MaximalNumberCouplingElements *= this->MaximalNumberCouplingElements;
  return TmpDimension;      
}



// evaluate Hilbert space dimension
//
// return value = Hilbert space dimension

long QuasiholeOnSphereWithSpinAndPairing::EvaluateHilbertSpaceDimension()
{
  long TmpDimension = 0;
  int TmpIndex = 0;
  int TmpComplementaryIndex;
  
  int TmpComplementaryNbrParticles;
  int TmpComplementaryLz;
  for (int TmpNbrParticles = this->NbrFermionsUpMin; TmpNbrParticles <= this->NbrFermionsUpMax; ++TmpNbrParticles)
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
		  TmpDimension += this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndex] * this->NbrQuasiholesPerNPerLzSingleLayer[TmpComplementaryIndex];  
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
  int SingleLayerIndex = this->GetLinearIndexSingleLayer(nbrParticlesUp, lzValueUp);
  int TmpStateIndex = this->FirstIndexWithNbrParticlesUpLzValueUp[SingleLayerIndex] + beta * this->NbrQuasiholesPerNPerLzSingleLayer[SingleLayerIndex] + alpha;
  
  return TmpStateIndex;
}

// find the values of beta in each layer for a given state
//
// index = state index
// alpha = reference on the value of alpha (up spin)
// beta = reference on the value of beta (down spsin)
void QuasiholeOnSphereWithSpinAndPairing::FindBetaIndices(int index, int& alpha, int& beta)
{
  int nbrParticlesUp = this->NbrFermionUpFullSpace[index];
  int lzValueUp = this->LzValueUpFullSpace[index];
  int TmpMinIndex = this->FirstIndexWithNbrParticlesUpLzValueUp[this->GetLinearIndexSingleLayer(nbrParticlesUp, lzValueUp)];
  int TmpNbr = this->NbrQuasiholesPerNPerLzSingleLayer[this->GetLinearIndexSingleLayer(nbrParticlesUp, lzValueUp)];
  if (TmpNbr != 0)
    {
      alpha = (index - TmpMinIndex) % TmpNbr;
      beta = (index - TmpMinIndex) / TmpNbr;
    }
  else
    {
      alpha = 0;
      beta = index - TmpMinIndex;
    }
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
  int lz = 2 * m - this->LzMax;
  
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index] - 1;
  int LzUp = this->LzValueUpFullSpace[index] - lz;  
  if ((NbrParticlesUp < this->NbrFermionsUpMin) || (abs(LzUp) > this->GetMaximalLzSingleLayer(NbrParticlesUp)))
    return 0;
  
  int TmpIndexUpLeft = this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUp);
  int TmpIndexUpRight = this->GetLinearIndexSingleLayer(NbrParticlesUp + 1, LzUp + lz);
  
  int NbrParticlesDown = NbrParticlesUp - this->TotalSpin;
  int LzDown = LzUp - this->TotalLz;
  if ((NbrParticlesDown < 0) || (abs(LzDown) > this->GetMaximalLzSingleLayer(NbrParticlesDown)))
    return 0;
  int TmpIndexDownLeft = this->GetLinearIndexSingleLayer(NbrParticlesDown, LzDown);
  int TmpIndexDownRight = this->GetLinearIndexSingleLayer(NbrParticlesDown + 1, LzDown + lz);
  
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  double Tmp1;
  
  int TmpNbrCouplingElements = 0;
  int TmpBetaLeftUp;
  int TmpBetaLeftDown;
  NumberCouplingElements = this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexUpLeft] * this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft];
  
  if (NumberCouplingElements == 0)
    return 0;
    
  for (int i = 0; i < NumberCouplingElements; ++i)
    {
      TmpBetaLeftUp = i / this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft];
      TmpBetaLeftDown = i % (this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft]);
      this->AnnihilationElementsOneLayer.GetMatrixElement((this->SingleLayerIndices[TmpIndexUpLeft] + TmpBetaLeftUp), (this->SingleLayerIndices[TmpIndexUpRight] + BetaUp),  Tmp);
      if (Tmp != 0.0)
	{
	  this->AnnihilationElementsOneLayer.GetMatrixElement((this->SingleLayerIndices[TmpIndexDownLeft] + TmpBetaLeftDown), (this->SingleLayerIndices[TmpIndexDownRight] + BetaDown), Tmp1);
	  if (Tmp1 != 0.0)
	    {
	      interactionElements[TmpNbrCouplingElements] = Tmp * Tmp1;
	      leftIndices[TmpNbrCouplingElements] = this->FindStateIndex(NbrParticlesUp, LzUp, TmpBetaLeftUp, TmpBetaLeftDown);
	      ++TmpNbrCouplingElements;
	    }	  
	} 
    }
  return TmpNbrCouplingElements;    
}


// apply a_u_m a_d_m to a given state
//
// index = index of the state on which the operator has to be applied
// m = index for destruction operators
// leftIndices = reference to an array containing the indices of the resulting states
// interactionElements = reference to an array containing the matrix elements 
// return value = number of left states that are connected to the initial state
int QuasiholeOnSphereWithSpinAndPairing::AduAdd (int index, int m, int*& leftIndices, double*& interactionElements)
{
  int NumberCouplingElements;
  int lz = 2 * m - this->LzMax;
  
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index] + 1;
  int LzUp = this->LzValueUpFullSpace[index] + lz;  
  if ((NbrParticlesUp > this->NbrFermionsUpMax) || (abs(LzUp) > this->GetMaximalLzSingleLayer(NbrParticlesUp)))
    return 0;
  
  int TmpIndexUpLeft = this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUp);
  int TmpIndexUpRight = this->GetLinearIndexSingleLayer(NbrParticlesUp - 1, LzUp - lz);
  
  int NbrParticlesDown = NbrParticlesUp - this->TotalSpin;
  int LzDown = LzUp - this->TotalLz;
  if ((NbrParticlesDown < 0) || (abs(LzDown) > this->GetMaximalLzSingleLayer(NbrParticlesDown)))
    return 0;
  int TmpIndexDownLeft = this->GetLinearIndexSingleLayer(NbrParticlesDown, LzDown);
  int TmpIndexDownRight = this->GetLinearIndexSingleLayer(NbrParticlesDown - 1, LzDown - lz);
  
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  double Tmp1;
  int TmpBetaLeftUp;
  int TmpBetaLeftDown;
  int TmpNbrCouplingElements = 0;
  
  NumberCouplingElements = this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexUpLeft] * this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft];
  
  if (NumberCouplingElements == 0)
    return 0;
  
  for (int i = 0; i < NumberCouplingElements; ++i)
    {
      TmpBetaLeftUp = i / this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft];
      TmpBetaLeftDown = i % (this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDownLeft]);
      this->AnnihilationElementsOneLayer.GetMatrixElement((this->SingleLayerIndices[TmpIndexUpRight] + BetaUp), (this->SingleLayerIndices[TmpIndexUpLeft] + TmpBetaLeftUp),  Tmp);
      if (Tmp != 0.0)
	{
	  this->AnnihilationElementsOneLayer.GetMatrixElement((this->SingleLayerIndices[TmpIndexDownRight] + BetaDown), (this->SingleLayerIndices[TmpIndexDownLeft] + TmpBetaLeftDown), Tmp1);
	  if (Tmp1 != 0.0)
	    {
	      interactionElements[TmpNbrCouplingElements] = Tmp * Tmp1;
	      leftIndices[TmpNbrCouplingElements] = this->FindStateIndex(NbrParticlesUp, LzUp, TmpBetaLeftUp, TmpBetaLeftDown);
	      ++TmpNbrCouplingElements;
	    }	  
	}
    }

  return TmpNbrCouplingElements;    
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
  int ShiftedOperatorLzValue = m;
  
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index];
  int NbrParticlesDown = NbrParticlesUp - this->TotalSpin;
  int LzUp = this->LzValueUpFullSpace[index];
  int LzDown = LzUp - this->TotalLz;
  int TmpIndexDown = this->GetLinearIndexSingleLayer(NbrParticlesDown, LzDown); 
  int NbrStatesUp = this->NbrQuasiholesPerNPerLzSingleLayer [this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUp)];
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  int TmpIndex = this->FindStateIndex(NbrParticlesUp, LzUp, BetaUp, 0);  
  NumberCouplingElements = 0;
   
  for (int i = 0; i < this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexDown]; ++i)
  {
    this->AdAElementsOneLayer[ShiftedOperatorLzValue].GetMatrixElement((this->SingleLayerIndices[TmpIndexDown] + i), (this->SingleLayerIndices[TmpIndexDown] + BetaDown), Tmp);
    if (Tmp != 0.0)
    {
      interactionElements[NumberCouplingElements] = Tmp;
      leftIndices[NumberCouplingElements] = TmpIndex;
      ++NumberCouplingElements;
    }
    TmpIndex += NbrStatesUp;    
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
  int ShiftedOperatorLzValue = m;
  
  int NbrParticlesUp = this->NbrFermionUpFullSpace[index];
  int LzUp = this->LzValueUpFullSpace[index];
  int TmpIndexUp = this->GetLinearIndexSingleLayer(NbrParticlesUp, LzUp);  
  
  int BetaUp;
  int BetaDown;  
  this->FindBetaIndices(index, BetaUp, BetaDown);
  
  double Tmp;
  int TmpIndex = this->FindStateIndex(NbrParticlesUp, LzUp, 0, BetaDown);  
  NumberCouplingElements = 0;
  for (int i = 0; i < this->NbrQuasiholesPerNPerLzSingleLayer[TmpIndexUp]; ++i)
    {
      this->AdAElementsOneLayer[ShiftedOperatorLzValue].GetMatrixElement((this->SingleLayerIndices[TmpIndexUp] + i), (this->SingleLayerIndices[TmpIndexUp] + BetaUp), Tmp);
      if (Tmp != 0.0)
	{
	  interactionElements[NumberCouplingElements] = Tmp;
	  leftIndices[NumberCouplingElements] = TmpIndex;
	  ++NumberCouplingElements;
	}    
      ++TmpIndex;
    }    
  return NumberCouplingElements;    
}


// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 
ostream& QuasiholeOnSphereWithSpinAndPairing::PrintState (ostream& Str, int state)
{
  int NbrFermionsUp = this->NbrFermionUpFullSpace[state];
  int NbrFermionsDown = NbrFermionsUp - this->TotalSpin;
  
  int lzValueUp = this->LzValueUpFullSpace[state];
  int lzValueDown = lzValueUp - this->TotalLz;
  
  int alpha;
  int beta;
  this->FindBetaIndices(state, alpha, beta);
   
  Str << "(N=" << NbrFermionsUp << ", lz=" << lzValueUp << ", " << alpha << ") (N=" << NbrFermionsDown << ", lz=" << lzValueDown << ", " << beta << ")";
  
  return Str;
}
