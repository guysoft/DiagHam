////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of abstract quantum Hall hamiltonian associated          //
//            to particles on a sphere with n-body interaction terms          //
//                                                                            //
//                        last modification : 22/09/2004                      //
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
#include "Hamiltonian/QHEHamiltonian/AbstractQHEOnSphereNBodyInteractionHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Complex.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// destructor
//

AbstractQHEOnSphereNBodyInteractionHamiltonian::~AbstractQHEOnSphereNBodyInteractionHamiltonian()
{
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereNBodyInteractionHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
										int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int* MIndices;
      int* NIndices;
      double* TmpInteraction;
      double Coefficient2;
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true)
	  {
	    for (int i = firstComponent; i < LastComponent; ++i)
	      {
		for (int j = this->MinSumIndices; j <= this->MaxSumIndices; ++j)
		  {
		    TmpInteraction = this->NBodyInteractionFactors[k][j];
		    int Lim = NbrSortedIndicesPerSum[k][j];
		    MIndices = this->SortedIndicesPerSum[k][j];
		    for (int i1 = 0; i1 < Lim; ++i1)
		      {
			Coefficient2 = vSource[i] * TmpInteraction[i1];
			NIndices = this->SortedIndicesPerSum[k][j];
			for (int i2 = 0; i2 < Lim; ++i2)
			  {
			    Index = TmpParticles->ProdAdProdA(i, MIndices, NIndices, k, Coefficient);
			    if (Index < Dim)
			      vDestination[Index] += Coefficient * TmpInteraction[i2] * Coefficient2;
			    NIndices += k;
			  }
			MIndices += k;
		      }
		  }
	      }
	  }	    
      for (int i = firstComponent; i < LastComponent; ++i)
	vDestination[i] += this->HamiltonianShift * vSource[i];
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      Coefficient = vSource[k];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      vDestination[k++] += this->HamiltonianShift * Coefficient;
	    }
	}
      else
	{
	  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  int Pos = firstComponent / this->FastMultiplicationStep; 
	  int PosMod = firstComponent % this->FastMultiplicationStep;
	  if (PosMod != 0)
	    {
	      ++Pos;
	      PosMod = this->FastMultiplicationStep - PosMod;
	    }
	  int l =  PosMod + firstComponent + this->PrecalculationShift;
	  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
	      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
	      Coefficient = vSource[l];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      vDestination[l] += this->HamiltonianShift * Coefficient;
	      l += this->FastMultiplicationStep;
	      ++Pos;
	    }
	  int Index;
	  int* MIndices;
	  int* NIndices;
	  double* TmpInteraction;
	  double Coefficient2;
	  firstComponent += this->PrecalculationShift;
	  LastComponent += this->PrecalculationShift;
	  for (l = 0; l < this->FastMultiplicationStep; ++l)
	    if (PosMod != l)
	      {	
		for (int k = 2; k <= this->MaxNBody; ++k)
		  if (this->NBodyFlags[k] == true)
		    {
		      for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			{
			  for (int j = this->MinSumIndices; j <= this->MaxSumIndices; ++j)
			    {
			      TmpInteraction = this->NBodyInteractionFactors[k][j];
			      int Lim = NbrSortedIndicesPerSum[k][j];
			      MIndices = this->SortedIndicesPerSum[k][j];
			      for (int i1 = 0; i1 < Lim; ++i1)
				{
				  Coefficient2 = vSource[i] * TmpInteraction[i1];
				  NIndices = this->SortedIndicesPerSum[k][j];
				  for (int i2 = 0; i2 < Lim; ++i2)
				    {
				      Index = TmpParticles->ProdAdProdA(i, MIndices, NIndices, k, Coefficient);
				      if (Index < Dim)
					vDestination[Index] += Coefficient * TmpInteraction[i2] * Coefficient2;
				      NIndices += k;
				    }
				  MIndices += k;
				}
			    }
			}		      
		    }
		for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
		  vDestination[i] += this->HamiltonianShift * vSource[i];					
	      }
	  delete TmpParticles;
	}
   }
  return vDestination;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long AbstractQHEOnSphereNBodyInteractionHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient;
  long Memory = 0;
  int* MIndices;
  int* NIndices;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int LastComponent = lastComponent + firstComponent;

/*  long NbrDifferent = 0;
  long Total = 0;
  long MaxNbrDifferent = 10000000;
  double* Different = new double [MaxNbrDifferent];*/

  for (int k = 2; k <= this->MaxNBody; ++k)
    if (this->NBodyFlags[k] == true)
      {
	for (int i = firstComponent; i < LastComponent; ++i)
	  {
	    for (int j = this->MinSumIndices; j <= this->MaxSumIndices; ++j)
	      {
		int Lim = NbrSortedIndicesPerSum[k][j];
		MIndices = this->SortedIndicesPerSum[k][j];
		for (int i1 = 0; i1 < Lim; ++i1)
		  {
		    NIndices = this->SortedIndicesPerSum[k][j];
		    for (int i2 = 0; i2 < Lim; ++i2)
		      {
			Index = TmpParticles->ProdAdProdA(i, MIndices, NIndices, k, Coefficient);
			if (Index < this->Particles->GetHilbertSpaceDimension())
			  {
			    ++Memory;
			    ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
/*			    ++Total;
			    if (NbrDifferent == 0)
			      {				
				Different[0] = Coefficient;
				++NbrDifferent;
			      }
			    else
			      {
				long MaxPos = NbrDifferent - 1;
				if (Coefficient <= Different[0])
				  {
				    if (Coefficient != Different[0])
				      {
					if (NbrDifferent < (MaxNbrDifferent - 1))
					  {
					    for (MaxPos = NbrDifferent - 1; MaxPos >= 0; --MaxPos)
					      Different[MaxPos + 1] = Different[MaxPos];				    
					    Different[0] = Coefficient;
					    ++NbrDifferent;
					  }
					else
					  {
					    MaxNbrDifferent *= 2;
					    double* TmpDifferent = new double [MaxNbrDifferent];
					    for (MaxPos = 0; MaxPos < NbrDifferent; ++MaxPos)
					      TmpDifferent[MaxPos + 1] = Different[MaxPos];				    
					    delete[] Different;
					    Different = TmpDifferent;
					    Different[0] = Coefficient;
					    ++NbrDifferent;
					  }
				      }
				  }
				else
				  if (Coefficient >= Different[MaxPos])
				    {
				      if (Coefficient != Different[MaxPos])
					{
					  if (NbrDifferent < (MaxNbrDifferent - 1))
					    {
					      Different[MaxPos + 1] = Coefficient;
					      ++NbrDifferent;
					    }
					  else
					    {
					      MaxNbrDifferent *= 2;
					      double* TmpDifferent = new double [MaxNbrDifferent];
					      for (MaxPos = 0; MaxPos < NbrDifferent; ++MaxPos)
						TmpDifferent[MaxPos] = Different[MaxPos];				    
					      delete[] Different;
					      Different = TmpDifferent;
					      Different[NbrDifferent] = Coefficient;
					      ++NbrDifferent;
					    }
					}
				    }
				  else
				    {
				      long MinPos = 0l;
				      long MedPos;
				      while ((MaxPos - MinPos) > 1)
					{
					  MedPos = (MaxPos + MinPos) >> 1;
					  if (Coefficient >= Different[MedPos])
					    {
					      MinPos = MedPos;
					    }
					  else
					    {
					      MaxPos = MedPos;
					    }
					}
				      if ((Coefficient != Different[MinPos]) && (Coefficient != Different[MaxPos]))
					{
					  ++NbrDifferent;
					  ++MinPos;
					  if (NbrDifferent < MaxNbrDifferent)
					    {
					      MaxPos = NbrDifferent - 1;
					      ++MinPos;
					      for (; MaxPos >= MinPos; --MaxPos)
						Different[MaxPos] = Different[MaxPos- 1];
					      Different[MinPos - 1] = Coefficient;
					    }
					  else
					    {
					      MaxNbrDifferent *= 2;
					      double* TmpDifferent = new double [MaxNbrDifferent];
					      for (MaxPos = 0; MaxPos < MinPos; ++MaxPos)
					      TmpDifferent[MaxPos] = Different[MaxPos];
					      TmpDifferent[MaxPos++] = Coefficient;
					      for (; MaxPos < NbrDifferent; ++MaxPos)
						TmpDifferent[MaxPos] = Different[MaxPos - 1];				    
					      delete[] Different;
					      Different = TmpDifferent;
					    }
					}
				    }
			      }*/
			  }
			NIndices += k;
		      }
		    MIndices += k;
		  }
	      }
	  }
      }
//  cout << "Total = " << Total << "   unique = " <<  NbrDifferent << "(" << (((double) NbrDifferent) * 100.0 / ((double) Total)) <<"%)" << endl;
  delete TmpParticles;

  return Memory;
}

// enable fast multiplication algorithm
//

void AbstractQHEOnSphereNBodyInteractionHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  int Index;
  double Coefficient;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int Pos;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficient = new double* [ReducedSpaceDimension];
  double* TmpInteraction;
  int* MIndices;
  int* NIndices;

  int TotalPos = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
    {
      this->InteractionPerComponentIndex[TotalPos] = new int [this->NbrInteractionPerComponent[TotalPos]];
      this->InteractionPerComponentCoefficient[TotalPos] = new double [this->NbrInteractionPerComponent[TotalPos]];      
      TmpIndexArray = this->InteractionPerComponentIndex[TotalPos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[TotalPos];
      Pos = 0;
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true)
	  {
	    for (int j = this->MinSumIndices; j <= this->MaxSumIndices; ++j)
	      {
		TmpInteraction = this->NBodyInteractionFactors[k][j];
		int Lim = this->NbrSortedIndicesPerSum[k][j];
		MIndices = this->SortedIndicesPerSum[k][j];
		for (int i1 = 0; i1 < Lim; ++i1)
		  {
		    NIndices = this->SortedIndicesPerSum[k][j];
		    for (int i2 = 0; i2 < Lim; ++i2)
		      {
			Index = this->Particles->ProdAdProdA(i, MIndices, NIndices, k, Coefficient);
			if (Index < this->Particles->GetHilbertSpaceDimension())
			  {
			    TmpIndexArray[Pos] = Index;
			    TmpCoefficientArray[Pos] = Coefficient * TmpInteraction[i1] *  TmpInteraction[i2];
			    ++Pos;
			  }
			NIndices += k;
		      }
		    MIndices += k;
		  }
	      }
	  }
      ++TotalPos;
    }
  this->FastMultiplicationFlag = true;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted

void AbstractQHEOnSphereNBodyInteractionHamiltonian::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int Pos;
  int Min = firstComponent / this->FastMultiplicationStep;
  int Max = lastComponent / this->FastMultiplicationStep;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  double* TmpInteraction;
  int* MIndices;
  int* NIndices;
 
  for (int i = Min; i < Max; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
      this->InteractionPerComponentCoefficient[i] = new double [this->NbrInteractionPerComponent[i]];      
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
      Pos = 0;
      for (int k = 2; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true)
	  {
	    for (int j = this->MinSumIndices; j <= this->MaxSumIndices; ++j)
	      {
		TmpInteraction = this->NBodyInteractionFactors[k][j];
		int Lim = NbrSortedIndicesPerSum[k][j];
		MIndices = this->SortedIndicesPerSum[k][j];
		for (int i1 = 0; i1 < Lim; ++i1)
		  {
		    NIndices = this->SortedIndicesPerSum[k][j];
		    for (int i2 = 0; i2 < Lim; ++i2)
		      {
			Index = TmpParticles->ProdAdProdA(i, MIndices, NIndices, k, Coefficient);
			if (Index < this->Particles->GetHilbertSpaceDimension())
			  {
			    TmpIndexArray[Pos] = Index;
			    TmpCoefficientArray[Pos] = Coefficient * TmpInteraction[i1] *  TmpInteraction[i2];
			    ++Pos;
			  }
			NIndices += k;
		      }
		    MIndices += k;
		  }
	      }
	  }
    }
  delete TmpParticles;
}

// get all indices needed to characterize a completly skew symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long  AbstractQHEOnSphereNBodyInteractionHamiltonian::GetAllSkewSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
										  int**& sortedIndicesPerSum)
{
  long** BinomialCoefficients = GetBinomialCoefficients(nbrValues);
  long NbrElements = BinomialCoefficients[nbrValues][nbrIndices];
  int** Indices = new int* [NbrElements];
  int* Sum = new int [NbrElements];
  int Min = nbrIndices - 1;
  int Max;
  int Step;
  int Pos = 0;
  for (int i = nbrValues - 1; i >= Min; --i)
    {
      Step = BinomialCoefficients[i][nbrIndices - 1];
      for (int j = 0; j < Step; ++j)
	{
	  Indices[Pos] = new int [nbrIndices];
	  Indices[Pos][0] = i;
	  Sum[Pos] = i;
	  ++Pos;
	}
    }
  for (int i = 1; i < nbrIndices; ++i)
    {
      int Pos = 0;
      Min = nbrIndices - i - 1;
      while (Pos < NbrElements)
	{
	  Max = Indices[Pos][i - 1] - 1;
	  for (; Max >= Min; --Max)
	    {
	      Step = BinomialCoefficients[Max][Min];
	      for (int j = 0; j < Step; ++j)
		{
		  Indices[Pos][i] = Max;
		  Sum[Pos] += Max;
		  ++Pos;
		}
	    }
	}
    }

  int MaxSum = (((nbrValues - 1) * nbrValues) - ((nbrIndices - 1) * (nbrIndices - 2)))/ 2;
  int MinSum = (nbrIndices * (nbrIndices - 1)) / 2;
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
  for (int i = 0; i <= MaxSum; ++i)
    nbrSortedIndicesPerSum[i] = 0;
  for (int i = 0; i < NbrElements; ++i)
    ++nbrSortedIndicesPerSum[Sum[i]];
  long* TmpPos = new long [MaxSum + 1];
  for (int i = MinSum; i <= MaxSum; ++i)
    {
      sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * nbrIndices];
      nbrSortedIndicesPerSum[i] = 0;
      TmpPos[i] = 0l;      
    }
  for (int i = 0; i < NbrElements; ++i)
    {   
      Pos = Sum[i];
      Max = nbrSortedIndicesPerSum[Pos];
      for (int j = 0; j < nbrIndices; ++j)
	{
	  sortedIndicesPerSum[Pos][TmpPos[Pos]] = Indices[i][j];
	  ++TmpPos[Pos];
	}
      ++nbrSortedIndicesPerSum[Pos];
      delete[] Indices[i];
    }
  delete[] TmpPos;
  delete[] Sum;
  delete[]Indices;
  for (int i = 0; i <= nbrValues; ++i)
    delete[] BinomialCoefficients[i];
  delete[] BinomialCoefficients;
  return NbrElements;
}

// get all indices needed to characterize a completly symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka inverse of the product of the factorial of the number 
//                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long AbstractQHEOnSphereNBodyInteractionHamiltonian::GetAllSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
									     int**& sortedIndicesPerSum,
									     double**& sortedIndicesPerSumSymmetryFactor)
{
  long** DimensionSymmetricGroup;
  if (nbrValues >= nbrIndices)
    DimensionSymmetricGroup = GetIrreducibleRepresentationDimensionSymmetricGroup(nbrValues);
  else
    DimensionSymmetricGroup = GetIrreducibleRepresentationDimensionSymmetricGroup(nbrIndices);
  long NbrElements = DimensionSymmetricGroup[nbrValues][nbrIndices];

  int** Indices = new int* [NbrElements];
  int* Sum = new int [NbrElements];
  int Max;
  int Step;
  int Pos = 0;
  int TmpNbrIndices;
  for (int i = nbrValues - 1; i >= 0; --i)
    {
      Step = DimensionSymmetricGroup[i + 1][nbrIndices - 1];
      for (int j = 0; j < Step; ++j)
	{
	  Indices[Pos] = new int [nbrIndices];
	  Indices[Pos][0] = i;
	  Sum[Pos] = i;
	  ++Pos;
	}
    }
  for (int i = 1; i < nbrIndices; ++i)
    {
      int Pos = 0;
      TmpNbrIndices = nbrIndices - i - 1;
      while (Pos < NbrElements)
	{
	  Max = Indices[Pos][i - 1];
	  Step = DimensionSymmetricGroup[Max + 1][TmpNbrIndices];
	  for (int j = 0; j < Step; ++j)
	    {
	      Indices[Pos][i] = Max;
	      Sum[Pos] += Max;
	      ++Pos;
	    }
	  --Max;
	  for (; Max >= 0; --Max)
	    {
	      Step = DimensionSymmetricGroup[Max + 1][TmpNbrIndices];
	      for (int j = 0; j < Step; ++j)
		{
		  Indices[Pos][i] = Max;
		  Sum[Pos] += Max;
		  ++Pos;
		}
	    }
	}
    }

  int MaxSum = (nbrValues - 1) * nbrIndices;
  long* TmpPos = new long [MaxSum + 1];
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
  sortedIndicesPerSumSymmetryFactor = new double* [MaxSum + 1];
  for (int i = 0; i <= MaxSum; ++i)
    nbrSortedIndicesPerSum[i] = 0;
  for (int i = 0; i < NbrElements; ++i)
    {
      ++nbrSortedIndicesPerSum[Sum[i]];
    }
  for (int i = 0; i <= MaxSum; ++i)
    {
      TmpPos[i] = 0l;
      sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * nbrIndices];
      sortedIndicesPerSumSymmetryFactor[i] = new double [nbrSortedIndicesPerSum[i]];
      nbrSortedIndicesPerSum[i] = 0;
    }
  for (int i = 0; i < NbrElements; ++i)
    {   
      Pos = Sum[i];
      Max = nbrSortedIndicesPerSum[Pos];
      int* TmpIndices = Indices[i];
      for (int j = 0; j < nbrIndices; ++j)
	{
	  sortedIndicesPerSum[Pos][TmpPos[Pos]] = TmpIndices[j];
	  ++TmpPos[Pos];
	}
      double& SymmetryFactor = sortedIndicesPerSumSymmetryFactor[Pos][Max];
      SymmetryFactor = 1.0;
      for (int j = 1; j < nbrIndices; ++j)
	{
	  int TmpSymmetryFactor = 1;
	  while ((j < nbrIndices) && (TmpIndices[j - 1] == TmpIndices[j]))
	    {
	      ++TmpSymmetryFactor;
	      ++j;
	    }
	  if (TmpSymmetryFactor != 1)
	    for (int k = 2; k <= TmpSymmetryFactor; ++k)
	      SymmetryFactor *= (double) k;
	}
      delete[] TmpIndices;
      SymmetryFactor = 1.0 / SymmetryFactor;
      ++nbrSortedIndicesPerSum[Pos];
    }

  delete[] TmpPos;
  delete[] Sum;
  delete[]Indices;
  for (int i = 0; i <= nbrValues; ++i)
    delete[] DimensionSymmetricGroup[i];
  delete[] DimensionSymmetricGroup;
  return NbrElements;
}

