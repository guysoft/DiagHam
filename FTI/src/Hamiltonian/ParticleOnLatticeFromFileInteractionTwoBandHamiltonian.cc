////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//          class of a two body interaction projected onto two bands          //
//         from an ASCII file providing the two body matrix elements          //
//                                                                            //
//                        last modification : 01/05/2020                      //
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
#include "Hamiltonian/ParticleOnLatticeFromFileInteractionTwoBandHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

ParticleOnLatticeFromFileInteractionTwoBandHamiltonian::ParticleOnLatticeFromFileInteractionTwoBandHamiltonian()
{
}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// matrixElementsInteractionFile = name of the ASCII file containing the matrix element for the generic two body interaction term
// tightBindingModel = pointer to the tight binding model
// flatBandFlag = use flat band model
// flatBandOneBodyGap = set the gap between the first band and the second band when using the flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeFromFileInteractionTwoBandHamiltonian::ParticleOnLatticeFromFileInteractionTwoBandHamiltonian(ParticleOnSphereWithSpin* particles, int nbrParticles,
													       int nbrSiteX, int nbrSiteY,
													       char* matrixElementsInteractionFile,
													       Abstract2DTightBindingModel* tightBindingModel, 
													       bool flatBandFlag, double flatBandOneBodyGap,
													       AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->HamiltonianShift = 0.0;
  this->TightBindingModel = tightBindingModel;
  this->FlatBand = flatBandFlag;
  this->FlatBandOneBodyGap = flatBandOneBodyGap;
  this->MatrixElementsInteractionFile = new char[strlen(matrixElementsInteractionFile) + 1];
  strcpy(this->MatrixElementsInteractionFile, matrixElementsInteractionFile);
  
  this->Architecture = architecture;
  this->Memory = memory;
  this->OneBodyInteractionFactorsupup = 0;
  this->OneBodyInteractionFactorsdowndown = 0;
  this->OneBodyInteractionFactorsupdown = 0;
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      if (TmpMemory < 1024)
	cout  << "fast = " <<  TmpMemory << "b ";
      else
	if (TmpMemory < (1 << 20))
	  cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	else
	  if (TmpMemory < (1 << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    {
	      cout  << "fast = " << (TmpMemory >> 30) << ".";
	      TmpMemory -= ((TmpMemory >> 30) << 30);
	      TmpMemory *= 100l;
	      TmpMemory >>= 30;
	      if (TmpMemory < 10l)
		cout << "0";
	      cout  << TmpMemory << " Gb ";
	    }
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeFromFileInteractionTwoBandHamiltonian::~ParticleOnLatticeFromFileInteractionTwoBandHamiltonian()
{
}
  
// evaluate all interaction factors
//   

void ParticleOnLatticeFromFileInteractionTwoBandHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0;
  ComplexMatrix* OneBodyBasis = new ComplexMatrix [this->TightBindingModel->GetNbrStatePerBand()];
  this->InteractionFactorsupup = 0;
  this->InteractionFactorsdowndown = 0;
  this->InteractionFactorsupdown = 0;
  this->InteractionFactorsupupupup = 0;
  this->InteractionFactorsupupdowndown = 0;
  this->InteractionFactorsdowndownupup = 0;
  this->InteractionFactorsdowndowndowndown = 0;
  this->InteractionFactorsupdownupup = 0;
  this->InteractionFactorsupdowndowndown = 0;
  this->InteractionFactorsupupupdown = 0;
  this->InteractionFactorsdowndownupdown = 0;
  this->InteractionFactorsupdownupdown = 0;

  MultiColumnASCIIFile TmpInteractionFile;
  if (TmpInteractionFile.Parse(this->MatrixElementsInteractionFile) == false)
    {
      TmpInteractionFile.DumpErrors(cout) << endl;
      exit(0);
    }
  if (TmpInteractionFile.GetNbrLines() == 0)
    {
      cout << this->MatrixElementsInteractionFile << " is an empty file" << endl;
      exit(0);
    }
  if (TmpInteractionFile.GetNbrColumns() < 13)
    {
      cout << this->MatrixElementsInteractionFile << " has a wrong number of column (has "
	   << TmpInteractionFile.GetNbrColumns() << ", should be at least 13)" << endl;
      exit(0);
    }
  int TmpNbrTwoBodyMatrixElements = TmpInteractionFile.GetNbrLines();
  cout << "nbr of two body matrix elements in " << this->MatrixElementsInteractionFile << " = " << TmpNbrTwoBodyMatrixElements << endl;

  int* TmpSigma1 = TmpInteractionFile.GetAsIntegerArray(0);
  int* TmpSigma2 = TmpInteractionFile.GetAsIntegerArray(3);
  int* TmpSigma3 = TmpInteractionFile.GetAsIntegerArray(6);
  int* TmpSigma4 = TmpInteractionFile.GetAsIntegerArray(9);
  int* TmpKx1 = TmpInteractionFile.GetAsIntegerArray(1);
  int* TmpKx2 = TmpInteractionFile.GetAsIntegerArray(4);
  int* TmpKx3 = TmpInteractionFile.GetAsIntegerArray(7);
  int* TmpKx4 = TmpInteractionFile.GetAsIntegerArray(10);
  int* TmpKy1 = TmpInteractionFile.GetAsIntegerArray(2);
  int* TmpKy2 = TmpInteractionFile.GetAsIntegerArray(5);
  int* TmpKy3 = TmpInteractionFile.GetAsIntegerArray(8);
  int* TmpKy4 = TmpInteractionFile.GetAsIntegerArray(11);
  double* TmpMatrixElements = TmpInteractionFile.GetAsDoubleArray(12);
  if (TmpMatrixElements == 0)
    {
      TmpInteractionFile.DumpErrors(cout) << endl;
      exit(0);
    }
  
  int* TmpLinearizedSumK = new int[TmpNbrTwoBodyMatrixElements];
  int* TmpLinearizedK1 = new int[TmpNbrTwoBodyMatrixElements];
  int* TmpLinearizedK2 = new int[TmpNbrTwoBodyMatrixElements];
  int* TmpLinearizedK3 = new int[TmpNbrTwoBodyMatrixElements];
  int* TmpLinearizedK4 = new int[TmpNbrTwoBodyMatrixElements];
  for (int i = 0; i < TmpNbrTwoBodyMatrixElements; ++i)
    {
      TmpLinearizedSumK[i] = this->TightBindingModel->GetLinearizedMomentumIndex((TmpKx1[i] + TmpKx2[i]) % this->NbrSiteX,
										 (TmpKy1[i] + TmpKy2[i]) % this->NbrSiteY);
      TmpLinearizedK1[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpKx1[i], TmpKy1[i]);
      TmpLinearizedK2[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpKx2[i], TmpKy2[i]);
      TmpLinearizedK3[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpKx3[i], TmpKy3[i]);
      TmpLinearizedK4[i] = this->TightBindingModel->GetLinearizedMomentumIndex(TmpKx4[i], TmpKy4[i]);
    }

  this->BandIndex1 = 0;
  this->BandIndex2 = 1;
  int* SigmaToBand = new int[2];
  SigmaToBand[0] = this->BandIndex1;
  SigmaToBand[1] = this->BandIndex2;

  if ((this->FlatBand == false) || (this->FlatBandOneBodyGap != 0.0))
    {
      this->OneBodyInteractionFactorsupup = new double [this->TightBindingModel->GetNbrStatePerBand()];
      this->OneBodyInteractionFactorsdowndown = new double [this->TightBindingModel->GetNbrStatePerBand()];
    }
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	int Index = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
	if (this->FlatBand == false)
	  {
	    this->OneBodyInteractionFactorsupup[Index] = this->TightBindingModel->GetEnergy(this->BandIndex1, Index);
	    this->OneBodyInteractionFactorsdowndown[Index] = this->TightBindingModel->GetEnergy(this->BandIndex2, Index);
	  }
	else
	  {
	    if (this->FlatBandOneBodyGap != 0.0)
	      {
		this->OneBodyInteractionFactorsupup[Index] = 0.0;
		this->OneBodyInteractionFactorsdowndown[Index] = this->FlatBandOneBodyGap;		
	      }
	  }
	OneBodyBasis[Index] =  this->TightBindingModel->GetOneBodyMatrix(Index);
      }

  int NbrInternalIndices = 2;
  this->InteractionFactorsSigma = new Complex***** [NbrInternalIndices];
  for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
    {
      this->InteractionFactorsSigma[sigma3] = new Complex****  [NbrInternalIndices];
      for (int sigma4 = sigma3; sigma4 < NbrInternalIndices; ++sigma4)
	{
	  this->InteractionFactorsSigma[sigma3][sigma4] = new Complex***[NbrInternalIndices];
	  for (int sigma1 = 0; sigma1 < NbrInternalIndices; ++sigma1)
	    {
	      this->InteractionFactorsSigma[sigma3][sigma4][sigma1] = new Complex**[NbrInternalIndices];
	    }
	}
    }

 
  this->NbrInterSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrInterSectorIndicesPerSum = new int[this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    this->NbrInterSectorIndicesPerSum[i] = 0;
  this->NbrIntraSectorSums = this->NbrSiteX * this->NbrSiteY;
  this->NbrIntraSectorIndicesPerSum = new int[this->NbrIntraSectorSums];
  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
    this->NbrIntraSectorIndicesPerSum[i] = 0;      

  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)
	  {
	    ++this->NbrInterSectorIndicesPerSum[this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, (ky1 + ky2) % this->NbrSiteY)];
	  }
  this->InterSectorIndicesPerSum = new int* [this->NbrInterSectorSums];
  for (int i = 0; i < this->NbrInterSectorSums; ++i)
    {
      if (this->NbrInterSectorIndicesPerSum[i] > 0)
	{
	  this->InterSectorIndicesPerSum[i] = new int[2 * this->NbrInterSectorIndicesPerSum[i]];      
	  this->NbrInterSectorIndicesPerSum[i] = 0;
	}
    }

  int*** TmpLinearizedKInterIndices = new int** [this->NbrInterSectorSums];
  for (int j = 0; j < this->NbrInterSectorSums; ++j)
    {
      TmpLinearizedKInterIndices[j] = new int* [this->NbrSiteX * this->NbrSiteY];
      for (int i = 0; i < (this->NbrSiteX * this->NbrSiteY); ++i)
	{
	  TmpLinearizedKInterIndices[j][i] = new int [this->NbrSiteX * this->NbrSiteY];
	}
    }
  int*** TmpLinearizedKIntraIndices = new int** [this->NbrIntraSectorSums];
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      TmpLinearizedKIntraIndices[j] = new int* [this->NbrSiteX * this->NbrSiteY];
      for (int i = 0; i < (this->NbrSiteX * this->NbrSiteY); ++i)
	{
	  TmpLinearizedKIntraIndices[j][i] = new int [this->NbrSiteX * this->NbrSiteY];
	}
    }
  
  for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
    for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
      for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2)    
	  {
	    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, (ky1 + ky2) % this->NbrSiteY);
	    this->InterSectorIndicesPerSum[TmpSum][this->NbrInterSectorIndicesPerSum[TmpSum] << 1] = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
	    this->InterSectorIndicesPerSum[TmpSum][1 + (this->NbrInterSectorIndicesPerSum[TmpSum] << 1)] = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
	    TmpLinearizedKInterIndices[TmpSum][this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1)][this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2)] = this->NbrInterSectorIndicesPerSum[TmpSum];
	    ++this->NbrInterSectorIndicesPerSum[TmpSum];    
	  }
 
  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 < Index2)
		  {
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, (ky1 + ky2) % this->NbrSiteY);
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      this->IntraSectorIndicesPerSum = new int* [this->NbrIntraSectorSums];
      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	{
	  if (this->NbrIntraSectorIndicesPerSum[i]  > 0)
	    {
	      this->IntraSectorIndicesPerSum[i] = new int[2 * this->NbrIntraSectorIndicesPerSum[i]];      
	      this->NbrIntraSectorIndicesPerSum[i] = 0;
	    }
	}
      for (int kx1 = 0; kx1 < this->NbrSiteX; ++kx1)
	for (int kx2 = 0; kx2 < this->NbrSiteX; ++kx2)
	  for (int ky1 = 0; ky1 < this->NbrSiteY; ++ky1)
	    for (int ky2 = 0; ky2 < this->NbrSiteY; ++ky2) 
	      {
		int Index1 = this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1);
		int Index2 = this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2);
		if (Index1 < Index2)
		  {
		    int TmpSum = this->TightBindingModel->GetLinearizedMomentumIndex((kx1 + kx2) % this->NbrSiteX, 
										     (ky1 + ky2) % this->NbrSiteY);
		    this->IntraSectorIndicesPerSum[TmpSum][this->NbrIntraSectorIndicesPerSum[TmpSum] << 1] = Index1;
		    this->IntraSectorIndicesPerSum[TmpSum][1 + (this->NbrIntraSectorIndicesPerSum[TmpSum] << 1)] = Index2;
		    TmpLinearizedKIntraIndices[TmpSum][this->TightBindingModel->GetLinearizedMomentumIndex(kx1, ky1)][this->TightBindingModel->GetLinearizedMomentumIndex(kx2, ky2)] = this->NbrIntraSectorIndicesPerSum[TmpSum];
		    ++this->NbrIntraSectorIndicesPerSum[TmpSum];    
		  }
	      }
      

      double TmpKx1 = 0.0;
      double TmpKx2 = 0.0;
      double TmpKx3 = 0.0;
      double TmpKx4 = 0.0;
      double TmpKy1 = 0.0;
      double TmpKy2 = 0.0;
      double TmpKy3 = 0.0;
      double TmpKy4 = 0.0;
      
      for (int sigma1 = 0; sigma1 < NbrInternalIndices; ++sigma1)
	{
	  for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
	    {
	      this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1] = new Complex*[this->NbrIntraSectorSums];
	      for (int j = 0; j < this->NbrIntraSectorSums; ++j)
		{
		  int Tmp = this->NbrIntraSectorIndicesPerSum[j] * this->NbrIntraSectorIndicesPerSum[j];
		  this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma1][j] = new Complex [Tmp];
		}
	    }
	}

      for (int sigma1 = 0; sigma1 < NbrInternalIndices; ++sigma1)
	{
	  for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
	    {
	      for (int sigma4 = sigma3 + 1; sigma4 < NbrInternalIndices; ++sigma4)
		{
		  this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1] = new Complex*[this->NbrIntraSectorSums];
		  for (int j = 0; j < this->NbrInterSectorSums; ++j)
		    {
		      this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma1][j] = new Complex [this->NbrIntraSectorIndicesPerSum[j] * this->NbrInterSectorIndicesPerSum[j]];
		    }
		}
	    }
	}

      for (int sigma1 = 0; sigma1 < NbrInternalIndices; ++sigma1)
	{
	  for (int sigma2 = sigma1 + 1; sigma2 < NbrInternalIndices; ++sigma2)
	    {
	      for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
		{
		  this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2] = new Complex*[this->NbrIntraSectorSums];
		  for (int j = 0; j < this->NbrInterSectorSums; ++j)
		    {
		      this->InteractionFactorsSigma[sigma3][sigma3][sigma1][sigma2][j] = new Complex [this->NbrInterSectorIndicesPerSum[j] * this->NbrIntraSectorIndicesPerSum[j]];

		    }
		}
	    }
	}

      for (int sigma1 = 0; sigma1 < NbrInternalIndices; ++sigma1)
	{
	  for (int sigma2 = sigma1 + 1; sigma2 < NbrInternalIndices; ++sigma2)
	    {
	      for (int sigma3 = 0; sigma3 < NbrInternalIndices; ++sigma3)
		{
		  for (int sigma4 = sigma3 + 1; sigma4 < NbrInternalIndices; ++sigma4)
		    {
		      this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2] = new Complex*[this->NbrIntraSectorSums];
		      for (int j = 0; j < this->NbrInterSectorSums; ++j)
			{
			  this->InteractionFactorsSigma[sigma3][sigma4][sigma1][sigma2][j] = new Complex [this->NbrInterSectorIndicesPerSum[j] * this->NbrInterSectorIndicesPerSum[j]];
			}
		    }
		}
	    }
	}
      
      for (int i = 0; i < TmpNbrTwoBodyMatrixElements; ++i)
	{
	  int TmpIndex = 0;
	  int TmpMaxIndexFactor = 0;
	  int TmpSumK = TmpLinearizedSumK[i];
	  int Sigma1 = TmpSigma1[i];
	  int Sigma2 = TmpSigma2[i];
	  int Sigma3 = TmpSigma3[i];
	  int Sigma4 = TmpSigma4[i];
	  int K1 = TmpLinearizedK1[i];
	  int K2 = TmpLinearizedK2[i];
	  int K3 = TmpLinearizedK3[i];
	  int K4 = TmpLinearizedK4[i];
	  double TmpSign = 1.0;
	  if (Sigma2 < Sigma1)
	    {
	      int Tmp = Sigma2;
	      Sigma2 = Sigma1;
	      Sigma1 = Tmp;
	      Tmp = K2;
	      K2 = K1;
	      K1 = Tmp;
	      TmpSign *= -1.0;
	    }
	  if (Sigma4 < Sigma3)
	    {
	      int Tmp = Sigma4;
	      Sigma4 = Sigma3;
	      Sigma3 = Tmp;
	      Tmp = K4;
	      K4 = K3;
	      K3 = Tmp;
	      TmpSign *= -1.0;
	    }
	  if (Sigma1 == Sigma2)
	    {
	      if (K1 > K2)
		{
		  int Tmp = K2;
		  K2 = K1;
		  K1 = Tmp;
		  TmpSign *= -1.0;
		}
	      if (K1 == K2)
		{
		  TmpSign = 0.0;
		}
	      else
		{
		  if (Sigma3 == Sigma4)
		    {
		      if (K3 > K4)
			{
			  int Tmp = K4;
			  K4 = K3;
			  K3 = Tmp;
			  TmpSign *= -1.0;
			}
		      if (K3 == K4)
			{
			  TmpSign = 0.0;
			}
		      else
			{
			  TmpIndex = ((this->NbrIntraSectorIndicesPerSum[TmpSumK] * TmpLinearizedKIntraIndices[TmpSumK][K3][K4])
				      + TmpLinearizedKIntraIndices[TmpSumK][K1][K2]);
			}
		    }
		  else
		    {
		      TmpIndex = ((this->NbrIntraSectorIndicesPerSum[TmpSumK] * TmpLinearizedKInterIndices[TmpSumK][K3][K4])
				  + TmpLinearizedKIntraIndices[TmpSumK][K1][K2]);
		    }
		}
	    }
	  else
	    {
	      if (Sigma3 == Sigma4)
		{
		  if (K3 > K4)
		    {
		      int Tmp = K4;
		      K4 = K3;
		      K3 = Tmp;
		      TmpSign *= -1.0;
		    }
		  if (K3 == K4)
		    {
		      TmpSign = 0.0;
		    }
		  else
		    {
		      TmpIndex = ((this->NbrInterSectorIndicesPerSum[TmpSumK] * TmpLinearizedKIntraIndices[TmpSumK][K3][K4])
				  + TmpLinearizedKInterIndices[TmpSumK][K1][K2]);
		    }
		}
	      else
		{
		  TmpIndex = ((this->NbrInterSectorIndicesPerSum[TmpSumK] * TmpLinearizedKInterIndices[TmpSumK][K3][K4])
			      + TmpLinearizedKInterIndices[TmpSumK][K1][K2]);
		}
	    }
	  if (TmpSign != 0.0)
	    {
	      this->InteractionFactorsSigma[Sigma1][Sigma2][Sigma3][Sigma4][TmpSumK][TmpIndex] += TmpSign * TmpMatrixElements[i];
	    }
	}
    }
  else
    {      
      // bosonic interaction
    }

  delete[] TmpSigma1;
  delete[] TmpSigma2;
  delete[] TmpSigma3;
  delete[] TmpSigma4;
  delete[] TmpKx1;
  delete[] TmpKx2;
  delete[] TmpKx3;
  delete[] TmpKx4;
  delete[] TmpKy1;
  delete[] TmpKy2;
  delete[] TmpKy3;
  delete[] TmpKy4;
  delete[] TmpMatrixElements;
  delete[] TmpLinearizedSumK;
  delete[] TmpLinearizedK1;
  delete[] TmpLinearizedK2;
  delete[] TmpLinearizedK3;
  delete[] TmpLinearizedK4;
  for (int j = 0; j < this->NbrInterSectorSums; ++j)
    {
      for (int i = 0; i < (this->NbrSiteX * this->NbrSiteY); ++i)
	{
	  delete[] TmpLinearizedKInterIndices[j][i];
	}
      delete[] TmpLinearizedKInterIndices[j];
    }
  delete[] TmpLinearizedKInterIndices;
  for (int j = 0; j < this->NbrIntraSectorSums; ++j)
    {
      for (int i = 0; i < (this->NbrSiteX * this->NbrSiteY); ++i)
	{
	  delete[] TmpLinearizedKIntraIndices[j][i];
	}
      delete[] TmpLinearizedKIntraIndices[j];
    }
  delete[] TmpLinearizedKIntraIndices;
  delete[] OneBodyBasis;
  delete[] SigmaToBand;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}

