////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                 Copyright (C) 2001-2004 Antoine Sterdyniak                 //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//          the 3-body interaction that produces the Haffnian state           //
//                                                                            //
//                        last modification : 23/07/2010                      //
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
#include "Hamiltonian/ParticleOnTorusHaffnianHamiltonian.h"
#include "Architecture/AbstractArchitecture.h"
#include "MathTools/ClebschGordanCoefficients.h"
  
#include <stdio.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


// default constructor
//

ParticleOnTorusHaffnianHamiltonian::ParticleOnTorusHaffnianHamiltonian()
{
}

// constructor from datas with three-body interaction and optional two-body and one-body
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusHaffnianHamiltonian::ParticleOnTorusHaffnianHamiltonian(ParticleOnTorus* particles, int nbrParticles, int lzmax, double ratio,
										       AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
										       char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = lzmax;
  this->NbrLzValue = this->LzMax + 1;
  this->NbrParticles = nbrParticles;
  this->MaxNBody = 3;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->NBodyFlags = new bool [this->MaxNBody + 1];
  this->NBodyInteractionFactors = new double** [this->MaxNBody + 1];
  this->NbrSortedIndicesPerSum = new int* [this->MaxNBody + 1];
  this->SortedIndicesPerSum = new int** [this->MaxNBody + 1];
  this->MinSumIndices = new int [this->MaxNBody + 1];
  this->MaxSumIndices = new int [this->MaxNBody + 1];
  this->NBodySign = new double[this->MaxNBody + 1];
  this->NbrNIndices = new long[this->MaxNBody + 1];
  this->NIndices = new int*[this->MaxNBody + 1];
  this->NbrMIndices = new long*[this->MaxNBody + 1];
  this->MIndices = new int**[this->MaxNBody + 1];
  this->MNNBodyInteractionFactors = new double**[this->MaxNBody + 1];
  this->OneBodyTermFlag = false;
  this->FullTwoBodyFlag = false;
  for (int k = 0; k <= this->MaxNBody; ++k)
    {
      this->MinSumIndices[k] = 1;
      this->MaxSumIndices[k] = 0;      
      this->NBodyFlags[k] = false;
      this->NBodySign[k] = 1.0;
      if ((this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic) && ((k & 1) == 0))
	this->NBodySign[k] = -1.0;
    }
  this->NBodyFlags[3] = true;
  this->Architecture = architecture;
  this->EvaluateInteractionFactors();
  this->HamiltonianShift = 0.0;
  this->L2Operator = 0;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->DiskStorageFlag = onDiskCacheFlag;
  this->Memory = memory;
  if (precalculationFileName == 0)
    {
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
		cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
	  if (this->DiskStorageFlag == false)
	    {
	      this->EnableFastMultiplication();
	    }
	  else
	    {
	      char* TmpFileName = this->Architecture->GetTemporaryFileName();
	      this->EnableFastMultiplicationWithDiskStorage(TmpFileName);	      
	      delete[] TmpFileName;
	    }
	}
      else
	{
	  this->FastMultiplicationFlag = false;
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnTorusHaffnianHamiltonian::~ParticleOnTorusHaffnianHamiltonian()
{
  for (int k = 1; k <= this->MaxNBody; ++k)
    if (this->NBodyFlags[k] == true)
      {
	for (int MinSum = this->MinSumIndices[k]; MinSum <= this->MaxSumIndices[k]; ++MinSum)
	  {
	    delete[] this->SortedIndicesPerSum[k][MinSum];
	    if (this->MNNBodyInteractionFactors == 0)
	      delete[] this->NBodyInteractionFactors[k][MinSum];	      
	  }
	delete[] this->NbrSortedIndicesPerSum[k];
	delete[] this->SortedIndicesPerSum[k];
	if (this->MNNBodyInteractionFactors == 0)
	  delete[] this->NBodyInteractionFactors[k];
	else
	  {
	    for (int i = 0; i < this->NbrNIndices[k]; ++i)
	      {
		delete[] this->MNNBodyInteractionFactors[k][i];		
		delete[] this->MIndices[k][i];
	      }
	    delete[] this->NIndices[k];
	    delete[] this->NbrMIndices[k];
	    delete[] this->MIndices[k];
	    delete[] this->MNNBodyInteractionFactors[k];
	  }
      }
  
  delete[] this->NbrNIndices;
  delete[] this->NIndices;
  delete[] this->NbrMIndices;
  delete[] this->MIndices;
  delete[] this->MNNBodyInteractionFactors;
  
  delete[] this->NBodyFlags;
  delete[] this->NBodyInteractionFactors;
  delete[] this->SortedIndicesPerSum;
  delete[] this->NbrSortedIndicesPerSum;
  delete[] this->MinSumIndices;
  delete[] this->MaxSumIndices;
  delete[] this->NBodySign;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* ParticleOnTorusHaffnianHamiltonian::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusHaffnianHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  if (this->FastMultiplicationFlag == true)
    {
      for (int i = 0; i < this->Particles->GetHilbertSpaceDimension(); ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
    }
  this->Particles = (ParticleOnTorus*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusHaffnianHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}

// evaluate all interaction factors
//   

void ParticleOnTorusHaffnianHamiltonian::EvaluateInteractionFactors()
{
  if (this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
    {
    }
  else
    {
      this->MinSumIndices[3] = 0;
      this->MaxSumIndices[3] = this->LzMax;
      double** SortedIndicesPerSumSymmetryFactor;
      GetAllSymmetricIndices(this->NbrLzValue, 3, this->NbrSortedIndicesPerSum[3], this->SortedIndicesPerSum[3],
			     SortedIndicesPerSumSymmetryFactor);
      
      
      long TmpNbrNIndices = 0;
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[3]; ++MinSum)
	TmpNbrNIndices += this->NbrSortedIndicesPerSum[3][MinSum];
      cout <<"TmpNbrNIndices = "<<TmpNbrNIndices<<endl;
      this->NbrNIndices[3] = TmpNbrNIndices;
      this->NIndices[3] = new int[TmpNbrNIndices * 3];
      this->NbrMIndices[3] = new long[TmpNbrNIndices];
      this->MIndices[3] = new int*[TmpNbrNIndices];
      this->MNNBodyInteractionFactors[3] = new double* [TmpNbrNIndices];
      TmpNbrNIndices = 0;
      int* TmpNIndices = this->NIndices[3];
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[3]; ++MinSum)
	{
	  int Lim = this->NbrSortedIndicesPerSum[3][MinSum];
	  double* TmpSymmetryFactors = SortedIndicesPerSumSymmetryFactor[MinSum];
	  int* TmpNIndices2 = this->SortedIndicesPerSum[3][MinSum];
	  
	  for (int i = 0; i < Lim; ++i)
	    {
	      this->NbrMIndices[3][TmpNbrNIndices] = Lim;		    
	      this->MIndices[3][TmpNbrNIndices] = new int [Lim * 3];
	      this->MNNBodyInteractionFactors[3][TmpNbrNIndices] = new double [Lim];
	      int* TmpMIndices = this->MIndices[3][TmpNbrNIndices];
	      int* TmpMIndices2 = this->SortedIndicesPerSum[3][MinSum];
	      double* TmpInteraction = this->MNNBodyInteractionFactors[3][TmpNbrNIndices];
	      for (int j = 0; j < Lim; ++j)
		{
		  double& TmpInteraction2 = TmpInteraction[j];
		  
		  TmpInteraction2 = 0.0;
		  if (TmpMIndices2[0] != TmpMIndices2[1])
		    {
		      if (TmpMIndices2[1] != TmpMIndices2[2])
			{
			  TmpInteraction2 += this->EvaluateSymmetricPermutationNIndicesInteractionFactor(TmpMIndices2[0],TmpMIndices2[1],TmpMIndices2[2],TmpNIndices2[0],TmpNIndices2[1],TmpNIndices2[2]);
			  TmpInteraction2 += this->EvaluateSymmetricPermutationNIndicesInteractionFactor(TmpMIndices2[0],TmpMIndices2[2],TmpMIndices2[1],TmpNIndices2[0],TmpNIndices2[1],TmpNIndices2[2]);
			  TmpInteraction2 += this->EvaluateSymmetricPermutationNIndicesInteractionFactor(TmpMIndices2[1],TmpMIndices2[0],TmpMIndices2[2],TmpNIndices2[0],TmpNIndices2[1],TmpNIndices2[2]);
			  TmpInteraction2 += this->EvaluateSymmetricPermutationNIndicesInteractionFactor(TmpMIndices2[1],TmpMIndices2[2],TmpMIndices2[0],TmpNIndices2[0],TmpNIndices2[1],TmpNIndices2[2]);
			  TmpInteraction2 += this->EvaluateSymmetricPermutationNIndicesInteractionFactor(TmpMIndices2[2],TmpMIndices2[0],TmpMIndices2[1],TmpNIndices2[0],TmpNIndices2[1],TmpNIndices2[2]);
			  TmpInteraction2 += this->EvaluateSymmetricPermutationNIndicesInteractionFactor(TmpMIndices2[2],TmpMIndices2[1],TmpMIndices2[0],TmpNIndices2[0],TmpNIndices2[1],TmpNIndices2[2]);
			}
		      else
			{
			  TmpInteraction2 += this->EvaluateSymmetricPermutationNIndicesInteractionFactor(TmpMIndices2[0],TmpMIndices2[1],TmpMIndices2[2],TmpNIndices2[0],TmpNIndices2[1],TmpNIndices2[2]);
			  TmpInteraction2 += this->EvaluateSymmetricPermutationNIndicesInteractionFactor(TmpMIndices2[1],TmpMIndices2[0],TmpMIndices2[2],TmpNIndices2[0],TmpNIndices2[1],TmpNIndices2[2]);
			  TmpInteraction2 += this->EvaluateSymmetricPermutationNIndicesInteractionFactor(TmpMIndices2[1],TmpMIndices2[2],TmpMIndices2[0],TmpNIndices2[0],TmpNIndices2[1],TmpNIndices2[2]);
			}
		    }
		  else
		    {
		      if (TmpMIndices2[1] != TmpMIndices2[2])
			{
			  TmpInteraction2 += this->EvaluateSymmetricPermutationNIndicesInteractionFactor(TmpMIndices2[0],TmpMIndices2[1],TmpMIndices2[2],TmpNIndices2[0],TmpNIndices2[1],TmpNIndices2[2]);
			  TmpInteraction2 += this->EvaluateSymmetricPermutationNIndicesInteractionFactor(TmpMIndices2[0],TmpMIndices2[2],TmpMIndices2[1],TmpNIndices2[0],TmpNIndices2[1],TmpNIndices2[2]);
			  TmpInteraction2 += this->EvaluateSymmetricPermutationNIndicesInteractionFactor(TmpMIndices2[2],TmpMIndices2[0],TmpMIndices2[1],TmpNIndices2[0],TmpNIndices2[1],TmpNIndices2[2]);
			}
		      else
			{
			  TmpInteraction2 += this->EvaluateSymmetricPermutationNIndicesInteractionFactor(TmpMIndices2[0],TmpMIndices2[1],TmpMIndices2[2],TmpNIndices2[0],TmpNIndices2[1],TmpNIndices2[2]);
			}
		    }
		  cout << TmpMIndices2[0] << " " << TmpMIndices2[1] << " " << TmpMIndices2[2] << " " << TmpNIndices2[0] << " " << TmpNIndices2[1] << " " << TmpNIndices2[2] << " : " << TmpInteraction2 << endl;
		  for (int l = 0; l < 3; ++l)
		    {
		      (*TmpMIndices) = (*TmpMIndices2);			
		      ++TmpMIndices;
		      ++TmpMIndices2;
		    }
		}
	      for (int j = 0; j < 3; ++j)
		{
		  (*TmpNIndices) = (*TmpNIndices2);			
		  ++TmpNIndices;
		  ++TmpNIndices2;
		}
	      ++TmpNbrNIndices;
	    }
	}
      for (int MinSum = 0; MinSum <= this->MaxSumIndices[3]; ++MinSum)
	{
	  delete[] SortedIndicesPerSumSymmetryFactor[MinSum];
	}
      delete[] SortedIndicesPerSumSymmetryFactor;
    }
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a+_m3 a_n1 a_n2 a_n3 coupling term
//
// m1 = first creation operator index
// m2 = second creation operator index
// m3 = third creation operator index
// n1 = first annihilation operator index
// n2 = second annihilation operator index
// n3 = thrid annihilation operator index
//
// return value = numerical coefficient

double ParticleOnTorusHaffnianHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int n1,int n2, int n3)
{
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  double PIOnM = M_PI / DoubleNbrLzValue ;
  double Q2Factor = 2.0 * PIOnM / DoubleNbrLzValue;//4.0 * M_PI * PIOnM;
  double ResultNy2 = 0.0;
  double ResultNy1 = 0.0;
  double ResultNx2 = 0.0;
  double ResultNx1 = 0.0;
  double Nx1;
  double Nx2;
  double Q1;
  double Q2;
  double Q3;
  double Q4;
  double Q5;
  double Q21;
  double Q22;
  double Q23;
  double IniNy1 = (double) (m1 - n3);
  double Ny1 = IniNy1;
  double IniNy2 = (double) (n1 - m3);
  double Ny2 = IniNy2;
  double PremFactor1 = ((2.0 * ((double) (m1 - n2))) - Ny2)* PIOnM;
  double PremFactor2 = ((2.0 * ((double) (n2 - m3))) - Ny1) * PIOnM;
  double Factor1 = PremFactor1;
  double Factor2 = PremFactor2;
  double Precision = 1.0;
  double Precision1 = 1.0;
  double Coefficient = 1.0;
  double Coefficient1 = 1.0;
  double Coefficient2 = 1.0;
  double Tmp;
  //  cout << "coef " << m1 << " "  << m2 << " "  << m3 << " "  << n1 << " "<< n2 << " "<< n3 << " : ";
  
  while ((fabs(ResultNy2) + fabs(Coefficient)) != fabs(ResultNy2))
    {
      Q1 = this->Ratio * Ny2 * Ny2;
      if (Ny2 != 0.0)
	Coefficient = exp(- PIOnM * Q1);
      else
	Coefficient = 1.0;
      ResultNy1 = 0.0;
      Ny1 = IniNy1;
      Factor2 = PremFactor2;
      Coefficient1 = 1.0;
      Nx1 = 0;
      Nx2 = 0;
      while ((fabs(ResultNy1) + fabs(Coefficient1)) != fabs(ResultNy1))
	{
	  Q2 = this->Ratio * Ny1 * (Ny1 - Ny2);
	  Q21 = (this->Ratio * Ny1 * Ny1);
	  Q22 = (this->Ratio * Ny2 * Ny2);
	  Q23 = (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
	  Q21 *= Q2Factor;
	  Q22 *= Q2Factor;
	  Q23 *= Q2Factor;
	  if ((Ny1 == 0.0)||(Ny1 == Ny2))
	    Coefficient1 = 1.0;
	  else
	    Coefficient1 = exp(- PIOnM * Q2);
	  Coefficient1 *= FourierTransformedInteraction(Q21, Q22, Q23);
	  
	  ResultNx2 = 1.0 ; // Nx1 = 0 Nx2 = 0
	  Precision1 = ResultNx2 ;
	  Nx2 = 1.0;
	  while ((fabs(ResultNx2) + Precision1) != fabs(ResultNx2)) // Nx1 = 0 Nx2!=0
	    {
	      Q5 = this->InvRatio * Nx2 * Nx2;
	      Q21 = (this->InvRatio * Nx1 * Nx1) + (this->Ratio * Ny1 * Ny1);
	      Q22 = (this->Ratio * Ny2 * Ny2);
	      Q23 = (this->InvRatio * Nx2 * Nx2) + (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
	      Q21 *= Q2Factor;
	      Q22 *= Q2Factor;
	      Q23 *= Q2Factor;
	      Precision1 = 2.0 * exp(- PIOnM * Q5);
	      Precision1 *= FourierTransformedInteraction(Q21, Q22, Q23);
	      ResultNx2 += Precision1 * cos (Nx2 * Factor2);
	      Nx2 += 1.0;
	    }
	  ResultNx1 = ResultNx2;
	  Nx1 = 1.0;
	  Coefficient2 = 1.0;
	  Nx2 = 0;
	  while ((fabs(ResultNx1) + fabs(Coefficient2)) != fabs(ResultNx1))
	    {
	      ResultNx2 =  2.0 * cos(Nx1 * Factor1); // Nx1 != 0 Nx2=0
	      Q3 = this->InvRatio * Nx1 * Nx1;
	      Q21 = (this->InvRatio * Nx1 * Nx1) + (this->Ratio * Ny1 * Ny1);
	      Q22 = (this->Ratio * Ny2 * Ny2);
	      Q23 = (this->InvRatio * Nx1 * Nx1) + (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
	      Q21 *= Q2Factor;
	      Q22 *= Q2Factor;
	      Q23 *= Q2Factor;
	      Coefficient2 = exp(- PIOnM * Q3);
	      Coefficient2 *= FourierTransformedInteraction(Q21, Q22, Q23);
	      Precision = 2.0;
	      Precision1 = 2.0;
	      Nx2 = 1.0;
	      while ((fabs(ResultNx2) + Precision + Precision1) != fabs(ResultNx2))// Nx1 != 0 Nx2!=0
		{
		  Q4 = this->InvRatio * Nx2 * (Nx2 - Nx1);
		  Q5 = this->InvRatio * Nx2 * (Nx2 + Nx1);
		  if (Nx1 == Nx2)
		    Precision = 2.0;
		  else
		    Precision = 2.0 * exp(- PIOnM * Q4);
		  Precision1 = 2.0 * exp(- PIOnM * Q5);
		  Q21 = (this->InvRatio * Nx1 * Nx1) + (this->Ratio * Ny1 * Ny1);
		  Q22 = (this->InvRatio * Nx2 * Nx2) + (this->Ratio * Ny2 * Ny2);
		  Q23 = (this->InvRatio * (Nx1 + Nx2) * (Nx1 + Nx2)) + (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
		  Q21 *= Q2Factor;
		  Q22 *= Q2Factor;
		  Q23 *= Q2Factor;
		  Tmp = FourierTransformedInteraction(Q21, Q22, Q23);
		  Precision *= Tmp;
		  Precision1 *= Tmp;
		  ResultNx2 += (Precision * cos (Nx1 * Factor1 + Nx2 * Factor2) + Precision1 * cos(Nx1 * Factor1 - Nx2 * Factor2));
		  Nx2 += 1.0;
		}
	      ResultNx1 += Coefficient2 * ResultNx2;
	      Nx1 += 1.0;
	    }
	  ResultNy1 += ResultNx1 * Coefficient1; 
	  Factor2 -= M_PI;
	  Ny1 += DoubleNbrLzValue;
	}
      ResultNy2 += ResultNy1 * Coefficient;
      Factor1 -= M_PI;
      Ny2 += DoubleNbrLzValue;
    }
  
  Ny2 = IniNy2 - DoubleNbrLzValue;
  Factor1 = PremFactor1 + M_PI;
  Factor2 = PremFactor2;
  Coefficient = 1.0;
  
  while ((fabs(ResultNy2) + fabs(Coefficient)) != fabs(ResultNy2))
    {
      Q1 = this->Ratio * Ny2 * Ny2;
      if (Ny2 != 0.0)
	Coefficient = exp(- PIOnM * Q1);
      else
	Coefficient = 1.0;
      ResultNy1 = 0.0;
      Ny1 = IniNy1;
      Factor2 = PremFactor2;
      Coefficient1 = 1.0;
      Nx1 = 0;
      Nx2 = 0;
      while ((fabs(ResultNy1) + fabs(Coefficient1)) != fabs(ResultNy1))
	{
	  Q2 = this->Ratio * Ny1 * (Ny1 - Ny2);
	  Q21 = (this->Ratio * Ny1 * Ny1);
	  Q22 = (this->Ratio * Ny2 * Ny2);
	  Q23 = (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
	  Q21 *= Q2Factor;
	  Q22 *= Q2Factor;
	  Q23 *= Q2Factor;
	  if ((Ny1 == 0.0)||(Ny1 == Ny2))
	    Coefficient1 = 1.0;
	  else
	    Coefficient1 = exp(- PIOnM * Q2);
	  Coefficient1 *= FourierTransformedInteraction(Q21, Q22, Q23);
	  ResultNx2 = 1.0 ; // Nx1 = 0 Nx2=0
	  Precision1 = ResultNx2 ;
	  Nx2 = 1.0;
	  while ((fabs(ResultNx2) + Precision1) != fabs(ResultNx2)) // Nx1 = 0 Nx2!=0
	    {
	      Q5 = this->InvRatio * Nx2 * Nx2;
	      Q21 = (this->Ratio * Ny1 * Ny1);
	      Q22 = (this->InvRatio * Nx2 * Nx2) + (this->Ratio * Ny2 * Ny2);
	      Q23 = (this->InvRatio * Nx2 * Nx2) + (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
	      Q21 *= Q2Factor;
	      Q22 *= Q2Factor;
	      Q23 *= Q2Factor;
	      Precision1 = 2.0 * exp(- PIOnM * Q5);
	      Precision1 *= FourierTransformedInteraction(Q21, Q22, Q23);
	      ResultNx2 += Precision1 * cos (Nx2 * Factor2);	      
	      Nx2 += 1.0;
	    }
	  ResultNx1 = ResultNx2;
	  Nx1 = 1.0;
	  Coefficient2 = 1.0;
	  while ((fabs(ResultNx1) + fabs(Coefficient2)) != fabs(ResultNx1))
	    {
	      ResultNx2 = 2.0 * cos (Nx1 * Factor1); // Nx1 != 0 Nx2=0
	      Q3 = this->InvRatio * Nx1 * Nx1;
	      Q21 = (this->InvRatio * Nx1 * Nx1) + (this->Ratio * Ny1 * Ny1);
	      Q22 = (this->Ratio * Ny2 * Ny2);
	      Q23 = (this->InvRatio * Nx1 * Nx1) + (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
	      Q21 *= Q2Factor;
	      Q22 *= Q2Factor;
	      Q23 *= Q2Factor;
	      Coefficient2 = exp(- PIOnM * Q3);
	      Coefficient2 *= FourierTransformedInteraction(Q21, Q22, Q23);
	      Precision = 2.0;
	      Precision1 = 2.0;
	      Nx2 = 1.0;
	      while ((fabs(ResultNx2) + Precision + Precision1) != fabs(ResultNx2))// Nx1 != 0 Nx2!=0
		{
		  Q4 = this->InvRatio * Nx2 * (Nx2 - Nx1);
		  Q5 = this->InvRatio * Nx2 * (Nx2 + Nx1);
		  if (Nx1 == Nx2)
		    Precision = 2.0;
		  else
		    Precision = 2.0 * exp(- PIOnM * Q4);
		  Precision1 = 2.0 * exp(- PIOnM * Q5);
		  Q21 = (this->InvRatio * Nx1 * Nx1) + (this->Ratio * Ny1 * Ny1);
		  Q22 = (this->InvRatio * Nx2 * Nx2) + (this->Ratio * Ny2 * Ny2);
		  Q23 = (this->InvRatio * (Nx1 + Nx2) * (Nx1 + Nx2)) + (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
		  Q21 *= Q2Factor;
		  Q22 *= Q2Factor;
		  Q23 *= Q2Factor;
		  Tmp = FourierTransformedInteraction(Q21, Q22, Q23);
		  Precision *= Tmp;
		  Precision1 *= Tmp;
		  ResultNx2 += (Precision * cos (Nx1 * Factor1 + Nx2 * Factor2) + Precision1 * cos(Nx1 * Factor1 - Nx2 * Factor2));
		  Nx2 += 1.0;
		}
	      ResultNx1 += Coefficient2 * ResultNx2;
	      Nx1 += 1.0;
	    }
	  ResultNy1 += ResultNx1 * Coefficient1;
	  Factor2 -= M_PI;
	  Ny1 += DoubleNbrLzValue;
	}
      ResultNy2 += ResultNy1 * Coefficient;
      Factor1 += M_PI;
      Ny2 -= DoubleNbrLzValue;
    }
  
  Ny2 = IniNy2;
  Factor1 = PremFactor1;
  Factor2 = PremFactor2 + M_PI;
  
  Coefficient = 1.0;
  while ((fabs(ResultNy2) + fabs(Coefficient)) != fabs(ResultNy2))
    {
      Q1 = this->Ratio * Ny2 * Ny2;
      if (Ny2 != 0.0)
	Coefficient = exp(- PIOnM * Q1);
      else
	Coefficient = 1.0;
      
      Ny1 = IniNy1 - DoubleNbrLzValue;
      Factor2 = PremFactor2 + M_PI;
      ResultNy1 = 0.0;
      Coefficient1 = 1.0;
      while ((fabs(ResultNy1) + fabs(Coefficient1)) != fabs(ResultNy1))
	{
	  Q2 = this->Ratio * Ny1 * (Ny1 - Ny2);
	  Q21 = (this->Ratio * Ny1 * Ny1);
	  Q22 = (this->Ratio * Ny2 * Ny2);
	  Q23 = (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
	  Q21 *= Q2Factor;
	  Q22 *= Q2Factor;
	  Q23 *= Q2Factor;
	  if ((Ny1 == 0.0)||(Ny1 == Ny2))
	    Coefficient1 = 1.0;
	  else
	    Coefficient1 = exp(- PIOnM * Q2);

	  ResultNx2 = FourierTransformedInteraction(Q21, Q22, Q23); // Nx1 = 0 Nx2=0
	  Precision1 = ResultNx2 ;
	  Nx2 = 1.0;
	  while ((fabs(ResultNx2) + Precision1) != fabs(ResultNx2)) // Nx1 = 0 Nx2!=0
	    {
	      Q5 = this->InvRatio * Nx2 * Nx2;
	      Q21 = (this->Ratio * Ny1 * Ny1);
	      Q22 = (this->InvRatio * Nx2 * Nx2) + (this->Ratio * Ny2 * Ny2);
	      Q23 = (this->InvRatio * Nx2 * Nx2) + (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
	      Q21 *= Q2Factor;
	      Q22 *= Q2Factor;
	      Q23 *= Q2Factor;
	      Precision1 = 2.0 * exp(- PIOnM * Q5) * FourierTransformedInteraction(Q21, Q22, Q23);
	      ResultNx2 += Precision1 * cos (Nx2 * Factor2);
	      Nx2 += 1.0;
	    }
	  ResultNx1 = ResultNx2;
	  Nx1 = 1.0;
	  Coefficient2 = 1.0;
	  while ((fabs(ResultNx1) + fabs(Coefficient2)) != fabs(ResultNx1))
	    {
	      Q21 = (this->InvRatio * Nx1 * Nx1) + (this->Ratio * Ny1 * Ny1);
	      Q22 = (this->Ratio * Ny2 * Ny2);
	      Q23 = (this->InvRatio * Nx1 * Nx1) + (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
	      Q21 *= Q2Factor;
	      Q22 *= Q2Factor;
	      Q23 *= Q2Factor;
	      ResultNx2 = 2.0 * cos (Nx1 * Factor1) * FourierTransformedInteraction(Q21, Q22, Q23); // Nx1 != 0 Nx2=0
	      Q3 = this->InvRatio * Nx1 * Nx1;
	      Coefficient2 = exp(- PIOnM * Q3);
	      Precision = 2.0;
	      Precision1 = 2.0;
	      Nx2 = 1.0;
	      while ((fabs(ResultNx2) + Precision + Precision1) != fabs(ResultNx2))// Nx1 != 0 Nx2!=0
		{
		  Q4 = this->InvRatio * Nx2 * (Nx2 - Nx1);
		  Q5 = this->InvRatio * Nx2 * (Nx2 + Nx1);
		  if (Nx1 == Nx2)
		    Precision = 2.0;
		  else
		    Precision = 2.0 * exp(- PIOnM * Q4);
		  Precision1 = 2.0 * exp(- PIOnM * Q5);
		  Q21 = (this->InvRatio * Nx1 * Nx1) + (this->Ratio * Ny1 * Ny1);
		  Q22 = (this->InvRatio * Nx2 * Nx2) + (this->Ratio * Ny2 * Ny2);
		  Q23 = (this->InvRatio * (Nx1 + Nx2) * (Nx1 + Nx2)) + (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
		  Q21 *= Q2Factor;
		  Q22 *= Q2Factor;
		  Q23 *= Q2Factor;
		  Tmp = FourierTransformedInteraction(Q21, Q22, Q23);
		  Precision *= Tmp;
		  Precision1 *= Tmp;
		  ResultNx2 += (Precision * cos (Nx1 * Factor1 + Nx2 * Factor2) + Precision1 * cos(Nx1 * Factor1 - Nx2 * Factor2));
		  Nx2 += 1.0;
		}
	      ResultNx1 += Coefficient2 * ResultNx2;
	      Nx1 += 1.0;
	    }
	  ResultNy1 += ResultNx1 * Coefficient1; 
	  Factor2 += M_PI;
	  Ny1 -= DoubleNbrLzValue;
	}
      ResultNy2 += ResultNy1 * Coefficient;
      Factor1 -= M_PI;
      Ny2 += DoubleNbrLzValue;
    }
  
  Ny2 = IniNy2 - DoubleNbrLzValue;
  Factor1 = PremFactor1 + M_PI;
  Factor2 = PremFactor2 + M_PI;
  
  Coefficient = 1.0;	
  while ((fabs(ResultNy2) + fabs(Coefficient)) != fabs(ResultNy2))
    {
      Q1 = this->Ratio * Ny2 * Ny2;
      if (Ny2 != 0.0)
	Coefficient = exp(- PIOnM * Q1);
      else
	Coefficient = 1.0;
      
      Ny1 = IniNy1 - DoubleNbrLzValue;
      Factor2 = PremFactor2 + M_PI;
      ResultNy1 = 0.0;
      Coefficient1 = 1.0;
      while ((fabs(ResultNy1) + fabs(Coefficient1)) != fabs(ResultNy1))
	{
	  Q21 = (this->Ratio * Ny1 * Ny1);
	  Q22 = (this->Ratio * Ny2 * Ny2);
	  Q23 = (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
	  Q21 *= Q2Factor;
	  Q22 *= Q2Factor;
	  Q23 *= Q2Factor;
	  Q2 = this->Ratio * Ny1 * (Ny1 - Ny2);
	  if ((Ny1 == 0.0)||(Ny1 == Ny2))
	    Coefficient1 = 1.0;
	  else
	    Coefficient1 = exp(- PIOnM * Q2);
	  
	  ResultNx2 = FourierTransformedInteraction(Q21, Q22, Q23); // Nx1 = 0 Nx2=0
	  Precision1 = ResultNx2 ;
	  Nx2 = 1.0;
	  while ((fabs(ResultNx2) + Precision1) != fabs(ResultNx2)) // Nx1 = 0 Nx2!=0
	    {
	      Q5 = this->InvRatio * Nx2 * Nx2;
	      Q21 = (this->Ratio * Ny1 * Ny1);
	      Q22 = (this->InvRatio * Nx2 * Nx2) + (this->Ratio * Ny2 * Ny2);
	      Q23 = (this->InvRatio * Nx2 * Nx2) + (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
	      Q21 *= Q2Factor;
	      Q22 *= Q2Factor;
	      Q23 *= Q2Factor;
	      Precision1 = 2.0 * exp(- PIOnM * Q5) * FourierTransformedInteraction(Q21, Q22, Q23);
	      ResultNx2 += Precision1 * cos (Nx2 * Factor2);
	      Nx2 += 1.0;
	    }
	  ResultNx1 = ResultNx2;
	  Nx1 = 1.0;
	  Coefficient2 = 1.0;
	  while ((fabs(ResultNx1) + fabs(Coefficient2)) != fabs(ResultNx1))
	    {
	      Q21 = (this->InvRatio * Nx1 * Nx1) + (this->Ratio * Ny1 * Ny1);
	      Q22 = (this->Ratio * Ny2 * Ny2);
	      Q23 = (this->InvRatio * Nx1 * Nx1) + (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
	      Q21 *= Q2Factor;
	      Q22 *= Q2Factor;
	      Q23 *= Q2Factor;
	      ResultNx2 = 2.0 * cos (Nx1 * Factor1) * FourierTransformedInteraction(Q21, Q22, Q23); // Nx1 != 0 Nx2=0
	      Q3 = this->InvRatio * Nx1 * Nx1;
	      Coefficient2 = exp(- PIOnM * Q3);
	      Precision = 2.0;
	      Precision1 = 2.0;
	      Nx2 = 1.0;
	      while ((fabs(ResultNx2) + Precision + Precision1) != fabs(ResultNx2))// Nx1 != 0 Nx2!=0
		{
		  Q4 = this->InvRatio * Nx2 * (Nx2 - Nx1);
		  Q5 = this->InvRatio * Nx2 * (Nx2 + Nx1);
		  Q21 = (this->InvRatio * Nx1 * Nx1) + (this->Ratio * Ny1 * Ny1);
		  Q22 = (this->InvRatio * Nx2 * Nx2) + (this->Ratio * Ny2 * Ny2);
		  Q23 = (this->InvRatio * (Nx1 + Nx2) * (Nx1 + Nx2)) + (this->Ratio * (Ny1 + Ny2) * (Ny1 + Ny2));
		  Q21 *= Q2Factor;
		  Q22 *= Q2Factor;
		  Q23 *= Q2Factor;
		  if (Nx1 == Nx2)
		    Precision = 2.0;
		  else
		    Precision = 2.0 * exp(- PIOnM * Q4);
		  Precision1 = 2.0 * exp(- PIOnM * Q5);		  
		  Tmp = FourierTransformedInteraction(Q21, Q22, Q23);
		  Precision *= Tmp;
		  Precision1 *= Tmp;
		  ResultNx2 += ((Precision * cos (Nx1 * Factor1 + Nx2 * Factor2) + Precision1 * cos(Nx1 * Factor1 - Nx2 * Factor2)));
		  Nx2 += 1.0;
		}
	      ResultNx1 += Coefficient2* ResultNx2;
	      Nx1 += 1.0;
	    }
	  ResultNy1 += ResultNx1 * Coefficient1; 
	  Factor2 += M_PI;
	  Ny1 -= DoubleNbrLzValue;
	}
      ResultNy2 += ResultNy1 * Coefficient;
      Factor1 += M_PI;
      Ny2 -= DoubleNbrLzValue;
     }
  return (ResultNy2 / (24.0 * (M_PI * DoubleNbrLzValue)*(M_PI * DoubleNbrLzValue)));
}

// evaluate the Fourier transformed interaction
//
// q21 = square value of the first momentum 
// q22 = square value of the second momentum 
// q23 = square value of the third momentum 
// return value = value of the Fourier transformed interaction 

double ParticleOnTorusHaffnianHamiltonian::FourierTransformedInteraction(double q21, double q22, double q23)
{
  double Tmp = q21 * q22 * (q21 + q22);
  Tmp += q21 * q23 * (q21 + q23);
  Tmp += q22 * q23 * (q22 + q23);
  Tmp -= q21 * q22 * q23 * (q21 + q22 + q23);
  return Tmp;
}
