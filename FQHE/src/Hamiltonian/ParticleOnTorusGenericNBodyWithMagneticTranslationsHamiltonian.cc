////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of hamiltonian associated to particles on a torus with         //
//                          generic n-body interaction                        //
//                                                                            //
//                        last modification : 24/01/2015                      //
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
#include "Hamiltonian/ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Architecture/AbstractArchitecture.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/Endian.h"

#include <iostream>
#include <algorithm>


// default constructor
//

ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian()
{
}

// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// lzmax = maximum Lz value reached by a particle in the state
// twoBodyDeltaStrength = strength of the additional two body delta interaction
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio,
																	 AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
																	 char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->NBodyValue = 2;//3;
  this->SqrNBodyValue = this->NBodyValue * this->NBodyValue;
  this->TwoBodyFlag = false;
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = true;
  this->OneBodyInteractionFactors = 0;
  this->Ratio = ratio;
  this->InvRatio = 1.0 / ratio;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;
  this->EvaluateExponentialFactors();
  this->HamiltonianShift = 0.0;
  this->EvaluateInteractionFactors();
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  cout << "fast memory = ";
	  PrintMemorySize(cout,TmpMemory)<<endl;
	  if (memory > 0)
	    {
	      this->EnableFastMultiplication();
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
}

// destructor
//

ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::~ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian()
{
}
  
// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particles = (ParticleOnTorusWithMagneticTranslations*) hilbertSpace;
  this->EvaluateInteractionFactors();
}
  
// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}

// evaluate all interaction factors
//   

void ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::EvaluateInteractionFactors()
{
  long TotalNbrInteractionFactors = 0l;
  this->GetIndices();
  if (this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
    {
      int NbrPermutations = 1;
      for (int i = 1; i <= this->NBodyValue; ++i)
	NbrPermutations *= i;
      int** Permutations = new int*[NbrPermutations]; 
      double* PermutationSign = new double[NbrPermutations]; 
      Permutations[0] = new int [this->NBodyValue];
      for (int i = 0; i < this->NBodyValue; ++i)
	Permutations[0][i] = i;
      PermutationSign[0] = 1.0;
      double TmpSign = 1.0;
      for (int i = 1; i < NbrPermutations; ++i)
	{
	  Permutations[i] = new int [this->NBodyValue];
	  for (int j = 0; j < this->NBodyValue; ++j)
	    Permutations[i][j] = Permutations[i - 1][j];
	  int* TmpArrayPerm = Permutations[i];
	  int Pos1 = this->NBodyValue - 1;
	  while (TmpArrayPerm[Pos1 - 1] >= TmpArrayPerm[Pos1])
	    --Pos1;
	  --Pos1;
	  int Pos2 = this->NBodyValue - 1;      
	  while (TmpArrayPerm[Pos2] <= TmpArrayPerm[Pos1])
	    --Pos2;
	  int TmpIndex = TmpArrayPerm[Pos1];
	  TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	  TmpArrayPerm[Pos2] = TmpIndex;
	  TmpSign *= -1.0;
	  Pos2 = this->NBodyValue - 1;   
	  Pos1++;
	  while (Pos1 < Pos2)
	    {
	      TmpIndex = TmpArrayPerm[Pos1];
	      TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	      TmpArrayPerm[Pos2] = TmpIndex;
	      ++Pos1;
	      --Pos2;
	      TmpSign *= -1.0;
	    }
	  PermutationSign[i] = TmpSign;
	}

      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      int* TmpNIndices =  new int [this->NBodyValue];
      int* TmpMIndices =  new int [this->NBodyValue];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  for (int k = 0; k < this->NBodyValue; ++k)
		    {
		      TmpNIndices[k]  = this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k];
		      TmpMIndices[k] = this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + k];
		    }
		  cout << this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices) << " <> ";
// 		  cout << this->EvaluateInteractionCoefficient(TmpMIndices[0], TmpMIndices[1], TmpMIndices[2], 
// 							       TmpNIndices[0], TmpNIndices[1], TmpNIndices[2]) << endl;
		  Complex TmpInteraction = 0.0;
		  for (int l1 = 0; l1 < NbrPermutations; ++l1)
		    {
		      int* TmpPerm1 = Permutations[l1];
		      for (int k = 0; k < this->NBodyValue; ++k)
			{
			  TmpNIndices[k]  = this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + TmpPerm1[k]];
			}
		      for (int l2 = 0; l2 < NbrPermutations; ++l2)
			{
			  int* TmpPerm2 = Permutations[l2];
			  for (int k = 0; k < this->NBodyValue; ++k)
			    {
			      TmpMIndices[k] = this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + TmpPerm2[k]];
			    }
			  Complex Tmp = this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices);
// 			  cout << i << " " << Index << " : " << TmpMIndices[0] << TmpMIndices[1] <<  TmpMIndices[2] << " | "
// 			       << TmpNIndices[0] << TmpNIndices[1] <<  TmpNIndices[2]
// 			       << " = " << Tmp << endl;
			  TmpInteraction += PermutationSign[l1] * PermutationSign[l2] * Tmp;
			}
		    }
// 		  cout << i << " " << Index << " : " << this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + 0]
// 		       << this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + 1]
// 		       << this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + 2] << " | "
// 		       << this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + 0]
// 		       << this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + 1]
// 		       << this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + 2] << " = "
// 		       << TmpInteraction << endl;
		  this->NBodyInteractionFactors[i][Index] = TmpInteraction;		  
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  else
    {
      this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
      int* TmpNIndices =  new int [this->NBodyValue];
      int* TmpMIndices =  new int [this->NBodyValue];
      for (int i = 0; i < this->NbrNBodySectorSums; ++i)
	{
	  this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
	  int Index = 0;
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = 0; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  for (int k = 0; k < this->NBodyValue; ++k)
		    {
		      TmpNIndices[k] = this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k];
		    }
		  Complex TmpInteraction = 0.0;
		  for (int k = 0; k < this->NBodyValue; ++k)
		    {
		      TmpMIndices[k] = this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + k];
		    }
		  TmpInteraction += this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices);
		  while (std::prev_permutation(TmpMIndices, TmpMIndices  + this->NBodyValue))
		    {
		      TmpInteraction += this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices);		    
		    }
		  while (std::prev_permutation(TmpNIndices, TmpNIndices  + this->NBodyValue))
		    {
		      for (int k = 0; k < this->NBodyValue; ++k)
			{
			  TmpMIndices[k] = this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + k];
			}
		      TmpInteraction += this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices);
		      while (std::prev_permutation(TmpMIndices, TmpMIndices  + this->NBodyValue))
			{
			  TmpInteraction += this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices);
			}
		    }
		  cout << i << " " << Index << " = " << TmpInteraction << endl;
		  this->NBodyInteractionFactors[i][Index] = TmpInteraction;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	}
    }
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
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

double ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int n1, int n2, int n3)
{
  double DoubleNbrLzValue = (double) this->NbrLzValue;
  double PIOnM = M_PI / DoubleNbrLzValue ;
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
  //cout << "coef " << m1 << " "  << m2 << " "  << m3 << " "  << n1 << " "<< n2 << " "<< n3 << endl;
  
//   if (this->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
//     {
//     }
//   else
    {
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
	  while ((fabs(ResultNy1) + fabs(Coefficient1)) != fabs(ResultNy1))
	    {
	      Q2 = this->Ratio * Ny1 * (Ny1 - Ny2);
	      if ((Ny1 == 0.0)||(Ny1 == Ny2))
		Coefficient1 = 1.0;
	      else
		Coefficient1 = exp(- PIOnM * Q2);
	      
	      
	      ResultNx2 = 1.0 ; // Nx1 = 0 Nx2 = 0
	      Precision1 = ResultNx2 ;
	      Nx2 = 1.0;
	      while ((fabs(ResultNx2) + Precision1) != fabs(ResultNx2)) // Nx1 = 0 Nx2!=0
		{
		  Q5 = this->InvRatio * Nx2 * Nx2;
		  Precision1 = 2.0 * exp(- PIOnM * Q5);
		  ResultNx2 += Precision1 * cos (Nx2 * Factor2);
		  Nx2 += 1.0;
		}
	      ResultNx1 = ResultNx2;
	      Nx1 = 1.0;
	      Coefficient2 = 1.0;
	      while ((fabs(ResultNx1) + fabs(Coefficient2)) != fabs(ResultNx1))
		{
		  ResultNx2 =  2.0 * cos(Nx1 * Factor1); // Nx1 != 0 Nx2=0
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
	  while ((fabs(ResultNy1) + fabs(Coefficient1)) != fabs(ResultNy1))
	    {
	      Q2 = this->Ratio * Ny1 * (Ny1 - Ny2);
	      if ((Ny1 == 0.0)||(Ny1 == Ny2))
		Coefficient1 = 1.0;
	      else
		Coefficient1 = exp(- PIOnM * Q2);
	      ResultNx2 = 1.0 ; // Nx1 = 0 Nx2=0
	      Precision1 = ResultNx2 ;
	      Nx2 = 1.0;
	      while ((fabs(ResultNx2) + Precision1) != fabs(ResultNx2)) // Nx1 = 0 Nx2!=0
		{
		  Q5 = this->InvRatio * Nx2 * Nx2;
		  Precision1 = 2.0 * exp(- PIOnM * Q5);
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
		      ResultNx2 += (Precision * cos (Nx1 * Factor1 + Nx2 * Factor2) + Precision1 * cos(Nx1 * Factor1 - Nx2 * Factor2));
		      Nx2 += 1.0;
		    }
		  ResultNx1 += Coefficient2* ResultNx2;
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
	      if ((Ny1 == 0.0)||(Ny1 == Ny2))
		Coefficient1 = 1.0;
	      else
		Coefficient1 = exp(- PIOnM * Q2);
	      
	      ResultNx2 = 1.0 ; // Nx1 = 0 Nx2=0
	      Precision1 = ResultNx2 ;
	      Nx2 = 1.0;
	      while ((fabs(ResultNx2) + Precision1) != fabs(ResultNx2)) // Nx1 = 0 Nx2!=0
		{
		  Q5 = this->InvRatio * Nx2 * Nx2;
		  Precision1 = 2.0 * exp(- PIOnM * Q5);
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
	      Q2 = this->Ratio * Ny1 * (Ny1 - Ny2);
	      if ((Ny1 == 0.0)||(Ny1 == Ny2))
		Coefficient1 = 1.0;
	      else
		Coefficient1 = exp(- PIOnM * Q2);
	      
	      ResultNx2 = 1.0 ; // Nx1 = 0 Nx2=0
	      Precision1 = ResultNx2 ;
	      Nx2 = 1.0;
	      while ((fabs(ResultNx2) + Precision1) != fabs(ResultNx2)) // Nx1 = 0 Nx2!=0
		{
		  Q5 = this->InvRatio * Nx2 * Nx2;
		  Precision1 = 2.0 * exp(- PIOnM * Q5);
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
		      ResultNx2 += (Precision * cos (Nx1 * Factor1 + Nx2 * Factor2) + Precision1 * cos(Nx1 * Factor1 - Nx2 * Factor2));
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
}
  
// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::EvaluateTwoBodyInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  return 0.0;
}

// evaluate the numerical coefficient  in front of the Prod a^+_mi Prod a+_n coupling term
//
// mIndices = array containing the creation operator indices
// nIndices = array containing the annihilation operator indices
// return value = numerical coefficient  

Complex ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(int* mIndices, int* nIndices)
{
  double* QxValues = new double [this->NBodyValue];
  double* QyValues = new double [this->NBodyValue];
  double* Q2Values = new double [this->NBodyValue];
  double* CosineCoffients = new double [this->NBodyValue];
  int Tmp;
  double QxFactor = sqrt(2.0 * M_PI / this->Ratio / ((double) this->MaxMomentum));
  double QyFactor = QxFactor * this->Ratio;
  for (int i = 0; i < this->NBodyValue; ++i)
    {
      QxValues[i] = 0.0;
      Tmp = (nIndices[i] - mIndices[i]);
      if (Tmp < 0)
	{
	  Tmp += this->MaxMomentum;
	}
      else
	{
	  if (Tmp >=  this->MaxMomentum)
	    {
	      Tmp -= this->MaxMomentum;
	    }
	}
      QyValues[i] = (double) Tmp;
      Q2Values[i] = (QyFactor * QyFactor * QyValues[i] * QyValues[i]);
      CosineCoffients[i] = M_PI * ((double) (nIndices[i] + mIndices[i] 
					     - nIndices[this->NBodyValue - 1] - mIndices[this->NBodyValue - 1])) /  ((double) this->MaxMomentum);
    }  
  QyFactor *= ((double) this->MaxMomentum);
  Complex Coefficient = this->RecursiveEvaluateInteractionCoefficient(QxValues, QyValues, Q2Values, CosineCoffients, QxFactor, QyFactor, 0);
  delete[] QxValues;
  delete[] QyValues;
  delete[] Q2Values;
  delete[] CosineCoffients;
//  cout << Coefficient << endl;
  return Coefficient;
}
  

Complex ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::RecursiveEvaluateInteractionCoefficient(double* qxValues, double* qyValues, double* q2Values, 
														    double* cosineCoffients, const double& qxFactor, const double& qyFactor,
														    int xPosition)
{
  if (xPosition < (this->NBodyValue - 1))
    {
      qxValues[xPosition] = 0.0;
      Complex Coefficient = this->RecursiveEvaluateInteractionCoefficient(qxValues, qyValues, q2Values, cosineCoffients, qxFactor, qyFactor, xPosition + 1);
      double CurrentQxValue = qxValues[xPosition];

      qxValues[xPosition] += qxFactor;
      Complex CurrentCoefficient = exp(-0.25 * (qxValues[xPosition] * qxValues[xPosition])) * this->RecursiveEvaluateInteractionCoefficient(qxValues, qyValues, q2Values, cosineCoffients, qxFactor, qyFactor, xPosition + 1);      
      while ((Norm(Coefficient) + Norm(CurrentCoefficient)) != Norm(Coefficient))
	{	  
	  Coefficient += Phase(cosineCoffients[xPosition] * qxValues[xPosition]) * CurrentCoefficient;
	  qxValues[xPosition] += qxFactor;
	  CurrentCoefficient = exp(-0.25 * (qxValues[xPosition] * qxValues[xPosition])) * this->RecursiveEvaluateInteractionCoefficient(qxValues, qyValues, q2Values, cosineCoffients, qxFactor, qyFactor, xPosition + 1);
	}
      qxValues[xPosition] = CurrentQxValue;
      qxValues[xPosition] - qxFactor;
      CurrentCoefficient = exp(-0.25 * (qxValues[xPosition] * qxValues[xPosition])) * this->RecursiveEvaluateInteractionCoefficient(qxValues, qyValues, q2Values, cosineCoffients, qxFactor, qyFactor, xPosition + 1);      
      while ((Norm(Coefficient) + Norm(CurrentCoefficient)) != Norm(Coefficient))
	{	  
	  Coefficient += Phase(cosineCoffients[xPosition] * qxValues[xPosition]) * CurrentCoefficient;
	  qxValues[xPosition] -= qxFactor;
	  CurrentCoefficient = exp(-0.25 * (qxValues[xPosition] * qxValues[xPosition])) * this->RecursiveEvaluateInteractionCoefficient(qxValues, qyValues, q2Values, cosineCoffients, qxFactor, qyFactor, xPosition + 1);
	}
      qxValues[xPosition] = CurrentQxValue;
      return Coefficient;
    }
  else
    {
      qxValues[xPosition] = 0.0;
      for (int k = 0; k < xPosition; ++k)
 	qxValues[xPosition] -= qxValues[k];	
      Complex Coefficient = exp(-0.25 * (qxValues[xPosition] * qxValues[xPosition])) * this->RecursiveEvaluateInteractionCoefficient2(qxValues, qyValues, q2Values, qyFactor, 0);
      return Coefficient;
    }
}


double ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::RecursiveEvaluateInteractionCoefficient2(double* qxValues, double* qyValues, double* q2Values, 
														      const double& qyFactor, int yPosition)
{
  if (yPosition < (this->NBodyValue - 1))
    {
      q2Values[yPosition] = (qxValues[yPosition] * qxValues[yPosition]) + (qyValues[yPosition] * qyValues[yPosition]);
      double Coefficient = exp(-0.25 * (qyValues[yPosition] * qyValues[yPosition])) * this->RecursiveEvaluateInteractionCoefficient2(qxValues, qyValues, q2Values, qyFactor, yPosition + 1);
      double CurrentQyValue = qyValues[yPosition];
      qyValues[yPosition] += qyFactor;
      q2Values[yPosition] = (qxValues[yPosition] * qxValues[yPosition]) + (qyValues[yPosition] * qyValues[yPosition]);
      double CurrentCoefficient = exp(-0.25 * (qyValues[yPosition] * qyValues[yPosition])) * this->RecursiveEvaluateInteractionCoefficient2(qxValues, qyValues, q2Values, qyFactor, yPosition + 1);      
      while ((fabs(Coefficient) + fabs(CurrentCoefficient)) != fabs(Coefficient))
	{	  
	  
	  Coefficient += CurrentCoefficient;
	  qyValues[yPosition] += qyFactor;
	  q2Values[yPosition] = (qxValues[yPosition] * qxValues[yPosition]) + (qyValues[yPosition] * qyValues[yPosition]);
	  CurrentCoefficient =  exp(-0.25 * (qyValues[yPosition] * qyValues[yPosition])) * this->RecursiveEvaluateInteractionCoefficient2(qxValues, qyValues, q2Values, qyFactor, yPosition + 1);
	}
      qyValues[yPosition] = CurrentQyValue;
      qyValues[yPosition] -= qyFactor;
      q2Values[yPosition] = (qxValues[yPosition] * qxValues[yPosition]) + (qyValues[yPosition] * qyValues[yPosition]);
      CurrentCoefficient = exp(-0.25 * (qyValues[yPosition] * qyValues[yPosition])) * this->RecursiveEvaluateInteractionCoefficient2(qxValues, qyValues, q2Values, qyFactor, yPosition + 1);      
      while ((fabs(Coefficient) + fabs(CurrentCoefficient)) != fabs(Coefficient))
	{	  
	  
	  Coefficient += CurrentCoefficient;
	  qyValues[yPosition] -= qyFactor;
	  q2Values[yPosition] = (qxValues[yPosition] * qxValues[yPosition]) + (qyValues[yPosition] * qyValues[yPosition]);
	  CurrentCoefficient =  exp(-0.25 * (qyValues[yPosition] * qyValues[yPosition])) * this->RecursiveEvaluateInteractionCoefficient2(qxValues, qyValues, q2Values, qyFactor, yPosition + 1);
	}
      qyValues[yPosition] = CurrentQyValue;
      return Coefficient;
    }
  else
    {
      qyValues[yPosition] = 0.0;
      for (int k = 0; k < yPosition; ++k)
	qyValues[yPosition] -= qyValues[k];	
      q2Values[yPosition] = (qxValues[yPosition] * qxValues[yPosition]) + (qyValues[yPosition] * qyValues[yPosition]);
      return exp(-0.25 * (qyValues[yPosition] * qyValues[yPosition])) * this->VFactor(q2Values);
    }
}


double ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::VFactor(double* q2Values)
{
  return 1.0;
}
