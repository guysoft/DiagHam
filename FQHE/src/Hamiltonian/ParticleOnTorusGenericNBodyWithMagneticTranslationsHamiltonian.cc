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
#include "GeneralTools/ArrayTools.h"

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
// maxMomentum = number of flux quanta
// xMomentum = relative angular momentum along x 
// ratio = torus aspect ratio (Lx/Ly)
// nbrNBody = type of interaction i.e. the number of density operators that are involved in the interaction
// interactionName = name of the interaction, will be use to generate the interaction matrix element output file name
// nbrMonomials = number of monomials in the Fourier transformed interaction
// monomialCoefficients = coefficients in front of each monomial in the Fourier transformed interaction
// monomialDescription = description of each monomial in the Fourier transformed interaction
// regenerateElementFlag = regenerate th interaction matrix elements instead of reading them from the harddrive
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// onDiskCacheFlag = flag to indicate if on-disk cache has to be used to store matrix elements
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian(ParticleOnTorusWithMagneticTranslations* particles, int nbrParticles, int maxMomentum, int xMomentum, double ratio,
															       int nbrNBody, char* interactionName, 
															       int nbrMonomials, double* monomialCoefficients, int** monomialDescription, 
															       bool regenerateElementFlag, AbstractArchitecture* architecture, long memory, bool onDiskCacheFlag, 
															       char* precalculationFileName)
{
  this->Particles = particles;
  this->LzMax = maxMomentum - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaxMomentum = maxMomentum;
  this->XMomentum = xMomentum;
  this->NbrParticles = nbrParticles;
  this->MomentumModulo = FindGCD(this->NbrParticles, this->MaxMomentum);
  this->NBodyValue = nbrNBody;
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
  char* MatrixElementFileName = new char [strlen(interactionName) + 256];
  sprintf (MatrixElementFileName, "%dbody_%s_2s_%d_ratio_%.14f.dat", this->NBodyValue, interactionName, this->MaxMomentum, this->Ratio);
  this->NbrMonomials = nbrMonomials;
  this->MonomialCoefficients = monomialCoefficients;
  this->MonomialDescription = monomialDescription;
  if ((regenerateElementFlag == true) || (!(IsFile(MatrixElementFileName))))
    {
      this->EvaluateInteractionFactors();
      this->WriteInteractionFactors(MatrixElementFileName);
    }
  else
    {
      this->ReadInteractionFactors(MatrixElementFileName);
    }
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
  this->QxValues = new double [this->NBodyValue];
  this->QyValues = new double [this->NBodyValue];
  this->Q2Values = new double [this->NBodyValue];
  this->CosineCoffients = new double [this->NBodyValue];

  int** LinearizedNBodySectorIndicesPerSum = new int* [this->NbrNBodySectorSums];
  for (int i = 0; i < this->NbrNBodySectorSums; ++i)
    {
      LinearizedNBodySectorIndicesPerSum[i] = new int [NbrNBodySectorIndicesPerSum[i]];
      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	{
	  int Tmp = 0;
	  for (int k = this->NBodyValue - 1; k >= 0; --k)
	    {
	      Tmp *= this->MaxMomentum;
	      Tmp += this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k];
	    }
	  LinearizedNBodySectorIndicesPerSum[i][j1] = Tmp;
	}
      SortArrayUpOrdering<int>(LinearizedNBodySectorIndicesPerSum[i], this->NbrNBodySectorIndicesPerSum[i]);
      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	{
	  int Tmp = LinearizedNBodySectorIndicesPerSum[i][j1];
	  for (int k = 0; k < this->NBodyValue; ++k)
	    { 
	      this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k] = Tmp % this->MaxMomentum;
	      Tmp /= this->MaxMomentum;
	    }	  
	}
    }
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
	      Index = j1 * this->NbrNBodySectorIndicesPerSum[i];
	      for (int j2 = 0; j2 <= j1; ++j2)
		{
		  for (int k = 0; k < this->NBodyValue; ++k)
		    {
		      TmpNIndices[k]  = this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k];
		      TmpMIndices[k] = this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + k];
		    }
// 		  if (this->FindCanonicalIndices(TmpMIndices, TmpNIndices, i, TmpCanonicalMIndices TmpCanonicalNIndices, TmpCanonicalSum) == true)
// 		    {
		  double TmpInteraction = 0.0;
		  long NbrOperations = 0l;
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
			  double Tmp = this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, NbrOperations);
			  TmpInteraction += PermutationSign[l1] * PermutationSign[l2] * Tmp;
			}
		    }
		  this->NBodyInteractionFactors[i][Index] = TmpInteraction;		  
		  for (int k = 0; k < this->NBodyValue; ++k)
		    cout << TmpMIndices[k] << " ";
		  cout << "| ";
		  for (int k = 0; k < this->NBodyValue; ++k)
		    cout << TmpNIndices[k] << " ";
		  cout << " = " << TmpInteraction << " (" << NbrOperations << ")"<< endl;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = j1 + 1; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  this->NBodyInteractionFactors[i][j1 * this->NbrNBodySectorIndicesPerSum[i] + j2] = this->NBodyInteractionFactors[i][j2 * this->NbrNBodySectorIndicesPerSum[i] + j1];
		  TotalNbrInteractionFactors++;		  
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
	      Index = j1 * this->NbrNBodySectorIndicesPerSum[i];
	      for (int j2 = 0; j2 <= j1; ++j2)// this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  for (int k = 0; k < this->NBodyValue; ++k)
		    {
		      TmpNIndices[k] = this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k];
		    }
		  double TmpInteraction = 0.0;
		  long NbrOperations = 0l;
		  for (int k = 0; k < this->NBodyValue; ++k)
		    {
		      TmpMIndices[k] = this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + k];
		    }
		  TmpInteraction += this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, NbrOperations);
		  while (std::prev_permutation(TmpMIndices, TmpMIndices  + this->NBodyValue))
		    {
		      TmpInteraction += this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, NbrOperations);		    
		    }
		  while (std::prev_permutation(TmpNIndices, TmpNIndices  + this->NBodyValue))
		    {
		      for (int k = 0; k < this->NBodyValue; ++k)
			{
			  TmpMIndices[k] = this->NBodySectorIndicesPerSum[i][(j2 * this->NBodyValue) + k];
			}
		      TmpInteraction += this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, NbrOperations);
		      while (std::prev_permutation(TmpMIndices, TmpMIndices  + this->NBodyValue))
			{
			  TmpInteraction += this->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices, NbrOperations);
			}
		    }
		  this->NBodyInteractionFactors[i][Index] = TmpInteraction;
		  for (int k = 0; k < this->NBodyValue; ++k)
		    cout << TmpMIndices[k] << " ";
		  cout << "| ";
		  for (int k = 0; k < this->NBodyValue; ++k)
		    cout << TmpNIndices[k] << " ";
		  cout << " = " << TmpInteraction << " (" << NbrOperations << ")"<< endl;
		  TotalNbrInteractionFactors++;
		  ++Index;
		}
	    }
	  for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	    {
	      for (int j2 = j1 + 1; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
		{
		  this->NBodyInteractionFactors[i][j1 * this->NbrNBodySectorIndicesPerSum[i] + j2] = this->NBodyInteractionFactors[i][j2 * this->NbrNBodySectorIndicesPerSum[i] + j1];
		  TotalNbrInteractionFactors++;		  
		}
	    }
	}
    }
  delete[] this->QxValues;
  delete[] this->QyValues;
  delete[] this->Q2Values;
  delete[] this->CosineCoffients;
  cout << "nbr interaction = " << TotalNbrInteractionFactors << endl;
  cout << "====================================" << endl;
}
	
  
// evaluate the numerical coefficient  in front of the Prod a^+_mi Prod a+_n coupling term
//
// mIndices = array containing the creation operator indices
// nIndices = array containing the annihilation operator indices
// return value = numerical coefficient  

double ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::EvaluateInteractionCoefficient(int* mIndices, int* nIndices, long& nbrOperations)
{
  int Tmp;
  double Prefactor = powl(this->MaxMomentum, -this->NBodyValue + 1.0);
  for (int i = 0; i < this->NBodyValue; ++i)
    {
      this->QxValues[i] = 0.0;
      this->QyValues[i] = (double) (nIndices[i] - mIndices[i]);
      this->CosineCoffients[i] = 2.0 * ((double) mIndices[i]);
    }  
  double CurrentPrecision;
  double Coefficient = Prefactor * this->RecursiveEvaluateInteractionCoefficient(0, 0.0, 0.0, 0.0, 0.0, CurrentPrecision, nbrOperations);
  return Coefficient;
}
  

double ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::RecursiveEvaluateInteractionCoefficient(int xPosition, double currentSumQx, double currentSumQy, double currentSumQ2, double currentSumPhase, double& currentPrecision, long& nbrOperations)
{
  if (xPosition < (this->NBodyValue - 1))
    {
      double TotalCoefficient  = 0.0;
      double Coefficient  = 1.0;
      int CurrentQy = this->QyValues[xPosition];
      currentPrecision = 0.0;
      double TmpPrecision;
      double Tmp;
      while ((fabs(Coefficient) + fabs(TotalCoefficient)) != fabs(TotalCoefficient))
	{	        
	  Tmp = (this->QyValues[xPosition] * this->QyValues[xPosition] * this->Ratio);
	  this->QxValues[xPosition] = 0.0;
	  this->Q2Values[xPosition] = Tmp;
	  Coefficient = this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + this->QxValues[xPosition], 
								      currentSumQy + this->QyValues[xPosition], 
								      currentSumQ2 + this->Q2Values[xPosition], 
								      currentSumPhase +  this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]), TmpPrecision, nbrOperations);
	  currentPrecision += TmpPrecision;
	  if (Coefficient == 0.0)
	    TmpPrecision = 1.0;
	  while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
	    {	  
	      ++this->QxValues[xPosition];
	      this->Q2Values[xPosition] = (this->QxValues[xPosition] * this->QxValues[xPosition] * this->InvRatio) + Tmp;
	      Coefficient += this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + this->QxValues[xPosition], 
									   currentSumQy + this->QyValues[xPosition], 
									   currentSumQ2 + this->Q2Values[xPosition], 
									   currentSumPhase +  this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]), TmpPrecision, nbrOperations);
	      currentPrecision += TmpPrecision;
	    }
	  this->QxValues[xPosition] = 0.0;
	  if (Coefficient == 0.0)
	    TmpPrecision = 1.0;
	  else
	    TmpPrecision = 2.0 * fabs(Coefficient);
	  while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
	    {	  
	      --this->QxValues[xPosition];
	      this->Q2Values[xPosition] = (this->QxValues[xPosition] * this->QxValues[xPosition] * this->InvRatio) + Tmp;
	      Coefficient += this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + this->QxValues[xPosition], 
									   currentSumQy + this->QyValues[xPosition], 
									   currentSumQ2 + this->Q2Values[xPosition], 
									   currentSumPhase +  this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]), TmpPrecision, nbrOperations);
	      currentPrecision += TmpPrecision;
	    }
	  this->QyValues[xPosition] +=  (double) this->MaxMomentum;
	  TotalCoefficient += Coefficient;
	}
      this->QyValues[xPosition] = CurrentQy -  (double) this->MaxMomentum;
      if (TotalCoefficient == 0.0)
	Coefficient = 1.0;
      else
	Coefficient = 2.0 * TotalCoefficient;
      while ((fabs(Coefficient) + fabs(TotalCoefficient)) != fabs(TotalCoefficient))
	{	        
	  this->QxValues[xPosition] = 0.0;
	  Tmp = (this->QyValues[xPosition] * this->QyValues[xPosition] * this->Ratio);
	  this->Q2Values[xPosition] = Tmp;
	  Coefficient = this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + this->QxValues[xPosition], 
								      currentSumQy + this->QyValues[xPosition], 
								      currentSumQ2 + this->Q2Values[xPosition], 
								      currentSumPhase +  this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]), TmpPrecision, nbrOperations);
	  currentPrecision += TmpPrecision;
	  if (Coefficient == 0.0)
	    TmpPrecision = 1.0;
	  while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
	    {	  
	      ++this->QxValues[xPosition];
	      this->Q2Values[xPosition] = (this->QxValues[xPosition] * this->QxValues[xPosition] * this->InvRatio) + Tmp;
	      Coefficient += this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + this->QxValues[xPosition], 
									   currentSumQy + this->QyValues[xPosition], 
									   currentSumQ2 + this->Q2Values[xPosition], 
									   currentSumPhase +  this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]), TmpPrecision, nbrOperations);
	      currentPrecision += TmpPrecision;
	    }
	  this->QxValues[xPosition] = 0.0;
	  if (Coefficient == 0.0)
	    TmpPrecision = 1.0;
	  else
	    TmpPrecision = 2.0 * fabs(Coefficient);
	  while ((fabs(Coefficient) + TmpPrecision) != fabs(Coefficient))
	    {	  
	      --this->QxValues[xPosition];
	      this->Q2Values[xPosition] = (this->QxValues[xPosition] * this->QxValues[xPosition] * this->InvRatio) + Tmp;
	      Coefficient += this->RecursiveEvaluateInteractionCoefficient(xPosition + 1, currentSumQx + this->QxValues[xPosition], 
									   currentSumQy + this->QyValues[xPosition], 
									   currentSumQ2 + this->Q2Values[xPosition], 
									   currentSumPhase +  this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]), TmpPrecision, nbrOperations);
	      currentPrecision += TmpPrecision;
	    }
	  this->QyValues[xPosition] -=  (double) this->MaxMomentum;
	  TotalCoefficient += Coefficient;
	}      
      this->QyValues[xPosition] = CurrentQy;
      return TotalCoefficient;
    }
  else
    {
      double TmpExponentialFactor = M_PI / ((double) this->MaxMomentum);
      this->QxValues[xPosition] = -currentSumQx;
      this->QyValues[xPosition] = -currentSumQy;  
      this->Q2Values[xPosition] = (this->QxValues[xPosition] * this->QxValues[xPosition] * this->InvRatio) + (this->QyValues[xPosition] * this->QyValues[xPosition] * this->Ratio);
      currentSumPhase += this->QxValues[xPosition] * (this->CosineCoffients[xPosition] + this->QyValues[xPosition]);
      currentPrecision = exp(- 0.5 * TmpExponentialFactor * (this->Q2Values[xPosition] + currentSumQ2)) * this->VFactor(this->Q2Values);
      ++nbrOperations;
      return (cos(TmpExponentialFactor * currentSumPhase) * currentPrecision);
// if we do not assume that it is invariant under {qx}<->{-qx}
//      return (Complex(TmpExponentialFactor * Sum2) * currentPrecision);
    }
}



double ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::VFactor(double* q2Values)
{
  double Tmp = 0.0;
  double Tmp2;
  for (int i = 0; i < this->NbrMonomials; ++i)
    {
      Tmp2 = this->MonomialCoefficients[i];
      for (int j = 0; j < (this->NBodyValue - 1); ++j)
	{
	  for (int k = 0; k < this->MonomialDescription[i][j]; ++k)
	    Tmp2 *= q2Values[j];
	}
      Tmp += Tmp2;
    }
  return Tmp;
  //  return (-q2Values[0] * q2Values[0] * q2Values[1]);
//  return 1.0;
}


// read the interaction matrix elements from disk
//
// fileName = name of the file where the interaction matrix elements are stored
// return value = true if no error occured

bool ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::ReadInteractionFactors(char* fileName)
{
  this->GetIndices();
  int** LinearizedNBodySectorIndicesPerSum = new int* [this->NbrNBodySectorSums];
  for (int i = 0; i < this->NbrNBodySectorSums; ++i)
    {
      LinearizedNBodySectorIndicesPerSum[i] = new int [NbrNBodySectorIndicesPerSum[i]];
      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	{
	  int Tmp = 0;
	  for (int k = this->NBodyValue - 1; k >= 0; --k)
	    {
	      Tmp *= this->MaxMomentum;
	      Tmp += this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k];
	    }
	  LinearizedNBodySectorIndicesPerSum[i][j1] = Tmp;
	}
      SortArrayUpOrdering<int>(LinearizedNBodySectorIndicesPerSum[i], this->NbrNBodySectorIndicesPerSum[i]);
      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	{
	  int Tmp = LinearizedNBodySectorIndicesPerSum[i][j1];
	  for (int k = 0; k < this->NBodyValue; ++k)
	    { 
	      this->NBodySectorIndicesPerSum[i][(j1 * this->NBodyValue) + k] = Tmp % this->MaxMomentum;
	      Tmp /= this->MaxMomentum;
	    }	  
	}
    }
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "cannot open " << fileName << endl;
      return false;
    }
  this->NBodyInteractionFactors = new Complex* [this->NbrNBodySectorSums];
  for (int i = 0; i < this->NbrNBodySectorSums; ++i)
    {
      this->NBodyInteractionFactors[i] = new Complex[this->NbrNBodySectorIndicesPerSum[i] * this->NbrNBodySectorIndicesPerSum[i]];
      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	{
	  int Index = j1 * this->NbrNBodySectorIndicesPerSum[i];
	  for (int j2 = 0; j2 <= j1; ++j2)
	    {
	      ReadLittleEndian(File, this->NBodyInteractionFactors[i][Index].Re);
	      ReadLittleEndian(File, this->NBodyInteractionFactors[i][Index].Im);
	      ++Index;
	    }
	}
      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	{
	  for (int j2 = j1 + 1; j2 < this->NbrNBodySectorIndicesPerSum[i]; ++j2)
	    {
	      this->NBodyInteractionFactors[i][j1 * this->NbrNBodySectorIndicesPerSum[i] + j2] = this->NBodyInteractionFactors[i][j2 * this->NbrNBodySectorIndicesPerSum[i] + j1];
	    }
	}
    }
  File.close();
  return true;
}

// write the interaction matrix elements from disk
//
// fileName = name of the file where the interaction matrix elements are stored
// return value = true if no error occured

bool ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian::WriteInteractionFactors(char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "cannot create " << fileName << endl;
      return false;
    }
  for (int i = 0; i < this->NbrNBodySectorSums; ++i)
    {
      for (int j1 = 0; j1 < this->NbrNBodySectorIndicesPerSum[i]; ++j1)
	{
	  int Index = j1 * this->NbrNBodySectorIndicesPerSum[i];
	  for (int j2 = 0; j2 <= j1; ++j2)
	    {
	      WriteLittleEndian(File, this->NBodyInteractionFactors[i][Index].Re);
	      WriteLittleEndian(File, this->NBodyInteractionFactors[i][Index].Im);
	      ++Index;
	    }
	}
    }
  File.close();
  return true;
}

