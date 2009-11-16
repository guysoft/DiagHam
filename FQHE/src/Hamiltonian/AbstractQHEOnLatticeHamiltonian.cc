////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//                class of quantum Hall Hamiltonian associated                //
//                          to particles on a lattice                         //
//                                                                            //
//                        last modification : 12/02/2008                      //
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
#include "Hamiltonian/AbstractQHEOnLatticeHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>
#include <fstream>
#include <climits>
#include <cstdlib>
#include <cstring>

using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;
using std::ostream;



// default constructor
//
AbstractQHEOnLatticeHamiltonian::AbstractQHEOnLatticeHamiltonian()
{
  this->NbrQ12Indices=0;
  this->NbrRealInteractionPerComponent=0;
  this->NbrComplexInteractionPerComponent=0;
  this->LoadBalancingArray=0;
  this->NbrBalancedTasks=0;
  this->FastMultiplicationStep=0;
  this->HermitianSymmetryFlag=false;

}

// destructor
//

AbstractQHEOnLatticeHamiltonian::~AbstractQHEOnLatticeHamiltonian()
{
  if (FastMultiplicationFlag)
    {
      delete [] this->NbrRealInteractionPerComponent;
      delete [] this->NbrComplexInteractionPerComponent;
      long MinIndex, MaxIndex;
      this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
      int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;

      int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
	++ReducedSpaceDimension;
      for (int i=0; i<ReducedSpaceDimension; ++i)
	{
	  delete [] this->InteractionPerComponentIndex[i];
	  delete [] this->InteractionPerComponentCoefficientIndex[i];
	}
      delete [] this->InteractionPerComponentIndex;
      delete [] this->InteractionPerComponentCoefficientIndex;
      this->RealInteractionCoefficients.Empty();
      this->ComplexInteractionCoefficients.Empty();
      this->FastMultiplicationFlag=false;
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void AbstractQHEOnLatticeHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
    if (NbrHoppingTerms>0)
    {
      delete [] this->HoppingTerms;
      delete [] this->KineticQi;
      delete [] this->KineticQf;
    }
  if (NbrInteractionFactors>0)
    {
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      delete [] this->Q3Value;
      delete [] this->Q4Value;
    }
  if (NbrQ12Indices>0)
    {
      for (int i=0; i<NbrQ12Indices; ++i)
	{
	  delete [] this->Q3PerQ12[i];
	  delete [] this->Q4PerQ12[i];
	}
      delete [] this->NbrQ34Values;
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      
    }
  if (NbrDiagonalInteractionFactors>0)
    {
      delete [] this->DiagonalInteractionFactors;
      delete [] this->DiagonalQValues;
    }
  this->Particles = (ParticleOnLattice*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// set flux density in units of flux quanta through the lattice
//
// nbrFluxQuanta = flux quantua piercing the lattice
void AbstractQHEOnLatticeHamiltonian::SetNbrFluxQuanta(int nbrFluxQuanta)
{
  this->NbrFluxQuanta=nbrFluxQuanta;
  this->FluxDensity=((double)nbrFluxQuanta)/NbrSites;
  this->Particles->SetNbrFluxQuanta(nbrFluxQuanta);
  if (NbrHoppingTerms>0)
    {
      delete [] this->HoppingTerms;
      delete [] this->KineticQi;
      delete [] this->KineticQf;
    }
  if (NbrInteractionFactors>0)
    {
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      delete [] this->Q3Value;
      delete [] this->Q4Value;
    }
  if (NbrQ12Indices>0)
    {
      for (int i=0; i<NbrQ12Indices; ++i)
	{
	  delete [] this->Q3PerQ12[i];
	  delete [] this->Q4PerQ12[i];
	}
      delete [] this->NbrQ34Values;
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      
    }
  if (NbrDiagonalInteractionFactors>0)
    {
      delete [] this->DiagonalInteractionFactors;
      delete [] this->DiagonalQValues;
    }  
  bool EnableFastCalculation=FastMultiplicationFlag;
  if (FastMultiplicationFlag)
    {
      delete [] this->NbrRealInteractionPerComponent;
      delete [] this->NbrComplexInteractionPerComponent;
      long MinIndex, MaxIndex;
      this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
      int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
      for (int i=0; i<EffectiveHilbertSpaceDimension; ++i)
	{
	  delete [] this->InteractionPerComponentIndex[i];
	  delete [] this->InteractionPerComponentCoefficientIndex[i];
	}
      delete [] this->InteractionPerComponentIndex;
      delete [] this->InteractionPerComponentCoefficientIndex;
      this->RealInteractionCoefficients.Empty();
      this->ComplexInteractionCoefficients.Empty();
      this->FastMultiplicationFlag=false;
    }    
  this->EvaluateInteractionFactors();
  if (this->LoadedPrecalculation)
    {
      cout << "Cannot re-enable fast calculation when precalculation data was loaded from file!"<<endl;
      cout << "Reverting to slow calculation"<<endl;
    }
  else if (EnableFastCalculation)
    {
      int TmpMemory = this->FastMultiplicationMemory(0);
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
      cout << endl;
      if (AllowedMemory > 0)
	{
	  this->EnableFastMultiplication();
	}
    }
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* AbstractQHEOnLatticeHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractQHEOnLatticeHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void AbstractQHEOnLatticeHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}

// ask if Hamiltonian implements hermitian symmetry operations
//
bool AbstractQHEOnLatticeHamiltonian::IsHermitian()
{
  return HermitianSymmetryFlag;
}

// ask if Hamiltonian implements conjugate methods
//
bool AbstractQHEOnLatticeHamiltonian::IsConjugate()
{
  return true;
}

// symmetrize interaction factors to enable hermitian matrix multiplication
// return = true upon success
bool AbstractQHEOnLatticeHamiltonian::HermitianSymmetrizeInteractionFactors()
{
  if (HermitianSymmetryFlag)
    return true;

  if (this->Particles->HaveOrder()==false)
    {
      cout << "Hamiltonian tried to use hermitian symmetry, but this is not implemented in HilbertSpace!"<<endl;
      HermitianSymmetryFlag=false;
      return false;
    }

  cout << "Using hermitian symmetry"<<endl;

  int *M = new int[2];
  int *N = new int[2];

  // single particle terms
  if (NbrHoppingTerms>0)
    {
      int TmpNbrHoppingTerms = 0;
      int *Flags = new int[this->NbrHoppingTerms];
      for (int j = 0; j < NbrHoppingTerms; ++j) 
	{
	  M[0] = this->KineticQi[j];
	  N[0] = this->KineticQf[j];
	  Flags[j] = this->Particles->CheckOrder(M, N, 1);
	  // cout << "M="<<M[0]<<", N="<<N[0]<<", order: "<<Flags[j]<<" element: "<<HoppingTerms[j]<<endl;
	  if (Flags[j]>0)
	    ++TmpNbrHoppingTerms;
	  else if (Flags[j]==0)
	    {
	      ++TmpNbrHoppingTerms;
	      HoppingTerms[j]*=0.5;
	    }
	}
      Complex *TmpHoppingTerms = new Complex[TmpNbrHoppingTerms];
      int *TmpQi = new int[TmpNbrHoppingTerms];
      int *TmpQf = new int[TmpNbrHoppingTerms];
      int Pos=0;
      for (int j = 0; j < this->NbrHoppingTerms; ++j)
	if (Flags[j]>=0)
	  {
	    TmpHoppingTerms[Pos]=this->HoppingTerms[j];
	    TmpQi[Pos]=this->KineticQi[j];
	    TmpQf[Pos]=this->KineticQf[j];
	    ++Pos;
	  }
      delete [] this->HoppingTerms;
      delete [] this->KineticQi;
      delete [] this->KineticQf;
      this->HoppingTerms = TmpHoppingTerms;
      this->KineticQi = TmpQi;
      this->KineticQf = TmpQf; 
      this->NbrHoppingTerms = TmpNbrHoppingTerms;
    }
  
  if (this->NbrQ12Indices == 0)
    {
      int TmpNbrInteractionFactors = 0;
      int *Flags = new int[NbrInteractionFactors];
      for (int j = 0; j < NbrInteractionFactors; ++j) 
	{
	  M[0] = this->Q1Value[j];
	  M[1] = this->Q2Value[j];
	  N[0] = this->Q3Value[j];
	  N[1] = this->Q4Value[j];
	  Flags[j] = this->Particles->CheckOrder (M, N, 2);
	  if (Flags[j]>0)
	    ++TmpNbrInteractionFactors;
	  else if (Flags[j]==0)
	    {
	      ++TmpNbrInteractionFactors;
	      this->InteractionFactors[j]*=0.5; // diagonal term: make up for double counting
	    }
	}
      Complex* TmpInteractionFactors = new Complex[TmpNbrInteractionFactors];
      int* TmpQ1Value = new int[TmpNbrInteractionFactors];
      int* TmpQ2Value = new int[TmpNbrInteractionFactors];
      int* TmpQ3Value = new int[TmpNbrInteractionFactors];
      int* TmpQ4Value = new int[TmpNbrInteractionFactors];
      int Pos=0;
      for (int j = 0; j < NbrInteractionFactors; ++j)
	{
	  if (Flags[j]>=0)
	    {
	      TmpInteractionFactors[Pos]=InteractionFactors[j];
	      TmpQ1Value[Pos]=Q1Value[j];
	      TmpQ2Value[Pos]=Q2Value[j];
	      TmpQ3Value[Pos]=Q3Value[j];
	      TmpQ4Value[Pos]=Q4Value[j];
	      ++Pos;
	    }
	}
      delete [] InteractionFactors;
      delete [] Q1Value;
      delete [] Q2Value;
      delete [] Q3Value;
      delete [] Q4Value;
      this->InteractionFactors = TmpInteractionFactors;
      this->NbrInteractionFactors = TmpNbrInteractionFactors;
      this->Q1Value = TmpQ1Value;
      this->Q2Value = TmpQ2Value;
      this->Q3Value = TmpQ3Value;
      this->Q4Value = TmpQ4Value;
      delete [] Flags;
    }
  else
    {
      int OldNbrQ34Values;
      int* OldQ3PerQ12;
      int* OldQ4PerQ12;
      int TmpNbrQ12Values = 0;
      int* Q12Flags = new int[this->NbrQ12Indices];
      int TmpNbrQ34Values;
      int* TmpQ3PerQ12;
      int* TmpQ4PerQ12;
      int* Q34Flags;
      // quick 5count of interaction factors
      int OldNbrInteractionFactors=0;
      for (int q12 = 0; q12 < this->NbrQ12Indices; ++q12)
	OldNbrInteractionFactors+=this->NbrQ34Values[q12];      
      int TmpNbrInteractionFactors=0;
      Complex *TmpInteractionFactors=new Complex[OldNbrInteractionFactors];
      int Pos=0;
      for (int q12 = 0; q12 < this->NbrQ12Indices; ++q12)
	{
	  M[0]=this->Q1Value[q12];
	  M[1]=this->Q2Value[q12];
	  OldNbrQ34Values = this->NbrQ34Values[q12];
	  OldQ3PerQ12 = this->Q3PerQ12[q12];
	  OldQ4PerQ12 = this->Q4PerQ12[q12];
	  TmpNbrQ34Values = 0;
	  Q34Flags = new int[OldNbrQ34Values];
	  for (int q34 = 0; q34 < OldNbrQ34Values; ++q34)
	    {
	      N[0]=OldQ3PerQ12[q34];
	      N[1]=OldQ4PerQ12[q34];
	      Q34Flags[q34] = this->Particles->CheckOrder(M, N, 2);
	      // cout << "Check Order ="<< Q3Flags[m3]<<endl;
	      if (Q34Flags[q34]>0)
		{
		  ++TmpNbrQ34Values;
		  TmpInteractionFactors[TmpNbrInteractionFactors++]=this->InteractionFactors[Pos];
		}
	      else if (Q34Flags[q34]==0)
		{
		  ++TmpNbrQ34Values;
		  TmpInteractionFactors[TmpNbrInteractionFactors++]=0.5*this->InteractionFactors[Pos];
		}
	      ++Pos;
	    }
	  if (TmpNbrQ34Values>0)
	    {
	      //cout << "Q1="<<M[0]<<", M2="<<M[1]<<": ";
	      ++TmpNbrQ12Values;
	      Q12Flags[q12]=1;
	      TmpQ3PerQ12 = new int[TmpNbrQ34Values];
	      TmpQ4PerQ12 = new int[TmpNbrQ34Values];
	      int Pos2=0;
	      for (int q34 = 0; q34 < OldNbrQ34Values; ++q34)
		if (Q34Flags[q34]>=0)
		  {
		    TmpQ3PerQ12[Pos2]=OldQ3PerQ12[q34];
		    TmpQ4PerQ12[Pos2]=OldQ4PerQ12[q34];
		    //    cout << " " << TmpQ3PerQ12[Pos2];
		    Pos2++;
		  }
	      //cout << endl;
	      delete [] OldQ3PerQ12;
	      delete [] OldQ4PerQ12;
	      this->Q3PerQ12[q12] = TmpQ3PerQ12;
	      this->Q4PerQ12[q12] = TmpQ4PerQ12;
	      this->NbrQ34Values[q12] = TmpNbrQ34Values;
 	    }
	  else
	    {
	      Q12Flags[q12]=-1;
	      delete [] OldQ3PerQ12;
	      delete [] OldQ4PerQ12;
	    }
	}
      if (this->NbrQ12Indices!=TmpNbrQ12Values)
	{
	  cout << "reducing unused Q1,Q2"<<endl;
	  int **NewQ3PerQ12=new int*[TmpNbrQ12Values];
	  int **NewQ4PerQ12=new int*[TmpNbrQ12Values];
	  int *NewNbrQ34Values=new int[TmpNbrQ12Values];
	  Pos = 0;
	  for (int q12 = 0; q12 < this->NbrQ12Indices; ++q12)
	    if (Q12Flags[q12]>0)
	      {
		NewQ3PerQ12[Pos]=this->Q3PerQ12[q12];
		NewQ4PerQ12[Pos]=this->Q4PerQ12[q12];
		NewNbrQ34Values[Pos]=this->NbrQ34Values[q12];
		++Pos;
	      }
	  delete [] this->Q3PerQ12;
	  delete [] this->Q4PerQ12;
	  delete [] this->NbrQ34Values;
	  this->Q3PerQ12=NewQ3PerQ12;
	  this->Q3PerQ12=NewQ3PerQ12;
	  this->NbrQ34Values=NewNbrQ34Values;
	}
      // reduce size of table InteractionFactors to match new size, and copy contents
      delete [] this->InteractionFactors;
      this->InteractionFactors = new Complex[TmpNbrInteractionFactors];
      for (int i=0; i<TmpNbrInteractionFactors; ++i)
	this->InteractionFactors[i]=TmpInteractionFactors[i];
      delete [] TmpInteractionFactors;
    }

  // diagonal terms are always the same... so we're done

  delete [] M;
  delete [] N;
  this->HermitianSymmetryFlag=true;
  return true;
}


  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractQHEOnLatticeHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Particles->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractQHEOnLatticeHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiply(ComplexVector& vSource,
								 ComplexVector& vDestination, 
								int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();  
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      double TmpInteractionRe,TmpInteractionIm;
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      if (NbrHoppingTerms>0)
	{
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteractionRe = this->HoppingTerms[j].Re;
	      TmpInteractionIm = this->HoppingTerms[j].Im;
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination.Re(Index) += Coefficient * (TmpInteractionRe*vSource[i].Re - TmpInteractionIm*vSource[i].Im);
		      vDestination.Im(Index) += Coefficient * (TmpInteractionRe*vSource[i].Im + TmpInteractionIm*vSource[i].Re);
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteractionRe = this->HoppingTerms[ReducedNbrHoppingTerms].Re;
	  TmpInteractionIm = this->HoppingTerms[ReducedNbrHoppingTerms].Im;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpParticles->AdA(i, qf, qi, Coefficient);
	      if (Index < Dim)
		{
		  vDestination.Re(Index) += Coefficient * (TmpInteractionRe*vSource[i].Re - TmpInteractionIm*vSource[i].Im);
		  vDestination.Im(Index) += Coefficient * (TmpInteractionRe*vSource[i].Im + TmpInteractionIm*vSource[i].Re);
		}
	      vDestination.Re(i) += this->HamiltonianShift * vSource[i].Re;
	      vDestination.Im(i) += this->HamiltonianShift * vSource[i].Im;
	    }         
	}
      // four-fermion interactions:
      if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
	      TmpInteractionRe = this->InteractionFactors[j].Re;
	      TmpInteractionIm = this->InteractionFactors[j].Im;
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination.Re(Index) += Coefficient * (TmpInteractionRe*vSource[i].Re - TmpInteractionIm*vSource[i].Im);
		      vDestination.Im(Index) += Coefficient * (TmpInteractionRe*vSource[i].Im + TmpInteractionIm*vSource[i].Re);
		    }
		}
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2, TmpRe, TmpIm;
	  int ProcessedNbrInteractionFactors;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      ProcessedNbrInteractionFactors = 0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = TmpParticles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      TmpRe = vSource[i].Re*Coefficient;
		      TmpIm = vSource[i].Im*Coefficient;
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = TmpParticles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpInteractionRe = this->InteractionFactors[ProcessedNbrInteractionFactors].Re;
			      TmpInteractionIm = this->InteractionFactors[ProcessedNbrInteractionFactors].Im;
			      vDestination.Re(Index) += Coefficient2 * (TmpRe*TmpInteractionRe-TmpIm*TmpInteractionIm);
			      vDestination.Im(Index) += Coefficient2 * (TmpRe*TmpInteractionIm+TmpIm*TmpInteractionRe);
			    }
			  ++ProcessedNbrInteractionFactors;
			}
		    }
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
	    }
	}	  

      // separated diagonal terms as these will be the general rule for contact interactions
      if (NbrDiagonalInteractionFactors>0)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient = TmpParticles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
							 DiagonalInteractionFactors, DiagonalQValues);
	      vDestination.Re(i) +=  Coefficient * vSource[i].Re;
	      vDestination.Im(i) +=  Coefficient * vSource[i].Im;
	    }
	}
      
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  unsigned short* TmpCoefficientIndexArray;
	  double TmpRe, TmpIm;
	  unsigned short TmpNbrRealInteraction;
	  unsigned short TmpNbrComplexInteraction;
	  Complex *TmpCPtr;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i, ++k)
	    {
	      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
	      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
	      TmpRe = vSource[k].Re;
	      TmpIm = vSource[k].Im;
	      int Pos=0;
	      for (; Pos < TmpNbrRealInteraction; ++Pos)
		{
		  vDestination.Re(TmpIndexArray[Pos]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpRe;
		  vDestination.Im(TmpIndexArray[Pos]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpIm;
		}
	      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos)
		{
		  TmpCPtr= &(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]]);
		  vDestination.Re(TmpIndexArray[Pos]) +=  TmpCPtr->Re*TmpRe-TmpCPtr->Im*TmpIm;
		  vDestination.Im(TmpIndexArray[Pos]) +=  TmpCPtr->Re*TmpIm+TmpCPtr->Im*TmpRe;		  
		}
	      vDestination.Re(k) += this->HamiltonianShift * TmpRe;
	      vDestination.Im(k) += this->HamiltonianShift * TmpIm;	      
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }

  //cout << "vDestination:" <<endl<<vDestination<<endl;
  
  return vDestination;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  int Index;
  double TmpInteractionRe,TmpInteractionIm;
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  int* TmpIndexArray;
  unsigned short* TmpCoefficientIndexArray;
  double TmpRe, TmpIm;
  Complex *TmpCPtr, TmpC;
  int TmpNbrRealInteraction;
  int TmpNbrComplexInteraction;
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
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      TmpRe = vSource[l].Re;
      TmpIm = vSource[l].Im;
      int Pos2=0;
      for (; Pos2 < TmpNbrRealInteraction; ++Pos2)
	{
	  vDestination.Re(TmpIndexArray[Pos2]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*TmpRe;
	  vDestination.Im(TmpIndexArray[Pos2]) +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*TmpIm;
	}
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos2)
	{
	  TmpCPtr= &(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]);
	  vDestination.Re(TmpIndexArray[Pos2]) +=  TmpCPtr->Re*TmpRe-TmpCPtr->Im*TmpIm;
	  vDestination.Im(TmpIndexArray[Pos2]) +=  TmpCPtr->Re*TmpIm+TmpCPtr->Im*TmpRe;		  
	}
      vDestination.Re(l) += this->HamiltonianShift * TmpRe;
      vDestination.Im(l) += this->HamiltonianShift * TmpIm;	      
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	
	if (NbrHoppingTerms>0)
	  {
	    // deal with kinetic energy terms first!      
	    int qi;
	    int qf;
	    int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	    for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	      {
		qi = this->KineticQi[j];
		qf = this->KineticQf[j];
		TmpInteractionRe = this->HoppingTerms[j].Re;
		TmpInteractionIm = this->HoppingTerms[j].Im;
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdA(i, qf, qi, Coefficient);
		    if (Index < Dim)
		      {
			vDestination.Re(Index) += Coefficient * (TmpInteractionRe*vSource[i].Re - TmpInteractionIm*vSource[i].Im);
			vDestination.Im(Index) += Coefficient * (TmpInteractionRe*vSource[i].Im + TmpInteractionIm*vSource[i].Re);
		      }
		  }
	      }
	    qi = this->KineticQi[ReducedNbrHoppingTerms];
	    qf = this->KineticQf[ReducedNbrHoppingTerms];
	    TmpInteractionRe = this->HoppingTerms[ReducedNbrHoppingTerms].Re;
	    TmpInteractionIm = this->HoppingTerms[ReducedNbrHoppingTerms].Im;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		Index = TmpParticles->AdA(i, qf, qi, Coefficient);
		if (Index < Dim)
		  {
		    vDestination.Re(Index) += Coefficient * (TmpInteractionRe*vSource[i].Re - TmpInteractionIm*vSource[i].Im);
		    vDestination.Im(Index) += Coefficient * (TmpInteractionRe*vSource[i].Im + TmpInteractionIm*vSource[i].Re);
		  }
		vDestination.Re(i) += this->HamiltonianShift * vSource[i].Re;
		vDestination.Im(i) += this->HamiltonianShift * vSource[i].Im;
	      }         
	  }
	// four-fermion interactions:
	if (this->NbrQ12Indices == 0) // full storage
	  { 	  
	    for (int j = 0; j < NbrInteractionFactors; ++j) 
	      {
		int q1 = this->Q1Value[j];
		int q2 = this->Q2Value[j];
		int q3 = this->Q3Value[j];
		int q4 = this->Q4Value[j];
		TmpInteractionRe = this->InteractionFactors[j].Re;
		TmpInteractionIm = this->InteractionFactors[j].Im;
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		    if (Index < Dim)
		      {
			vDestination.Re(Index) += Coefficient * (TmpInteractionRe*vSource[i].Re - TmpInteractionIm*vSource[i].Im);
			vDestination.Im(Index) += Coefficient * (TmpInteractionRe*vSource[i].Im + TmpInteractionIm*vSource[i].Re);
		      }
		  }
	      }
	  }
	else // intelligent storage
	  {
	    double Coefficient2, TmpRe, TmpIm;
	    int ProcessedNbrInteractionFactors;
	    int TmpNbrQ34Values;
	    int* TmpQ3Values;
	    int* TmpQ4Values;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		ProcessedNbrInteractionFactors = 0;
		for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		  {
		    Coefficient = TmpParticles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
		    if (Coefficient != 0.0)
		      {
			TmpRe = vSource[i].Re*Coefficient;
			TmpIm = vSource[i].Im*Coefficient;
			TmpNbrQ34Values = this->NbrQ34Values[i12];
			TmpQ3Values = this->Q3PerQ12[i12];
			TmpQ4Values = this->Q4PerQ12[i12];
			for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			  {
			    Index = TmpParticles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			    if (Index < Dim)
			      {
				TmpInteractionRe = this->InteractionFactors[ProcessedNbrInteractionFactors].Re;
				TmpInteractionIm = this->InteractionFactors[ProcessedNbrInteractionFactors].Im;
				vDestination.Re(Index) += Coefficient2 * (TmpRe*TmpInteractionRe-TmpIm*TmpInteractionIm);
				vDestination.Im(Index) += Coefficient2 * (TmpRe*TmpInteractionIm+TmpIm*TmpInteractionRe);
			      }
			    ++ProcessedNbrInteractionFactors;
			  }
		      }
		    else
		      ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		  }
	      }
	  }	  
  
	// separated diagonal terms as these will be the general rule for contact interactions
	if (NbrDiagonalInteractionFactors>0)
	  {
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		Coefficient = TmpParticles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
							   DiagonalInteractionFactors, DiagonalQValues);
		vDestination.Re(i) +=  Coefficient * vSource[i].Re;
		vDestination.Im(i) +=  Coefficient * vSource[i].Im;
	      }
	  }
      }
  delete TmpParticles;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using disk storage option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyDiskStorage(ComplexVector& vSource, ComplexVector& vDestination, 
									   int firstComponent, int nbrComponent)
  {
  cout << "Attention: AbstractQHEOnLatticeHamiltonian::LowLevelAddMultiplyDiskStorage must be defined" << endl;
  return vDestination;
}


// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();  
  double Coefficient;
  Complex TmpInteraction;
  Complex TmpCoefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      if (NbrHoppingTerms>0)
	{
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteraction = this->HoppingTerms[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      //cout << "element "<<qi<<"->"<<qf<<" on "<<i<<": "; 
	      Index = TmpParticles->AdA(i, qf, qi, Coefficient);
	      //cout << "target: "<<Index<<endl;
	      if (Index < Dim)
		{
		  TmpCoefficient = Coefficient * TmpInteraction;
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += this->HamiltonianShift * vSources[l][i];
	    }
	}
      // four-fermion interactions:
      if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
		    }
		}
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  int ProcessedNbrInteractionFactors;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	  Complex* TmpCoefficients = new Complex[nbrVectors];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      ProcessedNbrInteractionFactors = 0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = TmpParticles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      for (int l = 0; l < nbrVectors; ++l)
			TmpCoefficients[l] = vSources[l][i]*Coefficient;
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = TmpParticles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpCoefficient = this->InteractionFactors[ProcessedNbrInteractionFactors] * Coefficient2;
			      for (int l = 0; l < nbrVectors; ++l)
				vDestinations[l][Index] += TmpCoefficient * TmpCoefficients[l];
			    }
			  ++ProcessedNbrInteractionFactors;
			}
		    }
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
	    }
	  delete [] TmpCoefficients;
	}	  

      // separated diagonal terms as these will be the general rule for contact interactions
      if (NbrDiagonalInteractionFactors>0)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient = TmpParticles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
							 DiagonalInteractionFactors, DiagonalQValues);
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] +=  Coefficient * vSources[l][i];
	    }
	}
      
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  int Index;
	  double TmpRealCoefficient;
	  unsigned short* TmpCoefficientIndexArray;
	  unsigned short TmpNbrRealInteraction;
	  unsigned short TmpNbrComplexInteraction;
	  Complex* Coefficient2 = new Complex [nbrVectors];
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i, ++k)
	    {
	      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
	      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  Coefficient2[l] = vSources[l][k];
		  vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
		}
	      int Pos=0;
	      for (; Pos < TmpNbrRealInteraction; ++Pos)
		{
		  Index = TmpIndexArray[Pos];
		  TmpRealCoefficient = RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]];
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][Index] +=  TmpRealCoefficient * Coefficient2[l];
		}
	      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos)
		{
		  Index = TmpIndexArray[Pos];
		  TmpCoefficient = ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]];
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][Index] +=  TmpCoefficient * Coefficient2[l];
		}
	    }
	  delete [] Coefficient2;
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
											   int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply"<<endl;
  return vDestinations;
}

// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using disk storage option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyDiskStorage(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										   int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::LowLevelMultipleAddMultiplyDiskStorage!"<<endl;
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiply(ComplexVector& vSource,
									     ComplexVector& vDestination, 
									     int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();  
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      Complex TmpInteraction;
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      if (NbrHoppingTerms>0)
	{
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteraction = this->HoppingTerms[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    vDestination[i] += Coefficient * Conj(TmpInteraction) *vSource[Index];
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpParticles->AdA(i, qf, qi, Coefficient);
	      if (Index < Dim)
		vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    }
	}
      // four-fermion interactions:
      if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    vDestination[i] += Conj(Coefficient * TmpInteraction) * vSource[Index];
		}
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  Complex TmpSum;
	  int ProcessedNbrInteractionFactors;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpSum = 0.0;
	      ProcessedNbrInteractionFactors = 0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = TmpParticles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = TmpParticles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    TmpSum += (Coefficient*Coefficient2)*Conj(this->InteractionFactors[ProcessedNbrInteractionFactors])*vSource[Index];
			  ++ProcessedNbrInteractionFactors;
			}
		    }
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
	      vDestination[i]+=TmpSum;
	    }
	}	  

      // separated diagonal terms as these will be the general rule for contact interactions
      if (NbrDiagonalInteractionFactors>0)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient = TmpParticles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
							 DiagonalInteractionFactors, DiagonalQValues);
	      if (Coefficient != 0.0)
		vDestination[i] +=  Coefficient * vSource[i];
	    }
	}
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  unsigned short* TmpCoefficientIndexArray;
	  Complex TmpSum;
	  unsigned short TmpNbrRealInteraction;
	  unsigned short TmpNbrComplexInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i, ++k)
	    {
	      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
	      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
	      TmpSum=0.0;
	      int Pos=0;
	      for (; Pos < TmpNbrRealInteraction; ++Pos)
		TmpSum +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*vSource[TmpIndexArray[Pos]];
	      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos)
		TmpSum +=  Conj(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]])*vSource[TmpIndexArray[Pos]];
	      vDestination[k] += TmpSum + this->HamiltonianShift * vSource[k];
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestination;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  int Index;
  Complex TmpInteraction;
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  Complex TmpSum=0.0;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  int* TmpIndexArray;
  unsigned short* TmpCoefficientIndexArray;
  Complex TmpElement, TmpC;
  int TmpNbrRealInteraction;
  int TmpNbrComplexInteraction;
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
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      TmpSum=0.0;
      int Pos2=0;
      for (; Pos2 < TmpNbrRealInteraction; ++Pos2)
	TmpSum +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*vSource[TmpIndexArray[Pos2]];
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos2)
	TmpSum +=  Conj(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos2]])*vSource[TmpIndexArray[Pos2]];
      vDestination[l] += this->HamiltonianShift * TmpElement;
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	
	if (NbrHoppingTerms>0)
	  {
	    // deal with kinetic energy terms first!      
	    int qi;
	    int qf;
	    int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	    for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	      {
		qi = this->KineticQi[j];
		qf = this->KineticQf[j];
		TmpInteraction = this->HoppingTerms[j];
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdA(i, qf, qi, Coefficient);
		    if (Index < Dim)
		      vDestination[i] += Coefficient * Conj(TmpInteraction)*vSource[Index];
		  }
	      }
	    qi = this->KineticQi[ReducedNbrHoppingTerms];
	    qf = this->KineticQf[ReducedNbrHoppingTerms];
	    TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms];
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		Index = TmpParticles->AdA(i, qf, qi, Coefficient);
		if (Index < Dim)
		  {
		    vDestination[i] += Coefficient * Conj(TmpInteraction)*vSource[Index];
		  }
		vDestination[i] += this->HamiltonianShift * vSource[i];
	      }         
	  }
	// four-fermion interactions:
	if (this->NbrQ12Indices == 0) // full storage
	  { 	  
	    for (int j = 0; j < NbrInteractionFactors; ++j) 
	      {
		int q1 = this->Q1Value[j];
		int q2 = this->Q2Value[j];
		int q3 = this->Q3Value[j];
		int q4 = this->Q4Value[j];
		TmpInteraction = Conj(this->InteractionFactors[j]);
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		    if (Index < Dim)
		      vDestination[i] += Coefficient * TmpInteraction*vSource[Index];
		  }
	      }
	  }
	else // intelligent storage
	  {
	    double Coefficient2;
	    Complex TmpSum;
	    int ProcessedNbrInteractionFactors;
	    int TmpNbrQ34Values;
	    int* TmpQ3Values;
	    int* TmpQ4Values;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		TmpSum=0.0;
		ProcessedNbrInteractionFactors = 0;
		for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		  {
		    Coefficient = TmpParticles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
		    if (Coefficient != 0.0)
		      {
			TmpNbrQ34Values = this->NbrQ34Values[i12];
			TmpQ3Values = this->Q3PerQ12[i12];
			TmpQ4Values = this->Q4PerQ12[i12];
			for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			  {
			    Index = TmpParticles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			    if (Index < Dim)
			      {
				TmpSum += Coefficient * Coefficient2 *
				  Conj(this->InteractionFactors[ProcessedNbrInteractionFactors])*vSource[Index];
			      }
			    ++ProcessedNbrInteractionFactors;
			  }
		      }
		    else
		      ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		  }
	      }
	  }	  
  
	// separated diagonal terms as these will be the general rule for contact interactions
	if (NbrDiagonalInteractionFactors>0)
	  {
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		Coefficient = TmpParticles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
							   DiagonalInteractionFactors, DiagonalQValues);
		vDestination[i] +=  Coefficient * vSource[i];
	      }
	  }
      }
  delete TmpParticles;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using disk storage option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiplyDiskStorage(ComplexVector& vSource, ComplexVector& vDestination, 
									   int firstComponent, int nbrComponent)
  {
  cout << "Attention: AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelAddMultiplyDiskStorage must be defined" << endl;
  return vDestination;
}


// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();  
  double Coefficient;
  Complex TmpInteraction;
  Complex TmpCoefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      if (NbrHoppingTerms>0)
	{
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteraction = Conj(this->HoppingTerms[j]);
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][i] += TmpCoefficient * vSources[l][Index];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = Conj(this->HoppingTerms[ReducedNbrHoppingTerms]);
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      //cout << "element "<<qi<<"->"<<qf<<" on "<<i<<": "; 
	      Index = TmpParticles->AdA(i, qf, qi, Coefficient);
	      //cout << "target: "<<Index<<endl;
	      if (Index < Dim)
		{
		  TmpCoefficient = Coefficient * TmpInteraction;
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][i] += TmpCoefficient * vSources[l][Index];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += this->HamiltonianShift * vSources[l][i];
	    }
	}
      // four-fermion interactions:
      if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
	      TmpInteraction = Conj(this->InteractionFactors[j]);
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
		    }
		}
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  int ProcessedNbrInteractionFactors;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	  Complex* TmpSum = new Complex[nbrVectors];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      for (int l = 0; l < nbrVectors; ++l)
		TmpSum[l] = 0.0;
	      ProcessedNbrInteractionFactors = 0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = TmpParticles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = TmpParticles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpCoefficient = Conj(this->InteractionFactors[ProcessedNbrInteractionFactors])*(Coefficient * Coefficient2);
			      for (int l = 0; l < nbrVectors; ++l)
				TmpSum[l] += TmpCoefficient * vSources[l][Index];
			    }
			  ++ProcessedNbrInteractionFactors;
			}
		    }
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] +=TmpSum[l];
	    }
	  delete [] TmpSum;
	}	  

      // separated diagonal terms as these will be the general rule for contact interactions
      if (NbrDiagonalInteractionFactors>0)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient = TmpParticles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
							 DiagonalInteractionFactors, DiagonalQValues);
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] +=  Coefficient * vSources[l][i];
	    }
	}
      
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  int Index;
	  double TmpRealCoefficient;
	  unsigned short* TmpCoefficientIndexArray;
	  unsigned short TmpNbrRealInteraction;
	  unsigned short TmpNbrComplexInteraction;
	  Complex* TmpSum = new Complex [nbrVectors];
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i, ++k)
	    {
	      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
	      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
	      for (int l = 0; l < nbrVectors; ++l)
		TmpSum[l] = 0.0;
	      int Pos=0;
	      for (; Pos < TmpNbrRealInteraction; ++Pos)
		{
		  Index = TmpIndexArray[Pos];
		  TmpRealCoefficient = RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]];
		  for (int l = 0; l < nbrVectors; ++l)
		    TmpSum[l] +=  TmpRealCoefficient * vSources[l][Index];
		}
	      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos)
		{
		  Index = TmpIndexArray[Pos];
		  TmpCoefficient = Conj(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]]);
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][Index] +=  TmpCoefficient * vSources[l][Index];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		TmpSum[l] += TmpSum[l] + this->HamiltonianShift * vSources[l][k];
	    }
	  delete [] TmpSum;
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->ConjugateLowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
											   int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply"<<endl;
  return vDestinations;
}

// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using disk storage option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyDiskStorage(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										   int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::ConjugateLowLevelMultipleAddMultiplyDiskStorage!"<<endl;
  return vDestinations;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiply(ComplexVector& vSource,
								 ComplexVector& vDestination, 
								int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();  
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      Complex TmpInteraction;
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      if (NbrHoppingTerms>0)
	{
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteraction = this->HoppingTerms[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		      vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpParticles->AdA(i, qf, qi, Coefficient);
	      if (Index < Dim)
		{
		  vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		  vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
		}
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    }
	}
      // four-fermion interactions:
      if (this->NbrQ12Indices == 0) // full storage
	{
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += Coefficient * TmpInteraction*vSource[i];
		      vDestination[i] += Coefficient * Conj(TmpInteraction)*vSource[Index];
		    }
		}
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  Complex TmpSum;
	  Complex TmpCoefficient;
	  int ProcessedNbrInteractionFactors;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpSum=0.0;
	      ProcessedNbrInteractionFactors = 0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = TmpParticles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      TmpCoefficient = vSource[i];
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = TmpParticles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpInteraction = this->InteractionFactors[ProcessedNbrInteractionFactors];
			      Coefficient2*=Coefficient;
			      vDestination[Index] += Coefficient2 * TmpCoefficient * TmpInteraction;
			      TmpSum += Coefficient2*Conj(TmpInteraction)*vSource[Index];
			    }
			  ++ProcessedNbrInteractionFactors;
			}
		    }
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
	      vDestination[i] += TmpSum;
	    }
	}	  
      // separated diagonal terms as these will be the general rule for contact interactions
      if (NbrDiagonalInteractionFactors>0)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient = TmpParticles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
							 DiagonalInteractionFactors, DiagonalQValues);
	      vDestination[i] +=  Coefficient * vSource[i];
	    }
	}
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  unsigned short* TmpCoefficientIndexArray;
	  Complex TmpElement;
	  unsigned short TmpNbrRealInteraction;
	  unsigned short TmpNbrComplexInteraction;
	  Complex TmpSum;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i, ++k)
	    {
	      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
	      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
	      TmpElement = vSource[k];
	      TmpSum = 0.0;
	      int Pos=0;
	      for (; Pos < TmpNbrRealInteraction; ++Pos)
		{
		  vDestination[TmpIndexArray[Pos]] +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpElement;
		  TmpSum +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]] * vSource[TmpIndexArray[Pos]];
		}
	      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos)
		{
		  vDestination[TmpIndexArray[Pos]] +=  ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]]*TmpElement;
		  TmpSum +=  Conj(ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]]) * vSource[TmpIndexArray[Pos]];
		}
	      vDestination[k] += TmpSum + this->HamiltonianShift * TmpElement;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->HermitianLowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->HermitianLowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestination;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  int Index;
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
  int* TmpIndexArray;
  unsigned short* TmpCoefficientIndexArray;
  Complex TmpInteraction, TmpC;
  Complex TmpSum;
  int TmpNbrRealInteraction;
  int TmpNbrComplexInteraction;
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
      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[Pos];
      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[Pos];
      TmpSum=0.0;
      TmpC = vSource[l];
      int Pos2=0;
      for (; Pos2 < TmpNbrRealInteraction; ++Pos2)
	{
	  vDestination[TmpIndexArray[Pos2]] +=  RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*TmpC;
	  TmpSum += RealInteractionCoefficients[TmpCoefficientIndexArray[Pos2]]*vSource[TmpIndexArray[Pos2]];
	}
      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos2)
	{
	  TmpInteraction= ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos2]];
	  vDestination[TmpIndexArray[Pos2]] +=  TmpInteraction*TmpC;
	  TmpSum +=  Conj(TmpInteraction) * vSource[TmpIndexArray[Pos2]];
	}
      vDestination[l] += TmpSum + this->HamiltonianShift * TmpC;
      l += this->FastMultiplicationStep;
      ++Pos;
    }

  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	
	if (NbrHoppingTerms>0)
	  {
	    // deal with kinetic energy terms first!      
	    int qi;
	    int qf;
	    int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	    for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	      {
		qi = this->KineticQi[j];
		qf = this->KineticQf[j];
		TmpInteraction = this->HoppingTerms[j];
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdA(i, qf, qi, Coefficient);
		    if (Index < Dim)
		      {
			vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
			vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
		      }
		  }
	      }
	    qi = this->KineticQi[ReducedNbrHoppingTerms];
	    qf = this->KineticQf[ReducedNbrHoppingTerms];
	    TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms];
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		Index = TmpParticles->AdA(i, qf, qi, Coefficient);
		if (Index < Dim)
		  {
		    vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		    vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
		  }
		vDestination[i] += this->HamiltonianShift * vSource[i];
	      }         
	  }
	// four-fermion interactions:
	if (this->NbrQ12Indices == 0) // full storage
	  { 	  
	    for (int j = 0; j < NbrInteractionFactors; ++j) 
	      {
		int q1 = this->Q1Value[j];
		int q2 = this->Q2Value[j];
		int q3 = this->Q3Value[j];
		int q4 = this->Q4Value[j];
		TmpInteraction = this->InteractionFactors[j];
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		    if (Index < Dim)
		      {
			vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
			vDestination[i] += Coefficient * Conj(TmpInteraction) * vSource[Index];
		      }
		  }
	      }
	  }
	else // intelligent storage
	  {
	    double Coefficient2;
	    Complex TmpSum;
	    int ProcessedNbrInteractionFactors;
	    int TmpNbrQ34Values;
	    int* TmpQ3Values;
	    int* TmpQ4Values;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		TmpSum = 0.0;
		ProcessedNbrInteractionFactors = 0;
		for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		  {
		    Coefficient = TmpParticles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
		    if (Coefficient != 0.0)
		      {
			TmpC = vSource[i]*Coefficient;
			TmpNbrQ34Values = this->NbrQ34Values[i12];
			TmpQ3Values = this->Q3PerQ12[i12];
			TmpQ4Values = this->Q4PerQ12[i12];
			for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			  {
			    Index = TmpParticles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			    if (Index < Dim)
			      {
				TmpInteraction = this->InteractionFactors[ProcessedNbrInteractionFactors];
				vDestination[Index] += Coefficient2 * TmpInteraction * TmpC;
				TmpSum += (Coefficient * Coefficient2) * Conj(TmpInteraction) * vSource[Index];
			      }
			    ++ProcessedNbrInteractionFactors;
			  }
		      }
		    else
		      ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		  }
		vDestination[i] += TmpSum;
	      }
	  }	  
  
	// separated diagonal terms as these will be the general rule for contact interactions
	if (NbrDiagonalInteractionFactors>0)
	  {
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		Coefficient = TmpParticles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
							   DiagonalInteractionFactors, DiagonalQValues);
		vDestination[i] +=  Coefficient * vSource[i];
	      }
	  }
      }
  delete TmpParticles;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using disk storage option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiplyDiskStorage(ComplexVector& vSource, ComplexVector& vDestination, 
									   int firstComponent, int nbrComponent)
  {
  cout << "Attention: AbstractQHEOnLatticeHamiltonian::HermitianLowLevelAddMultiplyDiskStorage must be defined" << endl;
  return vDestination;
}


// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();  
  double Coefficient;
  Complex TmpInteraction;
  Complex TmpCoefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();
      if (NbrHoppingTerms>0)
	{
	  // deal with kinetic energy terms first!      
	  int qi;
	  int qf;
	  int ReducedNbrHoppingTerms = NbrHoppingTerms-1;
	  for (int j = 0; j < ReducedNbrHoppingTerms; ++j) 
	    {
	      qi = this->KineticQi[j];
	      qf = this->KineticQf[j];
	      TmpInteraction = this->HoppingTerms[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdA(i, qf, qi, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
		      TmpCoefficient.Conjugate();
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][i] += TmpCoefficient * vSources[l][Index];
		    }
		}
	    }
	  qi = this->KineticQi[ReducedNbrHoppingTerms];
	  qf = this->KineticQf[ReducedNbrHoppingTerms];
	  TmpInteraction = this->HoppingTerms[ReducedNbrHoppingTerms];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      //cout << "element "<<qi<<"->"<<qf<<" on "<<i<<": "; 
	      Index = TmpParticles->AdA(i, qf, qi, Coefficient);
	      //cout << "target: "<<Index<<endl;
	      if (Index < Dim)
		{
		  TmpCoefficient = Coefficient * TmpInteraction;
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
		  TmpCoefficient.Conjugate();
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][i] += TmpCoefficient * vSources[l][Index];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += this->HamiltonianShift * vSources[l][i];
	    }
	}
      // four-fermion interactions:
      if (this->NbrQ12Indices == 0) // full storage
	{ 	  
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      int q1 = this->Q1Value[j];
	      int q2 = this->Q2Value[j];
	      int q3 = this->Q3Value[j];
	      int q4 = this->Q4Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		  if (Index < Dim)
		    {
		      TmpCoefficient = Coefficient * TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += TmpCoefficient * vSources[l][i];
		      TmpCoefficient.Conjugate();
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][i] += TmpCoefficient * vSources[l][Index];
		    }
		}
	    }
	}
      else // intelligent storage
	{
	  double Coefficient2;
	  int ProcessedNbrInteractionFactors;
	  int TmpNbrQ34Values;
	  int* TmpQ3Values;
	  int* TmpQ4Values;
	  Complex* TmpCoefficients = new Complex[nbrVectors];
	  Complex* TmpSum = new Complex[nbrVectors];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      for (int l = 0; l < nbrVectors; ++l)
		TmpSum[l] = 0.0;
	      ProcessedNbrInteractionFactors = 0;
	      for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
		{
		  Coefficient = TmpParticles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
		  if (Coefficient != 0.0)
		    {
		      for (int l = 0; l < nbrVectors; ++l)
			TmpCoefficients[l] = vSources[l][i];
		      TmpNbrQ34Values = this->NbrQ34Values[i12];
		      TmpQ3Values = this->Q3PerQ12[i12];
		      TmpQ4Values = this->Q4PerQ12[i12];
		      for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
			{
			  Index = TmpParticles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpCoefficient = this->InteractionFactors[ProcessedNbrInteractionFactors] * (Coefficient * Coefficient2);
			      for (int l = 0; l < nbrVectors; ++l)
				vDestinations[l][Index] += TmpCoefficient * TmpCoefficients[l];
			      TmpCoefficient.Conjugate();
			      for (int l = 0; l < nbrVectors; ++l)
				TmpSum[l] += TmpCoefficient * vSources[l][Index];
			    }
			  ++ProcessedNbrInteractionFactors;
			}
		    }
		  else
		    ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += TmpSum[l];
	    }
	  delete [] TmpCoefficients;
	  delete [] TmpSum;
	}	  
      
      // separated diagonal terms as these will be the general rule for contact interactions
      if (NbrDiagonalInteractionFactors>0)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Coefficient = TmpParticles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
							 DiagonalInteractionFactors, DiagonalQValues);
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] +=  Coefficient * vSources[l][i];
	    }
	}
      
      delete TmpParticles;
    }
  else // fast calculation enabled
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  int Index;
	  double TmpRealCoefficient;
	  unsigned short* TmpCoefficientIndexArray;
	  unsigned short TmpNbrRealInteraction;
	  unsigned short TmpNbrComplexInteraction;
	  Complex* Coefficient2 = new Complex [nbrVectors];
	  Complex* TmpSum = new Complex[nbrVectors];
	  for (int l = 0; l < nbrVectors; ++l)
	    TmpSum[l] = 0.0;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i, ++k)
	    {
	      TmpNbrRealInteraction = this->NbrRealInteractionPerComponent[i];
	      TmpNbrComplexInteraction = this->NbrComplexInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientIndexArray = this->InteractionPerComponentCoefficientIndex[i];
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  TmpSum[l] = 0.0;
		  Coefficient2[l] = vSources[l][k];
		}
	      int Pos=0;
	      for (; Pos < TmpNbrRealInteraction; ++Pos)
		{
		  Index = TmpIndexArray[Pos];
		  TmpRealCoefficient = RealInteractionCoefficients[TmpCoefficientIndexArray[Pos]];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Index] +=  TmpRealCoefficient * Coefficient2[l];
		      TmpSum[l] += TmpRealCoefficient * vSources[l][Index];
		    }
		}
	      for (int j=0; j < TmpNbrComplexInteraction; ++j, ++Pos)
		{
		  Index = TmpIndexArray[Pos];
		  TmpCoefficient = ComplexInteractionCoefficients[TmpCoefficientIndexArray[Pos]];
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][Index] +=  TmpCoefficient * Coefficient2[l];
		  TmpCoefficient.Conjugate();
		  for (int l = 0; l < nbrVectors; ++l)
		    TmpSum[l] += TmpCoefficient * vSources[l][Index];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][k] += TmpSum[l] + this->HamiltonianShift * Coefficient2[l];
	    }
	  delete [] Coefficient2;
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->HermitianLowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	}
    }
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
											   int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply"<<endl;
  return vDestinations;
}

// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using disk storage option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyDiskStorage(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
										   int firstComponent, int nbrComponent)
{
  cout << "Calling non-defined function AbstractQHEOnLatticeHamiltonian::HermitianLowLevelMultipleAddMultiplyDiskStorage!"<<endl;
  return vDestinations;
}

 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> AbstractQHEOnLatticeHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> AbstractQHEOnLatticeHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// get the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
// nbrThreads = number of threads requested
// segmentIndices = array returning the reference to an array of the first index of each of the segments
//
bool AbstractQHEOnLatticeHamiltonian::GetLoadBalancing(int nbrTasks, long* &segmentIndices)
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;

  if (((NbrRealInteractionPerComponent!=0)||(NbrComplexInteractionPerComponent!=0))&&(this->FastMultiplicationStep!=0))
    {
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;

      if ((LoadBalancingArray==0)||(NbrBalancedTasks!=nbrTasks))
	{
	  if (LoadBalancingArray!=0)
	    delete [] LoadBalancingArray;
	  this->LoadBalancingArray = new long[nbrTasks+1];
	  this->NbrBalancedTasks=nbrTasks;
	  long TmpNbrElement=0;
	  for (int i=0; i<ReducedSpaceDimension; ++i)
	    TmpNbrElement+=NbrRealInteractionPerComponent[i]+NbrComplexInteractionPerComponent[i];
	  long TmpNbrPerSegment = TmpNbrElement/nbrTasks;
	  TmpNbrElement=0;
	  int Pos=0;
	  this->LoadBalancingArray[0]=MinIndex;
	  for (int i=0; i<ReducedSpaceDimension; ++i)
	    {
	      TmpNbrElement+=NbrRealInteractionPerComponent[i]+NbrComplexInteractionPerComponent[i];
	      if (TmpNbrElement>TmpNbrPerSegment)
		{
		  LoadBalancingArray[Pos+1]=MinIndex+i*this->FastMultiplicationStep;
		  TmpNbrElement=0;
		  ++Pos;
		}
	    }
	  LoadBalancingArray[nbrTasks]=MaxIndex+1;

	  cout << "LoadBalancingArray=["<<LoadBalancingArray[1]-LoadBalancingArray[0];
	  for (int i=2; i<=nbrTasks; ++i)
	    cout <<" "<<LoadBalancingArray[i]-LoadBalancingArray[i-1];
	  cout << "]"<< endl;
	}
    }
  else
    {
      if ((LoadBalancingArray==0)||(NbrBalancedTasks!=nbrTasks))
	{
	  if (LoadBalancingArray!=0)
	    delete [] LoadBalancingArray;
	  this->LoadBalancingArray = new long[nbrTasks+1];
	  
	  int Step = EffectiveHilbertSpaceDimension / nbrTasks;
	  this->LoadBalancingArray[0]=MinIndex;
	  for (int i=0; i<nbrTasks; ++i)
	    LoadBalancingArray[i]=MinIndex+i*Step;
	  LoadBalancingArray[nbrTasks]=MaxIndex+1;
	}
    }
  segmentIndices=LoadBalancingArray;
  return true;
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// if allowedMemory == 0, the value from preceding calls is used
// return value = amount of memory needed

long AbstractQHEOnLatticeHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  this->LoadedPrecalculation=false;
  if (allowedMemory>0)
    this->AllowedMemory = allowedMemory;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  this->NbrRealInteractionPerComponent = new unsigned short [EffectiveHilbertSpaceDimension];
  this->NbrComplexInteractionPerComponent = new unsigned short [EffectiveHilbertSpaceDimension];   
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    {
      this->NbrRealInteractionPerComponent[i] = 0x0;
      this->NbrComplexInteractionPerComponent[i] = 0x0;
    }
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start memory" << endl;

  QHEParticlePrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);
  long Memory = 0;
   
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    {
      Memory += this->NbrRealInteractionPerComponent[i];
      Memory += this->NbrComplexInteractionPerComponent[i];
    }
  
  cout << "nbr interaction = " << Memory << endl;
  
  // memory requirement, ignoring the actual storage size of the values of matrix
  // elements, which is assumed small (maybe need to add an estimate, at least)
  long TmpMemory = AllowedMemory - (2*sizeof (unsigned short) + sizeof (int*) + sizeof(unsigned short*)) * EffectiveHilbertSpaceDimension;
  cout << "of which can be stored: "<<(TmpMemory / ((int) (sizeof (int) + sizeof(unsigned short))))<<endl;
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(unsigned short)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(unsigned short)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if (this->Particles->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  // memory requirement, ignoring the actual storage size of the values of matrix
	  // elements, which is assumed small (maybe need to add an estimate, at least, again!)
	  TmpMemory = AllowedMemory - (2*sizeof (unsigned short) + sizeof (int*) + sizeof(unsigned short*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	    {
	      Memory += this->NbrRealInteractionPerComponent[i];
	      Memory += this->NbrComplexInteractionPerComponent[i];
	    }	  
	}
      
      Memory = ((2*sizeof (unsigned short) + sizeof (int*) + sizeof(unsigned short*)) * ReducedSpaceDimension) + (Memory * (sizeof (int) + sizeof(unsigned short)));
      
      if (this->DiskStorageFlag == false)
	{
	  int TotalReducedSpaceDimension = ReducedSpaceDimension;
	  unsigned short* TmpNbrRealInteractionPerComponent = new unsigned short [TotalReducedSpaceDimension];
	  unsigned short* TmpNbrComplexInteractionPerComponent = new unsigned short [TotalReducedSpaceDimension];	  
	  int Pos = 0;
	  for (int i = 0; i < ReducedSpaceDimension; ++i)
	    {
	      TmpNbrRealInteractionPerComponent[i] = this->NbrRealInteractionPerComponent[Pos];
	      TmpNbrComplexInteractionPerComponent[i] = this->NbrComplexInteractionPerComponent[Pos];
	      Pos += this->FastMultiplicationStep;
	    }
	  delete[] this->NbrRealInteractionPerComponent;
	  delete[] this->NbrComplexInteractionPerComponent;
	  this->NbrRealInteractionPerComponent = TmpNbrRealInteractionPerComponent;
	  this->NbrComplexInteractionPerComponent = TmpNbrComplexInteractionPerComponent;
	}
    }
  else
    {
      Memory = ((2*sizeof (unsigned short) + sizeof (int*) + sizeof(unsigned short*)) * EffectiveHilbertSpaceDimension) + (Memory * (sizeof (int) + sizeof(unsigned short)));
      this->FastMultiplicationStep = 1;
    }

  cout << "reduction factor=" << this->FastMultiplicationStep << endl;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  cout << "final Memory in bytes = " <<Memory<<endl;
  return Memory;    
}



// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that have to be precalcualted
// return value = number of non-zero matrix elements
//
long AbstractQHEOnLatticeHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Index;
  double Coefficient;
  long Memory = 0;
  // deal with kinetic energy terms first!      
  int qi;
  int qf;
  Complex TmpInteraction;
  ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();  
  for (int j = 0; j < NbrHoppingTerms; ++j)
    {
      qi = this->KineticQi[j];
      qf = this->KineticQf[j];
      TmpInteraction = this->HoppingTerms[j];
      if (fabs(TmpInteraction.Im)<LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpParticles->AdA(i, qf, qi, Coefficient);
	      if (Index < this->Particles->GetHilbertSpaceDimension())
		{
		  ++Memory;		
		  ++this->NbrRealInteractionPerComponent[i - this->PrecalculationShift];		  
		}
	    }
	}
      else
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpParticles->AdA(i, qf, qi, Coefficient);
	      if (Index < this->Particles->GetHilbertSpaceDimension())
		{
		  ++Memory;
		  ++this->NbrComplexInteractionPerComponent[i - this->PrecalculationShift];
		}
	    }
	}
    }
  
  // four-fermion interactions:
  if (this->NbrQ12Indices == 0) // full storage
    {
      for (int j = 0; j < NbrInteractionFactors; ++j) 
	{
	  int q1 = this->Q1Value[j];
	  int q2 = this->Q2Value[j];
	  int q3 = this->Q3Value[j];
	  int q4 = this->Q4Value[j];
	  TmpInteraction = this->InteractionFactors[j];
	  if (fabs(TmpInteraction.Im)<LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD)
	    {
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		  if (Index < this->Particles->GetHilbertSpaceDimension())
		    {
		      ++Memory;
		      ++this->NbrRealInteractionPerComponent[i - this->PrecalculationShift];
		    }
		}
	    }
	  else
	    {
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, q1, q2, q3, q4, Coefficient);
		  if (Index < this->Particles->GetHilbertSpaceDimension())
		    {
		      ++Memory;
		      ++this->NbrComplexInteractionPerComponent[i - this->PrecalculationShift];
		    }
		}
	    }	 	  
	}
    }
  else // intelligent storage
    {
      double Coefficient2;
      int ProcessedNbrInteractionFactors;
      int TmpNbrQ34Values;
      int* TmpQ3Values;
      int* TmpQ4Values;
      for (int i = firstComponent; i < LastComponent; ++i)
	{	  
	  ProcessedNbrInteractionFactors = 0;
	  for (int i12 = 0; i12 < this->NbrQ12Indices; ++i12)
	    {
	      Coefficient = TmpParticles->AA(i, this->Q1Value[i12], this->Q2Value[i12]);
	      if (Coefficient != 0.0)
		{
		  TmpNbrQ34Values = this->NbrQ34Values[i12];
		  TmpQ3Values = this->Q3PerQ12[i12];
		  TmpQ4Values = this->Q4PerQ12[i12];
		  for (int i34 = 0; i34 < TmpNbrQ34Values; ++i34)
		    {
		      Index = TmpParticles->AdAd(TmpQ3Values[i34], TmpQ4Values[i34], Coefficient2);
		      if (Index < this->Particles->GetHilbertSpaceDimension())
			{
			  ++Memory;
			  if (fabs(this->InteractionFactors[ProcessedNbrInteractionFactors].Im)<LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD)
			    ++this->NbrRealInteractionPerComponent[i - this->PrecalculationShift];
			  else
			    ++this->NbrComplexInteractionPerComponent[i - this->PrecalculationShift];
			  // cout << "4b - connecting :"<<Index<<", "<<i<<": "<<Coefficient<<"*"<<Coefficient2<<"*"<<this->InteractionFactors[ProcessedNbrInteractionFactors]<< " (q's=["<<this->Q1Value[i12]<<", "<<this->Q2Value[i12]<<", "<<TmpQ3Values[i34]<<", "<<TmpQ4Values[i34]<<"])"<<endl;
			}
		      ++ProcessedNbrInteractionFactors;
		    }
		}
	      else
		ProcessedNbrInteractionFactors += this->NbrQ34Values[i12];
	    }
	}
    }	  
  
  // separated diagonal terms as these will be the general rule for contact interactions
  if (NbrDiagonalInteractionFactors>0)
    {	  
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  Coefficient = TmpParticles->AdAdAADiagonal(i, NbrDiagonalInteractionFactors,
						     DiagonalInteractionFactors, DiagonalQValues);
	  if (fabs(Coefficient)>LATTICEHAMILTONIAN_IDENTICAL_ELEMENT_THRESHOLD)
	    {
	      ++Memory;
	      ++this->NbrRealInteractionPerComponent[i - this->PrecalculationShift];
	    }
	}
    }
  
  delete TmpParticles;

  return Memory;
}

// enable fast multiplication algorithm
//

void AbstractQHEOnLatticeHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficientIndex = new unsigned short* [ReducedSpaceDimension];

    // allocate all memory at the outset:
  for (int i = 0; i < ReducedSpaceDimension; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [this->NbrRealInteractionPerComponent[i]
						       +this->NbrComplexInteractionPerComponent[i]];
      this->InteractionPerComponentCoefficientIndex[i] = new unsigned short [this->NbrRealInteractionPerComponent[i]
									     +this->NbrComplexInteractionPerComponent[i]];
    }

  QHEParticlePrecalculationOperation Operation(this, false);
  Operation.ApplyOperation(this->Architecture);
      
  cout << "Nbr distinct matrix elements: "<<RealInteractionCoefficients.GetNbrElements()<<" real, "
       << ComplexInteractionCoefficients.GetNbrElements()<<" complex"<<endl;
   
  this->FastMultiplicationFlag = true;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}

// enable fast multiplication algorithm using on disk cache 
//
// fileName = prefix of the name of the file where temporary matrix elements will be stored

void AbstractQHEOnLatticeHamiltonian::EnableFastMultiplicationWithDiskStorage(char* fileName)
{
  cout << "Using non-defined function EnableFastMultiplicationWithDiskStorage!"<<endl;
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted

void AbstractQHEOnLatticeHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{
  int LastComponent = nbrComponent + firstComponent;
 ParticleOnLattice* TmpParticles = (ParticleOnLattice*) this->Particles->Clone();

  long TotalPos = ((firstComponent - this->PrecalculationShift - 1) / this->FastMultiplicationStep) + 1;
  int InitalPos = ((firstComponent - 1) / this->FastMultiplicationStep) + 1;
  InitalPos *= this->FastMultiplicationStep;
  for (int i = InitalPos; i < LastComponent; i += this->FastMultiplicationStep)
    {
      this->EvaluateFastMultiplicationComponent(TmpParticles, i, this->InteractionPerComponentIndex[TotalPos], 
						this->InteractionPerComponentCoefficientIndex[TotalPos], TotalPos);
    }
  
  delete TmpParticles;
}

// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs

bool AbstractQHEOnLatticeHamiltonian::SavePrecalculation (char* fileName)
{
  if (this->FastMultiplicationFlag)
    {
      ofstream File;
      File.open(fileName, ios::binary | ios::out);
      int Tmp = this->Particles->GetHilbertSpaceDimension();
      File.write((char*) &(Tmp), sizeof(int));
      File.write((char*) &(this->FastMultiplicationStep), sizeof(int));
      Tmp /= this->FastMultiplicationStep;
      if ((Tmp * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
	++Tmp;
      File.write((char*) this->NbrRealInteractionPerComponent, sizeof(unsigned short) * Tmp);
      File.write((char*) this->NbrComplexInteractionPerComponent, sizeof(unsigned short) * Tmp);
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * (this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i]));
	}
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentCoefficientIndex[i]), sizeof(unsigned short) * (this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i]));	  
	}
      RealInteractionCoefficients.WriteArray(File);
      ComplexInteractionCoefficients.WriteArray(File);
      File.close();
      return true;
    }
  else
    {
      return false;
    }
}

// load precalculations from a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be read
// return value = true if no error occurs

bool AbstractQHEOnLatticeHamiltonian::LoadPrecalculation (char* fileName)
{
  this->LoadedPrecalculation=true;
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  int Tmp;
  File.read((char*) &(Tmp), sizeof(int));
  if (Tmp != this->Particles->GetHilbertSpaceDimension())
    {
      File.close();
      return false;
    }
  File.read((char*) &(this->FastMultiplicationStep), sizeof(int));
  Tmp /= this->FastMultiplicationStep;
  if ((Tmp * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
    ++Tmp;
  this->NbrRealInteractionPerComponent = new unsigned short [Tmp];
  this->NbrComplexInteractionPerComponent = new unsigned short [Tmp];
  File.read((char*) this->NbrRealInteractionPerComponent, sizeof(unsigned short) * Tmp);
  File.read((char*) this->NbrComplexInteractionPerComponent, sizeof(unsigned short) * Tmp);

  this->InteractionPerComponentIndex = new int* [Tmp];
  this->InteractionPerComponentCoefficientIndex = new unsigned short* [Tmp];
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [(this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i])];
      File.read((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * (this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i]));	  
    }
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentCoefficientIndex[i]=new unsigned short[(this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i])];
      File.read((char*) (this->InteractionPerComponentCoefficientIndex[i]), sizeof(unsigned short) * (this->NbrRealInteractionPerComponent[i]+this->NbrComplexInteractionPerComponent[i]));	  
    }
  RealInteractionCoefficients.ReadArray(File);
  ComplexInteractionCoefficients.ReadArray(File);
   
  File.close();
  this->FastMultiplicationFlag = true;
  return true;
}

