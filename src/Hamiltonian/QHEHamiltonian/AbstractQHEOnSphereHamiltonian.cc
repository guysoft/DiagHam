////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                          to particles on a sphere with                     //
//                                                                            //
//                        last modification : 24/03/2003                      //
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
#include "Hamiltonian/QHEHamiltonian/AbstractQHEOnSphereHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Complex.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>
#include <fstream>


using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;
using std::ostream;


// destructor
//

AbstractQHEOnSphereHamiltonian::~AbstractQHEOnSphereHamiltonian()
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void AbstractQHEOnSphereHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnSphere*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* AbstractQHEOnSphereHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractQHEOnSphereHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void AbstractQHEOnSphereHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractQHEOnSphereHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
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

Complex AbstractQHEOnSphereHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

RealVector& AbstractQHEOnSphereHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination) 
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereHamiltonian::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
							       int firstComponent, int nbrComponent) 
{
  int LastComponent = firstComponent + nbrComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    vDestination[i] = 0.0;
  return this->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

RealVector& AbstractQHEOnSphereHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelAddMultiply(vSource, vDestination, 0, this->Particles->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      double TmpInteraction;
      int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
	{
	  m1 = this->M1Value[j];
	  m2 = this->M2Value[j];
	  m3 = this->M3Value[j];
	  TmpInteraction = this->InteractionFactors[j];
	  m4 = m1 + m2 - m3;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	      if (Index < Dim)
		vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	    }
	}
      m1 = this->M1Value[ReducedNbrInteractionFactors];
      m2 = this->M2Value[ReducedNbrInteractionFactors];
      m3 = this->M3Value[ReducedNbrInteractionFactors];
      TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
      m4 = m1 + m2 - m3;
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	  if (Index < Dim)
	    vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	  vDestination[i] += this->HamiltonianShift * vSource[i];
	}
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
	  int m1;
	  int m2;
	  int m3;
	  int m4;
	  double TmpInteraction;
	  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;
	  firstComponent += this->PrecalculationShift;
	  LastComponent += this->PrecalculationShift;
	  for (int k = 0; k < this->FastMultiplicationStep; ++k)
	    if (PosMod != k)
	      {		
		for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
		  {
		    m1 = this->M1Value[j];
		    m2 = this->M2Value[j];
		    m3 = this->M3Value[j];
		    TmpInteraction = this->InteractionFactors[j];
		    m4 = m1 + m2 - m3;
		    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		      {
			Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
			if (Index < Dim)
			  vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		      }
		  }
		m1 = this->M1Value[ReducedNbrInteractionFactors];
		m2 = this->M2Value[ReducedNbrInteractionFactors];
		m3 = this->M3Value[ReducedNbrInteractionFactors];
		TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
		m4 = m1 + m2 - m3;
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		    if (Index < Dim)
		      vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		    vDestination[i] += this->HamiltonianShift * vSource[i];
		  }
	      }
	  delete TmpParticles;
	}
   }
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

RealVector* AbstractQHEOnSphereHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      double TmpInteraction;
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{
	  m1 = this->M1Value[j];
	  m2 = this->M2Value[j];
	  m3 = this->M3Value[j];
	  TmpInteraction = this->InteractionFactors[j];
	  m4 = m1 + m2 - m3;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	      if (Index < Dim)
		{
		  Coefficient *= TmpInteraction;
		  for (int l = 0; l < nbrVectors; ++l)
		    vDestinations[l][Index] += Coefficient * vSources[l][i];
		}
	    }
	}
      for (int l = 0; l < nbrVectors; ++l)
	{
	  RealVector& TmpSourceVector = vSources[l];
	  RealVector& TmpDestinationVector = vDestinations[l];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	}
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  double* Coefficient2 = new double [nbrVectors];
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int j;
	  int Pos;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  Coefficient2[l] = vSources[l][k];
		  vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
		}
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		    }
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	}
      else
	{
	  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  double* Coefficient2 = new double [nbrVectors];
	  int j;
	  int TmpNbrInteraction;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  int Pos2;
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
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  Coefficient2[k] = vSources[k][l];
		  vDestinations[k][l] += this->HamiltonianShift * Coefficient2[k];
		}
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos2 = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      vDestinations[k][Pos2] += Coefficient  * Coefficient2[k];
		    }
		}
	      l += this->FastMultiplicationStep;
	      ++Pos;
	    }
	  int Index;
	  int m1;
	  int m2;
	  int m3;
	  int m4;
	  double TmpInteraction;
	  firstComponent += this->PrecalculationShift;
	  LastComponent += this->PrecalculationShift;
	  for (int k = 0; k < this->FastMultiplicationStep; ++k)
	    if (PosMod != k)
	      {		
		for (int j = 0; j < this->NbrInteractionFactors; ++j) 
		  {
		    m1 = this->M1Value[j];
		    m2 = this->M2Value[j];
		    m3 = this->M3Value[j];
		    TmpInteraction = this->InteractionFactors[j];
		    m4 = m1 + m2 - m3;
		    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		      {
			Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
			if (Index < Dim)
			  {
			    Coefficient *= TmpInteraction;
			    for (int l = 0; l < nbrVectors; ++l)
			      vDestinations[l][Index] += Coefficient * vSources[l][i];
			  }
		      }
		  }
		for (int l = 0; l < nbrVectors; ++l)
		  {
		    RealVector& TmpSourceVector = vSources[l];
		    RealVector& TmpDestinationVector = vDestinations[l];
		    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
		  }
	      }
	  delete[] Coefficient2;
	  delete TmpParticles;
	}
   }
  return vDestinations;
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractQHEOnSphereHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination) 
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnSphereHamiltonian::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							 int firstComponent, int nbrComponent)
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// return value = reference on vectorwhere result has been stored

ComplexVector& AbstractQHEOnSphereHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored
ComplexVector& AbstractQHEOnSphereHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								      int firstComponent, int nbrComponent)
{
  return vDestination;
}
 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> AbstractQHEOnSphereHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> AbstractQHEOnSphereHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long AbstractQHEOnSphereHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  
  this->NbrInteractionPerComponent = new int [EffectiveHilbertSpaceDimension];
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    this->NbrInteractionPerComponent[i] = 0;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;

  QHEParticlePrecalculationOperation Operation(this);
  this->Architecture->ExecuteOperation(&Operation);

  long Memory = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    Memory += this->NbrInteractionPerComponent[i];

  cout << "nbr interaction = " << Memory << endl;
  long TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(double*)) * EffectiveHilbertSpaceDimension;
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(double)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(double)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if (this->Particles->GetHilbertSpaceDimension() != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	    Memory += this->NbrInteractionPerComponent[i];
	}
      int* TmpNbrInteractionPerComponent = new int [ReducedSpaceDimension];
      for (int i = 0; i < ReducedSpaceDimension; ++i)
	TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[i * this->FastMultiplicationStep];
      delete[] this->NbrInteractionPerComponent;
      this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension) + (Memory * (sizeof (int) + sizeof(double)));
    }
  else
    {
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(double*)) * EffectiveHilbertSpaceDimension) + (Memory * (sizeof (int) + sizeof(double)));
      this->FastMultiplicationStep = 1;
    }

  cout << "reduction factor=" << this->FastMultiplicationStep << endl;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  return Memory;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long AbstractQHEOnSphereHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient;
  long Memory = 0;
  int m1;
  int m2;
  int m3;
  int m4;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int LastComponent = lastComponent + firstComponent;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{
	  m1 = this->M1Value[j];
	  m2 = this->M2Value[j];
	  m3 = this->M3Value[j];
	  m4 = m1 + m2 - m3;
	  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      ++Memory;
	      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
	    }
	}    
    }
  delete TmpParticles;

  return Memory;
}

// enable fast multiplication algorithm
//

void AbstractQHEOnSphereHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  int Index;
  double Coefficient;
  int m1;
  int m2;
  int m3;
  int m4;
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

  int TotalPos = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
    {
      this->InteractionPerComponentIndex[TotalPos] = new int [this->NbrInteractionPerComponent[TotalPos]];
      this->InteractionPerComponentCoefficient[TotalPos] = new double [this->NbrInteractionPerComponent[TotalPos]];      
      TmpIndexArray = this->InteractionPerComponentIndex[TotalPos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[TotalPos];
      Pos = 0;
      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{
	  m1 = this->M1Value[j];
	  m2 = this->M2Value[j];
	  m3 = this->M3Value[j];
	  m4 = m1 + m2 - m3;
	  Index = this->Particles->AdAdAA(i + this->PrecalculationShift, m1, m2, m3, m4, Coefficient);
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      TmpIndexArray[Pos] = Index;
	      TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
	      ++Pos;
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

void AbstractQHEOnSphereHamiltonian::PartialEnableFastMultiplication(int firstComponent, int lastComponent)
{
  int Index;
  double Coefficient;
  int m1;
  int m2;
  int m3;
  int m4;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int Pos;
  int Min = firstComponent / this->FastMultiplicationStep;
  int Max = lastComponent / this->FastMultiplicationStep;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
 
  for (int i = Min; i < Max; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
      this->InteractionPerComponentCoefficient[i] = new double [this->NbrInteractionPerComponent[i]];      
      TmpIndexArray = this->InteractionPerComponentIndex[i];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
      Pos = 0;
      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	{
	  m1 = this->M1Value[j];
	  m2 = this->M2Value[j];
	  m3 = this->M3Value[j];
	  m4 = m1 + m2 - m3;
	  Index = TmpParticles->AdAdAA(i * this->FastMultiplicationStep, m1, m2, m3, m4, Coefficient);
	  if (Index < this->Particles->GetHilbertSpaceDimension())
	    {
	      TmpIndexArray[Pos] = Index;
	      TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
	      ++Pos;
	    }
	}
    }
  delete TmpParticles;
}

// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs

bool AbstractQHEOnSphereHamiltonian::SavePrecalculation (char* fileName)
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
      File.write((char*) this->NbrInteractionPerComponent, sizeof(int) * Tmp);
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * this->NbrInteractionPerComponent[i]);	  
	}
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentCoefficient[i]), sizeof(double) * this->NbrInteractionPerComponent[i]);	  
	}
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

bool AbstractQHEOnSphereHamiltonian::LoadPrecalculation (char* fileName)
{
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
  this->NbrInteractionPerComponent = new int [Tmp];
  File.read((char*) this->NbrInteractionPerComponent, sizeof(int) * Tmp);
  this->InteractionPerComponentIndex = new int* [Tmp];
  this->InteractionPerComponentCoefficient = new double* [Tmp];
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
      File.read((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * this->NbrInteractionPerComponent[i]);	  
    }
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentCoefficient[i] = new double [this->NbrInteractionPerComponent[i]];
      File.read((char*) (this->InteractionPerComponentCoefficient[i]), sizeof(double) * this->NbrInteractionPerComponent[i]);	  
    }
  File.close();
  this->FastMultiplicationFlag = true;
  return true;
}

