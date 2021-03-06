////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of particle on sphere density-density operator           //
//                                                                            //
//                        last modification : 10/12/2002                      //
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
#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

#include <iostream>


using std::cout;
using std::endl;


// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationIndex1 = index of the leftmost creation operator
// creationIndex2 = index of the rightmost creation operator
// annihilationIndex1 = index of the leftmost annihilation operator
// annihilationIndex2 = index of the rightmost annihilation operator

ParticleOnSphereDensityDensityOperator::ParticleOnSphereDensityDensityOperator(ParticleOnSphere* particle, int creationIndex1, int creationIndex2,
									       int annihilationIndex1, int annihilationIndex2)
{
  this->Particle= (ParticleOnSphere*) (particle->Clone());
  this->CreationIndex1 = creationIndex1;
  this->CreationIndex2 = creationIndex2;
  this->AnnihilationIndex1 = annihilationIndex1;
  this->AnnihilationIndex2 = annihilationIndex2;
}

// copy constructor
//
// oper = operator to copy
  
ParticleOnSphereDensityDensityOperator::ParticleOnSphereDensityDensityOperator(ParticleOnSphereDensityDensityOperator& oper)
{
  this->Particle = (ParticleOnSphere*) (oper.Particle->Clone());
  this->CreationIndex1 = oper.CreationIndex1;
  this->CreationIndex2 = oper.CreationIndex2;
  this->AnnihilationIndex1 = oper.AnnihilationIndex1;
  this->AnnihilationIndex2 = oper.AnnihilationIndex2;
}

// destructor
//

ParticleOnSphereDensityDensityOperator::~ParticleOnSphereDensityDensityOperator()
{
  delete this->Particle;
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereDensityDensityOperator::Clone ()
{
  return new ParticleOnSphereDensityDensityOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereDensityDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphere*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereDensityDensityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereDensityDensityOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnSphereDensityDensityOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  if (this->Particle->GetHilbertSpaceDimension() > 0)
    {
      int Dim = firstComponent + nbrComponent;
      int FullDim = this->Particle->GetHilbertSpaceDimension();
      double Coefficient = 0.0;
      double Element = 0.0;
      int Index;
      for (int i = firstComponent; i < Dim; ++i)
	{
	  Index = this->Particle->AdAdAA(i, this->CreationIndex1, this->CreationIndex2, this->AnnihilationIndex1, this->AnnihilationIndex2, Coefficient);
	  if (Index != FullDim)
	    Element += V1[Index] * V2[i] * Coefficient;      
	}
      return Complex(Element);
    }
  else
    {
      long Dim = firstComponent + nbrComponent;
      long FullDim = this->Particle->GetLargeHilbertSpaceDimension();
      double Coefficient = 0.0;
      double Element = 0.0;
      long Index;
      for (long i = firstComponent; i < Dim; ++i)
	{
	  Index = this->Particle->AdAdAA(i, this->CreationIndex1, this->CreationIndex2, this->AnnihilationIndex1, this->AnnihilationIndex2, Coefficient);
	  if (Index != FullDim)
	    Element += V1[Index] * V2[i] * Coefficient;      
	}
      return Complex(Element);
    }
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereDensityDensityOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									int firstComponent, int nbrComponent)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  int Last = firstComponent + nbrComponent;;
  int Index;
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      Index = this->Particle->AdAdAA(i, this->CreationIndex1, this->CreationIndex2, this->AnnihilationIndex1, this->AnnihilationIndex2, Coefficient);
      if (Index != Dim)
	vDestination[Index] += vSource[i] * Coefficient;
    }
  return vDestination;
}
  

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnSphereDensityDensityOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  if (this->Particle->GetHilbertSpaceDimension() > 0)
    {
      int Dim = firstComponent + nbrComponent;
      int FullDim = this->Particle->GetHilbertSpaceDimension();
      double Coefficient = 0.0;
      Complex Element = 0.0;
      int Index;
      for (int i = firstComponent; i < Dim; ++i)
	{
	  Index = this->Particle->AdAdAA(i, this->CreationIndex1, this->CreationIndex2, this->AnnihilationIndex1, this->AnnihilationIndex2, Coefficient);
	  if (Index != FullDim)
	    Element += Conj(V1[Index]) * V2[i] * Coefficient;      
	}
      return Complex(Element);
    }
  else
    {
      long Dim = firstComponent + nbrComponent;
      long FullDim = this->Particle->GetLargeHilbertSpaceDimension();
      double Coefficient = 0.0;
      Complex Element = 0.0;
      long Index;
      for (long i = firstComponent; i < Dim; ++i)
	{
	  Index = this->Particle->AdAdAA(i, this->CreationIndex1, this->CreationIndex2, this->AnnihilationIndex1, this->AnnihilationIndex2, Coefficient);
	  if (Index != FullDim)
	    Element += Conj(V1[Index]) * V2[i] * Coefficient;      
	}
      return Complex(Element);
    }
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnSphereDensityDensityOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  int Last = firstComponent + nbrComponent;;
  int Index;
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      Index = this->Particle->AdAdAA(i, this->CreationIndex1, this->CreationIndex2, this->AnnihilationIndex1, this->AnnihilationIndex2, Coefficient);
      if (Index != Dim)
	vDestination[Index] += vSource[i] * Coefficient;
    }
  return vDestination;
}
  

