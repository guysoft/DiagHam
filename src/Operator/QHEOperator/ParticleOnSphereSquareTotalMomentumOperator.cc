////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of particle on sphere sqare of total momentum operator         //
//                                                                            //
//                        last modification : 06/03/2005                      //
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
#include "Operator/QHEOperator/ParticleOnSphereSquareTotalMomentumOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Complex.h"

  
// constructor from default datas
//
// particle = hilbert space associated to the particles
// totalLz = momentum total value
// lzMax = maximum Lz value reached by a fermion

ParticleOnSphereSquareTotalMomentumOperator::ParticleOnSphereSquareTotalMomentumOperator(ParticleOnSphere* particle, int totalLz, int lzMax)
{
  this->Particle = (ParticleOnSphere*) (particle->Clone());
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->Coefficients = RealMatrix(this->TotalLz + 1, this->TotalLz + 1);
  for (int i = 0; i <= this->TotalLz; ++i)
    {
      double TmpCoefficient = sqrt(0.5 * ((double) ((((this->TotalLz + 2) * this->TotalLz) - (((2 * i) - this->TotalLz) * (2 * i - this->TotalLz + 2))))));
      for (int j = 0; j <= this->TotalLz; ++j)
	{
     	  this->Coefficients(i, j) = TmpCoefficient;
	}
    }
  for (int i = 0; i <= this->TotalLz; ++i)
    {
      double TmpCoefficient = sqrt(0.5 * ((double) ((((this->TotalLz + 2) * this->TotalLz) - (((2 * i) - this->TotalLz) * (2 * i - this->TotalLz - 2))))));
      for (int j = 0; j <= this->TotalLz; ++j)
	{
     	  this->Coefficients(j, i) *= TmpCoefficient;
	}
    }
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnSphereSquareTotalMomentumOperator::ParticleOnSphereSquareTotalMomentumOperator(const ParticleOnSphereSquareTotalMomentumOperator& oper)
{
  this->Particle = (ParticleOnSphere*) (oper.Particle->Clone());
  this->TotalLz = oper.TotalLz;
  this->LzMax = oper.LzMax;
  this->Coefficients = oper.Coefficients;
}

// destructor
//

ParticleOnSphereSquareTotalMomentumOperator::~ParticleOnSphereSquareTotalMomentumOperator()
{
  delete this->Particle;
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereSquareTotalMomentumOperator::Clone ()
{
  return new ParticleOnSphereSquareTotalMomentumOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereSquareTotalMomentumOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphere*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereSquareTotalMomentumOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereSquareTotalMomentumOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereSquareTotalMomentumOperator::MatrixElement (RealVector& V1, RealVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Element = 0.0;
  int Index = 0;
  for (int i = 0; i < Dim; ++i)
    {
//      Index = this->Particle->ProdAdProdA(i, this->CreationIndices, this->AnnihilationIndices, this->NbrNBody, Coefficient);
      Element += V1[Index] * V2[i] * Coefficient;
    }
  return Complex(Element);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereSquareTotalMomentumOperator::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}
   
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereSquareTotalMomentumOperator::Multiply(RealVector& vSource, RealVector& vDestination, 
								  int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Index = 0;
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
//      Index = this->Particle->ProdAdProdA(i, this->CreationIndices, this->AnnihilationIndices, this->NbrNBody, Coefficient);
      vDestination[Index] = vSource[i] * Coefficient;
    }
  return vDestination;
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnSphereSquareTotalMomentumOperator::Multiply(ComplexVector& vSource, ComplexVector& vDestination, 
								     int firstComponent, int nbrComponent)
{
  return vDestination;
}


