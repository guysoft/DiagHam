////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                 class of particle on sphere density operator               //
//                                                                            //
//                        last modification : 11/12/2002                      //
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
#include "Operator/QHEOperator/ParticleOnSphereDensityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

  
using std::cout;
using std::endl;


// constructor from default datas
//
// particle = hilbert space associated to the particles
// index = index of the density operator

ParticleOnSphereDensityOperator::ParticleOnSphereDensityOperator(ParticleOnSphere* particle, int index)
{
  this->Particle= particle;
  this->OperatorIndex = index;
}

// destructor
//

ParticleOnSphereDensityOperator::~ParticleOnSphereDensityOperator()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereDensityOperator::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphere*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereDensityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereDensityOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereDensityOperator::MatrixElement (RealVector& V1, RealVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Element = 0.0;
  for (int i = 0; i < Dim; ++i)
    {
      Element += V1[i] * V2[i] * this->Particle->AdA(i, this->OperatorIndex);
    }
  return Complex(Element);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereDensityOperator::MatrixElement (ComplexVector& V1, ComplexVector& V2)
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

RealVector& ParticleOnSphereDensityOperator::Multiply(RealVector& vSource, RealVector& vDestination, 
						      int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  for (int i = firstComponent; i < Last; ++i)
    {
      vDestination[i] = vSource[i] * this->Particle->AdA(i, this->OperatorIndex);
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

ComplexVector& ParticleOnSphereDensityOperator::Multiply(ComplexVector& vSource, ComplexVector& vDestination, 
								int firstComponent, int nbrComponent)
{
  return vDestination;
}

