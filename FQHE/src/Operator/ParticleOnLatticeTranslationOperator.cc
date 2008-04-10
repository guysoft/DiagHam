////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of particle on lattice 1-body operator              //
//                                                                            //
//                        last modification : 09/04/2008                      //
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
#include "Operator/ParticleOnLatticeTranslationOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

  
// constructor from default datas
//
// particle = hilbert space associated to the particles
// rx = x-component of the desired translation
// ry = y-component of the desired translation
ParticleOnLatticeTranslationOperator::ParticleOnLatticeTranslationOperator(ParticleOnLattice* particle, int rx, int ry)
{
  this->Particle = (ParticleOnLattice*) (particle->Clone());
  this->Rx=rx;
  this->Ry=ry;  
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnLatticeTranslationOperator::ParticleOnLatticeTranslationOperator(const ParticleOnLatticeTranslationOperator& oper)
{
  this->Particle = (ParticleOnLattice*) (oper.Particle->Clone());
  this->Rx=oper.Rx;
  this->Ry=oper.Ry;
}

// destructor
//

ParticleOnLatticeTranslationOperator::~ParticleOnLatticeTranslationOperator()
{
  delete this->Particle;
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnLatticeTranslationOperator::Clone ()
{
  return new ParticleOnLatticeTranslationOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeTranslationOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnLattice*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeTranslationOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeTranslationOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}

// set components of translation vector
//
// rx = x-component of the desired translation
// ry = y-component of the desired translation
void ParticleOnLatticeTranslationOperator::SetTranslationComponents(int rx, int ry)
{
  this->Rx=rx;
  this->Ry=ry;  
}



  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnLatticeTranslationOperator::MatrixElement (RealVector& V1, RealVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  Complex TranslationPhase;  
  Complex Element = 0.0;
  int Index;
  for (int i = 0; i < Dim; ++i)
    {
      Index = this->Particle->TranslateState(i, this->Rx, this->Ry, TranslationPhase);      
      Element += V1[Index] * V2[i] * TranslationPhase;
    }
  return Element;
}


Complex ParticleOnLatticeTranslationOperator::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  Complex TranslationPhase;
  Complex Element;
  int Index;
  for (int i = 0; i < Dim; ++i)
    {
      Index = this->Particle->TranslateState(i, this->Rx, this->Ry, TranslationPhase);
      Element += Conj(V1[Index]) * V2[i] * TranslationPhase;
    }
  return Element; 
}
   
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnLatticeTranslationOperator::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Index;
  Complex TranslationPhase;
  for (int i = firstComponent; i < Last; ++i)
    {
      Index = this->Particle->TranslateState(i, this->Rx, this->Ry, TranslationPhase);
      vDestination[Index].Re = vSource[i].Re * TranslationPhase.Re - vSource[i].Im * TranslationPhase.Im;
      vDestination[Index].Im = vSource[i].Re * TranslationPhase.Im + vSource[i].Im * TranslationPhase.Re;
    }
  return vDestination;
}


