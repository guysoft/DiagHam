////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of QHE particle wave function evaluation operation         //
//                                                                            //
//                        last modification : 29/07/2004                      //
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
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"


// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// firstPass = flag to indicate if the operation has to be applied to the first pass of the precalculations

QHEParticleWaveFunctionOperation::QHEParticleWaveFunctionOperation (ParticleOnSphere* space, RealVector* state, RealVector* position, 
								    AbstractFunctionBasis* basis)
{
  this->FirstComponent = 0;
  this->NbrComponent = space->GetHilbertSpaceDimension();
  this->HilbertSpace = (ParticleOnSphere*) space->Clone();
  this->State = state;
  this->Position = position;
  this->OperationType = AbstractArchitectureOperation::QHEParticleWaveFunction;
  this->Basis = basis;
}

// copy constructor 
//
// operation = reference on operation to copy

QHEParticleWaveFunctionOperation::QHEParticleWaveFunctionOperation(const QHEParticleWaveFunctionOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->State = operation.State;
  this->HilbertSpace = (ParticleOnSphere*) operation.HilbertSpace->Clone();
  this->Position = operation.Position;
  this->OperationType = AbstractArchitectureOperation::QHEParticleWaveFunction;
  this->Basis = operation.Basis;
  this->Scalar = operation.Scalar;
}
  
// destructor
//

QHEParticleWaveFunctionOperation::~QHEParticleWaveFunctionOperation()
{
  delete this->HilbertSpace;
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void QHEParticleWaveFunctionOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* QHEParticleWaveFunctionOperation::Clone()
{
  return new QHEParticleWaveFunctionOperation (*this);
}
  
// apply operation
//
// return value = true if no error occurs

bool QHEParticleWaveFunctionOperation::ApplyOperation()
{
  this->Scalar = this->HilbertSpace->EvaluateWaveFunction(*(this->State), *(this->Position), *(this->Basis),
							  this->FirstComponent, this->NbrComponent);
  return true;
}

