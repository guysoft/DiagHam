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
// space = pointer to the Hilbert space to use
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// nextCoordinates = indicate which coordinates will be change during next time step (-1 if no time coherence has to be used)

QHEParticleWaveFunctionOperation::QHEParticleWaveFunctionOperation (ParticleOnSphere* space, RealVector* state, RealVector* position, 
								    AbstractFunctionBasis* basis, int nextCoordinates)
{
  this->FirstComponent = 0;
  this->NbrComponent = space->GetHilbertSpaceDimension();
  this->NextCoordinates = nextCoordinates;
  if (this->NextCoordinates != -1)
    space->InitializeWaveFunctionEvaluation(true);
  else
    space->InitializeWaveFunctionEvaluation(false);
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
  this->NextCoordinates = operation.NextCoordinates;
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
  if (this->NextCoordinates == -1)
    this->Scalar = this->HilbertSpace->EvaluateWaveFunction(*(this->State), *(this->Position), *(this->Basis),
							    this->FirstComponent, this->NbrComponent);
  else
    this->Scalar = this->HilbertSpace->EvaluateWaveFunctionWithTimeCoherence(*(this->State), *(this->Position), 
									     *(this->Basis), this->NextCoordinates, 
									     this->FirstComponent, this->NbrComponent);
  return true;
}

