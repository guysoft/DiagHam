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


#ifndef QHEPARTICLEWAVEFUNCTIONOPERATION_H
#define QHEPARTICLEWAVEFUNCTIONOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractScalarSumOperation.h"
#include "HilbertSpace/QHEHilbertSpace/ParticleOnSphere.h"


class AbstractFunctionBasis;
class RealVector;


class QHEParticleWaveFunctionOperation: public AbstractScalarSumOperation
{

 protected:

  // pointer to the Hilbert space
  ParticleOnSphere* HilbertSpace;
  // vector corresponding to the state in the Fock basis  
  RealVector* State;
  // vector whose components give coordinates of the point where the wave function has to be evaluated
  RealVector* Position;
  // one body real space basis to use
  AbstractFunctionBasis* Basis;


 public:
  
  // constructor 
  //
  // space = pointer to the Hilbert space to use
  // state = vector corresponding to the state in the Fock basis
  // position = vector whose components give coordinates of the point where the wave function has to be evaluated
  // basis = one body real space basis to use
  QHEParticleWaveFunctionOperation(ParticleOnSphere* space, RealVector* state, RealVector* position, AbstractFunctionBasis* basis);

  // copy constructor 
  //
  // operation = reference on operation to copy
  QHEParticleWaveFunctionOperation(const QHEParticleWaveFunctionOperation& operation);
  
  // destructor
  //
  ~QHEParticleWaveFunctionOperation();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // get dimension (i.e. Hilbert space dimension)
  //
  // return value = dimension
  int GetDimension ();

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // apply operation
  //
  // return value = true if no error occurs
  bool ApplyOperation();
  
};

// get dimension (i.e. Hilbert space dimension)
//
// return value = dimension

inline int QHEParticleWaveFunctionOperation::GetDimension ()
{
  return this->HilbertSpace->GetHilbertSpaceDimension();
}

#endif
