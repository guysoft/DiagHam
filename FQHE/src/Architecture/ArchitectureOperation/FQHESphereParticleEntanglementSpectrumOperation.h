////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of FQHE particle entanglement spectrum parallelization operation    //
//                                                                            //
//                        last modification : 15/12/2010                      //
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


#ifndef FQHESPHEREPARTICLEENTANGLEMENTSPECTRUMOPERATION_H
#define FQHESPHEREPARTICLEENTANGLEMENTSPECTRUMOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"
#include "HilbertSpace/ParticleOnSphere.h"


class AbstractFunctionBasis;
class RealVector;
class RealSymmetricMatrix;

class FQHESphereParticleEntanglementSpectrumOperation: public AbstractPrecalculationOperation
{

 protected:

  // fullSpace = pointer to the full Hilbert space to use
  ParticleOnSphere* FullSpace;
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  ParticleOnSphere* DestinationHilbertSpace;
  // pointer to the complementary Hilbert space (i.e. part B)
  ParticleOnSphere* ComplementaryHilbertSpace;
  // total system ground state
  RealVector GroundState;
  // reduced density matrix where result is stored
  RealSymmetricMatrix DensityMatrix;


  // a temporary array to store copies of operations in SMP mode
  FQHESphereParticleEntanglementSpectrumOperation** LocalOperations;
  // number of operation copies
  int NbrLocalOperations;

 public:
  
  // constructor 
  //
  // fullSpace = pointer to the full Hilbert space to use
  // destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
  // complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
  // groundState = reference on the total system ground state
  // densityMatrix = reference on the density matrix where result has to stored
  FQHESphereParticleEntanglementSpectrumOperation(ParticleOnSphere* fullSpace, ParticleOnSphere* destinationSpace, ParticleOnSphere* complementarySpace, RealVector& groundState, RealSymmetricMatrix& densityMatrix);

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereParticleEntanglementSpectrumOperation(const FQHESphereParticleEntanglementSpectrumOperation& operation);
  
  // destructor
  //
  ~FQHESphereParticleEntanglementSpectrumOperation();
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // get hilbert space dimension
  // 
  // return value = hilbert space dimension  
  int GetHilbertSpaceDimension ();

 protected:

  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
};

// get hilbert space dimension
// 
// return value = hilbert space dimension

inline int FQHESphereParticleEntanglementSpectrumOperation::GetHilbertSpaceDimension ()
{
  return this->HilbertSpace->GetHilbertSpaceDimension();
}

#endif
