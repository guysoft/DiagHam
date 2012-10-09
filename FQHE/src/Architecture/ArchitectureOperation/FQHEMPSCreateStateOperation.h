////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                       class of state creation from a MPS	              //
//                                                                            //
//                        last modification : 08/10/2012                      //
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


#ifndef FHQEMPSCREATESTATEOPERATION_H
#define FHQEMPSCREATESTATEOPERATION_H


#include "config.h"
#include "Vector/RealVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/FermionOnSpherePTruncated.h"



class RealVector;
class SparseComplexMatrix;


class FQHEMPSCreateStateOperation: public AbstractArchitectureOperation
{

 protected:
  
  // pointer to the Hilbert space
  FermionOnSpherePTruncated* Space;

  // vector where the MPS state will be stored
  RealVector* OutputState;

  // array that gives the B matrices 
  SparseComplexMatrix* BMatrices;

  // indicates the type of boundary conditions (-1 = trace, traceFlag >= 0 takes the final corresponding diagonal element)
  int TraceFlag;  

  // index of the first component
  long FirstComponent;
  // number of component to compute 
  long NbrComponent;
  
 public:
  
  // constructor 
  //
  // space = pointer to the Hilbert space
  // bMatrices = array that gives the B matrices 
  // state = pointer to the vector where the MPS state will be stored
  // traceFlag = indicates the type of boundary conditions (-1 = trace, traceFlag >= 0 takes the final corresponding diagonal element)
  FQHEMPSCreateStateOperation(FermionOnSpherePTruncated* space, SparseComplexMatrix* bMatrices, RealVector* state, int traceFlag);
  
  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHEMPSCreateStateOperation(const FQHEMPSCreateStateOperation & operation);
  
  // constructor from a master node information
  //
  // space= pointer to the HilbertSpace to use
  // architecture = pointer to the distributed architecture to use for communications
  FQHEMPSCreateStateOperation(FermionOnSpherePTruncated* space,  SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~FQHEMPSCreateStateOperation();  
  
  // set the output state 
  // 
  // state = pointer to the output state
  void SetOutputState (RealVector* state);
  
  // get the output state 
  // 
  // return value = pointer to the output state 
  Vector* GetOutputState ();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  virtual void SetIndicesRange (const long& firstComponent, const long& nbrComponent);

  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
 protected:
  
  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  
};

// get the output state 
// 
// return value = pointer to the output state 

inline Vector* FQHEMPSCreateStateOperation::GetOutputState ()
{
  return this->OutputState;
}

#endif
