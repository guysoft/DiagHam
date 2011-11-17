////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Antoine Sterdyniak                //
//                                                                            //
//                                                                            //
//                   class of U1U1 states symmetrization Operation	          //
//                                                                            //
//                        last modification : 03/03/2010                      //
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


#ifndef SYMMETRIZEU1U1OPERATION_H
#define SYMMETRIZEU1U1OPERATION_H


#include "config.h"
#include "Vector/RealVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/BosonOnSphereShort.h"



class RealVector;


class FQHESphereSymmetrizeU1U1StateOperation: public AbstractArchitectureOperation
{

 protected:
  
  // pointer to the HilbertSpace
  BosonOnSphereShort * FinalSpace;
	
  BosonOnSphereShort * LeftSpace;
	
  BosonOnSphereShort * RightSpace;
  
  // vector where the component are stored
  RealVector* LeftVector;
  
  // vector where the component are stored
  RealVector* RightVector;
  
  // vector where the result has to be stored
  RealVector* DestinationVector;
  
  // index of the first component
  long FirstComponent;
  
  long NbrComponent;
  
  bool UnnormalizedBasisFlag;
  
 public:
  
  // constructor 
  //
  // Space = pointer to the HilbertSpace to use
  // sourceVector = vector where the component are stored
  // destinationVector = vector where the result has to be stored
  FQHESphereSymmetrizeU1U1StateOperation( BosonOnSphereShort * finalSpace, BosonOnSphereShort * leftSpace , BosonOnSphereShort * rightSpace , RealVector* destinationVector, RealVector* leftVector, RealVector* rightVector,  bool unnormalizedBasisFlag);
  
  
  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereSymmetrizeU1U1StateOperation(const FQHESphereSymmetrizeU1U1StateOperation & operation);
  
  // constructor from a master node information
  //
  // Space= pointer to the HilbertSpace to use
  // architecture = pointer to the distributed architecture to use for communications
  FQHESphereSymmetrizeU1U1StateOperation( BosonOnSphereShort* space,  SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~FQHESphereSymmetrizeU1U1StateOperation();
  
  
  // set destination vector 
  // 
  // vector where the result has to be stored
  void SetDestinationVector (RealVector* DestinationVector);
  
  // get destination vector 
  // 
  // return value = pointer to destination vector
  Vector* GetDestinationVector ();
  
  //set FirstCoimponent
  //
  // firstComponent  
  void SetIndicesRange (const long& firstComponent, const long& nbrComponent);
  
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

// get destination vector 
// 
// return value = pointer to destination vector

inline Vector* FQHESphereSymmetrizeU1U1StateOperation::GetDestinationVector()
{
  return this->DestinationVector;
}

#endif
