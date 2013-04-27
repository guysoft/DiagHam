////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                       class of CFT evaluation for MPS	              //
//                                                                            //
//                        last modification : 22/04/2013                      //
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


#ifndef FHQEMPSEVALUATECFTOPERATION_H
#define FHQEMPSEVALUATECFTOPERATION_H


#include "config.h"
#include "Vector/RealVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"

#include "Matrix/LongRationalMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"



class RealVector;
class SparseRealMatrix;
class SparseComplexMatrix;
class FQHEMPSClustered2RMatrix;


class FQHEMPSEvaluateCFTOperation: public AbstractArchitectureOperation
{

 protected:

  // pointer to the MPS matrix 
  FQHEMPSClustered2RMatrix* MPSMatrix;
  
  // indicates if the overlap matrix has to be computed instead of the matrix elements
  bool OverlapMatrix;

  // level for the left state
  int LeftLevel;
  // level for the right state
  int RightLevel;

  // array that contains the Hilbert space for the partition at each level
  BosonOnDiskShort** U1BosonBasis;

  
  // array where the already computed overlap matrices are stored
  LongRationalMatrix* PreviousOverlapMatrices;
  // array where the already computed matrix element are stored
  LongRationalMatrix** PreviousMatrixElements;
  // number of entry of the PreviousMatrixElements first index
  int NbrLeftPreviousMatrixElements;
  // number of entry of the PreviousMatrixElements second index
  int NbrRightPreviousMatrixElements;
  
  // value of the central charge divided by 12
  LongRational CentralCharge12;
  // conformal weight of the left state at level 0
  LongRational WeightLeft;
  // conformal weight of the right state  at level 0
  LongRational WeightRight;
  // conformal weight of the field whose matrix elements have to be evaluated
  LongRational WeightMatrixElement;

  // index of the first component
  long FirstComponent;
  // number of component to compute 
  long NbrComponent;

  // matrix storing the CFT matrix elements
  LongRationalMatrix MatrixElements; 

 // number of part in the CFT calculation will be separated in MPI mode
  int NbrMPIStage;
  // number of part in the CFT calculation will be separated in SMP mode
  int NbrSMPStage;
  // array with size of SMP stages used to distribute work
  int* SMPStages; 

 public:
  
  // constructor 
  //
  // mPSMatrix = pointer to the MPS matrix 
  // u1BosonBasis = array that contains the Hilbert space for the partition at each level
  // leftLevel = level for the left state
  // rightLevel = level for the right state
  // centralCharge12 =value of the central charge divided by 12
  // weightLeft = conformal weight of the left state at level 0 
  // weightRight = weight of the right state  at level 0
  // WeightMatrixElement = weight of the field whose matrix elements have to be evaluated
  // previousMatrixElements = array where the already computed matrix element are stored
  // nbrLeftPreviousMatrixElements = number of entry of the PreviousMatrixElements first index
  // nbrRightPreviousMatrixElements = number of entry of the PreviousMatrixElements second index
  // nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
  // nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode
  FQHEMPSEvaluateCFTOperation(FQHEMPSClustered2RMatrix* mPSMatrix, BosonOnDiskShort** u1BosonBasis, int leftLevel, int rightLevel,
			      const LongRational& centralCharge12, const LongRational& weightLeft, const LongRational& weightRight, const LongRational& weightMatrixElement,
			      LongRationalMatrix** previousMatrixElements, int nbrLeftPreviousMatrixElements, int nbrRightPreviousMatrixElements,
			      int nbrMPIStage = 100, int nbrSMPStage = 100);
  
  // constructor to compute the CFT overlap matrix
  //
  // mPSMatrix = pointer to the MPS matrix 
  // u1BosonBasis = array that contains the Hilbert space for the partition at each level
  // leftLevel = level for the left or right state
  // centralCharge12 =value of the central charge divided by 12
  // weightLeft = conformal weight of the left or right state at level 0 
  // previousOverlapMatrices = array where the already computed overlap matrices are stored
  // nbrPreviousOverlapMatrices = number of entry of the PreviousMatrixElements first index
  // nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
  // nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode
  FQHEMPSEvaluateCFTOperation(FQHEMPSClustered2RMatrix* mPSMatrix, BosonOnDiskShort** u1BosonBasis, int leftLevel,
			      const LongRational& centralCharge12, const LongRational& weightLeft,
			      LongRationalMatrix* previousOverlapMatrices, int nbrPreviousOverlapMatrices,
			      int nbrMPIStage = 100, int nbrSMPStage = 100);

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHEMPSEvaluateCFTOperation(const FQHEMPSEvaluateCFTOperation & operation);
  
  // constructor from a master node information
  //
  // space= pointer to the HilbertSpace to use
  // architecture = pointer to the distributed architecture to use for communications
  //  FQHEMPSEvaluateCFTOperation(ParticleOnSphere* space,  SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~FQHEMPSEvaluateCFTOperation();  
  
  // get the matrix storing the CFT matrix elements
  // 
  // return value = reference on the matrix
  virtual LongRationalMatrix& GetMatrixElements ();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  virtual void SetIndicesRange (const long& firstComponent, const long& nbrComponent);

  // clone operation
  //
  // return value = pointer to cloned operation
  virtual AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  virtual bool RawApplyOperation();
  
  // apply operation for SMP using round robin scheduling
  //
  //  architecture = instance of architecture class
  // return value = true if no error occurs
  virtual bool ApplyOperationSMPRoundRobin(SMPArchitecture* architecture, int threadID);

 protected:
  
  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // apply operation for SimpleMPI architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture);
  

};

// get the matrix storing the CFT matrix elements
// 
// return value = reference on the matrix

inline LongRationalMatrix& FQHEMPSEvaluateCFTOperation::GetMatrixElements ()
{
  return this->MatrixElements;
}

#endif
