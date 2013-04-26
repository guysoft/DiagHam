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

#include "config.h"
#include "Architecture/ArchitectureOperation/FQHEMPSEvaluateCFTOperation.h"
#include "Vector/RealVector.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "Matrix/SparseComplexMatrix.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"



#include <iostream>
#include <sys/time.h>
#include <stdlib.h>

using std::cout;
using std::endl;


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


FQHEMPSEvaluateCFTOperation::FQHEMPSEvaluateCFTOperation(FQHEMPSClustered2RMatrix* mPSMatrix, BosonOnDiskShort** u1BosonBasis, int leftLevel, int rightLevel,
							 const LongRational& centralCharge12, const LongRational& weightLeft, 
							 const LongRational& weightRight, const LongRational& weightMatrixElement,
							 LongRationalMatrix** previousMatrixElements, int nbrLeftPreviousMatrixElements, int nbrRightPreviousMatrixElements,
							 int nbrMPIStage, int nbrSMPStage)
{
  this->MPSMatrix = mPSMatrix;
  this->U1BosonBasis = u1BosonBasis;
  this->LeftLevel = leftLevel;
  this->RightLevel = rightLevel;
  this->CentralCharge12 = centralCharge12;
  this->WeightLeft = weightLeft;
  this->WeightRight = weightRight;
  this->WeightMatrixElement = weightMatrixElement;
  this->PreviousMatrixElements = previousMatrixElements;
  this->NbrLeftPreviousMatrixElements = nbrLeftPreviousMatrixElements;
  this->NbrRightPreviousMatrixElements = nbrRightPreviousMatrixElements;
  this->FirstComponent = 0;
  this->NbrComponent = this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() * this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
  this->MatrixElements = LongRationalMatrix (this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension(), U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension(), true);
  this->OperationType = AbstractArchitectureOperation::FQHEMPSEvaluateCFTOperation;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
  this->SMPStages = new int[1]; 
}

// copy constructor 
//
// operation = reference on operation to copy

FQHEMPSEvaluateCFTOperation::FQHEMPSEvaluateCFTOperation(const FQHEMPSEvaluateCFTOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;

  this->MPSMatrix = operation.MPSMatrix;
  this->U1BosonBasis = operation.U1BosonBasis;
  this->LeftLevel = operation.LeftLevel;
  this->RightLevel = operation.RightLevel;
  this->CentralCharge12 = operation.CentralCharge12;
  this->WeightLeft = operation.WeightLeft;
  this->WeightRight = operation.WeightRight;
  this->WeightMatrixElement = operation.WeightMatrixElement;
  this->PreviousMatrixElements = operation.PreviousMatrixElements;
  this->NbrLeftPreviousMatrixElements = operation.NbrLeftPreviousMatrixElements;
  this->NbrRightPreviousMatrixElements = operation.NbrRightPreviousMatrixElements;
  this->MatrixElements = operation.MatrixElements;
  this->OperationType = AbstractArchitectureOperation::FQHEMPSEvaluateCFTOperation;
  this->NbrMPIStage = operation.NbrMPIStage;
  this->NbrSMPStage = operation.NbrSMPStage;
  this->SMPStages = operation.SMPStages;
}

// destructor
//

FQHEMPSEvaluateCFTOperation::~FQHEMPSEvaluateCFTOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHEMPSEvaluateCFTOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHEMPSEvaluateCFTOperation::Clone()
{
  return new FQHEMPSEvaluateCFTOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHEMPSEvaluateCFTOperation::RawApplyOperation()
{
  if (this->NbrComponent == 0)
    return true;
  int LastComponent = this->NbrComponent + this->FirstComponent;
  long* Partition = new long[this->U1BosonBasis[this->LeftLevel]->GetHilbertSpaceDimension() + this->U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension() + 1];
  unsigned long* TmpPartition = new unsigned long [this->MPSMatrix->GetTruncationLevel() + 2];
  unsigned long* TemporaryOccupationNumber = new unsigned long [this->MPSMatrix->GetTruncationLevel() + 2];
  for (int Index = this->FirstComponent; Index < LastComponent; ++Index)
    {
      int n = Index / U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
      int m = Index % U1BosonBasis[this->RightLevel]->GetHilbertSpaceDimension();
      int PartitionLength = 0;
      U1BosonBasis[this->LeftLevel]->GetOccupationNumber(n, TmpPartition);	    
      for (int k = 1; k <= this->LeftLevel; ++k)
	for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
	  {
	    ++PartitionLength;
	  }
      int Position = PartitionLength;
      PartitionLength = 0;
      for (int k = 1; k <= this->LeftLevel; ++k)
	for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
	  {
	    Partition[Position - PartitionLength - 1] = (long) k;
	    ++PartitionLength;
	  }
      U1BosonBasis[this->RightLevel]->GetOccupationNumber(m, TmpPartition);	    
      for (int k = 1; k <= this->RightLevel; ++k)
	for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
	  {
	    Partition[PartitionLength] = -(long) k;
	    ++PartitionLength;		  
	  }
      LongRational Tmp = this->MPSMatrix->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, this->CentralCharge12, 
									  this->WeightLeft, this->WeightRight, this->WeightMatrixElement,
									  this->PreviousMatrixElements, this->NbrLeftPreviousMatrixElements, this->NbrRightPreviousMatrixElements, 
									  this->U1BosonBasis, TemporaryOccupationNumber);
      this->MatrixElements.SetMatrixElement(n, m, Tmp);
    }
  delete[] TmpPartition;
  delete[] Partition;
  return true;
}



// apply operation for SMP using round robin scheduling
//
//  architecture = instance of architecture class
// return value = true if no error occurs

bool FQHEMPSEvaluateCFTOperation::ApplyOperationSMPRoundRobin(SMPArchitecture* architecture, int threadID)
{      
  int NbrComponents = this->NbrComponent;
  int FirstComponent = this->FirstComponent;
  
   int NbrStages = this->NbrSMPStage * architecture->GetNbrThreads();
   int StageIndex = 0;

   bool LockFlag = false;
   while (StageIndex < NbrStages) 
     {
       if (LockFlag == false)   
         {
           architecture->LockMutex();
           LockFlag = true;
         }
       StageIndex = this->SMPStages[0];
       if (StageIndex < NbrStages) 
         {         
           this->SMPStages[0]++;   
           architecture->UnLockMutex();
           LockFlag = false;
           this->SetIndicesRange(FirstComponent + this->GetRankChunkStart(NbrComponents, StageIndex,  NbrStages),  
				 this->GetRankChunkSize(NbrComponents, StageIndex,  NbrStages));
           this->RawApplyOperation();
           ++StageIndex;
         }
     }
   if (LockFlag == true)
     {
       architecture->UnLockMutex();
       LockFlag = false;
     }
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHEMPSEvaluateCFTOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent / architecture->GetNbrThreads();
  if (Step == 0)
    Step = this->NbrComponent;
  this->SMPStages[0] = 0;
  int TotalNbrComponent = this->FirstComponent + this->NbrComponent;
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHEMPSEvaluateCFTOperation** TmpOperations = new FQHEMPSEvaluateCFTOperation * [architecture->GetNbrThreads()];
  
   
  for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
    {
      TmpOperations[i] = (FQHEMPSEvaluateCFTOperation *) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  
  architecture->SendJobsRoundRobin();

//   for( int i = 0; i < architecture->GetNbrThreads(); i++)
//     {
//       TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
//       TmpFirstComponent += Step;
//       if (TmpFirstComponent >= TotalNbrComponent)
// 	{
// 	  Step = 0;
// 	  TmpFirstComponent = TotalNbrComponent;
// 	}
//       else
// 	{
// 	  if ((TmpFirstComponent + Step) >= TotalNbrComponent)
// 	    {
// 	      Step = TotalNbrComponent - TmpFirstComponent;	      
// 	    }
// 	}
//     }
//  architecture->SendJobs();

  for (int i = 0; i < architecture->GetNbrThreads(); i++)
    {
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  return true;
}

// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHEMPSEvaluateCFTOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__    
//   int Step = this->NbrComponent / architecture->GetNbrNodes();
//   int TmpFirstComponent = this->FirstComponent + (Step * architecture->GetNodeNbr());
//   int TmpNbrComponent = Step;
//   if ((architecture->GetNodeNbr() + 1) == architecture->GetNbrNodes())
//     {
//       TmpNbrComponent += this->NbrComponent % Step;
//     }
//   this->SetIndicesRange(TmpFirstComponent, TmpNbrComponent); 
//   switch (architecture->GetArchitectureID())
//     {	 
//     case AbstractArchitecture::MixedMPISMP:
//       this->ArchitectureDependentApplyOperation((SMPArchitecture*) (architecture->GetLocalArchitecture())); 
//       break;
//     default:
//       this->RawApplyOperation();
//       break;
//     }		
//   MPI::COMM_WORLD.Barrier();
//   if (this->OutputState != 0)
//     architecture->SumVector(*(this->OutputState));	
//   else
//     architecture->SumVector(*(this->ComplexOutputState));	
  return true;
#else
  return this->RawApplyOperation();
#endif
}
