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
#include "Architecture/ArchitectureOperation/FQHEMPSCreateStateOperation.h"
#include "Vector/RealVector.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "Matrix/SparseComplexMatrix.h"



#include <iostream>
#include <sys/time.h>
#include <stdlib.h>

using std::cout;
using std::endl;


// constructor 
//
// space = pointer to the Hilbert space
// bMatrices = array that gives the B matrices 
// state = pointer to the vector where the MPS state will be stored
// traceFlag = indicates the type of boundary conditions (-1 = trace, traceFlag >= 0 takes the final corresponding diagonal element)
// blockSize = indicates the size of the block for precalculations

FQHEMPSCreateStateOperation::FQHEMPSCreateStateOperation(ParticleOnSphere* space, SparseRealMatrix* bMatrices, RealVector* state, int traceFlag, int blockSize)
{
  this->FirstComponent = 0;
  this->NbrComponent = space->GetHilbertSpaceDimension();
  this->Space = (ParticleOnSphere*) space->Clone();
  this->OutputState = state;
  this->BMatrices = bMatrices;
  this->TraceFlag = traceFlag;  
  this->PrecalculationBlockSize = blockSize;
  this->OperationType = AbstractArchitectureOperation::FQHEMPSCreateStateOperation;
}


// copy constructor 
//
// operation = reference on operation to copy

FQHEMPSCreateStateOperation::FQHEMPSCreateStateOperation(const FQHEMPSCreateStateOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;

  this->Space = (ParticleOnSphere*) operation.Space->Clone();
  this->OutputState = operation.OutputState;
  this->BMatrices = operation.BMatrices;
  this->TraceFlag = operation.TraceFlag;  
  this->PrecalculationBlockSize = operation.PrecalculationBlockSize;
  this->OperationType = AbstractArchitectureOperation::FQHEMPSCreateStateOperation;	
}

// destructor
//

FQHEMPSCreateStateOperation::~FQHEMPSCreateStateOperation()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHEMPSCreateStateOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set the output state 
// 
// state = pointer to the output state

void FQHEMPSCreateStateOperation::SetOutputState (RealVector* state)
{
  this->OutputState = state;
}


// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHEMPSCreateStateOperation::Clone()
{
  return new FQHEMPSCreateStateOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHEMPSCreateStateOperation::RawApplyOperation()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);

  this->Space->CreateStateFromMPSDescription(this->BMatrices, *(this->OutputState), this->TraceFlag, (long) this->PrecalculationBlockSize, this->FirstComponent, this->NbrComponent);
  
  timeval TotalEndingTime;
  gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout << this->FirstComponent << " " <<  this->NbrComponent << " : " << Dt << "s" << endl;
  return true;
}



// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHEMPSCreateStateOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent / architecture->GetNbrThreads();
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHEMPSCreateStateOperation** TmpOperations = new FQHEMPSCreateStateOperation * [architecture->GetNbrThreads()];
  
   
  for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
    {
      TmpOperations[i] = (FQHEMPSCreateStateOperation *) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  
  for( int i = 0; i <  ReducedNbrThreads ; i++)
    {
      TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      TmpFirstComponent += Step;
    }

  TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->FirstComponent + this->NbrComponent - TmpFirstComponent);
 
  architecture->SendJobs();

  for (int i = 1; i < architecture->GetNbrThreads(); i++)
    {
      delete TmpOperations[i];
    }

  delete TmpOperations[0];
  delete[] TmpOperations;
  return true;
}

// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHEMPSCreateStateOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__    
  int Step = this->NbrComponent / architecture->GetNbrNodes();
  int TmpFirstComponent = this->FirstComponent + (Step * architecture->GetNodeNbr());
  int TmpNbrComponent = Step;
  if ((architecture->GetNodeNbr() + 1) == architecture->GetNbrNodes())
    {
      TmpNbrComponent += this->NbrComponent % Step;
    }
  this->SetIndicesRange(TmpFirstComponent, TmpNbrComponent); 
  switch (architecture->GetArchitectureID())
    {	 
    case AbstractArchitecture::MixedMPISMP:
      this->ArchitectureDependentApplyOperation((SMPArchitecture*) (architecture->GetLocalArchitecture())); 
      break;
    default:
      this->RawApplyOperation();
      break;
    }		
  MPI::COMM_WORLD.Barrier();
  architecture->SumVector(*(this->OutputState));	
      
  return true;
#else
  return this->RawApplyOperation();
#endif
}
