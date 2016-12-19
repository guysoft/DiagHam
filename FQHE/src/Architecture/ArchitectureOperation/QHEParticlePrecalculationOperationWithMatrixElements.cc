////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//          class of QHE particle hamiltonian precalculation operation        //
//                                                                            //
//                        last modification : 11/03/2003                      //
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
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperationWithMatrixElements.h"

#include <sys/time.h>

// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// firstPass = flag to indicate if the operation has to be applied to the first pass of the precalculations

QHEParticlePrecalculationOperationWithMatrixElements::QHEParticlePrecalculationOperationWithMatrixElements (AbstractQHEHamiltonian* hamiltonian, bool firstPass) :
  QHEParticlePrecalculationOperation (hamiltonian, firstPass)
{
  this->OperationType = AbstractArchitectureOperation::QHEParticlePrecalculationWithMatrixElements;
}

// copy constructor 
//
// operation = reference on operation to copy

QHEParticlePrecalculationOperationWithMatrixElements::QHEParticlePrecalculationOperationWithMatrixElements(const QHEParticlePrecalculationOperationWithMatrixElements& operation):
  QHEParticlePrecalculationOperation(operation)
{
  this->OperationType = AbstractArchitectureOperation::QHEParticlePrecalculationWithMatrixElements;
}
  
// destructor
//

QHEParticlePrecalculationOperationWithMatrixElements::~QHEParticlePrecalculationOperationWithMatrixElements()
{
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void QHEParticlePrecalculationOperationWithMatrixElements::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* QHEParticlePrecalculationOperationWithMatrixElements::Clone()
{
  return new QHEParticlePrecalculationOperationWithMatrixElements (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool QHEParticlePrecalculationOperationWithMatrixElements::RawApplyOperation()
{
  if (this->FirstPass ==  true)
    {
      this->RequiredMemory = this->Hamiltonian->PartialFastMultiplicationMemory(this->FirstComponent, this->NbrComponent, this->RealInteractionCoefficients, this->ComplexInteractionCoefficients);
    }
  else
    {
      this->Hamiltonian->PartialEnableFastMultiplication(this->FirstComponent, this->NbrComponent);
    }
  return true;
}


// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// mpiNodeNbr = provide the additional MPI node ID
// return value = true if no error occurs

bool QHEParticlePrecalculationOperationWithMatrixElements::ArchitectureDependentApplyOperation(SMPArchitecture* architecture, int mpiNodeNbr)
{
  long *SegmentIndices=0;
  int TmpNbrThreads = architecture->GetNbrThreads();
  if (Hamiltonian->GetLoadBalancing(TmpNbrThreads, SegmentIndices)==false)
    {
      SegmentIndices = new long[TmpNbrThreads+1];
      int Step = this->NbrComponent / TmpNbrThreads;
      SegmentIndices[0]=this->FirstComponent;
      for (int i=0; i<TmpNbrThreads; ++i)
	SegmentIndices[i]=this->FirstComponent+i*Step;
      SegmentIndices[TmpNbrThreads]=this->FirstComponent+this->NbrComponent;
    }
  QHEParticlePrecalculationOperationWithMatrixElements** TmpOperations = new QHEParticlePrecalculationOperationWithMatrixElements* [architecture->GetNbrThreads()];
  for (int i = 0; i < TmpNbrThreads; ++i)
    {
      TmpOperations[i] = (QHEParticlePrecalculationOperationWithMatrixElements*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
      TmpOperations[i]->SetIndicesRange(SegmentIndices[i], SegmentIndices[i+1]-SegmentIndices[i]);
    }
  architecture->SendJobs();
  
  if (this->FirstPass ==  true) //  merge matrix elements from threads
    {
      this->RealInteractionCoefficients = TmpOperations[0]->RealInteractionCoefficients;
      this->ComplexInteractionCoefficients = TmpOperations[0]->ComplexInteractionCoefficients;
      cout << "Merging matrix elements - this could be multi-threaded"<<endl;
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      gettimeofday (&(TotalStartingTime), 0);
      int StartTimeSecond = TotalStartingTime.tv_sec;      
      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
	{
	  this->RealInteractionCoefficients.MergeArray(TmpOperations[i]->RealInteractionCoefficients);
	  this->ComplexInteractionCoefficients.MergeArray(TmpOperations[i]->ComplexInteractionCoefficients);
	}
      gettimeofday (&(TotalEndingTime), 0);
      double Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);	          
      this->RealInteractionCoefficients.SortEntries();
      this->ComplexInteractionCoefficients.SortEntries();
      cout << "Done merging arrays in "<<Dt<<"s"<<endl;
    }
  for (int i = 0; i < architecture->GetNbrThreads(); ++i)
    {
      if (this->FirstPass ==  true)
	{
	  if (mpiNodeNbr>=0)
	    cout << "node "<<mpiNodeNbr<<" ";
	  cout << "thread "<<i<<" = "<<TmpOperations[i]->RequiredMemory<<endl;
	}
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  if (Hamiltonian->GetLoadBalancing(TmpNbrThreads, SegmentIndices)==false)
    delete [] SegmentIndices;
  return true;  
}


