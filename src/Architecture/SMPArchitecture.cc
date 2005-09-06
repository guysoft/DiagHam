////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                          class of SMP Architecture                         //
//                                                                            //
//                        last modification : 30/04/2002                      //
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
#include "Architecture/SMPArchitecture.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/MultipleVectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MatrixMatrixMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"
#include "Architecture/ArchitectureOperation/AbstractScalarSumOperation.h"

#ifdef __SMP__
#include <pthread.h>
#endif
#include <stdlib.h>
#include <iostream>


using std::cout;
using std::endl;



// thread for multiplicationmain function
//
// param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
// return value = unused pointer (null)
void* ThreadMain(void* param);

// thread for multiplicationmain function
//
// param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
// return value = unused pointer (null)
void* ThreadExecuteOperation(void* param);


// constructor
//
// nbrThreads = number of threads to run simultaneously (in principle, the number of processors that can be allocated)

SMPArchitecture::SMPArchitecture(int nbrThreads)
{
  this->ArchitectureID = AbstractArchitecture::SMP;
  this->NbrThreads = nbrThreads;
#ifdef __SMP__
  this->ThreadParameters = new ThreadMainParameter [this->NbrThreads];
  this->Threads = new pthread_t [this->NbrThreads];
  for (int i = 0; i < this->NbrThreads; ++i)
    {
      this->ThreadParameters[i].ThreadState = SMPArchitecture::Wait;
      this->ThreadParameters[i].ThreadID = i;
 /*     if (pthread_create (&(this->Threads[i]), 0, ThreadMain, (void*) &(this->ThreadParameters[i])))
	{
	  cout << "error, cannot create thread" << endl;
	  exit(1);
	}*/
    }
#endif
}
  
// destructor
//

SMPArchitecture::~SMPArchitecture()
{
#ifdef __SMP__
/*  for (int i = 0; i < this->NbrThreads; ++i)
    {
      if (this->ThreadParameters[i].ThreadState != SMPArchitecture::Dead)
	{
	  pthread_mutex_lock(this->ThreadParameters[i].LocalMutex);
	  this->ThreadParameters[i].ThreadState = SMPArchitecture::Exit;
	  pthread_mutex_unlock(this->ThreadParameters[i].LocalMutex);
	  pthread_cond_signal(this->ThreadParameters[i].LocalCondition);       
	}
    }
  void* ReturnValue;
  for (int i = 0; i < this->NbrThreads; ++i)
    {
      (void) pthread_join (this->Threads[i], &ReturnValue);
    }*/
#endif
}
  
// set the operation that has to be executed by a given thread
//
// operation = pointer to the operation
// index = thread index

void SMPArchitecture::SetThreadOperation(AbstractArchitectureOperation* operation, int index)
{
  this->ThreadParameters[index].Operation = operation;
}

// execute an architecture-dependent vector hamiltonian multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (VectorHamiltonianMultiplyOperation* operation)
{
  int Step = operation->GetDestinationVector()->GetVectorDimension() / this->NbrThreads;
  int FirstComponent = 0;
  int ReducedNbrThreads = this->NbrThreads - 1;
  operation->GetDestinationVector()->ClearVector();
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((VectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  this->ThreadParameters[ReducedNbrThreads].Operation = operation->Clone();
  ((VectorHamiltonianMultiplyOperation*) (this->ThreadParameters[ReducedNbrThreads].Operation))->SetIndicesRange(FirstComponent, 
												  operation->GetDestinationVector()->GetVectorDimension() - FirstComponent);  
  for (int i = 1; i < this->NbrThreads; ++i)
    {
      ((VectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->SetDestinationVector(operation->GetDestinationVector()->EmptyClone(true));
    }
  this->SendJobs();
  for (int i = 1; i < this->NbrThreads; ++i)
    {
      (*(operation->GetDestinationVector())) += (*(((VectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->GetDestinationVector()));
      delete ((VectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->GetDestinationVector();
      delete this->ThreadParameters[i].Operation;
    }
  return true;
}
  
// execute an architecture-dependent multiple vector hamiltonian multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (MultipleVectorHamiltonianMultiplyOperation* operation)
{
  bool RealFlag = false;
  int Dimension = 0;
  if (operation->GetDestinationComplexVectors() == 0)
    {
      RealFlag = true;
      Dimension = operation->GetDestinationRealVectors()[0].GetVectorDimension();
    }
  else
    {
      Dimension = operation->GetDestinationComplexVectors()[0].GetVectorDimension();
    }
  int Step = Dimension / this->NbrThreads;
  int FirstComponent = 0;
  int ReducedNbrThreads = this->NbrThreads - 1;
  int NbrVectors = operation->GetNbrVectors();
  if (RealFlag == true)
    {
      for (int i = 0; i < NbrVectors; ++i)
	operation->GetDestinationRealVectors()[i].ClearVector();
    }
  else
    {
      for (int i = 0; i < NbrVectors; ++i)
	operation->GetDestinationComplexVectors()[i].ClearVector();
    }
   for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((MultipleVectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  this->ThreadParameters[ReducedNbrThreads].Operation = operation->Clone();
  ((MultipleVectorHamiltonianMultiplyOperation*) (this->ThreadParameters[ReducedNbrThreads].Operation))->SetIndicesRange(FirstComponent, Dimension - FirstComponent);  
  for (int i = 1; i < this->NbrThreads; ++i)
    {
      if (RealFlag == true)
	((VectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->SetDestinationVector(operation->GetDestinationRealVectors()[0].EmptyCloneArray(NbrVectors, 
																				      true));
      else
	((VectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->SetDestinationVector(operation->GetDestinationComplexVectors()[0].EmptyCloneArray(NbrVectors, 
																					 true));

    }
  this->SendJobs();
  for (int i = 1; i < this->NbrThreads; ++i)
    {
      if (RealFlag == true)
	{
	  for (int j = 0; j < NbrVectors; ++j)
	    operation->GetDestinationRealVectors()[j] += ((MultipleVectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->GetDestinationRealVectors()[j];
	  delete[] ((MultipleVectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->GetDestinationRealVectors();
	}
      else
	{
	  for (int j = 0; j < NbrVectors; ++j)
	    operation->GetDestinationComplexVectors()[j] += ((MultipleVectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->GetDestinationComplexVectors()[j];
	  delete[] ((MultipleVectorHamiltonianMultiplyOperation*) (this->ThreadParameters[i].Operation))->GetDestinationComplexVectors();
	}
      delete this->ThreadParameters[i].Operation;
    }
  return true;
}

// execute an architecture-dependent vector abstract scalar sum operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully
  
bool SMPArchitecture::ExecuteOperation (AbstractScalarSumOperation* operation)
{
  int Step = operation->GetDimension() / this->NbrThreads;
  int FirstComponent = 0;
  int ReducedNbrThreads = this->NbrThreads - 1;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((AbstractScalarSumOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  this->ThreadParameters[ReducedNbrThreads].Operation = operation->Clone();
  ((AbstractScalarSumOperation*) (this->ThreadParameters[ReducedNbrThreads].Operation))->SetIndicesRange(FirstComponent,
													   operation->GetDimension() - FirstComponent);
  this->SendJobs();
  for (int i = 0; i < this->NbrThreads; ++i)
    {
      operation->GetScalar() += ((AbstractScalarSumOperation*) (this->ThreadParameters[i].Operation))->GetScalar();
      delete this->ThreadParameters[i].Operation;
    }
  return true;
}

// execute an architecture-dependent add real linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (AddRealLinearCombinationOperation* operation)
{
  int Step = operation->GetDestinationVector()->GetVectorDimension() / this->NbrThreads;
  int FirstComponent = 0;
  int ReducedNbrThreads = this->NbrThreads - 1;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((AddRealLinearCombinationOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  this->ThreadParameters[ReducedNbrThreads].Operation = operation->Clone();
  ((AddRealLinearCombinationOperation*) (this->ThreadParameters[ReducedNbrThreads].Operation))->SetIndicesRange(FirstComponent, 
														  operation->GetDestinationVector()
														  ->GetVectorDimension() - FirstComponent);  
  this->SendJobs();
  for (int i = 0; i < this->NbrThreads; ++i)
    {
      delete this->ThreadParameters[i].Operation;
    }
  return true;
}  

// execute an architecture-dependent add complex linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (AddComplexLinearCombinationOperation* operation)
{
  int Step = operation->GetDestinationVector()->GetVectorDimension() / this->NbrThreads;
  int FirstComponent = 0;
  int ReducedNbrThreads = this->NbrThreads - 1;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((AddComplexLinearCombinationOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  this->ThreadParameters[ReducedNbrThreads].Operation = operation->Clone();
  ((AddComplexLinearCombinationOperation*) (this->ThreadParameters[ReducedNbrThreads].Operation))->SetIndicesRange(FirstComponent, 
														  operation->GetDestinationVector()
														  ->GetVectorDimension() - FirstComponent);  
  this->SendJobs();
  for (int i = 0; i < this->NbrThreads; ++i)
    {
      delete this->ThreadParameters[i].Operation;
    }
  return true;
}  

// execute an architecture-dependent multiple real scalar product operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (MultipleRealScalarProductOperation* operation)
{
  int Step = operation->GetLeftVector()->GetVectorDimension() / this->NbrThreads;
  int FirstComponent = 0;
  int ReducedNbrThreads = this->NbrThreads - 1;
  this->ThreadParameters[0].Operation = operation;
  ((MultipleRealScalarProductOperation*) (this->ThreadParameters[0].Operation))->SetIndicesRange(FirstComponent, Step);
  FirstComponent += Step;
  for (int i = 1; i < ReducedNbrThreads; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((MultipleRealScalarProductOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      ((MultipleRealScalarProductOperation*) (this->ThreadParameters[i].Operation))->SetScalarProducts(new double [operation->GetNbrScalarProduct()]);
      FirstComponent += Step;            
    }
  this->ThreadParameters[ReducedNbrThreads].Operation = operation->Clone();
  ((MultipleRealScalarProductOperation*) (this->ThreadParameters[ReducedNbrThreads].Operation))->SetIndicesRange(FirstComponent, 
														   operation->GetLeftVector()->GetVectorDimension() - 
														   FirstComponent);  
  ((MultipleRealScalarProductOperation*) (this->ThreadParameters[ReducedNbrThreads].Operation))->SetScalarProducts(new double [operation->GetNbrScalarProduct()]);
  this->SendJobs();
  for (int i = 1; i < this->NbrThreads; ++i)
    {
      for (int j = 0; j < operation->GetNbrScalarProduct(); ++j)
	operation->GetScalarProducts()[j] += ((MultipleRealScalarProductOperation*) (this->ThreadParameters[i].Operation))->GetScalarProducts()[j];
      delete[] ((MultipleRealScalarProductOperation*) (this->ThreadParameters[i].Operation))->GetScalarProducts();
      delete this->ThreadParameters[i].Operation;
    }
  return true;
}  

// execute an architecture-dependent multiple complex scalar product operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (MultipleComplexScalarProductOperation* operation)
{
  int Step = operation->GetLeftVector()->GetVectorDimension() / this->NbrThreads;
  int FirstComponent = 0;
  int ReducedNbrThreads = this->NbrThreads - 1;
  this->ThreadParameters[0].Operation = operation;
  ((MultipleComplexScalarProductOperation*) (this->ThreadParameters[0].Operation))->SetIndicesRange(FirstComponent, Step);
  FirstComponent += Step;
  for (int i = 1; i < ReducedNbrThreads; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((MultipleComplexScalarProductOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      ((MultipleComplexScalarProductOperation*) (this->ThreadParameters[i].Operation))->SetScalarProducts(new Complex [operation->GetNbrScalarProduct()]);
      FirstComponent += Step;            
    }
  this->ThreadParameters[ReducedNbrThreads].Operation = operation->Clone();
  ((MultipleComplexScalarProductOperation*) (this->ThreadParameters[ReducedNbrThreads].Operation))->SetIndicesRange(FirstComponent, 
														   operation->GetLeftVector()->GetVectorDimension() - 
														   FirstComponent);  
  ((MultipleComplexScalarProductOperation*) (this->ThreadParameters[ReducedNbrThreads].Operation))->SetScalarProducts(new Complex [operation->GetNbrScalarProduct()]);
  this->SendJobs();
  for (int i = 1; i < this->NbrThreads; ++i)
    {
      for (int j = 0; j < operation->GetNbrScalarProduct(); ++j)
	operation->GetScalarProducts()[j] += ((MultipleComplexScalarProductOperation*) (this->ThreadParameters[i].Operation))->GetScalarProducts()[j];
      delete[] ((MultipleComplexScalarProductOperation*) (this->ThreadParameters[i].Operation))->GetScalarProducts();
      delete this->ThreadParameters[i].Operation;
    }
  return true;
}  

// execute an architecture-dependent matrix matrix multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (MatrixMatrixMultiplyOperation* operation)
{
  int Step = operation->GetDestinationMatrix()->GetNbrRow() / this->NbrThreads;
  int FirstComponent = 0;
  int ReducedNbrThreads = this->NbrThreads - 1;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((MatrixMatrixMultiplyOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  this->ThreadParameters[ReducedNbrThreads].Operation = operation->Clone();
  ((MatrixMatrixMultiplyOperation*) (this->ThreadParameters[ReducedNbrThreads].Operation))->SetIndicesRange(FirstComponent, 
													      operation->GetDestinationMatrix()->GetNbrRow() - 
													      FirstComponent);  
  this->SendJobs();
  for (int i = 1; i < this->NbrThreads; ++i)
    {
      delete this->ThreadParameters[i].Operation;
    }
  return true;
}
    
// execute an architecture-dependent abstract hamiltonian precalculation operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SMPArchitecture::ExecuteOperation (AbstractPrecalculationOperation* operation)
{
  int Step = operation->GetHilbertSpaceDimension() / this->NbrThreads;
  int FirstComponent = 0;
  int ReducedNbrThreads = this->NbrThreads - 1;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->ThreadParameters[i].Operation = operation->Clone();
      ((AbstractPrecalculationOperation*) (this->ThreadParameters[i].Operation))->SetIndicesRange(FirstComponent, Step);
      FirstComponent += Step;
    }
  this->ThreadParameters[ReducedNbrThreads].Operation = operation->Clone();
  ((AbstractPrecalculationOperation*) (this->ThreadParameters[ReducedNbrThreads].Operation))->
    SetIndicesRange(FirstComponent, operation->GetHilbertSpaceDimension() - FirstComponent);  
  this->SendJobs();
  return true;
}
    
// send jobs to threads
//

void SMPArchitecture::SendJobs ()
{
  int Flag = 0;
  void* ret;
#ifdef __SMP__
  pthread_mutex_t mut = PTHREAD_MUTEX_INITIALIZER;
#endif
  int ReducedNbrThreads =  this->NbrThreads - 1;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->ThreadParameters[i].ThreadID = i;
      this->ThreadParameters[i].Flag = &Flag;
#ifdef __SMP__
      this->ThreadParameters[i].mut = &mut;
#endif
    }
  this->ThreadParameters[ReducedNbrThreads].ThreadID = ReducedNbrThreads;
  this->ThreadParameters[ReducedNbrThreads].Flag = &Flag;
#ifdef __SMP__
  this->ThreadParameters[ReducedNbrThreads].mut = &mut;
#endif

#ifdef __SMP__
  pthread_t* Threads2 = new pthread_t [this->NbrThreads];
  for (int i = 0; i < this->NbrThreads; ++i)
    {
      if (pthread_create (&(Threads2[i]), 0, ThreadExecuteOperation, (void*) &(this->ThreadParameters[i])) )
	{
	  cout << "error, cannot create thread" << endl;
	  exit(1);
	}
    }
  for (int i = 0; i < this->NbrThreads; ++i)
    {
      (void) pthread_join (Threads2[i], &ret);
    }
#endif

/*  for (int i = 0; i < this->NbrThreads; ++i)
    {
//      pthread_mutex_lock(this->ThreadParameters[i].LocalMutex);
      this->ThreadParameters[i].ThreadState = SMPArchitecture::Execute;
//      pthread_mutex_unlock(this->ThreadParameters[i].LocalMutex);
      pthread_cond_signal(this->ThreadParameters[i].LocalCondition);       
    }
  for (int i = 0; i < this->NbrThreads; ++i)
    {
      while (this->ThreadParameters[i].ThreadState != SMPArchitecture::Accomplished)
	pthread_cond_wait (this->ThreadParameters[i].LocalCondition, this->ThreadParameters[i].LocalMutex);      
//      pthread_mutex_lock(this->ThreadParameters[i].LocalMutex);
      this->ThreadParameters[i].ThreadState = SMPArchitecture::Wait;
//      cout << "stop waiting thread " << i << endl;
//      pthread_mutex_unlock(this->ThreadParameters[i].LocalMutex);
    }*/
}
  
// thread for multiplicationmain function
//
// param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
// return value = unused pointer (null)

void* ThreadExecuteOperation(void* param)
{
#ifdef __SMP__
  ThreadMainParameter* LocalThreadParamater = (ThreadMainParameter*) param;
  LocalThreadParamater->Operation->ApplyOperation();
  pthread_mutex_lock(LocalThreadParamater->mut);
  (*(LocalThreadParamater->Flag)) = LocalThreadParamater->ThreadID;
  pthread_mutex_unlock(LocalThreadParamater->mut);
#endif
  return 0;
}

// main function for thread
//
// param = pointer to additional parameters, has to be cast into ThreadMainParameter pointer
// return value = unused pointer (null)

void* ThreadMain(void* param)
{
  ThreadMainParameter* LocalThreadParamater = (ThreadMainParameter*) param;
#ifdef __SMP__
  pthread_mutex_t LaunchJobMutex = PTHREAD_MUTEX_INITIALIZER;
  pthread_cond_t LaunchJobCondition = PTHREAD_COND_INITIALIZER;
  LocalThreadParamater->LocalMutex = &LaunchJobMutex;
  LocalThreadParamater->LocalCondition = &LaunchJobCondition;
  while (LocalThreadParamater->ThreadState != SMPArchitecture::Exit)
    {
      if (LocalThreadParamater->ThreadState == SMPArchitecture::Execute)
	{
//	  cout << "thread " << LocalThreadParamater->ThreadID << ": receiving new job " << LocalThreadParamater->Operation->GetOperationType() << endl;
	  LocalThreadParamater->Operation->ApplyOperation();
	  LocalThreadParamater->ThreadState = SMPArchitecture::Accomplished;
//	  cout << "thread " << LocalThreadParamater->ThreadID << ": job accomplished" << endl;	    
	}
      pthread_cond_signal(LocalThreadParamater->LocalCondition);
      pthread_mutex_lock(LocalThreadParamater->LocalMutex);
      while ((LocalThreadParamater->ThreadState == SMPArchitecture::Wait) || (LocalThreadParamater->ThreadState == SMPArchitecture::Accomplished))
	pthread_cond_wait (LocalThreadParamater->LocalCondition, LocalThreadParamater->LocalMutex);      
      pthread_mutex_unlock(LocalThreadParamater->LocalMutex);
    }  
#endif
//  cout << "exit thread " << LocalThreadParamater->ThreadID << endl;
  LocalThreadParamater->ThreadState = SMPArchitecture::Dead;
  return 0;
}
