////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of simple MPI Architecture                     //
//                                                                            //
//                        last modification : 17/05/2004                      //
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
#include "Architecture/SimpleMPIArchitecture.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include <sys/time.h>
#include <string.h>
#include <iostream>
#ifdef __MPI__
#include <mpi.h>
#endif


using std::cout;
using std::endl;


// constructor
//

SimpleMPIArchitecture::SimpleMPIArchitecture()
{
  this->PerformanceIndex = 1.0;
  this->ArchitectureID = AbstractArchitecture::SimpleMPI;
#ifdef __MPI__
  MPI::Init();
  this->NbrMPINodes = MPI::COMM_WORLD.Get_size();
  this->MPIRank = MPI::COMM_WORLD.Get_rank();
  this->ClusterPerformanceArray = new double [this->NbrMPINodes];

  if (this->MPIRank != 0)
    {
      this->MasterNodeFlag = false;
      MPI::COMM_WORLD.Send(&this->PerformanceIndex, 1, MPI::DOUBLE, 0, 1);
    }
  else
    {
      this->MasterNodeFlag = true;
      this->TotalPerformanceIndex = this->PerformanceIndex;
      this->ClusterPerformanceArray[0] = this->PerformanceIndex;
      for (int i = 1; i < this->NbrMPINodes; ++i)
	{
	  MPI::COMM_WORLD.Recv(&this->ClusterPerformanceArray[i], 1, MPI::DOUBLE, i, 1);	  
	  this->TotalPerformanceIndex += this->ClusterPerformanceArray[i];
	}
      for (int i = 0; i < this->NbrMPINodes; ++i)
	{
	  this->ClusterPerformanceArray[i] /= this->TotalPerformanceIndex;
	}
    }
  MPI::COMM_WORLD.Bcast(this->ClusterPerformanceArray, this->NbrMPINodes, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&this->TotalPerformanceIndex, 1, MPI::DOUBLE, 0);
#else
  this->MasterNodeFlag = true;
  this->NbrMPINodes = 1;
  this->MPIRank = 0;
  this->ClusterPerformanceArray = 0;
  this->TotalPerformanceIndex = this->PerformanceIndex;
#endif
  if (this->MasterNodeFlag == true)
    {
      cout << this->NbrMPINodes << " " << this->TotalPerformanceIndex << endl;
    }
  this->LocalArchitecture = new MonoProcessorArchitecture;
}
  
// destructor
//

SimpleMPIArchitecture::~SimpleMPIArchitecture()
{
#ifdef __MPI__
  MPI::Finalize();
#endif
  delete this->LocalArchitecture;
  if (this->ClusterPerformanceArray != 0)
    delete[] this->ClusterPerformanceArray;
}
  
// get typical range of indices on which the local architecture acts
//
// minIndex = reference on the minimum index on which the local architecture can act
// maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
//            architecture doesn't support this feature)

void SimpleMPIArchitecture::GetTypicalRange (long& minIndex, long& maxIndex)
{
  minIndex = this->MinimumIndex;
  maxIndex = this->MaximumIndex;
}
  
// set dimension of the Hilbert space on which the architecture has to work
// 
// dimension = dimension of the Hilbert space

void SimpleMPIArchitecture::SetDimension (long dimension)
{
  this->HilbertSpaceDimension = dimension;
  double Tmp = 0.0;
  for (int i = 0; i < this->MPIRank; ++i)
    Tmp += this->ClusterPerformanceArray[i];
  if (this->MPIRank == 0)
    {
      this->MinimumIndex = (long) 0;
    }
  else
    {
      cout << Tmp << endl;
      this->MinimumIndex = (long) (Tmp * ((double) dimension));
    }
  if (this->MPIRank == (this->NbrMPINodes - 1))
    {
      this->MaximumIndex = dimension - 1;
    }
  else
    {
      Tmp += this->ClusterPerformanceArray[this->MPIRank];      
      this->MaximumIndex = (long) (Tmp * ((double) dimension)) - (long) 1;
    }
//   this->MinimumIndex = (long) 0;
//   this->MaximumIndex = dimension - 1;
  cout << this->MPIRank << " " << this->MinimumIndex << " " << this->MaximumIndex << endl;
}

// request an operation to the slave nodes and wait till they are ready to get operation parameters
//
// operationType = operation ID
// return value = true if no error occured

bool SimpleMPIArchitecture::RequestOperation (int operationType)
{
#ifdef __MPI__
  if (this->MasterNodeFlag)
    {
      MPI::COMM_WORLD.Bcast(&operationType, 1, MPI::INT, 0);
      int NbrMPINodes = MPI::COMM_WORLD.Get_size();
      bool Flag = true;
      int Acknowledge = 0;
      for (int i = 1; i < NbrMPINodes; ++i)
	{
	  MPI::COMM_WORLD.Recv(&Acknowledge, 1, MPI::INT, i, 1);      
	  if ((Flag == true) && (Acknowledge == 0))
	    Flag = false;
	}
      return Flag;
    }
#endif
  return false;
}

// wait an operation request from the master node  (without sending acknowledge)
//
// operationType = reference on the integer where the operation ID will be stored
// return value = true until the free slave signal is sent or an error occurs

bool SimpleMPIArchitecture::WaitOperation (int& operationType)
{
#ifdef __MPI__
  if (this->MasterNodeFlag == false)
    {
      MPI::COMM_WORLD.Bcast(&operationType, 1, MPI::INT, 0);
      if (operationType == SimpleMPIArchitecture::FreeSlaveSignal)
	{
	  return false;
	}
      else
	{
	  return true;
	}
    }
#endif
  return false;
}


// send acknowledge to the master node 
//
// acknowledge = true to send a positive answer
// return value = true if no error occured

bool SimpleMPIArchitecture::SendAcknowledge (bool acknowledge)
{
#ifdef __MPI__
  if (!this->MasterNodeFlag)
    {
      int Acknowledge = 0;
      if (acknowledge == true)
	 Acknowledge = 1;
      MPI::COMM_WORLD.Send(&Acknowledge, 1, MPI::INT, 0, 1); 
      return true;
    }
#endif
  return false;
}

// broadcast an integer from master node to slave nodes
// 
// value = integer to broadcast
// return value = true if no error occured

bool SimpleMPIArchitecture::BroadcastToSlaves(int& value)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Bcast(&value, 1, MPI::INT, 0); 
  return true;
#else
  return false;
#endif  
}

// broadcast an integer array from master node to slave nodes
// 
// values = array of integesr to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::BroadcastToSlaves(int* values, int nbrValues)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Bcast(values, nbrValues, MPI::INT, 0); 
  return true;
#else
  return false;
#endif  
}

// broadcast a double from master node to slave nodes
// 
// value = integer to broadcast
// return value = true if no error occured

bool SimpleMPIArchitecture::BroadcastToSlaves(double& value)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Bcast(&value, 1, MPI::DOUBLE, 0); 
  return true;
#else
  return false;
#endif  
}

// broadcast a double array from master node to slave nodes
// 
// values = array of integesr to broadcast
// nbrValues = number of element in the array
// return value = true if no error occured

bool SimpleMPIArchitecture::BroadcastToSlaves(double* values, int nbrValues)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Bcast(values, nbrValues, MPI::DOUBLE, 0); 
  return true;
#else
  return false;
#endif  
}
  
// broadcast a vector on each slave node
//
// vector = pointer to the vector tobroadcast  (only usefull for the master node)
// return value = pointer to the broadcasted vector or null pointer if an error occured

Vector* SimpleMPIArchitecture::BroadcastVector(Vector* vector)
{
#ifdef __MPI__
  if ((this->MasterNodeFlag) && (vector != 0))
    {
      vector->BroadcastClone(MPI::COMM_WORLD, this->MPIRank);
      return vector;
    }
  else
    if (this->MasterNodeFlag == false)
      {
	Vector TmpVector;
	return TmpVector.BroadcastClone(MPI::COMM_WORLD, 0);
      }
#endif  
  return 0;
}

// broadcast a vector type and allocate a vector based on it on each slave node
//
// vector = pointer to the vector to be used as reference (only usefull for the master node)
// return value = pointer to the cloned vector or null pointer if an error occured

Vector* SimpleMPIArchitecture::BroadcastVectorType(Vector* vector)
{
#ifdef __MPI__
  if ((this->MasterNodeFlag) && (vector != 0))
    {
      vector->BroadcastEmptyClone(MPI::COMM_WORLD, this->MPIRank);
      return vector;
    }
  else
    if (this->MasterNodeFlag == false)
      {
	Vector TmpVector;
	return TmpVector.BroadcastEmptyClone(MPI::COMM_WORLD, 0);
      }
#endif  
  return 0;
}

// broadcast an array of vectors on each slave node
//
// nbrVectors = reference on the number of vectors to broadcast or get
// vector = pointer to the vector tobroadcast  (only usefull for the master node)
// return value =  pointer to the array of broadcasted vectors or null pointer if an error occured null pointer if an error occured

Vector** SimpleMPIArchitecture::BroadcastVectorArray(int& nbrVectors, Vector* vector)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Bcast(&nbrVectors, 1, MPI::INT, 0); 
  if ((this->MasterNodeFlag) && (vector != 0))
    {
      switch (vector->GetVectorType())
	{
	case Vector::RealDatas:
	  for (int i = 0; i < nbrVectors; ++i)
	    ((RealVector*) vector)[i].BroadcastClone(MPI::COMM_WORLD, this->MPIRank);
	  break;
	case Vector::ComplexDatas:	    
	  for (int i = 0; i < nbrVectors; ++i)
	    ((ComplexVector*) vector)[i].BroadcastClone(MPI::COMM_WORLD, this->MPIRank);
	  break;
	}
      return 0;
    }
  else
    if (this->MasterNodeFlag == false)
      {
	
	Vector** TmpVectorArray = new Vector*[nbrVectors];
	Vector TmpVector;
	for (int i = 0; i < nbrVectors; ++i)
	  TmpVectorArray[i] = TmpVector.BroadcastClone(MPI::COMM_WORLD, 0);
	return TmpVectorArray;
      }
#endif  
  return 0;
}

// broadcast vector type and allocate an array of vectors based on it on each slave node
//
// nbrVectors = reference on the number of vectors to broadcast or get
// vector = pointer to the vector to be used as reference (only usefull for the master node)
// return value =  pointer to the array of cloned vector or null pointer if an error occurednull pointer if an error occured

Vector** SimpleMPIArchitecture::BroadcastVectorTypeArray(int& nbrVectors, Vector* vector)
{
#ifdef __MPI__
  MPI::COMM_WORLD.Bcast(&nbrVectors, 1, MPI::INT, 0); 
  if ((this->MasterNodeFlag) && (vector != 0))
    {
      vector[0].BroadcastEmptyClone(MPI::COMM_WORLD, this->MPIRank);
      return 0;
    }
  else
    if (this->MasterNodeFlag == false)
      {
	Vector** TmpVectorArray = new Vector*[nbrVectors];
	Vector TmpVector;
	TmpVectorArray[0] = TmpVector.BroadcastEmptyClone(MPI::COMM_WORLD, 0);
	for (int i = 1; i < nbrVectors; ++i)
	  TmpVectorArray[i] = TmpVectorArray[0]->EmptyClone();
	return TmpVectorArray;
      }
#endif  
  return 0;
}

// get a temporary file name
//
// return value = string corresponding to a temporary file name

char* SimpleMPIArchitecture::GetTemporaryFileName()
{
  timeval Time;
  gettimeofday (&Time, 0);
  char* TmpString = new char [64];
  sprintf (TmpString, "diagam%d%d_node%d.tmp",(int)  Time.tv_sec, (int)  Time.tv_usec, this->MPIRank);
  return TmpString;
}
  



/*
// execute an architecture-dependent vector hamiltonian multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (VectorHamiltonianMultiplyOperation* operation)
{
#ifdef __MPI__
  operation->GetDestinationVector()->ClearVector();
  operation->SetIndicesRange(this->MinimumIndex, this->MaximumIndex - this->MinimumIndex + 1);
  operation->ApplyOperation();
  operation->GetDestinationVector()->SumVector(MPI::COMM_WORLD, 0);
  operation->GetDestinationVector()->BroadcastVector(MPI::COMM_WORLD, 0);
  return true;
#else
  return false;
#endif
}
  

// execute an architecture-dependent multiple vector hamiltonian multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (MultipleVectorHamiltonianMultiplyOperation* operation)
{
#ifdef __MPI__
  bool RealFlag = true;
  int NbrVectors = operation->GetNbrVectors();
   if (operation->GetDestinationComplexVectors() == 0)
    {
      RealFlag = true;
      for (int i = 0; i < NbrVectors; ++i)
	operation->GetDestinationRealVectors()[i].ClearVector();
    }
   else
     {
      for (int i = 0; i < NbrVectors; ++i)
	operation->GetDestinationComplexVectors()[i].ClearVector();
     }
   operation->SetIndicesRange(this->MinimumIndex, this->MaximumIndex - this->MinimumIndex + 1);
   operation->ApplyOperation();
   if (RealFlag == true)
     {
      for (int i = 0; i < NbrVectors; ++i)
	{
	  operation->GetDestinationRealVectors()[i].SumVector(MPI::COMM_WORLD, 0);
	  operation->GetDestinationRealVectors()[i].BroadcastVector(MPI::COMM_WORLD, 0);
	}
     }
   else
     {
       for (int i = 0; i < NbrVectors; ++i)
	 {
	   operation->GetDestinationComplexVectors()[i].SumVector(MPI::COMM_WORLD, 0);
	   operation->GetDestinationComplexVectors()[i].BroadcastVector(MPI::COMM_WORLD, 0);
	 }
     }
   return true;
#else
  return false;
#endif
}
  

// execute an architecture-dependent vector abstract scalar sum operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully
  
bool SimpleMPIArchitecture::ExecuteOperation (AbstractScalarSumOperation* operation)
{
#ifdef __MPI__
  operation->SetIndicesRange(this->MinimumIndex, this->MaximumIndex - this->MinimumIndex + 1);
  operation->ApplyOperation();
  Complex TmpResult;
  MPI::COMM_WORLD.Reduce(&(operation->GetScalar().Re), &(TmpResult.Re), 1, MPI::DOUBLE, MPI::SUM, 0);
  MPI::COMM_WORLD.Reduce(&(operation->GetScalar().Im), &(TmpResult.Im), 1, MPI::DOUBLE, MPI::SUM, 0);
  operation->GetScalar() = TmpResult;
  MPI::COMM_WORLD.Bcast(&(operation->GetScalar().Re), 1, MPI::DOUBLE, 0);  
  MPI::COMM_WORLD.Bcast(&(operation->GetScalar().Im), 1, MPI::DOUBLE, 0);  
  return true;
#else
  return false;
#endif
}

// execute an architecture-dependent add real linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (AddRealLinearCombinationOperation* operation)
{
#ifdef __MPI__
  operation->SetIndicesRange(this->MinimumIndex, this->MaximumIndex - this->MinimumIndex + 1);
  operation->ApplyOperation();
  for (int i = 0; i < this->NbrMPINodes; ++i)
    if (i == this->MPIRank)
      {
	operation->GetDestinationVector()->BroadcastPartialVector(MPI::COMM_WORLD, i, this->MinimumIndex, this->MaximumIndex - this->MinimumIndex + 1);
      }
    else
      {
	operation->GetDestinationVector()->BroadcastPartialVector(MPI::COMM_WORLD, i);
      }
  return true;
#else
  return false;
#endif
}  

// execute an architecture-dependent add complex linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (AddComplexLinearCombinationOperation* operation)
{
#ifdef __MPI__
  operation->SetIndicesRange(this->MinimumIndex, this->MaximumIndex - this->MinimumIndex + 1);
  operation->ApplyOperation();
  for (int i = 0; i < this->NbrMPINodes; ++i)
    if (i == this->MPIRank)
      {
	operation->GetDestinationVector()->BroadcastPartialVector(MPI::COMM_WORLD, i, this->MinimumIndex, this->MaximumIndex - this->MinimumIndex + 1);
      }
    else
      {
	operation->GetDestinationVector()->BroadcastPartialVector(MPI::COMM_WORLD, i);
      }
  return true;
#else
  return false;
#endif
}  

// execute an architecture-dependent multiple real scalar product operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (MultipleRealScalarProductOperation* operation)
{
#ifdef __MPI__
  operation->SetIndicesRange(this->MinimumIndex, this->MaximumIndex - this->MinimumIndex + 1);
  operation->ApplyOperation();
  double* TmpScalarProducts = 0;
  if (this->MPIRank == 0)
    TmpScalarProducts = new double [operation->GetNbrScalarProduct()];
  MPI::COMM_WORLD.Reduce(operation->GetScalarProducts(), TmpScalarProducts, operation->GetNbrScalarProduct(), MPI::DOUBLE, MPI::SUM, 0);
  if (this->MPIRank == 0)
    {
      for (int i = 0; i < operation->GetNbrScalarProduct(); ++i)
	operation->GetScalarProducts()[i] = TmpScalarProducts[i];
      delete[] TmpScalarProducts;
    }
  MPI::COMM_WORLD.Bcast(operation->GetScalarProducts(), operation->GetNbrScalarProduct(), MPI::DOUBLE, 0);  
  return true;
#else
  return false;
#endif
}  

// execute an architecture-dependent multiple complex scalar product operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (MultipleComplexScalarProductOperation* operation)
{
#ifdef __MPI__
  operation->SetIndicesRange(this->MinimumIndex, this->MaximumIndex - this->MinimumIndex + 1);
  operation->ApplyOperation();
  Complex* TmpScalarProducts = 0;
  if (this->MPIRank == 0)
    TmpScalarProducts = new Complex [operation->GetNbrScalarProduct()];
  MPI::COMM_WORLD.Reduce(operation->GetScalarProducts(), TmpScalarProducts, operation->GetNbrScalarProduct(), MPI::COMPLEX, MPI::SUM, 0);
  if (this->MPIRank == 0)
    {
      for (int i = 0; i < operation->GetNbrScalarProduct(); ++i)
	operation->GetScalarProducts()[i] = TmpScalarProducts[i];
      delete[] TmpScalarProducts;
    }
  MPI::COMM_WORLD.Bcast(operation->GetScalarProducts(), operation->GetNbrScalarProduct(), MPI::COMPLEX, 0);  
  return true;
#else
  return false;
#endif
}  

// execute an architecture-dependent matrix matrix multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (MatrixMatrixMultiplyOperation* operation)
{
  return operation->ApplyOperation(this->LocalArchitecture);
}
 
// execute an architecture-dependent abstract hamiltonian precalculation operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (AbstractPrecalculationOperation* operation)
{
#ifdef __MPI__
  operation->SetIndicesRange(this->MinimumIndex, this->MaximumIndex - this->MinimumIndex + 1);
  operation->ApplyOperation();
  return true;
#else
  return false;
#endif
}
*/
