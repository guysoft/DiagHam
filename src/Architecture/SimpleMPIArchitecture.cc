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
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MatrixMatrixMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"

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
#ifdef __MPI__
  MPI::Init();
  this->NbrMPINodes = MPI::COMM_WORLD.Get_size();
  this->MPIRank = MPI::COMM_WORLD.Get_rank();
  this->ClusterPerformanceArray = new double [this->NbrMPINodes];

  if (this->MPIRank != 0)
    {
      MPI::COMM_WORLD.Send(&this->PerformanceIndex, 1, MPI::DOUBLE, 0, 1);
    }
  else
    {
      this->TotalPerformanceIndex = this->PerformanceIndex;
      this->ClusterPerformanceArray[0] = this->PerformanceIndex;
      for (int i = 1; i < this->NbrMPINodes; ++i)
	{
	  MPI::COMM_WORLD.Recv(&this->ClusterPerformanceArray[i], 1, MPI::DOUBLE, i, 1);	  
	  this->TotalPerformanceIndex += this->PerformanceIndex;
	}
      for (int i = 0; i < this->NbrMPINodes; ++i)
	{
	  this->ClusterPerformanceArray[i] /= this->TotalPerformanceIndex;
	}
    }
  MPI::COMM_WORLD.Bcast(this->ClusterPerformanceArray, this->NbrMPINodes, MPI::DOUBLE, 0);
  MPI::COMM_WORLD.Bcast(&this->TotalPerformanceIndex, 1, MPI::DOUBLE, 0);
#else
  this->NbrMPINodes = 1;
  this->MPIRank = 0;
  this->ClusterPerformanceArray = 0;
  this->TotalPerformanceIndex = this->PerformanceIndex;
#endif
  if (this->MPIRank == 0)
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
  cout << this->MPIRank << " " << this->MinimumIndex << " " << this->MaximumIndex << endl;
}

// execute an architecture-dependent vector hamiltonian multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (VectorHamiltonianMultiplyOperation* operation)
{
  return this->LocalArchitecture->ExecuteOperation(operation);
}
  

// execute an architecture-dependent add real linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (AddRealLinearCombinationOperation* operation)
{
  return this->LocalArchitecture->ExecuteOperation(operation);
}  

// execute an architecture-dependent add complex linear combination operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (AddComplexLinearCombinationOperation* operation)
{
  return this->LocalArchitecture->ExecuteOperation(operation);
}  

// execute an architecture-dependent multiple real scalar product operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (MultipleRealScalarProductOperation* operation)
{
  return this->LocalArchitecture->ExecuteOperation(operation);
}  

// execute an architecture-dependent multiple complex scalar product operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (MultipleComplexScalarProductOperation* operation)
{
  return this->LocalArchitecture->ExecuteOperation(operation);
}  

// execute an architecture-dependent matrix matrix multiplication operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (MatrixMatrixMultiplyOperation* operation)
{
  return this->LocalArchitecture->ExecuteOperation(operation);
}
 
// execute an architecture-dependent abstract hamiltonian precalculation operation
//
// operation = pointer to the operation to execute
// return value = true if operation has been completed successfully

bool SimpleMPIArchitecture::ExecuteOperation (AbstractPrecalculationOperation* operation)
{
  return this->LocalArchitecture->ExecuteOperation(operation);
}
