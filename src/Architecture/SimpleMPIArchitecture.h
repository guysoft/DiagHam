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


#ifndef SIMPLEMPIARCHITECTURE_H
#define SIMPLEMPIARCHITECTURE_H


#include "config.h"
#include "Architecture/AbstractArchitecture.h"
#include "Vector/Vector.h"

#ifdef __MPI__
#include <mpi.h>
#endif


class AbstractArchitectureOperation;


class SimpleMPIArchitecture : public AbstractArchitecture
{

 private:

  // total number of MPI nodes
  int NbrMPINodes;
  // rank of the current MPI node
  int MPIRank;
  // flag to inidcate if local node is the master mode
  bool MasterNodeFlag;

  // pointer to the architecture used for local operation
  AbstractArchitecture* LocalArchitecture;

  // current node performance index
  double PerformanceIndex;
  // cluster total performance index
  double TotalPerformanceIndex;
  // array containing performance index of each node
  double* ClusterPerformanceArray;

  // minimum index on which the current MPI node can act
  long MinimumIndex;
  // maximum index on which the current MPI node can act
  long MaximumIndex;

 public:
  
  enum SimpleMPISignals{
    FreeSlaveSignal = 0x7fffffff,
    SynchronizeSignal = 0x6fffffff
  };

  // constructor
  //
  SimpleMPIArchitecture();
  
  // destructor
  //
  ~SimpleMPIArchitecture();
  
  // get typical range of indices on which the local architecture acts
  //
  // minIndex = reference on the minimum index on which the local architecture can act
  // maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
  //            architecture doesn't support this feature)
  void GetTypicalRange (long& minIndex, long& maxIndex);
  
  // set dimension of the Hilbert space on which the architecture has to work
  // 
  // dimension = dimension of the Hilbert space
  void SetDimension (long dimension);

  // indicate if the local node is the master node
  // 
  // return value = true if the local node is the master node
  bool IsMasterNode();
  
  // request an operation to the slave nodes 
  //
  // operationType = operation ID
  // return value = true if no error occured
  bool RequestOperation (int operationType);

  // wait an operation request from the master node (without sending acknowledge)
  //
  // operationType = reference on the integer where the operation ID will be stored
  // return value = true until the free slave signal is sent or an error occurs
  bool WaitOperation (int& operationType);

  // send acknowledge to the master node 
  //
  // acknowledge = true to send a positive answer
  // return value = true if no error occured
  bool SendAcknowledge (bool acknowledge = true);

  // broadcast an integer from master node to slave nodes
  // 
  // value = integer to broadcast
  // return value = true if no error occured
  bool BroadcastToSlaves(int& value);

  // broadcast an integer array from master node to slave nodes
  // 
  // values = array of integesr to broadcast
  // nbrValues = number of element in the array
  // return value = true if no error occured
  bool BroadcastToSlaves(int* values, int nbrValues);

  // broadcast a double from master node to slave nodes
  // 
  // value = integer to broadcast
  // return value = true if no error occured
  bool BroadcastToSlaves(double& value);

  // broadcast a double array from master node to slave nodes
  // 
  // values = array of integesr to broadcast
  // nbrValues = number of element in the array
  // return value = true if no error occured
  bool BroadcastToSlaves(double* values, int nbrValues);

  // broadcast a vector on each slave node
  //
  // vector = pointer to the vector tobroadcast  (only usefull for the master node)
  // return value = pointer to the broadcasted vector or null pointer if an error occured
  Vector* BroadcastVector(Vector* vector = 0);

  // broadcast a vector type and allocate a vector based on it on each slave node
  //
  // vector = pointer to the vector to be used as reference (only usefull for the master node)
  // return value = pointer to the cloned vector or null pointer if an error occured
  Vector* BroadcastVectorType(Vector* vector = 0);

  // broadcast an array of vectors on each slave node
  //
  // nbrVectors = reference on the number of vectors to broadcast or get
  // vector = pointer to the vector tobroadcast  (only usefull for the master node)
  // return value =  pointer to the array of broadcasted vectors or null pointer if an error occured null pointer if an error occured
  Vector** BroadcastVectorArray(int& nbrVectors, Vector* vector = 0);

  // broadcast vector type and allocate an array of vectors based on it on each slave node
  //
  // nbrVectors = reference on the number of vectors to broadcast or get
  // vector = pointer to the vector to be used as reference (only usefull for the master node)
  // return value =  pointer to the array of cloned vector or null pointer if an error occurednull pointer if an error occured
  Vector** BroadcastVectorTypeArray(int& nbrVectors, Vector* vector = 0);

  // add current vector to the one of the master nide
  // 
  // vector = reference on the vector to add (or the destination vector of the master node)
  // return value = reference on the vector
  Vector& SumVector(Vector& vector);

  // get a temporary file name
  //
  // return value = string corresponding to a temporary file name
  char* GetTemporaryFileName();

};

// indicate if the local node is the master node
// 
// return value = true if the local node is the master node

inline bool SimpleMPIArchitecture::IsMasterNode()
{
  return this->MasterNodeFlag;
}

// add current vector to the one of the master nide
// 
// vector = reference on the vector to add (or the destination vector of the master node)
// return value = reference on the vector

inline Vector& SimpleMPIArchitecture::SumVector(Vector& vector)
{
  return vector.SumVector(MPI::COMM_WORLD, 0);
}

#endif
