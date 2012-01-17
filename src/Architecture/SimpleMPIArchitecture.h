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
#include <math.h>

#ifdef __MPI__
#include <mpi.h>
#endif

using std::cout;
using std::endl;


class AbstractArchitectureOperation;


class SimpleMPIArchitecture : public AbstractArchitecture
{
  
 friend class FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation;

 protected:

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
  // minimum index on which MPI nodes can act
  long* MinimumIndices;
  // maximum index on which MPI nodes can act
  long* MaximumIndices;

  // name of the optional log file to allow code profiling on MPI architecture
  char* LogFile;
  // flag to indicate if the log file option is activated
  bool VerboseModeFlag;

 public:
  
  enum SimpleMPISignals{
    FreeSlaveSignal = 0x7fffffff,
    SynchronizeSignal = 0x6fffffff
  };

  // constructor
  //
  // logFile = name of the optional log file to allow code profiling on MPI architecture
  SimpleMPIArchitecture(char* logFile = 0);
  
  // destructor
  //
  ~SimpleMPIArchitecture();
  
  // get typical range of indices on which the local architecture acts
  //
  // minIndex = reference on the minimum index on which the local architecture can act
  // maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
  //            architecture doesn't support this feature)
  virtual void GetTypicalRange (long& minIndex, long& maxIndex);
  
  // get the ID of the node that handles a given index
  //
  // index = index to check
  // return value = corresponding node ID
  virtual int GetNodeIDFromIndex(long index);

  // set dimension of the Hilbert space on which the architecture has to work
  // 
  // dimension = dimension of the Hilbert space
  virtual void SetDimension (long dimension);

  // indicate if the local node is the master node
  // 
  // return value = true if the local node is the master node
  virtual bool IsMasterNode();
  
  // indicate how many slave nodes are available
  //
  // return value = number of slave nodes
  virtual int GetNbrSlaveNodes();
  
  // indicate how many nodes are available
  //
  // return value = number of nodes
  virtual int GetNbrNodes();
  
  // get the rank of the current node
  //
  // return value = rank of current node
  virtual int GetNodeNbr();

  // get the architecture used on the local MPI node
  //
  // return value = pointer to the local architecture
  virtual AbstractArchitecture* GetLocalArchitecture();

  // request an operation to the slave nodes 
  //
  // operationType = operation ID
  // return value = true if no error occured
  virtual bool RequestOperation (int operationType);

  // wait an operation request from the master node (without sending acknowledge)
  //
  // operationType = reference on the integer where the operation ID will be stored
  // return value = true until the free slave signal is sent or an error occurs
  virtual bool WaitOperation (int& operationType);

  // send acknowledge to the master node 
  //
  // acknowledge = true to send a positive answer
  // return value = true if no error occured
  virtual bool SendAcknowledge (bool acknowledge = true);

  // wait for a slave to send the done signal
  //
  // return value = id of the slave that send the done signal
  virtual int WaitAnySlave ();

  // send done signal to the master node 
  //
  // return value = true if no error occured
  virtual bool SendDone ();

  // broadcast an integer from master node to slave nodes
  // 
  // value = integer to broadcast
  // return value = true if no error occured
  virtual bool BroadcastToSlaves(int& value);

  // broadcast an integer array from master node to slave nodes
  // 
  // values = array of integesr to broadcast
  // nbrValues = number of element in the array
  // return value = true if no error occured
  virtual bool BroadcastToSlaves(int* values, int nbrValues);

  // send an integer array from master node to a given slave node
  // 
  // slaveID = slave ID
  // values = array of integesr to broadcast
  // nbrValues = number of element in the array
  // return value = true if no error occured
  virtual bool SendToSlaves(int slaveID, int* values, int nbrValues);

  // receive an integer array from master node to the current slave node
  // 
  // values = array of integesr to broadcast
  // nbrValues = number of element in the array
  // return value = true if no error occured
  virtual bool ReceiveFromMaster(int* values, int& nbrValues);

  // send an integer array from the current slave node to master node
  // 
  // values = array of integesr to broadcast
  // nbrValues = number of element in the array
  // return value = true if no error occured
  virtual bool SendToMaster(int* values, int& nbrValues);

  // receive an integer array from master node to the current slave node
  // 
  // slaveID = slave ID
  // values = array of integesr to broadcast
  // nbrValues = number of element in the array
  // return value = true if no error occured
  virtual bool ReceiveFromSlave(int slaveID, int* values, int& nbrValues);

  // send a double array from the current slave node to master node
  // 
  // values = array of integesr to broadcast
  // nbrValues = number of element in the array
  // return value = true if no error occured
  virtual bool SendToMaster(double* values, int& nbrValues);

  // receive a double array from master node to the current slave node
  // 
  // slaveID = slave ID
  // values = array of integesr to broadcast
  // nbrValues = number of element in the array
  // return value = true if no error occured
  virtual bool ReceiveFromSlave(int slaveID, double* values, int& nbrValues);

  // send a double complex array from the current slave node to master node
  // 
  // values = array of integesr to broadcast
  // nbrValues = number of element in the array
  // return value = true if no error occured
  virtual bool SendToMaster(doublecomplex* values, int& nbrValues);

  // receive a double complex array from master node to the current slave node
  // 
  // slaveID = slave ID
  // values = array of integesr to broadcast
  // nbrValues = number of element in the array
  // return value = true if no error occured
  virtual bool ReceiveFromSlave(int slaveID, doublecomplex* values, int& nbrValues);

  // broadcast a double from master node to slave nodes
  // 
  // value = integer to broadcast
  // return value = true if no error occured
  virtual bool BroadcastToSlaves(double& value);

  // broadcast a double array from master node to slave nodes
  // 
  // values = array of integesr to broadcast
  // nbrValues = number of element in the array
  // return value = true if no error occured
  virtual bool BroadcastToSlaves(double* values, int nbrValues);

  // broadcast a vector on each slave node
  //
  // vector = pointer to the vector tobroadcast  (only usefull for the master node)
  // return value = pointer to the broadcasted vector or null pointer if an error occured
  virtual Vector* BroadcastVector(Vector* vector = 0);

  // broadcast a vector from a node to the others 
  //
  // nodeID = id of the mode that broadcasts its vector
  // vector = vector to broadcast or to the vector where the content will be stored
  void BroadcastVector(int nodeID, Vector& vector);

  // scatter a vector upon each slave node
  //
  // vector = pointer to the vector to scatter  (only usefull for the master node)
  // return value = pointer to the broadcasted vector or null pointer if an error occured
  virtual Vector* ScatterVector(Vector* vector = 0);

  // broadcast a vector type and allocate a vector based on it on each slave node
  //
  // vector = pointer to the vector to be used as reference (only usefull for the master node)
  // return value = pointer to the cloned vector or null pointer if an error occured
  virtual Vector* BroadcastVectorType(Vector* vector = 0);

  // broadcast an array of vectors on each slave node
  //
  // nbrVectors = reference on the number of vectors to broadcast or get
  // vector = pointer to the vector tobroadcast  (only usefull for the master node)
  // return value =  pointer to the array of broadcasted vectors or null pointer if an error occured null pointer if an error occured
  virtual Vector** BroadcastVectorArray(int& nbrVectors, Vector* vector = 0);

  // broadcast vector type and allocate an array of vectors based on it on each slave node
  //
  // nbrVectors = reference on the number of vectors to broadcast or get
  // vector = pointer to the vector to be used as reference (only usefull for the master node)
  // return value =  pointer to the array of cloned vector or null pointer if an error occurednull pointer if an error occured
  virtual Vector** BroadcastVectorTypeArray(int& nbrVectors, Vector* vector = 0);

  // scatter an array of vectors upon each slave node
  //
  // nbrVectors = reference on the number of vectors to broadcast or get
  // vector = pointer to the vector to be used as reference (only usefull for the master node)
  // return value = pointer to the broadcasted vector or null pointer if an error occured
  virtual Vector** ScatterVectorArray(int& nbrVectors, Vector* vector = 0);

  // add current vector to the one of the master nide
  // 
  // vector = reference on the vector to add (or the destination vector of the master node)
  // return value = reference on the vector
  virtual Vector& SumVector(Vector& vector);

  // reassemble current vector into the one of the master node
  // 
  // vector = reference on the vector to add (or the destination vector of the master node)
  // return value = reference on the vector
  virtual Vector& ReassembleVector(Vector& vector);

  // get a temporary file name
  //
  // return value = string corresponding to a temporary file name
  virtual char* GetTemporaryFileName();

  // indicate if the log file option is activated
  //
  // return value = true if the option is activated
  virtual bool VerboseMode();

  // add an entry to the log file 
  //
  // message = string corresponding to entry to add to the log file
  // masterFlag = true if only the master node should add the entry
  // return value = true if no error occured
  virtual bool AddToLog(const char * message, bool masterFlag = false);

  // dump the log file into a string
  //
  // header = optional header to add before the log file
  // footer = optional footer to add at the end of the log file
  // return value = string or 0 if an error occured or log is not available
  virtual char* DumpLog(const char* header = 0, const char* footer = 0);
  
  // write vector in a file 
  //
  // vector = vector to write
  // fileName = name of the file where the vector has to be stored
  // return value = true if no error occurs
  virtual bool WriteVector(RealVector& vector, const char* fileName);
  
  // read a vector from a file but only on master
  //
  // vector = vector to read
  // fileName = name of the file where the vector is read from
  // return value = true if no error occurs
  virtual bool ReadVector(RealVector& vector, const char* fileName);

};
// indicate if the local node is the master node
// 
// return value = true if the local node is the master node

inline bool SimpleMPIArchitecture::IsMasterNode()
{
  return this->MasterNodeFlag;
}

// get the architecture used on the local MPI node
//
// return value = pointer to the local architecture

inline AbstractArchitecture* SimpleMPIArchitecture::GetLocalArchitecture()
{
  return this->LocalArchitecture;
}

// indicate if the log file option is activated
//
// return value = true if the option is activated

inline bool SimpleMPIArchitecture::VerboseMode()
{
  return this->VerboseModeFlag;
}

// add current vector to the one of the master nide
// 
// vector = reference on the vector to add (or the destination vector of the master node)
// return value = reference on the vector

inline Vector& SimpleMPIArchitecture::SumVector(Vector& vector)
{
#ifdef __MPI__
  return vector.SumVector(MPI::COMM_WORLD, 0);
#else
  return vector;
#endif
}

// reassemble current vector into the one of the master node
// 
// vector = reference on the vector to add (or the destination vector of the master node)
// return value = reference on the vector

inline Vector& SimpleMPIArchitecture::ReassembleVector(Vector& vector)
{
#ifdef __MPI__
  return vector.ReassembleVector(MPI::COMM_WORLD, 0);
#else
  return vector;
#endif
}

// indicate how many nodes are available
//
// return value = number of nodes

inline int SimpleMPIArchitecture::GetNbrSlaveNodes()
{
  return (this->NbrMPINodes - 1);
}

// indicate how many nodes there are
//
// return value = number of nodes

inline int SimpleMPIArchitecture::GetNbrNodes()
{
  return this->NbrMPINodes;
}

// indicate number of current node
//
// return value = number of the current node

inline int SimpleMPIArchitecture::GetNodeNbr()
{
  return this->MPIRank;
}



#endif
