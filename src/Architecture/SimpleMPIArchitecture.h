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
  // return value = true if np error occured
  bool RequestOperation (int operationType);

};

// indicate if the local node is the master node
// 
// return value = true if the local node is the master node

inline bool SimpleMPIArchitecture::IsMasterNode()
{
  return this->MasterNodeFlag;
}

#endif
