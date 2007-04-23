////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of mixed MPI - SMP Architecture                  //
//                                                                            //
//                        last modification : 23/04/2007                      //
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


#ifndef MIXEDMPISMPARCHITECTURE_H
#define MIXEDMPISMPARCHITECTURE_H


#include "config.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "Vector/Vector.h"

#ifdef __MPI__
#include <mpi.h>
#endif



class MixedMPISMPArchitecture : public SimpleMPIArchitecture
{

 protected:

  // number of cpu atteched to each MPI node
  int* NbrCPUPerNode;

  // array that conatins hostname of each MPI node (only relevant for the master node)
  char** NodeHostnames;

 public:
  
  // constructor
  //
  // clusterFileName = name of the file that describes the cluster, if none assume one cpu per MPI node. The file should be at least accessible by the master mode
  MixedMPISMPArchitecture(char* clusterFileName = 0);
  
  // destructor
  //
  ~MixedMPISMPArchitecture();
  
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

};


#endif
