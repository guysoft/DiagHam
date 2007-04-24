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

};


#endif
