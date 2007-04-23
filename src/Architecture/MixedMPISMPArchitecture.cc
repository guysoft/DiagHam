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


#include "config.h"
#include "Architecture/MixedMPISMPArchitecture.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Architecture/SMPArchitecture.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include <sys/time.h>
#include <string.h>
#include <iostream>
#include <unistd.h>
#ifdef __MPI__
#include <mpi.h>
#endif


using std::cout;
using std::endl;


// constructor
//
// clusterFileName = name of the file that describes the cluster, if none assume one cpu per MPI node. The file should be at least accessible by the master mode

MixedMPISMPArchitecture::MixedMPISMPArchitecture(char* clusterFileName)
{
  this->PerformanceIndex = 1.0;
  this->ArchitectureID = AbstractArchitecture::MixedMPISMP;
#ifdef __MPI__
   this->NbrCPUPerNode = new int [this->NbrMPINodes];

  char* TmpLocalHostname = new char [512];
  gethostname(TmpLocalHostname, 511);

  if (this->MPIRank != 0)
    {
      this->NodeHostnames = 0;
      int HostnameStringSize = strlen(TmpLocalHostname) + 1;
      MPI::COMM_WORLD.Send(&HostnameStringSize, 1, MPI::INT, 0, 1);      
      MPI::COMM_WORLD.Send(TmpLocalHostname, HostnameStringSize, MPI::CHAR, 0, 1);      
    }
  else
    {
      this->NodeHostnames = new char*[this->NbrMPINodes];
      int HostnameStringSize = strlen(TmpLocalHostname) + 1;
      this->NodeHostnames[0] = new char[HostnameStringSize];
      strncpy (this->NodeHostnames[0], TmpLocalHostname, HostnameStringSize);
      for (int i = 1; i < this->NbrMPINodes; ++i)
	{
	  MPI::COMM_WORLD.Recv(&HostnameStringSize, 1, MPI::INT, i, 1);
	  this->NodeHostnames[i] = new char[HostnameStringSize];
	  MPI::COMM_WORLD.Recv(this->NodeHostnames[i], HostnameStringSize, MPI::CHAR, i, 1);	  
	}
    }
  delete[] TmpLocalHostname;

  if (this->MPIRank == 0)
    {
      if (clusterFileName != 0)
	{
	  MultiColumnASCIIFile ClusterFile;
	  bool ErrorFlag = false;
	  if (ClusterFile.Parse(clusterFileName) == true)
	    {
	      if (ClusterFile.GetNbrColumns() < 2)
		ErrorFlag = true;
	      else
		{
		  int* TmpNbrCPUNode = ClusterFile.GetAsIntegerArray(1);
		  if (TmpNbrCPUNode != 0)
		    {
		      for (int i = 0; i < this->NbrMPINodes; ++i)
			{
			  this->NbrCPUPerNode[i] = 1;
			  for (int j = 0; j < ClusterFile.GetNbrLines(); ++j)
			    {
			      cout << this->NodeHostnames[i] << " " << ClusterFile(0,j) << endl;
			      if ((strcmp(this->NodeHostnames[i], ClusterFile(0,j)) == 0) || 
				  ((strncmp(this->NodeHostnames[i], ClusterFile(0,j), strlen(ClusterFile(0,j))) == 0) && (this->NodeHostnames[i][strlen(ClusterFile(0,j))] = '.')))
				{
				  this->NbrCPUPerNode[i] = TmpNbrCPUNode[j];
				  j = ClusterFile.GetNbrLines();
				}
			    }
			  delete[] TmpNbrCPUNode;
			}
		    }
		  else
		    ErrorFlag = true;
		}
	    }
	  if (ErrorFlag == true)
	    {
	      cout << "an error occured while opening " << clusterFileName << endl;
	      ClusterFile.DumpErrors(cout);
	      cout << "switching to one cpu per node mode" << endl;
	      clusterFileName = 0;
	    }
	}
      if (clusterFileName == 0)
	for (int i = 0; i < this->NbrMPINodes; ++i)
	  this->NbrCPUPerNode[i] = 1;
    }
  MPI::COMM_WORLD.Bcast(this->NbrCPUPerNode, this->NbrMPINodes, MPI::INT, 0);
  this->PerformanceIndex *= (double) this->NbrCPUPerNode[this->MPIRank];

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
  this->MasterNodeFlag = true;
  this->NbrMPINodes = 1;
  this->MPIRank = 0;
  this->ClusterPerformanceArray = 0;
  this->TotalPerformanceIndex = this->PerformanceIndex;
  this->NbrCPUPerNode = new int[1];
  this->NbrCPUPerNode[0] = 1;
#endif
  if (this->MasterNodeFlag == true)
    {
      cout << this->NbrMPINodes << " " << this->TotalPerformanceIndex << endl;
    }
  if (this->NbrCPUPerNode[this->MPIRank] > 1)
    {
      delete this->LocalArchitecture;
      this->LocalArchitecture = new SMPArchitecture(this->NbrCPUPerNode[this->MPIRank]);
    }
}
  
// destructor
//

MixedMPISMPArchitecture::~MixedMPISMPArchitecture()
{
  delete[] this->NbrCPUPerNode;
  if (this->NodeHostnames != 0)
    {
      for (int i = 0; i < this->NbrMPINodes; ++i)
	delete[] this->NodeHostnames[i];
      delete[] this->NodeHostnames;
    }
}
  
// get typical range of indices on which the local architecture acts
//
// minIndex = reference on the minimum index on which the local architecture can act
// maxIndex = reference on the maximum index on which the local architecture can act (= minIndex is the 
//            architecture doesn't support this feature)

void MixedMPISMPArchitecture::GetTypicalRange (long& minIndex, long& maxIndex)
{
  minIndex = this->MinimumIndex;
  maxIndex = this->MaximumIndex;
}
  
// set dimension of the Hilbert space on which the architecture has to work
// 
// dimension = dimension of the Hilbert space

void MixedMPISMPArchitecture::SetDimension (long dimension)
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

bool MixedMPISMPArchitecture::RequestOperation (int operationType)
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

bool MixedMPISMPArchitecture::WaitOperation (int& operationType)
{
#ifdef __MPI__
  if (this->MasterNodeFlag == false)
    {
      MPI::COMM_WORLD.Bcast(&operationType, 1, MPI::INT, 0);
      if (operationType == MixedMPISMPArchitecture::FreeSlaveSignal)
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

bool MixedMPISMPArchitecture::SendAcknowledge (bool acknowledge)
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

