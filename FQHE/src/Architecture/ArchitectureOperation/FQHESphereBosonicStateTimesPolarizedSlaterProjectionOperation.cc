////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of operator Kostka numbers computation operation         //
//                                                                            //
//                        last modification : 09/05/2010                      //
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
#include "Architecture/ArchitectureOperation/FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation.h"

#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/BosonOnSphereTwoLandauLevels.h"
#include "HilbertSpace/FermionOnSphereThreeLandauLevels.h"
#include "HilbertSpace/FermionOnSphereFourLandauLevels.h"



#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;

FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation(ParticleOnSphere * initialSpace, FermionOnSphere * fermionSpace, FermionOnSphereWithSpin * finalSpace, RealVector* bosonicVector, RealVector* outputVector,bool twoLandauLevels,int nbrStage)
{
  this->FirstComponent = 0;	
  this->NbrComponent = initialSpace->GetHilbertSpaceDimension();
  this->FermionSpace = fermionSpace;
  this->InitialSpace = initialSpace;
  this->FinalSpace = finalSpace;
  this->OperationType = AbstractArchitectureOperation::FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation;
  this->OutputVector =  outputVector;
  this->BosonicVector = bosonicVector;
  this->NbrStage = nbrStage;
  this->TwoLandauLevels = twoLandauLevels;
}

// copy constructor 
//
// operation = reference on operation to copy
FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation(const FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation & operation)	
{
  this->FirstComponent = operation.FirstComponent;	
  this->NbrComponent 	= operation.NbrComponent;
  this->FermionSpace = operation.FermionSpace;		
  this->InitialSpace = operation.InitialSpace;	
  this->FinalSpace = (FermionOnSphereWithSpin *) operation.FinalSpace->Clone();
  this->OperationType = AbstractArchitectureOperation::FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation;
  this->OutputVector =  operation.OutputVector;
  this->BosonicVector = operation.BosonicVector;
  this->TwoLandauLevels = operation.TwoLandauLevels;
}

// destructor
//

FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::~FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation()
{
}

// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}

// set destination vector 
// 
// vector where the result has to be stored

void FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::SetOutputVector (RealVector* outputVector)
{
  this->OutputVector = outputVector;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation*FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::Clone()
{
  return new FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation (*this);
}


// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::RawApplyOperation()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);
  
  if (this->TwoLandauLevels)
    ((BosonOnSphereTwoLandauLevels * )this->InitialSpace)->BosonicStateTimePolarizedSlaters(*this->BosonicVector, *this->OutputVector,this->FermionSpace , this->FinalSpace, this->FirstComponent ,this->NbrComponent);
  else
    ((BosonOnSphereShort *)this->InitialSpace)->BosonicStateTimePolarizedSlaters(*this->BosonicVector, *this->OutputVector,this->FermionSpace , this->FinalSpace, this->FirstComponent ,this->NbrComponent);
	
  timeval TotalEndingTime;
  gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout << this->FirstComponent << " " <<  this->NbrComponent << " : " << Dt << "s" << endl;
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = (int) this->NbrComponent / (this->NbrStage*architecture->GetNbrThreads());
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
	
  FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation** TmpOperations = new FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation* [architecture->GetNbrThreads()];
  for(int i = 0 ;i < architecture->GetNbrThreads() ;i++)
    {
      TmpOperations[i] = (FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      TmpOperations[i]->SetOutputVector((RealVector*)this->OutputVector->EmptyClone(true));
    }
  
    for(int p = 0 ; p <this->NbrStage-1 ; p++)
      {
	for (int i = 0; i < ReducedNbrThreads; ++i)
	  {	
	    TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
	    TmpFirstComponent += Step;
	  }
      TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->FirstComponent+(p+1)*Step*architecture->GetNbrThreads() - TmpFirstComponent);
      TmpFirstComponent = this->FirstComponent+(p+1)*Step*architecture->GetNbrThreads();
      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
	{
	  TmpOperations[i]->OutputVector->ClearVector();
	}
      architecture->SendJobs();
      cout << TmpFirstComponent << " /  " << this->NbrComponent << " (" << ((TmpFirstComponent * 100) / this->NbrComponent) << "%)                   \r";
      cout.flush();
      
      for (int i = 1; i < architecture->GetNbrThreads(); ++i)
	{
	  (*(this->OutputVector)) += (*(TmpOperations[i]->OutputVector));
	}
      }
    
    for (int i = 0; i < ReducedNbrThreads; ++i)
      {	
	TmpOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
	TmpFirstComponent += Step;
      }
    TmpOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->FirstComponent+this->NbrComponent - TmpFirstComponent);
    
    for (int i = 1; i < architecture->GetNbrThreads(); ++i)
      {
	TmpOperations[i]->OutputVector->ClearVector();
    }
    architecture->SendJobs();
    for (int i = 1; i < architecture->GetNbrThreads(); ++i)
      {
	(*(this->OutputVector)) += (*(TmpOperations[i]->OutputVector));
	delete TmpOperations[i]->OutputVector;
	delete TmpOperations[i];
      }
    delete TmpOperations[0];
    delete[] TmpOperations;
    return true;
}

// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__
  if (architecture->IsMasterNode())
    {
      this->OutputVector->ClearVector();
      RealVector* TmpVector = (RealVector*) this->OutputVector->EmptyClone(true);
      int Step = (int) this->NbrComponent / (architecture->GetNbrSlaveNodes());
      int TmpFirstComponent = this->FirstComponent;
      int TmpRange[2];
      int TmpNbrSlaves = architecture->GetNbrSlaveNodes();
      int TmpSlaveID = 0;
      for (int i = 0; (i < TmpNbrSlaves) && (Step >= 0); ++i)
	{
	  TmpRange[0] = TmpFirstComponent;
	   TmpRange[1] = Step;
	   TmpFirstComponent += Step;
	   if (TmpFirstComponent >= (this->FirstComponent +this->NbrComponent))
	     Step = this->FirstComponent + this->NbrComponent - TmpFirstComponent;
	}
      while ((Step > 0) && ((TmpSlaveID = architecture->WaitAnySlave())))
	{
	  architecture->SumVector(*(this->OutputVector));
	}
      while ((TmpNbrSlaves > 0) && ((TmpSlaveID = architecture->WaitAnySlave())))
	{
	  --TmpNbrSlaves;
	  architecture->SumVector(*(this->OutputVector));
	}
    }
  else
    {
      int TmpRange[2];
      int TmpNbrElement = 0;
      while ((architecture->ReceiveFromMaster(TmpRange, TmpNbrElement) == true) && (TmpRange[1] > 0))
	{
	  this->OutputVector->ClearVector();
	  this->SetIndicesRange(TmpRange[0], TmpRange[1]); 
	  if (architecture->GetLocalArchitecture()->GetArchitectureID() == AbstractArchitecture::SMP)
	    this->ArchitectureDependentApplyOperation((SMPArchitecture*) architecture->GetLocalArchitecture());
	  else
	    this->RawApplyOperation();
	  architecture->SendDone();
	  architecture->SumVector(*(this->OutputVector));
	}
    }
  return true;
#else
  return this->RawApplyOperation();
#endif
}
