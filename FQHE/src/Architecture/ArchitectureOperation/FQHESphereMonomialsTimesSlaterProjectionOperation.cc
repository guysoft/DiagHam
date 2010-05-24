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
#include "Architecture/ArchitectureOperation/FQHESphereMonomialsTimesSlaterProjectionOperation.h"

#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereTwoLandauLevels.h"



#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
class BosonOnSphereShort;
 
// constructor 
//
// space = pointer to the HilbertSpace
// fileName = file where the kostka Numbers will be stored

FQHESphereMonomialsTimesSlaterProjectionOperation::FQHESphereMonomialsTimesSlaterProjectionOperation(ParticleOnSphere* fermionSpace, ParticleOnSphere* lllSpace, ParticleOnSphere* finalSpace, 
										 RealVector* fermionVector, RealVector* lllVector, RealVector* outputVector, int resume, int nbrComponent, 
										 bool projection , int step, bool flag3LL, bool symmetry)
{
  this->FirstComponent = resume;
  if(nbrComponent==0)
    this->NbrComponent 	= lllSpace->GetHilbertSpaceDimension() - resume;
  else
    this->NbrComponent = nbrComponent;
  this->FermionSpace = fermionSpace;
  this->LLLSpace = lllSpace;
  this->FinalSpace = finalSpace;
  this->OperationType = AbstractArchitectureOperation::FQHESphereMonomialsTimesSlaterProjection;
  this->FermionVector = fermionVector;
  this->LLLVector = lllVector; 
  this->OutputVector= outputVector;
  this->Projection=projection;
  this->NbrStage=step;
  this->Flag3LL=flag3LL;
  this->Symmetry = symmetry;
  this->FinalStates = new unsigned long [finalSpace->GetHilbertSpaceDimension()];
  this->Weigth = new double [finalSpace->GetHilbertSpaceDimension()];
  if(lllSpace->GetParticleStatistic() == ParticleOnSphere::BosonicStatistic)
    this->BosonFlag=true;
  else
    this->BosonFlag=false;
}

// copy constructor 
//
// operation = reference on operation to copy

FQHESphereMonomialsTimesSlaterProjectionOperation::FQHESphereMonomialsTimesSlaterProjectionOperation(const FQHESphereMonomialsTimesSlaterProjectionOperation& operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->FermionSpace = operation.FermionSpace;
  this->LLLSpace = operation.LLLSpace;
  this->FinalSpace = (FermionOnSphere *) operation.FinalSpace->Clone();
  this->OperationType = AbstractArchitectureOperation::FQHESphereMonomialsTimesSlaterProjection;
  this->FermionVector = operation.FermionVector;
  this->LLLVector = operation.LLLVector;
  this->OutputVector = operation.OutputVector;
  this->Projection=operation.Projection;
  this->NbrStage=operation.NbrStage;
  this->Flag3LL=operation.Flag3LL;
  this->Symmetry=operation.Symmetry;
  this->FinalStates = new unsigned long [operation.FinalSpace->GetHilbertSpaceDimension()];
  this->Weigth = new double [operation.FinalSpace->GetHilbertSpaceDimension()];
  this->BosonFlag=operation.BosonFlag;
}
  
// destructor
//

FQHESphereMonomialsTimesSlaterProjectionOperation::~FQHESphereMonomialsTimesSlaterProjectionOperation()
{
  delete this->FinalStates;
  delete this->Weigth;
}
  
// set range of indices
// 
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHESphereMonomialsTimesSlaterProjectionOperation::SetIndicesRange (const int& firstComponent, const int& nbrComponent)
{
  this->FirstComponent = firstComponent;
  this->NbrComponent = nbrComponent;
}


// set destination vector 
// 
// vector where the result has to be stored

void FQHESphereMonomialsTimesSlaterProjectionOperation::SetOutputVector (RealVector* outputVector)
{
  this->OutputVector = outputVector;
}


// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHESphereMonomialsTimesSlaterProjectionOperation::Clone()
{
  return new FQHESphereMonomialsTimesSlaterProjectionOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereMonomialsTimesSlaterProjectionOperation::RawApplyOperation()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);
  if(this->BosonFlag==true)
    {
      if(this->Flag3LL==true)
	{
	  ((FermionOnSphereThreeLandauLevels *)this->FermionSpace)->BosonicStateTimeFermionicState(*(this->LLLVector), *(this->FermionVector), *(this->OutputVector),
												   this->FinalStates, this->Weigth, (BosonOnSphereShort*) this->LLLSpace, 
												   (FermionOnSphere*)this->FinalSpace,
												   this->FirstComponent, this->NbrComponent);
	}
      else
	{
	  if(this->Projection)
	    {
	      if(this->Symmetry)
		((FermionOnSphereTwoLandauLevels *)this->FermionSpace)->BosonicStateTimeFermionicStateSymmetric(*(this->LLLVector), *(this->FermionVector), *(this->OutputVector),this->FinalStates,this->Weigth,
														(BosonOnSphereShort*) this->LLLSpace, (FermionOnSphere*)this->FinalSpace,
														this->FirstComponent, this->NbrComponent);
	      else
		((FermionOnSphereTwoLandauLevels *)this->FermionSpace)->BosonicStateTimeFermionicState(*(this->LLLVector), *(this->FermionVector), *(this->OutputVector),this->FinalStates,this->Weigth, 
												       (BosonOnSphereShort*)this->LLLSpace,(FermionOnSphere*)this->FinalSpace,
												       this->FirstComponent, this->NbrComponent);
	    }
	  else
	    ((FermionOnSphereTwoLandauLevels *)this->FermionSpace)->BosonicStateTimeFermionicState(*(this->LLLVector), *(this->FermionVector), *(this->OutputVector),this->FinalStates,this->Weigth,
												   (BosonOnSphereShort*) this->LLLSpace, (FermionOnSphereTwoLandauLevels*)this->FinalSpace, 
												   this->FirstComponent,this->NbrComponent);
	}
    }
  else
    {
      if(this->Symmetry)
	((FermionOnSphereTwoLandauLevels *)this->FermionSpace)->LLLFermionicStateTimeFermionicStateSymmetric(*(this->LLLVector), *(this->FermionVector), *(this->OutputVector),
													     this->FinalStates,this->Weigth,(FermionOnSphere*)this->LLLSpace, 
													     (BosonOnSphereShort*) this->FinalSpace, this->FirstComponent, this->NbrComponent);
      else
	((FermionOnSphereTwoLandauLevels *)this->FermionSpace)->LLLFermionicStateTimeFermionicState(*(this->LLLVector), *(this->FermionVector), *(this->OutputVector),this->FinalStates,this->Weigth,(FermionOnSphere*)this->LLLSpace, (BosonOnSphereShort*)this->FinalSpace,this->FirstComponent,this->NbrComponent);
    }
  timeval TotalEndingTime;
  gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout << this->FirstComponent << " " <<  this->NbrComponent << " : " << Dt << "s" << endl;
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereMonomialsTimesSlaterProjectionOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = (int) this->NbrComponent / (this->NbrStage*architecture->GetNbrThreads());
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  char* OutputFileName = "temporary_projection_vector.vec";
  char* LogFile = "projection.dat";
  FQHESphereMonomialsTimesSlaterProjectionOperation** TmpOperations = new FQHESphereMonomialsTimesSlaterProjectionOperation* [architecture->GetNbrThreads()];
  for(int i = 0 ;i < architecture->GetNbrThreads() ;i++)
    {
      TmpOperations[i] = (FQHESphereMonomialsTimesSlaterProjectionOperation*) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  for (int i = 1; i < architecture->GetNbrThreads(); ++i)
    {
      TmpOperations[i]->SetOutputVector((RealVector*)this->OutputVector->EmptyClone(true));
    }
  for(int p = 0 ;p <this->NbrStage -1;p++)
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
      (*this->OutputVector).WriteVector(OutputFileName);
      ofstream File;
      File.open(LogFile, ios::binary | ios::out);
      File.precision(14);
      File << TmpFirstComponent;
      File.close();
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

bool FQHESphereMonomialsTimesSlaterProjectionOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__
   if (architecture->IsMasterNode())
     {
       this->OutputVector->ClearVector();
       RealVector* TmpVector = (RealVector*) this->OutputVector->EmptyClone(true);
       int Step = (int) this->NbrComponent / (this->NbrStage * architecture->GetNbrSlaveNodes());
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
