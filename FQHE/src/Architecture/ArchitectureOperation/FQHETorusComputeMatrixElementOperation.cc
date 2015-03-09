////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2002 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//            class of operations that compute the matrix elements            //
//                  of the many-body interaction on the torus                 // 
//                                                                            //
//                        last modification : 06/03/2015                      //
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
#include "Architecture/ArchitectureOperation/FQHETorusComputeMatrixElementOperation.h"


// constructor 
//
// hamiltonian = pointer to the generic n-body Hamiltonian
// nbrUniqueMatrixElements = number of unique matrix elements
// momentumSectorIndices = array that contains the momentum sector of each matrix element
// j1Indices = array that contains the creation indices of each matrix element
// j2Indices = array that contains the annihilation indices of each matrix element
// matrixElements = array where the matrix elements will be stored
// nbrMPIStage = number of stages in which the calculation has to be splitted in MPI mode
// nbrSMPStage = number of stages in which the calculation has to be splitted in SMP mode

FQHETorusComputeMatrixElementOperation::FQHETorusComputeMatrixElementOperation(ParticleOnTorusGenericNBodyWithMagneticTranslationsHamiltonian* hamiltonian, long nbrUniqueMatrixElements,
									       int* momentumSectorIndices, int* j1Indices, int* j2Indices, Complex* matrixElements,
									       int nbrMPIStage, int nbrSMPStage)
{
  this->Hamiltonian = hamiltonian;
  this->NbrUniqueMatrixElements = nbrUniqueMatrixElements;
  this->MomentumSectorIndices = momentumSectorIndices;
  this->J1Indices = j1Indices;
  this->J2Indices = j2Indices;
  this->MatrixElements = matrixElements;

  this->FirstComponent = 0;
  this->NbrComponent = (int) nbrUniqueMatrixElements;
  this->LargeFirstComponent = 0l;
  this->LargeNbrComponent = nbrUniqueMatrixElements;
  this->OperationType = AbstractArchitectureOperation::FQHETorusComputeMatrixElementOperation;
  this->NbrMPIStage = nbrMPIStage;
  this->NbrSMPStage = nbrSMPStage;
}
    
// copy constructor 
//
// operation = reference on operation to copy

FQHETorusComputeMatrixElementOperation::FQHETorusComputeMatrixElementOperation(const FQHETorusComputeMatrixElementOperation & operation)
{
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->LargeFirstComponent = operation.LargeFirstComponent;
  this->LargeNbrComponent = operation.LargeNbrComponent;
  this->OperationType = AbstractArchitectureOperation::FQHETorusComputeMatrixElementOperation;
  this->Hamiltonian = operation.Hamiltonian;
  this->NbrUniqueMatrixElements = operation.NbrUniqueMatrixElements;
  this->MomentumSectorIndices = operation.MomentumSectorIndices;
  this->J1Indices = operation.J1Indices;
  this->J2Indices = operation.J2Indices;
  this->MatrixElements = operation.MatrixElements;
  this->NbrMPIStage = operation.NbrMPIStage;
  this->NbrSMPStage = operation.NbrSMPStage;
  this->SMPStages = operation.SMPStages;
}
  
// destructor
//

FQHETorusComputeMatrixElementOperation::~FQHETorusComputeMatrixElementOperation()
{
}
  
// set range of indices
//
// firstComponent = index of the first component
// nbrComponent = number of component

void FQHETorusComputeMatrixElementOperation::SetIndicesRange (const long& firstComponent, const long& nbrComponent)
{
  this->LargeFirstComponent = firstComponent;
  this->LargeNbrComponent = nbrComponent;
  this->FirstComponent = (int) firstComponent;
  this->NbrComponent = (int) nbrComponent;
}

// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHETorusComputeMatrixElementOperation::Clone()
{
  return new FQHETorusComputeMatrixElementOperation (*this);
}

// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHETorusComputeMatrixElementOperation::RawApplyOperation()
{
  int LastIndex = this->FirstComponent + this->NbrComponent;
  int* TmpNIndices =  new int [this->Hamiltonian->NBodyValue];
  int* TmpMIndices =  new int [this->Hamiltonian->NBodyValue];
  if (this->Hamiltonian->Particles->GetParticleStatistic() == ParticleOnTorus::FermionicStatistic)
    {
      int NbrPermutations = 1;
      for (int i = 1; i <= this->Hamiltonian->NBodyValue; ++i)
	NbrPermutations *= i;
      int** Permutations = new int*[NbrPermutations]; 
      double* PermutationSign = new double[NbrPermutations]; 
      Permutations[0] = new int [this->Hamiltonian->NBodyValue];
      for (int i = 0; i < this->Hamiltonian->NBodyValue; ++i)
	Permutations[0][i] = i;
      PermutationSign[0] = 1.0;
      double TmpSign = 1.0;
      for (int i = 1; i < NbrPermutations; ++i)
	{
	  Permutations[i] = new int [this->Hamiltonian->NBodyValue];
	  for (int j = 0; j < this->Hamiltonian->NBodyValue; ++j)
	    Permutations[i][j] = Permutations[i - 1][j];
	  int* TmpArrayPerm = Permutations[i];
	  int Pos1 = this->Hamiltonian->NBodyValue - 1;
	  while (TmpArrayPerm[Pos1 - 1] >= TmpArrayPerm[Pos1])
	    --Pos1;
	  --Pos1;
	  int Pos2 = this->Hamiltonian->NBodyValue - 1;      
	  while (TmpArrayPerm[Pos2] <= TmpArrayPerm[Pos1])
	    --Pos2;
	  int TmpIndex = TmpArrayPerm[Pos1];
	  TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	  TmpArrayPerm[Pos2] = TmpIndex;
	  TmpSign *= -1.0;
	  Pos2 = this->Hamiltonian->NBodyValue - 1;   
	  Pos1++;
	  while (Pos1 < Pos2)
	    {
	      TmpIndex = TmpArrayPerm[Pos1];
	      TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	      TmpArrayPerm[Pos2] = TmpIndex;
	      ++Pos1;
	      --Pos2;
	      TmpSign *= -1.0;
	    }
	  PermutationSign[i] = TmpSign;
	}

      for (int Index = this->FirstComponent;  Index < LastIndex; ++Index)
	{
	  int J1 = this->J1Indices[Index];
	  int J2 = this->J2Indices[Index];
	  int MomentumSector = this->MomentumSectorIndices[Index];
	  double TmpInteraction = 0.0;
	  for (int l1 = 0; l1 < NbrPermutations; ++l1)
	    {
	      int* TmpPerm1 = Permutations[l1];
	      for (int k = 0; k < this->Hamiltonian->NBodyValue; ++k)
		{
		  TmpNIndices[k]  = this->Hamiltonian->NBodySectorIndicesPerSum[Index][(J1 * this->Hamiltonian->NBodyValue) + TmpPerm1[k]];
		}
	      for (int l2 = 0; l2 < NbrPermutations; ++l2)
		{
		  int* TmpPerm2 = Permutations[l2];
		  for (int k = 0; k < this->Hamiltonian->NBodyValue; ++k)
		    {
		      TmpMIndices[k] = this->Hamiltonian->NBodySectorIndicesPerSum[MomentumSector][(J2 * this->Hamiltonian->NBodyValue) + TmpPerm2[k]];
		    }
		  double Tmp = this->Hamiltonian->EvaluateInteractionCoefficient(TmpMIndices, TmpNIndices);
		  TmpInteraction += PermutationSign[l1] * PermutationSign[l2] * Tmp;
		}
	    }
	}
     delete[] PermutationSign;
      for (int i = 0; i < NbrPermutations; ++i)
	{
	  delete[] Permutations[i];
	}
      delete[] Permutations;
    }
  else
    {
    }
  delete[] TmpNIndices;
  delete[] TmpMIndices;
   return true;
}

// apply operation for SMP using round robin scheduling
//
//  architecture = instance of architecture class
// return value = true if no error occurs

bool FQHETorusComputeMatrixElementOperation::ApplyOperationSMPRoundRobin(SMPArchitecture* architecture, int threadID)
{
  int TmpNbrComponents = this->NbrComponent;
  int TmpFirstComponent = this->FirstComponent;
  
  int NbrStages = this->NbrSMPStage * architecture->GetNbrThreads();
  int StageIndex = 0;
   bool LockFlag = false;
   while (StageIndex < NbrStages) 
     {
       if (LockFlag == false)   
         {
           architecture->LockMutex();
           LockFlag = true;
         }
       StageIndex = this->SMPStages[0];
       if (StageIndex < NbrStages) 
         {         
           this->SMPStages[0]++;   
           architecture->UnLockMutex();
           LockFlag = false;
           this->SetIndicesRange(TmpFirstComponent + this->GetRankChunkStart(TmpNbrComponents, StageIndex,  NbrStages),  
				 this->GetRankChunkSize(TmpNbrComponents, StageIndex,  NbrStages));
           this->RawApplyOperation();
           ++StageIndex;
         }
     }
   if (LockFlag == true)
     {
       architecture->UnLockMutex();
       LockFlag = false;
     }
  return true;
}

// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHETorusComputeMatrixElementOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  int Step = this->NbrComponent / architecture->GetNbrThreads();
  if (Step == 0)
    Step = this->NbrComponent;
  this->SMPStages[0] = 0;
  int TotalNbrComponent = this->FirstComponent + this->NbrComponent;
  int TmpFirstComponent = this->FirstComponent;
  int ReducedNbrThreads = architecture->GetNbrThreads() - 1;
  FQHETorusComputeMatrixElementOperation** TmpOperations = new  FQHETorusComputeMatrixElementOperation* [architecture->GetNbrThreads()];
  
   
  for( int i = 0; i <  architecture->GetNbrThreads() ; i++)
    {
      TmpOperations[i] = (FQHETorusComputeMatrixElementOperation *) this->Clone();
      architecture->SetThreadOperation(TmpOperations[i], i);
    }
  
  architecture->SendJobsRoundRobin();

  for (int i = 0; i < architecture->GetNbrThreads(); i++)
    {
      delete TmpOperations[i];
    }
  delete[] TmpOperations;
  return true;
}
  
// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHETorusComputeMatrixElementOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __MPI__    
//   int Step = this->NbrComponent / architecture->GetNbrNodes();
//   int TmpFirstComponent = this->FirstComponent + (Step * architecture->GetNodeNbr());
//   int TmpNbrComponent = Step;
//   if ((architecture->GetNodeNbr() + 1) == architecture->GetNbrNodes())
//     {
//       TmpNbrComponent += this->NbrComponent % Step;
//     }
//   this->SetIndicesRange(TmpFirstComponent, TmpNbrComponent); 
//   switch (architecture->GetArchitectureID())
//     {	 
//     case AbstractArchitecture::MixedMPISMP:
//       this->ArchitectureDependentApplyOperation((SMPArchitecture*) (architecture->GetLocalArchitecture())); 
//       break;
//     default:
//       this->RawApplyOperation();
//       break;
//     }		
//   MPI::COMM_WORLD.Barrier();
//   if (this->OutputState != 0)
//     architecture->SumVector(*(this->OutputState));	
//   else
//     architecture->SumVector(*(this->ComplexOutputState));	
  return true;
#else
  return this->RawApplyOperation();
#endif
}
  
