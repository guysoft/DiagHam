////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of FQHE particle entanglement spectrum parallelization operation    //
//                                                                            //
//                        last modification : 15/12/2010                      //
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
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"


// constructor 
//
// space = pointer to the Hilbert space to use
// invAlpha = inverse of the Jack polynomial alpha coefficient
// rootPartition = root partition (in fermionic binary representation)
// indexArray = array where state indices are stored
// stateArray = array use to store computed state description
// componentArray = array where computed component numerical factors are stored
// rhoArray = rho factor associated to each state
// nbrComputedComponentArray = number of connected components associated to each state through the Jack generator
// fermionicFlag = true if we are dealing with fermions
// symmetricFlag = true if the state is Lz<->-Lz symmetric

FQHESphereParticleEntanglementSpectrumOperation::FQHESphereParticleEntanglementSpectrumOperation (ParticleOnSphere* fullSpace, ParticleOnSphere* destinationSpace, ParticleOnSphere* complementarySpace, RealVector& groundState, RealSymmetricMatrix& densityMatrix)
{
  this->FullSpace  = (ParticleOnSphere*) fullSpace->Clone();
  this->DestinationHilbertSpace = (ParticleOnSphere*) destinationSpace->Clone();
  this->ComplementaryHilbertSpace = (ParticleOnSphere*) complementarySpace->Clone();
  this->GroundState = groundState;
  this->DensityMatrix = densityMatrix;
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = this->ComplementaryHilbertSpace->GetLargeHilbertSpaceDimension();
  this->OperationType = AbstractArchitectureOperation::FQHESphereParticleEntanglementSpectrumOperation;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}

// copy constructor 
//
// operation = reference on operation to copy

FQHESphereParticleEntanglementSpectrumOperation::FQHESphereParticleEntanglementSpectrumOperation(const FQHESphereParticleEntanglementSpectrumOperation& operation)
{
  this->FullSpace  = (ParticleOnSphere*) operation.FullSpace->Clone();
  this->DestinationHilbertSpace = (ParticleOnSphere*) operation.DestinationHilbertSpace->Clone();
  this->ComplementaryHilbertSpace = (ParticleOnSphere*) operation.ComplementaryHilbertSpace->Clone();
  this->GroundState = operation.GroundState;
  this->DensityMatrix = operation.DensityMatrix;
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->LargeFirstComponent = operation.LargeFirstComponent;
  this->LargeNbrComponent = operation.LargeNbrComponent;
  this->OperationType = AbstractArchitectureOperation::FQHESphereParticleEntanglementSpectrumOperation;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}
  
// destructor
//

FQHESphereParticleEntanglementSpectrumOperation::~FQHESphereParticleEntanglementSpectrumOperation()
{
  if (this->LocalOperations != 0)
    {
      for (int i = 0; i < this->NbrLocalOperations; ++i)
	{
	  delete this->LocalOperations[i];
	}
      delete[] this->LocalOperations;
    }
  delete this->FullSpace;
  delete this->DestinationHilbertSpace;
  delete this->ComplementaryHilbertSpace;
}
  
// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* FQHESphereParticleEntanglementSpectrumOperation::Clone()
{
  return new FQHESphereParticleEntanglementSpectrumOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHESphereParticleEntanglementSpectrumOperation::RawApplyOperation()
{
  this->FullSpace->EvaluatePartialDensityMatrixParticlePartitionCore(this->LargeFirstComponent, this->LargeNbrComponent, this->ComplementaryHilbertSpace, this->DestinationHilbertSpace, this->GroundState, &this->DensityMatrix);
  return true;
}


// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHESphereParticleEntanglementSpectrumOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  long Step = this->LargeNbrComponent / ((long) architecture->GetNbrThreads());
  long TmpFirstComponent = this->LargeFirstComponent;
  if (this->LocalOperations == 0)
    {
      this->NbrLocalOperations = architecture->GetNbrThreads();
      this->LocalOperations = new FQHESphereParticleEntanglementSpectrumOperation* [this->NbrLocalOperations];
      for (int i = 0; i < this->NbrLocalOperations; ++i)
	this->LocalOperations[i] = (FQHESphereParticleEntanglementSpectrumOperation*) this->Clone();
    }
  int ReducedNbrThreads = this->NbrLocalOperations - 1;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->LocalOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
      this->LocalOperations[i]->DensityMatrix = RealSymmetricMatrix(this->DestinationHilbertSpace->GetHilbertSpaceDimension(), true);
      architecture->SetThreadOperation(this->LocalOperations[i], i);
      TmpFirstComponent += Step;
    }
  this->LocalOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->LargeNbrComponent + this->LargeFirstComponent - TmpFirstComponent);  
  architecture->SetThreadOperation(this->LocalOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();
  for (int i = 0; i < ReducedNbrThreads; ++i)
    this->LocalOperations[ReducedNbrThreads]->DensityMatrix += this->LocalOperations[i]->DensityMatrix;
  return true;
}
  
