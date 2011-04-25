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

#include <sys/time.h>


// constructor 
//
// fullSpace = pointer to the full Hilbert space to use
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored

FQHESphereParticleEntanglementSpectrumOperation::FQHESphereParticleEntanglementSpectrumOperation (ParticleOnSphere* fullSpace, ParticleOnSphere* destinationSpace, ParticleOnSphere* complementarySpace, RealVector& groundState, RealSymmetricMatrix& densityMatrix)
{
  this->FullSpace  = (ParticleOnSphere*) fullSpace->Clone();
  this->DestinationHilbertSpace = (ParticleOnSphere*) destinationSpace->Clone();
  this->ComplementaryHilbertSpace = (ParticleOnSphere*) complementarySpace->Clone();
  this->GroundState = groundState;
  this->DensityMatrix = densityMatrix;
  this->IncompleteBetaThetaTop = 0; 
  this->IncompleteBetaThetaBottom = 0; 
  this->PhiRange = 0.0;
  this->NbrNonZeroElements = 0l;
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = this->ComplementaryHilbertSpace->GetLargeHilbertSpaceDimension();
  this->OperationType = AbstractArchitectureOperation::FQHESphereParticleEntanglementSpectrumOperation;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}

// constructor 
//
// fullSpace = pointer to the full Hilbert space to use
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored

FQHESphereParticleEntanglementSpectrumOperation::FQHESphereParticleEntanglementSpectrumOperation(ParticleOnSphere* fullSpace, ParticleOnSphere* destinationSpace, ParticleOnSphere* complementarySpace, ComplexVector& groundState, HermitianMatrix& densityMatrix)
{
  this->FullSpace  = (ParticleOnSphere*) fullSpace->Clone();
  this->DestinationHilbertSpace = (ParticleOnSphere*) destinationSpace->Clone();
  this->ComplementaryHilbertSpace = (ParticleOnSphere*) complementarySpace->Clone();
  this->ComplexGroundState = groundState;
  this->ComplexDensityMatrix = densityMatrix;
  this->IncompleteBetaThetaTop = 0; 
  this->IncompleteBetaThetaBottom = 0; 
  this->PhiRange = 0.0;
  this->NbrNonZeroElements = 0l;
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = this->ComplementaryHilbertSpace->GetLargeHilbertSpaceDimension();
  this->OperationType = AbstractArchitectureOperation::FQHESphereParticleEntanglementSpectrumOperation;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}

// constructor 
//
// fullSpace = pointer to the full Hilbert space to use
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// incompleteBetaThetaTop = pointer to the array where the top part coefficients are stored
// incompleteBetaThetaBotton = pointer on the pointer to the array where the bottom part coefficients are stored
// phiRange = The angle traced in the \hat{phi} direction between the 2 longitudes defining the cut in degrees

FQHESphereParticleEntanglementSpectrumOperation::FQHESphereParticleEntanglementSpectrumOperation(ParticleOnSphere* fullSpace, ParticleOnSphere* destinationSpace, ParticleOnSphere* complementarySpace, RealVector& groundState, RealSymmetricMatrix& densityMatrix, double* incompleteBetaThetaBottom, double* incompleteBetaThetaTop, double phiRange)
{
  this->FullSpace  = (ParticleOnSphere*) fullSpace->Clone();
  this->DestinationHilbertSpace = (ParticleOnSphere*) destinationSpace->Clone();
  this->ComplementaryHilbertSpace = (ParticleOnSphere*) complementarySpace->Clone();
  this->GroundState = groundState;
  this->DensityMatrix = densityMatrix;
  this->IncompleteBetaThetaTop = incompleteBetaThetaTop; 
  this->IncompleteBetaThetaBottom = incompleteBetaThetaBottom; 
  this->PhiRange = phiRange;
  this->NbrNonZeroElements = 0l;
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
  this->ComplexGroundState = operation.ComplexGroundState;
  this->ComplexDensityMatrix = operation.ComplexDensityMatrix;
  this->IncompleteBetaThetaTop = operation.IncompleteBetaThetaTop; 
  this->IncompleteBetaThetaBottom = operation.IncompleteBetaThetaBottom; 
  this->PhiRange = operation.PhiRange;
  this->NbrNonZeroElements = operation.NbrNonZeroElements;
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
  if (this->IncompleteBetaThetaTop == 0)
    {
      if (this->ComplexGroundState.GetLargeVectorDimension() == 0l)
	{
	  this->NbrNonZeroElements = this->FullSpace->EvaluatePartialDensityMatrixParticlePartitionCore(this->LargeFirstComponent, this->LargeNbrComponent, this->ComplementaryHilbertSpace, this->DestinationHilbertSpace, this->GroundState, &this->DensityMatrix);
	}
      else
	{
	  this->NbrNonZeroElements = this->FullSpace->EvaluatePartialDensityMatrixParticlePartitionCore(this->LargeFirstComponent, this->LargeNbrComponent, this->ComplementaryHilbertSpace, this->DestinationHilbertSpace, this->ComplexGroundState, &this->ComplexDensityMatrix);
	}
    }
  else
    {
      this->NbrNonZeroElements = this->FullSpace->EvaluatePartialDensityMatrixRealSpacePartitionCore(this->LargeFirstComponent, this->LargeNbrComponent, this->ComplementaryHilbertSpace, this->DestinationHilbertSpace, this->GroundState, &this->DensityMatrix, this->IncompleteBetaThetaBottom, this->IncompleteBetaThetaTop, this->PhiRange);
    }
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
  if (this->ComplexGroundState.GetLargeVectorDimension() == 0l)
    {
      for (int i = 0; i < ReducedNbrThreads; ++i)
	{
	  this->LocalOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
	  this->LocalOperations[i]->DensityMatrix = RealSymmetricMatrix(this->DestinationHilbertSpace->GetHilbertSpaceDimension(), true);
	  architecture->SetThreadOperation(this->LocalOperations[i], i);
	  TmpFirstComponent += Step;
	}
    }
  else
    {
      for (int i = 0; i < ReducedNbrThreads; ++i)
	{
	  this->LocalOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
	  this->LocalOperations[i]->ComplexDensityMatrix = HermitianMatrix(this->DestinationHilbertSpace->GetHilbertSpaceDimension(), true);
	  architecture->SetThreadOperation(this->LocalOperations[i], i);
	  TmpFirstComponent += Step;
	}
    }
  this->LocalOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->LargeNbrComponent + this->LargeFirstComponent - TmpFirstComponent);  
  architecture->SetThreadOperation(this->LocalOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();
  if (this->ComplexGroundState.GetLargeVectorDimension() == 0l)
    {
      for (int i = 0; i < ReducedNbrThreads; ++i)
	{
	  this->LocalOperations[ReducedNbrThreads]->DensityMatrix += this->LocalOperations[i]->DensityMatrix;
	  this->LocalOperations[ReducedNbrThreads]->NbrNonZeroElements += this->LocalOperations[i]->NbrNonZeroElements;
	}
    }
  else
    {
      for (int i = 0; i < ReducedNbrThreads; ++i)
	{
	  this->LocalOperations[ReducedNbrThreads]->ComplexDensityMatrix += this->LocalOperations[i]->ComplexDensityMatrix;
	  this->LocalOperations[ReducedNbrThreads]->NbrNonZeroElements += this->LocalOperations[i]->NbrNonZeroElements;
	}
    }
  this->NbrNonZeroElements = this->LocalOperations[ReducedNbrThreads]->NbrNonZeroElements;
  return true;
}
  
