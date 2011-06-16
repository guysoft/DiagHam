////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//  class of FQHE particle entanglement spectrum parallelization operation    //
//                    for torus with magnetic translations                    //
//                                                                            //
//                        last modification : 13/06/2011                      //
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
#include "Architecture/ArchitectureOperation/FQHETorusParticleEntanglementSpectrumOperation.h"

#include <sys/time.h>


// constructor 
//
// fullSpace = pointer to the full Hilbert space to use
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored

FQHETorusParticleEntanglementSpectrumOperation::FQHETorusParticleEntanglementSpectrumOperation(ParticleOnTorusWithMagneticTranslations* fullSpace, ParticleOnTorusWithMagneticTranslations* destinationSpace, ParticleOnTorusWithMagneticTranslations* complementarySpace, ComplexVector& groundState, HermitianMatrix& densityMatrix)
{
  this->FullSpace  = (ParticleOnTorusWithMagneticTranslations*) fullSpace->Clone();
  this->DestinationHilbertSpace = (ParticleOnTorusWithMagneticTranslations*) destinationSpace->Clone();
  this->ComplementaryHilbertSpace = (ParticleOnTorusWithMagneticTranslations*) complementarySpace->Clone();
  this->ComplexGroundState = groundState;
  this->ComplexDensityMatrix = densityMatrix;
  this->NbrNonZeroElements = 0l;
  this->FirstComponent = 0;
  this->NbrComponent = 0;
  this->LargeFirstComponent = 0;
  this->LargeNbrComponent = this->ComplementaryHilbertSpace->GetLargeHilbertSpaceDimension();
  this->OperationType = AbstractArchitectureOperation::FQHETorusParticleEntanglementSpectrumOperation;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}

// copy constructor 
//
// operation = reference on operation to copy

FQHETorusParticleEntanglementSpectrumOperation::FQHETorusParticleEntanglementSpectrumOperation(const FQHETorusParticleEntanglementSpectrumOperation& operation)
{
  this->FullSpace  = (ParticleOnTorusWithMagneticTranslations*) operation.FullSpace->Clone();
  this->DestinationHilbertSpace = (ParticleOnTorusWithMagneticTranslations*) operation.DestinationHilbertSpace->Clone();
  this->ComplementaryHilbertSpace = (ParticleOnTorusWithMagneticTranslations*) operation.ComplementaryHilbertSpace->Clone();
  this->ComplexGroundState = operation.ComplexGroundState;
  this->ComplexDensityMatrix = operation.ComplexDensityMatrix;
  this->NbrNonZeroElements = operation.NbrNonZeroElements;
  this->FirstComponent = operation.FirstComponent;
  this->NbrComponent = operation.NbrComponent;
  this->LargeFirstComponent = operation.LargeFirstComponent;
  this->LargeNbrComponent = operation.LargeNbrComponent;
  this->OperationType = AbstractArchitectureOperation::FQHETorusParticleEntanglementSpectrumOperation;
  this->LocalOperations = 0;
  this->NbrLocalOperations = 0;
}
  
// destructor
//

FQHETorusParticleEntanglementSpectrumOperation::~FQHETorusParticleEntanglementSpectrumOperation()
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

AbstractArchitectureOperation* FQHETorusParticleEntanglementSpectrumOperation::Clone()
{
  return new FQHETorusParticleEntanglementSpectrumOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool FQHETorusParticleEntanglementSpectrumOperation::RawApplyOperation()
{
  this->NbrNonZeroElements = this->FullSpace->EvaluatePartialDensityMatrixParticlePartitionCore(this->LargeFirstComponent, this->LargeNbrComponent, this->ComplementaryHilbertSpace, this->DestinationHilbertSpace, this->ComplexGroundState, &this->ComplexDensityMatrix);
  return true;
}


// apply operation for SMP architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool FQHETorusParticleEntanglementSpectrumOperation::ArchitectureDependentApplyOperation(SMPArchitecture* architecture)
{
  long Step = this->LargeNbrComponent / ((long) architecture->GetNbrThreads());
  long TmpFirstComponent = this->LargeFirstComponent;
  if (this->LocalOperations == 0)
    {
      this->NbrLocalOperations = architecture->GetNbrThreads();
      this->LocalOperations = new FQHETorusParticleEntanglementSpectrumOperation* [this->NbrLocalOperations];
      for (int i = 0; i < this->NbrLocalOperations; ++i)
	this->LocalOperations[i] = (FQHETorusParticleEntanglementSpectrumOperation*) this->Clone();
    }
  int ReducedNbrThreads = this->NbrLocalOperations - 1;
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->LocalOperations[i]->SetIndicesRange(TmpFirstComponent, Step);
	  this->LocalOperations[i]->ComplexDensityMatrix = HermitianMatrix(this->DestinationHilbertSpace->GetHilbertSpaceDimension(), true);
	  architecture->SetThreadOperation(this->LocalOperations[i], i);
	  TmpFirstComponent += Step;
    }
  this->LocalOperations[ReducedNbrThreads]->SetIndicesRange(TmpFirstComponent, this->LargeNbrComponent + this->LargeFirstComponent - TmpFirstComponent);  
  architecture->SetThreadOperation(this->LocalOperations[ReducedNbrThreads], ReducedNbrThreads);
  architecture->SendJobs();
  for (int i = 0; i < ReducedNbrThreads; ++i)
    {
      this->LocalOperations[ReducedNbrThreads]->ComplexDensityMatrix += this->LocalOperations[i]->ComplexDensityMatrix;
      this->LocalOperations[ReducedNbrThreads]->NbrNonZeroElements += this->LocalOperations[i]->NbrNonZeroElements;
    }
  this->NbrNonZeroElements = this->LocalOperations[ReducedNbrThreads]->NbrNonZeroElements;
  return true;
}
  
