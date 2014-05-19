////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                         class of particle on a torus                       //
//                                                                            //
//                        last modification : 18/07/2002                      //
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
#include "HilbertSpace/ParticleOnTorus.h"
#include "Vector/ComplexVector.h"
#include "Architecture/ArchitectureOperation/FQHETorusApplyCNRotationOperation.h"


// virtual destructor
//

ParticleOnTorus::~ParticleOnTorus ()
{
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Ky sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// kySector = Ky sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

RealSymmetricMatrix ParticleOnTorus::EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kySector, RealVector& groundState)
{
  RealSymmetricMatrix TmpDensityMatrix;
  return TmpDensityMatrix;
}
// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Ky sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// kySector = Ky sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix ParticleOnTorus::EvaluatePartialDensityMatrixParticlePartition (int nbrBosonSector, int kySector, ComplexVector& groundState)
{
  HermitianMatrix TmpDensityMatrix;
  return TmpDensityMatrix;
}

// apply a magnetic translation along x to a given state
//
// index = state index 
// return value = translated state index

int ParticleOnTorus::ApplyXMagneticTranslation(int index)
{
  return -1;
}
  
// remove part of each Fock state, discarding component if the Fock state does not a given pattern
//
// inputVector = state to truncate
// reducedSpace = Hilbert space where the truncated state will lie
// pattern = array describing the pattern 
// patternSize = pattern size
// patternShift = indicate where the pattern has to be applied
// return value = trucated state

RealVector ParticleOnTorus::TruncateStateWithPatternConstraint(RealVector& inputVector, ParticleOnTorus* reducedSpace, int* pattern, int patternSize, int patternShift)
{
  RealVector Tmp;
  return Tmp;
}

// apply the C4 rotation to a given state assumin it is an eigenstate of both kx and ky
//
// inputState = reference on the state that has to be rotated
// inputSpace = Hilbert space associated to the input state
// architecture = pointer to the architecture
// clockwise = the rotation is done clockwise
// return value = rotated state

ComplexVector ParticleOnTorus::C4Rotation (ComplexVector& inputState, ParticleOnTorus* inputSpace, bool clockwise, AbstractArchitecture* architecture)
{
  ComplexVector OutputState (this->HilbertSpaceDimension, true);
  if (architecture != 0)
    {
      if (clockwise == false)
	{
	  FQHETorusApplyCNRotationOperation Operation(4, &inputState, &OutputState, inputSpace, this);
	  Operation.ApplyOperation(architecture);
	}
      else
	{
	  FQHETorusApplyCNRotationOperation Operation(-4, &inputState, &OutputState, inputSpace, this);
	  Operation.ApplyOperation(architecture);
	}
    }
  else
    {
      this->CoreC4Rotation(inputState, inputSpace, OutputState, 0 , this->HilbertSpaceDimension, clockwise);
    }
  return OutputState;
}

// core part of the C4 rotation
//
// inputState = reference on the state that has to be rotated
// inputSpace = Hilbert space associated to the input state
// outputState = reference on the rotated state
// minIndex = minimum index that has to be computed
// nbrIndices = number of indices that have to be computed
// clockwise = the rotation is done clockwise
// return value = reference on the rotated state

ComplexVector& ParticleOnTorus::CoreC4Rotation (ComplexVector& inputState, ParticleOnTorus* inputSpace, ComplexVector& outputState, int minIndex, int nbrIndices, bool clockwise)
{
  return outputState;
}
