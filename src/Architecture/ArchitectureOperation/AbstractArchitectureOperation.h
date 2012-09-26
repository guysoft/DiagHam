////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of Abstract architecture operation                 //
//                                                                            //
//                        last modification : 23/10/2002                      //
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


#ifndef ABSTRACTARCHITECTUREOPERATION_H
#define ABSTRACTARCHITECTUREOPERATION_H


#include "config.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"

class MonoProcessorArchitecture;


class AbstractArchitectureOperation
{

 protected:

  int OperationType;

 public:
  
  enum Operation
    {
      VectorHamiltonianMultiply = 0x1,
      AddRealLinearCombination = 0x2,
      MultipleRealScalarProduct = 0x4,      
      MatrixMatrixMultiply = 0x8,
      AddComplexLinearCombination = 0x10,
      MultipleComplexScalarProduct = 0x20,
      MultipleVectorHamiltonianMultiply = 0x40,
      HamiltonianFullDiagonalize = 0x41,
      Generic = 0x100,
      HamiltonianPrecalculation = 0x200,
      QHEOperation = 0x10000,
      SpinOperation = 0x20000,
      FTIOperation = 0x30000,
      QHEParticlePrecalculation = 0x200,
      NDMAPPrecalculation = 0x400,
      ScalarSum = 0x800,
      OperatorMatrixElement = 0x801,
      QHEParticleWaveFunction = 0x10800,
      MainTask = 0x1000,
      VectorOperatorMultiply=0x2000,
      QHEComplexVectorConversion=0x14000,
      FQHESphereJackGenerator = 0x18000,
      FQHESphereJastrowMultiplication = 0x18001,
      FQHESphereJastrowDivision = 0x18002,
      FQHESphereMonomialsTimesSlaterProjection = 0x18003,
      FQHESphereParticleEntanglementSpectrumOperation = 0x18004,
      FQHESphereMultipleMonomialsTimesSlaterProjection = 0x18005,
      FQHETorusParticleEntanglementSpectrumOperation = 0x18006,
      FQHELatticeParticleEntanglementSpectrumOperation = 0x18007,
      FQHESphereMonomialsProductOperation = 0x18017,
      FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation = 0x18018,
      FQHESphereSymmetrizeU1U1StateOperation = 0x18020,
      FQHESquareLatticeSymmetrizeU1U1StateOperation = 0x18021,
      FQHESphereBosonsWithSpinLandauLevelLiftOperation = 0x18022,
      FQHELatticeFourierTransformOperation = 0x18024,
      FTIComputeBandStructureOperation = 0x30001
    };
  
  // destructor
  //
  virtual ~AbstractArchitectureOperation();
  
  // clone operation
  //
  // return value = pointer to cloned operation
  virtual AbstractArchitectureOperation* Clone() = 0;
  
  // apply operation for a given architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ApplyOperation(AbstractArchitecture* architecture);
  
  // apply an SMP round robin operation 
  //
  // return value = true if no error occurs
  virtual bool ApplyOperationSMPRoundRobin(SMPArchitecture* architecture, int threadID);
  
  // get operation type
  //
  // return value = code corresponding to the operation
  int GetOperationType ();

  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  virtual bool RawApplyOperation() = 0;
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  
  
 protected:

  // apply operation for mono processor architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(MonoProcessorArchitecture* architecture);
  
  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // apply operation for simple MPI architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  virtual bool ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture);
  
};


// get operation type
//
// return value = code corresponding to the operation

inline int AbstractArchitectureOperation::GetOperationType ()
{
  return this->OperationType;
}

#endif
