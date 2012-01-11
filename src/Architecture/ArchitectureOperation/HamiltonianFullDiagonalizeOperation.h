////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of hamiltonian full diagonalization operation           //
//                                                                            //
//                        last modification : 06/01/2012                      //
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


#ifndef HAMILTONIANFULLDIAGONALIZEOPERATION_H
#define HAMILTONIANFULLDIAGONALIZEOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"


class AbstractHamiltonian;
class Vector;


class HamiltonianFullDiagonalizeOperation: public AbstractArchitectureOperation
{

 protected:

  // index of the first component
  int FirstComponent;
  // number of component 
  int NbrComponent;

  // pointer to the hamiltonian
  AbstractHamiltonian* Hamiltonian;

  // execution time measured in RawApply
  double ExecutionTime;

  // true if the hamiltonian is complex
  bool ComplexFlag;
  // true if the eigenstates have to be computed
  bool EigenstateFlag ;

 public:
  
  // constructor 
  //
  // hamiltonian = pointer to the hamiltonian to use
  // complexFlag = true if the hamiltonian is complex
  // eigenstateFlag = true if the eigenstates have to be computed
  HamiltonianFullDiagonalizeOperation(AbstractHamiltonian* hamiltonian, bool complexFlag, bool eigenstateFlag);

  // copy constructor 
  //
  // operation = reference on operation to copy
  HamiltonianFullDiagonalizeOperation(const HamiltonianFullDiagonalizeOperation& operation);

  // constructor from a master node information
  //
  // hamiltonian = pointer to the hamiltonian to use
  // architecture = pointer to the distributed architecture to use for communications
  HamiltonianFullDiagonalizeOperation(AbstractHamiltonian* hamiltonian, SimpleMPIArchitecture* architecture);
  
  // destructor
  //
  ~HamiltonianFullDiagonalizeOperation();
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();

 protected:

  // apply operation for SimpleMPI architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture);
  
};

#endif
