////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                 Copyright (C) 2001-2002 Antoine Sterdyniak                 //
//                                                                            //
//                                                                            //
//                       class of computation operation                       //
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


#ifndef FQHESPHERE2LLBOSONICTIMESPOLARIZEDSLATERPROJECTIONOPERATION_H
#define FQHESPHERE2LLBOSONICTIMESPOLARIZEDSLATERPROJECTIONOPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractArchitectureOperation.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/BosonOnSphereTwoLandauLevels.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"



class FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation : public AbstractArchitectureOperation
{
 protected:
  
  // index of the first component
  int FirstComponent;
  
  // number of component 
  int NbrComponent;
  
  // pointer to the HilbertSpace after product and projection if there is one
  FermionOnSphereWithSpin * FinalSpace;
  
  // pointer to the initial Bosonic HilbertSpace
  ParticleOnSphere * InitialSpace;
  
  // pointer to the initial Fermionic HilbertSpace
  FermionOnSphere * FermionSpace;
  
  // RealVector where the result will be store
  RealVector * OutputVector;
  
  // RealVector where the bosonic state is stored
  RealVector * BosonicVector;
  
  // number of part in which the initial bosonic vector will be separated
  int NbrStage;
  
  // true if 2LL are considered
  bool TwoLandauLevels;
  
 public:
  
  // constructor 
  //
  // Space = pointer to the HilbertSpace to use
  // fileName = name of the file where the kostka number will be store
  // nbrLL = number of Landau levels
  FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation(ParticleOnSphere * initialSpace, FermionOnSphere * fermionSpace, FermionOnSphereWithSpin * finalSpace, RealVector* bosonicVector, RealVector* outputVector, bool twoLandauLevels , int nbrStage = 20);
  
  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation(const FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation & operation);
  
  // destructor
  //
  ~FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation();
  
 protected:
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // set destination vector 
  // 
  // vector where the result has to be stored
  void SetOutputVector (RealVector* outputVector);
  
  // get destination vector 
  // 
  // return value = pointer to destination vector
  RealVector* GetDestinationVector ();
  
  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
  // apply operation for SMP architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SMPArchitecture* architecture);
  
  // apply operation for SimpleMPI architecture
  //
  // architecture = pointer to the architecture
  // return value = true if no error occurs
  bool ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture);
  
};


#endif	
