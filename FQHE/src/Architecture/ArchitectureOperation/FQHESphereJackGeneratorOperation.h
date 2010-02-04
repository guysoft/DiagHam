////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of FQHE Jack generator parallelization operation         //
//                                                                            //
//                        last modification : 02/02/2010                      //
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


#ifndef FQHESPHEREJACKGENERATOROPERATION_H
#define FQHESPHEREJACKGENERATOROPERATION_H


#include "config.h"
#include "Architecture/ArchitectureOperation/AbstractPrecalculationOperation.h"
#include "HilbertSpace/ParticleOnSphere.h"


class AbstractFunctionBasis;
class RealVector;


class FQHESphereJackGeneratorOperation: public AbstractPrecalculationOperation
{

 protected:

  // pointer to the Hilbert space
  ParticleOnSphere* HilbertSpace;

  // inverse of the Jack polynomial alpha coefficient
  double InvAlpha;
  // root partition (in fermionic binary representation)
  unsigned long RootPartition;

  // local shift to apply to access array elements
  long LocalShift;
  // array where state indices are stored
  long** IndexArray;
  // array use to store computed state description
  unsigned long** StateArray; 
  // array where computed component numerical factors are stored
  double** ComponentArray;
  // rho factor associated to each state
  double* RhoArray;
  // number of connected components associated to each state through the Jack generator
  int* NbrComputedComponentArray;

 public:
  
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
  FQHESphereJackGeneratorOperation(ParticleOnSphere* space, double invAlpha, unsigned long rootPartition, long** indexArray, unsigned long** stateArray, double** componentArray, double* rhoArray, int* nbrComputedComponentArray);

  // copy constructor 
  //
  // operation = reference on operation to copy
  FQHESphereJackGeneratorOperation(const FQHESphereJackGeneratorOperation& operation);
  
  // destructor
  //
  ~FQHESphereJackGeneratorOperation();
  
  // clone operation
  //
  // return value = pointer to cloned operation
  AbstractArchitectureOperation* Clone();
  
  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const int& firstComponent, const int& nbrComponent);

  // set range of indices
  // 
  // firstComponent = index of the first component
  // nbrComponent = number of component
  void SetIndicesRange (const long& firstComponent, const long& nbrComponent);

  // get hilbert space dimension
  // 
  // return value = hilbert space dimension  
  int GetHilbertSpaceDimension ();

 protected:

  // apply operation (architecture independent)
  //
  // return value = true if no error occurs
  bool RawApplyOperation();
  
};

// get hilbert space dimension
// 
// return value = hilbert space dimension

inline int FQHESphereJackGeneratorOperation::GetHilbertSpaceDimension ()
{
  return this->HilbertSpace->GetHilbertSpaceDimension();
}

#endif
