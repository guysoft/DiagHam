////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of quantum Hall hamiltonian                     //
//                                                                            //
//                        last modification : 03/07/2003                      //
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


#ifndef ABSTRACTQHEHAMILTONIAN_H
#define ABSTRACTQHEHAMILTONIAN_H


#include "config.h"
#include "Hamiltonian/AbstractHamiltonian.h"

#include <iostream>


using std::ostream;


class AbstractArchitecture;


class AbstractQHEHamiltonian : public AbstractHamiltonian
{

  friend class QHEParticlePrecalculationOperation;

 protected:
  
  // architecture used for precalculation
  AbstractArchitecture* Architecture;

 public:

  // destructor
  //
  virtual ~AbstractQHEHamiltonian();

  // save precalculations in a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be stored
  // return value = true if no error occurs
  virtual bool SavePrecalculation (char* fileName) = 0;

 protected:
 
  // evaluate all interaction factors
  //   
  virtual void EvaluateInteractionFactors() = 0;

  // test the amount of memory needed for fast multiplication algorithm
  //
  // allowedMemory = amount of memory that cam be allocated for fast multiplication
  // return value = amount of memory needed
  virtual int FastMultiplicationMemory(int allowedMemory) = 0;

  // test the amount of memory needed for fast multiplication algorithm (partial evaluation)
  //
  // firstComponent = index of the first component that has to be precalcualted
  // lastComponent  = index of the last component that has to be precalcualted
  // return value = number of non-zero matrix element
  virtual int PartialFastMultiplicationMemory(int firstComponent, int lastComponent) = 0;

  // enable fast multiplication algorithm
  //
  virtual void EnableFastMultiplication() = 0;

  // enable fast multiplication algorithm (partial evaluation)
  //
  // jobIndex = index of the job that proceeds part of the fast multiplication evaluation
  // nbrJob = number of jobs that proceed the fast multiplication evaluation
  virtual void PartialEnableFastMultiplication(int jobIndex, int nbrJob) = 0;

  // load precalculations from a file
  // 
  // fileName = pointer to a string containg the name of the file where precalculations have to be read
  // return value = true if no error occurs
  virtual bool LoadPrecalculation (char* fileName) = 0;

};

#endif
