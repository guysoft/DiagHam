////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of qhe on sphere main task                     //
//                                                                            //
//                        last modification : 10/06/2004                      //
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


#ifndef QHEONSPHEREMAINTASK_H
#define QHEONSPHEREMAINTASK_H


#include "config.h"

#include "MainTask/AbstractMainTask.h"


class AbstractQHEHamiltonian;
class AbstractHilbertSpace;
class OptionManager;


class QHEOnSphereMainTask: public AbstractMainTask
{

 protected:

  // name of the file where results have to be stored
  char* OutputFileName;

  // twice the total momentum value of the system
  int LValue;

  // energy shift applied to the hamiltonian
  double EnergyShift;

  // pointer to the current Hamiltonian
  AbstractQHEHamiltonian* Hamiltonian;
  // pointer to the current Hilbert space
  AbstractHilbertSpace* Space;

  // name of the file where hamiltonian precalculations have to be saved (null if no precalculation has to be saved)
  char* SavePrecalculationFileName;
  // maximum Hilbert space dimension for which full diagonalization is applied
  int FullDiagonalizationLimit;
  // enable Lanczos disk resume capabilities
  bool DiskFlag;
  // resume from disk datas
  bool ResumeFlag;
  // number of eigenvalues to evaluate 
  int NbrEigenvalue;
  // number of lanczos iteration (for the current run)
  int NbrIterLanczos;
  // maximum number of Lanczos iteration
  int MaxNbrIterLanczos;
  // maximum number of vector in RAM during Lanczos iteration
  int VectorMemory;

 public:

  // constructor
  //  
  // options = pointer to the options managers containing all running options
  // space = pointer to the current Hilbert space
  // hamiltonian = pointer to the current Hamiltonian
  // lValue = twice the total momentum value of the system
  // shift = energy shift that is applied to the hamiltonian
  // outputFileName = name of the file where results have to be stored
  QHEOnSphereMainTask(OptionManager* options, AbstractHilbertSpace* space, 
		      AbstractQHEHamiltonian* hamiltonian, int lValue, double shift, char* outputFileName);
  
  // destructor
  //  
  ~QHEOnSphereMainTask();
  
  // execute the main task
  // 
  // return value = 0 if no error occurs, else return error code
  int ExecuteMainTask();

};

#endif
