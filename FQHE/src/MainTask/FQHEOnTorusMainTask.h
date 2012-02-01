////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                        class of qhe on disk main task                      //
//                                                                            //
//                        last modification : 08/10/2009                      //
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


#ifndef FQHEONTORUSMAINTASK_H
#define FQHEONTORUSMAINTASK_H


#include "config.h"

#include "MainTask/QHEOnDiskMainTask.h"

#include <iostream>
using std::ofstream;

class AbstractQHEHamiltonian;
class AbstractHilbertSpace;
class OptionManager;


class FQHEOnTorusMainTask: public QHEOnDiskMainTask
{

 protected:

  // momentum along the x-axis
  int KxValue;
  // momentum along the y-axis
  int KyValue;
  // flag if only the Ky quantum number is fixed   
  bool KyOnlyFlag;

  // flag indicating whether multiplicity of state is known
  bool MultiplicityFlag;
  // value of multiplicity
  int Multiplicity;

 public:

  // constructor
  //  
  // options = pointer to the options managers containing all running options
  // space = pointer to the current Hilbert space
  // hamiltonian = pointer to the current Hamiltonian
  // kyValue = total momentum value of the system along the y-axis
  // shift = energy shift that is applied to the hamiltonian
  // outputFileName = name of the file where results have to be stored
  // firstRun = flag that indicates if it the first time the main task is used
  // eigenvectorFileName = prefix to add to the name of each file that will contain an eigenvector
  FQHEOnTorusMainTask(OptionManager* options, AbstractHilbertSpace* space, 
		      AbstractQHEHamiltonian* hamiltonian, int kyValue, double shift, char* outputFileName,
		      bool firstRun = true, char* eigenvectorFileName = 0);
  
  // destructor
  //  
  ~FQHEOnTorusMainTask();
  
  // execute the main task
  // 
  // return value = 0 if no error occurs, else return error code
  int ExecuteMainTask();

  // set a kx-value
  //
  // kxValue = kx value
  void SetKxValue(int kxValue);

  // set multiplicity of a given momentum sector
  //
  // multiplicity = sector multiplicity
  void SetMultiplicity(int multiplicity);

 protected:
  
  // write a line of output to the results file
  //
  // file = stream to write to
  // value = numerical value to be printed after columns for flux and momentum (if defined)
  // terminate = indicate if line should be terminated with endl
  void WriteResult(ofstream& file, double value, bool terminate=true);

  // do the Hamiltonian diagonalization in a given Hilbert subspace
  //
  // subspaceDescription = name of the file that contains the vector files used to describe the Hilbert subspace
  // file = reference on the output file stream where eigenvalues have to be stored
  void DiagonalizeInHilbertSubspace(char* subspaceDescription, ofstream& file);

};

#endif
