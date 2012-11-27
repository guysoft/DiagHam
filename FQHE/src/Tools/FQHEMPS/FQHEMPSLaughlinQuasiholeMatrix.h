////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of MPS matrix for the Laughlin state including quasiholes      //
//                                                                            //
//                        last modification : 17/11/2012                      //
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


#ifndef FQHEMPSLAUGHLINQUASIHOLEMATRIX_H
#define FQHEMPSLAUGHLINQUASIHOLEMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"
#include "MathTools/FactorialCoefficient.h"
#include "Tools/FQHEMPS/FQHEMPSLaughlinMatrix.h"


class FQHEMPSLaughlinQuasiholeMatrix : public FQHEMPSLaughlinMatrix
{


 protected:


 public:
  
  // default constructor 
  //
  FQHEMPSLaughlinQuasiholeMatrix();

  // constructor 
  //
  // laughlinIndex = power of the Laughlin part (i.e. 1/nu)
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSLaughlinQuasiholeMatrix(int laughlinIndex, int pLevel, int nbrBMatrices = 2, bool cylinderFlag = false, double kappa = 1.0);

  // constructor from stored B matrices
  //
  // laughlinIndex = power of the Laughlin part (i.e. 1/nu)
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSLaughlinQuasiholeMatrix(int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag = false, double kappa = 1.0);

  // destructor
  //
  ~FQHEMPSLaughlinQuasiholeMatrix();
  
  // get the B matrices corresponding to localized quasiholes
  //
  // nbrQuasiholes = number of quasiholes
  // quasiholePositions = quasihole positions
  // return value = array of nbrQuasiholes matrices corresponding to each quasihole
  virtual SparseComplexMatrix* GetQuasiholeMatrices(int nbrQuasiholes, Complex* quasiholePositions);

 protected:

  // create the B matrices for the laughlin state
  //
  virtual void CreateBMatrices ();

};

#endif
