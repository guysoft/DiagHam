////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of MPS matrix for the Laughlin state                //
//                                                                            //
//                        last modification : 30/10/2012                      //
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


#ifndef FQHEMPSLAUGHLINMATRIX_H
#define FQHEMPSLAUGHLINMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"
#include "MathTools/FactorialCoefficient.h"


class FQHEMPSLaughlinMatrix : public AbstractFQHEMPSMatrix
{

 protected:

  // index of the Laughlin (i.e. 1/nu) 
  int LaughlinIndex;

  // |P| level truncation
  int PLevel;

  // true if B_0 has to be normalized on the cylinder geometry
  bool CylinderFlag;
  // cylinder aspect ratio
  bool Kappa;

 public:
  
  // default constructor 
  //
  FQHEMPSLaughlinMatrix();

  // constructor 
  //
  // laughlinIndex = power of the Laughlin part (i.e. 1/nu)
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital)
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSLaughlinMatrix(int laughlinIndex, int pLevel, int nbrBMatrices = 1, bool cylinderFlag = false, double kappa = 1.0);

  // destructor
  //
  ~FQHEMPSLaughlinMatrix();
  
  // create the B matrices for the laughlin state
  //
  virtual void CreateBMatrices ();

 protected:

  // create the matrix element of the B matrix U(1) part
  //
  // chargeNumerator = numerator of the charge (in sqrt(q) unit)
  // chargeDenominator = denominator of the charge (in sqrt(q) unit)
  // partition1 = U(1) partition associated to the left state
  // p1Level = length of partition1
  // partition2 = U(1) partition associated to the left state
  // p1Level = length of partition2
  // coef = reference on a temporary factorial coefficient
  // return value = matrix element
  double CreateLaughlinAMatrixElement (int chargeNumerator, int chargeDenominator, 
				       unsigned long* partition1, unsigned long* partition2, 
				       int p1Level, int p2Level, FactorialCoefficient& coef);

};

#endif
