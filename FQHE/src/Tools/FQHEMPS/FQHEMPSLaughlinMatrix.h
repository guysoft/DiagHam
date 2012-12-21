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
  
  // number of charge indices
  int NbrNValue;
  // first linearized index for each truncation level
  int* TotalStartingIndexPerPLevel;
  // number of linearized index per truncation level
  int* NbrIndicesPerPLevel;

  // true if B_0 has to be normalized on the cylinder geometry
  bool CylinderFlag;
  // cylinder aspect ratio
  double Kappa;

 public:
  
  // default constructor 
  //
  FQHEMPSLaughlinMatrix();

  // constructor 
  //
  // laughlinIndex = power of the Laughlin part (i.e. 1/nu)
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSLaughlinMatrix(int laughlinIndex, int pLevel, int nbrBMatrices = 2, bool cylinderFlag = false, double kappa = 1.0);

  // constructor from stored B matrices
  //
  // laughlinIndex = power of the Laughlin part (i.e. 1/nu)
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSLaughlinMatrix(int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag = false, double kappa = 1.0);

  // destructor
  //
  ~FQHEMPSLaughlinMatrix();
  
  // create the B matrices for the laughlin state
  //
  virtual void CreateBMatrices ();

  // extract a block with fixed quantum numbers of a given matrix written the MPS basis
  //
  // matrix = reference on the matrix
  // pLevel1 = tuncation level of the block left indices
  // q1 = charge index of the block left indices
  // pLevel1 = tuncation level of the block right indices
  // q2 = charge index of the block left indices
  // return value = block corresponding to the quantum numbers
  virtual SparseRealMatrix ExtractBlock(SparseRealMatrix& matrix, int pLevel1, int q1, int pLevel2, int q2);

  // get the range for the bond index when fixing the tuncation level and the charge index
  //
  // pLevel = tuncation level of the block
  // qValue = charge index of the block
  // return value = range for the bond index with fixed tuncation level and charge index
  virtual int GetBondIndexRange(int pLevel, int qValue);

  // get the bond index for a fixed truncation level and the charge index 
  //
  // localIndex = bond index in the pLevel and qValue restricted range
  // pLevel = tuncation level of the block
  // qValue = charge index of the block
  // return value = bond index in the full bond index range
  virtual int GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue);

  // get the charge index range
  // 
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index
  virtual void GetChargeIndexRange (int& minQ, int& maxQ);

 protected:

  // load the specific informations from the file header
  // 
  // file = reference on the input file stream
  // return value = true if no error occurred  
  virtual bool LoadHeader (ifstream& file);

  // save the specific informations to the file header 
  // 
  // file = reference on the output file stream
  // return value = true if no error occurred  
  virtual bool SaveHeader (ofstream& file);

  // compute the linearized index of the B matrix for the Laughlin states
  //
  // charge = charge index
  // chargedPartitionIndex =index of the partition in the charge sector
  // nbrCharges = total number of charge indices
  // globalIndexShift = index of the first state at the considered level
  // return value = linearized index
  int GetMatrixIndex(int charge, int chargedPartitionIndex, int nbrCharges, int globalIndexShift);

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

// compute the linearized index of the B matrix for the Laughlin states
//
// charge = charge index
// chargedPartitionIndex =index of the partition in the charge sector
// nbrCharges = total number of charge indices
// globalIndexShift = index of the first state at the considered level
// return value = linearized index

inline int FQHEMPSLaughlinMatrix::GetMatrixIndex(int charge, int chargedPartitionIndex, 
						 int nbrCharges, int globalIndexShift)
{
  return ((chargedPartitionIndex * nbrCharges + charge) + globalIndexShift);
}

#endif
