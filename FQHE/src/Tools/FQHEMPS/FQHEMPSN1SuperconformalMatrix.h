////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix for the N=1 superconformal states           //
//                                                                            //
//                        last modification : 11/03/2013                      //
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


#ifndef FQHEMPSN1SUPERCONFORMALMATRIX_H
#define FQHEMPSN1SUPERCONFORMALMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"
#include "Vector/LongRationalVector.h"


class LongRationalMatrix;
class BosonOnDiskShort;


class FQHEMPSN1SuperconformalMatrix : public FQHEMPSClustered2RMatrix
{

 protected:


  // effective P level truncation that has to be applied to the descendant of the primary field
  int EffectivePLevel;

 public:
  
  // default constructor 
  //
  FQHEMPSN1SuperconformalMatrix();

  // constructor 
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSN1SuperconformalMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices = 2, bool cylinderFlag = false, double kappa = 1.0);

  // constructor from a file describing the state
  //
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // fileName = name of the file that contains the state description
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSN1SuperconformalMatrix(int pLevel, int nbrBMatrices, char* fileName, bool cylinderFlag = false, double kappa = 1.0);

  // constructor from stored B matrices
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSN1SuperconformalMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag = false, double kappa = 1.0);

  // destructor
  //
  ~FQHEMPSN1SuperconformalMatrix();
  
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
//  virtual SparseRealMatrix ExtractBlock(SparseRealMatrix& matrix, int pLevel1, int q1, int pLevel2, int q2);

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

  // compute the scalar product matrices of the Virasoro descendant
  // 
  // partition = partition that desribes the product of Virasoro generators involved in the scalar product
  // partitionLength = partition length
  // position = position in partition starting from which all the indices are negative
  // centralCharge12 = reference on the value of the central charge divided by 12
  // invCentralCharge3 = reference on the value of three divided by the central charge
  // weight = weight of the primary field that is considered
  // return value = scalar product
  virtual LongRational ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
							       LongRational& centralCharge12, LongRational& invCentralCharge3, LongRational& weight);
  
  // compute the scalar product matrices of the Virasoro descendant, using information from previous levels
  // 
  // partition = partition that desribes the product of Virasoro generators involved in the scalar product
  // partitionLength = partition length
  // position = position in partition starting from which all the indices are negative
  // centralCharge12 = reference on the value of the central charge divided by 12
  // invCentralCharge3 = reference on the value of three divided by the central charge
  // weight = weight of the primary field that is considered
  // precomputedScalarProduct = matrices where scalar product matrix elements computed for previous levels are stored
  // precomputedScalarProductMaxPLevel = maxixum P level that can be accessed through precomputedScalarProduct
  // basis = basis that related the partitions to their index
  // return value = scalar product  
  virtual LongRational ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
							       LongRational& centralCharge12, LongRational& invCentralCharge3, LongRational& weight,
							       LongRationalMatrix* precomputedScalarProduct, int precomputedScalarProductMaxPLevel, 
							       BosonOnDiskShort** basis);

  // compute the matrix elements of any primary field in the Virasoro descendant basis
  // 
  // partition = partition that desribes the product of Virasoro generators involved in the scalar product
  // partitionLength = partition length
  // descendantPosition = location of the primary field
  // position = position in partition starting from which all the indices are negative
  // centralCharge12 = reference on the value of the central charge divided by 12
  // weight1 = weight of the primary field that is considered for the left state
  // weight2 = weight of the primary field that is considered for the right state
  // weight = weight of the primary field whose matrix elements are computed
  // return value = matrix element
  virtual LongRational ComputeDescendantMatrixElement (long* partition, int partitionLength, int descendantPosition, int position, 
						       LongRational& centralCharge12, LongRational& weight1, LongRational& weight2, 
						       LongRational& weight);
  
  // compute the matrix elements of any primary field in the Virasoro descendant basis
  // 
  // partition = partition that desribes the product of Virasoro generators involved in the scalar product
  // partitionLength = partition length
  // descendantPosition = location of the primary field
  // position = position in partition starting from which all the indices are negative
  // centralCharge12 = reference on the value of the central charge divided by 12
  // weight1 = weight of the primary field that is considered for the left state
  // weight2 = weight of the primary field that is considered for the right state
  // weight = weight of the primary field whose matrix elements are computed
  // precomputedDescendantMatrixElement = matrices where matrix elements computed for previous levels are stored
  // precomputedDescendantMatrixElementMaxLeftPLevel = maxixum P level that can be accessed through precomputedDescendantMatrixElement for the left entry
  // precomputedDescendantMatrixElementMaxRightPLevel = maxixum P level that can be accessed through precomputedDescendantMatrixElement for the right entry
  // basis = basis that related the partitions to their index
  // return value = matrix element
  virtual LongRational ComputeDescendantMatrixElement (long* partition, int partitionLength, int descendantPosition, int position, 
						       LongRational& centralCharge12, LongRational& weight1, LongRational& weight2, 
						       LongRational& weight, LongRationalMatrix** precomputedDescendantMatrixElement,
						       int precomputedDescendantMatrixElementMaxLeftPLevel, 
						       int precomputedDescendantMatrixElementMaxRightPLevel, 
						       BosonOnDiskShort** basis);


};

  

#endif
