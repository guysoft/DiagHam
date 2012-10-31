////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix for the clustered (k=2,r) states            //
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


#ifndef FQHEMPSCLUSTERED2RMATRIX_H
#define FQHEMPSCLUSTERED2RMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSLaughlinMatrix.h"
#include "Vector/LongRationalVector.h"


class FQHEMPSClustered2RMatrix : public FQHEMPSLaughlinMatrix
{

 protected:

  // r index (i.e. clustered (k=2,r) states) 
  int RIndex;

 public:
  
  // default constructor 
  //
  FQHEMPSClustered2RMatrix();

  // constructor 
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices = 2, bool cylinderFlag = false, double kappa = 1.0);

  // constructor from stored B matrices
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag = false, double kappa = 1.0);

  // destructor
  //
  ~FQHEMPSClustered2RMatrix();
  
  // create the B matrices for the laughlin state
  //
  virtual void CreateBMatrices ();

 protected:

  // compute the scalar product matrices of the Virasoro descendant
  // 
  // partition = partition that desribes the product of Virasoro generators involved in the scalar product
  // partitionLength = partition length
  // position = position in partition starting from which all the indices are negative
  // centralCharge12 = reference on the value of the central charge divided by 12
  // weight = weight of the primary field that is considered
  // return value = scalar product
  LongRational ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, LongRational& centralCharge12, LongRational& weight);
  

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
  LongRational ComputeDescendantMatrixElement (long* partition, int partitionLength, int descendantPosition, int position, 
					       LongRational& centralCharge12, LongRational& weight1, LongRational& weight2, 
					       LongRational& weight);

  // compute the linearized index of the B matrix for the (k=2,r) clustered states
  //
  // charge = charge index
  // chargedPartitionIndex =index of the partition in the charge sector
  // nbrCharges = total number of charge indices
  // chargeSectorDimension =  total number of partitions in the charge sector
  // fieldIndex = field index (0 for the identity, 1 for psi)
  // neutralPartitionIndex = index of the state in the neutral sector (linearly independant basis)
  // nbrIdentityDescendant = number of linearly independent descendant of the identity at the current level
  // globalIndexShift = index of the first state at the considered level
  // return value = linearized index
  int Get2RMatrixIndex(int charge, int chargedPartitionIndex, int nbrCharges, int chargeSectorDimension, 
		       int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int globalIndexShift);


};

// compute the linearized index of the B matrix for the (k=2,r) clustered states
//
// charge = charge index
// chargedPartitionIndex =index of the partition in the charge sector
// nbrCharges = total number of charge indices
// chargeSectorDimension =  total number of partitions in the charge sector
// fieldIndex = field index (0 for the identity, 1 for psi)
// neutralPartitionIndex = index of the state in the neutral sector (linearly independant basis)
// nbrIdentityDescendant = number of linearly independent descendant of the identity at the current level
// globalIndexShift = index of the first state at the considered level
// return value = linearized index

inline int FQHEMPSClustered2RMatrix::Get2RMatrixIndex(int charge, int chargedPartitionIndex, 
						      int nbrCharges, int chargeSectorDimension, 
						      int fieldIndex, int neutralPartitionIndex, 
						      int nbrIdentityDescendant, int globalIndexShift)
{
  return (((((nbrIdentityDescendant * fieldIndex) + neutralPartitionIndex) * chargeSectorDimension + chargedPartitionIndex) * nbrCharges + charge) + globalIndexShift);
}

#endif
