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


class LongRationalMatrix;
class BosonOnDiskShort;


class FQHEMPSClustered2RMatrix : public FQHEMPSLaughlinMatrix
{

 protected:

  // r index (i.e. clustered (k=2,r) states) 
  int RIndex;


  // indicate the first linearized index for fixed of the total p-level and the U(1) p-level
  int** StartingIndexPerPLevel;
  // number of non-zero vectors for the identity descendants at a given level 
  int* IdentityBasisDimension;
  // number of non-zero vectors for the psi  descendants at a given level 
  int* PsiBasisDimension;
  // number of U(1) mode at a given level 
  int* U1BasisDimension;		  

  // a temporary array to store a partition in the occupation number basis
  unsigned long* TemporaryOccupationNumber;

  // value of the central charge
  LongRational CentralCharge;

  // conformal weight of the identity (or the sigma field for the quasihole sector)
  LongRational WeightIdentity;
  // conformal weight of the psi field (or the phi field for the quasihole sector)
  LongRational WeightPsi;

  // matrix element of <Psi|Psi(1)|0> (or <phi|Psi(1)|sigma>  for the quasihole sector)
  double MatrixElementNormalization;
  // square of matrix element of <Psi|Psi(1)|0> (or <phi|Psi(1)|sigma>  for the quasihole sector)
  LongRational SquareMatrixElementNormalization;
  // conformal weight of the primary field used in the matrix element
  LongRational WeightPrimaryFieldMatrixElement;

  // degeneracy of the transfer matrix largest eigenvalue
  int TransferMatrixDegeneracy;

  // name describing the B matrices 
  char* BMatrixOutputName;

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

  // constructor 
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // cftDirectory = path to the directory where all the pure CFT matrices are stored
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool cylinderFlag = false, double kappa = 1.0);

  // constructor from a file describing the state
  //
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // fileName = name of the file that contains the state description
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSClustered2RMatrix(int pLevel, int nbrBMatrices, char* fileName, bool cylinderFlag = false, double kappa = 1.0);

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
  
  // get the name describing the B matrices 
  // 
  // return value = name 
  virtual char* GetName ();

  // get the filling factor of the state associated the B matrices 
  // 
  // numerator = reference on the filling factor numerator
  // denominator = reference on the filling factor denominator
  virtual void GetFillingFactor(int& numerator, int& denominator);

  // get the degeneracy of the transfer matrix largest eigenvalue
  // 
  // return value = degeneracy 
  virtual int GetTransferMatrixLargestEigenvalueDegeneracy();

  // create the B matrices for the laughlin state
  //
  // cftDirectory = an optional path to the directory where all the CFT matrices are stored
  virtual void CreateBMatrices (char* cftDirectory = 0);

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

  // get the boundary indices of the MPS representation
  //
  // rowIndex = matrix row index
  // columnIndex = matrix column index
  // padding = assume that the state has the estra padding
  virtual void GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding = false);

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
  // weight = weight of the primary field that is considered
  // return value = scalar product
  virtual LongRational ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, LongRational& centralCharge12, LongRational& weight);
  
  // compute the scalar product matrices of the Virasoro descendant, using information from previous levels
  // 
  // partition = partition that desribes the product of Virasoro generators involved in the scalar product
  // partitionLength = partition length
  // position = position in partition starting from which all the indices are negative
  // centralCharge12 = reference on the value of the central charge divided by 12
  // weight = weight of the primary field that is considered
  // precomputedScalarProduct = matrices where scalar product matrix elements computed for previous levels are stored
  // precomputedScalarProductMaxPLevel = maxixum P level that can be accessed through precomputedScalarProduct
  // basis = basis that related the partitions to their index
  // return value = scalar product  
  virtual LongRational ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
							       LongRational& centralCharge12, LongRational& weight,
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
  virtual int Get2RMatrixIndex(int charge, int chargedPartitionIndex, int nbrCharges, int chargeSectorDimension, 
			       int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int globalIndexShift);

  // compute the linearized index of a block of the matrices with fixed charge and p-level values
  //
  // chargedPartitionIndex = index of the partition in the charge sector
  // nbrCharges = total number of partitions at the current level
  // chargeSectorDimension =  total number of partitions in the charge sector
  // fieldIndex = field index (0 for the identity, 1 for psi)
  // neutralPartitionIndex = index of the state in the neutral sector (linearly independant basis)
  // nbrIdentityDescendant = number of linearly independent descendant of the identity at the current level
  // globalIndexShift = index of the first state at the considered level 
  // return value = linearized index
  virtual int Get2RReducedMatrixIndex(int chargedPartitionIndex, int chargeSectorDimension, 
				      int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int globalIndexShift);

};

  
// get the name describing the B matrices 
// 
// return value = name 

inline char* FQHEMPSClustered2RMatrix::GetName ()
{
  return this->BMatrixOutputName;
}

// get the degeneracy of the transfer matrix largest eigenvalue
// 
// return value = degeneracy 

inline int FQHEMPSClustered2RMatrix::GetTransferMatrixLargestEigenvalueDegeneracy()
{
  return this->TransferMatrixDegeneracy;
}

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

// compute the linearized index of a block of the matrices with fixed charge and p-level values
//
// chargedPartitionIndex = index of the partition in the charge sector
// nbrCharges = total number of partitions at the current level
// chargeSectorDimension =  total number of partitions in the charge sector
// fieldIndex = field index (0 for the identity, 1 for psi)
// neutralPartitionIndex = index of the state in the neutral sector (linearly independant basis)
// nbrIdentityDescendant = number of linearly independent descendant of the identity at the current level
// globalIndexShift = index of the first state at the considered level 
// return value = linearized index

inline int FQHEMPSClustered2RMatrix::Get2RReducedMatrixIndex(int chargedPartitionIndex, int chargeSectorDimension, 
							     int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int globalIndexShift)
{
  return (globalIndexShift + (chargedPartitionIndex + chargeSectorDimension * ((fieldIndex * nbrIdentityDescendant) + neutralPartitionIndex)));
}


#endif
