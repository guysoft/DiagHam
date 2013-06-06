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
class RealMatrix;
class RealSymmetricMatrix;
class BosonOnDiskShort;
class AbstractArchitecture;

class FQHEMPSClustered2RMatrix : public FQHEMPSLaughlinMatrix
{


  friend class FQHEMPSEvaluateCFTOperation;

 protected:

  // r index (i.e. clustered (k=2,r) states) 
  int RIndex;

  // number of CFT sectors
  int NbrCFTSectors;

  // indicate the first linearized index for fixed of the total p-level and the U(1) p-level
  int** StartingIndexPerPLevel;
  // number of non-zero vectors for the identity descendants at a given level 
  int* IdentityBasisDimension;
  // number of non-zero vectors for the psi  descendants at a given level 
  int* PsiBasisDimension;
  // number of U(1) mode at a given level 
  int* U1BasisDimension;		  

  // number of N (i.e. charge) values per  truncation level
  int** NbrNValuesPerPLevelCFTSector;
  // initial N (i.e. charge) value per  truncation level
  int** NInitialValuePerPLevelCFTSector;
  // last N (i.e. charge) value per  truncation level
  int** NLastValuePerPLevelCFTSector;

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

  // use arbitrary precision numbers for all the CFT calculations
  bool UseRationalFlag;
  
  // first linearized index for each truncation level, CFT sector, Q sector
  int*** StartingIndexPerPLevelCFTSectorQValue;
  // number of linearized indices for each truncation level, CFT sector, Q sector
  int*** NbrIndexPerPLevelCFTSectorQValue;
  // first linearized index for each truncation level, CFT sector, Q sector and U(1) sector
  int**** StartingIndexPerPLevelCFTSectorQValueU1Sector;
  // number of linearized indices for each truncation level, CFT sector, Q sector and U(1) sector
  int**** NbrIndexPerPLevelCFTSectorQValueU1Sector;
  // dimensions of each CFT sector at a given level
  int** NeutralSectorDimension;

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
  // useRational = use arbitrary precision numbers for all the CFT calculations
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // architecture = architecture to use for precalculation
  FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices = 2, bool useRational = true, 
			   bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, AbstractArchitecture* architecture = 0);

  // constructor 
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // cftDirectory = path to the directory where all the pure CFT matrices are stored
  // useRational = use arbitrary precision numbers for all the CFT calculations
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // architecture = architecture to use for precalculation
  FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool useRational = true,
			   bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, AbstractArchitecture* architecture = 0);

  // constructor from a file describing the state
  //
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // fileName = name of the file that contains the state description
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // architecture = architecture to use for precalculation
  FQHEMPSClustered2RMatrix(int pLevel, int nbrBMatrices, char* fileName, bool trimChargeIndices = false, bool cylinderFlag = false, 
			   double kappa = 1.0, AbstractArchitecture* architecture = 0);

  // constructor from stored B matrices
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSClustered2RMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0);

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
  // architecture = architecture to use for precalculation
  virtual void CreateBMatrices (char* cftDirectory, AbstractArchitecture* architecture);

  // get the number of CFT sectors invloved on the MPS
  //
  // return value = number of CFT sectors
  virtual int GetNbrCFTSectors();

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

  // get the range for the bond index when fixing the tuncation level, charge and CFT sector index
  //
  // pLevel = tuncation level of the block
  // qValue = charge index of the block
  // cftSector = CFT sector index of the block
  // return value = range for the bond index with fixed tuncation level and charge index
  virtual int GetBondIndexRange(int pLevel, int qValue, int cftSector);

  // get the bond index for a fixed truncation level and the charge index 
  //
  // localIndex = bond index in the pLevel and qValue restricted range
  // pLevel = tuncation level of the block
  // qValue = charge index of the block
  // return value = bond index in the full bond index range
  virtual int GetBondIndexWithFixedChargeAndPLevel(int localIndex, int pLevel, int qValue);

  // get the bond index for a fixed truncation level, charge and CFT sector index
  //
  // localIndex = bond index in the pLevel and qValue restricted range
  // pLevel = tuncation level of the block
  // qValue = charge index of the block
  // cftSector = CFT sector index of the block
  // return value = bond index in the full bond index range
  virtual int GetBondIndexWithFixedChargePLevelCFTSector(int localIndex, int pLevel, int qValue, int cftSector);

  // get the charge index range
  // 
  // pLevel = tuncation level
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index
  virtual void GetChargeIndexRange (int pLevel, int& minQ, int& maxQ);

  // get the charge index range at a given truncation level and in a given CFT sector
  // 
  // pLevel = tuncation level
  // cftSector = CFT sector
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index
  virtual void GetChargeIndexRange (int pLevel, int cftSector, int& minQ, int& maxQ);

  // get the boundary indices of the MPS representation
  //
  // rowIndex = matrix row index
  // columnIndex = matrix column index
  // padding = assume that the state has the estra padding
  virtual void GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding = false);

  // compute the charge index range at a given truncation level
  // 
  // pLevel = tuncation level
  // cftSector = CFT sector
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index
  virtual void ComputeChargeIndexRange(int pLevel, int cftSector, int& minQ, int& maxQ);

  // get the number of particles that fit the root configuration once the number of flux quanta is fixed
  // 
  // nbrFluxQuanta = number of flux quanta
  // padding = assume that the state has the extra padding
  // return value = number of partciles
  virtual int GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding = false);

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
  // temporaryOccupationNumber = local temporary to store the occupation numbers 
  // return value = scalar product  
  virtual LongRational ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
							       LongRational& centralCharge12, LongRational& weight,
							       LongRationalMatrix* precomputedScalarProduct, int precomputedScalarProductMaxPLevel, 
							       BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber);

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
  // temporaryOccupationNumber = local temporary to store the occupation numbers 
  // return value = scalar product
  virtual double ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, 
							 double& centralCharge12, double& weight,
							 RealSymmetricMatrix* precomputedScalarProduct, int precomputedScalarProductMaxPLevel, 
							 BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber);

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
  // temporaryOccupationNumber = local temporary to store the occupation numbers 
  // return value = matrix element
  virtual LongRational ComputeDescendantMatrixElement (long* partition, int partitionLength, int descendantPosition, int position, 
						       LongRational& centralCharge12, LongRational& weight1, LongRational& weight2, 
						       LongRational& weight, LongRationalMatrix** precomputedDescendantMatrixElement,
						       int precomputedDescendantMatrixElementMaxLeftPLevel, 
						       int precomputedDescendantMatrixElementMaxRightPLevel, 
						       BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber);


  // compute the matrix elements of any primary field in the Virasoro descendant basis, using double numbers instead of long rational
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
  // temporaryOccupationNumber = local temporary to store the occupation numbers 
  // return value = matrix element  
  virtual double ComputeDescendantMatrixElement (long* partition, int partitionLength, 
						 int descendantPosition, int position,
						 double& centralCharge12, double& weight1, 
						 double& weight2, double& weight,
						 RealMatrix** precomputedDescendantMatrixElement, 
						 int precomputedDescendantMatrixElementMaxLeftPLevel, 
						 int precomputedDescendantMatrixElementMaxRightPLevel, 
						 BosonOnDiskShort** basis, unsigned long* temporaryOccupationNumber);

  // compute the various arrays required to convert from quantum numbers and local indices to a global linearized index
  //
  // return value = dimension of the B matrix
  virtual int ComputeLinearizedIndexArrays();

  // compute the linearized index of the B matrix for the (k=2,r) clustered states
  //
  // pLevel = current total level
  // cftSector = CFT sector
  // qValue = charge index (i.e. Q)
  // chargeSectorLevel = level for the charge sector
  // chargeSectorIndex = index of the charge sector
  // cftSectorIndex = index within the current CFT sector and corresponding level (pLevel - chargeSectorLevel)
  // return value = linearized index
  virtual int Get2RMatrixIndexV2(int pLevel, int cftSector, int qValue, 
				 int chargeSectorLevel, int chargeSectorIndex, int cftSectorIndex);

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
// pLevel = current total level
// cftSector = CFT sector
// qValue = charge index (i.e. Q)
// chargeSectorLevel = level for the charge sector
// chargeSectorIndex = index of the charge sector
// cftSectorIndex = index within the current CFT sector and corresponding level (pLevel - chargeSectorLevel)
// return value = linearized index

inline int FQHEMPSClustered2RMatrix::Get2RMatrixIndexV2(int pLevel, int cftSector, int qValue, 
							int chargeSectorLevel, int chargeSectorIndex, int cftSectorIndex)
{
  return (this->StartingIndexPerPLevelCFTSectorQValueU1Sector[pLevel][cftSector][qValue - this->NInitialValuePerPLevelCFTSector[pLevel][cftSector]][chargeSectorLevel] 
	  + (this->NeutralSectorDimension[cftSector][pLevel - chargeSectorLevel] * chargeSectorIndex) + cftSectorIndex);
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

// get the number of CFT sectors invloved on the MPS
//
// return value = number of CFT sectors

inline int  FQHEMPSClustered2RMatrix::GetNbrCFTSectors()
{
  return this->NbrCFTSectors;
}


#endif
