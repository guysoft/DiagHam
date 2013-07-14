////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract MPS matrix for the FQHE                 //
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


#ifndef ABSTRACTFQHEMPSMATRIX_H
#define ABSTRACTFQHEMPSMATRIX_H


#include "config.h"
#include "MathTools/Complex.h" 

#include <fstream>

class SparseRealMatrix;
class SparseComplexMatrix;


using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;

class AbstractFQHEMPSMatrix
{

 protected:

  // number of B matrices
  int NbrBMatrices;

  // arrar that describes all the physical indices
  unsigned long* PhysicalIndices;

  // array where the B matrices are stored (for real matrices)
  SparseRealMatrix* RealBMatrices;

  // array where the B matrices are stored (for complex matrices)
  SparseComplexMatrix* ComplexBMatrices;

  // array where the B matrices for quasiholes are stored
  SparseRealMatrix* QuasiholeBMatrices;

 public:
  
  // default constructor 
  //
  AbstractFQHEMPSMatrix();

  // destructor
  //
  ~AbstractFQHEMPSMatrix();
  
  // save the matrices 
  // 
  // fileName = name of the file where the matrices have to be stored
  // return value = true if no error occurred  
  virtual bool SaveMatrices (char* fileName);

  // load the matrices 
  // 
  // fileName = name of the file where the matrices are stored
  // return value = true if no error occurred  
  virtual bool LoadMatrices (char* fileName);

  // get the number of B matrices
  //
  // return value = number of B matrices
  virtual int GetNbrMatrices();

  // get the array where the matrices are stored
  //
  // return value = pointer to the array
  virtual SparseRealMatrix* GetMatrices();

  // get the edge matrix for localized quasiholes, with normal ordering
  //
  // nbrQuasiholes = number of quasiholes
  // quasiholePositions = quasihole positions (for cylinder, positions have to be expressed in perimeter units)
  // return value = pointer to the edge matrix
  virtual SparseComplexMatrix* GetQuasiholeMatrices(int nbrQuasiholes, Complex* quasiholePositions);
  
  // get the name describing the B matrices 
  // 
  // return value = name 
  virtual char* GetName();

  // get the filling factor of the state associated the B matrices 
  // 
  // numerator = reference on the filling factor numerator
  // denominator = reference on the filling factor denominator
  virtual void GetFillingFactor(int& numerator, int& denominator);

  // get the number of particles that fit the root configuration once the number of flux quanta is fixed
  // 
  // nbrFluxQuanta = number of flux quanta
  // padding = assume that the state has the extra padding
  // return value = number of partciles
  virtual int GetMatrixNaturalNbrParticles(int nbrFluxQuanta, bool padding = false);

  // get the degeneracy of the transfer matrix largest eigenvalue
  // 
  // return value = degeneracy 
  virtual int GetTransferMatrixLargestEigenvalueDegeneracy();

  // get the MPS truncation level
  //
  // return value = truncation level
  virtual int GetTruncationLevel();

  // get the number of CFT sectors invloved on the MPS
  //
  // return value = number of CFT sectors
  virtual int GetNbrCFTSectors();

  // get the Q sector shift for a given CFT sector compared to the x=0 CFT sector
  //
  // cftSector = index of the CFT sector
  // return value = Q sector shift
  virtual int GetQValueCFTSectorShift(int cftSector);

  // extract a block with fixed quantum numbers of a given matrix written the MPS basis
  //
  // matrix = reference on the matrix
  // pLevel1 = truncation level of the block left indices
  // q1 = charge index of the block left indices
  // pLevel1 = truncation level of the block right indices
  // q2 = charge index of the block left indices
  // return value = block corresponding to the quantum numbers
  virtual SparseRealMatrix ExtractBlock(SparseRealMatrix& matrix, int pLevel1, int q1, int pLevel2, int q2);

  // extract a block with fixed quantum numbers of a given matrix written the MPS basis
  //
  // matrix = reference on the matrix
  // pLevel1 = truncation level of the block left indices
  // cftSector1 = CFT sector of the blck left indices
  // q1 = charge index of the block left indices
  // pLevel1 = truncation level of the block right indices
  // cftSector2 = CFT sector of the blck right indices
  // q2 = charge index of the block left indices
  // return value = block corresponding to the quantum numbers
  virtual SparseRealMatrix ExtractBlock(SparseRealMatrix& matrix, int pLevel1, int cftSector1, int q1, int pLevel2, int cftSector2, int q2);

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

  // get the charge index range at a given truncation level
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

  // compute the global charge index range at a given truncation level
  // 
  // pLevel = tuncation level
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index

  virtual void ComputeGlobalChargeIndexRange(int pLevel, int& minQ, int& maxQ);

  // compute P, N from the linearized index of the B matrix for the Laughlin states
  //
  // index = linearized index
  // charge = charge index
  // chargedPartitionIndex =index of the partition in the charge sector
  virtual void GetPNFromMatrixIndex(int index, int& charge, int& chargedPartitionIndex);

  // compute the level and the charge index of a given matrix index
  //
  // index = matrix index
  // pLevel = reference on the level
  // qValue = reference on the charge index
  virtual void GetChargeAndPLevelFromMatrixIndex(int index, int& pLevel, int& qValue);

  // get the boundary indices of the MPS representation
  //
  // rowIndex = matrix row index
  // columnIndex = matrix column index
  // padding = assume that the state has the extra padding
  virtual void GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding = false);

  // get the array of physical indices
  //
  // return value  = array of physical indices
  virtual unsigned long* GetPhysicalIndices();

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

};

// get the number of B matrices
//
// return value = number of B matrices

inline int AbstractFQHEMPSMatrix::GetNbrMatrices()
{
  return this->NbrBMatrices;
}

// get the array where the matrices are stored
//
// return value = pointer to the array

inline SparseRealMatrix* AbstractFQHEMPSMatrix::GetMatrices()
{
  return this->RealBMatrices;
}

// compute P, N from the linearized index of the B matrix for the Laughlin states
//
// index = linearized index
// charge = charge index
// chargedPartitionIndex =index of the partition in the charge sector

inline void AbstractFQHEMPSMatrix::GetPNFromMatrixIndex(int index, int& charge, int& chargedPartitionIndex)
{
  cout << "Dummy. " << endl;
}

// compute the level and the charge index of a given matrix index
//
// index = matrix index
// pLevel = reference on the level
// qValue = reference on the charge index

inline void AbstractFQHEMPSMatrix::GetChargeAndPLevelFromMatrixIndex(int index, int& pLevel, int& qValue)
{
  pLevel = -1;
  qValue = -1;
}

// get the degeneracy of the transfer matrix largest eigenvalue
// 
// return value = degeneracy 

inline int AbstractFQHEMPSMatrix::GetTransferMatrixLargestEigenvalueDegeneracy()
{
  return 1;
}

// get the MPS truncation level
//
// return value = truncation level

inline int AbstractFQHEMPSMatrix::GetTruncationLevel()
{
  return 0;
}

// get the number of CFT sectors invloved on the MPS
//
// return value = number of CFT sectors

inline int AbstractFQHEMPSMatrix::GetNbrCFTSectors()
{
  return 1;
}

// get the array of physical indices
//
// return value  = array of physical indices

inline unsigned long* AbstractFQHEMPSMatrix::GetPhysicalIndices()
{
  return this->PhysicalIndices;
}

#endif
