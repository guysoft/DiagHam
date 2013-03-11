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

  // get the B matrices corresponding to localized quasiholes
  //
  // nbrQuasiholes = number of quasiholes
  // quasiholePositions = quasihole positions
  // return value = array of nbrQuasiholes matrices corresponding to each quasihole
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

  // get the degeneracy of the transfer matrix largest eigenvalue
  // 
  // return value = degeneracy 
  virtual int GetTransferMatrixLargestEigenvalueDegeneracy();

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

#endif
