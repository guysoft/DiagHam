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

};

// get the array where the matrices are stored
//
// return value = pointer to the array

inline SparseRealMatrix* AbstractFQHEMPSMatrix::GetMatrices()
{
  return this->RealBMatrices;
}

#endif
