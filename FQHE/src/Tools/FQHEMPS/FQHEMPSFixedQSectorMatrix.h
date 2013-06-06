////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//      class of MPS matrix built from a fixed charge sector of another MPS   //
//                                                                            //
//                        last modification : 19/03/2013                      //
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


#ifndef FQHEMPSFIXEDQSECTORMATRIX_H
#define FQHEMPSFIXEDQSECTORMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"


class FQHEMPSFixedQSectorMatrix : public AbstractFQHEMPSMatrix
{

 protected:

  // pointer to the full MPS matrix without a fixed charge sector
  AbstractFQHEMPSMatrix* MPSMatrix;

  // number of CFT sectors
  int NbrCFTSectors;

  int BMatrixGroupSize;
  // periodicity of the charge sector
  int QPeriodicity;
  // Q sector that has to be selected
  int QSector;
  // CFT sector that has to be selected
  int CFTSector;

  // degeneracy of the transfer matrix largest eigenvalue
  int TransferMatrixLargestEigenvalueDegeneracy;

  // |P| truncation level
  int PLevel;

  // number of N (i.e. charge) values per  truncation level
  int* NbrNValuesPerPLevel;
  // initial N (i.e. charge) value per  truncation level
  int* NInitialValuePerPLevel;
  // last N (i.e. charge) value per  truncation level
  int* NLastValuePerPLevel;
  // first linearized index for each truncation level and charge sector
  int** TotalStartingIndexPerPLevelAndQSector;
  // number of linearized index per truncation level and charge sector
  int** NbrIndicesPerPLevelAndQSector;
  int*** GlobalIndexMapper;

 public:
  
  // default constructor 
  //
  FQHEMPSFixedQSectorMatrix();

  // constructor from two MPS matrices (the number of B matrices has to be identical for all of them)
  //
  // matrix = MPS matrix
  // qPeriodicity = periodicity of the charge sector, if set to zero, guess it from the filling factor
  // qSector = Q sector that has to be selected (from 0 to qPeriodicity-1)
  FQHEMPSFixedQSectorMatrix(AbstractFQHEMPSMatrix* matrix, int qPeriodicity = 0, int qSector = 0);

  // destructor
  //
  ~FQHEMPSFixedQSectorMatrix();
  
  // create the B matrices for the block state
  //
  virtual void CreateBMatrices ();

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

  // get the MPS truncation level
  //
  // return value = truncation level
  virtual int GetTruncationLevel();

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
  // pLevel = tuncation level
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index
  virtual void GetChargeIndexRange (int pLevel, int& minQ, int& maxQ);

  
  // get the boundary indices of the MPS representation
  //
  // rowIndex = matrix row index
  // columnIndex = matrix column index
  // padding = assume that the state has the estra padding
  virtual void GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding = false);

};

// get the filling factor of the state associated the B matrices 
// 
// numerator = reference on the filling factor numerator
// denominator = reference on the filling factor denominator

inline void FQHEMPSFixedQSectorMatrix::GetFillingFactor(int& numerator, int& denominator)
{
  this->MPSMatrix->GetFillingFactor(numerator, denominator);
}

// get the degeneracy of the transfer matrix largest eigenvalue
// 
// return value = degeneracy 

inline int FQHEMPSFixedQSectorMatrix::GetTransferMatrixLargestEigenvalueDegeneracy()
{
  return this->TransferMatrixLargestEigenvalueDegeneracy;;
}

// get the MPS truncation level
//
// return value = truncation level

inline int FQHEMPSFixedQSectorMatrix::GetTruncationLevel()
{
  return this->PLevel;
}

#endif
