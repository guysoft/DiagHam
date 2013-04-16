////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of MPS matrix for the Read-Rezayi k=3 state             //
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


#ifndef FQHEMPSREADREZAYI3MATRIX_H
#define FQHEMPSREADREZAYI3MATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"


class FQHEMPSReadRezayi3Matrix : public FQHEMPSClustered2RMatrix
{

 protected:


 public:
  
  // default constructor 
  //
  FQHEMPSReadRezayi3Matrix();

  // constructor 
  //
  // laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic RR state)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital)
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSReadRezayi3Matrix(int laughlinIndex, int pLevel, int nbrBMatrices = 1, bool cylinderFlag = false, double kappa = 1.0);

  // constructor from stored B matrices
  //
  // laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic RR state)  
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSReadRezayi3Matrix(int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag = false, double kappa = 1.0);

  // destructor
  //
  ~FQHEMPSReadRezayi3Matrix();
  
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
  virtual void CreateBMatrices ();

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

  // compute the linearized index of the B matrix for the (k=2,r) clustered states
  //
  // charge = charge index
  // chargedPartitionIndex =index of the partition in the charge sector
  // nbrCharges = total number of charge indices
  // chargeSectorDimension =  total number of partitions in the charge sector
  // fieldIndex = field index (0 for the identity, 1 for psi_{+1}, 3 for psi_{-1}, 5 for W_-3)
  // neutralPartitionIndex = index of the state in the neutral sector (linearly independant basis)
  // nbrIdentityDescendant = number of linearly independent descendants of the identity at the current level
  // nbrPsiDescendant = number of linearly independent descendants of the psi_{+1} at the current level
  // globalIndexShift = index of the first state at the considered level
  // return value = linearized index
  int GetReadRezayiK3MatrixIndex(int charge, int chargedPartitionIndex, int nbrCharges, int chargeSectorDimension, 
				 int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int nbrPsiDescendant, 
				 int globalIndexShift);


};


// get the degeneracy of the transfer matrix largest eigenvalue
// 
// return value = degeneracy 

inline int FQHEMPSReadRezayi3Matrix::GetTransferMatrixLargestEigenvalueDegeneracy()
{
  return 5;
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
// nbrPsiDescendant = number of linearly independent descendants of the psi_{+1} at the current level
// globalIndexShift = index of the first state at the considered level
// return value = linearized index

inline int FQHEMPSReadRezayi3Matrix::GetReadRezayiK3MatrixIndex(int charge, int chargedPartitionIndex, 
								int nbrCharges, int chargeSectorDimension, 
								int fieldIndex, int neutralPartitionIndex, 
								int nbrIdentityDescendant, int nbrPsiDescendant, 
								int globalIndexShift)
{
  return ((((((nbrIdentityDescendant * (fieldIndex & 1)) + (nbrPsiDescendant * (fieldIndex >> 1))) + neutralPartitionIndex) * chargeSectorDimension + chargedPartitionIndex) * nbrCharges + charge) + globalIndexShift);
}

#endif
