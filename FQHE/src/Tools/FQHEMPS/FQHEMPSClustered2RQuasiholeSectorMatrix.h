////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of MPS matrix for the clustered (k=2,r) states            //
//                         in their quasihole sector                          //
//                                                                            //
//                        last modification : 11/02/2013                      //
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


#ifndef FQHEMPSCLUSTERED2RQUASIHOLESECTORMATRIX_H
#define FQHEMPSCLUSTERED2RQUASIHOLESECTORMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSClustered2RMatrix.h"
#include "Vector/LongRationalVector.h"


class LongRationalMatrix;
class BosonOnDiskShort;


class FQHEMPSClustered2RQuasiholeSectorMatrix : public FQHEMPSClustered2RMatrix
{

 protected:

  // conformal weights of the sigma field
  LongRational WeightSigma;
  // conformal weights of the phi field (if r is even, then the sigma field has to be identical to the phi field)
  LongRational WeightPhi;

  // flag that indicates that the field is self conjugate
  bool SelfDualFlag;

  // matrix element of <phi|Psi(1)|sigma>
//  double MatrixElementNormalization;
  // conformal weight of the primary field used in the matrix element
//  LongRational WeightPrimaryFieldMatrixElement;

 public:
  
  // default constructor 
  //
  FQHEMPSClustered2RQuasiholeSectorMatrix();

  // constructor 
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // architecture = architecture to use for precalculation
  FQHEMPSClustered2RQuasiholeSectorMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices = 2, bool cylinderFlag = false, double kappa = 1.0, 
					  AbstractArchitecture* architecture = 0);

  // constructor 
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // cftDirectory = path to the directory where all the pure CFT matrices are stored
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // architecture = architecture to use for precalculation
  FQHEMPSClustered2RQuasiholeSectorMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool cylinderFlag = false, double kappa = 1.0, 
					  AbstractArchitecture* architecture = 0);

  // constructor from a file describing the state
  //
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
  // fileName = name of the file that contains the state description
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // architecture = architecture to use for precalculation
  FQHEMPSClustered2RQuasiholeSectorMatrix(int pLevel, int nbrBMatrices, char* fileName, bool cylinderFlag = false, double kappa = 1.0, 
					  AbstractArchitecture* architecture = 0);

  // constructor from stored B matrices
  //
  // rindex = r index (i.e. clustered (k=2,r) states) 
  // laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)  
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  FQHEMPSClustered2RQuasiholeSectorMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag = false, double kappa = 1.0);

  // destructor
  //
  ~FQHEMPSClustered2RQuasiholeSectorMatrix();
  
  // create the B matrices for the laughlin state
  //
  // cftDirectory = an optional path to the directory where all the CFT matrices are stored
  // architecture = architecture to use for precalculation
  virtual void CreateBMatrices (char* cftDirectory = 0, AbstractArchitecture* architecture = 0);

 protected:

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

inline int FQHEMPSClustered2RQuasiholeSectorMatrix::Get2RMatrixIndex(int charge, int chargedPartitionIndex, 
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

inline int FQHEMPSClustered2RQuasiholeSectorMatrix::Get2RReducedMatrixIndex(int chargedPartitionIndex, int chargeSectorDimension, 
									    int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int globalIndexShift)
{
  return (globalIndexShift + (chargedPartitionIndex + chargeSectorDimension * ((fieldIndex * nbrIdentityDescendant) + neutralPartitionIndex)));
}


#endif
