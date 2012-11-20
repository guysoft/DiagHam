////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of MPS matrix for the Laughlin state including quasiholes      //
//                                                                            //
//                        last modification : 17/11/2012                      //
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


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSLaughlinQuasiholeMatrix.h"
#include "Matrix/SparseComplexMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"

#include "GeneralTools/Endian.h"

#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor 
//

FQHEMPSLaughlinQuasiholeMatrix::FQHEMPSLaughlinQuasiholeMatrix()
{
}

// constructor 
//
// laughlinIndex = power of the Laughlin part (i.e. 1/nu)
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSLaughlinQuasiholeMatrix::FQHEMPSLaughlinQuasiholeMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, 
							       bool cylinderFlag, double kappa)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->QuasiholeBMatrices = new SparseRealMatrix [1];
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->CreateBMatrices();
}

// constructor from stored B matrices
//
// laughlinIndex = power of the Laughlin part (i.e. 1/nu)
// pLevel = |P| level truncation
// fileName = name of the file that contains the B matrices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSLaughlinQuasiholeMatrix::FQHEMPSLaughlinQuasiholeMatrix(int laughlinIndex, int pLevel, char* fileName, 
							       bool cylinderFlag, double kappa)
{
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->LoadMatrices(fileName);
}

// destructor
//

FQHEMPSLaughlinQuasiholeMatrix::~FQHEMPSLaughlinQuasiholeMatrix()
{
  delete[] this->TotalStartingIndexPerPLevel;
  delete[] this->NbrIndicesPerPLevel;
}

// get the B matrices corresponding to localized quasiholes
//
// nbrQuasiholes = number of quasiholes
// quasiholePositions = quasihole positions
// return value = array of nbrQuasiholes matrices corresponding to each quasihole

SparseComplexMatrix* FQHEMPSLaughlinQuasiholeMatrix::GetQuasiholeBMatrices(int nbrQuasiholes, Complex* quasiholePositions)
{
  SparseComplexMatrix* TmpQuasiholeBMatrices = new SparseComplexMatrix[nbrQuasiholes];
  for (int i = 0; i < nbrQuasiholes; ++i)
    TmpQuasiholeBMatrices[i] = SparseComplexMatrix(this->QuasiholeBMatrices[0]);
  return TmpQuasiholeBMatrices;
}

  
// create the B matrices for the laughlin state
//

void FQHEMPSLaughlinQuasiholeMatrix::CreateBMatrices ()
{
  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  SparseRealMatrix* BMatrices = new SparseRealMatrix[this->NbrBMatrices];
  SparseComplexMatrix* BQuasiholeMatrices = new SparseComplexMatrix[1];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort(i, i, this->PLevel + 1);
    }
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel[0] = 0;
  this->NbrNValue = ((2 * this->PLevel) + this->LaughlinIndex);
  int NValueShift = this->NbrNValue - 1;
  this->NbrIndicesPerPLevel[0] = U1BosonBasis[0]->GetHilbertSpaceDimension() * this->NbrNValue;
  for (int i = 1; i <= this->PLevel; ++i)
    {
      this->TotalStartingIndexPerPLevel[i] = this->TotalStartingIndexPerPLevel[i - 1] + this->NbrIndicesPerPLevel[i - 1];
      this->NbrIndicesPerPLevel[i] = U1BosonBasis[i]->GetHilbertSpaceDimension()  * this->NbrNValue;
    }
  int MatrixSize = this->NbrIndicesPerPLevel[this->PLevel] + this->TotalStartingIndexPerPLevel[this->PLevel];
  cout << "B matrix size = " << MatrixSize << endl;
  BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace = U1BosonBasis[i];
      for (int j = 1; j < this->NbrNValue; ++j)
	{
	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
	    {
	      double Tmp = 1.0;
              if (this->CylinderFlag)
		Tmp *= exp(-this->Kappa * this->Kappa * (((double) i)
							 + ((j - 1.0 - 0.5 * NValueShift) * (j - 1.0 - 0.5 * NValueShift) / (4.0 * this->LaughlinIndex))
							 + (((j - 0.5 * NValueShift) * (j - 0.5 * NValueShift)) / (4.0 * this->LaughlinIndex))));
	      BMatrices[0].SetMatrixElement(this->GetMatrixIndex(j - 1, k, this->NbrNValue, this->TotalStartingIndexPerPLevel[i]),
					    this->GetMatrixIndex(j, k, this->NbrNValue, this->TotalStartingIndexPerPLevel[i]), Tmp);
	    }
	}
    }
   
  for (int m = 1; m < this->NbrBMatrices; ++m)
    BMatrices[m] = SparseRealMatrix(MatrixSize, MatrixSize);

  unsigned long* Partition1 = new unsigned long [this->PLevel + 2];
  unsigned long* Partition2 = new unsigned long [this->PLevel + 2];
  FactorialCoefficient Coef;

  for (int i = 0; i <= this->PLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace1 = U1BosonBasis[i];
      for (int j = 0; j <= this->PLevel; ++j)
	{
	  BosonOnDiskShort* TmpSpace2 = U1BosonBasis[j];
	  int N2 = (2 * (j - i) - this->LaughlinIndex + 1 + NValueShift) / 2;
	  int N1 = N2 + (this->LaughlinIndex - 1);
	  for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
	    {
	      TmpSpace1->GetOccupationNumber(k1, Partition1);
	      for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
		{
		  TmpSpace2->GetOccupationNumber(k2, Partition2);
		  for (int m = 1; m < this->NbrBMatrices; ++m)
		    {
		      double Tmp = this->CreateLaughlinAMatrixElement(this->LaughlinIndex * m * m, 1, Partition1, Partition2, i, j, Coef);
		      if (this->CylinderFlag)
			Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (((double) (i + j))
								       + ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift)  / (2.0 * this->LaughlinIndex))
								       + (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift))  / (2.0 * this->LaughlinIndex))));
		      BMatrices[m].SetMatrixElement(this->GetMatrixIndex(N1, k1, this->NbrNValue, this->TotalStartingIndexPerPLevel[i]), 
						    this->GetMatrixIndex(N2, k2, this->NbrNValue, this->TotalStartingIndexPerPLevel[j]), Tmp);
		    }
		}
	    }
	}
    }

  delete[] Partition1;
  delete[] Partition2;

  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      this->RealBMatrices[i] = BMatrices[i];
    }
  delete[] BMatrices;
  for (int i = 0; i <= this->PLevel; ++i)
    delete U1BosonBasis[i];
  delete[] U1BosonBasis;
}

