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

SparseComplexMatrix* FQHEMPSLaughlinQuasiholeMatrix::GetQuasiholeMatrices(int nbrQuasiholes, Complex* quasiholePositions)
{
  SparseComplexMatrix* TmpQuasiholeBMatrices = new SparseComplexMatrix[nbrQuasiholes];
  for (int i = 0; i < nbrQuasiholes; ++i)
    TmpQuasiholeBMatrices[i] = SparseComplexMatrix(this->QuasiholeBMatrices[0]);
  Complex* TmpCoefficient = new Complex[nbrQuasiholes];
  for (int p = 0; p <= this->PLevel; ++p)
    {
      int NbrParitions1 = this->NbrIndicesPerPLevel[p] / this->NbrNValue;
      for (int q = 0; q <= this->PLevel; ++q)
	{
	  int NbrParitions2 = this->NbrIndicesPerPLevel[q] / this->NbrNValue;	  
	  for (int N2 = 0; N2 < (this->NbrNValue - 1); ++N2)
	    {
	      int N1 = N2 + 1;
	      for (int i = 0; i < nbrQuasiholes; ++i)
		{
		  if (Norm(quasiholePositions[i]) > MACHINE_PRECISION)
		    TmpCoefficient[i] = pow (quasiholePositions[i], ((double) (2 * N1 - this->NbrNValue + 2)) / (20 * ((double) this->LaughlinIndex)) + ((double) (p - q)));
		  else
		    {
		      cout << N1 << " " << this->NbrNValue << " " << (2 * N1 - this->NbrNValue + 2) << endl;
		      if ((2 * N1 - this->NbrNValue - 1) == (2 * this->LaughlinIndex * (q - p)))
			{
			  TmpCoefficient[i] = 1.0;
			}
		      else
			{
			  TmpCoefficient[i] = 0.0;
			}
		    }
		}
	      for (int k1 = 0; k1 < NbrParitions1; ++k1)
		for (int k2 = 0; k2 < NbrParitions2; ++k2)
		  {
		    double Tmp;
		    int Index1 = this->GetMatrixIndex(N1, k1, this->NbrNValue, this->TotalStartingIndexPerPLevel[p]);
		    int Index2 = this->GetMatrixIndex(N2, k2, this->NbrNValue, this->TotalStartingIndexPerPLevel[q]);
		    this->QuasiholeBMatrices[0].GetMatrixElement(Index1, Index2, Tmp);
		    for (int i = 0; i < nbrQuasiholes; ++i)
		      {
			TmpQuasiholeBMatrices[i].SetMatrixElement(Index1, Index2, Tmp * TmpCoefficient[i]);
		      }
		  }
	    }
	}
    }
  delete[] TmpCoefficient;
  for (int i = 0; i < nbrQuasiholes; ++i)
    {
      cout << " final matrix " << i << " : ";
      TmpQuasiholeBMatrices[i].PrintNonZero(cout) << endl;
      cout << "--------------------------------" << endl;
    }
  return TmpQuasiholeBMatrices;
}

  
// create the B matrices for the laughlin state
//

void FQHEMPSLaughlinQuasiholeMatrix::CreateBMatrices ()
{
  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
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
  this->RealBMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize);
  this->QuasiholeBMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize);
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
	      this->RealBMatrices[0].SetMatrixElement(this->GetMatrixIndex(j - 1, k, this->NbrNValue, this->TotalStartingIndexPerPLevel[i]),
					    this->GetMatrixIndex(j, k, this->NbrNValue, this->TotalStartingIndexPerPLevel[i]), Tmp);
	    }
	}
    }
   
  for (int m = 1; m < this->NbrBMatrices; ++m)
    this->RealBMatrices[m] = SparseRealMatrix(MatrixSize, MatrixSize);
  this->QuasiholeBMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize);

  unsigned long* Partition1 = new unsigned long [this->PLevel + 2];
  unsigned long* Partition2 = new unsigned long [this->PLevel + 2];
  FactorialCoefficient Coef;

  for (int i = 0; i <= this->PLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace1 = U1BosonBasis[i];
      for (int j = 0; j <= this->PLevel; ++j)
	{
	  BosonOnDiskShort* TmpSpace2 = U1BosonBasis[j];
	  for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
	    {
	      TmpSpace1->GetOccupationNumber(k1, Partition1);
	      for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
		{
		  TmpSpace2->GetOccupationNumber(k2, Partition2);
		  int N2 = (2 * (j - i) - this->LaughlinIndex + 1 + NValueShift) / 2;
		  int N1 = N2 + (this->LaughlinIndex - 1);
		  double Tmp = this->CreateLaughlinAMatrixElement(this->LaughlinIndex, 1, Partition1, Partition2, i, j, Coef);
		  if (this->CylinderFlag)
		    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (((double) (i + j))
								   + ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift)  / (2.0 * this->LaughlinIndex))
								   + (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift))  / (2.0 * this->LaughlinIndex))));
		  this->RealBMatrices[1].SetMatrixElement(this->GetMatrixIndex(N1, k1, this->NbrNValue, this->TotalStartingIndexPerPLevel[i]), 
						this->GetMatrixIndex(N2, k2, this->NbrNValue, this->TotalStartingIndexPerPLevel[j]), Tmp);
		  Tmp = this->CreateLaughlinAMatrixElement(1, this->LaughlinIndex, Partition1, Partition2, i, j, Coef);
		  for (int N2 = 0; N2 < (this->NbrNValue - 1); ++N2)
		    {
		      int N1 = N2 + 1;
		      double Tmp2 = Tmp;
		      if (this->CylinderFlag)
			Tmp2 *= exp(-0.5 * this->Kappa * this->Kappa * (((double) (i + j))
									+ ((N1 - 0.5 * NValueShift) * (N1 - 0.5 * NValueShift)  / (2.0 * this->LaughlinIndex))
									+ (((N2 - 0.5 * NValueShift) * (N2 - 0.5 * NValueShift))  / (2.0 * this->LaughlinIndex))));
		      this->QuasiholeBMatrices[0].SetMatrixElement(this->GetMatrixIndex(N1, k1, this->NbrNValue, this->TotalStartingIndexPerPLevel[i]), 
								   this->GetMatrixIndex(N2, k2, this->NbrNValue, this->TotalStartingIndexPerPLevel[j]), Tmp2);
		    }		    
		}
	    }
	}
    }

  delete[] Partition1;
  delete[] Partition2;

  for (int i = 0; i <= this->PLevel; ++i)
    delete U1BosonBasis[i];
  delete[] U1BosonBasis;
}

