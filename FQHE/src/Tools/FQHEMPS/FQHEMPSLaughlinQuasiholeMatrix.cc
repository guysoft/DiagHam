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

FQHEMPSLaughlinQuasiholeMatrix::FQHEMPSLaughlinQuasiholeMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, bool trimChargeIndices, 
							       bool cylinderFlag, double kappa)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->QuasiholeBMatrices = 0;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->NbrNValue = ((2 * this->PLevel) + this->LaughlinIndex);
  this->NValueGlobalShift = this->PLevel;
  this->UniformChargeIndexRange = !trimChargeIndices;
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

FQHEMPSLaughlinQuasiholeMatrix::FQHEMPSLaughlinQuasiholeMatrix(int laughlinIndex, int pLevel, char* fileName, bool trimChargeIndices, 
							       bool cylinderFlag, double kappa)
{
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->NbrNValue = ((2 * this->PLevel) + this->LaughlinIndex);
  this->NValueGlobalShift = this->PLevel;
  this->UniformChargeIndexRange = !trimChargeIndices;
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

// get the edge matrix for localized quasiholes, with normal ordering
//
// nbrQuasiholes = number of quasiholes
// quasiholePositions = quasihole positions in unit of magnetic length
// return value = pointer to the edge matrix

SparseComplexMatrix* FQHEMPSLaughlinQuasiholeMatrix::GetQuasiholeMatrices(int nbrQuasiholes, Complex* quasiholePositions)
{
    SparseComplexMatrix* QuasiholeMatrix = new SparseComplexMatrix(1, this->RealBMatrices[0].GetNbrRow());

    BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
    for (int p = 0; p <= this->PLevel; ++p)
        U1BosonBasis[p] = new BosonOnDiskShort(p, p, this->PLevel + 1);
    unsigned long* Monomial = new unsigned long[this->PLevel];

    FactorialCoefficient Coef;
    unsigned long* Occupation = new unsigned long [this->PLevel + 2];
    unsigned long* ZeroOccupation = new unsigned long [this->PLevel + 2];
    for (int i = 0; i < this->PLevel + 2; ++i)
        ZeroOccupation[i] = 0;

    int N = this->NValueGlobalShift + this->LaughlinIndex - 1 - nbrQuasiholes;
    for (int p = 0; p <= this->PLevel; ++p)
    {
        if ((N < this->NInitialValuePerPLevel[p]) || (N > this->NLastValuePerPLevel[p]))
            continue;
        BosonOnDiskShort* Space = U1BosonBasis[p];
        for (int k = 0; k < Space->GetHilbertSpaceDimension(); ++k)
        {
            Space->GetOccupationNumber(k, Occupation);
            Complex coeff = this->CreateLaughlinAMatrixElement(1, this->LaughlinIndex, ZeroOccupation, Occupation, 0, p, Coef);

            Space->GetMonomial(k, Monomial);
            for (int i = 0; (i < p) && (Monomial[i] > 0); ++i) // partition ENDS at the first zero in Monomial
            {
                Complex sum = 0;
                if (this->CylinderFlag)
                {
                    for (int a = 0; a < nbrQuasiholes; ++a)
                        sum += exp(this->Kappa * quasiholePositions[a] * Monomial[i]);
                }
                else
                {
                    for (int a = 0; a < nbrQuasiholes; ++a)
                        sum += pow(quasiholePositions[a], ((double)Monomial[i])); // don't add minus sign to unsigned long...
                }
                coeff *= sum;
            }

            QuasiholeMatrix->SetMatrixElement(0, this->GetMatrixIndex(p, k, N), coeff);
        }
    }

    // cout << "--- quasihole matrix" << endl;
    // QuasiholeMatrix->PrintNonZero(cout);
    // cout << "---" <<endl;
    cout << "quasihole matrix size = " << QuasiholeMatrix->GetNbrRow() << "x" << QuasiholeMatrix->GetNbrColumn() << endl;
    for (int p = 0; p <= this->PLevel; ++p)
        delete U1BosonBasis[p];
    delete[] U1BosonBasis;
    delete[] Monomial;
    delete[] Occupation;
    delete[] ZeroOccupation;
    return QuasiholeMatrix;
}

// get the boundary indices of the MPS representation
//
// rowIndex = matrix row index
// columnIndex = matrix column index
// padding = assume that the state has the estra padding

void FQHEMPSLaughlinQuasiholeMatrix::GetMatrixBoundaryIndices(int& rowIndex, int& columnIndex, bool padding)
{
    int MinQ;
    int MaxQ;
    this->GetChargeIndexRange(0, MinQ, MaxQ);
    if (padding == true)
        cout << "padding is not supported!" << endl;
    rowIndex = 0;
    columnIndex = this->NValueGlobalShift - MinQ;
}

