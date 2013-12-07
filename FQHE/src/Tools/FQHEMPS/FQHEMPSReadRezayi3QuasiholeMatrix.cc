////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                           class author: Yang-Le Wu                         //
//                                                                            //
//      class of MPS matrix for quasihole in the Read-Rezayi k=3 state        //
//                                                                            //
//                        last modification : 07/12/2013                      //
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
#include "Tools/FQHEMPS/FQHEMPSReadRezayi3QuasiholeMatrix.h"
#include "GeneralTools/ConfigurationParser.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/LongRationalMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "Architecture/ArchitectureOperation/FQHEMPSEvaluateCFTOperation.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/Endian.h"

#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;


// constructor
//
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// cftDirectory = path to the directory where all the pure CFT matrices are stored
// useRational = use arbitrary precision numbers for all the CFT calculations
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation

FQHEMPSReadRezayi3QuasiholeMatrix::FQHEMPSReadRezayi3QuasiholeMatrix(int pLevel, int nbrBMatrices, char* cftDirectory,
        bool useRational, bool trimChargeIndices, bool cylinderFlag, double kappa, AbstractArchitecture* architecture)
{
    this->NbrBMatrices = nbrBMatrices;
    this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
    this->LaughlinIndex = 0; // nonsense
    this->PLevel = pLevel;
    this->CylinderFlag = cylinderFlag;
    this->UseRationalFlag = useRational;
    this->UniformChargeIndexRange = !trimChargeIndices;
    this->Kappa = kappa;

    this->CentralCharge = LongRational(4l, 5l);

    this->BMatrixOutputName = new char[256];
    sprintf(this->BMatrixOutputName, "readrezayi3_vacqh");
    this->CreateBMatrices(cftDirectory, architecture);
}


// destructor
//

FQHEMPSReadRezayi3QuasiholeMatrix::~FQHEMPSReadRezayi3QuasiholeMatrix()
{
}

// create the B matrices for the laughlin state
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// architecture = architecture to use for precalculation

void FQHEMPSReadRezayi3QuasiholeMatrix::CreateBMatrices(char* cftDirectory, AbstractArchitecture* architecture)
{
    int nbrSectors = 6;
    //                          0         1     2      3          4       5
    char* sectorNames[6] = {"identity", "psi", "w", "epsilon", "sigma", "phi"};
    LongRational weights[6] = {LongRational(0l, 1l), LongRational(2l, 3l), LongRational(3l, 1l),
        LongRational(2l, 5l), LongRational(1l, 15l), LongRational(7l, 5l)};

    int nbrChannels = 15;
    int* leftSectors = new int[nbrChannels];
    int* rightSectors = new int[nbrChannels];
    double* globalFactors = new double[nbrChannels];

    double sqrtc = sqrt(0.5) * exp(0.25 * (lgamma(0.2) + 3 * lgamma(0.6) - lgamma(0.8) - 3 * lgamma(0.4)));

    // <σ₁|σ₁|1>
    leftSectors  [0] = 4;
    rightSectors [0] = 0;
    globalFactors[0] = 1.0;

    // <ε|σ₁|ψ₁>
    leftSectors  [1] = 3;
    rightSectors [1] = 1;
    globalFactors[1] = sqrt(2.0 / 3);

    // <φ|σ₁|ψ₁>
    leftSectors  [2] = 5;
    rightSectors [2] = 1;
    globalFactors[2] = - sqrt(7.0 / 2) / 3;

    // <σ₂|σ₁|ψ₂>
    leftSectors  [3] = 4;
    rightSectors [3] = 1;
    globalFactors[3] = 1 / sqrt(3.0);

    // <σ₁|σ₁|W>
    leftSectors  [4] = 4;
    rightSectors [4] = 2;
    globalFactors[4] = 1.0; // FIXME

    // <ψ₂|σ₁|ε>
    leftSectors  [5] = 1;
    rightSectors [5] = 3;
    globalFactors[5] = sqrt(2.0 / 3);

    // <σ₁|σ₁|ε>
    leftSectors  [6] = 4;
    rightSectors [6] = 3;
    globalFactors[6] = sqrtc;

    // <ψ₁|σ₁|σ₁>
    leftSectors  [7] = 1;
    rightSectors [7] = 4;
    globalFactors[7] = 1 / sqrt(3.0);

    // <σ₂|σ₁|σ₁>
    leftSectors  [8] = 4;
    rightSectors [8] = 4;
    globalFactors[8] = sqrt(2) * sqrtc;

    // <1|σ₁|σ₂>
    leftSectors  [9] = 0;
    rightSectors [9] = 4;
    globalFactors[9] = 1.0;

    // <W|σ₁|σ₂>
    leftSectors  [10] = 2;
    rightSectors [10] = 4;
    globalFactors[10] = 1.0; // FIXME

    // <ε|σ₁|σ₂>
    leftSectors  [11] = 3;
    rightSectors [11] = 4;
    globalFactors[11] = sqrtc;

    // <φ|σ₁|σ₂>
    leftSectors  [12] = 5;
    rightSectors [12] = 4;
    globalFactors[12] = sqrtc / sqrt(21.0);

    // <ψ₂|σ₁|φ>
    leftSectors  [13] = 1;
    rightSectors [13] = 5;
    globalFactors[13] = 1.0; // FIXME

    // <σ₁|σ₁|φ>
    leftSectors  [14] = 4;
    rightSectors [14] = 5;
    globalFactors[14] = 1.0; // FIXME

    RealMatrix*** Matrices = this->ComputeMatrixElements(cftDirectory, architecture, "sigma", LongRational(1l, 15l),
            nbrSectors, sectorNames, weights, nbrChannels, leftSectors, rightSectors, globalFactors);

    char* TmpFileName = new char[512];
    for (int i = 0; i <= this->PLevel; ++i)
    {
        for (int j = 0; j <= this->PLevel; ++j)
        {
            for (int c = 0; c < nbrChannels; ++c)
            {
                sprintf(TmpFileName, "cft_%s_finalmatrixelement_%s%s_level_%d_%d.dat", this->BMatrixOutputName, sectorNames[leftSectors[c]], sectorNames[rightSectors[c]], i, j);
                Matrices[c][i][j].WriteMatrix(TmpFileName);
            }
        }
    }
    delete[] TmpFileName;

    delete[] leftSectors;
    delete[] rightSectors;
    delete[] globalFactors;

    for (int c = 0; c < nbrChannels; ++c)
    {
        for (int i = 0; i <= this->PLevel; ++i)
            delete[] Matrices[c][i];
        delete[] Matrices[c];
    }
    delete[] Matrices;
}

bool FQHEMPSReadRezayi3QuasiholeMatrix::SaveMatrices(char* fileName)
{
    cout << "there's no B matrices to save" << endl;
    return true;
}


