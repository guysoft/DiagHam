////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//     class of MPS matrix for quasihole in the clustered (k=2,r) states      //
//                                                                            //
//                        last modification : 25/06/2013                      //
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
#include "Tools/FQHEMPS/FQHEMPSClustered2RQuasiholeMatrix.h"
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


// default constructor
//

FQHEMPSClustered2RQuasiholeMatrix::FQHEMPSClustered2RQuasiholeMatrix()
{
}

// constructor
//
// rindex = r index (i.e. clustered (k=2,r) states)
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// useRational = use arbitrary precision numbers for all the CFT calculations
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation

FQHEMPSClustered2RQuasiholeMatrix::FQHEMPSClustered2RQuasiholeMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, bool useRational,
        bool trimChargeIndices, bool cylinderFlag, double kappa, AbstractArchitecture* architecture)
{
    this->NbrBMatrices = nbrBMatrices;
    this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
    this->RIndex = rIndex;
    this->LaughlinIndex = laughlinIndex;
    this->PLevel = pLevel;
    this->CylinderFlag = cylinderFlag;
    this->UseRationalFlag = useRational;
    this->UniformChargeIndexRange = !trimChargeIndices;
    this->Kappa = kappa;
    this->CentralCharge = LongRational ((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
    this->SelfDualFlag = ((this->RIndex & 1) == 0);
    this->WeightIdentity = LongRational(0l, 1l);
    this->WeightPsi = LongRational(this->RIndex, 4l);
    this->WeightSigma = LongRational(5l - (2l * this->RIndex), 4l * (this->RIndex + 2l));
    this->WeightPhi = LongRational((this->RIndex - 1l) * (this->RIndex - 1l), 4l * (this->RIndex + 2l));

    this->WeightPrimaryFieldMatrixElement = WeightSigma;
    this->MatrixElementNormalization = 1.0 / M_SQRT2;
    this->SquareMatrixElementNormalization = LongRational(1, 2);

    this->TransferMatrixDegeneracy = this->RIndex + 2;
    this->NbrCFTSectors = 2;
    if (this->SelfDualFlag == true)
    {
        this->TransferMatrixDegeneracy /= 2;
        this->NbrCFTSectors = 1;
    }
    this->BMatrixOutputName = new char[256];
    sprintf(this->BMatrixOutputName, "clustered_k_2_vacqh_r_%d", this->RIndex);
    this->CreateBMatrices(0, architecture);
}

// constructor
//
// rindex = r index (i.e. clustered (k=2,r) states)
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// cftDirectory = path to the directory where all the pure CFT matrices are stored
// useRational = use arbitrary precision numbers for all the CFT calculations
// trimChargeIndices = trim the charge indices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation

FQHEMPSClustered2RQuasiholeMatrix::FQHEMPSClustered2RQuasiholeMatrix(int rIndex, int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory,
        bool useRational, bool trimChargeIndices,
        bool cylinderFlag, double kappa, AbstractArchitecture* architecture)
{
    this->NbrBMatrices = nbrBMatrices;
    this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
    this->RIndex = rIndex;
    this->LaughlinIndex = laughlinIndex;
    this->PLevel = pLevel;
    this->CylinderFlag = cylinderFlag;
    this->Kappa = kappa;
    this->CentralCharge = LongRational ((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
    this->SelfDualFlag = ((this->RIndex & 1) == 0);
    this->WeightIdentity = LongRational(0l, 1l);
    this->WeightPsi = LongRational(this->RIndex, 4l);
    this->WeightSigma = LongRational(5l - (2l * this->RIndex), 4l * (this->RIndex + 2l));
    this->WeightPhi = LongRational((this->RIndex - 1l) * (this->RIndex - 1l), 4l * (this->RIndex + 2l));

    this->WeightPrimaryFieldMatrixElement = WeightSigma;
    this->MatrixElementNormalization = 1.0 / M_SQRT2;
    this->SquareMatrixElementNormalization = LongRational(1, 2);

    this->UseRationalFlag = useRational;
    this->UniformChargeIndexRange = !trimChargeIndices;
    this->TransferMatrixDegeneracy = this->RIndex + 2;
    this->NbrCFTSectors = 2;
    if (this->SelfDualFlag == true)
    {
        this->TransferMatrixDegeneracy /= 2;
        this->NbrCFTSectors = 1;
    }
    this->BMatrixOutputName = new char[256];
    sprintf(this->BMatrixOutputName, "clustered_k_2_vacqh_r_%d", this->RIndex);
    this->CreateBMatrices(cftDirectory, architecture);
}

// constructor from a file describing the state
//
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital + 1)
// fileName = name of the file that contains the state description
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation

FQHEMPSClustered2RQuasiholeMatrix::FQHEMPSClustered2RQuasiholeMatrix(int pLevel, int nbrBMatrices, char* fileName, bool cylinderFlag, double kappa,
        AbstractArchitecture* architecture)
{
    this->NbrBMatrices = nbrBMatrices;
    this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
    this->CylinderFlag = cylinderFlag;
    this->Kappa = kappa;
    this->PLevel = pLevel;

    ConfigurationParser StateDefinition;
    if (StateDefinition.Parse(fileName) == false)
    {
        StateDefinition.DumpErrors(cout) << endl;
    }
    else
    {
        bool ErrorFlag = true;
        ErrorFlag = StateDefinition.GetAsSingleInteger("RIndex", this->RIndex);
        ErrorFlag = StateDefinition.GetAsSingleInteger("LaughlinIndex", this->LaughlinIndex);
        ErrorFlag = StateDefinition.GetAsSingleLongRational("WeightIdentity", this->WeightIdentity);
        ErrorFlag = StateDefinition.GetAsSingleLongRational("WeightPsi", this->WeightPsi);
        ErrorFlag = StateDefinition.GetAsSingleLongRational("WeightSigma", this->WeightSigma);
        ErrorFlag = StateDefinition.GetAsSingleLongRational("WeightPhi", this->WeightPhi);
        ErrorFlag = StateDefinition.GetAsBoolean("SelfDual", this->SelfDualFlag);
        ErrorFlag = StateDefinition.GetAsSingleLongRational("CentralCharge", this->CentralCharge);
        if (StateDefinition["PsiSquareMatrixElement"] != 0)
        {
            ErrorFlag = StateDefinition.GetAsSingleLongRational("PsiSquareMatrixElement", this->SquareMatrixElementNormalization);
        }
        else
        {
            this->SquareMatrixElementNormalization = LongRational(1, 2);
        }
        this->MatrixElementNormalization = sqrt(fabs(this->SquareMatrixElementNormalization.GetNumericalValue()));
        if (StateDefinition["EMatrixDegeneracy"] != 0)
        {
            ErrorFlag = StateDefinition.GetAsSingleInteger("EMatrixDegeneracy", this->TransferMatrixDegeneracy);
        }
        else
        {
            switch (this->RIndex)
            {
                case 2:
                    this->TransferMatrixDegeneracy = 2;
                    break;
                case 3:
                    this->TransferMatrixDegeneracy = 5;
                    break;
                case 6:
                    this->TransferMatrixDegeneracy = 5;
                    break;
            }
        }
        if (StateDefinition["Name"] != 0)
        {
            this->BMatrixOutputName = new char[strlen(StateDefinition["Name"]) + 1];
            strcpy(this->BMatrixOutputName, StateDefinition["Name"]);
        }
        else
        {
            this->BMatrixOutputName = new char[256];
            sprintf(this->BMatrixOutputName, "clustered_k_2_r_%d", this->RIndex);
        }
        this->NbrCFTSectors = 2;
        if (this->SelfDualFlag == true)
        {
            this->NbrCFTSectors = 1;
        }
        if (ErrorFlag == true)
        {
            this->WeightPrimaryFieldMatrixElement = this->WeightSigma;
            this->CreateBMatrices(StateDefinition["CFTMatrixDirectory"], architecture);
        }
    }
}


// constructor from stored B matrices
//
// rindex = r index (i.e. clustered (k=2,r) states)
// laughlinIndex = power of the Laughlin part (i.e.  laughlinIndex=2 for the fermionic MR at nu=1/2)
// pLevel = |P| level truncation
// fileName = name of the file that contains the B matrices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSClustered2RQuasiholeMatrix::FQHEMPSClustered2RQuasiholeMatrix(int rIndex, int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag, double kappa)
{
    this->RIndex = rIndex;
    this->LaughlinIndex = laughlinIndex;
    this->PLevel = pLevel;
    this->CylinderFlag = cylinderFlag;
    this->Kappa = kappa;
    this->LoadMatrices(fileName);
    this->CentralCharge = LongRational ((this->RIndex + 2l) - (2l * (this->RIndex - 1l) * (this->RIndex - 1l)), this->RIndex + 2l);
    this->SelfDualFlag = ((this->RIndex & 1) == 0);
    this->WeightIdentity = LongRational(0l, 1l);
    this->WeightPsi = LongRational(this->RIndex, 4l);
    this->WeightSigma = LongRational(5l - (2l * this->RIndex), 4l * (this->RIndex + 2l));
    this->WeightPhi = LongRational((this->RIndex - 1l) * (this->RIndex - 1l), 4l * (this->RIndex + 2l));

    this->WeightPrimaryFieldMatrixElement = WeightSigma;
    this->MatrixElementNormalization = 1.0 / M_SQRT2;
    this->SquareMatrixElementNormalization = LongRational(1, 2);

    this->TransferMatrixDegeneracy = this->RIndex + 2;
    this->BMatrixOutputName = new char[256];
    sprintf(this->BMatrixOutputName, "clustered_k_2_vacqh_r_%d", this->RIndex);
}

// destructor
//

FQHEMPSClustered2RQuasiholeMatrix::~FQHEMPSClustered2RQuasiholeMatrix()
{
}

// create the B matrices for the laughlin state
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// architecture = architecture to use for precalculation

void FQHEMPSClustered2RQuasiholeMatrix::CreateBMatrices(char* cftDirectory, AbstractArchitecture* architecture)
{
    int nbrSectors = 4;
    char* sectorNames[4] = {"identity", "psi", "sigma", "phi"};
    LongRational weights[4] = {this->WeightIdentity, this->WeightPsi, this->WeightSigma, this->WeightPhi};
    if (this->SelfDualFlag)
        nbrSectors = 3;

    int nbrChannels = (this->SelfDualFlag) ? 4 : 6;
    int* leftSectors = new int[nbrChannels];
    int* rightSectors = new int[nbrChannels];
    double* globalFactors = new double[nbrChannels];

    if (this->SelfDualFlag)
    {
        // <σ|σ|1>
        leftSectors[0] = 2;
        rightSectors[0] = 0;
        globalFactors[0] = 1.0;

        // <σ|σ|ψ>
        leftSectors[1] = 2;
        rightSectors[1] = 1;
        globalFactors[1] = 1.0 / M_SQRT2;

        // <1|σ|σ>
        leftSectors[2] = 0;
        rightSectors[2] = 2;
        globalFactors[2] = 1.0;

        // <ψ|σ|σ>
        leftSectors[3] = 1;
        rightSectors[3] = 2;
        globalFactors[3] = 1.0 / M_SQRT2;
    }
    else
    {
        // <σ|σ|1>
        leftSectors[0] = 2;
        rightSectors[0] = 0;
        globalFactors[0] = 1.0;

        // <φ|σ|ψ>
        leftSectors[1] = 3;
        rightSectors[1] = 1;
        globalFactors[1] = 1.0 / M_SQRT2; // FIXME

        // <1|σ|σ>
        leftSectors[2] = 0;
        rightSectors[2] = 2;
        globalFactors[2] = 1.0;

        // <φ|σ|σ>
        leftSectors[3] = 3;
        rightSectors[3] = 2;
        globalFactors[3] = 1.0; // FIXME

        // <σ|σ|φ>
        leftSectors[4] = 2;
        rightSectors[4] = 3;
        globalFactors[4] = 1.0; // FIXME

        // <ψ|σ|φ>
        leftSectors[5] = 1;
        rightSectors[5] = 3;
        globalFactors[5] = 1.0 / M_SQRT2; // FIXME
    }

    RealMatrix*** Matrices = this->ComputeMatrixElements(cftDirectory, architecture, "sigma", this->WeightSigma,
            nbrSectors, sectorNames, weights,
            nbrChannels, leftSectors, rightSectors, globalFactors);

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

bool FQHEMPSClustered2RQuasiholeMatrix::SaveMatrices(char* fileName)
{
    cout << "there's no B matrices to save" << endl;
    return true;
}

