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
    LongRational CentralCharge12 (this->CentralCharge);
    cout << "central charge = " << CentralCharge12 << endl;
    CentralCharge12 /= 12l;
    double CentralCharge12Numerical = CentralCharge12.GetNumericalValue();

    double WeightIdentityNumerical = this->WeightIdentity.GetNumericalValue();
    double WeightPsiNumerical = this->WeightPsi.GetNumericalValue();
    double WeightSigmaNumerical = this->WeightSigma.GetNumericalValue();
    double WeightPhiNumerical = this->WeightPhi.GetNumericalValue();

    BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
    LongRational** RationalMultiplicityFactor = new LongRational*[this->PLevel + 1];
    double** MultiplicityFactor = new double*[this->PLevel + 1];
    unsigned long* TmpPartition = new unsigned long [this->PLevel + 2];
    for (int i = 0; i <= this->PLevel; ++i)
    {
        U1BosonBasis[i] = new BosonOnDiskShort(i, i, this->PLevel + 1);
        RationalMultiplicityFactor[i] = new LongRational[U1BosonBasis[i]->GetHilbertSpaceDimension()];
        MultiplicityFactor[i] = new double[U1BosonBasis[i]->GetHilbertSpaceDimension()];
        for (int j = 0; j < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++j)
        {
            U1BosonBasis[i]->GetOccupationNumber(j, TmpPartition);
            RationalMultiplicityFactor[i][j] = 1l;
            MultiplicityFactor[i][j] = 1.0;
            for (int k = 1; k <= i; ++k)
            {
                if (TmpPartition[k] > 1ul)
                {
                    RationalMultiplicityFactor[i][j].FactorialDivide(TmpPartition[k]);
                    double Tmp = 1.0;
                    for (unsigned long l = 2l; l <= TmpPartition[k]; ++l)
                        Tmp *=  (double) l;
                    MultiplicityFactor[i][j] /= Tmp;
                }
            }
        }
    }
    delete[] TmpPartition;

    RealSymmetricMatrix* ScalarProductIdentity = new RealSymmetricMatrix[this->PLevel + 1];
    RealSymmetricMatrix* ScalarProductPsi = new RealSymmetricMatrix[this->PLevel + 1];
    RealSymmetricMatrix* ScalarProductSigma = new RealSymmetricMatrix[this->PLevel + 1];
    RealSymmetricMatrix* ScalarProductPhi = new RealSymmetricMatrix[this->PLevel + 1];
    LongRationalMatrix* RationalScalarProductIdentity = new LongRationalMatrix[this->PLevel + 1];
    LongRationalMatrix* RationalScalarProductPsi = new LongRationalMatrix[this->PLevel + 1];
    LongRationalMatrix* RationalScalarProductSigma = new LongRationalMatrix[this->PLevel + 1];
    LongRationalMatrix* RationalScalarProductPhi = new LongRationalMatrix[this->PLevel + 1];
    RealMatrix* OrthogonalBasisIdentityLeft = new RealMatrix[this->PLevel + 1];
    RealMatrix* OrthogonalBasisIdentityRight = new RealMatrix[this->PLevel + 1];
    RealMatrix* OrthogonalBasisPsiLeft = new RealMatrix[this->PLevel + 1];
    RealMatrix* OrthogonalBasisPsiRight = new RealMatrix[this->PLevel + 1];
    RealMatrix* OrthogonalBasisSigmaLeft = new RealMatrix[this->PLevel + 1];
    RealMatrix* OrthogonalBasisSigmaRight = new RealMatrix[this->PLevel + 1];
    RealMatrix* OrthogonalBasisPhiLeft = new RealMatrix[this->PLevel + 1];
    RealMatrix* OrthogonalBasisPhiRight = new RealMatrix[this->PLevel + 1];

    cout << "computing vacuum sector scalar products" << endl;
    cout << "weight: " << this->WeightIdentity << " " << this->WeightPsi << endl;
    char* TmpScalarProductIdentityFileName = 0;
    char* TmpScalarProductPsiFileName = 0;
    if (cftDirectory != 0)
    {
        TmpScalarProductIdentityFileName = new char[512 + strlen(cftDirectory)];
        TmpScalarProductPsiFileName = new char[512 + strlen(cftDirectory)];
    }
    for (int i = 0; i <= this->PLevel; ++i)
    {
        cout << "Level = " <<  i << endl;
        if (cftDirectory != 0)
        {
            if (this->UseRationalFlag == true)
            {
                sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_scalarproducts_identity_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
                sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_scalarproducts_psi_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
            }
            else
            {
                sprintf (TmpScalarProductIdentityFileName, "%s/cft_%s_num_scalarproducts_identity_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
                sprintf (TmpScalarProductPsiFileName, "%s/cft_%s_num_scalarproducts_psi_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
            }
        }
        this->ComputeFullScalarProductMatrix(cftDirectory, TmpScalarProductIdentityFileName, architecture, RationalScalarProductIdentity, ScalarProductIdentity, i, U1BosonBasis,
                CentralCharge12, CentralCharge12Numerical, this->WeightIdentity, WeightIdentityNumerical, "identity",
                OrthogonalBasisIdentityLeft, OrthogonalBasisIdentityRight, RationalMultiplicityFactor, MultiplicityFactor);
        this->ComputeFullScalarProductMatrix(cftDirectory, TmpScalarProductPsiFileName, architecture, RationalScalarProductPsi, ScalarProductPsi, i, U1BosonBasis,
                CentralCharge12, CentralCharge12Numerical, this->WeightPsi, WeightPsiNumerical, "psi",
                OrthogonalBasisPsiLeft, OrthogonalBasisPsiRight, RationalMultiplicityFactor, MultiplicityFactor);
        cout << "---------------------------------" << endl;
    }
    delete[] TmpScalarProductIdentityFileName;
    delete[] TmpScalarProductPsiFileName;
    this->RescaleFullScalarProductMatrix(RationalScalarProductIdentity, ScalarProductIdentity, RationalMultiplicityFactor, MultiplicityFactor);
    this->RescaleFullScalarProductMatrix(RationalScalarProductPsi, ScalarProductPsi, RationalMultiplicityFactor, MultiplicityFactor);

    cout << "computing quasihole sector scalar products" << endl;
    cout << "weight: " << this->WeightSigma << " " << this->WeightPhi << endl;
    char* TmpScalarProductSigmaFileName = 0;
    char* TmpScalarProductPhiFileName = 0;
    if (cftDirectory != 0)
    {
        TmpScalarProductSigmaFileName = new char[512 + strlen(cftDirectory)];
        TmpScalarProductPhiFileName = new char[512 + strlen(cftDirectory)];
    }
    for (int i = 0; i <= this->PLevel; ++i)
    {
        cout << "Level = " <<  i << endl;
        if (cftDirectory != 0)
        {
            if (this->UseRationalFlag == true)
            {
                sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_scalarproducts_sigma_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
                if (this->SelfDualFlag == false)
                    sprintf (TmpScalarProductPhiFileName, "%s/cft_%s_scalarproducts_phi_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
            }
            else
            {
                sprintf (TmpScalarProductSigmaFileName, "%s/cft_%s_num_scalarproducts_sigma_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
                if (this->SelfDualFlag == false)
                    sprintf (TmpScalarProductPhiFileName, "%s/cft_%s_num_scalarproducts_phi_level_%d.dat", cftDirectory, this->BMatrixOutputName, i);
            }
        }
        this->ComputeFullScalarProductMatrix(cftDirectory, TmpScalarProductSigmaFileName, architecture, RationalScalarProductSigma, ScalarProductSigma, i, U1BosonBasis,
                CentralCharge12, CentralCharge12Numerical, this->WeightSigma, WeightSigmaNumerical, "sigma",
                OrthogonalBasisSigmaLeft, OrthogonalBasisSigmaRight, RationalMultiplicityFactor, MultiplicityFactor);
        if (this->SelfDualFlag == false)
        {
            this->ComputeFullScalarProductMatrix(cftDirectory, TmpScalarProductPhiFileName, architecture, RationalScalarProductPhi, ScalarProductPhi, i, U1BosonBasis,
                    CentralCharge12, CentralCharge12Numerical, this->WeightPhi, WeightPhiNumerical, "phi",
                    OrthogonalBasisPhiLeft, OrthogonalBasisPhiRight, RationalMultiplicityFactor, MultiplicityFactor);
        }

        cout << "---------------------------------" << endl;
    }
    delete[] TmpScalarProductSigmaFileName;
    delete[] TmpScalarProductPhiFileName;
    this->RescaleFullScalarProductMatrix(RationalScalarProductSigma, ScalarProductSigma, RationalMultiplicityFactor, MultiplicityFactor);
    if (this->SelfDualFlag == false)
        this->RescaleFullScalarProductMatrix(RationalScalarProductPhi, ScalarProductPhi, RationalMultiplicityFactor, MultiplicityFactor);

    cout << "computing Sigma matrix elements" << endl;
    RealMatrix** MatrixSigmaIdentitySigma = new RealMatrix*[this->PLevel + 1];
    RealMatrix** MatrixSigmaSigmaIdentity = new RealMatrix*[this->PLevel + 1];
    RealMatrix** MatrixSigmaPsiSigma = new RealMatrix*[this->PLevel + 1];
    RealMatrix** MatrixSigmaSigmaPsi = new RealMatrix*[this->PLevel + 1];
    LongRationalMatrix** RationalMatrixSigmaIdentitySigma = new LongRationalMatrix*[this->PLevel + 1];
    LongRationalMatrix** RationalMatrixSigmaSigmaIdentity = new LongRationalMatrix*[this->PLevel + 1];
    LongRationalMatrix** RationalMatrixSigmaPsiSigma = new LongRationalMatrix*[this->PLevel + 1];
    LongRationalMatrix** RationalMatrixSigmaSigmaPsi = new LongRationalMatrix*[this->PLevel + 1];

    for (int i = 0; i <= this->PLevel; ++i)
    {
        MatrixSigmaIdentitySigma[i] = new RealMatrix[this->PLevel + 1];
        MatrixSigmaSigmaIdentity[i] = new RealMatrix[this->PLevel + 1];
        MatrixSigmaPsiSigma[i] = new RealMatrix[this->PLevel + 1];
        MatrixSigmaSigmaPsi[i] = new RealMatrix[this->PLevel + 1];
        RationalMatrixSigmaIdentitySigma[i] = new LongRationalMatrix[this->PLevel + 1];
        RationalMatrixSigmaSigmaIdentity[i] = new LongRationalMatrix[this->PLevel + 1];
        RationalMatrixSigmaPsiSigma[i] = new LongRationalMatrix[this->PLevel + 1];
        RationalMatrixSigmaSigmaPsi[i] = new LongRationalMatrix[this->PLevel + 1];
    }

    char* TmpMatrixSigmaIdentitySigmaFileName = 0;
    char* TmpMatrixSigmaSigmaIdentityFileName = 0;
    char* TmpMatrixSigmaPsiSigmaFileName = 0;
    char* TmpMatrixSigmaSigmaPsiFileName = 0;
    if (cftDirectory != 0)
    {
        TmpMatrixSigmaIdentitySigmaFileName = new char[512 + strlen(cftDirectory)];;
        TmpMatrixSigmaSigmaIdentityFileName = new char[512 + strlen(cftDirectory)];;
        TmpMatrixSigmaPsiSigmaFileName = new char[512 + strlen(cftDirectory)];;
        TmpMatrixSigmaSigmaPsiFileName = new char[512 + strlen(cftDirectory)];;
    }

    // FIXME: not sure how to handle the non-self-dual case...
    if (this->SelfDualFlag == false)
    {
        cout << "only the self-dual case is implemented." << endl;
        exit(1);
    }

    char* TmpFilenameHeader = 0;
    if (cftDirectory != 0)
    {
        TmpFilenameHeader = new char[512 + strlen(cftDirectory)];
        if (this->UseRationalFlag == true)
            sprintf(TmpFilenameHeader, "%s/cft_%s_matrixelement", cftDirectory, this->BMatrixOutputName);
        else
            sprintf(TmpFilenameHeader, "%s/cft_%s_num_matrixelement", cftDirectory, this->BMatrixOutputName);
    }
    for (int j = 0; j <= this->PLevel; ++j)
    {
        for (int i = 0; i <= this->PLevel; ++i)
        {
            cout << "Levels = " <<  i << " " << j << endl;
            if (cftDirectory != 0)
            {
                sprintf(TmpMatrixSigmaIdentitySigmaFileName, "%s_identitysigma_level_%d_%d.dat", TmpFilenameHeader, i, j);
                sprintf(TmpMatrixSigmaSigmaIdentityFileName, "%s_sigmaidentity_level_%d_%d.dat", TmpFilenameHeader, i, j);
                sprintf(TmpMatrixSigmaPsiSigmaFileName, "%s_psisigma_level_%d_%d.dat", TmpFilenameHeader, i, j);
                sprintf(TmpMatrixSigmaSigmaPsiFileName, "%s_sigmapsi_level_%d_%d.dat", TmpFilenameHeader, i, j);
            }

            this->ComputeFullMatrixElements(cftDirectory, TmpMatrixSigmaIdentitySigmaFileName, architecture,
                    RationalMatrixSigmaIdentitySigma, MatrixSigmaIdentitySigma, i, j, U1BosonBasis,
                    CentralCharge12, CentralCharge12Numerical,
                    this->WeightIdentity, WeightIdentityNumerical,
                    this->WeightSigma, WeightSigmaNumerical,
                    this->WeightSigma, WeightSigmaNumerical);
            this->ComputeFullMatrixElements(cftDirectory, TmpMatrixSigmaSigmaIdentityFileName, architecture,
                    RationalMatrixSigmaSigmaIdentity, MatrixSigmaSigmaIdentity, i, j, U1BosonBasis,
                    CentralCharge12, CentralCharge12Numerical,
                    this->WeightSigma, WeightSigmaNumerical,
                    this->WeightIdentity, WeightIdentityNumerical,
                    this->WeightSigma, WeightSigmaNumerical);

            this->ComputeFullMatrixElements(cftDirectory, TmpMatrixSigmaPsiSigmaFileName, architecture,
                    RationalMatrixSigmaPsiSigma, MatrixSigmaPsiSigma, i, j, U1BosonBasis,
                    CentralCharge12, CentralCharge12Numerical,
                    this->WeightPsi, WeightPsiNumerical,
                    this->WeightSigma, WeightSigmaNumerical,
                    this->WeightSigma, WeightSigmaNumerical);
            this->ComputeFullMatrixElements(cftDirectory, TmpMatrixSigmaSigmaPsiFileName, architecture,
                    RationalMatrixSigmaSigmaPsi, MatrixSigmaSigmaPsi, i, j, U1BosonBasis,
                    CentralCharge12, CentralCharge12Numerical,
                    this->WeightSigma, WeightSigmaNumerical,
                    this->WeightPsi, WeightPsiNumerical,
                    this->WeightSigma, WeightSigmaNumerical);
        }
    }
    this->RescaleFullMatrixElements(RationalMatrixSigmaIdentitySigma, MatrixSigmaIdentitySigma, RationalMultiplicityFactor, MultiplicityFactor, 1.0);
    this->RescaleFullMatrixElements(RationalMatrixSigmaSigmaIdentity, MatrixSigmaSigmaIdentity, RationalMultiplicityFactor, MultiplicityFactor, 1.0);
    this->RescaleFullMatrixElements(RationalMatrixSigmaPsiSigma, MatrixSigmaPsiSigma, RationalMultiplicityFactor, MultiplicityFactor, 1.0 / M_SQRT2);
    this->RescaleFullMatrixElements(RationalMatrixSigmaSigmaPsi, MatrixSigmaSigmaPsi, RationalMultiplicityFactor, MultiplicityFactor, 1.0 / M_SQRT2);


    if (cftDirectory != 0)
        sprintf(TmpFilenameHeader, "%s/cft_%s_finalmatrixelement", cftDirectory, this->BMatrixOutputName);
    else
    {
        TmpFilenameHeader = new char[512];
        TmpMatrixSigmaIdentitySigmaFileName = new char[512];
        TmpMatrixSigmaSigmaIdentityFileName = new char[512];
        TmpMatrixSigmaPsiSigmaFileName = new char[512];
        TmpMatrixSigmaSigmaPsiFileName = new char[512];
        sprintf(TmpFilenameHeader, "cft_%s_finalmatrixelement", this->BMatrixOutputName);
    }
    for (int i = 0; i <= this->PLevel; ++i)
    {
        OrthogonalBasisIdentityLeft[i].Transpose();
        OrthogonalBasisSigmaLeft[i].Transpose();
        OrthogonalBasisPsiLeft[i].Transpose();
        for (int j = 0; j <= this->PLevel; ++j)
        {
            sprintf(TmpMatrixSigmaIdentitySigmaFileName, "%s_identitysigma_level_%d_%d.dat", TmpFilenameHeader, i, j);
            sprintf(TmpMatrixSigmaSigmaIdentityFileName, "%s_sigmaidentity_level_%d_%d.dat", TmpFilenameHeader, i, j);
            sprintf(TmpMatrixSigmaPsiSigmaFileName, "%s_psisigma_level_%d_%d.dat", TmpFilenameHeader, i, j);
            sprintf(TmpMatrixSigmaSigmaPsiFileName, "%s_sigmapsi_level_%d_%d.dat", TmpFilenameHeader, i, j);
            ((OrthogonalBasisIdentityLeft[i] * MatrixSigmaIdentitySigma[i][j]) * OrthogonalBasisSigmaRight[j]).WriteMatrix(TmpMatrixSigmaIdentitySigmaFileName);
            ((OrthogonalBasisSigmaLeft[i] * MatrixSigmaSigmaIdentity[i][j]) * OrthogonalBasisIdentityRight[j]).WriteMatrix(TmpMatrixSigmaSigmaIdentityFileName);
            ((OrthogonalBasisPsiLeft[i] * MatrixSigmaPsiSigma[i][j]) * OrthogonalBasisSigmaRight[j]).WriteMatrix(TmpMatrixSigmaPsiSigmaFileName);
            ((OrthogonalBasisSigmaLeft[i] * MatrixSigmaSigmaPsi[i][j]) * OrthogonalBasisPsiRight[j]).WriteMatrix(TmpMatrixSigmaSigmaPsiFileName);
        }
    }
    delete[] TmpFilenameHeader;
    delete[] TmpMatrixSigmaIdentitySigmaFileName;
    delete[] TmpMatrixSigmaSigmaIdentityFileName;
    delete[] TmpMatrixSigmaPsiSigmaFileName;
    delete[] TmpMatrixSigmaSigmaPsiFileName;

    delete[] ScalarProductIdentity;
    delete[] ScalarProductPsi;
    delete[] ScalarProductSigma;
    delete[] ScalarProductPhi;
    for (int i = 0; i <= this->PLevel; ++i)
    {
        delete U1BosonBasis[i];
        delete[] RationalMultiplicityFactor[i];
        delete[] MultiplicityFactor[i];
        delete[] MatrixSigmaIdentitySigma[i];
        delete[] MatrixSigmaSigmaIdentity[i];
        delete[] MatrixSigmaPsiSigma[i];
        delete[] MatrixSigmaSigmaPsi[i];
        delete[] RationalMatrixSigmaIdentitySigma[i];
        delete[] RationalMatrixSigmaSigmaIdentity[i];
        delete[] RationalMatrixSigmaPsiSigma[i];
        delete[] RationalMatrixSigmaSigmaPsi[i];
    }
    delete[] U1BosonBasis;
    delete[] OrthogonalBasisIdentityLeft;
    delete[] OrthogonalBasisPsiLeft;
    delete[] OrthogonalBasisSigmaLeft;
    delete[] OrthogonalBasisPhiLeft;
    delete[] OrthogonalBasisIdentityRight;
    delete[] OrthogonalBasisPsiRight;
    delete[] OrthogonalBasisSigmaRight;
    delete[] OrthogonalBasisPhiRight;
    delete[] RationalMultiplicityFactor;
    delete[] MultiplicityFactor;
    delete[] MatrixSigmaIdentitySigma;
    delete[] MatrixSigmaSigmaIdentity;
    delete[] MatrixSigmaPsiSigma;
    delete[] MatrixSigmaSigmaPsi;
    delete[] RationalMatrixSigmaIdentitySigma;
    delete[] RationalMatrixSigmaSigmaIdentity;
    delete[] RationalMatrixSigmaPsiSigma;
    delete[] RationalMatrixSigmaSigmaPsi;
}

bool FQHEMPSClustered2RQuasiholeMatrix::SaveMatrices(char* fileName)
{
    cout << "there's no B matrices to save" << endl;
    return true;
}
