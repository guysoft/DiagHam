////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of MPS matrix for the Read-Rezayi k=3 state             //
//                           in its quasihole sector                          //
//                                                                            //
//                        last modification : 22/02/2013                      //
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
#include "Tools/FQHEMPS/FQHEMPSReadRezayi3QuasiholeSectorMatrix.h"
#include "Matrix/SparseRealMatrix.h"
#include "Matrix/LongRationalMatrix.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "Architecture/ArchitectureOperation/FQHEMPSEvaluateCFTOperation.h"

#include "GeneralTools/FilenameTools.h"



// constructor 
//
// laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital)
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation

FQHEMPSReadRezayi3QuasiholeSectorMatrix::FQHEMPSReadRezayi3QuasiholeSectorMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, bool cylinderFlag, double kappa, 
										 AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = 2;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->CreateBMatrices(0, architecture);
}

// constructor 
//
// laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic MR at nu=1/2)  
// pLevel = |P| level truncation
// nbrBMatrices = number of B matrices to compute (max occupation per orbital)
// cftDirectory = path to the directory where all the pure CFT matrices are stored
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// architecture = architecture to use for precalculation

FQHEMPSReadRezayi3QuasiholeSectorMatrix::FQHEMPSReadRezayi3QuasiholeSectorMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool cylinderFlag, double kappa, 
										 AbstractArchitecture* architecture)
{
  this->NbrBMatrices = nbrBMatrices;
  this->RealBMatrices = new SparseRealMatrix [this->NbrBMatrices];
  this->RIndex = 2;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->CreateBMatrices(cftDirectory, architecture);
}

// constructor from stored B matrices
//
// laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic RR state)  
// pLevel = |P| level truncation
// fileName = name of the file that contains the B matrices
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

FQHEMPSReadRezayi3QuasiholeSectorMatrix::FQHEMPSReadRezayi3QuasiholeSectorMatrix(int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag, double kappa)
{
  this->RIndex = 2;
  this->LaughlinIndex = laughlinIndex;
  this->PLevel = pLevel;
  this->CylinderFlag = cylinderFlag;
  this->Kappa = kappa;
  this->LoadMatrices(fileName);
}


// destructor
//

FQHEMPSReadRezayi3QuasiholeSectorMatrix::~FQHEMPSReadRezayi3QuasiholeSectorMatrix()
{
}
  
// get the name describing the B matrices 
// 
// return value = name 

char* FQHEMPSReadRezayi3QuasiholeSectorMatrix::GetName()
{
  char* TmpName = new char[24];
  sprintf (TmpName, "readrezayi3_qh", this->LaughlinIndex);
  return TmpName;
}

// create the B matrices for the laughlin state
//
// cftDirectory = an optional path to the directory where all the CFT matrices are stored
// architecture = architecture to use for precalculation

void FQHEMPSReadRezayi3QuasiholeSectorMatrix::CreateBMatrices (char* cftDirectory, AbstractArchitecture* architecture)
{
  LongRational CentralCharge (4l, 5l);
  cout << "central charge = " << CentralCharge << endl;
  LongRational CentralCharge12(CentralCharge);
  CentralCharge12 /= 12l;
  LongRational WeightPsi (2l, 3l);
  LongRational WeightSigma (1l, 15l);
  LongRational WeightPhi (7l, 5l);
  LongRational WeightEpsilon (2l, 5l);
  double WeightSigmaNumerical = WeightSigma.GetNumericalValue();
  double WeightPhiNumerical = WeightPhi.GetNumericalValue();
  double WeightEpsilonNumerical = WeightEpsilon.GetNumericalValue();
  long* Partition = new long[2 * (this->PLevel + 1)];
  unsigned long* TmpPartition = new unsigned long [this->PLevel + 2];
  this->TemporaryOccupationNumber = new unsigned long [this->PLevel + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductSigma = new RealSymmetricMatrix[this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductPhi = new RealSymmetricMatrix[this->PLevel + 1];
  RealSymmetricMatrix* ScalarProductEpsilon = new RealSymmetricMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductSigma = new LongRationalMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductPhi = new LongRationalMatrix[this->PLevel + 1];
  LongRationalMatrix* RationalScalarProductEpsilon = new LongRationalMatrix[this->PLevel + 1];
  RealMatrix** MatrixPsi01 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi10 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi11 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi12 = new RealMatrix*[this->PLevel + 1];
  RealMatrix** MatrixPsi21 = new RealMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi01 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi10 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi11 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi12 = new LongRationalMatrix*[this->PLevel + 1];
  LongRationalMatrix** RationalMatrixPsi21 = new LongRationalMatrix*[this->PLevel + 1];
  RealMatrix* OrthogonalBasisSigmaLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPhiLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisEpsilonLeft = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisSigmaRight = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisPhiRight = new RealMatrix[this->PLevel + 1];
  RealMatrix* OrthogonalBasisEpsilonRight = new RealMatrix[this->PLevel + 1];

  for (int i = 0; i <= this->PLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort (i, i, this->PLevel + 1);
      MatrixPsi01[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi10[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi11[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi21[i] = new RealMatrix[this->PLevel + 1];
      MatrixPsi12[i] = new RealMatrix[this->PLevel + 1];
      RationalMatrixPsi01[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi10[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi11[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi21[i] = new LongRationalMatrix[this->PLevel + 1];
      RationalMatrixPsi12[i] = new LongRationalMatrix[this->PLevel + 1];
    }
  
  char* TmpScalarProductSigmaFileName = 0; 
  char* TmpScalarProductEpsilonFileName = 0;
  char* TmpScalarProductPhiFileName = 0;
  char* TmpMatrixElementSigmaEpsilonFileName = 0;
  char* TmpMatrixElementEpsilonSigmaFileName = 0;
  char* TmpMatrixElementSigmaSigmaFileName = 0;
  char* TmpMatrixElementPhiSigmaFileName = 0;
  char* TmpMatrixElementSigmaPhiFileName = 0;
  if (cftDirectory != 0)
    {
      TmpScalarProductSigmaFileName = new char[512 + strlen(cftDirectory)];
      TmpScalarProductEpsilonFileName = new char[512 + strlen(cftDirectory)];
      TmpScalarProductPhiFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementSigmaEpsilonFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementEpsilonSigmaFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementSigmaSigmaFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementPhiSigmaFileName = new char[512 + strlen(cftDirectory)];
      TmpMatrixElementSigmaPhiFileName = new char[512 + strlen(cftDirectory)];
    }
  for (int i = 0; i <= this->PLevel; ++i)
    {
      cout << "Level = " <<  i << endl;
      RationalScalarProductSigma[i] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      RationalScalarProductEpsilon[i] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      if ((1 + i) <= this->PLevel)
	RationalScalarProductPhi[1 + i] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      if (i < 1)
	RationalScalarProductPhi[i] = LongRationalMatrix();
      if (cftDirectory != 0)
	{
	  sprintf (TmpScalarProductSigmaFileName, "%s/cft_readrezayi3_scalarproducts_sigma_level_%d.dat", cftDirectory, i);
	  sprintf (TmpScalarProductEpsilonFileName, "%s/cft_readrezayi3_scalarproducts_epsilon_level_%d.dat", cftDirectory, i);	  
	  if ((1 + i) <= this->PLevel)
	    sprintf (TmpScalarProductPhiFileName, "%s/cft_readrezayi3_scalarproducts_phi_level_%d.dat", cftDirectory, (i + 1));
	}
       if ((cftDirectory != 0) && (IsFile(TmpScalarProductSigmaFileName)) && (IsFile(TmpScalarProductEpsilonFileName))
	   && (((1 + i) > this->PLevel) || (IsFile(TmpScalarProductPhiFileName))))
	{
	  RationalScalarProductSigma[i].ReadMatrix(TmpScalarProductSigmaFileName);
	  RationalScalarProductEpsilon[i].ReadMatrix(TmpScalarProductEpsilonFileName);
	  if ((1 + i) <= this->PLevel)
	    RationalScalarProductPhi[3 + i].ReadMatrix(TmpScalarProductPhiFileName);
	}
       else
	 {
	  if (architecture == 0)
	    {
	      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
		for (int m = n; m < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++m)
		  {
		    int PartitionLength = 0;
		    U1BosonBasis[i]->GetOccupationNumber(n, TmpPartition);	    
		    for (int k = 1; k <= i; ++k)
		      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			{
			  ++PartitionLength;
			}
		    int Position = PartitionLength;
		    PartitionLength = 0;
		    for (int k = 1; k <= i; ++k)
		      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			{
			  Partition[Position - PartitionLength - 1] = (long) k;
			  ++PartitionLength;
			}
		    U1BosonBasis[i]->GetOccupationNumber(m, TmpPartition);	    
		    for (int k = 1; k <= i; ++k)
		      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			{
			  Partition[PartitionLength] = -(long) k;
			  ++PartitionLength;		  
			}
		    LongRational Tmp = this->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightSigma,
										     RationalScalarProductSigma, i - 1, U1BosonBasis, this->TemporaryOccupationNumber);
		    RationalScalarProductSigma[i].SetMatrixElement(m, n, Tmp);
		    if (n != m)
		      {
			RationalScalarProductSigma[i].SetMatrixElement(n, m, Tmp);	      
		      }
		    if ((1 + i) <= this->PLevel)
		      {
			Tmp = this->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightPhi);
			RationalScalarProductPhi[1 + i].SetMatrixElement(m, n, Tmp);
			if (n != m)
			  {
			    RationalScalarProductPhi[1 + i].SetMatrixElement(n, m, Tmp);	      
			  }
		      }
		    Tmp = this->ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightEpsilon,
									RationalScalarProductEpsilon, i - 1, U1BosonBasis, this->TemporaryOccupationNumber);
		    RationalScalarProductEpsilon[i].SetMatrixElement(m, n, Tmp);
		    if (n != m)
		      {
			RationalScalarProductEpsilon[i].SetMatrixElement(n, m, Tmp);	      
		      }
		  }
	      if (cftDirectory != 0)
		{
		  RationalScalarProductSigma[i].WriteMatrix(TmpScalarProductSigmaFileName);
		  RationalScalarProductEpsilon[i].WriteMatrix(TmpScalarProductEpsilonFileName);
		  if ((1 + i) <= this->PLevel)
		    RationalScalarProductPhi[1 + i].WriteMatrix(TmpScalarProductPhiFileName);
		}
	    }
	  else
	    {
	      if ((cftDirectory != 0) && (IsFile(TmpScalarProductSigmaFileName)))
		{		
		  RationalScalarProductSigma[i].ReadMatrix(TmpScalarProductSigmaFileName);
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, CentralCharge12, 
							 WeightSigma,
							 RationalScalarProductSigma,  i- 1);
		  Operation1.ApplyOperation(architecture);
		  RationalScalarProductSigma[i] = Operation1.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalScalarProductSigma[i].WriteMatrix(TmpScalarProductSigmaFileName);
		    }
		}
	      if ((cftDirectory != 0) && (IsFile(TmpScalarProductEpsilonFileName)))
		{
		  RationalScalarProductEpsilon[i].ReadMatrix(TmpScalarProductEpsilonFileName);
		}
	      else
		{
		  FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, CentralCharge12, 
							 WeightEpsilon,
							 RationalScalarProductEpsilon,  i - 1);
		  Operation2.ApplyOperation(architecture);
		  RationalScalarProductEpsilon[i] = Operation2.GetRationalMatrixElements();
		  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
		    {
		      RationalScalarProductEpsilon[i].WriteMatrix(TmpScalarProductEpsilonFileName);
		    }
		}
	      if ((1 + i) <= this->PLevel)
		{
		  if ((cftDirectory != 0) && (IsFile(TmpScalarProductPhiFileName)))
		    {
		      RationalScalarProductPhi[1 + i].ReadMatrix(TmpScalarProductPhiFileName);
		    }
		  else
		    {
		      FQHEMPSEvaluateCFTOperation Operation2(this, U1BosonBasis, i, CentralCharge12, 
							     WeightPhi,
							     RationalScalarProductPhi + 1,  i - 1);
		      Operation2.ApplyOperation(architecture);
		      RationalScalarProductPhi[1 + i] = Operation2.GetRationalMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  RationalScalarProductPhi[1 + i].WriteMatrix(TmpScalarProductPhiFileName);
			}
		    }
		}
	    }
	  }
      ScalarProductSigma[i] = RationalScalarProductSigma[i];      
      ScalarProductEpsilon[i] = RationalScalarProductEpsilon[i];      
      if ((1 + i) <= this->PLevel)
	ScalarProductPhi[1 + i] = RationalScalarProductPhi[1 + i];
      if (i < 1)
	ScalarProductPhi[i] = RealSymmetricMatrix();
      RealSymmetricMatrix TmpMatrix;
      TmpMatrix.Copy(ScalarProductSigma[i]);
      RealMatrix TmpBasis(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension());
      TmpBasis.SetToIdentity();
      RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
      TmpMatrix.LapackDiagonalize(TmpDiag, TmpBasis);
#else
      TmpMatrix.Diagonalize(TmpDiag, TmpBasis);
#endif
      double Error = 0.0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) > Error)
	  Error = fabs(TmpDiag(n, n));
      Error *= 1e-14;
      if (Error < 1e-14)
	Error = 1e-14;
      int Count  = 0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) < Error)
	  ++Count;
      cout << "nbr of null vectors sigma sector = " << Count << " (" << (U1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;
      if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	{
	  OrthogonalBasisSigmaLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  OrthogonalBasisSigmaRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisSigmaLeft[i][Count].Copy(TmpBasis[n]);
		OrthogonalBasisSigmaRight[i][Count].Copy(TmpBasis[n]);
		if (TmpDiag(n, n) > 0)
		  {
		    OrthogonalBasisSigmaLeft[i][Count] /=  sqrt(TmpDiag(n, n));
		    OrthogonalBasisSigmaRight[i][Count] /=  sqrt(TmpDiag(n, n));
		  }
		else
		  {
		    OrthogonalBasisSigmaLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
		    OrthogonalBasisSigmaRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		  }
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisSigmaLeft[i] = RealMatrix();
	  OrthogonalBasisSigmaRight[i] = RealMatrix();
	}
      if ((1 + i) <= this->PLevel)
	{	  
	  TmpMatrix.Copy(ScalarProductPhi[1 + i]);
	  TmpBasis.SetToIdentity();
#ifdef __LAPACK__
	  TmpMatrix.LapackDiagonalize(TmpDiag, TmpBasis);
#else
	  TmpMatrix.Diagonalize(TmpDiag, TmpBasis);
#endif
	  Error = 0.0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      Error = fabs(TmpDiag(n, n));
	  Error *= 1e-14;
	  if (Error < 1e-14)
	    Error = 1e-14;
	  Count  = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) < Error)
	      ++Count;
	  cout << "nbr of null vectors phi sector = " << Count << " (" << (U1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;
	  if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	    {
	      OrthogonalBasisPhiLeft[1 + i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	      OrthogonalBasisPhiRight[1 + i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	      Count = 0;
	      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
		if (fabs(TmpDiag(n, n)) > Error)
		  {
		    OrthogonalBasisPhiLeft[1 + i][Count].Copy(TmpBasis[n]);
		    OrthogonalBasisPhiRight[1 + i][Count].Copy(TmpBasis[n]);
		    if (TmpDiag(n, n) > 0)
		      {
			OrthogonalBasisPhiLeft[1 + i][Count] /=  sqrt(TmpDiag(n, n));
			OrthogonalBasisPhiRight[1 + i][Count] /=  sqrt(TmpDiag(n, n));
		      }
		    else
		      {
			OrthogonalBasisPhiLeft[1 + i][Count] /=  sqrt(-TmpDiag(n, n));
			OrthogonalBasisPhiRight[1 + i][Count] /=  -sqrt(-TmpDiag(n, n));
		      }
		    ++Count;
	      }
	    }
	  else
	    {
	      OrthogonalBasisPhiLeft[1 + i] = RealMatrix();
	      OrthogonalBasisPhiRight[1 + i] = RealMatrix();
	    }
	}
      if (i < 1)
	{
	  OrthogonalBasisPhiLeft[i] = RealMatrix();
	  OrthogonalBasisPhiRight[i] = RealMatrix();
	}

      TmpMatrix.Copy(ScalarProductEpsilon[i]);
      TmpBasis.SetToIdentity();
#ifdef __LAPACK__
      TmpMatrix.LapackDiagonalize(TmpDiag, TmpBasis);
#else
      TmpMatrix.Diagonalize(TmpDiag, TmpBasis);
#endif
      Error = 0.0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) > Error)
	  Error = fabs(TmpDiag(n, n));
      Error *= 1e-14;
      if (Error < 1e-14)
	Error = 1e-14;
      Count  = 0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) < Error)
	  ++Count;
      cout << "nbr of null vectors Epsilon sector = " << Count << " (" << (U1BosonBasis[i]->GetHilbertSpaceDimension() - Count) << " non null vectors)" << endl;
      if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	{
	  OrthogonalBasisEpsilonLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	  OrthogonalBasisEpsilonRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisEpsilonLeft[i][Count].Copy(TmpBasis[n]);
		OrthogonalBasisEpsilonRight[i][Count].Copy(TmpBasis[n]);
		if (TmpDiag(n, n) > 0)
		  {
		    OrthogonalBasisEpsilonLeft[i][Count] /=  sqrt(TmpDiag(n, n));
		    OrthogonalBasisEpsilonRight[i][Count] /=  sqrt(TmpDiag(n, n));
		  }
		else
		  {
		    OrthogonalBasisEpsilonLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
		    OrthogonalBasisEpsilonRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		  }
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisEpsilonLeft[i] = RealMatrix();
	  OrthogonalBasisEpsilonRight[i] = RealMatrix();
	}
      cout << "---------------------------------" << endl;
    }

  cout << "computing Psi matrix elements" << endl;
  LongRational Weight (WeightPsi);
  for (int j = 0; j <= this->PLevel; ++j)
    {
      for (int i = 0; i <= this->PLevel; ++i)
	{
	  RationalMatrixPsi01[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  RationalMatrixPsi10[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  RationalMatrixPsi11[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  if ((1 + j) <= this->PLevel)
	    {	  
	      RationalMatrixPsi12[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	    }
	  if ((1 + i) <= this->PLevel)
	    {	  
	      RationalMatrixPsi21[i][j] = LongRationalMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	    }
	  cout << "Levels = " <<  i << " " << j << endl;
	  if (cftDirectory != 0)
	    {
	      sprintf (TmpMatrixElementSigmaEpsilonFileName, "%s/cft_readrezayi3_matrixelement_sigmaepsilon_level_%d_%d.dat", cftDirectory, i, j);
	      sprintf (TmpMatrixElementEpsilonSigmaFileName, "%s/cft_readrezayi3_matrixelement_epsilonsigma_level_%d_%d.dat", cftDirectory, i, j);
	      sprintf (TmpMatrixElementSigmaSigmaFileName, "%s/cft_readrezayi3_matrixelement_sigmasigma_level_%d_%d.dat", cftDirectory, i, j);
	      if ((1 + i) <= this->PLevel)
		sprintf (TmpMatrixElementPhiSigmaFileName, "%s/cft_readrezayi3_matrixelement_phisigma_level_%d_%d.dat", cftDirectory, (1 + i), j);
	      if ((1 + j) <= this->PLevel)
		sprintf (TmpMatrixElementSigmaPhiFileName, "%s/cft_readrezayi3_matrixelement_sigmaphi_level_%d_%d.dat", cftDirectory, i, (1 + j));
	    }
	  if ((cftDirectory != 0) && (IsFile(TmpMatrixElementSigmaEpsilonFileName)) && 
	      (IsFile(TmpMatrixElementEpsilonSigmaFileName)) && (IsFile(TmpMatrixElementSigmaSigmaFileName)) && 
	      (((1 + i) > this->PLevel) || IsFile(TmpMatrixElementPhiSigmaFileName)) &&
	      (((1 + j) > this->PLevel) || IsFile(TmpMatrixElementSigmaPhiFileName)))
	    {
	      RationalMatrixPsi01[i][j].ReadMatrix(TmpMatrixElementEpsilonSigmaFileName);
	      RationalMatrixPsi10[i][j].ReadMatrix(TmpMatrixElementSigmaEpsilonFileName);
	      RationalMatrixPsi11[i][j].ReadMatrix(TmpMatrixElementSigmaSigmaFileName);
	      if ((1 + i) <= this->PLevel)
		RationalMatrixPsi21[i][j].ReadMatrix(TmpMatrixElementPhiSigmaFileName);
	      if ((1 + j) <= this->PLevel)
		RationalMatrixPsi12[i][j].ReadMatrix(TmpMatrixElementSigmaPhiFileName);
	    }
	  else
	    {
	      if (architecture == 0)
		{
		  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
		    for (int m = 0; m < U1BosonBasis[j]->GetHilbertSpaceDimension(); ++m)
		      {
			int PartitionLength = 0;
			U1BosonBasis[i]->GetOccupationNumber(n, TmpPartition);	    
			for (int k = 1; k <= i; ++k)
			  for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			    {
			      ++PartitionLength;
			    }
			int Position = PartitionLength;
			PartitionLength = 0;
			for (int k = 1; k <= i; ++k)
			  for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			    {
			      Partition[Position - PartitionLength - 1] = (long) k;
			      ++PartitionLength;
			    }
			U1BosonBasis[j]->GetOccupationNumber(m, TmpPartition);	    
			for (int k = 1; k <= j; ++k)
			  for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
			    {
			      Partition[PartitionLength] = -(long) k;
			      ++PartitionLength;		  
			    }
			LongRational Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightEpsilon, WeightSigma, Weight,
										 RationalMatrixPsi01, i - 1, j, U1BosonBasis, this->TemporaryOccupationNumber);
			RationalMatrixPsi01[i][j].SetMatrixElement(n, m, Tmp);
			Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightSigma, WeightEpsilon, Weight,
								    RationalMatrixPsi10, i - 1, j, U1BosonBasis, this->TemporaryOccupationNumber);
			RationalMatrixPsi10[i][j].SetMatrixElement(n, m, Tmp);
			Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightSigma, WeightSigma, Weight,
								    RationalMatrixPsi11, i - 1, j, U1BosonBasis, this->TemporaryOccupationNumber);
			RationalMatrixPsi11[i][j].SetMatrixElement(n, m, Tmp);
			if ((1 + j) <= this->PLevel)
			  {	  
			    Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightSigma, WeightPhi, Weight,
									RationalMatrixPsi12, i - 1, j, U1BosonBasis, this->TemporaryOccupationNumber);
			    RationalMatrixPsi12[i][j].SetMatrixElement(n, m, Tmp);
			  }
			if ((1 + i) <= this->PLevel)
			  {	  
			    Tmp = this->ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightPhi, WeightSigma, Weight,
									RationalMatrixPsi21, i - 1, j, U1BosonBasis, this->TemporaryOccupationNumber);
			    RationalMatrixPsi21[i][j].SetMatrixElement(n, m, Tmp);
			  }
		      }
		  if (cftDirectory != 0)
		    {
		      RationalMatrixPsi01[i][j].WriteMatrix(TmpMatrixElementEpsilonSigmaFileName);
		      RationalMatrixPsi10[i][j].WriteMatrix(TmpMatrixElementSigmaEpsilonFileName);
		      RationalMatrixPsi11[i][j].WriteMatrix(TmpMatrixElementSigmaSigmaFileName);
		      if ((3 + i) <= this->PLevel)
			RationalMatrixPsi21[i][j].WriteMatrix(TmpMatrixElementPhiSigmaFileName);
		      if ((3 + j) <= this->PLevel)
			RationalMatrixPsi12[i][j].WriteMatrix(TmpMatrixElementSigmaPhiFileName);
		    }
		}
	      else
		{
		  if ((cftDirectory != 0) && (IsFile(TmpMatrixElementEpsilonSigmaFileName)))
		    {
		      RationalMatrixPsi01[i][j].ReadMatrix(TmpMatrixElementEpsilonSigmaFileName);
		    }
		  else
		    {
		      FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, j, CentralCharge12, 
							     WeightEpsilon, WeightSigma, Weight,
							     RationalMatrixPsi01,  i - 1, j);
		      Operation1.ApplyOperation(architecture);
		      RationalMatrixPsi01[i][j] = Operation1.GetRationalMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  RationalMatrixPsi01[i][j].WriteMatrix(TmpMatrixElementEpsilonSigmaFileName);
			}
		    }
		  if ((cftDirectory != 0) && (IsFile(TmpMatrixElementSigmaEpsilonFileName)))
		    {
		      RationalMatrixPsi10[i][j].ReadMatrix(TmpMatrixElementSigmaEpsilonFileName);
		    }
		  else
		    {
		      FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, j, CentralCharge12, 
							     WeightSigma, WeightEpsilon, Weight,
							     RationalMatrixPsi10,  i - 1, j);
		      Operation1.ApplyOperation(architecture);
		      RationalMatrixPsi10[i][j] = Operation1.GetRationalMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  RationalMatrixPsi10[i][j].WriteMatrix(TmpMatrixElementSigmaEpsilonFileName);
			}
		    }
		  if ((cftDirectory != 0) && (IsFile(TmpMatrixElementSigmaSigmaFileName)))
		    {
		      RationalMatrixPsi11[i][j].ReadMatrix(TmpMatrixElementSigmaSigmaFileName);
		    }
		  else
		    {
		      FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, j, CentralCharge12, 
							     WeightSigma, WeightSigma, Weight,
							     RationalMatrixPsi11,  i - 1, j);
		      Operation1.ApplyOperation(architecture);
		      RationalMatrixPsi11[i][j] = Operation1.GetRationalMatrixElements();
		      if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			{
			  RationalMatrixPsi11[i][j].WriteMatrix(TmpMatrixElementSigmaSigmaFileName);
			}
		    }
		  if ((1 + j) <= this->PLevel)
		    {	  
		      if ((cftDirectory != 0) && (IsFile(TmpMatrixElementSigmaPhiFileName)))
			{
			  RationalMatrixPsi12[i][j].ReadMatrix(TmpMatrixElementSigmaPhiFileName);
			}
		      else
			{
			  FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, j, CentralCharge12, 
								 WeightSigma, WeightPhi, Weight,
								 RationalMatrixPsi12,  i - 1, j);
			  Operation1.ApplyOperation(architecture);
			  RationalMatrixPsi12[i][j] = Operation1.GetRationalMatrixElements();
			  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			    {
			      RationalMatrixPsi12[i][j].WriteMatrix(TmpMatrixElementSigmaPhiFileName);
			    }
			}
		    }
		  if ((1 + i) <= this->PLevel)
		    {	  
		      if ((cftDirectory != 0) && (IsFile(TmpMatrixElementPhiSigmaFileName)))
			{
			  RationalMatrixPsi21[i][j].ReadMatrix(TmpMatrixElementPhiSigmaFileName);
			}
		      else
			{
			  FQHEMPSEvaluateCFTOperation Operation1(this, U1BosonBasis, i, j, CentralCharge12, 
								 WeightPhi, WeightSigma, Weight,
								 RationalMatrixPsi21,  i - 1, j);
			  Operation1.ApplyOperation(architecture);
			  RationalMatrixPsi21[i][j] = Operation1.GetRationalMatrixElements();
			  if ((cftDirectory != 0) && (architecture->CanWriteOnDisk() == true))
			    {
			      RationalMatrixPsi21[i][j].WriteMatrix(TmpMatrixElementPhiSigmaFileName);
			    }
			}
		    }
		}
	    }
	  MatrixPsi01[i][j] = RationalMatrixPsi01[i][j];
	  MatrixPsi01[i][j] *= sqrt(2.0 / 3.0);
	  MatrixPsi10[i][j] = RationalMatrixPsi10[i][j];
	  MatrixPsi10[i][j] *= sqrt(2.0 / 3.0);
	  MatrixPsi11[i][j] = RationalMatrixPsi11[i][j];
	  MatrixPsi11[i][j] *= 1.0 / sqrt(3.0);

	  if ((1 + j) <= this->PLevel)
	    {	
	      MatrixPsi12[i][1 + j] = RationalMatrixPsi12[i][j];
	      MatrixPsi12[i][1 + j] *= sqrt(7.0 / 2.0) / 3.0;
	    }
	  if ((1 + i) <= this->PLevel)
	    {	  
	      MatrixPsi21[1 + i][j] = RationalMatrixPsi21[i][j];
	      MatrixPsi21[1 + i][j] *= -sqrt(7.0 / 2.0) / 3.0;
	    }
	}
    }

  SparseRealMatrix* BMatrices = new SparseRealMatrix[this->NbrBMatrices];

  int** StartingIndexPerPLevel = new int* [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel = new int [this->PLevel + 1];
  this->NbrIndicesPerPLevel = new int [this->PLevel + 1];
  this->TotalStartingIndexPerPLevel[0] = 0;
  StartingIndexPerPLevel[0] = new int [1];
  StartingIndexPerPLevel[0][0] = 0;

  
  int NValueShift;
  int QValue;
  int QValueDenominator;
  QValue = 5;
  QValueDenominator = 3;
  this->NbrNValue = QValueDenominator * (2 * this->PLevel) + 5;
  NValueShift = QValueDenominator * this->PLevel;

  this->NbrNValuesPerPLevel = new int [this->PLevel + 1];
  this->NInitialValuePerPLevel = new int [this->PLevel + 1];
  this->NLastValuePerPLevel = new int [this->PLevel + 1];     
  for (int i = 0; i <= this->PLevel; ++i)
    {
      this->ComputeChargeIndexRange(i, this->NInitialValuePerPLevel[i], this->NLastValuePerPLevel[i]);
      this->NbrNValuesPerPLevel[i] =  this->NLastValuePerPLevel[i] - this->NInitialValuePerPLevel[i] + 1;
    }
     
  this->NbrIndicesPerPLevel[0] = (U1BosonBasis[0]->GetHilbertSpaceDimension() * (OrthogonalBasisSigmaLeft[0].GetNbrColumn() 
										 + OrthogonalBasisSigmaLeft[0].GetNbrColumn() 
										 + OrthogonalBasisEpsilonLeft[0].GetNbrColumn() 
										 + OrthogonalBasisPhiLeft[0].GetNbrColumn())) * this->NbrNValue;
  for (int i = 1; i <= this->PLevel; ++i)
    {
      this->TotalStartingIndexPerPLevel[i] = this->TotalStartingIndexPerPLevel[i - 1] + this->NbrIndicesPerPLevel[i - 1];
      StartingIndexPerPLevel[i] = new int [i + 1];      
      StartingIndexPerPLevel[i][0] = this->TotalStartingIndexPerPLevel[i];
      int Tmp = 0;
      int Tmp2;
      for (int j = 0; j < i; ++j)
	{
	  Tmp2 = U1BosonBasis[j]->GetHilbertSpaceDimension() * (OrthogonalBasisSigmaLeft[i - j].GetNbrColumn() 
								+ OrthogonalBasisSigmaLeft[i - j].GetNbrColumn() 
								+ OrthogonalBasisEpsilonLeft[i - j].GetNbrColumn() 
								+ OrthogonalBasisPhiLeft[i - j].GetNbrColumn()) * this->NbrNValue;
	  StartingIndexPerPLevel[i][j + 1] = Tmp2 + StartingIndexPerPLevel[i][j];
	  Tmp += Tmp2;
	}
      Tmp += U1BosonBasis[i]->GetHilbertSpaceDimension() * (OrthogonalBasisSigmaLeft[0].GetNbrColumn() 
							    + OrthogonalBasisSigmaLeft[0].GetNbrColumn()
							    + OrthogonalBasisEpsilonLeft[0].GetNbrColumn() 
							    + OrthogonalBasisPhiLeft[0].GetNbrColumn()) * this->NbrNValue;
      this->NbrIndicesPerPLevel[i] =  Tmp;
    }
  
  int MatrixSize = this->NbrIndicesPerPLevel[this->PLevel] + this->TotalStartingIndexPerPLevel[this->PLevel];
  cout << "B matrix size = " << MatrixSize << endl;
  int* TmpNbrElementPerRow = new int[MatrixSize];
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;


  // B^[0]  matrix evaluation
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisSigmaLeft = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiLeft = OrthogonalBasisPhiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilonLeft = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisSigmaRight = OrthogonalBasisSigmaRight[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilonRight = OrthogonalBasisEpsilonRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiRight = OrthogonalBasisPhiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductSigma = ScalarProductSigma[i - p];
	  RealSymmetricMatrix& TmpScalarProductPhi = ScalarProductPhi[i - p];
	  RealSymmetricMatrix& TmpScalarProductEpsilon = ScalarProductEpsilon[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      for (int j = QValueDenominator; j < this->NbrNValue; ++j)
		{
		  // Epsilon
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  ++TmpNbrElementPerRow[this->GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(), TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p])];
			}
		    }
		  // sigma1 and sigma2
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  ++TmpNbrElementPerRow[this->GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(), TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p])];
			  ++TmpNbrElementPerRow[this->GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 3, NeutralIndex1, TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(), TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p])];
			}
		    }
		}
	    }
	}
    }

  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= (i - 1); ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p - 1];
	  RealMatrix& TmpOrthogonalBasisSigmaLeft = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilonLeft = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiLeft = OrthogonalBasisPhiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiRight = OrthogonalBasisPhiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductPhi = ScalarProductPhi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      for (int j = QValueDenominator; j < this->NbrNValue; ++j)
		{
		  // Phi
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  ++TmpNbrElementPerRow[this->GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 5, NeutralIndex1, TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(), TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p])];
			}
		    }
		}
	    }
	}
    }
  BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisSigmaLeft = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiLeft = OrthogonalBasisPhiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilonLeft = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisSigmaRight = OrthogonalBasisSigmaRight[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilonRight = OrthogonalBasisEpsilonRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiRight = OrthogonalBasisPhiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductSigma = ScalarProductSigma[i - p];
	  RealSymmetricMatrix& TmpScalarProductPhi = ScalarProductPhi[i - p];
	  RealSymmetricMatrix& TmpScalarProductEpsilon = ScalarProductEpsilon[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      for (int j = QValueDenominator; j < this->NbrNValue; ++j)
		{
		  // Epsilon
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  double Tmp = 0.0;
			  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
			    {
			      double Tmp1 = 0.0;			      
			      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				{
				  Tmp1 += TmpScalarProductEpsilon(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisEpsilonRight(NeutralIndex4, NeutralIndex2);				  
				}
			      Tmp += TmpOrthogonalBasisEpsilonLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
			    }
			  if (this->CylinderFlag)
			    Tmp *= exp(-this->Kappa * this->Kappa * (WeightEpsilonNumerical +  ((double) i)
								     + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								     + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
			  BMatrices[0].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(), TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							this->GetReadRezayiK3MatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(), TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			}
		    }
		  // sigma1 and sigma2
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigmaLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  double Tmp = 0.0;
			  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
			    {
			      double Tmp1 = 0.0;			      
			      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				{
				  Tmp1 += TmpScalarProductSigma(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigmaRight(NeutralIndex4, NeutralIndex2);				  
				}
				  Tmp += TmpOrthogonalBasisSigmaLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
			    }
			  if (this->CylinderFlag)
			    Tmp *= exp(-this->Kappa * this->Kappa * (WeightSigmaNumerical +  ((double) i)
								     + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								     + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
			  BMatrices[0].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(), TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							this->GetReadRezayiK3MatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(), TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			  BMatrices[0].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 3, NeutralIndex1, TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(), TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							this->GetReadRezayiK3MatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 3, NeutralIndex2, TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(), TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			}
		    }
		}
	    }
	}
    }

  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= (i - 1); ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p - 1];
	  RealMatrix& TmpOrthogonalBasisSigmaLeft = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilonLeft = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiLeft = OrthogonalBasisPhiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhiRight = OrthogonalBasisPhiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductPhi = ScalarProductPhi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      for (int j = QValueDenominator; j < this->NbrNValue; ++j)
		{
		  // Phi
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhiLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  double Tmp = 0.0;
			  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
			    {
			      double Tmp1 = 0.0;			      
			      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				{
				  Tmp1 += TmpScalarProductPhi(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPhiRight(NeutralIndex4, NeutralIndex2);
				}
			      Tmp += TmpOrthogonalBasisPhiLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
			    }
			  if (this->CylinderFlag)
			    Tmp *= exp(-this->Kappa * this->Kappa * (WeightEpsilonNumerical +  ((double) i)
								     + (((j - QValueDenominator) - NValueShift) * ((j - QValueDenominator) - NValueShift) / (4.0 * QValue * QValueDenominator))
								     + (((j - NValueShift) * (j - NValueShift)) / (4.0 * QValue*  QValueDenominator))));
			  BMatrices[0].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 5, NeutralIndex1, TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(), TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							this->GetReadRezayiK3MatrixIndex(j, ChargedIndex, this->NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 5, NeutralIndex2, TmpOrthogonalBasisEpsilonLeft.GetNbrColumn(), TmpOrthogonalBasisSigmaLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			}
		    }
		}
	    }
	}
    }

  // B^[1]  matrix evaluation
  for (int i = 0; i < MatrixSize; ++i)
    TmpNbrElementPerRow[i] = 0;
  FactorialCoefficient Coef;
  unsigned long* Partition1 = new unsigned long [this->PLevel + 2];
  unsigned long* Partition2 = new unsigned long [this->PLevel + 2];
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		  RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		  RealMatrix& TmpMatrixPsi01 = MatrixPsi01[i - p][j - q];
		  RealMatrix& TmpMatrixPsi10 = MatrixPsi10[i - p][j - q];	
		  RealMatrix& TmpMatrixPsi11 = MatrixPsi11[i - p][j - q];
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  N2 = QValueDenominator * (j - i) + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisEpsilon1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
				{
				  ++TmpNbrElementPerRow[this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisEpsilon1.GetNbrColumn(), TmpOrthogonalBasisSigma1.GetNbrColumn(), StartingIndexPerPLevel[i][p])];
				}

			    }
			  N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisEpsilon2.GetNbrColumn(); ++NeutralIndex2)
				{
				  ++TmpNbrElementPerRow[this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisEpsilon1.GetNbrColumn(), TmpOrthogonalBasisSigma1.GetNbrColumn(), StartingIndexPerPLevel[i][p])];
				}
			    }

			  N2 = QValueDenominator * (j - i) + 1 + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
				{
				  ++TmpNbrElementPerRow[this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 3, NeutralIndex1, TmpOrthogonalBasisEpsilon1.GetNbrColumn(), TmpOrthogonalBasisSigma1.GetNbrColumn(), StartingIndexPerPLevel[i][p])];
				}
			    }
			}
		    }
		}	      
	    }
	}
    }
  
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= (i - 1); ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p - 1];
	  RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		  RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		  RealMatrix& TmpMatrixPsi21 = MatrixPsi21[i - p][j - q];	
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  N2 = QValueDenominator * (j - i) + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhi1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
				{
				  ++TmpNbrElementPerRow[this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 5, NeutralIndex1, TmpOrthogonalBasisEpsilon1.GetNbrColumn(), TmpOrthogonalBasisSigma1.GetNbrColumn(), StartingIndexPerPLevel[i][p])];
				}

			    }
			}
		    }
		}	      
	    }
	}
    }
  
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= (j - 1); ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q - 1];
		  RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		  RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		  RealMatrix& TmpMatrixPsi12 = MatrixPsi12[i - p][j - q];
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhi2.GetNbrColumn(); ++NeutralIndex2)
				{
				  ++TmpNbrElementPerRow[this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisEpsilon1.GetNbrColumn(), TmpOrthogonalBasisSigma1.GetNbrColumn(), StartingIndexPerPLevel[i][p])];
				}
			    }
			}
		    }
		}	      
	    }
	}
    }

  BMatrices[1] = SparseRealMatrix(MatrixSize, MatrixSize, TmpNbrElementPerRow);
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		  RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		  RealMatrix& TmpMatrixPsi01 = MatrixPsi01[i - p][j - q];
		  RealMatrix& TmpMatrixPsi10 = MatrixPsi10[i - p][j - q];	
		  RealMatrix& TmpMatrixPsi11 = MatrixPsi11[i - p][j - q];
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  N2 = QValueDenominator * (j - i) + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisEpsilon1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisEpsilon1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  if (this->CylinderFlag)
				    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightEpsilonNumerical + WeightSigmaNumerical + ((double) (i + j))
										   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
										   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
				  BMatrices[1].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisEpsilon1.GetNbrColumn(), TmpOrthogonalBasisSigma1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								this->GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 3, NeutralIndex2, TmpOrthogonalBasisEpsilon2.GetNbrColumn(), TmpOrthogonalBasisSigma2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}

			    }
			  N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisEpsilon2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi10(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisEpsilon2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  if (this->CylinderFlag)
				    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightEpsilonNumerical + WeightSigmaNumerical + ((double) (i + j))
										   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
										   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
				  BMatrices[1].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisEpsilon1.GetNbrColumn(), TmpOrthogonalBasisSigma1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								this->GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisEpsilon2.GetNbrColumn(), TmpOrthogonalBasisSigma2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}
			    }

			  N2 = QValueDenominator * (j - i) + 1 + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi11(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  if (this->CylinderFlag)
				    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightSigmaNumerical + WeightSigmaNumerical + ((double) (i + j))
										   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
										   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
				  BMatrices[1].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 3, NeutralIndex1, TmpOrthogonalBasisEpsilon1.GetNbrColumn(), TmpOrthogonalBasisSigma1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								this->GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisEpsilon2.GetNbrColumn(), TmpOrthogonalBasisSigma2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}
			    }
			}
		    }
		}	      
	    }
	}
    }
  
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= (i - 1); ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p - 1];
	  RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		  RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		  RealMatrix& TmpMatrixPsi21 = MatrixPsi21[i - p][j - q];	
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  N2 = QValueDenominator * (j - i) + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPhi1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisSigma2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi21(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisSigma2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisPhi1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  if (this->CylinderFlag)
				    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightEpsilonNumerical + WeightSigmaNumerical + ((double) (i + j))
										   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
										   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
				  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  BMatrices[1].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 5, NeutralIndex1, TmpOrthogonalBasisEpsilon1.GetNbrColumn(), TmpOrthogonalBasisSigma1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								this->GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 3, NeutralIndex2, TmpOrthogonalBasisEpsilon2.GetNbrColumn(), TmpOrthogonalBasisSigma2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}

			    }
			}
		    }
		}	      
	    }
	}
    }
  
  for (int i = 0; i <= this->PLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisEpsilon1 = OrthogonalBasisEpsilonLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisSigma1 = OrthogonalBasisSigmaLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPhi1 = OrthogonalBasisPhiLeft[i - p];
	  for (int j = 0; j <= this->PLevel; ++j)
	    {
	      for (int q = 0; q <= (j - 1); ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q - 1];
		  RealMatrix& TmpOrthogonalBasisEpsilon2 = OrthogonalBasisEpsilonRight[j - q];
		  RealMatrix& TmpOrthogonalBasisSigma2 = OrthogonalBasisSigmaRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPhi2 = OrthogonalBasisPhiRight[j - q];
		  RealMatrix& TmpMatrixPsi12 = MatrixPsi12[i - p][j - q];
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisSigma1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPhi2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi12(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPhi2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisSigma1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  if (this->CylinderFlag)
				    Tmp *= exp(-0.5 * this->Kappa * this->Kappa * (WeightEpsilonNumerical + WeightSigmaNumerical + ((double) (i + j))
										   + ((N1 - NValueShift) * (N1 - NValueShift) / (2.0 * QValue * QValueDenominator))
										   + (((N2 - NValueShift) * (N2 - NValueShift)) / (2.0 * QValue * QValueDenominator))));
				  Tmp *= this->CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  BMatrices[1].SetMatrixElement(this->GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, this->NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisEpsilon1.GetNbrColumn(), TmpOrthogonalBasisSigma1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								this->GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, this->NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 5, NeutralIndex2, TmpOrthogonalBasisEpsilon2.GetNbrColumn(), TmpOrthogonalBasisSigma2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}
			    }
			}
		    }
		}	      
	    }
	}
    }

  for (int i = 0; i < this->NbrBMatrices; ++i)
    {
      this->RealBMatrices[i] = BMatrices[i];
    }
  delete[] BMatrices;

  delete[] ScalarProductSigma;
  delete[] ScalarProductPhi;
  delete[] ScalarProductEpsilon;
  delete[] RationalScalarProductSigma;
  delete[] RationalScalarProductPhi;
  delete[] RationalScalarProductEpsilon;
  for (int i = 0; i <= this->PLevel; ++i)
    {
      delete[] MatrixPsi01[i];
      delete[] MatrixPsi10[i];
      delete[] MatrixPsi21[i];
      delete[] MatrixPsi12[i];
      delete[] MatrixPsi11[i];
      delete[] RationalMatrixPsi01[i];
      delete[] RationalMatrixPsi10[i];
      delete[] RationalMatrixPsi21[i];
      delete[] RationalMatrixPsi12[i];
      delete[] RationalMatrixPsi11[i];
      delete U1BosonBasis[i];
    }
  delete[] TmpNbrElementPerRow;
  delete[] U1BosonBasis;
  delete[] MatrixPsi01;
  delete[] MatrixPsi10;
  delete[] MatrixPsi21;
  delete[] MatrixPsi12;
  delete[] MatrixPsi11;
  delete[] RationalMatrixPsi01;
  delete[] RationalMatrixPsi10;
  delete[] RationalMatrixPsi21;
  delete[] RationalMatrixPsi12;
  delete[] RationalMatrixPsi11;
  delete[] OrthogonalBasisSigmaLeft;
  delete[] OrthogonalBasisEpsilonLeft;
  delete[] OrthogonalBasisPhiLeft;
  delete[] OrthogonalBasisSigmaRight;
  delete[] OrthogonalBasisEpsilonRight;
  delete[] OrthogonalBasisPhiRight;
}

