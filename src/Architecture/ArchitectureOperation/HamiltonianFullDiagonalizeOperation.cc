////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of hamiltonian full diagonalization operation           //
//                                                                            //
//                        last modification : 06/01/2012                      //
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
#include "Architecture/ArchitectureOperation/HamiltonianFullDiagonalizeOperation.h"
#include "Hamiltonian/AbstractHamiltonian.h"
#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"
#include "Matrix/RealMatrix.h"

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>


using std::cout;
using std::endl;


#ifdef __SCALAPACK__

// binding to the BLACS SL_init function
//
extern "C" void FORTRAN_NAME(sl_init)(const int* context, const int* nprow, const int* npcol);

// binding to the BLACS BLACS_GRIDINFO function
//
extern "C" void FORTRAN_NAME(blacs_gridinfo)(const int* context, const int* nprow, const int* npcol, const int* myrow, const int* mycol);

// binding to the BLACS BLACS_GRIDEXIT function
//
extern "C" void FORTRAN_NAME(blacs_gridexit)(const int* context); 

// binding to the BLACS ZGEBS2D function
//
extern "C" void FORTRAN_NAME(zgebs2d)(const int* context, const char* scope, const char* top, const int* m, const int* n, const doublecomplex* matrix, const int* lda);

// binding to the BLACS ZGEBR2D function
//
extern "C" void FORTRAN_NAME(zgebr2d)(const int* context, const char* scope, const char* top, const int* m, const int* n, const doublecomplex* matrix, const int* lda, const int* myrow, const int* mycol);

// binding to the SCALAPACK DESCINIT function
//
extern "C" void FORTRAN_NAME(descinit)(const int* desc, const int* m, const int* n, const int* mblock, const int* nblock, 
				       const int* rsrc, const int* csrc, const int* context, const int* mxllda, const int* info); 

// binding to the SCALAPACK NUMROC function
//
extern "C" int FORTRAN_NAME(numroc)(const int* n, const int* nblock, const int* iproc, const int* isrcproc, const int* nprocs); 

// binding to the SCALAPACK PDLAMCH function
//
extern "C" double FORTRAN_NAME(pdlamch)(const int* context, const char* cmach);

// binding to the SCALAPACK PZELSET function
//
extern "C" double FORTRAN_NAME(pzelset) (const doublecomplex* localMatrix, const int* rowIndex, const int* columnIndex, const int* desc, const doublecomplex* element);

// binding to the SCALAPACK PDELSET function
//
extern "C" double FORTRAN_NAME(pdelset) (const double* localMatrix, const int* rowIndex, const int* columnIndex, const int* desc, const double* element);


// binding to the SCALAPACK PZHETRD function
//
// extern "C" void FORTRAN_NAME(pzhetrd) (const char* uplo, const int* dimension, const doublecomplex* matrix, const int* localStartingRowIndex, , const int* localStartingColumnIndex, const int* desc, const doublecomplex* diagonal, const doublecomplex* offdiagonal, const doublecomplex* tau, const int* work, const double* lwork, const int* info);

// binding to the SCALAPACK PZHEEVX function
//
extern "C" void FORTRAN_NAME(pzheevx) (const char* jobz, const char* range, const char* uplo, 
				       const int* dimension, const doublecomplex* matrix, 
				       const int* localStartingRowIndex, const int* localStartingColumnIndex, const int* desc, 
				       const double* spectrumLowerBound, const double* spectrumUpperBound, 
				       const int* spectrumLowerIndex, const int* spectrumUpperIndex, 
				       const double* absoluteErrorTolerance, const int* nbrFoundEigenvalues,
				       const int* nbrFoundEigenstates, const double* eigenvalues, 
				       const double* orthogonalizationFactor, const doublecomplex* eigenstates, 
				       const int* localRowEigenstateIndex, const int* localColumnEigenstateIndex, const int* descEigenstateMatrix, 
				       const doublecomplex* work, const int* lwork, 
				       const double* rwork, const int* lrwork, 
				       const int* iwork, const int* liwork, 
				       const int* iFail, const int* iCluster, const double* gap, const int* info);

// binding to the SCALAPACK PZHEEV function
//
extern "C" void FORTRAN_NAME(pzheev) (const char* jobz, const char* uplo, 
				      const int* dimension, const doublecomplex* matrix, 
				      const int* localStartingRowIndex, const int* localStartingColumnIndex, const int* desc, 
				      const double* eigenvalues, const doublecomplex* eigenstates, 
				      const int* localRowEigenstateIndex, const int* localColumnEigenstateIndex, const int* descEigenstateMatrix, 
				      const doublecomplex* work, const int* lwork, 
				      const double* rwork, const int* lrwork, 
				      const int* info);


// binding to the SCALAPACK PDSYEVX function
//
extern "C" void FORTRAN_NAME(pdsyevx) (const char* jobz, const char* range, const char* uplo, 
				       const int* dimension, const double* matrix, 
				       const int* localStartingRowIndex, const int* localStartingColumnIndex, const int* desc, 
				       const double* spectrumLowerBound, const double* spectrumUpperBound, 
				       const int* spectrumLowerIndex, const int* spectrumUpperIndex, 
				       const double* absoluteErrorTolerance, const int* nbrFoundEigenvalues,
				       const int* nbrFoundEigenstates, const double* eigenvalues, 
				       const double* orthogonalizationFactor, const double* eigenstates, 
				       const int* localRowEigenstateIndex, const int* localColumnEigenstateIndex, const int* descEigenstateMatrix, 
				       const double* work, const int* lwork, 
				       const int* iwork, const int* liwork, 
				       const int* iFail, const int* iCluster, const double* gap, const int* info);


// binding to the SCALAPACK PDSYEV function
//
extern "C" void FORTRAN_NAME(pdsyev) (const char* jobz, const char* uplo, 
				      const int* dimension, const double* matrix, 
				      const int* localStartingRowIndex, const int* localStartingColumnIndex, const int* desc, 
				      const double* eigenvalues, const double* eigenstates, 
				      const int* localRowEigenstateIndex, const int* localColumnEigenstateIndex, const int* descEigenstateMatrix, 
				      const double* work, const int* lwork, 
				      const int* info);


#endif


// constructor 
//
// hamiltonian = pointer to the hamiltonian to use
// complexFlag = true if the hamiltonian is complex
// eigenstateFlag = true if the eigenstates have to be computed
// nbrEigenstates = number of eigenstates that have to be computed (<=0 if all eigenstates have to be computed)

HamiltonianFullDiagonalizeOperation::HamiltonianFullDiagonalizeOperation (AbstractHamiltonian* hamiltonian, bool complexFlag, bool eigenstateFlag, int nbrEigenstates)
{
  this->Hamiltonian = hamiltonian;
  this->OperationType = AbstractArchitectureOperation::HamiltonianFullDiagonalize;
  this->ComplexFlag = complexFlag;
  this->EigenstateFlag = eigenstateFlag;
  this->NbrEigenstates = nbrEigenstates;
  if (this->NbrEigenstates <= 0)
    this->NbrEigenstates = this->Hamiltonian->GetHilbertSpaceDimension();
}

// copy constructor 
//
// operation = reference on operation to copy

HamiltonianFullDiagonalizeOperation::HamiltonianFullDiagonalizeOperation(const HamiltonianFullDiagonalizeOperation& operation)
{
  this->Hamiltonian = operation.Hamiltonian;
  this->OperationType = AbstractArchitectureOperation::HamiltonianFullDiagonalize;
  this->ComplexFlag = operation.ComplexFlag;
  this->EigenstateFlag = operation.EigenstateFlag;
  this->NbrEigenstates = operation.NbrEigenstates;
}
  
// constructor from a master node information
//
// hamiltonian = pointer to the hamiltonian to use
// architecture = pointer to the distributed architecture to use for communications

HamiltonianFullDiagonalizeOperation::HamiltonianFullDiagonalizeOperation(AbstractHamiltonian* hamiltonian, SimpleMPIArchitecture* architecture)
{
  this->Hamiltonian = hamiltonian;
  this->OperationType = AbstractArchitectureOperation::HamiltonianFullDiagonalize;
  int TmpFlag = 0;
  architecture->BroadcastToSlaves(TmpFlag);
  if (TmpFlag == 1)
    this->ComplexFlag = true;
  else
    this->ComplexFlag = false;
  TmpFlag = 0;
  architecture->BroadcastToSlaves(TmpFlag);
  if (TmpFlag == 1)
    this->EigenstateFlag = true;
  else
    this->EigenstateFlag = false;
  TmpFlag = 0;
  architecture->BroadcastToSlaves(this->NbrEigenstates);
}
  
// destructor
//

HamiltonianFullDiagonalizeOperation::~HamiltonianFullDiagonalizeOperation()
{
}
  
// clone operation
//
// return value = pointer to cloned operation

AbstractArchitectureOperation* HamiltonianFullDiagonalizeOperation::Clone()
{
  return new HamiltonianFullDiagonalizeOperation (*this);
}
  
// apply operation (architecture independent)
//
// return value = true if no error occurs

bool HamiltonianFullDiagonalizeOperation::RawApplyOperation()
{
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  gettimeofday (&(TotalStartingTime2), 0);

  gettimeofday (&(TotalEndingTime2), 0);
  this->ExecutionTime = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  return true;
}

// apply operation for SimpleMPI architecture
//
// architecture = pointer to the architecture
// return value = true if no error occurs

bool HamiltonianFullDiagonalizeOperation::ArchitectureDependentApplyOperation(SimpleMPIArchitecture* architecture)
{
#ifdef __SCALAPACK__
  if (architecture->IsMasterNode())
    {
      if (architecture->RequestOperation(this->OperationType) == false)
	{
	  return false;
	}
       
      int TmpFlag = 0;
      if (this->ComplexFlag == true)
	TmpFlag = 1;
      architecture->BroadcastToSlaves(TmpFlag);
      TmpFlag = 0;
      if (this->EigenstateFlag == true)
	TmpFlag = 1;
      architecture->BroadcastToSlaves(TmpFlag);
      architecture->BroadcastToSlaves(this->NbrEigenstates);
    }


  int Context;
  int NbrNodePerColumn = architecture->GetNbrNodes();
  int NbrNodePerRow = (int) sqrt((double) NbrNodePerColumn);
  while ((NbrNodePerRow > 1) && ((NbrNodePerColumn % NbrNodePerRow) != 0))
    --NbrNodePerRow;
  NbrNodePerColumn /= NbrNodePerRow;
  FORTRAN_NAME(sl_init) (&Context, &NbrNodePerRow, &NbrNodePerColumn);
  
  int LocalNodeRow;
  int LocalNodeColumn;
  FORTRAN_NAME(blacs_gridinfo) (&Context, &NbrNodePerRow, &NbrNodePerColumn, &LocalNodeRow, &LocalNodeColumn);
  cout << "NbrNodePerRow=" << NbrNodePerRow << " NbrNodePerColumn=" << NbrNodePerColumn 
       << "   LocalNodeRow=" << LocalNodeRow << "   LocalNodeColumn=" << LocalNodeColumn << endl; 
  
  int* Desc = new int[9];
  int TmpGlobalNbrRow = this->Hamiltonian->GetHilbertSpaceDimension();
  int TmpGlobalNbrColumn = this->Hamiltonian->GetHilbertSpaceDimension();
  int NbrRowPerBlock = 64;
  int NbrColumnPerBlock = 64;
  int Information = 0;
  int TmpZero = 0;
  int LocalLeadingDimensionRow = FORTRAN_NAME(numroc) (&TmpGlobalNbrRow, &NbrRowPerBlock, &LocalNodeRow, &TmpZero, &NbrNodePerRow);
  int LocalLeadingDimensionColumn = FORTRAN_NAME(numroc) (&TmpGlobalNbrColumn, &NbrColumnPerBlock, &LocalNodeColumn, &TmpZero, &NbrNodePerColumn);
  FORTRAN_NAME(descinit) (Desc, &TmpGlobalNbrRow, &TmpGlobalNbrColumn, &NbrRowPerBlock, &NbrColumnPerBlock,
			  &TmpZero, &TmpZero, &Context, &LocalLeadingDimensionRow, &Information);
  
  
  int* DescEingenstateMatrix = new int[9];
  FORTRAN_NAME(descinit) (DescEingenstateMatrix, &TmpGlobalNbrRow, &TmpGlobalNbrColumn, &NbrRowPerBlock, &NbrColumnPerBlock,
			  &TmpZero, &TmpZero, &Context, &LocalLeadingDimensionRow, &Information);
   
  
  cout << LocalNodeRow << "," << LocalNodeColumn << "  LocalLeadingDimensionRow = " << LocalLeadingDimensionRow << endl;
  
  const char* DoublePrecisionMachineParameterIndex = "U";
  double UnderflowThreshold = FORTRAN_NAME(pdlamch) (&Context, DoublePrecisionMachineParameterIndex);  

  if (this->ComplexFlag == true)
    {
      HermitianMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension(), true);
      this->Hamiltonian->GetHamiltonian(HRep);
      doublecomplex* LocalScalapackMatrix = new doublecomplex[LocalLeadingDimensionRow * LocalLeadingDimensionColumn];
      
      Complex Tmp;
      doublecomplex TmpElement;
      for (int j = 1; j <= TmpGlobalNbrRow; ++j)
	{
	  for (int i = 1; i <= TmpGlobalNbrRow; ++i)
	    {
	      HRep.GetMatrixElement(i - 1, j - 1, Tmp);
	      TmpElement.r = Tmp.Re;
	      TmpElement.i = Tmp.Im;
	      FORTRAN_NAME(pzelset) (LocalScalapackMatrix, &i, &j, Desc, &TmpElement);
	    }
	}
      
      Information = 0; 
      const char* JobZ = "N";
      if (this->EigenstateFlag == true)
	JobZ = "V";
      const char* Range = "A";
      const char* UpperLower = "U";
      int LocalStartingRowIndex = 1;
      int LocalStartingColumnIndex = 1;
      double SpectrumLowerBound = 0.0;
      double SpectrumUpperBound = 0.0;
      int SpectrumLowerIndex = 0;
      int SpectrumUpperIndex = 0;
      double AbsoluteErrorTolerance = UnderflowThreshold;
      int NbrFoundEigenvalues = 0;
      int NbrFoundEigenstates = 0;
      double* Eigenvalues = new double [TmpGlobalNbrRow];
      double OrthogonalizationFactor = -1.0;
      doublecomplex* Eigenstates = 0;
      if (this->EigenstateFlag == true)
	Eigenstates = new doublecomplex [LocalLeadingDimensionRow * LocalLeadingDimensionColumn];
      int LocalRowEigenstateIndex = 1;
      int LocalColumnEigenstateIndex = 1;
      int* IFail = new int[TmpGlobalNbrRow];
      int* ICluster = new int [2 * NbrNodePerRow * NbrNodePerColumn];
      double* Gap = new double[NbrNodePerRow * NbrNodePerColumn];
      doublecomplex* ScalapackWorkingArea = new doublecomplex[10];
      int ScalapackWorkingAreaSize = -1;
      double* ScalapackRWorkingArea = new  double[10];
      int ScalapackRWorkingAreaSize = -1; 
      int* ScalapackIWorkingArea = new int[10];
      int ScalapackIWorkingAreaSize = -1;
  
      for (int i = 0; i < TmpGlobalNbrRow; ++i)
	Eigenvalues[i] = 0.0;
      
//       FORTRAN_NAME(pzheevx)(JobZ, Range, UpperLower, 
// 			    &TmpGlobalNbrRow, LocalScalapackMatrix, 
// 			    &LocalStartingRowIndex, &LocalStartingColumnIndex, Desc, 
// 			    &SpectrumLowerBound, &SpectrumUpperBound, &SpectrumLowerIndex, &SpectrumUpperIndex,
// 			    &AbsoluteErrorTolerance, &NbrFoundEigenvalues, 
// 			    &NbrFoundEigenstates, Eigenvalues,
// 			    &OrthogonalizationFactor, Eigenstates, 
// 			    &LocalRowEigenstateIndex, &LocalColumnEigenstateIndex, DescEingenstateMatrix,
// 			    ScalapackWorkingArea, &ScalapackWorkingAreaSize, 
// 			    ScalapackRWorkingArea, &ScalapackRWorkingAreaSize, 
// 			    ScalapackIWorkingArea, &ScalapackIWorkingAreaSize, 
// 			    IFail, ICluster, Gap, &Information);  
      
//       ScalapackWorkingAreaSize = (int) ScalapackWorkingArea[0].r;
//       ScalapackRWorkingAreaSize = (int) ScalapackRWorkingArea[0];
//       ScalapackIWorkingAreaSize = ScalapackIWorkingArea[0];
//       delete[] ScalapackWorkingArea;
//       delete[] ScalapackRWorkingArea;
//       delete[] ScalapackIWorkingArea;
//       ScalapackWorkingArea = new doublecomplex[ScalapackWorkingAreaSize];
//       ScalapackRWorkingArea = new double [ScalapackRWorkingAreaSize];
//       ScalapackIWorkingArea = new int [ScalapackIWorkingAreaSize];
      
//       FORTRAN_NAME(pzheevx)(JobZ, Range, UpperLower, 
// 			    &TmpGlobalNbrRow, LocalScalapackMatrix, 
// 			    &LocalStartingRowIndex, &LocalStartingColumnIndex, Desc, 
// 			    &SpectrumLowerBound, &SpectrumUpperBound, &SpectrumLowerIndex, &SpectrumUpperIndex,
// 			    &AbsoluteErrorTolerance, &NbrFoundEigenvalues, 
// 			    &NbrFoundEigenstates, Eigenvalues,
// 			    &OrthogonalizationFactor, Eigenstates, 
// 			    &LocalRowEigenstateIndex, &LocalColumnEigenstateIndex, DescEingenstateMatrix,
// 			    ScalapackWorkingArea, &ScalapackWorkingAreaSize, 
// 			    ScalapackRWorkingArea, &ScalapackRWorkingAreaSize, 
// 			    ScalapackIWorkingArea, &ScalapackIWorkingAreaSize, 
// 			    IFail, ICluster, Gap, &Information); 


      if (architecture->IsMasterNode())
	{
	  cout << "starting diagonalization" << endl;
	}
      timeval TotalStartingTime;
      gettimeofday (&TotalStartingTime, 0);

      FORTRAN_NAME(pzheev)(JobZ, UpperLower, 
			   &TmpGlobalNbrRow, LocalScalapackMatrix, 
			   &LocalStartingRowIndex, &LocalStartingColumnIndex, Desc, 
			   Eigenvalues,
			   Eigenstates, 
			   &LocalRowEigenstateIndex, &LocalColumnEigenstateIndex, DescEingenstateMatrix,
			   ScalapackWorkingArea, &ScalapackWorkingAreaSize, 
			   ScalapackRWorkingArea, &ScalapackRWorkingAreaSize, 
			   &Information);  
      
      ScalapackWorkingAreaSize = (int) ScalapackWorkingArea[0].r;
      ScalapackRWorkingAreaSize = (int) ScalapackRWorkingArea[0];
      delete[] ScalapackWorkingArea;
      delete[] ScalapackRWorkingArea;
      ScalapackWorkingArea = new doublecomplex[ScalapackWorkingAreaSize];
      ScalapackRWorkingArea = new double [ScalapackRWorkingAreaSize];
      
      FORTRAN_NAME(pzheev)(JobZ, UpperLower, 
			   &TmpGlobalNbrRow, LocalScalapackMatrix, 
			   &LocalStartingRowIndex, &LocalStartingColumnIndex, Desc, 
			   Eigenvalues,
			   Eigenstates, 
			   &LocalRowEigenstateIndex, &LocalColumnEigenstateIndex, DescEingenstateMatrix,
			   ScalapackWorkingArea, &ScalapackWorkingAreaSize, 
			   ScalapackRWorkingArea, &ScalapackRWorkingAreaSize, 
			   &Information);  

      timeval TotalEndingTime;
      gettimeofday (&TotalEndingTime, 0);
      
      if (architecture->IsMasterNode())
	{
	  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) + 
			(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));		      
	  cout << "diagonalization done in " << Dt << "s" << endl;
	}

      NbrFoundEigenvalues = TmpGlobalNbrRow;
      
      if (architecture->IsMasterNode())
 	{
 	  this->DiagonalizedMatrix = RealDiagonalMatrix (Eigenvalues, this->Hamiltonian->GetHilbertSpaceDimension());
	  if (this->EigenstateFlag == true)
	    {
	      this->ComplexEigenstates = ComplexMatrix (this->Hamiltonian->GetHilbertSpaceDimension(), this->NbrEigenstates, true);
	    }
 	}
      else
	{
	  if (this->EigenstateFlag == true)
	    {
	    }
	  delete[] Eigenvalues;
	}
      if (this->EigenstateFlag == true)
	delete[] Eigenstates;
      delete[] ScalapackWorkingArea;
      delete[] ScalapackRWorkingArea;
      delete[] LocalScalapackMatrix;
    }
  else
    {
      RealSymmetricMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension(), true);
      this->Hamiltonian->GetHamiltonian(HRep);
      
      double* LocalScalapackMatrix = new double[LocalLeadingDimensionRow * LocalLeadingDimensionColumn];
      
      double Tmp;
      for (int j = 1; j <= TmpGlobalNbrRow; ++j)
	{
	  for (int i = 1; i <= TmpGlobalNbrRow; ++i)
	    {
	      HRep.GetMatrixElement(i - 1, j - 1, Tmp);
	      FORTRAN_NAME(pdelset) (LocalScalapackMatrix, &i, &j, Desc, &Tmp);
	    }
	}

      Information = 0; 
      const char* JobZ = "N";
      if (this->EigenstateFlag == true)
	JobZ = "V";

      const char* Range = "A";
      const char* UpperLower = "U";
      int LocalStartingRowIndex = 1;
      int LocalStartingColumnIndex = 1;
      double SpectrumLowerBound = 0.0;
      double SpectrumUpperBound = 0.0;
      int SpectrumLowerIndex = 0;
      int SpectrumUpperIndex = 0;
      double AbsoluteErrorTolerance = UnderflowThreshold;
      int NbrFoundEigenvalues = 0;
      int NbrFoundEigenstates = 0;
      double* Eigenvalues = new double [TmpGlobalNbrRow];
      double OrthogonalizationFactor = -1.0;
      double* Eigenstates = 0;
      if (this->EigenstateFlag == true)
	Eigenstates = new double [LocalLeadingDimensionRow * LocalLeadingDimensionColumn];
      int LocalRowEigenstateIndex = 1;
      int LocalColumnEigenstateIndex = 1;
      int* IFail = new int[TmpGlobalNbrRow];
      int* ICluster = new int [2 * NbrNodePerRow * NbrNodePerColumn];
      double* Gap = new double[NbrNodePerRow * NbrNodePerColumn];
      double* ScalapackWorkingArea = new double[10];
      int ScalapackWorkingAreaSize = -1;
      int* ScalapackIWorkingArea = new int[10];
      int ScalapackIWorkingAreaSize = -1;
  
      for (int i = 0; i < TmpGlobalNbrRow; ++i)
	Eigenvalues[i] = 0.0;
      
//       FORTRAN_NAME(pdsyevx)(JobZ, Range, UpperLower, 
// 			    &TmpGlobalNbrRow, LocalScalapackMatrix, 
// 			    &LocalStartingRowIndex, &LocalStartingColumnIndex, Desc, 
// 			    &SpectrumLowerBound, &SpectrumUpperBound, &SpectrumLowerIndex, &SpectrumUpperIndex,
// 			    &AbsoluteErrorTolerance, &NbrFoundEigenvalues, 
// 			    &NbrFoundEigenstates, Eigenvalues,
// 			    &OrthogonalizationFactor, Eigenstates, 
// 			    &LocalRowEigenstateIndex, &LocalColumnEigenstateIndex, DescEingenstateMatrix,
// 			    ScalapackWorkingArea, &ScalapackWorkingAreaSize, 
// 			    ScalapackIWorkingArea, &ScalapackIWorkingAreaSize, 
// 			    IFail, ICluster, Gap, &Information);  
      
//       ScalapackWorkingAreaSize = (int) ScalapackWorkingArea[0];
//       ScalapackIWorkingAreaSize = ScalapackIWorkingArea[0];
//       delete[] ScalapackWorkingArea;
//       delete[] ScalapackIWorkingArea;
//       ScalapackWorkingArea = new double[ScalapackWorkingAreaSize];
//       ScalapackIWorkingArea = new int [ScalapackIWorkingAreaSize];

//       cout << "ScalapackWorkingAreaSize = " << ScalapackWorkingAreaSize << "  ScalapackIWorkingAreaSize = " << ScalapackIWorkingAreaSize << endl;
 
//       FORTRAN_NAME(pdsyevx)(JobZ, Range, UpperLower, 
// 			    &TmpGlobalNbrRow, LocalScalapackMatrix, 
// 			    &LocalStartingRowIndex, &LocalStartingColumnIndex, Desc, 
// 			    &SpectrumLowerBound, &SpectrumUpperBound, &SpectrumLowerIndex, &SpectrumUpperIndex,
// 			    &AbsoluteErrorTolerance, &NbrFoundEigenvalues, 
// 			    &NbrFoundEigenstates, Eigenvalues,
// 			    &OrthogonalizationFactor, Eigenstates, 
// 			    &LocalRowEigenstateIndex, &LocalColumnEigenstateIndex, DescEingenstateMatrix,
// 			    ScalapackWorkingArea, &ScalapackWorkingAreaSize, 
// 			    ScalapackIWorkingArea, &ScalapackIWorkingAreaSize, 
// 			    IFail, ICluster, Gap, &Information); 
      
      FORTRAN_NAME(pdsyev)(JobZ, UpperLower, 
			   &TmpGlobalNbrRow, LocalScalapackMatrix, 
			   &LocalStartingRowIndex, &LocalStartingColumnIndex, Desc, 
			   Eigenvalues,
			   Eigenstates, 
			   &LocalRowEigenstateIndex, &LocalColumnEigenstateIndex, DescEingenstateMatrix,
			   ScalapackWorkingArea, &ScalapackWorkingAreaSize, 
			   &Information); 
      
      ScalapackWorkingAreaSize = (int) ScalapackWorkingArea[0];
      delete[] ScalapackWorkingArea;
      ScalapackWorkingArea = new double[ScalapackWorkingAreaSize];

      cout << "ScalapackWorkingAreaSize = " << ScalapackWorkingAreaSize << "  ScalapackIWorkingAreaSize = " << ScalapackIWorkingAreaSize << endl;
 
      FORTRAN_NAME(pdsyev)(JobZ, UpperLower, 
			   &TmpGlobalNbrRow, LocalScalapackMatrix, 
			   &LocalStartingRowIndex, &LocalStartingColumnIndex, Desc, 
			   Eigenvalues,
			   Eigenstates, 
			   &LocalRowEigenstateIndex, &LocalColumnEigenstateIndex, DescEingenstateMatrix,
			   ScalapackWorkingArea, &ScalapackWorkingAreaSize, 
			   &Information); 

      NbrFoundEigenvalues = TmpGlobalNbrRow;
      NbrFoundEigenstates = TmpGlobalNbrRow;

      cout << "Information = " << Information << endl;
      if (architecture->IsMasterNode())
	{
	  cout << "NbrFoundEigenvalues = " << NbrFoundEigenvalues << endl;
	  for (int i = 0; i < NbrFoundEigenvalues; ++i)
	    {
	      cout << i << " : " << Eigenvalues[i] << endl;
	    }      
// 	  for (int i = 0; i < NbrFoundEigenstates; ++i)
// 	    {
// 	      cout << "eigenstate " << i << " : ";
// 	      for (int j = 0; j < LocalLeadingDimensionRow; ++j)
// 		cout << Eigenstates[j + (i * LocalLeadingDimensionRow)] << " ";
// 	      cout << endl;
// 	    }
	  cout << "------------------------" << endl;
	  RealDiagonalMatrix TmpDiag (this->Hamiltonian->GetHilbertSpaceDimension());
	  RealMatrix Q(this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
	  HRep.LapackDiagonalize(TmpDiag, Q);
	  for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
	    cout << j << " : " << TmpDiag[j] << endl;
// 	  for (int i = 0; i < this->Hamiltonian->GetHilbertSpaceDimension(); ++i)
// 	    {
// 	      cout << "eigenstate " << i << " : ";
// 	      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension(); ++j)
// 		cout << Q[i][j] << " ";
// 	      cout << endl;
// 	    }
	}
    }

  return true;
#else
  return this->RawApplyOperation();
#endif
}

  
