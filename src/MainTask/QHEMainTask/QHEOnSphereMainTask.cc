////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                       class of qhe on sphere main task                     //
//                                                                            //
//                        last modification : 10/06/2004                      //
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
#include "MainTask/QHEMainTask/QHEOnSphereMainTask.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Hamiltonian/QHEHamiltonian/AbstractQHEHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithGroundState.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedBlockLanczosAlgorithm.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithGroundStateDiskStorage.h"

#include "Options/OptionManager.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Architecture/ArchitectureOperation/ArchitectureBaseOperationManager.h"

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// constructor
//  
// options = pointer to the options managers containing all running options
// space = pointer to the current Hilbert space
// hamiltonian = pointer to the current Hamiltonian
// lValue = twice the total momentum value of the system
// shift = energy shift that is applied to the hamiltonian
// outputFileName = name of the file where results have to be stored
// firstRun = flag that indicates if it the first time the main task is used
// eigenvectorFileName = prefix to add to the name of each file that will contain an eigenvector

QHEOnSphereMainTask::QHEOnSphereMainTask(OptionManager* options, AbstractHilbertSpace* space, 
					 AbstractQHEHamiltonian* hamiltonian, int lValue, double shift, char* outputFileName,
					 bool firstRun, char* eigenvectorFileName)
{
  this->OutputFileName = new char [strlen(outputFileName) + 1];
  strncpy(this->OutputFileName, outputFileName, strlen(outputFileName));
  this->OutputFileName[strlen(outputFileName)] = '\0';
  if (eigenvectorFileName == 0)
    {
      this->EigenvectorFileName = 0;
    }
  else
    {
      this->EigenvectorFileName = new char [strlen(eigenvectorFileName) + 1];
      strncpy(this->EigenvectorFileName, eigenvectorFileName, strlen(eigenvectorFileName));
      this->EigenvectorFileName[strlen(eigenvectorFileName)] = '\0';
    }
  this->Hamiltonian = hamiltonian;
  this->Space = space;
  this->LValue = lValue;
  this->EnergyShift = shift;
  this->ResumeFlag = ((BooleanOption*) (*options)["resume"])->GetBoolean();
  this->DiskFlag = ((BooleanOption*) (*options)["disk"])->GetBoolean();
  this->MaxNbrIterLanczos = ((SingleIntegerOption*) (*options)["iter-max"])->GetInteger();
  this->NbrIterLanczos = ((SingleIntegerOption*) (*options)["nbr-iter"])->GetInteger();
  this->NbrEigenvalue = ((SingleIntegerOption*) (*options)["nbr-eigen"])->GetInteger();
  if (this->NbrEigenvalue > this->Space->GetHilbertSpaceDimension())
    {
      this->NbrEigenvalue = this->Space->GetHilbertSpaceDimension();
    }
  this->FullDiagonalizationLimit = ((SingleIntegerOption*) (*options)["full-diag"])->GetInteger();
  this->BlockLanczosFlag = false;
  if ((*options)["block-lanczos"] != 0)
    {
      this->BlockLanczosFlag = ((BooleanOption*) (*options)["block-lanczos"])->GetBoolean();
    }
  this->SizeBlockLanczos = 1;
  if ((*options)["block-size"] != 0)
    {
      this->SizeBlockLanczos = ((SingleIntegerOption*) (*options)["block-size"])->GetInteger();
    }
  this->VectorMemory = ((SingleIntegerOption*) (*options)["nbr-vector"])->GetInteger();
  this->SavePrecalculationFileName = ((SingleStringOption*) (*options)["save-precalculation"])->GetString();
  this->FullReorthogonalizationFlag = ((BooleanOption*) (*options)["force-reorthogonalize"])->GetBoolean();
  this->EvaluateEigenvectors = ((BooleanOption*) (*options)["eigenstate"])->GetBoolean();
  this->EigenvectorConvergence = ((BooleanOption*) (*options)["eigenstate-convergence"])->GetBoolean();
  if ((*options)["show-itertime"] != 0)
    {
      this->ShowIterationTime = ((BooleanOption*) (*options)["show-itertime"])->GetBoolean();
    }
  else
    this->ShowIterationTime = false;
  if ((*options)["initial-vector"] != 0)
    {
      this->InitialVectorFileName = ((SingleStringOption*) (*options)["initial-vector"])->GetString();
    }
  else
    {
      this->InitialVectorFileName = 0;
    }
  if ((*options)["partial-lanczos"] != 0)
    {
      this->PartialLanczos = ((BooleanOption*) (*options)["partial-lanczos"])->GetBoolean();
    }
  else
    {
      this->PartialLanczos = false;
    }
  if ((*options)["use-lapack"] != 0)
    {
      this->LapackFlag = ((BooleanOption*) (*options)["use-lapack"])->GetBoolean();
    }
  else
    {
      this->LapackFlag = false;
    }
  if ((*options)["limit-time"] != 0)
    {
      this->MaximumAllowedTime = (((SingleIntegerOption*) (*options)["limit-time"])->GetInteger());
    }
  else
    {
      this->MaximumAllowedTime = 0;
    }
  this->FirstRun = firstRun;
}  
 
// destructor
//  

QHEOnSphereMainTask::~QHEOnSphereMainTask()
{
  delete[] this->OutputFileName;
  if (this->EigenvectorFileName != 0)
    delete[] this->EigenvectorFileName;
}
  
// set architecture binded to the task
// 
// architecture = pointer to the architecture to use

void QHEOnSphereMainTask::SetArchitecture(AbstractArchitecture* architecture)
{
  this->Architecture = architecture;
  if ((this->Architecture->GetArchitectureID() & AbstractArchitecture::WithCommunicator) != 0)
    if (this->OperationManagers.GetNbrElement() == 0)
      {
	this->OperationManagers += new ArchitectureBaseOperationManager((SimpleMPIArchitecture*) this->Architecture, this->Hamiltonian);
      }
}

// execute the main task
// 
// return value = 0 if no error occurs, else return error code

int QHEOnSphereMainTask::ExecuteMainTask()
{
  ofstream File;
  if (this->FirstRun == true)
    {
      File.open(this->OutputFileName, ios::binary | ios::out);
      this->FirstRun = false;
    }
  else
    {
      File.open(this->OutputFileName, ios::binary | ios::out | ios::app);
    }
  File.precision(14);
  cout.precision(14);
  cout << "----------------------------------------------------------------" << endl;
  cout << " LzTotal = " << this->LValue << endl;
  cout << " Hilbert space dimension = " << this->Space->GetHilbertSpaceDimension() << endl;
  if (this->SavePrecalculationFileName != 0)
    {
      this->Hamiltonian->SavePrecalculation(this->SavePrecalculationFileName);
    }
  
  if (this->Hamiltonian->GetHilbertSpaceDimension() < this->FullDiagonalizationLimit)
    {
      RealSymmetricMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension());
      this->Hamiltonian->GetHamiltonian(HRep);
//      cout << HRep << endl;
      if (this->Hamiltonian->GetHilbertSpaceDimension() > 1)
	{
#ifdef __LAPACK__
	  if (this->LapackFlag == true)
	    {
	      RealDiagonalMatrix TmpDiag (this->Hamiltonian->GetHilbertSpaceDimension());
	      if (this->EvaluateEigenvectors == false)
		{
		  HRep.LapackDiagonalize(TmpDiag);
		}
	      else
		{
		  RealMatrix Q(this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
		  HRep.LapackDiagonalize(TmpDiag, Q);
		  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
		  RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
		  for (int j = 0; j < this->NbrEigenvalue; ++j)
		    {
		      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
		      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
		      Q[j].WriteVector(TmpVectorName);
		      cout << ((TmpEigenvector * Q[j]) - this->EnergyShift) << " " << endl;		  
		    }	      
		  cout << endl;
		  delete[] TmpVectorName;
		}
	      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
		{
		  File << (this->LValue/ 2) << " " << (TmpDiag[j] - this->EnergyShift) << endl;
		}
	    }
	  else
	    {
#endif
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (this->Hamiltonian->GetHilbertSpaceDimension());
	      if (this->EvaluateEigenvectors == false)
		{
		  HRep.Householder(TmpTriDiag, 1e-7);
		  TmpTriDiag.Diagonalize();
		  TmpTriDiag.SortMatrixUpOrder();
		}
	      else
		{
		  RealMatrix Q(this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
		  HRep.Householder(TmpTriDiag, 1e-7, Q);
		  TmpTriDiag.Diagonalize(Q);
		  TmpTriDiag.SortMatrixUpOrder(Q);
		  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
		  RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
		  for (int j = 0; j < this->NbrEigenvalue; ++j)
		    {
		      this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
		      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, j);
		      Q[j].WriteVector(TmpVectorName);
		      cout << ((TmpEigenvector * Q[j]) - this->EnergyShift) << " " << endl;		  
		    }	      
		  cout << endl;
		  delete[] TmpVectorName;
		}
	      for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
		{
		  File << (this->LValue/ 2) << " " << (TmpTriDiag.DiagonalElement(j) - this->EnergyShift) << endl;
		}
#ifdef __LAPACK__
	    }
#endif
	}
      else
	{
	  File << (this->LValue/ 2) << " " << (HRep(0, 0)  - this->EnergyShift) << endl;
	}
    }
  else
    {
      AbstractLanczosAlgorithm* Lanczos;
      if ((this->NbrEigenvalue == 1) && (this->FullReorthogonalizationFlag == false))
	{
	  if (this->DiskFlag == false)
	    if (this->EvaluateEigenvectors == true)
	      Lanczos = new BasicLanczosAlgorithmWithGroundState(this->Architecture, this->MaxNbrIterLanczos);
	    else
	      Lanczos = new BasicLanczosAlgorithm(this->Architecture, this->NbrEigenvalue, this->MaxNbrIterLanczos);
	  else
	    if (this->EvaluateEigenvectors == true)
	      Lanczos = new BasicLanczosAlgorithmWithGroundStateDiskStorage(this->Architecture, this->NbrIterLanczos, this->MaxNbrIterLanczos);
	    else
	      Lanczos = new BasicLanczosAlgorithmWithDiskStorage(this->Architecture, this->NbrEigenvalue, this->MaxNbrIterLanczos);
	}
      else
	{
	  if (this->DiskFlag == false)
	    {
	      if (this->BlockLanczosFlag == true)
		Lanczos = new FullReorthogonalizedBlockLanczosAlgorithm (this->Architecture, this->NbrEigenvalue, this->SizeBlockLanczos, this->MaxNbrIterLanczos);
	      else
		Lanczos = new FullReorthogonalizedLanczosAlgorithm (this->Architecture, this->NbrEigenvalue, this->MaxNbrIterLanczos);
	    }
	  else
	    Lanczos = new FullReorthogonalizedLanczosAlgorithmWithDiskStorage (this->Architecture, this->NbrEigenvalue, this->VectorMemory, this->MaxNbrIterLanczos);
	}
      double GroundStateEnergy;
      double Precision = 1.0;
      double PreviousLowest = 1e50;
      double Lowest = PreviousLowest;
      int CurrentNbrIterLanczos = 0;
      Lanczos->SetHamiltonian(this->Hamiltonian);
      if ((this->DiskFlag == true) && (this->ResumeFlag == true))
	Lanczos->ResumeLanczosAlgorithm();
      else
	{
	  if (this->InitialVectorFileName == 0)
	    {
	      Lanczos->InitializeLanczosAlgorithm();
	    }
	  else
	    {	      
	      RealVector InitialVector;
	      InitialVector.ReadVector(this->InitialVectorFileName);
	      Lanczos->InitializeLanczosAlgorithm(InitialVector);
	    }
	}
      cout << "Run Lanczos Algorithm" << endl;
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      timeval TotalCurrentTime;
      double Dt;
      gettimeofday (&(TotalStartingTime), 0);
      int StartTimeSecond = TotalStartingTime.tv_sec;
      if (this->ResumeFlag == false)
	{
	  Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
	  CurrentNbrIterLanczos = NbrEigenvalue + 3;
	}
      RealTriDiagonalSymmetricMatrix TmpMatrix;
      gettimeofday (&(TotalCurrentTime), 0); 
      int CurrentTimeSecond = TotalCurrentTime.tv_sec;
      while ((Lanczos->TestConvergence() == false) && (((this->DiskFlag == true) && (((this->MaximumAllowedTime == 0) && (CurrentNbrIterLanczos < this->NbrIterLanczos)) || 
										     ((this->MaximumAllowedTime > 0) && (this->MaximumAllowedTime > (CurrentTimeSecond - StartTimeSecond))))) ||
						       ((this->DiskFlag == false) && ((this->PartialLanczos == false) && (CurrentNbrIterLanczos < this->MaxNbrIterLanczos)) ||
							((this->PartialLanczos == true) && (CurrentNbrIterLanczos < this->NbrIterLanczos)))))
	{
	  if (this->BlockLanczosFlag == true)
	    CurrentNbrIterLanczos += this->SizeBlockLanczos;
	  else
	    ++CurrentNbrIterLanczos;
	  Lanczos->RunLanczosAlgorithm(1);
	  TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	  TmpMatrix.SortMatrixUpOrder();
	  Lowest = TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->EnergyShift;
	  Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	  PreviousLowest = Lowest; 
	  cout << (TmpMatrix.DiagonalElement(0) - this->EnergyShift) << " " << Lowest << " " << Precision << " ";
	  gettimeofday (&(TotalEndingTime), 0);
	  CurrentTimeSecond = TotalEndingTime.tv_sec;
	  if (this->ShowIterationTime == true)
	    {
	      Dt = (double) (TotalEndingTime.tv_sec - TotalCurrentTime.tv_sec) + 
		((TotalEndingTime.tv_usec - TotalCurrentTime.tv_usec) / 1000000.0);		      
	      cout << "(" << Dt << " s for step " << CurrentNbrIterLanczos << ")";
	      TotalCurrentTime.tv_usec = TotalEndingTime.tv_usec;
	      TotalCurrentTime.tv_sec = TotalEndingTime.tv_sec;
	    }
	  cout << endl;
	}
      if ((Lanczos->TestConvergence() == true) && (CurrentNbrIterLanczos == 0))
	{
	  TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	  TmpMatrix.SortMatrixUpOrder();
	}
      if (CurrentNbrIterLanczos >= this->MaxNbrIterLanczos)
	{
	  cout << "too much Lanczos iterations" << endl;
	  File << "too much Lanczos iterations" << endl;
	  File.close();
	  return 1;
	}
      GroundStateEnergy = Lowest;
      cout << endl;
      cout << (TmpMatrix.DiagonalElement(0) - this->EnergyShift) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	   << CurrentNbrIterLanczos << endl;
      for (int i = 0; i < this->NbrEigenvalue; ++i)
	{
	  cout << (TmpMatrix.DiagonalElement(i) - this->EnergyShift) << " ";
	  File << (this->LValue/ 2) << " " << (TmpMatrix.DiagonalElement(i) - this->EnergyShift) << endl;
	}
      cout << endl;
      if ((this->EvaluateEigenvectors == true) && (((this->DiskFlag == true) && (((this->MaximumAllowedTime == 0) && (CurrentNbrIterLanczos < this->NbrIterLanczos)) || 
										 ((this->MaximumAllowedTime > 0) && (this->MaximumAllowedTime > (CurrentTimeSecond - StartTimeSecond))))) ||
						   ((this->DiskFlag == false) && ((this->PartialLanczos == false) && (CurrentNbrIterLanczos < this->MaxNbrIterLanczos)) ||
						    ((this->PartialLanczos == true) && (CurrentNbrIterLanczos < this->NbrIterLanczos)))))
	{
	  RealVector* Eigenvectors = (RealVector*) Lanczos->GetEigenstates(this->NbrEigenvalue);
	  if (Eigenvectors != 0)
	    {
	      RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
	      if ((this->EigenvectorConvergence == true) && ((this->PartialLanczos == false) || (CurrentNbrIterLanczos <= this->NbrIterLanczos)))
		{
		  VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(Eigenvectors[this->NbrEigenvalue - 1]), &TmpEigenvector);
		  Operation1.ApplyOperation(this->Architecture);
		  double Scalar = TmpEigenvector * Eigenvectors[this->NbrEigenvalue - 1];
		  Precision = fabs((Scalar - TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1)) / TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1));
		  while (Precision > 1e-7)
		    {
		      ++CurrentNbrIterLanczos;
		      Lanczos->RunLanczosAlgorithm(1);
		      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
		      TmpMatrix.SortMatrixUpOrder();
		      Lowest = TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->EnergyShift;
		      delete[] Eigenvectors;
		      Eigenvectors = (RealVector*) Lanczos->GetEigenstates(this->NbrEigenvalue);
		      VectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &(Eigenvectors[this->NbrEigenvalue - 1]), &TmpEigenvector);
		      Operation1.ApplyOperation(this->Architecture);
		      Scalar = TmpEigenvector * Eigenvectors[this->NbrEigenvalue - 1];
		      Scalar = TmpEigenvector * Eigenvectors[this->NbrEigenvalue - 1];
		      Precision = fabs((Scalar - TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1)) / TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1));		  
		      cout << (TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->EnergyShift) << " " << (Scalar - this->EnergyShift) << " " 
			   << Precision << " ";
		      if (this->ShowIterationTime == true)
			{
			  gettimeofday (&(TotalEndingTime), 0);
			  Dt = (double) (TotalEndingTime.tv_sec - TotalCurrentTime.tv_sec) + 
			    ((TotalEndingTime.tv_usec - TotalCurrentTime.tv_usec) / 1000000.0);		      
			  cout << "(" << Dt << " s)";
			  TotalCurrentTime.tv_usec = TotalEndingTime.tv_usec;
			  TotalCurrentTime.tv_sec = TotalEndingTime.tv_sec;
			}
		      cout << endl;
		    }
		}
	      char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 32];
	      for (int i = 0; i < this->NbrEigenvalue; ++i)
		{
		  VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(Eigenvectors[i]), &TmpEigenvector);
		  Operation1.ApplyOperation(this->Architecture);
		  cout << ((TmpEigenvector * Eigenvectors[i]) - this->EnergyShift) << " ";	
		  if ((this->PartialLanczos == false) || (CurrentNbrIterLanczos < this->NbrIterLanczos))
		    {	  
		      sprintf (TmpVectorName, "%s.%d.vec", this->EigenvectorFileName, i);
		    }
		  else
		    {
		      sprintf (TmpVectorName, "%s.%d.part.vec", this->EigenvectorFileName, i);		  
		    }
		  Eigenvectors[i].WriteVector(TmpVectorName);
		}
	      cout << endl;
	      delete[] TmpVectorName;
	      delete[] Eigenvectors;
	    }
	  else
	    {
	      cout << "eigenvectors can't be computed" << endl;
	    }
	}
      gettimeofday (&(TotalEndingTime), 0);
      cout << "------------------------------------------------------------------" << endl << endl;;
      Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
      cout << "time = " << Dt << endl;
      delete Lanczos;
    }
  cout << "----------------------------------------------------------------" << endl;
  File.close(); 
  return 0;
}

