////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Nicolas Regnault                  //
//                                                                            //
//                       class author: Gunnar Moeller                         //
//                                                                            //
//                      class of a generic real main task                     //
//                                                                            //
//                       last modification : 03/08/2009                       //
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

#include "MainTask/GenericRealMainTask.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Hamiltonian/AbstractHamiltonian.h"

#include "LanczosAlgorithm/LanczosManager.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"

#include "Options/Options.h"

#include "Architecture/ArchitectureOperation/ArchitectureBaseOperationManager.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/GenericSignalHandler.h"


#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <cstring>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;



// constructor
//  
// options = pointer to the options managers containing all running options
// space = pointer to the current Hilbert space
// hamiltonian = pointer to the current Hamiltonian
// subspaceStr = string to prepend in output file for each eigenvalue in this subspace
// subspaceLegend = legend indicating contents of subspaceStr to include in output file
// shift = energy shift that is applied to the hamiltonian
// outputFileName = name of the file where results have to be stored
// firstRun = flag that indicates if it the first time the main task is used
// eigenvectorFileName = prefix to add to the name of each file that will contain an eigenvector

GenericRealMainTask::GenericRealMainTask(OptionManager* options, AbstractHilbertSpace* space, LanczosManager* lanczos, 
					 AbstractHamiltonian* hamiltonian, char *subspaceStr, char *subspaceLegend,
					 double shift, char* outputFileName, bool firstRun, char* eigenvectorFileName)
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
  this->AlgorithmManager = lanczos;
  this->EnergyShift = shift;
  this->SubspaceStr = new char[strlen(subspaceStr)+1];
  strcpy(this->SubspaceStr,subspaceStr);
  this->SubspaceLegend = new char[strlen(subspaceLegend)+1];
  strcpy(this->SubspaceLegend,subspaceLegend);
  this->ResumeFlag = options->GetBoolean("resume");
  this->DiskFlag = options->GetBoolean("disk");
  this->MaxNbrIterLanczos = options->GetInteger("iter-max");
  this->NbrIterLanczos = options->GetInteger("nbr-iter");
  this->NbrEigenvalue = options->GetInteger("nbr-eigen");
  if (this->NbrEigenvalue > this->Space->GetHilbertSpaceDimension())
    {
      this->NbrEigenvalue = this->Space->GetHilbertSpaceDimension();
    }
  this->FullDiagonalizationLimit = options->GetInteger("full-diag");
  this->BlockLanczosFlag = false;
  if ((*options)["block-lanczos"] != 0)
    {
      this->BlockLanczosFlag = options->GetBoolean("block-lanczos");
    }
  this->SizeBlockLanczos = 1;
  if ((*options)["block-size"] != 0)
    {
      this->SizeBlockLanczos = options->GetInteger("block-size");
    }
  this->VectorMemory = options->GetInteger("nbr-vector");
  if ((*options)["save-precalculation"] != 0)
    {
      this->SavePrecalculationFileName = options->GetString("save-precalculation");
    }
  else
    this->SavePrecalculationFileName = 0;
  this->FullReorthogonalizationFlag = options->GetBoolean("force-reorthogonalize");
  this->EvaluateEigenvectors = options->GetBoolean("eigenstate");
  this->EigenvectorConvergence = options->GetBoolean("eigenstate-convergence");
  if ((*options)["show-itertime"] != 0)
    {
      this->ShowIterationTime = options->GetBoolean("show-itertime");
    }
  else
    this->ShowIterationTime = false;
  if ((*options)["initial-vector"] != 0)
    {
      this->InitialVectorFileName = options->GetString("initial-vector");
    }
  else
    {
      this->InitialVectorFileName = 0;
    }
  if ((*options)["initial-blockvectors"] != 0)
    {
      this->InitialBlockVectorFileName = options->GetString("initial-blockvectors");
    }
  else
    {
      this->InitialBlockVectorFileName = 0;
    }
  if ((*options)["partial-lanczos"] != 0)
    {
      this->PartialLanczos = options->GetBoolean("partial-lanczos");
    }
  else
    {
      this->PartialLanczos = false;
    }
  if ((*options)["use-lapack"] != 0)
    {
      this->LapackFlag = options->GetBoolean("use-lapack");
    }
  else
    {
      this->LapackFlag = false;
    }
  if ((*options)["limit-time"] != 0)
    {
      this->MaximumAllowedTime = (options->GetInteger("limit-time"));
    }
  else
    {
      this->MaximumAllowedTime = 0;
    }
  if ((((*options)["use-hilbert"]) != 0) && (options->GetString("use-hilbert") != 0))
    {
      this->ReducedHilbertSpaceDescription = options->GetString("use-hilbert");
    }
  else
    {
      this->ReducedHilbertSpaceDescription = 0;
    }
  if ((*options)["get-hvalue"] != 0)
    {
      this->ComputeEnergyFlag = options->GetBoolean("get-hvalue");
    }
  else
    {
      this->ComputeEnergyFlag = false;
    }
  if (((*options)["show-hamiltonian"] != 0) && (options->GetBoolean("show-hamiltonian") == true))
    {
      RealSymmetricMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension());
      this->Hamiltonian->GetHamiltonian(HRep);
      cout << HRep << endl;
    }  
  if (((*options)["lanczos-precision"] != 0) && (options->GetDouble("lanczos-precision") > 0))
    {
      this->LanczosPrecision = options->GetDouble("lanczos-precision");
    }
  else
    {
      this->LanczosPrecision = 0.0;
    }
  if (((*options)["fast-disk"] != 0) && (this->NbrEigenvalue == 1) && (this->EvaluateEigenvectors == true))
    {
      this->FastDiskFlag = options->GetBoolean("fast-disk");
      if ((*options)["resume-fastdisk"] != 0)
	{
	  this->ResumeFastDiskFlag = options->GetBoolean("resume-fastdisk");
	}
    }
  else
    {
      this->FastDiskFlag = false;
      this->ResumeFastDiskFlag = false;
    }
  this->PartialEigenstateFlag = 0;
  if (((*options)["partial-eigenstate"] != 0) && (this->EvaluateEigenvectors == true))
    {
      this->PartialEigenstateFlag = options->GetInteger("partial-eigenstate");
    }
  this->FirstRun = firstRun;
}
 
// destructor
//  

GenericRealMainTask::~GenericRealMainTask()
{
  delete[] this->OutputFileName;
  if (this->EigenvectorFileName != 0)
    delete[] this->EigenvectorFileName;
  delete [] this->SubspaceStr;
  delete [] this->SubspaceLegend;
}
  
// set architecture bound to the task
// 
// architecture = pointer to the architecture to use

void GenericRealMainTask::SetArchitecture(AbstractArchitecture* architecture)
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

int GenericRealMainTask::ExecuteMainTask()
{
  ofstream File;
  if (this->FirstRun == true)
    {
      File.open(this->OutputFileName, ios::binary | ios::out);
      this->FirstRun = false;
      File << "# "<<SubspaceLegend;
      File <<" E";
      if ((this->EvaluateEigenvectors == true) && (this->ComputeEnergyFlag == true))
	File << " <H>";
      File << endl;
    }
  else
    {
      File.open(this->OutputFileName, ios::binary | ios::out | ios::app);
    }
  File.precision(14);
  cout.precision(14);
  cout << "----------------------------------------------------------------" << endl;
  cout << " Hilbert space dimension = " << this->Space->GetHilbertSpaceDimension() << endl;
  if (this->SavePrecalculationFileName != 0)
    {
      this->Hamiltonian->SavePrecalculation(this->SavePrecalculationFileName);
    }
  if (this->Hamiltonian->GetHilbertSpaceDimension() == 0)
    return 0;
  if (this->Hamiltonian->GetHilbertSpaceDimension() < this->FullDiagonalizationLimit)
    {
      RealSymmetricMatrix HRep (this->Hamiltonian->GetHilbertSpaceDimension(), true);
      this->Hamiltonian->GetHamiltonian(HRep);
      if (this->Hamiltonian->GetHilbertSpaceDimension() > 1)
	{
#ifdef __LAPACK__
	  if (this->LapackFlag == true)
	    {
	      RealDiagonalMatrix TmpDiag (this->Hamiltonian->GetHilbertSpaceDimension());
	      if (this->EvaluateEigenvectors == false)
		{
		  HRep.LapackDiagonalize(TmpDiag);
		  for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
		    this->WriteResult(File, TmpDiag[j] - this->EnergyShift);
		}
	      else
		{
		  RealMatrix Q(this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
		  HRep.LapackDiagonalize(TmpDiag, Q);
		  if (this->EvaluateEigenvectors == true)
		    {
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
		      this->WriteResult(File,TmpDiag[j] - this->EnergyShift, false);
		      if (this->ComputeEnergyFlag == true)
			{
			  RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
			  this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			  File << " " << ((TmpEigenvector * Q[j]) - this->EnergyShift);
			}
		      File << endl;
		    }
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
		  for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; ++j)
		    this->WriteResult(File, TmpTriDiag.DiagonalElement(j) - this->EnergyShift);
		}
	      else
		{
		  RealMatrix Q(this->Hamiltonian->GetHilbertSpaceDimension(), this->Hamiltonian->GetHilbertSpaceDimension());
		  HRep.Householder(TmpTriDiag, 1e-7, Q);
		  TmpTriDiag.Diagonalize(Q);
		  TmpTriDiag.SortMatrixUpOrder(Q);
		  if (this->EvaluateEigenvectors == true)
		    {
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
		      this->WriteResult(File, TmpTriDiag.DiagonalElement(j) - this->EnergyShift, false);
		      if (this->ComputeEnergyFlag == true)
			{
			  RealVector TmpEigenvector(this->Hamiltonian->GetHilbertSpaceDimension());
			  this->Hamiltonian->LowLevelMultiply(Q[j], TmpEigenvector);
			  File << " " << ((TmpEigenvector * Q[j]) - this->EnergyShift);
			}
		      File << endl;
		    }
		}
#ifdef __LAPACK__
	    }
#endif
	}
      else
	{
	  this->WriteResult(File, HRep(0, 0)  - this->EnergyShift, false);
	  if (this->ComputeEnergyFlag == true)
	    File << " " << (HRep(0, 0)  - this->EnergyShift) ;
	  File << endl;	      
	  if (this->EvaluateEigenvectors)
	    {
	      char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 16];
	      RealVector TmpEigenvector(1);
	      TmpEigenvector[0]=1.0;
	      sprintf (TmpVectorName, "%s.0.vec", this->EigenvectorFileName);
	      TmpEigenvector.WriteVector(TmpVectorName);
	      delete [] TmpVectorName;
	    }
	}
    }
  else
    {
      AbstractLanczosAlgorithm* Lanczos = AlgorithmManager->GetLanczosAlgorithm(this->Architecture, this->EvaluateEigenvectors);
      if (this->LanczosPrecision != 0.0)
	Lanczos->SetEigenvaluePrecision(this->LanczosPrecision);
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
	  if (this->BlockLanczosFlag == false)
	    {
	      if (this->InitialVectorFileName == 0)
		Lanczos->InitializeLanczosAlgorithm();
	      else
		{	   
		  RealVector InitialVector;
		  InitialVector.ReadVector(this->InitialVectorFileName);
		  Lanczos->InitializeLanczosAlgorithm(InitialVector);
		}
	    }
	  else
	    {
	      if (this->InitialBlockVectorFileName == 0)
		Lanczos->InitializeLanczosAlgorithm();
	      else
		{
		  int TmpNbrInitialVectors;
		  ConfigurationParser InitialVectorDescription;
		  if (InitialVectorDescription.Parse(this->InitialBlockVectorFileName) == false)
		    {
		      InitialVectorDescription.DumpErrors(cout) << endl;
		    }
		  else
		    {
		      char** VectorFileNames;
		      if (InitialVectorDescription.GetAsStringArray("InitialVectors", ' ', VectorFileNames, TmpNbrInitialVectors) == false)
			{
			  cout << "Vectors are not defined or have a wrong value in " << this->InitialBlockVectorFileName << endl;
			}
		      else
			{
			  RealVector* InitialVectors = new RealVector[TmpNbrInitialVectors];
			  for (int i = 0; i < TmpNbrInitialVectors; ++i)
			    {
			      InitialVectors[i].ReadVector(VectorFileNames[i]);
			    }
			  Lanczos->InitializeLanczosAlgorithm(InitialVectors, TmpNbrInitialVectors);		  
			  delete[] InitialVectors;
			}
		    }
		}
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
      GenericSignalHandler Usr1Handler(SIGUSR1);
      while ((Lanczos->TestConvergence() == false) && (((this->DiskFlag == true) && (((this->MaximumAllowedTime == 0) && (CurrentNbrIterLanczos < this->NbrIterLanczos)) || 
										     ((this->MaximumAllowedTime > 0) && (this->MaximumAllowedTime > (CurrentTimeSecond - StartTimeSecond))))) ||
						       ((this->DiskFlag == false) && ((this->PartialLanczos == false) && (CurrentNbrIterLanczos < this->MaxNbrIterLanczos)) ||
							((this->PartialLanczos == true) && (CurrentNbrIterLanczos < this->NbrIterLanczos)))))
	{
	  if (this->BlockLanczosFlag == true)
	    CurrentNbrIterLanczos += this->SizeBlockLanczos;
	  else
	    ++CurrentNbrIterLanczos;
	  Usr1Handler.StartToDeferSignal();
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

	  if (Usr1Handler.HavePendingSignal())
	    {
	      cout << "Terminating Lanczos iteration on user signal"<<endl;
	      File << "# Lanczos terminated at step "<<CurrentNbrIterLanczos<<" with precision "<<Precision<<endl;
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      for (int i = 0; i < this->NbrEigenvalue; ++i)
		{
		  cout << (TmpMatrix.DiagonalElement(i) - this->EnergyShift) << " ";
		  this->WriteResult(File, TmpMatrix.DiagonalElement(i) - this->EnergyShift, true);
		}
	      cout << endl;
	    }
	  
	  if ( ((Usr1Handler.HavePendingSignal()) && (this->EvaluateEigenvectors == true))
	       || ((this->PartialEigenstateFlag > 0) && ((CurrentNbrIterLanczos % (this->PartialEigenstateFlag * this->SizeBlockLanczos)) == 0)))
	    {
	      RealVector* Eigenvectors = (RealVector*) Lanczos->GetEigenstates(this->NbrEigenvalue);
	      if (Eigenvectors != 0)
		{
		  char* TmpVectorName = new char [strlen(this->EigenvectorFileName) + 32];
		  for (int i = 0; i < this->NbrEigenvalue; ++i)
		    {
		      sprintf (TmpVectorName, "%s.%d.part.%d.vec", this->EigenvectorFileName, i, CurrentNbrIterLanczos);		  
		      Eigenvectors[i].WriteVector(TmpVectorName);
		    }
		  delete[] TmpVectorName;
		  delete[] Eigenvectors;
		}
	      else
		{
		  cout << "eigenvectors can't be computed" << endl;
		}
	    }
	  Usr1Handler.ProcessDeferredSignal();
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
	  if  (this->ComputeEnergyFlag == false)
	    WriteResult(File, TmpMatrix.DiagonalElement(i) - this->EnergyShift);
	}
      cout << endl;
      if ((this->EvaluateEigenvectors == true) && 
	  (((this->DiskFlag == true) && (((this->MaximumAllowedTime == 0) && (CurrentNbrIterLanczos < this->NbrIterLanczos)) || 
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
		  if (this->ComputeEnergyFlag == true)
		    WriteResult(File, TmpMatrix.DiagonalElement(i) - this->EnergyShift, false);
		  VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(Eigenvectors[i]), &TmpEigenvector);
		  Operation1.ApplyOperation(this->Architecture);
		  cout << ((TmpEigenvector * Eigenvectors[i]) - this->EnergyShift) << " ";	
		  if (this->EvaluateEigenvectors == true)
		    {
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
		  if (this->ComputeEnergyFlag == true)
		    {
		      File << " " << ((TmpEigenvector * Eigenvectors[i]) - this->EnergyShift);
		    }
		  if (this->ComputeEnergyFlag == true)
		    File << endl;
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
      AlgorithmManager->FreeLanczosAlgorithm();
    }
  cout << "----------------------------------------------------------------" << endl;
  File.close(); 
  return 0;
}

// add optiongroup with options related to this module to the given OptionManager
//
void GenericRealMainTask::AddOptionGroup(OptionManager *optionManager)
{
  this->AlgorithmManager->AddOptionGroup(optionManager);
}


// write a line of output to the results file
//
// file = stream to write to
// value = numerical value to be printed after columns for flux and momentum (if defined)
// terminate = indicate if line should be terminated with endl
void GenericRealMainTask::WriteResult(ofstream& file, double value, bool terminate)
{
  file << SubspaceStr <<" ";
  file << value;
  // cout <<"wrote value "<<value<<" to file"<<endl;
  if (terminate)
    file << endl;
}
