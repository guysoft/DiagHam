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

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/AbstractHilbertSpace.h"
#include "Hamiltonian/QHEHamiltonian/AbstractQHEHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithmWithDiskStorage.h"

#include "Options/OptionManager.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>


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

QHEOnSphereMainTask::QHEOnSphereMainTask(OptionManager* options, AbstractHilbertSpace* space, 
					 AbstractQHEHamiltonian* hamiltonian, int lValue, double shift, char* outputFileName)
{
  this->OutputFileName = new char [strlen(outputFileName) + 1];
  strncpy(this->OutputFileName, outputFileName, strlen(outputFileName));
  this->OutputFileName[strlen(outputFileName)] = '\0';
  this->Hamiltonian = hamiltonian;
  this->Space = space;
  this->LValue = lValue;
  this->EnergyShift = shift;
  this->ResumeFlag = ((BooleanOption*) (*options)["resume"])->GetBoolean();
  this->DiskFlag = ((BooleanOption*) (*options)["disk"])->GetBoolean();
  this->MaxNbrIterLanczos = ((SingleIntegerOption*) (*options)["iter-max"])->GetInteger();
  this->NbrIterLanczos = ((SingleIntegerOption*) (*options)["nbr-iter"])->GetInteger();
  this->NbrEigenvalue = ((SingleIntegerOption*) (*options)["nbr-eigen"])->GetInteger();
  this->FullDiagonalizationLimit = ((SingleIntegerOption*) (*options)["full-diag"])->GetInteger();
  this->VectorMemory = ((SingleIntegerOption*) (*options)["nbr-vector"])->GetInteger();
  this->SavePrecalculationFileName = ((SingleStringOption*) (*options)["save-precalculation"])->GetString();
}  
 
// destructor
//  

QHEOnSphereMainTask::~QHEOnSphereMainTask()
{
  delete[] this->OutputFileName;
}
  
// execute the main task
// 
// return value = 0 if no error occurs, else return error code

int QHEOnSphereMainTask::ExecuteMainTask()
{
  ofstream File;
  File.open(this->OutputFileName, ios::binary | ios::out | ios::app);
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
      if (this->Hamiltonian->GetHilbertSpaceDimension() > 1)
	{
	  RealTriDiagonalSymmetricMatrix TmpTriDiag (this->Hamiltonian->GetHilbertSpaceDimension());
	  HRep.Householder(TmpTriDiag, 1e-7);
	  TmpTriDiag.Diagonalize();
	  TmpTriDiag.SortMatrixUpOrder();
	  for (int j = 0; j < this->Hamiltonian->GetHilbertSpaceDimension() ; j++)
	    {
	      File << (this->LValue/ 2) << " " << (TmpTriDiag.DiagonalElement(j) - this->EnergyShift) << endl;
	    }
	}
      else
	{
	  File << (this->LValue/ 2) << " " << (HRep(0, 0)  - this->EnergyShift) << endl;
	}
    }
  else
    {
      AbstractLanczosAlgorithm* Lanczos;
      if (this->NbrEigenvalue == 1)
	{
	  if (this->DiskFlag == false)
	    Lanczos = new BasicLanczosAlgorithm(this->Architecture, this->NbrEigenvalue, this->MaxNbrIterLanczos);
	  else
	    Lanczos = new BasicLanczosAlgorithmWithDiskStorage(this->Architecture, this->NbrEigenvalue, this->MaxNbrIterLanczos);
	}
      else
	{
	  if (this->DiskFlag == false)
	    Lanczos = new FullReorthogonalizedLanczosAlgorithm (this->Architecture, this->NbrEigenvalue, this->MaxNbrIterLanczos);
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
	Lanczos->InitializeLanczosAlgorithm();
      cout << "Run Lanczos Algorithm" << endl;
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      double Dt;
      gettimeofday (&(TotalStartingTime), 0);
      if (ResumeFlag == false)
	{
	  Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
	  CurrentNbrIterLanczos = NbrEigenvalue + 3;
	}
      RealTriDiagonalSymmetricMatrix TmpMatrix;
      while ((Lanczos->TestConvergence() == false) && (((this->DiskFlag == true) && (CurrentNbrIterLanczos < this->NbrIterLanczos)) ||
						       ((this->DiskFlag == false) && (CurrentNbrIterLanczos < this->MaxNbrIterLanczos))))
	{
	  ++CurrentNbrIterLanczos;
	  Lanczos->RunLanczosAlgorithm(1);
	  TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	  TmpMatrix.SortMatrixUpOrder();
	  Lowest = TmpMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->EnergyShift;
	  Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	  PreviousLowest = Lowest; 
	  cout << (TmpMatrix.DiagonalElement(0) - this->EnergyShift) << " " << Lowest << " " << Precision << " "<< endl;
	}
      if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
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
      for (int i = 0; i <= this->NbrEigenvalue; ++i)
	{
	  cout << (TmpMatrix.DiagonalElement(i) - this->EnergyShift) << " ";
	  File << (this->LValue/ 2) << " " << (TmpMatrix.DiagonalElement(i) - this->EnergyShift) << endl;
	}
      cout << endl;
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

