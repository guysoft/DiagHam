////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                                                                            //
//           class of full reorthogonalized block Lanczos algorithm           //
//                (with full re-orthogonalization at each step)               //
//                                                                            //
//                        last modification : 10/03/2005                      //
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


#include "LanczosAlgorithm/FullReorthogonalizedBlockLanczosAlgorithm.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Matrix/RealMatrix.h"

#include <stdlib.h>
#include <iostream>


using std::cout;
using std::endl;


// default constructor
//
// architecture = architecture to use for matrix operations
// nbrEigenvalue = number of wanted eigenvalues (rounded to the upper multiple of blockSize)
// blockSize = size of the block used for the block Lanczos algorithm
// maxIter = an approximation of maximal number of iteration (rounded to the upper multiple of blockSize)

FullReorthogonalizedBlockLanczosAlgorithm::FullReorthogonalizedBlockLanczosAlgorithm(AbstractArchitecture* architecture, int nbrEigenvalue, int blockSize, int maxIter) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->BlockSize = blockSize;
  this->MaximumNumberIteration = maxIter;
  this->NbrEigenvalue = nbrEigenvalue;
  if ((this->NbrEigenvalue % this->BlockSize) != 0)
    this->NbrEigenvalue = (this->NbrEigenvalue / this->BlockSize) + this->BlockSize;
  if ((this->MaximumNumberIteration % this->BlockSize) != 0)
    this->MaximumNumberIteration = (this->MaximumNumberIteration / this->BlockSize) + this->BlockSize;
  this->LanczosVectors = new RealVector [this->MaximumNumberIteration];

  if (this->MaximumNumberIteration > 0)
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
      this->ReducedMatrix = RealSymmetricMatrix(this->MaximumNumberIteration, true);
      this->TemporaryReducedMatrix = RealSymmetricMatrix(this->MaximumNumberIteration, true);
    }
  else
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->DiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->ReducedMatrix = RealSymmetricMatrix();
      this->TemporaryReducedMatrix = RealSymmetricMatrix();
    }
  this->Architecture = architecture;
  this->Flag.Initialize();
  this->PreviousLastWantedEigenvalue = 0.0;
  this->EigenvaluePrecision = MACHINE_PRECISION;
  this->EigenvectorPrecision = 0.0;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

FullReorthogonalizedBlockLanczosAlgorithm::FullReorthogonalizedBlockLanczosAlgorithm(const FullReorthogonalizedBlockLanczosAlgorithm& algorithm) 
{
  this->LanczosVectors = algorithm.LanczosVectors;
  this->Index = algorithm.Index;

  this->MaximumNumberIteration = algorithm.MaximumNumberIteration;
  this->NbrEigenvalue = algorithm.NbrEigenvalue;
  this->BlockSize = algorithm.algorithm;

  this->Hamiltonian = algorithm.Hamiltonian;
  this->Architecture = algorithm.Architecture;

  this->PreviousLastWantedEigenvalue = algorithm.PreviousLastWantedEigenvalue;
  this->EigenvaluePrecision = algorithm.EigenvaluePrecision;
  this->EigenvectorPrecision = algorithm.EigenvectorPrecision;
 
  this->ReducedMatrix = algorithm.ReducedMatrix;
  this->TemporaryReducedMatrix = algorithm.TemporaryReducedMatrix;
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->DiagonalizedMatrix = algorithm.DiagonalizedMatrix;

  this->Flag = algorithm.Flag;
}

// destructor
//

FullReorthogonalizedBlockLanczosAlgorithm::~FullReorthogonalizedBlockLanczosAlgorithm() 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->LanczosVectors;
    }
}

// initialize Lanczos algorithm with a random vector
//

void FullReorthogonalizedBlockLanczosAlgorithm::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->LanczosVectors[0] = RealVector (Dimension);
  for (int i = 0; i < Dimension; i++)
    this->LanczosVectors[0][i] = (drand48() - 0.5) * 2.0;
  this->LanczosVectors[0] /= this->LanczosVectors[0].Norm();
  double* TmpCoef = new double [this->NbrEigenvalue];
  for (int j = 1; j < this->BlockSize; ++j)
    {
      this->LanczosVectors[j] = RealVector (Dimension);
      for (int i = 0; i < Dimension; ++i)
	this->LanczosVectors[j][i] = (drand48() - 0.5) * 2.0;
      for (int i = 0; i < j; ++i)
	TmpCoef[i] = this->LanczosVectors[j] * this->LanczosVectors[i];
      for (int i = 0; i < j; ++i)
	this->LanczosVectors[j].AddLinearCombination(TmpCoef[i], this->LanczosVectors[i]);
      double TmpNorm = this->LanczosVectors[j].Norm();
      if (TmpNorm < this->DeflationPrecision)
	--j;
      else
	this->LanczosVectors[j] /= TmpNorm;
    }
  delete[] TmpCoef;
  this->Index = 0;
  this->ReducedMatrix.Resize(0, 0);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void FullReorthogonalizedBlockLanczosAlgorithm::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->LanczosVectors[0] = vector;
  double* TmpCoef = new double [this->NbrEigenvalue];
  for (int j = 1; j < this->NbrEigenvalue; ++j)
    {
      this->LanczosVectors[j] = RealVector (Dimension);
      for (int i = 0; i < Dimension; ++i)
	this->LanczosVectors[j][i] = (drand48() - 0.5) * 2.0;
      for (int i = 0; i < j; ++i)
	TmpCoef[i] = this->LanczosVectors[j] * this->LanczosVectors[i];
      for (int i = 0; i < j; ++i)
	  this->LanczosVectors[j].AddLinearCombination(TmpCoef[i], this->LanczosVectors[i]);
      double TmpNorm = this->LanczosVectors[j].Norm();
      if (TmpNorm < this->DeflationPrecision)
	--j;
      else
	this->LanczosVectors[j] /= TmpNorm;
    }
  delete[] TmpCoef;
  this->Index = 0;
  this->ReducedMatrix.Resize(0, 0);
}

// get last produced vector
//
// return value = reference on lest produced vector

Vector& FullReorthogonalizedBlockLanczosAlgorithm::GetGroundState()
{
  RealVector* Eigenstates = new RealVector [nbrEigenstates];
  RealMatrix TmpEigenvector (this->ReducedMatrix.GetNbrRow(), this->ReducedMatrix.GetNbrRow(), true);
  for (int i = 0; i < this->ReducedMatrix.GetNbrRow(); ++i)
    TmpEigenvector(i, i) = 1.0;

  RealTriDiagonalSymmetricMatrix SortedDiagonalizedMatrix (this->ReducedMatrix.GetNbrRow());
  this->TemporaryReducedMatrix.Copy(this->ReducedMatrix);
  this->TemporaryReducedMatrix.Householder(SortedDiagonalizedMatrix, 1e-7, TmpEigenstates);
  SortedDiagonalizedMatrix.Diagonalize(TmpEigenstates);
  SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);

  double* TmpCoefficents = new double [this->ReducedMatrix.GetNbrRow()];
  for (int j = 0; j < this->ReducedMatrix.GetNbrRow(); ++j)
    TmpCoefficents[j] = TmpEigenvector(j, 0);
  this->GroundState = RealVector (this->Hamiltonian->GetHilbertSpaceDimension());
  this->GroundState.Copy(this->LanczosVectors[0], TmpEigenvector(0, 0));
  AddRealLinearCombinationOperation Operation (&(this->GroundState), &(this->LanczosVectors[1]), this->ReducedMatrix.GetNbrRow() - 1, &(TmpCoefficents[1]));
  this->Architecture->ExecuteOperation(&Operation);
  this->GroundState /= this->GroundState.Norm();
  delete[] TmpCoefficents;

  return this->GroundState;
}

// get the n first eigenstates
//
// nbrEigenstates = number of needed eigenstates
// return value = array containing the eigenstates

Vector* FullReorthogonalizedBlockLanczosAlgorithm::GetEigenstates(int nbrEigenstates)
{
  RealVector* Eigenstates = new RealVector [nbrEigenstates];
  RealMatrix TmpEigenvector (this->ReducedMatrix.GetNbrRow(), this->ReducedMatrix.GetNbrRow(), true);
  for (int i = 0; i < this->ReducedMatrix.GetNbrRow(); ++i)
    TmpEigenvector(i, i) = 1.0;

  RealTriDiagonalSymmetricMatrix SortedDiagonalizedMatrix (this->ReducedMatrix.GetNbrRow());
  this->TemporaryReducedMatrix.Copy(this->ReducedMatrix);
  this->TemporaryReducedMatrix.Householder(SortedDiagonalizedMatrix, 1e-7, TmpEigenstates);
  SortedDiagonalizedMatrix.Diagonalize(TmpEigenstates);
  SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);

  double* TmpCoefficents = new double [this->ReducedMatrix.GetNbrRow()];
  for (int i = 0; i < nbrEigenstates; ++i) 
    {
      for (int j = 0; j < this->ReducedMatrix.GetNbrRow(); ++j)
	TmpCoefficents[j] = TmpEigenvector(j, i);
      Eigenstates[i] = RealVector (this->Hamiltonian->GetHilbertSpaceDimension());
      Eigenstates[i].Copy(this->LanczosVectors[0], TmpEigenvector(0, i));
      AddRealLinearCombinationOperation Operation (&(Eigenstates[i]), &(this->LanczosVectors[1]), this->ReducedMatrix.GetNbrRow() - 1, &(TmpCoefficents[1]));
      this->Architecture->ExecuteOperation(&Operation);
      Eigenstates[i] /= Eigenstates[i].Norm();
    }
  delete[] TmpCoefficents;
  return Eigenstates;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void FullReorthogonalizedBlockLanczosAlgorithm::RunLanczosAlgorithm (int nbrIter) 
{
  int Dimension;
  if (this->Index == 0)
    {
      Dimension = this->ReducedMatrix.GetNbrRow() + nbrIter;
      if (nbrIter < 2)
	Dimension = this->ReducedMatrix.GetNbrRow() + 2;
      this->ReducedMatrix.Resize(Dimension, Dimension);

       for (int k = 0; k < this->BlockSize; ++k)
	{
	  this->LanczosVectors[k + this->BlockSize] = RealVector(this->Hamiltonian->GetHilbertSpaceDimension());
	  VectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, &(this->LanczosVectors[k]), &(this->LanczosVectors[k + this->BlockSize]));
	  this->Architecture->ExecuteOperation(&Operation);
	}

      int Lim = 2 * this->BlockSize;
      for (int i = 0; i < this->BlockSize; ++i)
	{
	  MultipleRealScalarProductOperation Operation (&(this->LanczosVectors[i + this->BlockSize]), this->LanczosVectors,   
							i + 1, this->TemporaryCoefficient);
	  this->Architecture->ExecuteOperation(&Operation);
	  for (int j = 0; j < i; ++j)
	    {
	      this->ReducedMatrix(i, j) = this->TemporaryCoefficient[j];
	    }
	}
      
      for (int i = 0; i < this->BlockSize; ++i)
	{
	  for (int j = 0; j < this->BlockSize; ++j)
	    this->TemporaryCoefficient[j] = -this->ReducedMatrix(i, j);
	  AddRealLinearCombinationOperation Operation (&(this->LanczosVectors[this->BlockSize + i]), this->LanczosVectors, this->BlockSize, this->TemporaryCoefficient);
	  this->Architecture->ExecuteOperation(&Operation);
	}
	  
      double TmpNorm;
      for (int k = this->BlockDimension[1]; k < Lim; ++k)
	{
	  TmpNorm = this->LanczosVectors[k].Norm();
	  if (TmpNorm > this->DeflationPrecision)
	    this->LanczosVectors[k] =/ TmpNorm; 
	else
	  {
	    --Lim;
	    RealVector TmpVector = this->LanczosVectors[Lim];
	    this->LanczosVectors[Lim] = this->LanczosVectors[k]; 
	    this->LanczosVectors[k] = TmpVector; 
	    --this->CurrentBlockSize;
	  }
	}

      this->BlockDimension[1] = this->CurrentBlockSize;

      this->BlockPosition[2] = this->BlockDimension[1] + this->BlockPosition[1];
      Lim = 3 * this->BlockSize;      
      for (int k =  2 * this->BlockSize; k < Lim; ++k)
	{
	  this->LanczosVectors[k + this->BlockDimension[1]] = RealVector(this->Hamiltonian->GetHilbertSpaceDimension());
	  VectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, &(this->LanczosVectors[k]), &(this->LanczosVectors[k + this->BlockDimension[1]]));
	  this->Architecture->ExecuteOperation(&Operation);
	}

      Lim = this->BlockPosition[2] + this->CurrentBlockSize;
      for (int j = 2 * this->BlockSize; j < Lim; ++j)
	{
	  for (int i = 0; i < (2 * this->BlockSize); ++i)
	    this->ReducedMatrix(i, j) = this->LanczosVectors[i] * this->LanczosVectors[j];
	}
      this->Index = 2;

    }
  else
    {
      Dimension = this->ReducedMatrix.GetNbrRow() + nbrIter;
      this->ReducedMatrix.Resize(Dimension, Dimension);
    }
  for (; nbrIter >= 0; --nbrIter)
    {

      int Lim = this->BlockPosition[this->Index] + this->CurrentBlockSize;
      for (int i = this->BlockPosition[this->Index]; i < Lim; ++i)
	{
	  int j = 0;
	  for (; j < this->BlockPosition[Index - 2]; ++j)
	    this->TemporaryCoefficient[j] = -this->ReducedMatrix(i, j);
	  for (; j < this->BlockPosition[Index]; ++j)
	    this->TemporaryCoefficient[j] = -this->ReducedMatrix(i, j);
	  AddRealLinearCombinationOperation Operation (&(this->LanczosVectors[i]), this->LanczosVectors, this->BlockPosition[Index], this->TemporaryCoefficient);
	  this->Architecture->ExecuteOperation(&Operation);
	}

      double VectorNorm = this->LanczosVectors[i].Norm();

      double TmpNorm;
      for (int k = this->BlockPosition[this->Index]; k < Lim; ++k)
	{
	  TmpNorm = this->LanczosVectors[k].Norm();
	  if (TmpNorm > this->DeflationPrecision)
	  this->LanczosVectors[k] =/ TmpNorm; 
	  else
	    {
	      --Lim;
	      RealVector TmpVector = this->LanczosVectors[Lim];
	      this->LanczosVectors[Lim] = this->LanczosVectors[k]; 
	      this->LanczosVectors[k] = TmpVector; 
	      --this->CurrentBlockSize;
	    }
	}
      this->BlockDimension[this->Index] = this->CurrentBlockSize;
      this->Index++;
      this->BlockPosition[this->Index] = this->BlockPosition[this->Index - 1] + this->CurrentBlockSize;
      this->BlockDimension[this->Index] = this->CurrentBlockSize;

      for (int k = this->BlockPosition[this->Index - 1]; k < this->BlockPosition[this->Index] ; ++k)
	{
	  this->LanczosVectors[k + this->CurrentBlockSize] = RealVector(this->Hamiltonian->GetHilbertSpaceDimension());
	  VectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, &(this->LanczosVectors[k]), &(this->LanczosVectors[k + this->CurrentBlockSize]));
	  this->Architecture->ExecuteOperation(&Operation);
	}

      Lim = this->BlockPosition[this->Index] + this->CurrentBlockSize;
      for (int j = this->BlockPosition[this->Index]; j < Lim; ++j)
	{
	  for (int i = this->BlockPosition[this->Index - 1]; i < this->BlockPosition[this->Index]; ++i)
	    this->ReducedMatrix(i, j) = this->LanczosVectors[i] * this->LanczosVectors[j];
	  for (int i = j - this->BlockPosition[2]; i < this->BlockPosition[this->Index]; ++i)
	    this->ReducedMatrix(i, j) = this->LanczosVectors[i] * this->LanczosVectors[j];	    
	}

    }
  if (this->PreviousLastWantedEigenvalue != 0.0)
    {
      this->PreviousLastWantedEigenvalue = this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1);
      this->Diagonalize();
      this->DiagonalizedMatrix.SortMatrixUpOrder();
    }
  else
    {
      this->Diagonalize();
      this->DiagonalizedMatrix.SortMatrixUpOrder();
      this->PreviousLastWantedEigenvalue = 2.0 * this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1);
    }
}

  
// diagonalize tridiagonalized matrix and find ground state energy
//

void FullReorthogonalizedBlockLanczosAlgorithm::Diagonalize () 
{
  int Dimension = this->ReducedMatrix.GetNbrRow();
  this->TemporaryReducedMatrix.Copy(this->ReducedMatrix);
  this->TemporaryReducedMatrix.Householder(this->DiagonalizedMatrix, 1e-7);
  this->DiagonalizedMatrix.Diagonalize();
  this->GroundStateEnergy = this->DiagonalizedMatrix.DiagonalElement(0);
  for (int DiagPos = 1; DiagPos < Dimension; DiagPos++)
    if (this->DiagonalizedMatrix.DiagonalElement(DiagPos) < this->GroundStateEnergy)
      this->GroundStateEnergy = this->DiagonalizedMatrix.DiagonalElement(DiagPos);  
  return;
}

// test if convergence has been reached
//
// return value = true if convergence has been reached

bool FullReorthogonalizedBlockLanczosAlgorithm::TestConvergence ()
{
  cout << this->PreviousLastWantedEigenvalue << " " << this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1) << " " << this->EigenvaluePrecision<< endl;
  if (((this->ReducedMatrix.GetNbrRow() > this->NbrEigenvalue) && 
      (fabs(this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1) - this->PreviousLastWantedEigenvalue) < 
       (this->EigenvaluePrecision * fabs(this->DiagonalizedMatrix.DiagonalElement(this->NbrEigenvalue - 1))))) || 
      (this->CurrentBlockSize <= 0))
    return true;
  else
    return false;
}

