////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2003 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of basic  Arnoldi algorithm                     //
//                         for non symmetric matrices                         //
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


#include "LanczosAlgorithm/BasicArnoldiAlgorithm.h"
#include "Vector/ComplexVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"

#include <stdlib.h>
#include <iostream>


using std::cout;
using std::endl;


// default constructor
//
// architecture = architecture to use for matrix operations
// nbrEigenvalue = number of wanted eigenvalues
// maxIter = an approximation of maximal number of iteration
// strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 

BasicArnoldiAlgorithm::BasicArnoldiAlgorithm(AbstractArchitecture* architecture, int nbrEigenvalue, int maxIter, 
											 bool strongConvergence) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->MaximumNumberIteration = maxIter;
  this->NbrEigenvalue = nbrEigenvalue;
  this->ArnoldiVectors = new ComplexVector [this->MaximumNumberIteration];
  this->TemporaryCoefficients = new Complex [this->MaximumNumberIteration];
  if (maxIter > 0)
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
      this->ComplexDiagonalizedMatrix = ComplexDiagonalMatrix(this->MaximumNumberIteration, true);
    }
  else
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->ComplexDiagonalizedMatrix = ComplexDiagonalMatrix();
    }
  this->Architecture = architecture;
  this->Flag.Initialize();
  this->StrongConvergenceFlag = strongConvergence;
  this->PreviousLastWantedEigenvalue = 0.0;
  this->PreviousWantedEigenvalues = new double [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    this->PreviousWantedEigenvalues[i] = 0.0;
  this->EigenvaluePrecision = MACHINE_PRECISION;
  this->EigenvectorPrecision = 0.0;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

BasicArnoldiAlgorithm::BasicArnoldiAlgorithm(const BasicArnoldiAlgorithm& algorithm) 
{
  this->Index = algorithm.Index;
  this->MaximumNumberIteration = algorithm.MaximumNumberIteration;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->ArnoldiVectors = new ComplexVector [this->MaximumNumberIteration];
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Flag = algorithm.Flag;
  this->Architecture = algorithm.Architecture;
  this->NbrEigenvalue = algorithm.NbrEigenvalue;
  this->PreviousLastWantedEigenvalue = algorithm.PreviousLastWantedEigenvalue;
  this->EigenvaluePrecision = algorithm.EigenvaluePrecision;
  this->EigenvectorPrecision = algorithm.EigenvectorPrecision;
  this->StrongConvergenceFlag = algorithm.StrongConvergenceFlag;
  this->PreviousWantedEigenvalues = new double [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    this->PreviousWantedEigenvalues[i] = 0.0;
  this->TemporaryCoefficients = algorithm.TemporaryCoefficients;
}

// destructor
//

BasicArnoldiAlgorithm::~BasicArnoldiAlgorithm() 
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->ArnoldiVectors;
      delete[]  this->TemporaryCoefficients;
    }
  delete[] this->PreviousWantedEigenvalues;
}

// initialize Lanczos algorithm with a random vector
//

void BasicArnoldiAlgorithm::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->ArnoldiVectors[0] = ComplexVector (Dimension);
  this->ArnoldiVectors[1] = ComplexVector (Dimension);
  this->ArnoldiVectors[2] = ComplexVector (Dimension);
  for (int i = 0; i < Dimension; i++)
    {
      this->ArnoldiVectors[0].Re(i) = (rand() - 32767) * 0.5;
      this->ArnoldiVectors[0].Im(i) = (rand() - 32767) * 0.5;
    }
  this->ArnoldiVectors[0] /= this->ArnoldiVectors[0].Norm();
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void BasicArnoldiAlgorithm::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->ArnoldiVectors[0] = vector;
  this->ArnoldiVectors[1] = ComplexVector (Dimension);
  this->ArnoldiVectors[2] = ComplexVector (Dimension);
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}

// get last produced vector
//
// return value = reference on lest produced vector

Vector& BasicArnoldiAlgorithm::GetGroundState()
{
  return (*(this->GetEigenstates(1)));
}

// get the n first eigenstates
//
// nbrEigenstates = number of needed eigenstates
// return value = array containing the eigenstates

Vector* BasicArnoldiAlgorithm::GetEigenstates(int nbrEigenstates)
{
  ComplexVector* Eigenstates = new ComplexVector [nbrEigenstates];
  ComplexMatrix TmpEigenvector (this->ReducedMatrix.GetNbrRow(), this->ReducedMatrix.GetNbrRow(), true);
  for (int i = 0; i < this->TridiagonalizedMatrix.GetNbrRow(); ++i)
    TmpEigenvector(i, i) = 1.0;

  ComplexDiagonalMatrix  SortedDiagonalizedMatrix (this->ReducedMatrix.GetNbrRow());
  this->TemporaryReducedMatrix.Copy(this->ReducedMatrix);
#ifdef __LAPACK__
  ComplexDiagonalMatrix TmpDiag (SortedDiagonalizedMatrix.GetNbrColumn());
  this->TemporaryReducedMatrix.LapackDiagonalize(TmpDiag, TmpEigenvector);
  for (int i = 0; i < SortedDiagonalizedMatrix.GetNbrColumn(); ++i)
    SortedDiagonalizedMatrix[i] = TmpDiag[i];
#else
  cout << "lapack is required for BasicArnoldiAlgorithm" << endl;
#endif
  SortedDiagonalizedMatrix.SortMatrixUpOrder(TmpEigenvector);

//   Complex* TmpCoefficents = new Complex [this->BlockSize];
//   for (int i = 0; i < nbrEigenstates; ++i)
//     {
//       Eigenstates[i].Copy(this->InitialStates[0], TmpEigenvector(0, i));
//       for (int j = 1; j < this->BlockSize; ++j)
// 	TmpCoefficents[j - 1] = TmpEigenvector(j, i);	  
//       AddRealLinearCombinationOperation Operation (&(Eigenstates[i]), &(this->InitialStates[1]), this->BlockSize - 1,  TmpCoefficents);
//       Operation.ApplyOperation(this->Architecture);
//     }
//   for (int i = 0; i < this->BlockSize; ++i)
//     this->LanczosVectors[i].Copy(this->InitialStates[i]);
//   MultipleVectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, this->LanczosVectors, &(this->LanczosVectors[this->BlockSize]), this->BlockSize);
//   Operation.ApplyOperation(this->Architecture);
//   for (int i = 0; i < this->BlockSize; ++i)
//     {
//       MultipleRealScalarProductOperation Operation2 (&(this->LanczosVectors[i + this->BlockSize]), this->LanczosVectors,   
// 						     i + 1, this->TemporaryCoefficients);
//       Operation2.ApplyOperation(this->Architecture);
//       for (int j = 0; j <= i; ++j)
// 	{
// 	  this->ReducedMatrix(i, j) = this->TemporaryCoefficients[j];
// 	}
//     }
//   for (int i = 0; i < this->BlockSize; ++i)
//     {
//       for (int j = 0; j < this->BlockSize; ++j)
// 	this->TemporaryCoefficients[j] = -this->ReducedMatrix(i, j);
//       AddRealLinearCombinationOperation Operation2 (&(this->LanczosVectors[this->BlockSize + i]), this->LanczosVectors, 
// 						    this->BlockSize, this->TemporaryCoefficients);
//       Operation2.ApplyOperation(this->Architecture);
//     }
  
//   this->ReorthogonalizeVectors(&(this->LanczosVectors[this->BlockSize]), this->BlockSize, this->ReducedMatrix, 0, this->BlockSize);
//   for (int k = 0; k < nbrEigenstates; ++k)
//     {
//       for (int j = 0; j < this->BlockSize; ++j)
// 	TmpCoefficents[j] = TmpEigenvector(j, k);	  
//       AddRealLinearCombinationOperation Operation5 (&(Eigenstates[k]), &(this->LanczosVectors[this->BlockSize]), this->BlockSize,  TmpCoefficents);
//       Operation5.ApplyOperation(this->Architecture);
//     }       
//   for (int i = 1; i < this->Index; ++i)
//     {
//       this->ReducedMatrix.Resize((i + 1), (i + 1);
//       int NewVectorPosition = i * this->BlockSize;
//       MultipleVectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &(this->LanczosVectors[this->BlockSize]), 
// 							     &(this->LanczosVectors[2 * this->BlockSize]), this->BlockSize);
//       Operation2.ApplyOperation(this->Architecture);
//       for (int k = 0; k < this->BlockSize; ++k)
// 	{
// 	  MultipleRealScalarProductOperation Operation (&(this->LanczosVectors[k + 2]), 
// 							&(this->LanczosVectors[this->BlockSize]),   
// 							k + 1, this->TemporaryCoefficients);
// 	  Operation.ApplyOperation(this->Architecture);
// 	  for (int j = 0; j <= k; ++j)
// 	    {
// 	      this->ReducedMatrix(NewVectorPosition + k, NewVectorPosition + j) = this->TemporaryCoefficients[j];
// 	    }
// 	}
//       int Lim = (i - 1) * this->BlockSize;
//       for (int j = 0; j < this->BlockSize; ++j)
// 	{
// 	  for (int k = j; k < (2 * this->BlockSize); ++k)
// 	    this->TemporaryCoefficients[k - j] = -this->ReducedMatrix(Lim + this->BlockSize + j, Lim + k);
// 	  AddRealLinearCombinationOperation Operation2 (&(this->LanczosVectors[j + (2 * this->BlockSize)]), &(this->LanczosVectors[j]), 
// 							2 * this->BlockSize - j,
// 							this->TemporaryCoefficients);	  
// 	  Operation2.ApplyOperation(this->Architecture);
// 	}
//       this->ReorthogonalizeVectors(&(this->LanczosVectors[2 * this->BlockSize]), this->BlockSize, this->ReducedMatrix, 
// 				   NewVectorPosition - this->BlockSize, NewVectorPosition);  
//       for (int k = 0; k < nbrEigenstates; ++k)
// 	{
// 	  for (int j = 0; j < this->BlockSize; ++j)
// 	    TmpCoefficents[j] = TmpEigenvector(((i + 1) * this->BlockSize) + j, k);	  
// 	  AddRealLinearCombinationOperation Operation5 (&(Eigenstates[k]), &(this->LanczosVectors[2 * this->BlockSize]), this->BlockSize,  TmpCoefficents);
// 	  Operation5.ApplyOperation(this->Architecture);
// 	}       	  
//       for (int k = 0; k < this->BlockSize; ++k)
// 	{
// 	  RealVector TmpVector = this->LanczosVectors[k];
// 	  this->LanczosVectors[k] = this->LanczosVectors[k + this->BlockSize];
// 	  this->LanczosVectors[k + this->BlockSize] = this->LanczosVectors[k + (2 * this->BlockSize)];
// 	  this->LanczosVectors[k + (2 * this->BlockSize)] = TmpVector;
// 	}
//       cout << i << "/" << this->Index  << "           \r";
//       cout.flush();
//     }
//   for (int i = 0; i < this->BlockSize; ++i)
//     Eigenstates[i] /= Eigenstates[i].Norm();
//  cout << endl;
//  delete[] TmpCoefficents;
  return Eigenstates;
}

// run current Lanczos algorithm (continue from previous results if Lanczos algorithm has already been run)
//
// nbrIter = number of iteration to do 

void BasicArnoldiAlgorithm::RunLanczosAlgorithm (int nbrIter) 
{
  int Dimension;
  if (this->Index == 0)
    {
      if (nbrIter < 2)
	nbrIter = 2;
      Dimension = nbrIter;
      this->ReducedMatrix.Resize(Dimension, Dimension);

      VectorHamiltonianMultiplyOperation Operation1 (this->Hamiltonian, &(this->ArnoldiVectors[0]), &(this->ArnoldiVectors[1]));
      Operation1.ApplyOperation(this->Architecture);
      this->ReducedMatrix[0][0] = (this->ArnoldiVectors[0] * this->ArnoldiVectors[1]);
      this->ArnoldiVectors[1].AddLinearCombination(-this->ReducedMatrix[0][0], this->ArnoldiVectors[0]);
      this->ReducedMatrix[0][1] = this->ArnoldiVectors[1].Norm(); 
      this->ArnoldiVectors[1] /= this->ReducedMatrix[0][1]; 
      VectorHamiltonianMultiplyOperation Operation2 (this->Hamiltonian, &(this->ArnoldiVectors[1]), &(this->ArnoldiVectors[2]));
      Operation2.ApplyOperation(this->Architecture);
      this->ReducedMatrix[1][0] = (this->ArnoldiVectors[0] * this->ArnoldiVectors[2]);
      this->ReducedMatrix[1][1] = (this->ArnoldiVectors[1] * this->ArnoldiVectors[2]);
    }
  else
    {
      Dimension = this->ReducedMatrix.GetNbrRow() + nbrIter;
      this->ReducedMatrix.Resize(Dimension, Dimension);
    }
  for (int i = this->Index + 2; i < Dimension; ++i)
    {
      cout << i << " " << Dimension << endl;
      for (int k = 0; k < i; ++k)
	{
	  this->ReducedMatrix.GetMatrixElement(k, i, this->TemporaryCoefficients[k]);
	  this->TemporaryCoefficients[k] *= -1.0;
	}
      AddComplexLinearCombinationOperation Operation2 (&(this->ArnoldiVectors[i]), this->ArnoldiVectors, i - 1, this->TemporaryCoefficients);	  
      Operation2.ApplyOperation(this->Architecture);
      double VectorNorm = this->ArnoldiVectors[i].Norm();
      this->ReducedMatrix.SetMatrixElement(i, i - 1, VectorNorm);
      if (VectorNorm < 1e-5)
	{
	  cout << "subspace !!! " << i << endl;
	}
      this->ArnoldiVectors[i] /= VectorNorm;
      this->Index++;
      this->ArnoldiVectors[i + 1] = ComplexVector(this->Hamiltonian->GetHilbertSpaceDimension());
      VectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, &(this->ArnoldiVectors[i]), &(this->ArnoldiVectors[i + 1]));
      Operation.ApplyOperation(this->Architecture);
      MultipleComplexScalarProductOperation Operation3 (&(this->ArnoldiVectors[i + 1]), this->ArnoldiVectors,
							i, this->TemporaryCoefficients);
      Operation3.ApplyOperation(this->Architecture);
      for (int j = 0; j < i; ++j)
	{
	  this->ReducedMatrix.SetMatrixElement(j, i, this->TemporaryCoefficients[j]);
	}
    }
  if (this->PreviousLastWantedEigenvalue != 0.0)
    {
      this->PreviousLastWantedEigenvalue = Norm(this->ComplexDiagonalizedMatrix[this->NbrEigenvalue - 1]);
      for (int i = 0; i < this->NbrEigenvalue; ++i)
	this->PreviousWantedEigenvalues[i] = Norm(this->ComplexDiagonalizedMatrix[i]);
      this->Diagonalize();
      this->ComplexDiagonalizedMatrix.SortMatrixUpOrder();
    }
  else
    {
      this->Diagonalize();
      this->ComplexDiagonalizedMatrix.SortMatrixUpOrder();
      this->PreviousLastWantedEigenvalue = 2.0 * Norm(this->ComplexDiagonalizedMatrix[this->NbrEigenvalue - 1]);
      for (int i = 0; i < this->NbrEigenvalue; ++i)
	this->PreviousWantedEigenvalues[i] = 2.0 * Norm(this->ComplexDiagonalizedMatrix[i]);
    }
}

  
// test if convergence has been reached
//
// return value = true if convergence has been reached

bool BasicArnoldiAlgorithm::TestConvergence ()
{
  if (this->ReducedMatrix.GetNbrRow() > this->NbrEigenvalue)
    {
      if (this->StrongConvergenceFlag == true)
	{
	  for (int i = this->NbrEigenvalue - 1; i >= 0; --i)
	    {
	      if (Norm(this->ComplexDiagonalizedMatrix[i] - this->PreviousWantedEigenvalues[i]) > 
		  (this->EigenvaluePrecision * Norm(this->ComplexDiagonalizedMatrix[i])))
		{ 
		  if (Norm(ComplexDiagonalizedMatrix[i]) > MACHINE_PRECISION)
		    return false;
		  else
		    if (Norm(this->PreviousWantedEigenvalues[i]) > MACHINE_PRECISION)
		      return false;
		}
	    }
	  return true;
	}
      else
	if ((Norm(this->ComplexDiagonalizedMatrix[this->NbrEigenvalue - 1]) - this->PreviousLastWantedEigenvalue) < 
	    (this->EigenvaluePrecision * Norm(this->ComplexDiagonalizedMatrix[this->NbrEigenvalue - 1])))
	  {
	    return true;
	  }
	else
	  return false;
    }
  return false;
}

// diagonalize tridiagonalized matrix and find ground state energy
//

void BasicArnoldiAlgorithm::Diagonalize () 
{
  int Dimension = this->ReducedMatrix.GetNbrRow();
  this->TemporaryReducedMatrix.Copy(this->ReducedMatrix);
#ifdef __LAPACK__
  ComplexDiagonalMatrix TmpDiag (this->TemporaryReducedMatrix.GetNbrColumn());
  this->TemporaryReducedMatrix.LapackDiagonalize(TmpDiag);
  this->ComplexDiagonalizedMatrix.Resize(this->TemporaryReducedMatrix.GetNbrColumn(), this->TemporaryReducedMatrix.GetNbrColumn());
  for (int i = 0; i < this->TemporaryReducedMatrix.GetNbrColumn(); ++i)
    this->ComplexDiagonalizedMatrix[i] = TmpDiag[i];
#else
  cout << "error, LAPACK is required for BasicArnoldiAlgorithm" << endl;
#endif
  this->GroundStateEnergy = Norm(this->ComplexDiagonalizedMatrix[0]);
  for (int DiagPos = 1; DiagPos < Dimension; DiagPos++)
    if (Norm(this->ComplexDiagonalizedMatrix[DiagPos]) < this->GroundStateEnergy)
      this->GroundStateEnergy = Norm(this->ComplexDiagonalizedMatrix[DiagPos]);  
  return;
}

