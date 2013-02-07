////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2003 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//               class of implicetly restarted Arnoldi algorithm              //
//                         for non symmetric matrices                         //
//                                                                            //
//                        last modification : 06/02/2013                      //
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


#include "LanczosAlgorithm/ImplicitlyRestartedArnoldiAlgorithm.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/AddRealLinearCombinationOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"
#include "Architecture/ArchitectureOperation/MultipleRealScalarProductOperation.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexDiagonalMatrix.h"
#include "Matrix/RealUpperTriangularMatrix.h"

#include <stdlib.h>
#include <iostream>


using std::cout;
using std::endl;


// default constructor
//

ImplicitlyRestartedArnoldiAlgorithm::ImplicitlyRestartedArnoldiAlgorithm()
{ 
}

// default constructor
//
// architecture = architecture to use for matrix operations
// nbrEigenvalue = number of wanted eigenvalues
// maxNbrVectors = maximum number of vectors that can be stored in memory before restarting the Arnoldi algorithm
// nbrKeptVectors = number of vectors that are kept when restarting the Arnoldi algorithm (can't be lower than nbrEigenvalue)
// maxIter = an approximation of maximal number of iteration
// highEnergy = true if the higher energy part of the spectrum has to be computed instead of the lower energy part
// leftFlag= compute left eigenvalues/eigenvectors instead of right eigenvalues/eigenvectors
// strongConvergence = flag indicating if the convergence test has to be done on the latest wanted eigenvalue (false) or all the wanted eigenvalue (true) 

ImplicitlyRestartedArnoldiAlgorithm::ImplicitlyRestartedArnoldiAlgorithm(AbstractArchitecture* architecture, int nbrEigenvalue, int maxNbrVectors, int nbrKeptVectors,
									 int maxIter, bool highEnergy, bool leftFlag, bool strongConvergence) 
{
  this->Index = 0;
  this->Hamiltonian = 0;
  this->MaximumNumberIteration = maxIter;
  this->NbrEigenvalue = nbrEigenvalue;
  this->MaxNbrVectors = maxNbrVectors;
  this->NbrKeptVectors = nbrKeptVectors;
  if (this->NbrKeptVectors < this->NbrEigenvalue)
    this->NbrKeptVectors = this->NbrEigenvalue;
  if (this->NbrKeptVectors < 2)
    this->NbrKeptVectors = 2;
  this->TemporaryCoefficients = new double [this->MaximumNumberIteration];
  this->ArnoldiVectors = new RealVector[this->MaxNbrVectors];
  this->ArnoldiVectorMatrix = RealMatrix();
  if (maxIter > 0)
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix(this->MaximumNumberIteration, true);
      this->ComplexDiagonalizedMatrix = ComplexDiagonalMatrix(this->MaximumNumberIteration, true);
      this->ReducedMatrix = RealUpperHessenbergMatrix (this->MaximumNumberIteration, true);
    }
  else
    {
      this->TridiagonalizedMatrix = RealTriDiagonalSymmetricMatrix();
      this->ComplexDiagonalizedMatrix = ComplexDiagonalMatrix();
      this->ReducedMatrix = RealUpperHessenbergMatrix();
   }
  this->Architecture = architecture;
  this->Flag.Initialize();
  this->StrongConvergenceFlag = strongConvergence;
  this->PreviousLastWantedEigenvalue = 0.0;
  this->ComplexPreviousWantedEigenvalues = new Complex [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    this->ComplexPreviousWantedEigenvalues[i] = 0.0;
  this->EigenvaluePrecision = MACHINE_PRECISION;
  this->EigenvectorPrecision = 0.0;
  this->HighEnergyFlag = highEnergy;
  this->LeftFlag = leftFlag;
}

// copy constructor
//
// algorithm = algorithm from which new one will be created

ImplicitlyRestartedArnoldiAlgorithm::ImplicitlyRestartedArnoldiAlgorithm(const ImplicitlyRestartedArnoldiAlgorithm& algorithm) 
{
  this->Index = algorithm.Index;
  this->MaximumNumberIteration = algorithm.MaximumNumberIteration;
  this->Hamiltonian = algorithm.Hamiltonian;
  this->ArnoldiVectors = new RealVector[this->MaxNbrVectors];
  this->ArnoldiVectorMatrix = algorithm.ArnoldiVectorMatrix;
  for (int i = 0; i < this->MaxNbrVectors; ++i)
    this->ArnoldiVectors[i] = this->ArnoldiVectorMatrix[i];
  this->TridiagonalizedMatrix = algorithm.TridiagonalizedMatrix;
  this->Flag = algorithm.Flag;
  this->Architecture = algorithm.Architecture;
  this->NbrEigenvalue = algorithm.NbrEigenvalue;
  this->MaxNbrVectors = algorithm.MaxNbrVectors;
  this->NbrKeptVectors = algorithm.NbrKeptVectors;
  this->PreviousLastWantedEigenvalue = algorithm.PreviousLastWantedEigenvalue;
  this->EigenvaluePrecision = algorithm.EigenvaluePrecision;
  this->EigenvectorPrecision = algorithm.EigenvectorPrecision;
  this->StrongConvergenceFlag = algorithm.StrongConvergenceFlag;
  this->ComplexDiagonalizedMatrix = algorithm.ComplexDiagonalizedMatrix;
  this->ReducedMatrix = algorithm.ReducedMatrix;
  this->ComplexPreviousWantedEigenvalues = new Complex [this->NbrEigenvalue];
  for (int i = 0; i < this->NbrEigenvalue; ++i)
    this->ComplexPreviousWantedEigenvalues[i] = 0.0;
  this->TemporaryCoefficients = algorithm.TemporaryCoefficients;
  this->HighEnergyFlag = algorithm.HighEnergyFlag;
  this->LeftFlag = algorithm.LeftFlag;
}

// destructor
//

ImplicitlyRestartedArnoldiAlgorithm::~ImplicitlyRestartedArnoldiAlgorithm() 
{
}

// initialize Lanczos algorithm with a random vector
//

void ImplicitlyRestartedArnoldiAlgorithm::InitializeLanczosAlgorithm() 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->ArnoldiVectorMatrix = RealMatrix (Dimension, this->MaxNbrVectors);
  for (int i = 0; i < this->MaxNbrVectors; ++i)
    this->ArnoldiVectors[i] = this->ArnoldiVectorMatrix[i];
  for (int i = 0; i < Dimension; i++)
    {
      this->ArnoldiVectorMatrix[0][i] = (rand() - 32767) * 0.5;
    }
  this->ArnoldiVectorMatrix[0] /= this->ArnoldiVectorMatrix[0].Norm();
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}
  
// initialize Lanczos algorithm with a given vector
//
// vector = reference to the vector used as first step vector

void ImplicitlyRestartedArnoldiAlgorithm::InitializeLanczosAlgorithm(const Vector& vector) 
{
  int Dimension = this->Hamiltonian->GetHilbertSpaceDimension();
  this->ArnoldiVectorMatrix = RealMatrix (Dimension, this->MaxNbrVectors);
  for (int i = 0; i < this->MaxNbrVectors; ++i)
    this->ArnoldiVectors[i] = this->ArnoldiVectorMatrix[i];
  this->ArnoldiVectorMatrix[0] = vector;
  this->Index = 0;
  this->TridiagonalizedMatrix.Resize(0, 0);
}

// restart the Arnoldi algorithm if needed
//

void ImplicitlyRestartedArnoldiAlgorithm::RestartAlgorithm()
{
  if (this->Index == 2)
    {
#ifdef __LAPACK__
      cout << "restarting Arnoldi algorithm" << endl;
      RealUpperTriangularMatrix TmpR(this->ReducedMatrix.GetNbrRow());
      RealMatrix TmpQ (this->ReducedMatrix.GetNbrRow(), this->ReducedMatrix.GetNbrColumn());
      Complex Determinant = 1.0;
      for (int i = 0; i < this->ComplexDiagonalizedMatrix.GetNbrColumn(); ++i)
	Determinant *= this->ComplexDiagonalizedMatrix[i];
      this->ReducedMatrix.LapackQRFactorization(TmpR, TmpQ);
      double Determinant2 = 1.0;
      for (int i = 0; i < TmpR.GetNbrColumn(); ++i)
	{
	  double Tmp;
	  TmpR.GetMatrixElement(i, i, Tmp); 
	  cout << i << " " << Tmp << endl;
	  Determinant2 *= Tmp;
	}
      cout << "det " << Determinant << " " << Determinant2 << endl;
      RealMatrix Tmp (TmpR);
      RealMatrix Tmp2 = TmpQ * Tmp;
      cout <<  this->ReducedMatrix << endl;
      cout << "R" << endl; 
      cout << TmpR << endl;
      cout << "Q" << endl;
      cout << TmpQ << endl;
      cout << "QR" << endl;
      cout << Tmp2 << endl;
      exit(0);
      this->Index = this->NbrKeptVectors;
#else
  cout << "error, LAPACK is required for ImplicitlyRestartedArnoldiAlgorithm" << endl;
#endif
    }
}
