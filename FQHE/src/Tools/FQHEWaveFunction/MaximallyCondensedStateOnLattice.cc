#include "MaximallyCondensedStateOnLattice.h"

#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/MultipleComplexScalarProductOperation.h"

#include "Operator/ParticleOnLatticeOneBodyOperator.h"

#include "MathTools/RandomNumber/NumRecRandomGenerator.h"

#include "Tools/NewUnconstrainedOptimizsation.h"

#include <cstdlib>
#include <sys/time.h>
#include <iostream>
using std::cout;
using std::endl;


// constructor for contact interactions on a square lattice
//
// nbrStates = number of quantum states
// states = state vectors to superpose
// space = Hilbert-space of states
// lx = Lx dimension
// ly = Ly dimension
// sublattices = number of sublattices
// randomGenerator = external random number generator
MaximallyCondensedStateOnLattice::MaximallyCondensedStateOnLattice(AbstractArchitecture *architecture, int nbrStates, ComplexVector *states, ParticleOnLattice* space, int lx, int ly, int sublattices, AbstractRandomNumberGenerator *randomGenerator)
{
  this->Architecture=architecture;
  this->NbrVectors = nbrStates;
  this->Vectors = states;
  this->Space = space;
  // lattice dimensions
  this->Lx = lx;
  this->Ly = ly;
  this->NbrSubLattices = sublattices;
  this->LastMaximumEV = 0.0;
  this->SphereParametrization= NSphereParameters(nbrStates,true);
  this->VariationalParameters.Resize(this->SphereParametrization.GetNbrParameters());

  this->DensityMatrixDimension = Lx * Ly * NbrSubLattices;

  // storage for density matrices
  this->DiagonalDensityMatrices = new HermitianMatrix[NbrVectors];
  this->NbrOffDiagonal = (this->NbrVectors * (this->NbrVectors - 1));
  this->OffDiagonalDensityMatrices = new ComplexMatrix[NbrOffDiagonal];
  for (int i=0; i<NbrVectors; ++i)
    this->DiagonalDensityMatrices[i].Resize(DensityMatrixDimension,DensityMatrixDimension);
  for (int i=0; i<NbrOffDiagonal; ++i)
    this->OffDiagonalDensityMatrices[i].Resize(DensityMatrixDimension,DensityMatrixDimension);
  this->CurrentDensityMatrix.Resize(DensityMatrixDimension, DensityMatrixDimension);
  this->CurrentHermitianMatrix.Resize(DensityMatrixDimension, DensityMatrixDimension);
  
  // matrices as temporary space for calculations
  this->M.Resize(DensityMatrixDimension,DensityMatrixDimension);
  this->Q.Resize(DensityMatrixDimension,DensityMatrixDimension); 

  if (randomGenerator!=NULL)
    {
      this->RandomNumbers = randomGenerator;
      ExternalGenerator=true;
    }
  else
    {
      timeval RandomTime;
      gettimeofday (&(RandomTime), 0);
      this->RandomNumbers = new NumRecRandomGenerator(RandomTime.tv_sec);
      ExternalGenerator=false;
    }
}

  
 

// destructor
//
MaximallyCondensedStateOnLattice::~MaximallyCondensedStateOnLattice()
{
  if (NbrVectors>0)
    {
      delete [] DiagonalDensityMatrices;
      delete [] OffDiagonalDensityMatrices;
    }
  if (ExternalGenerator==false)
    delete this->RandomNumbers;
}

// get the parameters of the Many-Body state that was last calculated
// return = state
ComplexVector & MaximallyCondensedStateOnLattice::GetVariationalParameters()
{
  return this->SphereParametrization.GetComplexCoordinates();
}

// set trial parameters
void MaximallyCondensedStateOnLattice::SetVariationalParameters(RealVector &variationalParameters)
{
  this->SphereParametrization.SetParameters(&(variationalParameters[0]));
}

// get the wavefunction corresponding to the current parameters
// return = complex vector of local amplitudes and phases
ComplexVector MaximallyCondensedStateOnLattice::GetWaveFunction()
{
  ComplexVector TmpParameters = this->SphereParametrization.GetComplexCoordinates();
  ComplexVector Result;
  Result.Copy(Vectors[0],TmpParameters[0]);
  for (int i=1; i<NbrVectors; ++i)
    Result.AddLinearCombination(TmpParameters[i],Vectors[i]);
  return ComplexVector(Result);
}
  
  
// optimize wavefunction starting from present settings of VariationalParameters
// tolerance = final tolerance on the variational parameters
// maxIter = maximal number of function evaluations
//
double MaximallyCondensedStateOnLattice::Optimize(double tolerance, int maxIter)
{
  double InitialStepSize=1.0;
  int EffectiveNbrVariationalParameters = SphereParametrization.GetNbrParameters();
  cout << "Starting Optimization ";
  this->NbrEvaluations=0;
  int NbrPoints = 2 * EffectiveNbrVariationalParameters + 1;
  int rnf;
  double Result;
  double *Work = new double[(NbrPoints+13)*(NbrPoints+EffectiveNbrVariationalParameters)
			    + 3*EffectiveNbrVariationalParameters*(EffectiveNbrVariationalParameters+3)/2 + 12];
  // passing parameter vector to optimizer as vector indexed from 1, not 0:
  double *x = &(this->VariationalParameters[0]);
  double (MaximallyCondensedStateOnLattice::*TargetFunction)(int, double*)=&MaximallyCondensedStateOnLattice::EvaluateCondensateFraction;
  MaximallyCondensedStateOnLattice *TargetObject=this;
  Result = NewUOA::newuoa(EffectiveNbrVariationalParameters, NbrPoints, x, InitialStepSize,
			  tolerance, &rnf, maxIter, Work, TargetObject, TargetFunction);
  cout << endl << "total: "<<NbrEvaluations<< " evaluations"<<endl;
  delete [] Work;
  return Result;
}

double MaximallyCondensedStateOnLattice::SimplexOptimize(double targetSize, int maxIter, double initialStep)
{
  cout << "Attention: MaximallyCondensedStateOnLattice::SimplexOptimize not implemented"<<endl;
  return 0.0;
}



// target function for optimizer routine:
double MaximallyCondensedStateOnLattice::EvaluateCondensateFraction(int nbrParameters, double *x)
{
  for (int i=0; i<this->SphereParametrization.GetNbrParameters(); ++i)
    if (this->VariationalParameters[i]!=x[i])
      this->VariationalParameters[i]=x[i];
  this->SetVariationalParameters(this->VariationalParameters);
  Complex TmpC;
  Complex TmpC2;
  this->CurrentDensityMatrix.ResizeAndClean(DensityMatrixDimension,DensityMatrixDimension);
  this->ResultingParameters = SphereParametrization.GetComplexCoordinates();
  for (int n=0; n<NbrVectors; ++n)
    {
      for (int m=0; m<n; ++m)
	CurrentDensityMatrix.AddLinearCombination(this->ResultingParameters[n]*Conj(this->ResultingParameters[m]),
					   OffDiagonalDensityMatrices[(NbrVectors-1)*n+m]);
      CurrentDensityMatrix.AddLinearCombination(SqrNorm(this->ResultingParameters[n]),DiagonalDensityMatrices[n]);
      for (int m=n+1; m<NbrVectors; ++m)
	CurrentDensityMatrix.AddLinearCombination(this->ResultingParameters[n]*Conj(this->ResultingParameters[m]),
					   OffDiagonalDensityMatrices[(NbrVectors-1)*n+m-1]);
    }
  for (int i=0; i<DensityMatrixDimension; ++i)
    {
      CurrentDensityMatrix.GetMatrixElement(i,i,TmpC);
      if (fabs(TmpC.Im)>1e-13)
	{
	  cout << "Error: Imaginary part on diagonal in element ("<<i<<","<<i<<"): "<< TmpC <<" !"<<endl;
	  exit(1);
	}
      CurrentHermitianMatrix.SetMatrixElement(i,i,TmpC.Re);
      for (int j=i+1; j<DensityMatrixDimension; ++j)
	{
	  CurrentDensityMatrix.GetMatrixElement(i,j,TmpC);
	  CurrentDensityMatrix.GetMatrixElement(j,i,TmpC2);
	  if (Norm(TmpC-Conj(TmpC2))>1e-13)
	    {
	      cout << "Error: Matrix not hermitian in elements ("<<i<<","<<j<<"): "<< TmpC << " vs "<< TmpC2 <<" !"<<endl;
	      exit(1);
	    }
	  CurrentHermitianMatrix.SetMatrixElement(i,j,TmpC);
	}
    }
  CurrentHermitianMatrix.Diagonalize(M, Q, 1e-12, 1000);
  this->LastMaximumEV = M[DensityMatrixDimension-1];
  return this->LastMaximumEV;
}

// evaluate all interaction factors
//   
void MaximallyCondensedStateOnLattice::EvaluateDensityMatrices()
{
  ComplexVector TargetVector;
  Complex Tmp;
  Complex *ScalarProducts = new Complex[NbrVectors];
  int CreationIndex, AnnihilationIndex;
  TargetVector.Resize(Vectors[0].GetVectorDimension());
  ParticleOnLatticeOneBodyOperator *DensityOperator= new ParticleOnLatticeOneBodyOperator(this->Space);
  for (int n=0; n<NbrVectors; ++n)
    {
      for (int CreationX=0; CreationX<this->Lx; ++CreationX)
	for (int CreationY=0; CreationY<this->Ly; ++CreationY)
	  for (int CreationSub=0; CreationSub<this->NbrSubLattices; ++CreationSub)
	    {
	      CreationIndex = Space->EncodeQuantumNumber(CreationX, CreationY, CreationSub, Tmp);	
	      for (int AnnihilationX=0; AnnihilationX<this->Lx; ++AnnihilationX)
		for (int AnnihilationY=0; AnnihilationY<this->Ly; ++AnnihilationY)
		  for (int AnnihilationSub=0; AnnihilationSub<this->NbrSubLattices; ++AnnihilationSub)
		    {
		      AnnihilationIndex = Space->EncodeQuantumNumber(AnnihilationX, AnnihilationY, AnnihilationSub, Tmp);
		      DensityOperator->SetCreationAnnihilationIndex(CreationIndex,AnnihilationIndex);
		      // calculate possible matrix elements in subspace of vectors
		      VectorOperatorMultiplyOperation Operation(DensityOperator,&(Vectors[n]),&TargetVector);
		      Operation.ApplyOperation(this->Architecture);
		      MultipleComplexScalarProductOperation Operation2(&TargetVector, Vectors, NbrVectors, ScalarProducts);
		      Operation2.ApplyOperation(this->Architecture);
		      for (int m=0; m<n; ++m)
			OffDiagonalDensityMatrices[(NbrVectors-1)*n+m].
			  SetMatrixElement(CreationIndex, AnnihilationIndex, Conj(ScalarProducts[m]));
		      DiagonalDensityMatrices[n].SetMatrixElement(CreationIndex, AnnihilationIndex, Conj(ScalarProducts[n]));
		      for (int m=n+1; m<NbrVectors; ++m)
			OffDiagonalDensityMatrices[(NbrVectors-1)*n+m-1].
			  SetMatrixElement(CreationIndex, AnnihilationIndex, Conj(ScalarProducts[m]));
		    }
	    }
    }
  delete [] ScalarProducts;
  delete DensityOperator;
}
