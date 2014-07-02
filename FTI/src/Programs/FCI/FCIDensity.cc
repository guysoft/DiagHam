#include "Vector/RealVector.h"
#include "Matrix/ComplexMatrix.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Operator/ParticleOnSphereDensityOperator.h"
#include "FunctionBasis/ParticleOnDiskFunctionBasis.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


// get the Hilbert space and the vector state form the input file name
//
// inputState = input file name
// nbrParticles = reference on the number of particles 
// nbrSitesX = reference on the number of site along the x direction
// nbrSitesY = reference on the number of site along the y direction
// kxMomentum = reference on the momentum along the x direction
// kyMomentum = reference on the momentum along the y direction
// statistics = reference on the statistic flag
// space = reference on the pointer to the Hilbert space
// state = reference on the state vector
// return value = true if no error occured
bool FQHECheckerboardLatticeModelDensityGetHilbertSpace(char* inputState, int& nbrParticles, int& nbrSitesX, int& nbrSitesY,
							int& kxMomentum, int& kyMomentum, bool& statistics,
							ParticleOnSphere*& space, ComplexVector& state);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHECheckerboardLatticeModelDensity" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PlotOptionGroup = new OptionGroup ("plot options");  
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += PlotOptionGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "input-states", "use a file to describe the state as a linear combination");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (overriding the one found in the vector file name if greater than 0)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0 - 0.5 * M_SQRT2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tpp", "second next nearest neighbor hoping amplitude", 0.5 * (M_SQRT2 - 1.0));
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "coefficients-only", "only compute the one body coefficients that are requested to evaluate the density profile", false);

  (*PrecalculationGroup) += new SingleStringOption  ('\n', "use-precomputed", "use precomputed matrix elements to do the plot");
  (*PlotOptionGroup) += new SingleStringOption ('\n', "output", "output file name (default output name replace the .vec extension of the input file with .rho.dat)", 0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECheckerboardLatticeModelDensity -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (Manager.GetString("input-states") == 0)
    {
      cout << "FQHECheckerboardLatticeModelDensity requires an input state" << endl;
      return -1;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey"); 
  bool Statistics = false;
  bool CoefficientOnlyFlag = Manager.GetBoolean("coefficients-only");
  char* OutputName = Manager.GetString("output");
  double NNHoping = Manager.GetDouble("t1");
  double NextNNHoping = Manager.GetDouble("t2");
  double SecondNextNNHoping = Manager.GetDouble("tpp");
  double MuS = Manager.GetDouble("mu-s");

  MultiColumnASCIIFile InputVectors;
  if (InputVectors.Parse(Manager.GetString("input-states")) == false)
    {
      InputVectors.DumpErrors(cout) << endl;
      return -1;
    }

  int CurrentKx = 0;
  int CurrentKy = 0;
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(InputVectors(0, 0), NbrParticles, NbrSitesX, NbrSitesY, CurrentKx, CurrentKy, Statistics) == false)
    {
      return -1;      
    }
  for (int i = 1; i < InputVectors.GetNbrLines(); ++i)
    {
      int TmpNbrParticles = 0;
      int TmpNbrSitesX = 0;
      int TmpNbrSitesY = 0;
      int TmpCurrentKx = 0;
      int TmpCurrentKy = 0;
      if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(InputVectors(0, i), TmpNbrParticles, TmpNbrSitesX, TmpNbrSitesY, TmpCurrentKx, TmpCurrentKy, Statistics) == false)
	{
	  return -1;      
	}
      if (NbrParticles != TmpNbrParticles)
	{
	  cout << InputVectors(0, i) << " and " << InputVectors(0, 0) << " have different numbers of particles" << endl;
	  return -1;
	}
      if (TmpNbrSitesX != NbrSitesX)
	{
	  cout << InputVectors(0, i) << " and " << InputVectors(0, 0) << " have different number of sites in the x direction" << endl;
	  return -1;
	}
      if (TmpNbrSitesY != NbrSitesY)
	{
	  cout << InputVectors(0, i) << " and " << InputVectors(0, 0) << " have different number of sites in the y direction" << endl;
	  return -1;
	}
    }
  Complex* Coefficients = 0;
  if (InputVectors(1, 0) != 0)
    {
      Coefficients = InputVectors.GetAsComplexArray(1);
    }
  else
    {
      if (CoefficientOnlyFlag == true)
	{
	  Coefficients = new Complex [InputVectors.GetNbrLines()];
	  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
	    {
	      Coefficients[i] = 1.0;
	    }
	}
      else
	{
	  cout << "no coefficients defined in " << Manager.GetString("input-states") << endl;
	  return -1; 
	}
    }
  
  //  ParticleOnDiskFunctionBasis Basis(TotalLz);

  int ForceMaxMomentum = NbrSitesX * NbrSitesY - 1;      
  //  double Scale = 2.0 * sqrt((double) TotalLz);
  ComplexMatrix PrecalculatedValues(ForceMaxMomentum + 1, ForceMaxMomentum + 1, true);
  ComplexMatrix** RawPrecalculatedValues = new ComplexMatrix*[InputVectors.GetNbrLines()];
  for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
    {
      RawPrecalculatedValues[i] = new ComplexMatrix[InputVectors.GetNbrLines()];
      for (int j = i; j < InputVectors.GetNbrLines(); ++j)
	RawPrecalculatedValues[i][j] = ComplexMatrix(ForceMaxMomentum + 1, ForceMaxMomentum + 1, true);
    }
  
  if (Manager.GetString("use-precomputed") == 0)
    {	  
      for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
	{
	  ParticleOnSphere* LeftSpace = 0;
	  ComplexVector LeftState;
	  int LeftKxMomentum = 0;
	  int LeftKyMomentum = 0;
	  if (FQHECheckerboardLatticeModelDensityGetHilbertSpace(InputVectors(0, i), NbrParticles, NbrSitesX, NbrSitesY, 
								 LeftKxMomentum, LeftKyMomentum, Statistics, 
								 LeftSpace, LeftState) == false)
	    return -1;
	  for (int m = 0; m <= ForceMaxMomentum; ++m)
	    {
	      int TotalKx = m / NbrSitesY;
	      int TotalKy = m % NbrSitesY;
	      ParticleOnSphereDensityOperator Operator (LeftSpace, m);
	      OperatorMatrixElementOperation Operation(&Operator, LeftState, LeftState);
	      Operation.ApplyOperation(Architecture.GetArchitecture());
	      //	      Complex Tmp = Operation.GetScalar();
	      Complex Tmp = Operator.MatrixElement(LeftState, LeftState);
	      RawPrecalculatedValues[i][i].AddToMatrixElement(m, m, Tmp);
	      Tmp *= (Conj(Coefficients[i]) * Coefficients[i]);
	      PrecalculatedValues.AddToMatrixElement(m, m, Tmp);
	    }
	  for (int j = i + 1; j < InputVectors.GetNbrLines(); ++j)
	    {
	      ParticleOnSphere* RightSpace = 0;
	      ComplexVector RightState;
	      int RightKxMomentum = 0;
	      int RightKyMomentum = 0;
	      if (FQHECheckerboardLatticeModelDensityGetHilbertSpace(InputVectors(0, j), NbrParticles, NbrSitesX, NbrSitesY, 
								     RightKxMomentum, RightKyMomentum, Statistics, 
								     RightSpace, RightState) == false)
		return -1;
	      RightSpace->SetTargetSpace(LeftSpace);
	      for (int m = 0; m <= ForceMaxMomentum; ++m)
		{
		  for (int n = 0; n <= ForceMaxMomentum; ++n)
		    {
		      int TotalKx1 = m / NbrSitesY;
		      int TotalKy1 = m % NbrSitesY;
		      int TotalKx2 = n / NbrSitesY;
		      int TotalKy2 = n % NbrSitesY;
		      if ((((RightKxMomentum - TotalKx2) % NbrSitesX) == ((LeftKxMomentum - TotalKx1) % NbrSitesX))
			  && (((RightKyMomentum - TotalKy2) % NbrSitesY) == ((LeftKyMomentum - TotalKy1) % NbrSitesY)))
			{
			  ParticleOnSphereDensityOperator Operator (RightSpace, m, n);
			  OperatorMatrixElementOperation Operation(&Operator, LeftState, RightState);
			  Operation.ApplyOperation(Architecture.GetArchitecture());
			  //Complex Tmp = Operation.GetScalar();
			  Complex Tmp = -Operator.MatrixElement(LeftState, RightState);
			  RawPrecalculatedValues[i][j].AddToMatrixElement(m, n, Tmp);
			  Tmp *=  (Conj(Coefficients[i]) * Coefficients[j]);
			  PrecalculatedValues.AddToMatrixElement(m, n, Tmp);
			  PrecalculatedValues.AddToMatrixElement(n, m, Conj(Tmp));
			}
		    }
		}
	      delete RightSpace;
	    }
	}
    }
  else
    {
      MultiColumnASCIIFile PrecomputedElements;
      if (PrecomputedElements.Parse(Manager.GetString("use-precomputed")) == false)
	{
	  PrecomputedElements.DumpErrors(cout) << endl;
	  return -1;
	}
      int* StateIndexLeft = PrecomputedElements.GetAsIntegerArray(0);
      int* StateIndexRight = PrecomputedElements.GetAsIntegerArray(1);
      int* OperatorIndexLeft = PrecomputedElements.GetAsIntegerArray(2);
      int* OperatorIndexRight = PrecomputedElements.GetAsIntegerArray(3);	  
      Complex* OneBodyCoefficients = PrecomputedElements.GetAsComplexArray(4);	  
      for (int i = 0; i < PrecomputedElements.GetNbrLines(); ++i)
	{
	  RawPrecalculatedValues[StateIndexLeft[i]][StateIndexRight[i]].AddToMatrixElement(OperatorIndexLeft[i], OperatorIndexRight[i], OneBodyCoefficients[i]);
	}
      for (int i = 0; i < InputVectors.GetNbrLines(); ++i)
	{
	  for (int m = 0; m <= ForceMaxMomentum; ++m)
	    {
	      Complex Tmp = RawPrecalculatedValues[i][i][m][m];
	      Tmp *= (Coefficients[i] * Coefficients[i]);
	      PrecalculatedValues.AddToMatrixElement(m, m, Tmp);
	    }
	  for (int j = i + 1; j < InputVectors.GetNbrLines(); ++j)
	    {
	      for (int m = 0; m <= ForceMaxMomentum; ++m)
		{
		  for (int n = 0; n <= ForceMaxMomentum; ++n)
		    {
		      Complex Tmp = RawPrecalculatedValues[i][j][m][n];
		      Tmp *= (Coefficients[i] * Coefficients[j]);
		      PrecalculatedValues.AddToMatrixElement(m, n, Tmp);
		      PrecalculatedValues.AddToMatrixElement(n, m, Conj(Tmp));
		    }
		}
	    }
	}
    }

  ComplexMatrix* OneBodyBasis = new ComplexMatrix [NbrSitesX * NbrSitesY];
  double GammaX = 0.0;
  double GammaY = 0.0;
  for (int kx = 0; kx < NbrSitesX; ++kx)
    for (int ky = 0; ky < NbrSitesY; ++ky)
      {
	int Index = (kx * NbrSitesY) + ky;
	Complex B1 = 4.0 * NNHoping * Complex (cos (1.0 * M_PI * (((double) kx) + GammaX) / ((double) NbrSitesX)) * cos (1.0 * M_PI * (((double) ky) + GammaY) / ((double) NbrSitesY)) * cos(M_PI * 0.25), 
						     sin (1.0 * M_PI * (((double) kx) + GammaX) / ((double) NbrSitesX)) * sin (1.0 * M_PI * (((double) ky) + GammaY) / ((double) NbrSitesY)) * sin(M_PI * 0.25));
	double d1 = 4.0 * SecondNextNNHoping * cos (2.0 * M_PI * (((double) kx) + GammaX) / ((double) NbrSitesX)) * cos (2.0 * M_PI * (((double) ky) + GammaY) / ((double) NbrSitesY));
	double d3 =  MuS + (2.0 * NextNNHoping * (cos (2.0 * M_PI * (((double) kx) + GammaX) / ((double) NbrSitesX)) 
							      - cos (2.0 * M_PI * (((double) ky) + GammaY) / ((double) NbrSitesY))));
	HermitianMatrix TmpOneBobyHamiltonian(2, true);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 0, d1 + d3);
	TmpOneBobyHamiltonian.SetMatrixElement(0, 1, B1);
	TmpOneBobyHamiltonian.SetMatrixElement(1, 1, d1 - d3);
	ComplexMatrix TmpMatrix(2, 2, true);
	TmpMatrix[0][0] = 1.0;
	TmpMatrix[1][1] = 1.0;
	RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	TmpOneBobyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
	TmpOneBobyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif   
	OneBodyBasis[Index] = TmpMatrix;	
      }

  ofstream File;
  File.precision(14);
  if (OutputName == 0)
    OutputName = ReplaceExtensionToFileName(Manager.GetString("input-states"), "dat", "rho.dat");
  File.open(OutputName, ios::binary | ios::out);
  File << "# density coefficients for " << Manager.GetString("input-states") << endl;
  if (CoefficientOnlyFlag == false)
    {
      File << "#" << endl << "# m  n  c_{m,n}" << endl;
      for (int i = 0; i <= ForceMaxMomentum; ++i)
	for (int j = 0; j <= ForceMaxMomentum; ++j)
	  File << "# " << i << " " << j << " " << PrecalculatedValues[i][j] << endl;
    }
  else
    {
      File << "#" << endl << "# state_index_left state_index_right m  n  c_{m,n}" << endl;
      for (int m = 0; m < InputVectors.GetNbrLines(); ++m)
	for (int n = m; n < InputVectors.GetNbrLines(); ++n)
	  for (int i = 0; i <= ForceMaxMomentum; ++i)
	    for (int j = 0; j <= ForceMaxMomentum; ++j)
	      {
		if ((RawPrecalculatedValues[m][n][i][j].Re != 0.0) || (RawPrecalculatedValues[m][n][i][j].Im != 0.0))
		  File << m << " " << n << " " << i << " " << j << " " << RawPrecalculatedValues[m][n][i][j] << endl;
	      }
    }
  
  if (CoefficientOnlyFlag == false)
    {
      Complex Sum = 0.0;
      double Normalization = 1.0 / (((double) NbrSitesX) * ((double) NbrSitesY));      
      for (int i = 0; i < NbrSitesX; ++i)
	{
	  for (int j = 0; j < NbrSitesY; ++j)
	    {
	      Complex DensityA = 0.0;
	      for (int m = 0; m <= ForceMaxMomentum; ++m)
		{
		  for (int n = 0; n <= ForceMaxMomentum; ++n)
		    {
		      int TotalKx1 = m / NbrSitesY;
		      int TotalKy1 = m % NbrSitesY;
		      int TotalKx2 = n / NbrSitesY;
		      int TotalKy2 = n % NbrSitesY;
		      DensityA += Phase(-2.0 * M_PI * ((double) ((TotalKx1 - TotalKx2) * i)) / ((double) NbrSitesX) 
					- 2.0 * M_PI * ((double) ((TotalKy1 - TotalKy2) * j)) / ((double) NbrSitesY)) * Normalization * Conj(OneBodyBasis[m][0][0]) * OneBodyBasis[n][0][0] * PrecalculatedValues[m][n];
		    }
		}
	      Sum += DensityA;
	      File << i << " " << j << " " << DensityA.Re << endl;
	    }
	  File << endl;
	  for (int j = 0; j < NbrSitesY; ++j)
	    {
	      Complex DensityB = 0.0;
	      for (int m = 0; m <= ForceMaxMomentum; ++m)
		{
		  for (int n = 0; n <= ForceMaxMomentum; ++n)
		    {
		      int TotalKx1 = m / NbrSitesY;
		      int TotalKy1 = m % NbrSitesY;
		      int TotalKx2 = n / NbrSitesY;
		      int TotalKy2 = n % NbrSitesY;
		      DensityB += Phase(-2.0 * M_PI * (((double) (TotalKx1 - TotalKx2)) * (0.5 + (double) i)) / ((double) NbrSitesX) 
					- 2.0 * M_PI * (((double) (TotalKy1 - TotalKy2)) * (0.5 + (double) j)) / ((double) NbrSitesY)) * Normalization * OneBodyBasis[m][0][1] * Conj(OneBodyBasis[n][0][1]) * PrecalculatedValues[m][n];
		    }
		}
	      Sum += DensityB;
	      File << (((double) i) + 0.5) << " " << (((double) j) + 0.5) << " " << DensityB.Re << endl;	    
	    }
	  File << endl;
	}
      cout << "Sum = " << Sum << endl;
    }
  File.close();
}

// get the Hilbert space and the vector state form the input file name
//
// inputState = input file name
// nbrParticles = reference on the number of particles 
// nbrSitesX = reference on the number of site along the x direction
// nbrSiteY = reference on the number of site along the y direction
// kxMomentum = reference on the momentum along the x direction
// kyMomentum = reference on the momentum along the y direction
// statistics = reference on the statistic flag
// space = reference on the pointer to the Hilbert space
// state = reference on the state vector
// return value = true if no error occured

bool FQHECheckerboardLatticeModelDensityGetHilbertSpace(char* inputState, int& nbrParticles, int& nbrSitesX, int& nbrSitesY,
							int& kxMomentum, int& kyMomentum, bool& statistics,
							ParticleOnSphere*& space, ComplexVector& state)
{
  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(inputState, nbrParticles, nbrSitesX, nbrSitesY, kxMomentum, kyMomentum, statistics) == false)
    {
      return false;      
    }
  if (state.ReadVector (inputState) == false)
    {
      cout << "can't open vector file " << inputState << endl;
      return false;      
    }

  if (statistics == true)
    {
#ifdef __64_BITS__
      if ((nbrSitesX * nbrSitesY) <= 63)
#else
	if ((nbrSitesX * nbrSitesY) <= 31)
#endif
	  {
	    space = new FermionOnSquareLatticeMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, kxMomentum, kyMomentum);
	  }
	else
	  {
	    space = new FermionOnSquareLatticeMomentumSpaceLong (nbrParticles, nbrSitesX, nbrSitesY, kxMomentum, kyMomentum);
	  }
    }
  else
    {
      space = new BosonOnSquareLatticeMomentumSpace (nbrParticles, nbrSitesX, nbrSitesY, kxMomentum, kyMomentum);
    }
  if (space->GetLargeHilbertSpaceDimension() != state.GetLargeVectorDimension())
    {
      cout << "dimension mismatch between the state (" << state.GetLargeVectorDimension() << ") and the Hilbert space (" << space->GetLargeHilbertSpaceDimension() << ")" << endl;
      return false;
    }
  return true;
}
