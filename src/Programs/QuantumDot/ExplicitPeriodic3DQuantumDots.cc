#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Complex.h"

#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "HilbertSpace/QuantumDotHilbertSpace/Periodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/XYReflexionSymmetricPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/ImpairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/ImpairXPairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PairXPairYPeriodic3DOneParticle.h"

#include "Hamiltonian/ExplicitHamiltonian.h"
#include "Hamiltonian/QuantumDotHamiltonian/PeriodicQuantumDots3DHamiltonian.h"
#include "Hamiltonian/QuantumDotHamiltonian/XYReflexionSymmetricPeriodic3DHamiltonian.h"

#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"

#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/QuantumDot/Potential/ThreeDConstantCellPotential.h"
#include "Tools/QuantumDot/Potential/DotEmbeddedWellThreeDConstantCellPotential.h"

#include <iostream>
#include <stdlib.h>
#include <fstream.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>

using std::cout;
using std::endl;
using std::ostream;
using std::ios;
using std::ofstream;

bool EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray);

int main(int argc, char** argv)
{
  cout.precision(14);
  // some running options and help
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption EigenstateOption ('e', "eigenstate", "evaluate eigenstates", false);
  SingleIntegerOption IterationOption ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 50);
  SingleIntegerOption NumberXValueOption ('M', "M-cell", "number of cells in the x direction", 160);
  SingleIntegerOption NumberYValueOption ('N', "N-cell", "number of cells in the y direction", 160);
  SingleIntegerOption NumberZValueOption ('H', "H-cell", "number of cells in the z direction", 21);
  SingleIntegerOption UnderBarrierValueOption ('\n', "barrier", "number of cells in the well barrier", 2);
  SingleIntegerOption BelowWettingLayerValueOption ('\n', "below", "number of cells between well barrier and wetting layer", 2);
  SingleIntegerOption WettingWidthOption ('\n', "wetting", "number of cells in wetting layer", 1);
  SingleIntegerOption BaseRadiusOption ('\n', "base", "base radius in cell unit", 18);
  SingleIntegerOption DotHeightOption ('\n', "height", "height of dot in cell unit", 7);
  SingleIntegerOption TopRadiusOption ('\n', "top", "top radius in cell unit", 6);
  SingleDoubleOption CellXSizeOption ('X', "cell-xsize", "cell size in the x direction in Angstrom", 5.65);
  SingleDoubleOption CellYSizeOption ('Y', "cell-ysize", "cell size in the y direction in Angstrom", 5.65);
  SingleDoubleOption CellZSizeOption ('Z', "cell-zsize", "cell size in the z direction in Angstrom", 5.65);
  SingleDoubleOption XMassOption ('\n', "mu-x", "electron effective mass in x direction (in vacuum electron mass unit)", 0.07);
  SingleDoubleOption YMassOption ('\n', "mu-y", "electron effective mass in y direction (in vacuum electron mass unit)", 0.07);
  SingleDoubleOption ZMassOption ('\n', "mu-z", "electron effective mass in z direction (in vacuum electron mass unit)", 0.07);
  SingleDoubleOption WellPotentialOption ('\n', "well", "potential in the well", 1.079);
  SingleDoubleOption DotPotentialOption ('\n', "dot", "potential in the dot", -0.73);
  BooleanOption DiskOption ('d', "disk", "enable disk resume capabilities", false);
  BooleanOption ResumeOption ('r', "resume", "resume from disk datas", false);
  SingleIntegerOption VectorMemoryOption ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 400);
  SingleIntegerOption NbrIterationOption ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 60);
  BooleanOption PairXOption ('\n', "pairX", "pair function in X direction", false);
  BooleanOption PairYOption ('\n', "pairY", "pair function in Y direciton", false);

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &EigenstateOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &NumberXValueOption;
  OptionList += &NumberYValueOption;
  OptionList += &NumberZValueOption;
  OptionList += &UnderBarrierValueOption;
  OptionList += &BelowWettingLayerValueOption;
  OptionList += &WettingWidthOption;
  OptionList += &BaseRadiusOption;
  OptionList += &DotHeightOption;
  OptionList += &TopRadiusOption;
  OptionList += &CellXSizeOption;
  OptionList += &CellYSizeOption;
  OptionList += &CellZSizeOption;
  OptionList += &XMassOption;
  OptionList += &YMassOption;
  OptionList += &ZMassOption;
  OptionList += &WellPotentialOption;
  OptionList += &DotPotentialOption;
  OptionList += &VectorMemoryOption;
  OptionList += &DiskOption;
  OptionList += &ResumeOption;
  OptionList += &NbrIterationOption;
  OptionList += &PairXOption;
  OptionList += &PairYOption;

  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type ExplicitMatrixExample -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }

  bool SMPFlag = SMPOption.GetBoolean();
  bool EigenstateFlag = EigenstateOption.GetBoolean();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  int M = NumberXValueOption.GetInteger();
  int N = NumberYValueOption.GetInteger();
  int H = NumberZValueOption.GetInteger();
  int UnderBarrier = UnderBarrierValueOption.GetInteger();
  int BelowWettingLayer = BelowWettingLayerValueOption.GetInteger();
  int WettingWidth = WettingWidthOption.GetInteger();
  int BaseRadius = BaseRadiusOption.GetInteger();
  int DotHeight = DotHeightOption.GetInteger();
  int TopRadius = TopRadiusOption.GetInteger();
  double Lx = CellXSizeOption.GetDouble();
  double Ly = CellYSizeOption.GetDouble();
  double Lz = CellZSizeOption.GetDouble();
  double Mux = XMassOption.GetDouble();
  double Muy = YMassOption.GetDouble();
  double Muz = ZMassOption.GetDouble();
  double WellPotential = WellPotentialOption.GetDouble();
  double DotPotential = DotPotentialOption.GetDouble();
  int VectorMemory = VectorMemoryOption.GetInteger();
  bool ResumeFlag = ResumeOption.GetBoolean();
  bool DiskFlag = DiskOption.GetBoolean();
  int NbrIterLanczos = NbrIterationOption.GetInteger();
  bool PairX = PairXOption.GetBoolean();
  bool PairY = PairYOption.GetBoolean();   

  // DotEmbeddedWellThreeDConstantCellPotential(int numberX, int numberY, int numberZ, int underBarrier, int belowWettingLayer, int wettingWidth, int baseRadius, int dotHeight, int topRadius)
  DotEmbeddedWellThreeDConstantCellPotential* potential = new DotEmbeddedWellThreeDConstantCellPotential(M, N, H, UnderBarrier, BelowWettingLayer, WettingWidth, BaseRadius, DotHeight, TopRadius);

  // ConstructPotential(double wellPotential, double dotPotential)
  //potential->ConstructPotential(WellPotential, DotPotential);
  potential->LoadPotential("DotPotential.txt");

  Periodic3DOneParticle* Space = new Periodic3DOneParticle(M / 2 + 1, -M / 4, N / 2 + 1, -N / 4, H, -H / 2);

  int NbrCellX = M, NbrCellY = N, NbrCellZ = H; 
  int NbrStateX = Space->GetNbrStateX(), NbrStateY = Space->GetNbrStateY(), NbrStateZ = Space->GetNbrStateZ();  
  int LowerImpulsionX = Space->GetLowerImpulsionX(), LowerImpulsionY = Space->GetLowerImpulsionY(), LowerImpulsionZ = Space->GetLowerImpulsionZ();

  double PERIODIC_HAMILTONIAN_FACTOR = 150.4;
 
  double** RealWaveFunctionOverlapX; double** ImaginaryWaveFunctionOverlapX;
  double** RealWaveFunctionOverlapY; double** ImaginaryWaveFunctionOverlapY; 
  double** RealWaveFunctionOverlapZ; double** ImaginaryWaveFunctionOverlapZ;
  double XSize = NbrCellX * Lx, YSize = NbrCellY * Ly, ZSize = NbrCellZ * Lz;
  int Dimension = NbrStateX * NbrStateY * NbrStateZ;

  if (!EvaluateWaveFunctionOverlap(NbrCellX, NbrStateX, RealWaveFunctionOverlapX, ImaginaryWaveFunctionOverlapX))
    cout << "Error in evaluation of function overlap in X direction. Stop!" << endl;  
  if (!EvaluateWaveFunctionOverlap(NbrCellY, NbrStateY, RealWaveFunctionOverlapY, ImaginaryWaveFunctionOverlapY))
    cout << "Error in evaluation of function overlap in Y direction. Stop!" << endl;
  if (!EvaluateWaveFunctionOverlap(NbrCellZ, NbrStateZ, RealWaveFunctionOverlapZ, ImaginaryWaveFunctionOverlapZ))
    cout << "Error in evaluation of function overlap in Z direction. Stop!" << endl;

  double InvXFactor = PERIODIC_HAMILTONIAN_FACTOR / (Mux * XSize * XSize);
  double InvYFactor = PERIODIC_HAMILTONIAN_FACTOR / (Muy * YSize * YSize);
  double InvZFactor = PERIODIC_HAMILTONIAN_FACTOR / (Muz * ZSize * ZSize);
  
  double* KineticElements = new double[Dimension];

  double FactorX = 0.0, FactorY = 0.0;
  int TotalIndex1 = 0;
  for (int i = 0; i < NbrStateX; ++i)
    {
      FactorX = double((i + LowerImpulsionX) * (i + LowerImpulsionX)) * InvXFactor;
      for (int j = 0; j < NbrStateY; ++j)
	{
	  FactorY = double((j + LowerImpulsionY) * (j + LowerImpulsionY)) * InvYFactor + FactorX;
	  for (int k = 0; k < NbrStateZ; ++k)
	    {	      
	      KineticElements[TotalIndex1] = FactorY + double((k + LowerImpulsionZ) * (k + LowerImpulsionZ)) * InvZFactor;	      
	      ++TotalIndex1;
	    }
	}
    }

  int LengthX = (NbrStateX - 1) * 2 + 1; int LengthY = (NbrStateY - 1) * 2 + 1; int LengthZ = (NbrStateZ - 1) * 2 + 1;

  double*** TmpReal = new double** [LengthX];
  double*** TmpImaginary = new double** [LengthX];

  double TmpRe, TmpIm;
  double TmpRe2, TmpIm2;
  double* TmpRealWaveFunctionOverlapX;
  double* TmpImaginaryWaveFunctionOverlapX;
  double* TmpRealWaveFunctionOverlapY;
  double* TmpImaginaryWaveFunctionOverlapY;
  double* TmpRealPrecalculatedHamiltonian;
  double* TmpImaginaryPrecalculatedHamiltonian;

  for (int m = 0; m < LengthX; ++m)
    {
      TmpReal[m] = new double* [LengthY];
      TmpImaginary[m] = new double* [LengthY];
      TmpRealWaveFunctionOverlapX = RealWaveFunctionOverlapX[m];
      TmpImaginaryWaveFunctionOverlapX = ImaginaryWaveFunctionOverlapX[m];	      	  
      for (int n = 0; n < LengthY; ++n)
	{	  
	  TmpReal[m][n] = new double [NbrCellZ];
	  TmpImaginary[m][n] = new double [NbrCellZ];
	  TmpRealWaveFunctionOverlapY = RealWaveFunctionOverlapY[n];
	  TmpImaginaryWaveFunctionOverlapY = ImaginaryWaveFunctionOverlapY[n];	  
	  TmpRealPrecalculatedHamiltonian = TmpReal[m][n];
	  TmpImaginaryPrecalculatedHamiltonian = TmpImaginary[m][n];		  
	  for (int CellZ = 0; CellZ < NbrCellZ; ++CellZ)
	    {
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (int CellY = 0; CellY < NbrCellY; ++CellY)
		{
		  TmpRe2 = TmpRealWaveFunctionOverlapY[CellY];
		  TmpIm2 = TmpImaginaryWaveFunctionOverlapY[CellY];
		  for (int CellX = 0; CellX < NbrCellX; ++CellX)
		    {		      
		      TmpRe += potential->GetPotential(CellX, CellY, CellZ) * (TmpRealWaveFunctionOverlapX[CellX] * TmpRe2 - TmpImaginaryWaveFunctionOverlapX[CellX] * TmpIm2);
		      TmpIm += potential->GetPotential(CellX, CellY, CellZ) * (TmpRealWaveFunctionOverlapX[CellX] * TmpIm2 + TmpImaginaryWaveFunctionOverlapX[CellX] * TmpRe2);		      
		    }
		}
	      TmpRealPrecalculatedHamiltonian[CellZ] = TmpRe;  
	      TmpImaginaryPrecalculatedHamiltonian[CellZ] = TmpIm;  
	    }
	}
    }

  double*** RealPrecalculatedHamiltonian = new double** [LengthX];
  double*** ImaginaryPrecalculatedHamiltonian = new double** [LengthX];
  double* TmpRealWaveFunctionOverlapZ;
  double* TmpImaginaryWaveFunctionOverlapZ;
  for (int m = 0; m < LengthX; ++m)
    {
      RealPrecalculatedHamiltonian[m] = new double* [LengthY];      
      ImaginaryPrecalculatedHamiltonian[m] = new double* [LengthY]; 
      for (int n = 0; n < LengthY; ++n)
	{
	  RealPrecalculatedHamiltonian[m][n] = new double [LengthZ];      
	  ImaginaryPrecalculatedHamiltonian[m][n] = new double [LengthZ]; 
	  TmpRealPrecalculatedHamiltonian = TmpReal[m][n];
	  TmpImaginaryPrecalculatedHamiltonian = TmpImaginary[m][n];
	  for (int p = 0; p < LengthZ; ++p)
	    {
	      TmpRealWaveFunctionOverlapZ = RealWaveFunctionOverlapZ[p];
	      TmpImaginaryWaveFunctionOverlapZ = ImaginaryWaveFunctionOverlapZ[p];
	      TmpRe = 0.0; TmpIm = 0.0;
	      for (int CellZ = 0; CellZ < NbrCellZ; ++CellZ)
		{
		  TmpRe += (TmpRealPrecalculatedHamiltonian[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ] - TmpImaginaryPrecalculatedHamiltonian[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ]);
		  TmpIm += (TmpRealPrecalculatedHamiltonian[CellZ] * TmpImaginaryWaveFunctionOverlapZ[CellZ] + TmpImaginaryPrecalculatedHamiltonian[CellZ] * TmpRealWaveFunctionOverlapZ[CellZ]);
		}
	      RealPrecalculatedHamiltonian[m][n][p] = TmpRe;
	      ImaginaryPrecalculatedHamiltonian[m][n][p] = TmpIm;
	    }
	}
    }
  delete[] TmpReal; delete[] TmpImaginary;

  HermitianMatrix HamiltonianRepresentation (Dimension);
  
  int m1, m2, n1, n2, p1, p2;
  int IndexX, IndexY, IndexZ;
  int TmpIndex = 0; int** TotalIndex = new int* [NbrStateX];
  for (m1 = 0; m1 < NbrStateX; ++m1) 
    {
      TotalIndex[m1] = new int [NbrStateY];
      for (n1 = 0; n1 < NbrStateY; ++n1)	
	{
	  TotalIndex[m1][n1] = (m1 * NbrStateY + n1) * NbrStateZ;
	  for (p1 = 0; p1 < NbrStateZ; ++p1)
	    {	      
	      HamiltonianRepresentation.SetMatrixElement(TmpIndex, TmpIndex, KineticElements[TmpIndex]);
	      ++TmpIndex;
	    }
	}
    }

  int OriginX = NbrStateX - 1; int OriginY = NbrStateY - 1; int OriginZ = NbrStateZ - 1;
  int Index1, Index2;
  Complex Tmp;
  int* TmpTotalIndex1; int* TmpTotalIndex2;
  for (m1 = 0; m1 < NbrStateX; ++m1)
    {
      for (n1 = 0; n1 < NbrStateY; ++n1)
	{
	  Index1 = TotalIndex[m1][n1];
	  for (p1 = 0; p1 < NbrStateZ; ++p1)
	    {	      
	      for (m2 = 0; m2 < NbrStateX; ++m2)
		{		  
		  IndexX = -m1 + m2 + OriginX;		  
		  for (n2 = 0; n2 < NbrStateY; ++n2)
		    {
		      IndexY = -n1 + n2 + OriginY;
		      Index2 = TotalIndex[m2][n2];
		      for (p2 = 0; p2 < NbrStateZ; ++p2)
			{
			  IndexZ = -p1 + p2 + OriginZ;
			  Tmp = Complex(RealPrecalculatedHamiltonian[IndexX][IndexY][IndexZ], ImaginaryPrecalculatedHamiltonian[IndexX][IndexY][IndexZ]);
			  HamiltonianRepresentation.SetMatrixElement(Index1, Index2, Tmp);
			  ++Index2;
			}

		    }
		}
	      ++Index1;
	    }
	}
    }


  // RealVector* Eigenstates = 0;
  double* Eigenvalues = new double [NbrEigenvalue];

  // find the eigenvalues (and eigenvectors if needed)
  
  // architecture type (i.e. 1 CPU or multi CPU)
  AbstractArchitecture* Architecture;
  if (SMPFlag == true)
    Architecture = new SMPArchitecture(2);
  else
    Architecture = new MonoProcessorArchitecture;
  
  double Precision;
  double PreviousLowest;
  double Lowest;
  int CurrentNbrIterLanczos;

  // type of lanczos algorithm (with or without reorthogonalization)
  // BasicLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
  FullReorthogonalizedComplexLanczosAlgorithm Lanczos(Architecture, NbrEigenvalue, MaxNbrIterLanczos);      
     
  ExplicitHamiltonian Hamiltonian(Space, &HamiltonianRepresentation);
  
  
  // initialization of lanczos algorithm
  Precision = 1.0;
  PreviousLowest = 1e50;
  Lowest = PreviousLowest;
  CurrentNbrIterLanczos = NbrEigenvalue + 3;
  Lanczos.SetHamiltonian(&Hamiltonian);
  Lanczos.InitializeLanczosAlgorithm();
  Lanczos.RunLanczosAlgorithm(NbrEigenvalue + 2);
  RealTriDiagonalSymmetricMatrix TmpMatrix;
  
  // run Lancos algorithm up to desired precision on the n-th eigenvalues
  while ((Precision > 1e-14) && (CurrentNbrIterLanczos++ < MaxNbrIterLanczos))
    {
      Lanczos.RunLanczosAlgorithm(1);
      TmpMatrix.Copy(Lanczos.GetDiagonalizedMatrix());
      TmpMatrix.SortMatrixUpOrder();
      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
      cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << endl;
      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
      PreviousLowest = Lowest; 
    }      
  if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
    {
      cout << "too much Lanczos iterations" << endl;
      exit(0);
    }
  
  // store eigenvalues      
  for (int i = 0; i < NbrEigenvalue; ++i)    
    {
      Eigenvalues[i] = TmpMatrix.DiagonalElement(i);
      cout << Eigenvalues[i] << '\t';
    }

  return 0;
}

// evaluate the wave function overlap
//
// nbrStep = number of steps in the given direction
// nbrState = number of states chosen for this direction
// realArray = 2D array containing the real elements of the overlap
// imaginaryArray = 2D array containing the imaginary elements of the overlap

bool EvaluateWaveFunctionOverlap(int nbrStep, int nbrState, double** &realArray, double** &imaginaryArray)
{
  double Diff = 0.0;
  double Tmp = 0.0;
  double Tmp1 = 1.0 / double (nbrStep);
  int Length = (nbrState - 1) * 2 + 1;
  realArray = new double* [Length];
  imaginaryArray = new double* [Length];  
  int Origin = nbrState - 1;
  for (int delta = 0; delta < Length; ++delta)
    {
      realArray[delta] = new double [nbrStep];
      imaginaryArray[delta] = new double [nbrStep];
      if (delta != Origin)
	{
	  Diff = 2.0 * M_PI * double (delta - Origin);
	  Tmp = Diff / nbrStep;	
	  Diff = 1.0 / Diff;	
	  for (int i = 0; i < nbrStep; ++i)
	    {
	      realArray[delta][i] = Diff * (sin(Tmp * (i + 1)) - sin(Tmp * i));
	      imaginaryArray[delta][i] = Diff * (cos(Tmp * (i + 1)) - cos(Tmp * i));
	    }
	}
      else
	for (int i = 0; i < nbrStep; ++i)
	  {
	    realArray[delta][i] = Tmp1;
	    imaginaryArray[delta][i] = 0.0;
	  }	
    }
  return true;
}
