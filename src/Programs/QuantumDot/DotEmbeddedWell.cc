#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "HilbertSpace/QuantumDotHilbertSpace/Periodic3DOneParticle.h"

#include "Hamiltonian/QuantumDotHamiltonian/PeriodicQuantumDots3DHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
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


int main(int argc, char** argv)
{
  cout.precision(14);
  // some running options and help
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption VerboseOption ('v', "verbose", "verbose mode", false);
  BooleanOption EigenstateOption ('e', "eigenstate", "evaluate eigenstates", false);
  SingleIntegerOption IterationOption ('i', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 100);
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

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &VerboseOption;
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
  bool VerboseFlag = VerboseOption.GetBoolean();
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

  // DotEmbeddedWellThreeDConstantCellPotential(int numberX, int numberY, int numberZ, int underBarrier, int belowWettingLayer, int wettingWidth, int baseRadius, int dotHeight, int topRadius)
  DotEmbeddedWellThreeDConstantCellPotential* potential = new DotEmbeddedWellThreeDConstantCellPotential(M, N, H, UnderBarrier, BelowWettingLayer, WettingWidth, BaseRadius, DotHeight, TopRadius);

  // ConstructPotential(double wellPotential, double dotPotential)
  potential->ConstructPotential(WellPotential, DotPotential);

  // define Hilbert space
  Periodic3DOneParticle Space(M / 2, M / 4, N / 2, N / 4, H, H / 2);
  timeval PrecalculationStartingTime;
  timeval PrecalculationEndingTime;
  gettimeofday (&(PrecalculationStartingTime), 0);

  PeriodicQuantumDots3DHamiltonian Hamiltonian(&Space, Lx * ((double) M), Ly * ((double) N),  Lz * ((double) H), Mux, Muy, Muz, M, N, H, potential);

  gettimeofday (&(PrecalculationEndingTime), 0);
  double Dt = (double) (PrecalculationEndingTime.tv_sec - PrecalculationStartingTime.tv_sec) +
    ((PrecalculationEndingTime.tv_usec - PrecalculationStartingTime.tv_usec) / 1000000.0);
  cout << "precalculation time = " << Dt << endl;

  ComplexVector* Eigenstates = 0;
  double* Eigenvalues = 0;

  // architecture type (i.e. 1 CPU or multi CPU)
  AbstractArchitecture* Architecture;
  if (SMPFlag == true)
    Architecture = new SMPArchitecture(2);
  else
	Architecture = new MonoProcessorArchitecture;
  
  double HamiltonianShift = - (150.4 * ((1.0 / (Lx * Lx * Mux)) + (1.0 / (Ly * Ly * Muy)) + (1.0 / (Lz * Lz * Muz))));
  Hamiltonian.ShiftHamiltonian (HamiltonianShift);
  cout << "Shift:  " << HamiltonianShift << endl;
  // type of lanczos algorithm (with or without reorthogonalization)
  gettimeofday (&(PrecalculationStartingTime), 0);
  FullReorthogonalizedComplexLanczosAlgorithm Lanczos(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
  
  // initialization of lanczos algorithm
  double Precision = 1.0;
  double PreviousLowest = 1e50;
  double Lowest = PreviousLowest;
  int CurrentNbrIterLanczos = NbrEigenvalue + 3;
  Lanczos.SetHamiltonian(&Hamiltonian);
  Lanczos.InitializeLanczosAlgorithm();
  cout << "Run Lanczos Algorithm" << endl;
  Lanczos.RunLanczosAlgorithm(NbrEigenvalue + 2);
  RealTriDiagonalSymmetricMatrix TmpMatrix;
  
  // run Lancos algorithm up to desired precision on the n-th eigenvalues
  while (Lanczos.TestConvergence() == false)
    {
      ++CurrentNbrIterLanczos;
      Lanczos.RunLanczosAlgorithm(1);
      TmpMatrix.Copy(Lanczos.GetDiagonalizedMatrix());
      TmpMatrix.SortMatrixUpOrder();
      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
      PreviousLowest = Lowest;
      if (VerboseFlag == true)
	cout << CurrentNbrIterLanczos << '\t' <<  TmpMatrix.DiagonalElement(0) << '\t' << Lowest << endl;
    }
  if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
    {
      cout << "too much Lanczos iterations" << endl;
      exit(0);
    }
  
  // store eigenvalues
  Eigenvalues = new double [NbrEigenvalue];
  for (int i = 0; i < NbrEigenvalue; ++i)
    {
      Eigenvalues[i] = (TmpMatrix.DiagonalElement(i) - HamiltonianShift);
      cout << Eigenvalues[i] << '\t';
    }

  //compute eigenstates
  if (EigenstateFlag == true)
    Eigenstates = (ComplexVector*) Lanczos.GetEigenstates(NbrEigenvalue);
  
  gettimeofday (&(PrecalculationEndingTime), 0);
  Dt = (double) (PrecalculationEndingTime.tv_sec - PrecalculationStartingTime.tv_sec) +
    ((PrecalculationEndingTime.tv_usec - PrecalculationStartingTime.tv_usec) / 1000000.0);

  // insert here your code using the eigenvalues and the eigenvectors
  if (EigenstateFlag == true)
    {
      ofstream OutputFile;
      OutputFile.precision(14);
      OutputFile.open("eigenvalues", ios::binary | ios::out);
      for (int i = 0; i < NbrEigenvalue; ++i)
	OutputFile << Eigenvalues[i] << " ";
      OutputFile << endl;
      OutputFile.close();
    }

  if ((EigenstateFlag == true) && (Eigenstates != 0))
    {
      char* TmpFileName = new char[256];
      for (int i = 0; i < NbrEigenvalue; ++i)
	{
	  sprintf  (TmpFileName, "eigenvector.%d", i);
	  Eigenstates[i].WriteAsciiVector(TmpFileName);
	}
      delete[] TmpFileName;
    }
  
  return 0;
}
