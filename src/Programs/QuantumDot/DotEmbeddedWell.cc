#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "HilbertSpace/QuantumDotHilbertSpace/Periodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/XYReflexionSymmetricPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/ImpairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/ImpairXPairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PairXImpairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/PairXPairYPeriodic3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/InPlaneReflexionSymmetricPeriodic3DOneParticle.h"

#include "Hamiltonian/QuantumDotHamiltonian/PeriodicQuantumDots3DHamiltonian.h"
#include "Hamiltonian/QuantumDotHamiltonian/ReflexionSymmetricPeriodic3DHamiltonian.h"
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


int main(int argc, char** argv)
{  
  cout.precision(14);
  // some running options and help
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption EigenstateOption ('e', "eigenstate", "evaluate eigenstates", false);
  SingleIntegerOption IterationOption ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 50);
  SingleIntegerOption NumberXValueOption ('M', "M-cell", "number of cells in the x direction", 161);
  SingleIntegerOption NumberYValueOption ('N', "N-cell", "number of cells in the y direction", 161);
  SingleIntegerOption NumberZValueOption ('H', "H-cell", "number of cells in the z direction", 21);
  SingleIntegerOption UnderBarrierValueOption ('\n', "barrier", "number of cells in the well barrier", 2);
  SingleIntegerOption BelowWettingLayerValueOption ('\n', "below", "number of cells between well barrier and wetting layer", 2);
  SingleIntegerOption WettingWidthOption ('\n', "wetting", "number of cells in wetting layer", 1);
  SingleIntegerOption BaseRadiusOption ('\n', "base", "base radius in cell unit", 18);
  SingleIntegerOption DotHeightOption ('\n', "height", "height of dot in cell unit", 3);
  SingleIntegerOption TopRadiusOption ('\n', "top", "top radius in cell unit", 13);
  SingleDoubleOption AnisotropyOption('a', "anisotropy", "anisotropy factor", 0.0);
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
  OptionList += &AnisotropyOption;
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
  double Anisotropy = AnisotropyOption.GetDouble();
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
  potential->ConstructPotential(WellPotential, DotPotential, Anisotropy);
  //potential->LoadPotential("DotPotential.txt");

  //InPlaneReflexionSymmetricPeriodic3DOneParticle* Space = new InPlaneReflexionSymmetricPeriodic3DOneParticle(PairX, M / 4, (N / 4) * 2 + 1, -N / 4, H, -H / 2);  
  // define Hilbert space  
  
  XYReflexionSymmetricPeriodic3DOneParticle GeneralSpace(M / 4, N / 4, H, -H / 2);
  XYReflexionSymmetricPeriodic3DOneParticle* Space;
  if (PairX)
    if (PairY)
      Space = new PairXPairYPeriodic3DOneParticle(GeneralSpace);     
    else
      Space = new PairXImpairYPeriodic3DOneParticle(GeneralSpace); 
  else
     if (PairY)
      Space = new ImpairXPairYPeriodic3DOneParticle(GeneralSpace);     
    else
      Space = new ImpairXImpairYPeriodic3DOneParticle(GeneralSpace);        
  //Periodic3DOneParticle* Space = new Periodic3DOneParticle((M / 4) * 2 + 1, -M / 4, (N / 4) * 2 + 1, -N / 4, H, -H / 2);
  
  timeval PrecalculationStartingTime;
  timeval PrecalculationEndingTime;
  gettimeofday (&(PrecalculationStartingTime), 0);
  
  //cout << "General space dimension: " << GeneralSpace.GetHilbertSpaceDimension() << endl;
  cout << "Sample size in cell unit: " << M << '\t' << N << '\t' << H << endl;
  cout << "Hilbert space dimensions: " << Space->GetNbrStateX() << '\t' << Space->GetNbrStateY() << '\t' << Space->GetNbrStateZ() << endl;
  cout << "Minimal impulsions:       " << Space->GetLowerImpulsionX() << '\t' << Space->GetLowerImpulsionY() << '\t' << Space->GetLowerImpulsionZ() << endl;

  //ReflexionSymmetricPeriodic3DHamiltonian Hamiltonian(Space, PairX, Lx * ((double) M), Ly * ((double) N),  Lz * ((double) H), Mux, Muy, Muz, M, N, H, potential);  
  XYReflexionSymmetricPeriodic3DHamiltonian Hamiltonian(Space, PairX, PairY, Lx * ((double) M), Ly * ((double) N),  Lz * ((double) H), Mux, Muy, Muz, M, N, H, potential);
  //PeriodicQuantumDots3DHamiltonian Hamiltonian(Space, Lx * ((double) M), Ly * ((double) N),  Lz * ((double) H), Mux, Muy, Muz, M, N, H, potential);

  cout << endl;
  gettimeofday (&(PrecalculationEndingTime), 0);
  double Dt = (double) (PrecalculationEndingTime.tv_sec - PrecalculationStartingTime.tv_sec) +
    ((PrecalculationEndingTime.tv_usec - PrecalculationStartingTime.tv_usec) / 1000000.0);
  cout << "Precalculation time = " << Dt << endl;

  ComplexVector* Eigenstates = 0;
  double* Eigenvalues = 0;

  // architecture type (i.e. 1 CPU or multi CPU)
  AbstractArchitecture* Architecture;
  if (SMPFlag == true)
    Architecture = new SMPArchitecture(2);
  else
    Architecture = new MonoProcessorArchitecture;

  cout << "----------------------------------------------------------------" << endl;

  double HamiltonianShift = -Hamiltonian.MaxPartialDiagonalElement();
  Hamiltonian.ShiftHamiltonian (HamiltonianShift);
  cout << "Hamiltonian shift =  " << HamiltonianShift << endl;
  gettimeofday (&(PrecalculationStartingTime), 0);

  // type of lanczos algorithm (with or without reorthogonalization)
  AbstractLanczosAlgorithm* Lanczos;
  if (DiskFlag == false)
    Lanczos = new FullReorthogonalizedComplexLanczosAlgorithm(Architecture, NbrEigenvalue, MaxNbrIterLanczos);   
  else
    Lanczos = new FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage(Architecture, NbrEigenvalue, VectorMemory, MaxNbrIterLanczos);
  cout << "Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl; 

  // initialization of lanczos algorithm
  double Precision = 1.0;
  double PreviousLowest = 1e50;
  double Lowest = PreviousLowest;
  int CurrentNbrIterLanczos = 0;
  Lanczos->SetHamiltonian(&Hamiltonian);
  if ((DiskFlag == true) && (ResumeFlag == true))
    Lanczos->ResumeLanczosAlgorithm();
  else
    Lanczos->InitializeLanczosAlgorithm();
  cout << "------------------- Run Lanczos Algorithm ---------------------" << endl;
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&(TotalStartingTime), 0);
  if (ResumeFlag == false)
    {
      Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
      CurrentNbrIterLanczos = NbrEigenvalue + 3;
      if ((DiskFlag == true) && (CurrentNbrIterLanczos >= NbrIterLanczos))
        {
          NbrIterLanczos = CurrentNbrIterLanczos + 1;
        }
    }
  RealTriDiagonalSymmetricMatrix TmpMatrix;
  while ((Lanczos->TestConvergence() == false) &&  (((DiskFlag == true) && (CurrentNbrIterLanczos < NbrIterLanczos)) ||
                                                    ((DiskFlag == false) && (CurrentNbrIterLanczos < MaxNbrIterLanczos))))
    {
      ++CurrentNbrIterLanczos;
      Lanczos->RunLanczosAlgorithm(1);
      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
      TmpMatrix.SortMatrixUpOrder();
      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
      PreviousLowest = Lowest;
      cout << CurrentNbrIterLanczos << "\t" <<  TmpMatrix.DiagonalElement(0) - HamiltonianShift << "\t\t" << Lowest - HamiltonianShift << "\t\t" << Precision << endl;
    }
  if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
    {
      cout << "too much Lanczos iterations" << endl;
      exit(0);
    }
  cout << "------------------ Actual eigenvalues ------------------" << endl;
  // store eigenvalues
  Eigenvalues = new double [NbrEigenvalue];
  for (int i = 0; i < NbrEigenvalue; ++i)
    {
      Eigenvalues[i] = (TmpMatrix.DiagonalElement(i) - HamiltonianShift);
      cout << Eigenvalues[i] << '\t';
    }
  cout << endl;
  if (Lanczos->TestConvergence())
    {
      //compute eigenstates
      if (EigenstateFlag == true)
        Eigenstates = (ComplexVector*) Lanczos->GetEigenstates(NbrEigenvalue);
      
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
            
      cout << "----------------- End of calculation ---------------------" << endl;      
      cout << "     ==========  CALCULATION IS FINALIZED  =========  " << endl;
      cout << "Sample size in cell unit: " << M << '\t' << N << '\t' << H << endl;
      cout << "Hilbert space dimensions: " << Space->GetNbrStateX() << '\t' << Space->GetNbrStateY() << '\t' << Space->GetNbrStateZ() << endl;
      cout << "Minimal impulsions:       " << Space->GetLowerImpulsionX() << '\t' << Space->GetLowerImpulsionY() << '\t' << Space->GetLowerImpulsionZ() << endl;
    }
  gettimeofday (&(TotalEndingTime), 0);
  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);  
  cout << endl << "Total time = " << Dt << endl;
  delete Lanczos;
  
  return 0;
}
