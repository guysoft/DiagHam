#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/QuantumDotHilbertSpace/VerticalPeriodicParticleInMagneticField.h"

#include "Hamiltonian/QuantumDotHamiltonian/CylindricalHamiltonianInMagneticField.h"

#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"

#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/QuantumDot/Potential/ThreeDConstantCylinderPotential.h"
#include "Tools/QuantumDot/Potential/QuantumDotThreeDConstantCylinderPotential.h"

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
 
  // QuantumDotThreeDConstantCylinderPotential(double belowHeight, double wettingWidth, int nbrCylinderDot, double dotHeight, double baseRadius, double topRadius, double aboveHeight);
  cout.precision(14);
  // some running options and help  
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption EigenstateOption ('e', "eigenstate", "evaluate eigenstates", false);
  SingleIntegerOption IterationOption ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 50);
  SingleIntegerOption NumberRValueOption ('R', "R-states", "number of states in plane", 100);
  SingleIntegerOption NumberZValueOption ('Z', "Z-states", "number of cells in z direction", 100);
  SingleIntegerOption NumberMValueOption ('m', "momentum", "quantum number of kinetic in z direction", 0);

  SingleDoubleOption BelowValueOption ('\n', "below", "width of the layer below the wetting layer (in Angstrom unit)", 0.0);
  SingleDoubleOption WettingWidthOption ('\n', "wetting", "width of the wetting layer (in Angstrom unit)", 8.0);
  SingleIntegerOption NumberDotOption('\n', "dot", "number of uniformly high layer in the dot", 3);
  SingleDoubleOption BaseRadiusOption ('\n', "base", "base radius in Angstrom unit", 100);
  SingleDoubleOption DotHeightOption ('\n', "height", "height of dot in Angstrom unit", 18);
  SingleDoubleOption TopRadiusOption ('\n', "top", "top radius in Anstrom unit", 74);
  SingleDoubleOption AboveValueOption ('\n', "above", "width of the layer above the dot layer (in Angstrom unit)", 100.0);

  SingleDoubleOption RMassOption ('\n', "mu-r", "electron effective mass in plane (in vacuum electron mass unit)", 0.07);
  SingleDoubleOption ZMassOption ('\n', "mu-z", "electron effective mass in z direction (in vacuum electron mass unit)", 0.07);
  SingleDoubleOption DotPotentialOption ('\n', "dot", "potential in the dot", -0.4);
  SingleDoubleOption MagneticFieldOption ('b', "magnetic", "magnetic field in Z direction (in Tesla unit)", 30);
  SingleDoubleOption WaveVectorOption('w', "wave", "wave vector of Bloch function in Z direction (in 1/Angstrom unit)", 0.0);
  BooleanOption DiskOption ('d', "disk", "enable disk resume capabilities", false);
  BooleanOption ResumeOption ('r', "resume", "resume from disk datas", false);
  SingleIntegerOption VectorMemoryOption ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 400);
  SingleIntegerOption NbrIterationOption ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 60);


  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &EigenstateOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &NumberRValueOption;
  OptionList += &NumberZValueOption;
  OptionList += &NumberMValueOption;
  OptionList += &BelowValueOption;
  OptionList += &WettingWidthOption;
  OptionList += &NumberDotOption;
  OptionList += &BaseRadiusOption;
  OptionList += &DotHeightOption;
  OptionList += &TopRadiusOption;
  OptionList += &AboveValueOption;
  OptionList += &RMassOption;
  OptionList += &ZMassOption;
  OptionList += &DotPotentialOption;
  OptionList += &MagneticFieldOption;
  OptionList += &WaveVectorOption;
  OptionList += &VectorMemoryOption;
  OptionList += &DiskOption;
  OptionList += &ResumeOption;
  OptionList += &NbrIterationOption;

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
  int NbrStateR = NumberRValueOption.GetInteger();
  int NbrStateZ = NumberZValueOption.GetInteger();
  int NumberM = NumberMValueOption.GetInteger();
  double Below = BelowValueOption.GetDouble();
  double WettingWidth = WettingWidthOption.GetDouble();
  double BaseRadius = BaseRadiusOption.GetDouble();
  double DotHeight = DotHeightOption.GetDouble();
  int DotNbr = NumberDotOption.GetInteger();
  double TopRadius = TopRadiusOption.GetDouble();
  double Above = AboveValueOption.GetDouble();
  double Mur = RMassOption.GetDouble();
  double Muz = ZMassOption.GetDouble();
  double DotPotential = DotPotentialOption.GetDouble();
  double MagneticField = MagneticFieldOption.GetDouble();
  double WaveVector= WaveVectorOption.GetDouble();
  int VectorMemory = VectorMemoryOption.GetInteger();
  bool ResumeFlag = ResumeOption.GetBoolean();
  bool DiskFlag = DiskOption.GetBoolean();
  int NbrIterLanczos = NbrIterationOption.GetInteger();
 
  // QuantumDotThreeDConstantCylinderPotential(double belowHeight, double wettingWidth, int nbrCylinderDot, double dotHeight, double baseRadius, double topRadius, double aboveHeight);
  QuantumDotThreeDConstantCylinderPotential* potential = new QuantumDotThreeDConstantCylinderPotential(Below, WettingWidth, DotNbr, DotHeight, BaseRadius, TopRadius, Above);
  // void ConstructPotential(double dotPotential);
  potential->ConstructPotential(DotPotential);

  // define Hilbert space
  // VerticalPeriodicParticleInMagneticField(int nbrStateR, int nbrStateZ, int lowerImpulsionZ);
  VerticalPeriodicParticleInMagneticField* Space = new VerticalPeriodicParticleInMagneticField(NumberM, NbrStateR, NbrStateZ, -NbrStateZ / 2);

  timeval PrecalculationStartingTime;
  timeval PrecalculationEndingTime;
  gettimeofday (&(PrecalculationStartingTime), 0);
  
  //cout << "General space dimension: " << GeneralSpace.GetHilbertSpaceDimension() << endl;
  cout << "Hilbert space Component: " << Space->GetNbrStateR() << '\t' << Space->GetNbrStateZ() << endl;
  cout << "Minimal impulsions:       " << Space->GetLowerImpulsionZ() << endl;
  cout << "Hilbert space dimension: " << Space->GetHilbertSpaceDimension() << endl;
  
  // CylindricalHamiltonianInMagneticField(VerticalPeriodicParticleInMagneticField* space, double mur, double muz, double bz, double waveVector, ThreeDConstantCylinderPotential* PotentialInput);
  CylindricalHamiltonianInMagneticField Hamiltonian(Space, Mur, Muz, MagneticField, WaveVector, potential);
  
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

  double HamiltonianShift = Hamiltonian.MaxPartialDiagonalElement();
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
          OutputFile.open("eigenvalues", ios::binary | ios::out | ios::app);
	  OutputFile << MagneticField << " ";
          for (int i = 0; i < NbrEigenvalue; ++i)
            OutputFile << Eigenvalues[i] << " ";
          OutputFile << endl;
          OutputFile.close();
        }
      /*
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
      */      
      cout << "----------------- End of calculation ---------------------" << endl;      
      cout << "     ==========  CALCULATION IS FINALIZED  =========  " << endl;
    }
  gettimeofday (&(TotalEndingTime), 0);
  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);  
  cout << endl << "Total time = " << Dt << endl;
  delete Lanczos;
  
  return 0;
}
