#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"

#include "Hamiltonian/ExplicitHamiltonian.h"
#include "HilbertSpace/UndescribedHilbertSpace.h"
#include "HilbertSpace/QuantumDotHilbertSpace/Confined3DOneParticle.h"
#include "HilbertSpace/QuantumDotHilbertSpace/Periodic3DOneParticle.h"

#include "Hamiltonian/QuantumDotHamiltonian/PeriodicQuantumDots3DHamiltonian.h"
#include "Hamiltonian/QuantumDotHamiltonian/NewPeriodicQuantumDots3DHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/QuantumDot/Potential/ThreeDPotential.h"
#include "Tools/QuantumDot/Potential/PeriodicPyramidQuantumDot.h"

#include "Tools/QuantumDot/Spectra/PeriodicSpectra.h"

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
  BooleanOption LanczosOption ('l', "lanczos", "enable lanczos diagonalization algorithm", true);
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption VerboseOption ('v', "verbose", "verbose mode", false);
  BooleanOption EigenstateOption ('e', "eigenstate", "evaluate eigenstates", true);
  SingleIntegerOption IterationOption ('i', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 6);
  SingleIntegerOption MValueOption ('M', "M-cell", "number of cells in the x direction", 61);
  SingleIntegerOption NValueOption ('N', "N-cell", "number of cells in the y direction", 61);
  SingleIntegerOption HValueOption ('H', "H-cell", "number of cells in the z direction", 60);
  SingleIntegerOption LeftSizeOption ('\n', "left-size", "size of the leftmost part in the z direction with constant null potential (in cell unit)", 18);
  SingleIntegerOption RightSizeOption ('\n', "right-size", "size of the rightmost part in the z direction with constant potential (in cell unit)", 26);
  SingleDoubleOption CellXSizeOption ('X', "cell-xsize", "cell size in the x direction in Angstrom", 2.97);
  SingleDoubleOption CellYSizeOption ('Y', "cell-ysize", "cell size in the y direction in Angstrom", 2.97);
  SingleDoubleOption CellZSizeOption ('Z', "cell-zsize", "cell size in the z direction in Angstrom", 2.64);  
  SingleDoubleOption XMassOption ('\n', "mu-x", "electron effective mass in x direction (in vacuum electron mass unit)", 0.5045);
  SingleDoubleOption YMassOption ('\n', "mu-y", "electron effective mass in y direction (in vacuum electron mass unit)", 0.5045);
  SingleDoubleOption ZMassOption ('\n', "mu-z", "electron effective mass in z direction (in vacuum electron mass unit)", 1.1);
  SingleStringOption CoefficientFileNameOption('\n', "coefficients", "name of the file where interaction coeffcients are stored", 
					       "/home/regnault/development/DMRG/DiagHam/potentiel_10_10_10_2");
  BooleanOption CarrierTypeOption('c', "carrier", "carrier type, true for hole, false for electron", true);

  List<AbstractOption*> OptionList;
  OptionList += &LanczosOption;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &VerboseOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &EigenstateOption;
  OptionList += &MValueOption;
  OptionList += &NValueOption;
  OptionList += &HValueOption;
  OptionList += &CellXSizeOption;
  OptionList += &CellYSizeOption;
  OptionList += &CellZSizeOption;
  OptionList += &XMassOption;
  OptionList += &YMassOption;
  OptionList += &ZMassOption;
  OptionList += &CoefficientFileNameOption;
  OptionList += &LeftSizeOption;
  OptionList += &RightSizeOption;
  OptionList += &CarrierTypeOption; 

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

  bool LanczosFlag = LanczosOption.GetBoolean();
  bool SMPFlag = SMPOption.GetBoolean();
  bool VerboseFlag = VerboseOption.GetBoolean();
  bool EigenstateFlag = EigenstateOption.GetBoolean();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  char* CoefficientFileName = CoefficientFileNameOption.GetString();
  int M = MValueOption.GetInteger();
  int N = NValueOption.GetInteger();
  int H = HValueOption.GetInteger();
  double Lx = CellXSizeOption.GetDouble();
  double Ly = CellYSizeOption.GetDouble();
  double Lz = CellZSizeOption.GetDouble();
  double Mux = XMassOption.GetDouble();
  double Muy = YMassOption.GetDouble();
  double Muz = ZMassOption.GetDouble();
  int LeftSize = LeftSizeOption.GetInteger();
  int RightSize = RightSizeOption.GetInteger();
  bool Carrier = CarrierTypeOption.GetBoolean();
 
  // **** PROBABILITIES ****
  double p = 0.175;

  // *** Dot geometry ****
  int Rb = 20, Rt = 5, w = 4;

   // **** Offset ****
  double Offset = 0.0;

  // *** Electric field (absolute value) ****;
  double Piezo = 0.03;

  // Pour la simulation du puits quantique InGaAs / GaAs
  Offset = 0.713; Piezo = 0.0; Mux = 0.0715; Muy = 0.0715; Muz = 0.0715;     
  Rb = 0; Rt = 0; w = 4; p = 0.5;

  /*
  if (Carrier)
    {
      Offset = 0.9; Piezo = -0.03;
      Mux = 0.5045; Muy = 0.5045; Muz = 1.1;
    }
  else
    {
      Offset = 1.8; Piezo = 0.03;
      Mux = 0.166; Muy = 0.166; Muz = 0.184;
    }
  */
  // PeriodicPyramidQuantumDotPotential(int NbrCellX, int NbrCellY, int NbrCellZ, double Lz, int u, int a, int rb, int rt, int w, double offset, double concentration, double piezofield, bool scratch, char* logfile)
  PeriodicPyramidQuantumDotPotential* potential = new PeriodicPyramidQuantumDotPotential(M, N, H, Lz, LeftSize, RightSize, Rb, Rt, w, Offset, p, Piezo, Carrier, "DotInput.txt");

  ofstream pfile("DotPotential.txt");
  pfile.precision(14);
  potential->PrintPotential(pfile);
  pfile.close();

  ofstream diagram("Diagram.txt");
  potential->PrintDiagram(diagram);
  diagram.close();  

  // define Hilbert space
  Periodic3DOneParticle Space(M, M / 2, N, N / 2, H / 2 + 1, H / 4);
  timeval PrecalculationStartingTime;
  timeval PrecalculationEndingTime;
  gettimeofday (&(PrecalculationStartingTime), 0);
  
  //  QuantumDots3DHamiltonian(Confined3DOneParticle* space, double xSize, double ySize, double zSize, double mux, double muy, double muz, int nbrCellX, int nbrCellY, int nbrCellZ, ThreeDPotential* PotentialInput, int memory = -1);
 
  NewPeriodicQuantumDots3DHamiltonian Hamiltonian(&Space, Lx * ((double) M), Ly * ((double) N),  Lz * ((double) H), Mux, Muy, Muz, M, N, H, potential);

  gettimeofday (&(PrecalculationEndingTime), 0);
  double Dt = (double) (PrecalculationEndingTime.tv_sec - PrecalculationStartingTime.tv_sec) +
    ((PrecalculationEndingTime.tv_usec - PrecalculationStartingTime.tv_usec) / 1000000.0);
  cout << "precalculation time = " << Dt << endl;

  ofstream Input;
  Input.open("DotInput.txt", ios::out | ios::app);
  Input << "Lattice constants: X = " << Lx << ", Y = " << Ly << ", Z = " << Lz << '\n';
  Input << "Effective masses: Mx = " << Mux << ", My = " << Muy << ", Mz = " << Muz << '\n';
  Input << "Precalculation time: " << Dt << '\n';  
 
  ComplexVector* Eigenstates = 0;
  double* Eigenvalues = 0;
  
  if ((LanczosFlag == false) || (Space.GetHilbertSpaceDimension() < 300))
    {
      if (Hamiltonian.GetHilbertSpaceDimension() > 1)
	{
	  // diagonalize the hamiltonian
	  RealSymmetricMatrix HRep (Hamiltonian.GetHilbertSpaceDimension());
	  Hamiltonian.GetHamiltonian(HRep);
	  RealTriDiagonalSymmetricMatrix TmpTriDiag (Space.GetHilbertSpaceDimension());

	  cout << "start diagonalization..." << endl;
	  if (EigenstateFlag == false)
	    {
	      HRep.Householder(TmpTriDiag, MACHINE_PRECISION);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	    }
	  else
	    {
	      RealMatrix TmpEigenvectors (Space.GetHilbertSpaceDimension(), Space.GetHilbertSpaceDimension());
	      ((RealSymmetricMatrix*) Hamiltonian.GetHamiltonian())->Householder(TmpTriDiag, MACHINE_PRECISION, TmpEigenvectors);
	      TmpTriDiag.Diagonalize(TmpEigenvectors);
	      TmpTriDiag.SortMatrixUpOrder(TmpEigenvectors);
	      Eigenstates = new ComplexVector [NbrEigenvalue];
	      for (int i = 0; i < NbrEigenvalue; ++i)
		{
		  Eigenstates[i] = TmpEigenvectors[i];
		}
	    }

	  // store eigenvalues
	  int Max = Hamiltonian.GetHilbertSpaceDimension();
	  NbrEigenvalue = Max;
	  if (Max > NbrEigenvalue)
	    Max = NbrEigenvalue;
	  Eigenvalues = new double [NbrEigenvalue];
	  for (int j = 0; j < Max ; j++)
	    {
	      Eigenvalues[j]= TmpTriDiag.DiagonalElement(j);
	    }
	}
      else
	{
	  cout << (*(Hamiltonian.GetHamiltonian()))(0, 0) << endl;
	}
    }
  else
    {

      // architecture type (i.e. 1 CPU or multi CPU)
      AbstractArchitecture* Architecture;
      if (SMPFlag == true)
	Architecture = new SMPArchitecture(2);
      else
	Architecture = new MonoProcessorArchitecture;

      double HamiltonianShift = - (150.4 * ((1.0 / (Lx * Lx * Mux)) + (1.0 / (Ly * Ly * Muy)) + (1.0 / (Lz * Lz * Muz))));
      Hamiltonian.ShiftHamiltonian (HamiltonianShift);
      cout << "Décalage:  " << HamiltonianShift << endl;
      // type of lanczos algorithm (with or without reorthogonalization)
      // BasicLanczosAlgorithm Lanczos(Architecture, MaxNbrIterLanczos);
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
	  //cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << " "<< endl;
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
      cout << "diagonalisation time = " << Dt << endl;
    }

  Input << "Diagonalization time: " << Dt << '\n';
  Input << "To verify: M, N, H, under, above = " << M << ", " << N << ", " << H << ", " << LeftSize << ", " << RightSize << endl;
  Input.close(); 

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
