#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/QHEHilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/QHEHilbertSpace/FermionOnTorusWithMagneticTranslations.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnTorusCoulombHamiltonian.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian.h"

#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  SingleIntegerOption SMPNbrProcessorOption ('\n', "processors", "number of processors to use in SMP mode", 2);
  SingleIntegerOption IterationOption ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrIterationOption ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 40);
  BooleanOption GroundOption ('g', "ground", "restrict to the largest subspace");
  SingleIntegerOption NbrFermionOption ('p', "nbr-particles", "number of particles", 6);
  SingleIntegerOption MaxMomentumOption ('l', "max-momentum", "maximum momentum for a single particle", 9);
  SingleIntegerOption XMomentumOption ('x', "x-momentum", "constraint on the total momentum in the x direction (negative if none)", -1);
  SingleIntegerOption YMomentumOption ('y', "y-momentum", "constraint on the total momentum in the y direction (negative if none)", -1);
  SingleDoubleOption XRatioOption ('r', "ratio", "ratio between lengths along the x and y directions (-1 if has to be taken equal to nbr-particles/4)", -1);
  SingleIntegerOption MaxFullDiagonalizationOption ('f', "max-full", "maximum hilbert space size allowed to use full diagonalization", 300);
  BooleanOption DiskOption ('\n', "disk", "enable disk resume capabilities", false);
  BooleanOption ResumeOption ('\n', "resume", "resume from disk datas", false);
  SingleIntegerOption VectorMemoryOption ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  SingleStringOption SavePrecalculationOption ('\n', "save-precalculation", "save precalculation in a file",0);
  SingleStringOption LoadPrecalculationOption ('\n', "load-precalculation", "load precalculation from a file",0);

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &GroundOption;
  OptionList += &SMPNbrProcessorOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &NbrFermionOption;
  OptionList += &MaxMomentumOption;
  OptionList += &XMomentumOption;
  OptionList += &YMomentumOption;
  OptionList += &XRatioOption;
  OptionList += &MaxFullDiagonalizationOption;
  OptionList += &NbrIterationOption;
  OptionList += &VectorMemoryOption;
  OptionList += &DiskOption;
  OptionList += &ResumeOption;
  OptionList += &LoadPrecalculationOption;
  OptionList += &SavePrecalculationOption;
  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsTorus -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }


  bool GroundFlag = GroundOption.GetBoolean();
  bool SMPFlag = SMPOption.GetBoolean();
  int NbrProcessor = SMPNbrProcessorOption.GetInteger();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrIterLanczos = NbrIterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  int NbrFermions = NbrFermionOption.GetInteger();
  int MaxMomentum = MaxMomentumOption.GetInteger();
  int XMomentum = XMomentumOption.GetInteger();
  int YMomentum = YMomentumOption.GetInteger();
  int MaxFullDiagonalization = MaxFullDiagonalizationOption.GetInteger();
  double XRatio = NbrFermions / 4.0;
  if (XRatioOption.GetDouble() > 0)
    {
       XRatio = XRatioOption.GetDouble();
    }
  bool ResumeFlag = ResumeOption.GetBoolean();
  bool DiskFlag = DiskOption.GetBoolean();
  int VectorMemory = VectorMemoryOption.GetInteger();
  char* LoadPrecalculationFileName = LoadPrecalculationOption.GetString();
  char* SavePrecalculationFileName = SavePrecalculationOption.GetString();

  int InvNu = 3;
//  int MaxMomentum = InvNu * NbrFermions;
  int L = 0;
  double GroundStateEnergy = 0.0;

  char* OutputNameLz = new char [512];
  sprintf (OutputNameLz, "fermions_torus_coulomb_n_%d_2s_%d_ratio_%f.dat", NbrFermions, MaxMomentum, XRatio);
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);


  
  AbstractArchitecture* Architecture = 0;
  if (SMPFlag == false)
    Architecture = new MonoProcessorArchitecture;
  else
    Architecture = new SMPArchitecture(NbrProcessor);

  int MomentumModulo = FindGCD(NbrFermions, MaxMomentum);
  int XMaxMomentum = (MomentumModulo - 1);
  if (XMomentum < 0)
    XMomentum = 0;
  else
    XMaxMomentum = XMomentum;
//  int YMaxMomentum = (MomentumModulo - 1);
  int YMaxMomentum = (MaxMomentum - 1);
  if (YMomentum < 0)
    YMomentum = 0;
  else
    YMaxMomentum = YMomentum; 

  for (; XMomentum <= XMaxMomentum; ++XMomentum)
    for (int YMomentum2 = YMomentum; YMomentum2 <= YMaxMomentum; ++YMomentum2)
      {     
	cout << "----------------------------------------------------------------" << endl;
	cout << " Ratio = " << XRatio << endl;
//	FermionOnTorus TotalSpace (NbrFermions, MaxMomentum, y);
	FermionOnTorusWithMagneticTranslations TotalSpace (NbrFermions, MaxMomentum, XMomentum, YMomentum2);
	cout << " Total Hilbert space dimension = " << TotalSpace.GetHilbertSpaceDimension() << endl;
	//      cout << "momentum = " << Momentum << endl;
	cout << "momentum = (" << XMomentum << "," << YMomentum2 << ")" << endl;
//	cout << "momentum = (" << y << ")" << endl;
/*	for (int i = 0; i < TotalSpace.GetHilbertSpaceDimension(); ++i)
	  {
	    cout << i << " = ";
	    TotalSpace.PrintState(cout, i) << endl;
	  }
	cout << endl << endl;*/

/*      for (int i = 0; i < TotalSpace.GetHilbertSpaceDimension(); ++i)
	{
	  cout << "---------------------------------------------" << endl;
	  cout << i << " = " << endl;;
	  for (int m1 = 0; m1 < MaxMomentum; ++m1)
	    for (int m2 = 0; m2 < m1; ++m2)
	      for (int m3 = 0; m3 < MaxMomentum; ++m3)
		{
		  int m4 = m1 + m2 - m3;
		  if (m4 < 0)
		    m4 += MaxMomentum;
		  else
		    if (m4 >= MaxMomentum)
		      m4 -= MaxMomentum;
		  if (m3 > m4)
		    {
		      double Coefficient = 0.0;
		      int NbrTranslations = 0;
		      TotalSpace.AdAdAA(i, m1, m2, m3, m4, Coefficient, NbrTranslations);
		    }
		}
		}*/
	  
      AbstractArchitecture* Architecture = 0;
      if (SMPFlag == false)
	Architecture = new MonoProcessorArchitecture;
      else
	Architecture = new SMPArchitecture(NbrProcessor);
//      AbstractHamiltonian* Hamiltonian = new ParticleOnTorusCoulombHamiltonian (&TotalSpace, NbrFermions, MaxMomentum, XRatio);
      AbstractHamiltonian* Hamiltonian = new ParticleOnTorusCoulombWithMagneticTranslationsHamiltonian (&TotalSpace, NbrFermions, MaxMomentum, XMomentum, XRatio,
													Architecture);
      if (Hamiltonian->GetHilbertSpaceDimension() < MaxFullDiagonalization)
	{
	  HermitianMatrix HRep2 (Hamiltonian->GetHilbertSpaceDimension());
	  Hamiltonian->GetHamiltonian(HRep2);
	  Complex Zero;
//	  HRep2.SetMatrixElement(0, 1, Zero);
//	  HRep2.SetMatrixElement(1, 5, Zero);
//	  cout << HRep2 << endl;	  
	  RealSymmetricMatrix HRep (HRep2.ConvertToSymmetricMatrix());
	  if (Hamiltonian->GetHilbertSpaceDimension() > 1)
	    {
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian->GetHilbertSpaceDimension(), true);
	      HRep.Householder(TmpTriDiag, 1e-7);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	      if (L == 0)
		GroundStateEnergy = TmpTriDiag.DiagonalElement(0);
	      for (int j = 0; j < Hamiltonian->GetHilbertSpaceDimension() ; j++)
		{
		  File << XMomentum << " " << YMomentum2 << " " << TmpTriDiag.DiagonalElement(2 * j) << endl;
		  cout << XMomentum << " " << YMomentum2 << " " << TmpTriDiag.DiagonalElement(2 * j) << endl;
		}
	      cout << endl;
	    }
	  else
	    {
	      File << XMomentum << " " << YMomentum2 << " " << HRep(0, 0) << endl;
	    }
	}
      else
	{
	  int MaxNbrIterLanczos = 4000;
	  AbstractLanczosAlgorithm* Lanczos;
	  if (NbrEigenvalue == 1)
	    {
	      if (DiskFlag == false)
		Lanczos = new ComplexBasicLanczosAlgorithm(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	      else
		Lanczos = new ComplexBasicLanczosAlgorithmWithDiskStorage(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	    }
	  else
	    {
	      if (DiskFlag == false)
		Lanczos = new FullReorthogonalizedComplexLanczosAlgorithm (Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	      else
		Lanczos = new FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage (Architecture, NbrEigenvalue, VectorMemory, MaxNbrIterLanczos);
	    }
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = NbrEigenvalue + 3;
	  Lanczos->SetHamiltonian(Hamiltonian);
	  if ((DiskFlag == true) && (ResumeFlag == true))
	    Lanczos->ResumeLanczosAlgorithm();
	  else
	    Lanczos->InitializeLanczosAlgorithm();
	  cout << "Run Lanczos Algorithm" << endl;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  double Dt;
	  gettimeofday (&(TotalStartingTime), 0);
	  if (ResumeFlag == false)
	    {
	      Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
	      CurrentNbrIterLanczos = NbrEigenvalue + 3;
	    }
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  while ((Lanczos->TestConvergence() == false) && (((DiskFlag == true) && (CurrentNbrIterLanczos < NbrIterLanczos)) ||
							   ((DiskFlag == false) && (CurrentNbrIterLanczos < MaxNbrIterLanczos))))
	    {
	      Lanczos->RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest; 
	      cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << " "<< endl;
	      ++CurrentNbrIterLanczos;
	    }
	  if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
	    {
	      cout << "too much Lanczos iterations" << endl;
	      File << "too much Lanczos iterations" << endl;
	      File.close();
	      exit(0);
	    }
	  GroundStateEnergy = Lowest;
	  cout << endl;
	  cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	       << CurrentNbrIterLanczos << endl;
	  for (int i = 0; i <= NbrEigenvalue; ++i)
	    {
	      cout << TmpMatrix.DiagonalElement(i) << " ";
	      File << XMomentum << " " << YMomentum2 << " " << TmpMatrix.DiagonalElement(i) << endl;
	    }
	  cout << endl;
	  gettimeofday (&(TotalEndingTime), 0);
	  cout << "------------------------------------------------------------------" << endl << endl;;
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	}
      cout << "----------------------------------------------------------------" << endl;
      cout << " ground state energy = " << GroundStateEnergy << endl;
      cout << " energy per particle in the ground state = " << (GroundStateEnergy / (double) NbrFermions) << endl;
      delete Hamiltonian;

    }
  File.close();

  return 0;
}
