#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/QHEHilbertSpace/FermionOnDisk.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnDiskLaplacianDeltaHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithmWithDiskStorage.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ios;

double HamiltonianEvaluateInteractionCoefficient(int m1, int m2, int m3, int m4);

int main(int argc, char** argv)
{
  cout.precision(14);
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  SingleIntegerOption SMPNbrProcessorOption ('\n', "processors", "number of processors to use in SMP mode", 2);
  SingleIntegerOption NbrFermionOption ('p', "nbr-particles", "number of particles", 8);
  SingleIntegerOption IterationOption ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 30);
  SingleIntegerOption LzMaxOption ('l', "maximum-momentum", "maximum single particle momentum to study", 10, true, 1);
  SingleIntegerOption LzMinOption ('\n', "minimum-momentum", "minimum single particle momentum to study", 0, true, 1);
  SingleIntegerOption MemoryOption ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  SingleStringOption SavePrecalculationOption ('\n', "save-precalculation", "save precalculation in a file",0);
  SingleStringOption LoadPrecalculationOption ('\n', "load-precalculation", "load precalculation from a file",0);

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &SMPNbrProcessorOption;
  OptionList += &NbrFermionOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &LzMaxOption;
  OptionList += &LzMinOption;
  OptionList += &MemoryOption;
  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsDiskLaplacianDelta -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }
  
  
  bool SMPFlag = SMPOption.GetBoolean();
  int NbrProcessor = SMPNbrProcessorOption.GetInteger();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  int NbrFermions = NbrFermionOption.GetInteger();
  int MMin = LzMinOption.GetInteger();
  int MMax = LzMaxOption.GetInteger();
  if (MMin < (((NbrFermions - 1) * (NbrFermions)) / 2))
    MMin = (((NbrFermions - 1) * (NbrFermions)) / 2);
  if (MMax < MMin)
    MMax = MMin;
  long Memory = MemoryOption.GetInteger() << 20;
  char* LoadPrecalculationFileName = LoadPrecalculationOption.GetString();
  char* SavePrecalculationFileName = SavePrecalculationOption.GetString();

  char* OutputName = new char [1024];
  sprintf (OutputName, "fermions_disk_laplaciandelta_n_%d_l_%d.dat", NbrFermions, MMax);
  ofstream File;
  File.open(OutputName, ios::binary | ios::out);

  for (int  L = MMin; L <= MMax; ++L)
    {
      AbstractArchitecture* Architecture = 0;
      if (SMPFlag == false)
	Architecture = new MonoProcessorArchitecture;
      else
	Architecture = new SMPArchitecture(NbrProcessor);
      int MaxMomentum = L - (((NbrFermions - 1) * (NbrFermions - 2)) / 2);
      cout << "MaxMomentum=" << MaxMomentum << endl;
      int m4 = 1;
      for (int m1 = 0; m1 <= MaxMomentum; ++m1)
	for (int m2 = 2; m2 < m1; ++m2)
	  for (int m3 = 2; m3 <= 2; ++m3)
	    {
	      cout << m1 << " " << m2 << " " << HamiltonianEvaluateInteractionCoefficient(m1, m2, m3, m4) << endl;
	    }
/*	  for (int m3 = 0; m3 <= MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if ((m4 >= 0) && (m3 > m4))
		HamiltonianEvaluateInteractionCoefficient(m1, m2, m3, m4);
	    }*/
      return 0;
      FermionOnDisk Space(NbrFermions, L);
      cout << "Nbr fermions = " << NbrFermions << "    L = " << L << "    Dimension = " << Space.GetHilbertSpaceDimension() << endl;
      ParticleOnDiskLaplacianDeltaHamiltonian* Hamiltonian = new ParticleOnDiskLaplacianDeltaHamiltonian(&Space, NbrFermions, Architecture, Memory, LoadPrecalculationFileName);
      if (SavePrecalculationFileName != 0)
	{
	  Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
	}
      if (Hamiltonian->GetHilbertSpaceDimension() < 500)
	{
	  RealSymmetricMatrix HRep (Hamiltonian->GetHilbertSpaceDimension());
	  Hamiltonian->GetHamiltonian(HRep);
	  if (Hamiltonian->GetHilbertSpaceDimension() > 1)
	    {
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian->GetHilbertSpaceDimension());
	      HRep.Householder(TmpTriDiag, 1e-7);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	      for (int j = 0; j < Hamiltonian->GetHilbertSpaceDimension(); j++)
		{
		  File << L << " " << TmpTriDiag.DiagonalElement(j) << endl;
		}
	    }
	  else
	    {
	      double TmpVal = HRep(0, 0);
	      File << L << " " << TmpVal << endl;
	    }
	}
      else
	{
	  double Shift = -10.0;
	  Hamiltonian->ShiftHamiltonian(Shift);
	  AbstractLanczosAlgorithm* Lanczos;
	  if (NbrEigenvalue == 1)
	    {
	       Lanczos = new BasicLanczosAlgorithm(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	    }
	  else
	    {
	      Lanczos = new FullReorthogonalizedLanczosAlgorithm (Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	    }
	  double Precision = 1.0;
	  double PreviousLowest = 1e50;
	  double Lowest = PreviousLowest;
	  int CurrentNbrIterLanczos = 0;
	  Lanczos->SetHamiltonian(Hamiltonian);
	  Lanczos->InitializeLanczosAlgorithm();
	  cout << "Run Lanczos Algorithm" << endl;
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  double Dt;
	  gettimeofday (&(TotalStartingTime), 0);
	  Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
	  CurrentNbrIterLanczos = NbrEigenvalue + 3;
	  RealTriDiagonalSymmetricMatrix TmpMatrix;
	  while ((Lanczos->TestConvergence() == false) &&  (CurrentNbrIterLanczos < MaxNbrIterLanczos))
	    {
	      ++CurrentNbrIterLanczos;
	      Lanczos->RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest; 
	      cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << " "<< endl;
	    }
	  if (CurrentNbrIterLanczos >= MaxNbrIterLanczos)
	    {
	      cout << "too much Lanczos iterations" << endl;
	      File << "too much Lanczos iterations" << endl;
	      File.close();
	      exit(0);
	    }
	  for (int i = 0; i <= NbrEigenvalue; ++i)
	    {
	      cout << (TmpMatrix.DiagonalElement(i) - Shift) << " ";
	      File << L << " " << (TmpMatrix.DiagonalElement(i) - Shift) << endl;
	    }
	  cout << endl;
	  gettimeofday (&(TotalEndingTime), 0);
	  cout << "------------------------------------------------------------------" << endl << endl;;
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	  delete Lanczos;
	}
      cout << endl;
      cout << "//////////////////////////////////////////////////////" << endl;
    }
  File.close();
  return 0;
}


double HamiltonianEvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  if ((m1 == m2) || (m3 == m4))
    return 0.0;
  FactorialCoefficient Coef;
  Coef.SetToOne();
  if (m2 > 1)
    {
      Coef.PartialFactorialMultiply(m1 + 1, m1 + m2 - 1);
      Coef.FactorialDivide(m2);
    }
  else
    {
      if (m2 == 0)
	Coef /= m1;	
    }
  if (m4 > 1)
    {
      Coef.PartialFactorialMultiply(m3 + 1, m3 + m4 - 1);
      Coef.FactorialDivide(m4);
    }
  else
    {
      if (m4 == 0)
	Coef /= m3;	
    }
  Coef.Power2Divide(2 * (m1 + m2));
  double Val1 = (sqrt(Coef.GetNumericalValue()) * ((double) ((m2 - m1) * (m3 - m4)))/ M_PI);

  FactorialCoefficient Coef2;
  Coef2.SetToOne();
  if (m2 > 1)
    {
      Coef2.PartialFactorialMultiply(m1 + 1, m1 + m2 - 1);
      Coef2.FactorialDivide(m2);
    }
  else
    {
      if (m2 == 0)
	Coef2 /= m1;	
    }
  Coef2.Power2Divide(m1 + m2);
  double Val2 = sqrt (Coef2.GetNumericalValue());
  return Val2;
/*  Coef2.SetToOne();
  if (m4 > 1)
    {
      Coef2.PartialFactorialMultiply(m3 + 1, m3 + m4 - 1);
      Coef2.FactorialDivide(m4);
    }
  else
    {
      if (m4 == 0)
	Coef2 /= m3;	
    }
  Coef2.Power2Divide(m3 + m4);
  Val2 *= sqrt (Coef2.GetNumericalValue()) * (((double) ((m2 - m1) * (m3 - m4)))/ M_PI);
  if (fabs(Val1 - Val2) > (1e-14 * fabs(Val1)))
    {
      cout << "error ";
    }
  cout << m1 << " "  << m2 << " "  << m3 << " "  << m4 << " " << Val1 << " " << Val2 << endl;
  return Val1;*/
}

