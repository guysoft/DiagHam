#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/RealVector.h"
#include "HilbertSpace/SpinHilbertSpace/Spin1ChainWithTranslations.h"
#include "HilbertSpace/SpinHilbertSpace/Spin1Chain.h"
#include "Hamiltonian/SpinHamiltonian/PeriodicSpinChainHamiltonian.h"
#include "Hamiltonian/SpinHamiltonian/SpinChainHamiltonianWithTranslations.h"
#include "GeneralTools/List.h"
#include "GeneralTools/ListIterator.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

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


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);
 
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  SingleIntegerOption NbrSpinOption ('\0', "nbr-spin", "number of spins", 8);
  SingleIntegerOption SMPNbrProcessorOption ('\n', "processors", "number of processors to use in SMP mode", 2);
  SingleIntegerOption IterationOption ('\n', "iter-max", "maximum number of lanczos iteration (including resume run)", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 1);
  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &SMPNbrProcessorOption;
  OptionList += &IterationOption;
  OptionList += &NbrSpinOption;
  OptionList += &NbrEigenvaluesOption;
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
  int NbrSpin = NbrSpinOption.GetInteger();
  bool SMPFlag = SMPOption.GetBoolean();
  int NbrProcessor = SMPNbrProcessorOption.GetInteger();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
/*  for (int i = 0; i < NbrSpinPerChain; ++i)
    {
      Spin1ChainWithTranslations Chain(NbrSpin, i, 0, 10000000, 10000000);
      cout << "--------------------------------------------" << endl;
      cout << "momentum = " << i << endl;
      cout << "dimension = " << Chain.GetHilbertSpaceDimension() << endl;
//      for (int i = 0; i < Chain.GetHilbertSpaceDimension(); ++i)
//	Chain.PrintState(cout, i) << endl;
    }*/
  double* CouplingConstants = new double[NbrSpin + 1];
  for (int i = 0; i <= NbrSpin; ++i)
    CouplingConstants[i] = 1.0;
//  Spin1Chain Space(NbrSpin, 0, 1000000);
  for (int k = 0; k < NbrSpin; ++k)
    {
//      cout << "momentum = " << k << endl;
  Spin1ChainWithTranslations Space(NbrSpin, k, 0, 10000000, 10000000);
  AbstractArchitecture* Architecture = 0;
  if (SMPFlag == false)
    Architecture = new MonoProcessorArchitecture;
  else
    Architecture = new SMPArchitecture(NbrProcessor);
//  PeriodicSpinChainHamiltonian Hamiltonian(&Space, NbrSpin, CouplingConstants);
  SpinChainHamiltonianWithTranslations Hamiltonian(&Space, NbrSpin, 1.0);
  if (Hamiltonian.GetHilbertSpaceDimension() < 200)
    {
      HermitianMatrix HRep2 (Hamiltonian.GetHilbertSpaceDimension());
      Hamiltonian.GetHamiltonian(HRep2);
//      cout << HRep2 << endl;
      RealSymmetricMatrix HRep (HRep2.ConvertToSymmetricMatrix());
/*      RealSymmetricMatrix HRep (Hamiltonian.GetHilbertSpaceDimension());
      Hamiltonian.GetHamiltonian(HRep);*/
      if (Hamiltonian.GetHilbertSpaceDimension() > 1)
	{
	  RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian.GetHilbertSpaceDimension());
	  HRep.Householder(TmpTriDiag, 1e-7);
	  TmpTriDiag.Diagonalize();
	  TmpTriDiag.SortMatrixUpOrder();
	  for (int j = 0; j < HRep.GetNbrRow() ; j += 2)
	    {
	     cout << TmpTriDiag.DiagonalElement(j) << " ";
	    }
	  cout << endl;
	}
      else
	{
	  cout << HRep(0, 0) << endl;
	}
    }
  else
    {
      int MaxNbrIterLanczos = 4000;
      AbstractLanczosAlgorithm* Lanczos;
      if (NbrEigenvalue == 1)
	{
//	  Lanczos = new BasicLanczosAlgorithm(Architecture, NbrEigenvalue, MaxNbrIterLanczos);	
	  Lanczos = new ComplexBasicLanczosAlgorithm(Architecture, NbrEigenvalue, MaxNbrIterLanczos);	
	}
      else
	{
//	  Lanczos = new FullReorthogonalizedLanczosAlgorithm (Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	  Lanczos = new FullReorthogonalizedComplexLanczosAlgorithm(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	}
      double Precision = 1.0;
      double PreviousLowest = 1e50;
      double Lowest = PreviousLowest;
      int CurrentNbrIterLanczos = 0;
      Lanczos->SetHamiltonian(&Hamiltonian);
      Lanczos->InitializeLanczosAlgorithm();
      cout << "Run Lanczos Algorithm" << endl;
      timeval TotalStartingTime;
      timeval TotalEndingTime;
      double Dt;
      gettimeofday (&(TotalStartingTime), 0);
      Lanczos->RunLanczosAlgorithm(NbrEigenvalue + 2);
      CurrentNbrIterLanczos = NbrEigenvalue + 3;
      RealTriDiagonalSymmetricMatrix TmpMatrix;
      while (Lanczos->TestConvergence() == false)      
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
	  exit(0);
	}
      cout << endl;
      cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	 << CurrentNbrIterLanczos << endl;
      for (int i = 0; i < NbrEigenvalue; ++i)
	{
	  cout << TmpMatrix.DiagonalElement(i) << " ";
	}
      cout << endl;
    }
    }
  return 0;
}

