#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/QHEHilbertSpace/TrappedBosons.h"
#include "Hamiltonian/QHEHamiltonian/TrappedBosonHamiltonian.h"

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

int main(int argc, char** argv)
{
  cout.precision(14);
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  SingleIntegerOption SMPNbrProcessorOption ('\n', "processors", "number of processors to use in SMP mode", 2);
  SingleIntegerOption NbrBosonOption ('p', "nbr-particles", "number of particles", 8);
  SingleIntegerOption IterationOption ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 30);
  SingleIntegerOption LzMaxOption ('l', "maximum-momentum", "maximum single particle momentum to study", 10, true, 1);
  SingleIntegerOption LzMinOption ('\n', "minimum-momentum", "minimum single particle momentum to study", 1, true, 1);

  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &SMPNbrProcessorOption;
  OptionList += &NbrBosonOption;
  OptionList += &IterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &LzMaxOption;
  OptionList += &LzMinOption;
  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsDiskDelta -h" << endl;
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
  int NbrBosons = NbrBosonOption.GetInteger();
  int MMin = LzMinOption.GetInteger();
  int MMax = LzMaxOption.GetInteger();
  if (MMax < MMin)
    MMax = MMin;

  char* OutputName = new char [1024];
  sprintf (OutputName, "bosons_disk_delta_n_%d_l_%d.dat", NbrBosons, MMax);
  ofstream File;
  File.open(OutputName, ios::binary | ios::out);

  for (int  L = MMin; L <= MMax; ++L)
    {
      AbstractArchitecture* Architecture = 0;
      if (SMPFlag == false)
	Architecture = new MonoProcessorArchitecture;
      else
	Architecture = new SMPArchitecture(NbrProcessor);
      TrappedBosons Space (NbrBosons, L);
      cout << "Nbr bosons = " << NbrBosons << "    L = " << L << "    Dimension = " << Space.GetHilbertSpaceDimension() << endl;
      TrappedBosonHamiltonian* Hamiltonian = new TrappedBosonHamiltonian(&Space, L);
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
	      cout << TmpMatrix.DiagonalElement(i) << " ";
	      File << (int) (L / 2) << " " << (TmpMatrix.DiagonalElement(i)) << endl;
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

