#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/QHEHilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/QHEHilbertSpace/FermionOnSphere.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereDeltaHamiltonian.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereDeltaModifiedHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithmWithDiskStorage.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "MathTools/ClebschGordanCoefficients.h"
#include "BitmapPicture/AbstractBitmapPicture.h"

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


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption DiskOption ('d', "disk", "enable disk resume capabilities", false);
  BooleanOption ResumeOption ('r', "resume", "resume from disk datas", false);
  SingleIntegerOption SMPNbrProcessorOption ('\n', "processors", "number of processors to use in SMP mode", 2);
  SingleIntegerOption IterationOption ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrIterationOption ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 40);
  SingleIntegerOption LzMaxOption ('l', "lzmax", "twice the maximum momentum for a single particle", 14);
  SingleIntegerOption MinTotalLzOption ('\n', "min-lz", "minimum total lz value to evaluate (negative if it is 0 or 1 depending on the parameters)", -1);
  SingleIntegerOption MaxTotalLzOption ('\n', "max-lz", "maximum total lz value to evaluate (negative if all values have to be evaluated)", -1);
  SingleIntegerOption NbrBosonOption ('p', "nbr-particles", "number of particles", 8);
  SingleIntegerOption MemoryOption ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  SingleIntegerOption VectorMemoryOption ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  SingleStringOption SavePrecalculationOption ('\n', "save-precalculation", "save precalculation in a file",0);
  SingleStringOption LoadPrecalculationOption ('\n', "load-precalculation", "load precalculation from a file",0);
  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &SMPNbrProcessorOption;
  OptionList += &IterationOption;
  OptionList += &NbrIterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &NbrBosonOption;
  OptionList += &LzMaxOption;
  OptionList += &MinTotalLzOption;
  OptionList += &MaxTotalLzOption;
  OptionList += &MemoryOption;
  OptionList += &VectorMemoryOption;
  OptionList += &DiskOption;
  OptionList += &ResumeOption;
  OptionList += &LoadPrecalculationOption;
  OptionList += &SavePrecalculationOption;
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
  bool ResumeFlag = ResumeOption.GetBoolean();
  bool DiskFlag = DiskOption.GetBoolean();
  bool SMPFlag = SMPOption.GetBoolean();
  int NbrProcessor = SMPNbrProcessorOption.GetInteger();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrIterLanczos = NbrIterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  int NbrBosons = NbrBosonOption.GetInteger();
  int LzMax = LzMaxOption.GetInteger();
  int MinTotalLz = MinTotalLzOption.GetInteger();
  int MaxTotalLz = MaxTotalLzOption.GetInteger();
  int Memory = MemoryOption.GetInteger() << 20;
  int VectorMemory = VectorMemoryOption.GetInteger();
  char* LoadPrecalculationFileName = LoadPrecalculationOption.GetString();
  char* SavePrecalculationFileName = SavePrecalculationOption.GetString();

  int InvNu = 2;
  double GroundStateEnergy = 0.0;
  int Shift = 0;
  char* OutputNameLz = new char [256];
  sprintf (OutputNameLz, "bosons_delta_n_%d_2s_%d_lz.dat", NbrBosons, LzMax);
  char* OutputNameL = "bosons_l.dat";
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);
  int Max = (LzMax * NbrBosons);
  if (MaxTotalLz >= 0)
    Max = MaxTotalLz;
  int TotalSize = 0;
  double** Eigenvalues = new double* [2 * Max + 1];
  int* Dimensions = new int [2 * Max + 1];
  int  L = MinTotalLz;
  if (L < 0)
    {
      if ((abs(Max) & 1) != 0)
	L = 1;
      else
	L = 0;
    }
  for (; L <= Max; L += 2)
    {
      cout << "----------------------------------------------------------------" << endl;
      cout << " LzTotal = " << L << endl;
      BosonOnSphere Space (NbrBosons, L, LzMax);
//      FermionOnSphere Space (NbrBosons, L, LzMax);
      cout << " Hilbert space dimension = " << Space.GetHilbertSpaceDimension() << endl;
      TotalSize += Space.GetHilbertSpaceDimension();
/*     for (int i = 0; i < Space.GetHilbertSpaceDimension(); ++i)
	{
	  cout << i << " = ";
	  Space.PrintState(cout, i) << endl;
	}*/
      AbstractArchitecture* Architecture = 0;
      if (SMPFlag == false)
	Architecture = new MonoProcessorArchitecture;
      else
	Architecture = new SMPArchitecture(NbrProcessor);
      ParticleOnSphereDeltaHamiltonian* Hamiltonian = new ParticleOnSphereDeltaHamiltonian(&Space, NbrBosons, LzMax, Architecture, Memory, LoadPrecalculationFileName);
      if (SavePrecalculationFileName != 0)
	{
	  Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
	}

      if (Hamiltonian->GetHilbertSpaceDimension() < 300)
	{
//	  Dimensions[L >> 1] = Hamiltonian->GetHilbertSpaceDimension();
//	  Eigenvalues[L >> 1] = new double [Hamiltonian->GetHilbertSpaceDimension()];
	  RealSymmetricMatrix HRep (Hamiltonian->GetHilbertSpaceDimension());
	  Hamiltonian->GetHamiltonian(HRep);
/*	  AbstractBitmapPicture* TmpPic = Hamiltonian->GetHamiltonianColorPicture(1e-10);
	  char* TmpPicName = new char [256];
	  sprintf (TmpPicName, "pic_lz_%d.tga", L);
	  TmpPic->SavePicture(TmpPicName);*/
//	  cout << HRep << endl;
	  if (Hamiltonian->GetHilbertSpaceDimension() > 1)
	    {
	      RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian->GetHilbertSpaceDimension());
	      HRep.Householder(TmpTriDiag, 1e-7);
	      TmpTriDiag.Diagonalize();
	      TmpTriDiag.SortMatrixUpOrder();
	      if (L == 0)
		GroundStateEnergy = TmpTriDiag.DiagonalElement(0);
	      //	  cout << "eigenvalues : " << endl;
	      for (int j = 0; j < Hamiltonian->GetHilbertSpaceDimension() ; j++)
		{
//		  Eigenvalues[L >> 1][j] = TmpTriDiag.DiagonalElement(j);
//		  cout << TmpTriDiag.DiagonalElement(j) << " ";
		  File << (L / 2) << " " << TmpTriDiag.DiagonalElement(j) << endl;
// (TmpTriDiag.DiagonalElement(j) - GroundStateEnergy) << endl;
		}
	      cout << endl;
	    }
	  else
	    {
//	      Eigenvalues[L >> 1][0] = HRep(0, 0);
//	      cout << HRep(0, 0) << endl;
	      //	      GroundStateEnergy = HRep(0, 0);
	      File << (L / 2) << " " << HRep(0, 0) << endl;// - GroundStateEnergy) / (4 * M_PI)) << endl;
	    }
	}
      else
	{
	  int MaxNbrIterLanczos = 4000;
	  AbstractLanczosAlgorithm* Lanczos;
	  if (NbrEigenvalue == 1)
	    {
	      if (DiskFlag == false)
		Lanczos = new BasicLanczosAlgorithm(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	      else
		Lanczos = new BasicLanczosAlgorithmWithDiskStorage(Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	    }
	  else
	    {
	      if (DiskFlag == false)
		Lanczos = new FullReorthogonalizedLanczosAlgorithm (Architecture, NbrEigenvalue, MaxNbrIterLanczos);
	      else
		Lanczos = new FullReorthogonalizedLanczosAlgorithmWithDiskStorage (Architecture, NbrEigenvalue, VectorMemory, MaxNbrIterLanczos);
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
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue);
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
	      File << (L / 2) << " " << TmpMatrix.DiagonalElement(i) << endl;
	    }
	  cout << endl;
	  gettimeofday (&(TotalEndingTime), 0);
	  cout << "------------------------------------------------------------------" << endl << endl;;
	  Dt = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
	    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
	  cout << "time = " << Dt << endl;
	}
      cout << "----------------------------------------------------------------" << endl;
      cout << " Total Hilbert space dimension = " << TotalSize << endl;
      cout << " ground state energy = " << GroundStateEnergy << endl;
      cout << " energy per particle in the ground state = " << (GroundStateEnergy / (double) NbrBosons) << endl;
      delete Hamiltonian;
    }
  File.close();

  return 0;
}
