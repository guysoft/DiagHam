#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/QHEHilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/QHEHilbertSpace/FermionOnSphere.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereCoulombHamiltonian.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereLaplacianDeltaHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithmWithDiskStorage.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include "MathTools/ClebschGordanCoefficients.h"

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
  BooleanOption HelpOption ('h', "help", "display this help");
  BooleanOption SMPOption ('S', "SMP", "enable SMP mode");
  BooleanOption GroundOption ('g', "ground", "restrict to the largest subspace");
  SingleIntegerOption SMPNbrProcessorOption ('\n', "processors", "number of processors to use in SMP mode", 2);
  SingleIntegerOption NbrIterationOption ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  SingleIntegerOption IterationOption ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  SingleIntegerOption NbrEigenvaluesOption ('n', "nbr-eigen", "number of eigenvalues", 30);
  SingleIntegerOption LzMaxOption ('l', "lzmax", "twice the maximum momentum for a single particle", 21);
  SingleIntegerOption InitialLzOption ('\n', "initial-lz", "twice the inital momentum projection for the system", 
					-1);
  SingleIntegerOption NbrFermionOption ('p', "nbr-particles", "number of particles", 8);
  SingleIntegerOption MemoryOption ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  SingleIntegerOption NbrLzOption ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  BooleanOption DiskOption ('d', "disk", "enable disk resume capabilities", true);
  BooleanOption ResumeOption ('r', "resume", "resume from disk datas", false);
  SingleIntegerOption VectorMemoryOption ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  SingleStringOption SavePrecalculationOption ('\n', "save-precalculation", "save precalculation in a file",0);
  SingleStringOption LoadPrecalculationOption ('\n', "load-precalculation", "load precalculation from a file",0);
  List<AbstractOption*> OptionList;
  OptionList += &HelpOption;
  OptionList += &SMPOption;
  OptionList += &GroundOption;
  OptionList += &SMPNbrProcessorOption;
  OptionList += &IterationOption;
  OptionList += &NbrIterationOption;
  OptionList += &NbrEigenvaluesOption;
  OptionList += &NbrFermionOption;
  OptionList += &LzMaxOption;
  OptionList += &MemoryOption;
  OptionList += &InitialLzOption;
  OptionList += &NbrLzOption;
  OptionList += &VectorMemoryOption;
  OptionList += &DiskOption;
  OptionList += &ResumeOption;
  OptionList += &LoadPrecalculationOption;
  OptionList += &SavePrecalculationOption;
  if (ProceedOptions(argv, argc, OptionList) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsLaplacianDelta -h" << endl;
      return -1;
    }
  if (HelpOption.GetBoolean() == true)
    {
      DisplayHelp (OptionList, cout);
      return 0;
    }

  bool ResumeFlag = ResumeOption.GetBoolean();
  bool DiskFlag = DiskOption.GetBoolean();
  bool GroundFlag = GroundOption.GetBoolean();
  bool SMPFlag = SMPOption.GetBoolean();
  int NbrProcessor = SMPNbrProcessorOption.GetInteger();
  int MaxNbrIterLanczos = IterationOption.GetInteger();
  int NbrIterLanczos = NbrIterationOption.GetInteger();
  int NbrEigenvalue = NbrEigenvaluesOption.GetInteger();
  int NbrFermions = NbrFermionOption.GetInteger();
  int LzMax = LzMaxOption.GetInteger();
  long Memory = MemoryOption.GetInteger() << 20;
  int InitialLz = InitialLzOption.GetInteger();
  int NbrLz = NbrLzOption.GetInteger();
  int VectorMemory = VectorMemoryOption.GetInteger();
  char* LoadPrecalculationFileName = LoadPrecalculationOption.GetString();
  char* SavePrecalculationFileName = SavePrecalculationOption.GetString();

  int InvNu = 2;
  double GroundStateEnergy = 0.0;
  int Shift = 0;
  char* OutputNameLz = new char [256];
  sprintf (OutputNameLz, "fermions_laplaciandelta_n_%d_2s_%d_lz.dat", NbrFermions, LzMax);
  char* OutputNameL = "fermions_l.dat";
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);
  int Max = ((LzMax - NbrFermions + 1) * NbrFermions);
  int TotalSize = 0;
  double** Eigenvalues = new double* [2 * Max + 1];
  int* Dimensions = new int [2 * Max + 1];

  int  L = 0;
  if ((abs(Max) & 1) != 0)
     L = 1;
  if (InitialLz >= 0)
    {
      L = InitialLz;
      if ((abs(Max) & 1) != 0)
	L |= 1;
      else
	L &= ~0x1;
    }
  if (GroundFlag == true)
      Max = L;
  else
    {
      if (NbrLz > 0)
	{
	  Max = L + (2 * (NbrLz - 1));
	}
    }
  for (; L <= Max; L += 2)
    {
      cout << "----------------------------------------------------------------" << endl;
      cout << " LzTotal = " << L << endl;
      FermionOnSphere Space (NbrFermions, L, LzMax);
      /*      for (int j = 0; j < 100; ++j)
	for (int i = 0; i < Space.GetHilbertSpaceDimension(); ++i)
	  Space.FindStateIndex(Space.StateDescription[i], Space.StateLzMax[i]);
      for (int i = 0; i <= LzMax; ++i)
	cout << i << " = " << Space.LookUpTableShift[i] << endl;
	//	Space.PrintState(cout, i) << endl;
	return 0;*/
      cout << " Hilbert space dimension = " << Space.GetHilbertSpaceDimension() << endl;
      TotalSize += Space.GetHilbertSpaceDimension();
      AbstractArchitecture* Architecture = 0;
      if (SMPFlag == false)
	Architecture = new MonoProcessorArchitecture;
      else
	Architecture = new SMPArchitecture(NbrProcessor);
      AbstractQHEOnSphereHamiltonian* Hamiltonian;
      Hamiltonian = new ParticleOnSphereLaplacianDeltaHamiltonian(&Space, NbrFermions, LzMax, Architecture, Memory, LoadPrecalculationFileName);
//      Hamiltonian = new ParticleOnSphereCoulombHamiltonian(&Space, NbrFermions, LzMax, Architecture, Memory, LoadPrecalculationFileName);
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
	      if (L == 0)
		GroundStateEnergy = TmpTriDiag.DiagonalElement(0);
	      for (int j = 0; j < Hamiltonian->GetHilbertSpaceDimension() ; j++)
		{
		  File << (L / 2) << " " << TmpTriDiag.DiagonalElement(j) << endl;
		}
	      cout << endl;
	    }
	  else
	    {
	      File << (L / 2) << " " << HRep(0, 0) << endl;
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
	  int CurrentNbrIterLanczos = 0;
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
	  cout << "check" << endl;
	  while ((Lanczos->TestConvergence() == false) &&  (((DiskFlag == true) && (CurrentNbrIterLanczos < NbrIterLanczos)) ||
							    ((DiskFlag == false) && (CurrentNbrIterLanczos < MaxNbrIterLanczos))))
	    {
	      cout << "check in" << endl;
	      ++CurrentNbrIterLanczos;
	      Lanczos->RunLanczosAlgorithm(1);
	      TmpMatrix.Copy(Lanczos->GetDiagonalizedMatrix());
	      TmpMatrix.SortMatrixUpOrder();
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1);
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest; 
	      cout << TmpMatrix.DiagonalElement(0) << " " << Lowest << " " << Precision << " "<< endl;
	      cout << "check out" << endl;
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
      cout << "----------------------------------------------------------------" << endl;
      cout << " Total Hilbert space dimension = " << TotalSize << endl;
      cout << " ground state energy = " << GroundStateEnergy << endl;
      cout << " energy per particle in the ground state = " << (GroundStateEnergy / (double) NbrFermions) << endl;
      delete Hamiltonian;
    }
  File.close();

  return 0;
}
