#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/QHEHilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/QHEHilbertSpace/FermionOnSphere.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereCoulombHamiltonian.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereLaplacianDeltaHamiltonian.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereCoulombLaplacianDeltaHamiltonian.h"

#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedLanczosAlgorithmWithDiskStorage.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"
#include "Architecture/SimpleMPIArchitecture.h"

#include "MathTools/ClebschGordanCoefficients.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
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

  OptionManager Manager ("QHEFermionsLaplacianDelta" , "0.01");
  OptionGroup* ParallelizationGroup  = new OptionGroup ("parallelization options");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  Manager += SystemGroup;
  Manager += ParallelizationGroup;
  Manager += LanczosGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 8);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 21);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "add-coulomb", "add coefficent in front of the coulomb pseudo-potentials (pure laplacian delta if 0)", 0.0);
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");
  (*ParallelizationGroup) += new BooleanOption  ('S', "SMP", "enable SMP mode");
  (*ParallelizationGroup) += new SingleIntegerOption  ('\n', "processors", "number of processors to use in SMP mode", 2);
#ifdef __MPI__
  (*ParallelizationGroup) += new BooleanOption  ('\n', "mpi", "enable MPI mode");  
#endif
  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "full-diag", 
					     "maximum Hilbert space dimension for which full diagonalization is applied", 500, true, 100);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup)  += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsLaplacianDelta -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  bool ResumeFlag = ((BooleanOption*) Manager["help"])->GetBoolean();
  bool DiskFlag = ((BooleanOption*) Manager["disk"])->GetBoolean();
  bool GroundFlag = ((BooleanOption*) Manager["ground"])->GetBoolean();
  bool SMPFlag = ((BooleanOption*) Manager["SMP"])->GetBoolean();
#ifdef __MPI__
  bool MPIFlag = ((BooleanOption*) Manager["mpi"])->GetBoolean();
#else
  bool MPIFlag = false;
#endif
  int NbrProcessor = ((SingleIntegerOption*) Manager["processors"])->GetInteger();
  int MaxNbrIterLanczos = ((SingleIntegerOption*) Manager["iter-max"])->GetInteger();
  int NbrIterLanczos = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  int NbrEigenvalue = ((SingleIntegerOption*) Manager["nbr-eigen"])->GetInteger();
  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int FullDiagonalizationLimit = ((SingleIntegerOption*) Manager["full-diag"])->GetInteger();
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) Manager["fast-search"])->GetInteger()) << 20;
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();
  int VectorMemory = ((SingleIntegerOption*) Manager["nbr-vector"])->GetInteger();
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();
  char* SavePrecalculationFileName = ((SingleStringOption*) Manager["save-precalculation"])->GetString();
  double CoulombFactor = ((SingleDoubleOption*) Manager["add-coulomb"])->GetDouble();

  int InvNu = 2;
  double GroundStateEnergy = 0.0;
  double Shift = -10.0;
  char* OutputNameLz = new char [256];
  if (CoulombFactor == 0.0)
    sprintf (OutputNameLz, "fermions_laplaciandelta_n_%d_2s_%d_lz.dat", NbrFermions, LzMax);
  else
    sprintf (OutputNameLz, "fermions_coulomblaplaciandelta_n_%d_2s_%d_%f_lz.dat", NbrFermions, LzMax, CoulombFactor);

  char* OutputNameL = "fermions_l.dat";
  ofstream File;
  File.open(OutputNameLz, ios::binary | ios::out);
  File.precision(14);
  int Max = ((LzMax - NbrFermions + 1) * NbrFermions);
  int TotalSize = 0;
  double** Eigenvalues = new double* [2 * Max + 1];
  int* Dimensions = new int [2 * Max + 1];

  AbstractArchitecture* Architecture = 0;
  if (SMPFlag == false)
    if (MPIFlag == false)
      Architecture = new MonoProcessorArchitecture;
    else
      Architecture = new SimpleMPIArchitecture;
  else
    Architecture = new SMPArchitecture(NbrProcessor);

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
      FermionOnSphere Space (NbrFermions, L, LzMax, MemorySpace);
      cout << " Hilbert space dimension = " << Space.GetHilbertSpaceDimension() << endl;
      TotalSize += Space.GetHilbertSpaceDimension();
      Architecture->SetDimension(Space.GetHilbertSpaceDimension());
      AbstractQHEOnSphereHamiltonian* Hamiltonian;
      if (CoulombFactor == 0.0)
	Hamiltonian = new ParticleOnSphereLaplacianDeltaHamiltonian(&Space, NbrFermions, LzMax, Architecture, Memory, LoadPrecalculationFileName);
      else
	Hamiltonian = new ParticleOnSphereCoulombLaplacianDeltaHamiltonian(&Space, NbrFermions, LzMax, CoulombFactor, Architecture, Memory, LoadPrecalculationFileName);
      Hamiltonian->ShiftHamiltonian(Shift);
      if (SavePrecalculationFileName != 0)
	{
	  Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
	}
      if (Hamiltonian->GetHilbertSpaceDimension() < FullDiagonalizationLimit)
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
		GroundStateEnergy = (TmpTriDiag.DiagonalElement(0) - Shift);
	      for (int j = 0; j < Hamiltonian->GetHilbertSpaceDimension() ; j++)
		{
		  File << (L / 2) << " " << (TmpTriDiag.DiagonalElement(j) - Shift) << endl;
		}
	      cout << endl;
	    }
	  else
	    {
	      File << (L / 2) << " " << (HRep(0, 0) - Shift) << endl;
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
	      Lowest = TmpMatrix.DiagonalElement(NbrEigenvalue - 1) - Shift;
	      Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	      PreviousLowest = Lowest; 
	      cout << (TmpMatrix.DiagonalElement(0) - Shift) << " " << Lowest << " " << Precision << " "<< endl;
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
	  cout << (TmpMatrix.DiagonalElement(0) - Shift) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	       << CurrentNbrIterLanczos << endl;
	  for (int i = 0; i <= NbrEigenvalue; ++i)
	    {
	      cout << TmpMatrix.DiagonalElement(i) << " ";
	      File << (int) (L / 2) << " " << (TmpMatrix.DiagonalElement(i) - Shift) << endl;
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
  delete Architecture;

  return 0;
}
