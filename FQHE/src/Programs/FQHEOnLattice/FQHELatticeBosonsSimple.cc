#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"

#include "Vector/ComplexVector.h"

#include "HilbertSpace/BosonOnLattice.h"

#include "Hamiltonian/ParticleOnLatticeDeltaHamiltonian.h"

#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Options/Options.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEBosonsLatticeSimple" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 8);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 5);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 1);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice (-1=all)", -1);
  (*SystemGroup) += new SingleDoubleOption  ('u', "contactU", "prefactor U of the contact interaction (kinetic term ~ 1)", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('d', "deltaPotential", "Introduce a delta-potential at the origin", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "negative-hopping", "reverse sign of hopping terms", false);


  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "full-diag", 
					       "maximum Hilbert space dimension for which full diagonalization is applied", 
					       500, true, 100);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup) += new BooleanOption  ('\n', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", 
					 false);
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);  

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);

  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsDelta -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }




  //int MaxNbrIterLanczos = ((SingleIntegerOption*) Manager["iter-max"])->GetInteger();
  int NbrIterLanczos = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  int NbrEigenvalue = ((SingleIntegerOption*) Manager["nbr-eigen"])->GetInteger();
   
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");
  int NbrSites = Lx*Ly;  
  bool ReverseHopping = Manager.GetBoolean("negative-hopping");
  double ContactU = Manager.GetDouble("contactU");
  double Delta = Manager.GetDouble("deltaPotential");

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  
  int MaxFullDiagonalization = ((SingleIntegerOption*) Manager["full-diag"])->GetInteger();
  
  bool ResumeFlag = ((BooleanOption*) Manager["resume"])->GetBoolean();
  bool DiskFlag = ((BooleanOption*) Manager["disk"])->GetBoolean();
  int VectorMemory = ((SingleIntegerOption*) Manager["nbr-vector"])->GetInteger();
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();
  //char* SavePrecalculationFileName = ((SingleStringOption*) Manager["save-precalculation"])->GetString();

  int L = 0;
  double GroundStateEnergy = 0.0;
  //bool FirstRun = true;

  int NbrFluxValues = 1;
  if (NbrFluxQuanta == -1)
    {
      NbrFluxQuanta = 0;
      NbrFluxValues = NbrSites;
    }


  char* OutputName;
  if ( (OutputName = Manager.GetString("output-file")) == NULL)
    {
      OutputName = new char [256];
      char reverseHoppingString[4]="";
      char deltaString[20]="";
      if (ReverseHopping)
	sprintf(reverseHoppingString,"rh_");
      if (Delta!=0.0)
	sprintf(deltaString,"d_%g_",Delta);  
      if (NbrFluxValues == 1)
	sprintf (OutputName, "bosons_lattice_n_%d_x_%d_y_%d_u_%g_%s%sq_%d.dat", NbrBosons, Lx, Ly, ContactU, reverseHoppingString, deltaString, NbrFluxQuanta);
      else
	sprintf (OutputName, "bosons_lattice_n_%d_x_%d_y_%d_u_%g_%s%sq.dat", NbrBosons, Lx, Ly, ContactU, reverseHoppingString, deltaString);
    }
  BosonOnLattice TotalSpace= BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
  
  Architecture.GetArchitecture()->SetDimension(TotalSpace.GetHilbertSpaceDimension());
  
  AbstractQHEOnLatticeHamiltonian* Hamiltonian;
  Hamiltonian = new ParticleOnLatticeDeltaHamiltonian(&TotalSpace, NbrBosons, Lx, Ly, NbrFluxQuanta, ContactU,
						      ReverseHopping, Delta, Architecture.GetArchitecture(), Memory, LoadPrecalculationFileName);

  ofstream File;
  File.open(OutputName, ios::binary | ios::out);
  File.precision(14);


  {     
	cout << "----------------------------------------------------------------" << endl;
//	FermionOnTorus TotalSpace (NbrFermions, MaxMomentum, y);
	cout << " Total Hilbert space dimension = " << TotalSpace.GetHilbertSpaceDimension() << endl;
	//      cout << "momentum = " << Momentum << endl;
	cout << "flux = " << NbrFluxQuanta << endl;
//	cout << "momentum = (" << y << ")" << endl;
/*	for (int i = 0; i < TotalSpace.GetHilbertSpaceDimension(); ++i)
	  {
	    cout << i << " = ";
	    TotalSpace.PrintState(cout, i) << endl;
	  }
	cout << endl << endl;
	exit(0);*/
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
	Architecture.GetArchitecture()->SetDimension(TotalSpace.GetHilbertSpaceDimension());
	
	
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
		  File << NbrFluxQuanta << " " << TmpTriDiag.DiagonalElement(2 * j) << endl;
		  cout << NbrFluxQuanta << " " << TmpTriDiag.DiagonalElement(2 * j) << endl;
		}
	      cout << endl;
	    }
	  else
	    {
	      File << NbrFluxQuanta << " " << HRep(0, 0) << endl;
	    }
	}
      else
	{
	  int MaxNbrIterLanczos = 4000;
	  AbstractLanczosAlgorithm* Lanczos;
	  if (NbrEigenvalue == 1)
	    {
	      if (DiskFlag == false)
		Lanczos = new ComplexBasicLanczosAlgorithm(Architecture.GetArchitecture(), NbrEigenvalue, MaxNbrIterLanczos);
	      else
		Lanczos = new ComplexBasicLanczosAlgorithmWithDiskStorage(Architecture.GetArchitecture(), NbrEigenvalue, MaxNbrIterLanczos);
	    }
	  else
	    {
	      if (DiskFlag == false)
		Lanczos = new FullReorthogonalizedComplexLanczosAlgorithm (Architecture.GetArchitecture(), NbrEigenvalue, MaxNbrIterLanczos);
	      else
		Lanczos = new FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage (Architecture.GetArchitecture(), NbrEigenvalue, VectorMemory, 
											  MaxNbrIterLanczos);
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
	      File << NbrFluxQuanta << " " << TmpMatrix.DiagonalElement(i) << endl;
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
      cout << " ground state energy = " << GroundStateEnergy << endl;
      cout << " energy per particle in the ground state = " << (GroundStateEnergy / (double) NbrBosons) << endl;
      delete Hamiltonian;

    }
  File.close();

  return 0;
}
