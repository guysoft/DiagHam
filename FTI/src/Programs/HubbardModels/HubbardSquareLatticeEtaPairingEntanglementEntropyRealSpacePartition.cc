#include "Options/Options.h"

#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelSimpleSquareLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"

#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/FactorialCoefficient.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;



// extract the correlation matrix for a given region out of the correlation matrix for a bigger region
//
// correlationMatrix = correlation matrix for the bigger region
// sourceNbrSitesX = number of sites along x for the bigger region
// sourceNbrSitesY = number of sites along y for the bigger region
// targetNbrSitesX = number of sites along x for the smaller region
// targetNbrSitesY = number of sites along y for the smaller region
HermitianMatrix EtaPairaingEntanglementEntropyExtractCorrelationMatrix(HermitianMatrix& correlationMatrix, int sourceNbrSitesX, int sourceNbrSitesY, 
								       int targetNbrSitesX, int targetNbrSitesY);


int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("HubbardSquareLatticeEtaPairingEntanglementEntropyRealSpacePartition" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-pairs", "number of pairs", 0);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of unit cells along the x direction", 4);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of unit cells along the y direction", 4);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrsitex-a", "number of unit cells along the x direction for the part A", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrsitey-a", "number of unit cells along the y direction for the part A", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-nbrsitexa", "maximum number of unit cells along the x direction for the part A (equal to --nbrsitex-a if negative)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-nbrsiteya", "maximum number of unit cells along the y direction for the part A (equal to --nbrsitey-a if negative)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nearbyeta-x", "x distance of the broken pair when generating a nearby eta pairing state", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nearbyeta-y", "y distance of the broken pair when generating a nearby eta pairing state", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-nonvacuum", "apply the eta^+ operators to a non-vacuum state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-particles", "number of particles for the non-vacuum state", 0);
  (*SystemGroup) += new BooleanOption ('\n', "use-fermisea", "apply the eta^+ operators to the Fermi sea");
  (*SystemGroup) += new SingleStringOption  ('\n', "nonvacuum-file", "provide the description of the non-vacuum state as a two-column ASCII file");
  (*SystemGroup) += new BooleanOption ('\n', "use-random", "generate a random non-vacuum state");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "run-id", "add an additional run id to the file name when using the --use-random option", 0);
  
  (*SystemGroup) += new BooleanOption  ('\n', "show-nonvacuum", "show the non-vacuum state in the momentum basis");
  (*SystemGroup) += new BooleanOption  ('\n', "show-time", "show time required for each operation");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardSquareLatticeEtaPairingEntanglementEntropyRealSpacePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  bool ShowTimeFlag = Manager.GetBoolean("show-time");
  int NbrPairs = Manager.GetInteger("nbr-pairs"); 

  int XMomentum = 0;
  int YMomentum = 0;
  int NbrSitesX = 0;
  int NbrSitesY = 0;
  int NbrSites = 0; 
  int TotalSz = 0;
  bool Statistics = true;


  NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  NbrSitesY = Manager.GetInteger("nbr-sitey"); 
  NbrSites = NbrSitesX * NbrSitesY; 

  if ((NbrSitesX & 1) != 0)
    {
      cout << "error, eta pairing states require an even number of sites in the x direction" << endl;
      return 0;
    }
  if ((NbrSitesY & 1) != 0)
    {
      cout << "error, eta pairing states require an even number of sites in the y direction" << endl;
      return 0;
    }

  if ((NbrPairs & 1) != 0)
    {
      XMomentum = NbrSitesX >> 1;
      YMomentum = NbrSitesY >> 1;          
    }
  int NbrSitesXA = Manager.GetInteger("nbrsitex-a"); 
  int NbrSitesYA = Manager.GetInteger("nbrsitey-a"); 
  int NbrSitesA = NbrSitesXA * NbrSitesYA;
  int MaxNbrSitesXA = Manager.GetInteger("max-nbrsitexa"); 
  if (MaxNbrSitesXA < NbrSitesXA)
    {
      MaxNbrSitesXA = NbrSitesXA;
    }
  int MaxNbrSitesYA = Manager.GetInteger("max-nbrsiteya"); 
  if (MaxNbrSitesYA < NbrSitesYA)
    {
      MaxNbrSitesYA = NbrSitesYA;
    }
  int MaxNbrSitesA = MaxNbrSitesXA * MaxNbrSitesYA;

  
  Abstract2DTightBindingModel* TightBindingModel;
  TightBindingModel = new TightBindingModelSimpleSquareLattice (NbrSitesX, NbrSitesY, 1.0, 0.0, 0.0, 0.0,
								Architecture.GetArchitecture(), true);

  double* TightBindingModelEnergies = 0;
  int* TightBindingModelLinearizedMomenta = 0;
  TightBindingModel->GetEnergies(TightBindingModelEnergies, TightBindingModelLinearizedMomenta, 0);

  int* VacuumOneBodyLinearizedMomenta = 0; 
  double VacuumTotalEnergy = 0.0;
  int VacuumXMomentum = 0;
  int VacuumYMomentum = 0;
  int TmpMomentumX;
  int TmpMomentumY;
  int VacuumNbrParticles = Manager.GetInteger("nbr-particles"); 
  int VacuumTotalSz = VacuumNbrParticles;

  if (Manager.GetBoolean("use-fermisea") == true)
    { 
      VacuumOneBodyLinearizedMomenta = new int[VacuumNbrParticles];
      for (int i = 0; i < VacuumNbrParticles; ++i)
	{
	  TightBindingModel->GetLinearizedMomentumIndex(TightBindingModelLinearizedMomenta[i], TmpMomentumX, TmpMomentumY);
	  VacuumXMomentum += TmpMomentumX;
	  VacuumYMomentum += TmpMomentumY;
	  VacuumOneBodyLinearizedMomenta[i] = TightBindingModel->GetLinearizedMomentumIndex(TmpMomentumX, TmpMomentumY);
	  VacuumTotalEnergy += TightBindingModel->GetEnergy(0, VacuumOneBodyLinearizedMomenta[i]);
	}
    }
  else
    {
     if (Manager.GetString("nonvacuum-file") != 0)
       {
	 MultiColumnASCIIFile NonVacuumFile;
	 if (NonVacuumFile.Parse(Manager.GetString("nonvacuum-file")) == false)
	   {
	     NonVacuumFile.DumpErrors(cout);
	     return -1;
	   }
	 VacuumNbrParticles = NonVacuumFile.GetNbrLines();
	 VacuumTotalSz = VacuumNbrParticles;
	 VacuumOneBodyLinearizedMomenta = new int[VacuumNbrParticles];
	 int* TmpXMomenta = NonVacuumFile.GetAsIntegerArray(0);
	 int* TmpYMomenta = NonVacuumFile.GetAsIntegerArray(1);
	 for (int i = 0; i < VacuumNbrParticles; ++i)
	   {
	     VacuumOneBodyLinearizedMomenta[i] = TightBindingModel->GetLinearizedMomentumIndex(TmpXMomenta[i], TmpYMomenta[i]);
	     VacuumTotalEnergy += TightBindingModel->GetEnergy(0, VacuumOneBodyLinearizedMomenta[i]);
	     VacuumXMomentum += TmpXMomenta[i];
	     VacuumYMomentum += TmpYMomenta[i];	     
	   }
       }
     else
       {
	 if (Manager.GetBoolean("use-random") == true)
	   {
	     VacuumOneBodyLinearizedMomenta = new int[VacuumNbrParticles];
	     AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (0);
	     RandomNumber->UseTimeSeed();
	     for (int i = 0; i < VacuumNbrParticles; ++i)
	       {
		 VacuumOneBodyLinearizedMomenta[i] = TightBindingModelLinearizedMomenta[i];
	       }
	     for (int i = VacuumNbrParticles; i < NbrSites; ++i)
	       {
		 int Tmp = (int) (RandomNumber->GetRealRandomNumber() * (((double) i) + 0.0001));
		 if (Tmp < VacuumNbrParticles)
		   VacuumOneBodyLinearizedMomenta[Tmp] = TightBindingModelLinearizedMomenta[i];
	       }
	     for (int i = 0; i < VacuumNbrParticles; ++i)
	       {
		 VacuumTotalEnergy += TightBindingModel->GetEnergy(0, VacuumOneBodyLinearizedMomenta[i]);
		 TightBindingModel->GetLinearizedMomentumIndex(VacuumOneBodyLinearizedMomenta[i], TmpMomentumX, TmpMomentumY);
		 VacuumXMomentum += TmpMomentumX;
		 VacuumYMomentum += TmpMomentumY;	     
	       }
	   }
       }
    }

  if (Manager.GetBoolean("show-nonvacuum") == true)
    {
      cout << "using states : " << endl;
      for (int i = 0; i < VacuumNbrParticles; ++i)
	{
	  TightBindingModel->GetLinearizedMomentumIndex(VacuumOneBodyLinearizedMomenta[i], TmpMomentumX, TmpMomentumY);
	  cout << "(" << TmpMomentumX << ", " << TmpMomentumY << ")" << endl;
	}
    }
  
  XMomentum += VacuumXMomentum;
  YMomentum += VacuumYMomentum;      
  VacuumXMomentum %= NbrSitesX;
  VacuumYMomentum %= NbrSitesY;
  XMomentum %= NbrSitesX;
  YMomentum %= NbrSitesY;

  int NbrParticles = VacuumNbrParticles; 
  NbrParticles += 2 * NbrPairs;
  TotalSz += VacuumTotalSz;
  char* StatisticPrefix = new char [64];
  sprintf (StatisticPrefix, "fermions_hubbard");
  char* FilePrefix = new char [256];
  if ((Manager.GetInteger("nearbyeta-x") == 0) && (Manager.GetInteger("nearbyeta-y") == 0))
    {
      sprintf (FilePrefix, "%s_square_etapairing_nbrpairs_%ld_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, Manager.GetInteger("nbr-pairs"), NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
    }
  else
    {
      sprintf (FilePrefix, "%s_square_nearbyetapairing_nbrpairs_%ld_alphax_%ld_alphay_%ld_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, Manager.GetInteger("nbr-pairs"), Manager.GetInteger("nearbyeta-x"), 
	       Manager.GetInteger("nearbyeta-y"), NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
    }
  char* EntropyFileName = new char [512];
  if (Manager.GetBoolean("use-random") == true)
    {
      sprintf(EntropyFileName, "%s_sz_%d_runid_%ld_dxa_%d_ya_%d.ent", FilePrefix, TotalSz, Manager.GetInteger("run-id"), 
	      MaxNbrSitesXA, MaxNbrSitesYA);
    }
  else
    {
      sprintf(EntropyFileName, "%s_sz_%d_xa_%d_ya_%d.ent", FilePrefix, TotalSz, MaxNbrSitesXA, MaxNbrSitesYA);
    }
  if (Manager.GetString("nonvacuum-file") == 0)
    {
      char* NonVacuumFileName = new char [512];
      if (Manager.GetBoolean("use-random") == true)
	{
	  sprintf(NonVacuumFileName, "%s_sz_%d_runid_%ld_nonvacuum.dat", FilePrefix, TotalSz, Manager.GetInteger("run-id"));
	}
      else
	{
	  sprintf(NonVacuumFileName, "%s_sz_%d_nonvacuum.dat", FilePrefix, TotalSz);
	}
      ofstream NonVacuumFile;
      NonVacuumFile.open(NonVacuumFileName, ios::binary | ios::out);
      NonVacuumFile.precision(14);
      NonVacuumFile << "# Non-vacuum state total momentum along x = " << VacuumXMomentum << endl;
      NonVacuumFile << "# Non-vacuum state total momentum along y = " << VacuumYMomentum << endl;
      NonVacuumFile << "# Non-vacuum state total energy = " << VacuumTotalEnergy << endl;
      for (int i = 0; i < VacuumNbrParticles; ++i)
	{
	  TightBindingModel->GetLinearizedMomentumIndex(VacuumOneBodyLinearizedMomenta[i], TmpMomentumX, TmpMomentumY);
	  NonVacuumFile << TmpMomentumX << " " << TmpMomentumY << endl;
	}
      NonVacuumFile.close();
    }

  ofstream File;
  File.open(EntropyFileName, ios::binary | ios::out);
  File.precision(14);

  File << "# Non-vacuum state total momentum along x = " << VacuumXMomentum << endl;
  File << "# Non-vacuum state total momentum along y = " << VacuumYMomentum << endl;
  File << "# Non-vacuum state total energy = " << VacuumTotalEnergy << endl;
  cout << "Non-vacuum state total momentum along x = " << VacuumXMomentum << endl;
  cout << "Non-vacuum state total momentum along y = " << VacuumYMomentum << endl;
  cout << "Non-vacuum state total energy = " << VacuumTotalEnergy << endl;

  timeval TotalStartingTime;
  timeval TotalEndingTime;
  if (ShowTimeFlag == true)
    {
      gettimeofday (&(TotalStartingTime), 0);
    }
  HermitianMatrix EntanglementHamiltonian = TightBindingModel->EvaluateFullTwoPointCorrelationFunction(MaxNbrSitesXA, MaxNbrSitesYA, VacuumOneBodyLinearizedMomenta, VacuumNbrParticles, 0);
  if (ShowTimeFlag == true)
    {
      gettimeofday (&(TotalEndingTime), 0);
      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
      cout << "correlation matrix evaluated in " << Dt << "s" << endl;
    }
  for (; NbrSitesXA <= MaxNbrSitesXA; ++NbrSitesXA)
    {
      NbrSitesYA = NbrSitesXA;
      cout << "computing entropy of a " << NbrSitesXA << "x" << NbrSitesYA << " patch" << endl;
      int TotalNbrSitesA = NbrSitesXA * NbrSitesYA;
      if (ShowTimeFlag == true)
	{
	  gettimeofday (&(TotalStartingTime), 0);
	}
      RealDiagonalMatrix VacuumOneBodyEntanglementEnergies(TotalNbrSitesA, true);
      HermitianMatrix TmpEntanglementHamiltonian = EtaPairaingEntanglementEntropyExtractCorrelationMatrix(EntanglementHamiltonian, MaxNbrSitesXA, MaxNbrSitesYA, NbrSitesXA, NbrSitesYA);
#ifdef __LAPACK__
      TmpEntanglementHamiltonian.LapackDiagonalize(VacuumOneBodyEntanglementEnergies);
#else
      TmpEntanglementHamiltonian.Diagonalize(VacuumOneBodyEntanglementEnergies);
#endif
      VacuumOneBodyEntanglementEnergies.SortMatrixUpOrder();
      int MinOneBodyEntanglementEnergyIndex = 0;
      int MaxOneBodyEntanglementEnergyIndex = VacuumOneBodyEntanglementEnergies.GetNbrRow() - 1;  
      while ((MinOneBodyEntanglementEnergyIndex <= MaxOneBodyEntanglementEnergyIndex) && (VacuumOneBodyEntanglementEnergies[MinOneBodyEntanglementEnergyIndex] <= 0.0))
	++MinOneBodyEntanglementEnergyIndex;
      while ((MinOneBodyEntanglementEnergyIndex <= MaxOneBodyEntanglementEnergyIndex) && (VacuumOneBodyEntanglementEnergies[MaxOneBodyEntanglementEnergyIndex] >= 1.0))
	--MaxOneBodyEntanglementEnergyIndex;
      if (ShowTimeFlag == true)
	{
	  gettimeofday (&(TotalEndingTime), 0);
	  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	  cout << "diagonalization done in " << Dt << "s" << endl;
	}
      
      double EntanglementEntropy = 0.0;
      int NbrVacuumOneBodyEntanglementTrimmedEnergies = MaxOneBodyEntanglementEnergyIndex - MinOneBodyEntanglementEnergyIndex + 1;
      int NbrRejectedOneBodyEntropies = VacuumOneBodyEntanglementEnergies.GetNbrRow() - NbrVacuumOneBodyEntanglementTrimmedEnergies;
      double* VacuumOneBodyEntanglementTrimmedEnergies = new double[NbrVacuumOneBodyEntanglementTrimmedEnergies];
      for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies; ++i)
	{
	  VacuumOneBodyEntanglementTrimmedEnergies[i] = VacuumOneBodyEntanglementEnergies[i + MinOneBodyEntanglementEnergyIndex];
	}
      if (NbrPairs == 0)
	{
	  for (int i = 0; i < NbrVacuumOneBodyEntanglementTrimmedEnergies; ++i)
	    {
	      EntanglementEntropy -= VacuumOneBodyEntanglementTrimmedEnergies[i] * log (VacuumOneBodyEntanglementTrimmedEnergies[i]);
	      EntanglementEntropy -= (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i]) * log (1.0 - VacuumOneBodyEntanglementTrimmedEnergies[i]);
	    }
	}
      else
	{
	  int MaxNbrParticlesA = NbrSitesA;
	  if (MaxNbrParticlesA > NbrParticles)
	    MaxNbrParticlesA = NbrParticles;
	  FactorialCoefficient TmpCoefficient;
	  
	  for (int i = 0; i <= MaxNbrParticlesA; ++i)
	    {
	      //GetEntanglementEntropyPerNbrParticlesA(VacuumOneBodyEntanglementTrimmedEnergies, NbrVacuumOneBodyEntanglementTrimmedEnergies);
	      int MaxSumIndex = NbrSitesA - i;
	      if (MaxSumIndex > NbrPairs)
		MaxSumIndex = NbrPairs;
	      double Tmp = 0.0;
	      for (int j = 0; j < MaxSumIndex; ++j)
		{
		  TmpCoefficient.SetToOne();
		  TmpCoefficient.BinomialMultiply(NbrSitesA - i, j);
		  TmpCoefficient.BinomialMultiply(NbrSites - NbrSitesA - VacuumNbrParticles + i, NbrPairs - j);
		  TmpCoefficient.BinomialDivide(NbrSites - VacuumNbrParticles, NbrPairs);
		  double Tmp2 = TmpCoefficient.GetNumericalValue();
		  Tmp -= Tmp2 * log (Tmp2);
		}
	      cout << i << " " << Tmp << endl;
	    }
	  
	}
      cout << "Entanglement entropy = " << EntanglementEntropy << endl;
      cout << "Nbr Rejected one-body entanglement energies = " << NbrRejectedOneBodyEntropies << " / " << TotalNbrSitesA << endl;
      File << NbrSitesXA << " " << NbrSitesYA << " " << EntanglementEntropy << endl;
      delete[] VacuumOneBodyEntanglementTrimmedEnergies;
    }
  File.close();
  return 0;
}


// extract the correlation matrix for a given region out of the correlation matrix for a bigger region
//
// correlationMatrix = correlation matrix for the bigger region
// sourceNbrSitesX = number of sites along x for the bigger region
// sourceNbrSitesY = number of sites along y for the bigger region
// targetNbrSitesX = number of sites along x for the smaller region
// targetNbrSitesY = number of sites along y for the smaller region

HermitianMatrix EtaPairaingEntanglementEntropyExtractCorrelationMatrix(HermitianMatrix& correlationMatrix, int sourceNbrSitesX, int sourceNbrSitesY, 
								       int targetNbrSitesX, int targetNbrSitesY)
{
  int TargetTotalNbrSites = targetNbrSitesX * targetNbrSitesY;
  HermitianMatrix TmpMatrix (TargetTotalNbrSites, true);
  Complex Tmp;
  for (int i = 0; i < TargetTotalNbrSites; ++i)
    {
      int TmpTargetX1 = i / targetNbrSitesY;
      int TmpTargetY1 = i % targetNbrSitesY;
      int TmpSourceIndex1 = (TmpTargetX1 * sourceNbrSitesY) + TmpTargetY1;
      for (int j = i; j < TargetTotalNbrSites; ++j)
	{
	  int TmpTargetX2 = j / targetNbrSitesY;
	  int TmpTargetY2 = j % targetNbrSitesY;
	  int TmpSourceIndex2 = (TmpTargetX2 * sourceNbrSitesY) + TmpTargetY2;
	  correlationMatrix.GetMatrixElement(TmpSourceIndex1, TmpSourceIndex2, Tmp);
	  TmpMatrix.SetMatrixElement(i, j, Tmp);
	}
    }
  return TmpMatrix;
}

// double GetEntanglementEntropyPerNbrParticlesA(double* oneBodyEntanglementTrimmedEnergies, int nbrOneBodyEntanglementTrimmedEnergies, int nbrParticlesA)
// {
  
// }
