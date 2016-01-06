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

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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


// apply the time evolution operator to a state
//
// tau = amount of time to evolve
// hamiltonian = pointer to the Hamiltonian
// inputVector = reference on the input state
// outputVector = reference on the output state (allocation should be done outside this method)
// architecture = pointer to the architecture
// convergenceError = error below which the nom of the time evolved state should be considered equal to one
// return value = error on the time evolved state norm (the outputVector state is normalized to one)
double ApplyTimeEvolution(double tau, AbstractHamiltonian* hamiltonian, ComplexVector& inputVector, ComplexVector& outputVector, 
			  AbstractArchitecture* architecture, double convergenceError = MACHINE_PRECISION);


int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("HubbardSquareLatticeEtaPairingStateTimeEvolution" , "0.01");
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

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of unit cells along the x direction", 2);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of unit cells along the y direction", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use the Sz <-> -Sz symmetry");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive on-site (Hubbard) potential strength", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "nn-t", "nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "nnn-t", "next nearest neighbor hopping amplitude", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1); 
  (*SystemGroup) += new BooleanOption  ('\n', "no-evolution", "do not perform any time evolution and just store the eta pairing state");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "friendlyshow-hamiltonian", "show matrix representation of the hamiltonian, displaying only non-zero matrix elements");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleDoubleOption ('\n', "testhermitian-error", "error threshold when testing hermiticy (0 for machine accuracy)", 0.0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HubbardSquareLatticeEtaPairingStateTimeEvolution -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey"); 
  if ((NbrParticles & 1) != 0)
    {
      cout << "error, eta pairing states require an even number of particles" << endl;
      return 0;
    }
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
  int NbrSites = NbrSitesX * NbrSitesY; 
  bool SzSymmetryFlag = Manager.GetBoolean("szsymmetrized-basis");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  int XMomentum = 0;
  int YMomentum = 0;
  if ((NbrParticles & 2) != 0)
    {
      XMomentum = NbrSitesX >> 1;
      YMomentum = NbrSitesY >> 1;      
    }
  int TotalSz = 0;
  int SzParitySector = 1;

  char* StatisticPrefix = new char [64];
  if (SzSymmetryFlag == false)
    {
      sprintf (StatisticPrefix, "fermions_hubbard");
    }
  else
    {
      sprintf (StatisticPrefix, "fermions_hubbard_szsym");
    }
    
  

  char* FilePrefix = new char [256];
  sprintf (FilePrefix, "%s_square_etapairing_x_%d_y_%d_n_%d_ns_%d", StatisticPrefix, NbrSitesX, NbrSitesY, NbrParticles, NbrSites);
  
  char* FileParameterString = new char [256];
  sprintf (FileParameterString, "t_%.6f_tp_%.6f", Manager.GetDouble("nn-t"), Manager.GetDouble("nnn-t"));

  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetDouble("u-potential") == 0.0)
    sprintf(EigenvalueOutputFile, "%s_%s_sz_%d.dat", FilePrefix, FileParameterString, TotalSz);
  else
    sprintf(EigenvalueOutputFile, "%s_%s_u_%f_sz_%d.dat", FilePrefix, FileParameterString, Manager.GetDouble("u-potential"), TotalSz);

  Abstract2DTightBindingModel* TightBindingModel;
  TightBindingModel = new TightBindingModelSimpleSquareLattice (NbrSitesX, NbrSitesY, Manager.GetDouble("nn-t"), Manager.GetDouble("nnn-t"), 0.0, 0.0,
								Architecture.GetArchitecture(), true);

  RealSymmetricMatrix DensityDensityInteractionupup(NbrSites, true);
  RealSymmetricMatrix DensityDensityInteractiondowndown(NbrSites, true);
  RealSymmetricMatrix DensityDensityInteractionupdown(NbrSites, true);
  if (Manager.GetDouble("u-potential") != 0.0)
    {
      double UPotential = Manager.GetDouble("u-potential");
      for (int i = 0; i < NbrSites; ++i)
	{
	  DensityDensityInteractionupdown.SetMatrixElement(i, i, UPotential);
	}
    }

  FermionOnLatticeWithSpinRealSpaceAnd2DTranslation* Space = 0;
  AbstractHamiltonian* Hamiltonian = 0;
  if (SzSymmetryFlag == false)
    {
      Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, 0, NbrSites, XMomentum, NbrSitesX,
								     YMomentum, NbrSitesY);
    }
  else
    {
      bool MinusParitySector = true;
      Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (NbrParticles, 0, NbrSites, MinusParitySector, 
									       XMomentum, NbrSitesX,
									       YMomentum, NbrSitesY);
    }
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

  ComplexVector EtaPairingState = Space->GenerateEtaPairingState();

  HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
  Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian(Space, NbrParticles, NbrSites,XMomentum, NbrSitesX,
										  YMomentum, NbrSitesY, TightBindingMatrix, TightBindingMatrix,
										  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
										  DensityDensityInteractionupdown, 
										  Architecture.GetArchitecture(), Memory);
  if (Manager.GetBoolean("no-evolution"))  
    {
      char* EigenstateOutputFile = new char [512];
      if (SzSymmetryFlag == false)
	{
	  sprintf (EigenstateOutputFile, "%s_kx_%d_ky_%d_sz_%d.0.vec", FilePrefix, XMomentum, YMomentum, TotalSz);
	}
      else
	{
	  sprintf (EigenstateOutputFile, "%s_szp_%d_kx_%d_ky_%d_sz_%d.0.vec", FilePrefix, SzParitySector, XMomentum, YMomentum, TotalSz);
	}
      ComplexVector TmpVector (EtaPairingState.GetVectorDimension(), true);
      VectorHamiltonianMultiplyOperation Operation1 (Hamiltonian, &EtaPairingState, &TmpVector);
      Operation1.ApplyOperation(Architecture.GetArchitecture());
      cout << "check eta pairing state E=" << (TmpVector * EtaPairingState) 
	   << " var(E)" << ((TmpVector * TmpVector) -  ((TmpVector * EtaPairingState) * (TmpVector * EtaPairingState))) << endl;
      if (EtaPairingState.WriteVector(EigenstateOutputFile) == false)
	{
	  cout << "error while writing " << EigenstateOutputFile << endl;
	}
      delete[] EigenstateOutputFile;
    }
  else
    {
      ComplexVector TmpVector (EtaPairingState.GetVectorDimension(), true);
      double Tau = 1.0;
      double EvolutionError = ApplyTimeEvolution(Tau, Hamiltonian, EtaPairingState, TmpVector, Architecture.GetArchitecture());
    }

  delete Hamiltonian;
  delete Space;
  return 0;
}

// apply the time evolution operator to a state
//
// tau = amount of time to evolve
// hamiltonian = pointer to the Hamiltonian
// inputVector = reference on the input state
// outputVector = reference on the output state (allocation should be done outside this method)
// architecture = pointer to the architecture
// convergenceError = error below which the nom of the time evolved state should be considered equal to one
// return value = error on the time evolved state norm (the outputVector state is normalized to one)

double ApplyTimeEvolution(double tau, AbstractHamiltonian* hamiltonian, ComplexVector& inputVector, ComplexVector& outputVector, 
			  AbstractArchitecture* architecture, double convergenceError)
{  
  Complex TmpFactor (0.0, -tau);
  outputVector.Copy(inputVector);
  ComplexVector TmpVector (inputVector.GetLargeVectorDimension());
  ComplexVector TmpVector2 (inputVector.GetLargeVectorDimension());
  VectorHamiltonianMultiplyOperation Operation1 (hamiltonian, &inputVector, &TmpVector);
  Operation1.ApplyOperation(architecture);
  outputVector.AddLinearCombination(TmpFactor, TmpVector);
  double TmpNorm = outputVector.Norm();
  int Index = 2;
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  timeval TmpStartingTime;
  timeval TmpEndingTime;
  gettimeofday (&(TotalStartingTime), 0);
  while (fabs(1.0 - TmpNorm) > convergenceError)
    {
      TmpFactor *= Complex(0.0 , -tau / ((double) Index));
      timeval TmpStartingTime;
      gettimeofday (&(TmpStartingTime), 0);
      VectorHamiltonianMultiplyOperation Operation2 (hamiltonian, &TmpVector, &TmpVector2);
      Operation2.ApplyOperation(architecture);
      outputVector.AddLinearCombination(TmpFactor, TmpVector2);
      ComplexVector TmpVector3 = TmpVector;
      TmpVector = TmpVector2;
      TmpVector2 = TmpVector3;
      TmpNorm = outputVector.Norm();
      gettimeofday (&(TmpEndingTime), 0);
      double Dt = (double) ((TmpEndingTime.tv_sec - TmpStartingTime.tv_sec) + 
			    ((TmpEndingTime.tv_usec - TmpStartingTime.tv_usec) / 1000000.0));	
      cout << Index << " " << TmpNorm << " " << fabs(1.0 - TmpNorm) << " (" <<  Dt << "s)" << endl;
      ++Index;
    }
   gettimeofday (&(TotalEndingTime), 0);
   double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
			 ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));	
   cout << "time evolution of tau=" << tau << " done in " << Dt << "s" << endl;
   outputVector /= TmpNorm;
  return fabs(1.0 - TmpNorm);
}
