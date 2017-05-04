#include "Options/Options.h"

#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSU4SpinMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU4SpinMomentumSpace.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong.h"
#include "HilbertSpace/BosonOnLatticeRealSpace.h"
#include "HilbertSpace/BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslationLong.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslationLong.h"


#include "Hamiltonian/ParticleOnLatticeHofstadterSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeHofstadterSingleBandGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeTwoBandHofstadterHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeFourBandHofstadterHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian.h"

#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"
#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Architecture/ArchitectureOperation/AddComplexLinearCombinationOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FCIHofstadterModelTimeEvolution" , "0.01");
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

  (*SystemGroup) += new SingleStringOption('\n', "initial-state", "name of the file containing the initial vector upon which e^{-iHt} acts");  
  (*SystemGroup) += new BooleanOption  ('\n', "complex", "initial vector is a complex vector");
  (*SystemGroup) += new BooleanOption  ('\n', "compute-energy", "compute the energy of each time-evolved vector");
  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-final", "final value of the onsite potential strength", -1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tau", "time where the final value of U is reached", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "final-time", "time where the code stops", 10.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-steps", "number of points to evaluate", 1);

  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  
  (*SystemGroup) += new BooleanOption  ('\n', "landau-x", "Use Landau gauge along the x-axis within unit cell");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrParticles; 
  int NbrCellX; 
  int NbrCellY;
  int UnitCellX; 
  int UnitCellY;
  int FluxPerCell;
  int MinNbrSinglets;
  int Sz;
  int Kx,Ky;
  int SzSymmetry;
  bool StatisticFlag;
  bool GutzwillerFlag;
  bool TranslationFlag;
  bool UsingConstraintMinNbrSinglets;
  char Axis ='y';
  bool SzSymmetryFlag;
  if (Manager.GetBoolean("landau-x"))
    Axis ='x';
  
  
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  if (Manager.GetString("initial-state") == 0)
    {
      cout << " Error: an initial state must be provided" << endl;
      return -1;
    }
  char* StateFileName = Manager.GetString("initial-state");
  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }
  
  if (FTIHofstadterdModelWith2DTranslationFindSystemInfoFromVectorFileName(StateFileName, NbrParticles, Kx, Ky, FluxPerCell , NbrCellX,  NbrCellY, UnitCellX, UnitCellY,  StatisticFlag, GutzwillerFlag, TranslationFlag) == false )
    {
      return -1;
    }

  if (FTIHofstadterModelWithSzFindSystemInfoFromVectorFileName(StateFileName, Sz, SzSymmetry, MinNbrSinglets, SzSymmetryFlag,UsingConstraintMinNbrSinglets) == false)
    {
      return -1;
    }

  Abstract2DTightBindingModel *TightBindingModel;
  
  TightBindingModel = new TightBindingModelHofstadterSquare(NbrCellX, NbrCellY, UnitCellX, UnitCellY, FluxPerCell, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), true, false);
  
  
  ParticleOnSphere* Space = 0;
  AbstractQHEHamiltonian* Hamiltonian = 0;
  
  int NbrSites = TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand();
  
  if (TranslationFlag)
    {
      if ((SzSymmetryFlag) && (Sz == 0))
	{				      
#ifdef __64_BITS__
	  if ( NbrSites < 31)
#else
	    if (NbrSites < 15)
#endif
	      {
		if (UsingConstraintMinNbrSinglets == false)
		  {
		    Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslation (NbrParticles, Sz,NbrSites ,(SzSymmetry == -1), Kx, NbrCellX, Ky,  NbrCellY );	  
		  }
		else
		  {
		    Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, MinNbrSinglets, Sz, NbrSites, (SzSymmetry == -1),   Kx, NbrCellX, Ky,  NbrCellY, 10000000);
		  } 
	      }
	    else
	      {
		if (UsingConstraintMinNbrSinglets == false)
		  {
		    Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationLong (NbrParticles, Sz, NbrSites, (SzSymmetry == -1),Kx, NbrCellX, Ky,  NbrCellY );	  
		  }
		else
		  {
		    Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, MinNbrSinglets, Sz, NbrSites , (SzSymmetry == -1),   Kx, NbrCellX, Ky,  NbrCellY);
		  }
	      }
	}
      else
	{
#ifdef __64_BITS__
	  if (NbrSites < 31)
#else
	    if ( NbrSites < 15)
#endif
	      {
		if (UsingConstraintMinNbrSinglets == false)
		  {
		    Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation (NbrParticles, Sz, NbrSites, Kx, NbrCellX, Ky,  NbrCellY);	  
		  }
		else
		  {
		    Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets (NbrParticles, MinNbrSinglets, Sz,  NbrSites , Kx, NbrCellX, Ky,  NbrCellY, 10000000ul);
		  }
	      }
	    else
	      {
		if (UsingConstraintMinNbrSinglets == false)
		  {
		    Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationLong (NbrParticles, Sz,  NbrSites, Kx, NbrCellX, Ky,  NbrCellY);	  
		  }
		else
		  {
		    Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSingletsLong (NbrParticles, MinNbrSinglets, Sz,  NbrSites, Kx, NbrCellX, Ky,  NbrCellY, 10000000ul);
		  }
	      }
	}
    }
  else
    {
      if ((SzSymmetryFlag) && (Sz == 0))
	{
	  Space = new FermionOnLatticeWithSpinSzSymmetryRealSpace (NbrParticles, Sz,  NbrSites, (SzSymmetry == -1));	  
	}
      else
	{
	  Space = new FermionOnLatticeWithSpinRealSpace (NbrParticles, Sz, NbrSites,10000000);	  
	}
    }

  

  ComplexVector TmpInitialState (Space->GetHilbertSpaceDimension());
  if (Manager.GetBoolean("complex") == false)
    {
      RealVector InputState;
      if (InputState.ReadVector(StateFileName) == false)
	{
	  cout << "error while reading " << StateFileName << endl;
	  return -1;
	}
      if (InputState.GetVectorDimension() != Space->GetHilbertSpaceDimension())
	{
	  cout << "error: vector and Hilbert-space have unequal dimensions " << InputState.GetVectorDimension() << " "<< Space->GetHilbertSpaceDimension() << endl;
	  return -1;
	}
      TmpInitialState = InputState;
    }
  else
    {
      ComplexVector InputState;
      if (InputState.ReadVector(StateFileName) == false)
	{
	  cout << "error while reading " << StateFileName << endl;
	  return -1;
	}
      if (InputState.GetVectorDimension() != Space->GetHilbertSpaceDimension())
	{
	  cout << "error: vector and Hilbert-space have unequal dimensions " << InputState.GetVectorDimension() << " "<< Space->GetHilbertSpaceDimension() << endl;
	  return -1;
	}
      TmpInitialState = InputState;
    }
  

  RealSymmetricMatrix DensityDensityInteractionupup(NbrSites, true);
  RealSymmetricMatrix DensityDensityInteractiondowndown(NbrSites, true);
  RealSymmetricMatrix DensityDensityInteractionupdown(NbrSites, true);
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  HermitianMatrix TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();
      

  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-final", "final value of the onsite potential strength", -1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tau", "time where the final value of U is reached", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "final-time", "time where the code stops", 10.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-steps", "number of points to evaluate", 1);
  
  int NbrSteps = Manager.GetInteger("nbr-steps");
  double UFinal = Manager.GetDouble("u-final");
  double Tau = Manager.GetDouble("tau");
  double FinalTime = Manager.GetDouble("final-time");
  double TimeStep = FinalTime/ NbrSteps;
  double UPotential;
  
  char* OutputNamePrefix = new char [512];
  strcpy(OutputNamePrefix,StateFileName);
  
  RemoveExtensionFromFileName(OutputNamePrefix,".vec");
  char * ParameterString = new char [512];
  sprintf(ParameterString,"uf_%f_tau_%f_tf_%f_nstep_%d",UFinal,Tau,FinalTime,NbrSteps);
  
  for(int i = 0 ; i < NbrSteps; i++)
    {
      double UPotential =  UFinal/ Tau * TimeStep * i;
      if ( fabs(UPotential) > fabs(UFinal) ) 
	UPotential = UFinal ;
      for (int p = 0; p < NbrSites; ++p)
	{
	  DensityDensityInteractionupdown.SetMatrixElement(p, p, UPotential);
	}
      
      if (TranslationFlag)
	{
	  Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceAnd2DTranslationHamiltonian((ParticleOnSphereWithSpin  *)  Space, NbrParticles, NbrSites, Kx, NbrCellX, Ky,  NbrCellY,
											  TightBindingMatrix, TightBindingMatrix,
											  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
											  DensityDensityInteractionupdown, 
											  Architecture.GetArchitecture(), Memory);
	}
      else
	{
	  Hamiltonian = new ParticleOnLatticeWithSpinRealSpaceHamiltonian((ParticleOnSphereWithSpin  *)  Space, NbrParticles, NbrSites,
									  TightBindingMatrix, TightBindingMatrix,
									  DensityDensityInteractionupup, DensityDensityInteractiondowndown, 
									  DensityDensityInteractionupdown, 
									  Architecture.GetArchitecture(), Memory);
	}

      double Norm;
      int TmpExpansionOrder;
      ComplexVector TmpState (Space->GetHilbertSpaceDimension()) ;
      ComplexVector TmpState1 (Space->GetHilbertSpaceDimension()) ;
      Complex TmpCoefficient;
      
      TmpState.Copy(TmpInitialState);
      Norm = TmpState.Norm();
      double TmpNorm = 1.0;
      TmpExpansionOrder = 0;
      TmpCoefficient = 1.0;
      cout << "Computing state " << (i + 1) << "/" <<  NbrSteps << " at t = " << (TimeStep * i) <<" with Interaction "<< UPotential  <<endl;
      while ( ((fabs(TmpNorm) > 1e-8 ) || (TmpExpansionOrder < 1)) && (TmpExpansionOrder <= Manager.GetInteger("iter-max")))
	{
	  TmpExpansionOrder += 1;
	  TmpCoefficient = -TmpCoefficient * TimeStep * Complex(0.0, 1.0) / ((double) TmpExpansionOrder);
	  VectorHamiltonianMultiplyOperation Operation (Hamiltonian, (&TmpState), (&TmpState1));
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  TmpState.Copy(TmpState1);
	  TmpNorm = sqrt(TmpCoefficient.Re*TmpCoefficient.Re + TmpCoefficient.Im*TmpCoefficient.Im) * TmpState.Norm();
	  AddComplexLinearCombinationOperation Operation1 (&TmpInitialState, &TmpState1, 1, &TmpCoefficient);
	  Operation1.ApplyOperation(Architecture.GetArchitecture());
	  Norm = TmpInitialState.Norm();
	  cout << "Norm = " << Norm << " +/- " << TmpNorm << " for step " << TmpExpansionOrder << endl;      
	} 

      char* OutputName = new char [strlen(OutputNamePrefix)+ strlen(ParameterString) + 16];
      sprintf (OutputName, "%s_%s.%d.vec",OutputNamePrefix,ParameterString,i);
      TmpInitialState.WriteVector(OutputName);      
      delete [] OutputName;
      delete  Hamiltonian;
    }  
  delete[] OutputNamePrefix;
  delete[]  ParameterString;
  delete TightBindingModel;
  return 0;
}

