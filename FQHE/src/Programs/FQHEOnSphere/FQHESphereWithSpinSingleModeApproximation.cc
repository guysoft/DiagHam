#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSpinAllSz.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"

#include "Operator/ParticleOnSphereWithSpinDensityOperator.h"

#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithm.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithDiskStorage.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundState.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithEigenstates.h"
#include "LanczosAlgorithm/ComplexBasicLanczosAlgorithmWithGroundStateFastDisk.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithm.h"
#include "LanczosAlgorithm/FullReorthogonalizedComplexLanczosAlgorithmWithDiskStorage.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <cstring> 


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereSingleModeApproximation" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  Manager += ToolsGroup;

  (*SystemGroup) += new SingleStringOption('\n', "ground-state", "name of the file containing the ground state vector upon which rho_k acts");
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistics (default value))");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the total sz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new BooleanOption  ('A', "all-sz", "assume a hilbert space including all sz values");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "lz-boost", "Lz momentum that has to be transfer via the SMA", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "spin-indices", "0 for c_{up,M}^+ c_{up,M+boost}, 1 for c_{up,M}^+ c_{down,M+boost}, 2 for c_{down,M}^+ c_{up,M+boost}, 3 for c_{down,M}^+ c_{down,M+boost}");
  (*SystemGroup) += new BooleanOption ('\n', "compute-bilinears", "compute the action of all the bilinear operators on the ground state");
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "sma");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereSingleModeApproximation -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = Manager.GetInteger("total-lz");
  int TotalSz = Manager.GetInteger("total-sz");
  int LzBoost = Manager.GetInteger("lz-boost");
  int SpinIndices = Manager.GetInteger("spin-indices");
  bool Statistics = true;
  if (Manager.GetBoolean("all-sz"))
    TotalSz = -1;
  char* InputStateFile = new char [strlen(Manager.GetString("ground-state")) + 1];
  strcpy(InputStateFile, Manager.GetString("ground-state"));

  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(InputStateFile, NbrParticles, LzMax, TotalLz, 
							   TotalSz, Statistics) == false)
    {
      return -1;
    }

  int Parity = TotalLz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the same parity" << endl;
      return -1;           
    }
  if ((Manager.GetBoolean("all-sz")==false) && ((NbrParticles&1) != (TotalSz&1)))
    {
      cout << "Sz and NbrParticles must have the same parity" << endl;
      return -1;
    }

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
 
  char* OutputNamePrefix = new char [512];
  unsigned long MemorySpace = 9ul << 20;
  ParticleOnSphereWithSpin* InputSpace = 0;
  ParticleOnSphereWithSpin* OutputSpace = 0;
  int ResultingTotalLz = TotalLz + (2 * LzBoost);
  int ResultingTotalSz;
  if ((SpinIndices == 0) || (SpinIndices == 3))
    ResultingTotalSz = TotalSz;
  else if (SpinIndices == 1) //up-down
    ResultingTotalSz += 2;
  else if (SpinIndices == 2) //down-up
    ResultingTotalSz -= 2;
  else
    cout << "SpinIndices invalid value." << endl;

  if (Statistics == true)
    {
      if (Manager.GetBoolean("all-sz")==false)
	{
#ifdef __64_BITS__
	  if (LzMax <= 31)
	    {
	      InputSpace = new FermionOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz, MemorySpace);
	      OutputSpace = new FermionOnSphereWithSpin(NbrParticles, ResultingTotalLz, LzMax, ResultingTotalSz, MemorySpace);
	    }
	  else
	    {
	      cout << "Fermions with Spin not defined yet for LzMax > 31"<<endl;
	      exit(-1);
	    }
#else
	  if (LzMax <= 15)
	    {
	      InputSpace = new FermionOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz, MemorySpace);
	      OutputSpace = new FermionOnSphereWithSpin(NbrParticles, ResultingTotalLz, LzMax, ResultingTotalSz, MemorySpace);
	    }
	  else
	    {
	      cout << "Fermions with Spin not defined yet for LzMax > 15, consider using a 64 bit machine!"<<endl;
	      exit(-1);
	    }
#endif
	}
      else
	{
	  InputSpace = new FermionOnSphereWithSpinAllSz (NbrParticles, TotalLz, LzMax, MemorySpace);
	  OutputSpace = new FermionOnSphereWithSpinAllSz (NbrParticles, ResultingTotalLz, LzMax, MemorySpace);
	}
      sprintf (OutputNamePrefix, "fermions_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, ResultingTotalSz, ResultingTotalLz);
    }
  else
    {
      if (Manager.GetBoolean("all-sz")==false)
	{
	  InputSpace = new BosonOnSphereWithSpin(NbrParticles, TotalLz, LzMax, TotalSz);
	  OutputSpace = new BosonOnSphereWithSpin(NbrParticles, ResultingTotalLz, LzMax, ResultingTotalSz);
	}
      else
	{
	  InputSpace = new BosonOnSphereWithSpinAllSz(NbrParticles, TotalLz, LzMax, MemorySpace);
	  OutputSpace = new BosonOnSphereWithSpinAllSz(NbrParticles, ResultingTotalLz, LzMax, MemorySpace);
	}
      sprintf (OutputNamePrefix, "bosons_sphere_su2_%s_n_%d_2s_%d_sz_%d_lz_%d", Manager.GetString("interaction-name"), NbrParticles, LzMax, ResultingTotalSz, ResultingTotalLz);
    }

  InputSpace->SetTargetSpace(OutputSpace);
  Architecture.GetArchitecture()->SetDimension(InputSpace->GetHilbertSpaceDimension());

  char* StateFileName = Manager.GetString("ground-state");
  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }

  RealVector InputState;
  if (InputState.ReadVector(StateFileName) == false)
    {
      cout << "error while reading " << StateFileName << endl;
      return -1;
    }
  if (InputState.GetVectorDimension() != InputSpace->GetHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions " << InputState.GetVectorDimension() << " "<< InputSpace->GetHilbertSpaceDimension() << endl;
      return -1;
    }
 
  if (Manager.GetBoolean("compute-bilinears"))
    {
      RealVector TmpState(OutputSpace->GetHilbertSpaceDimension());
      int MinLzValue = 0;
      int MaxLzValue = LzMax;
      if (LzBoost >= 0)
	{
	  MaxLzValue = LzMax - LzBoost;	  
	}
      else
	{
	  MinLzValue = -LzBoost;	  
	}
      for (int m = MinLzValue; m <= MaxLzValue; ++m)
	{
	  cout << "computing c^+_"<< (m + LzBoost) << " c_" << m << " |Psi>" << endl;

          ParticleOnSphereWithSpinDensityOperator* TmpOperator;
	  if (SpinIndices == 0) //up-up
	      TmpOperator = new ParticleOnSphereWithSpinDensityOperator(InputSpace, m + LzBoost, 1, m, 1);
	  else if (SpinIndices == 1) //up-down
 	      TmpOperator = new ParticleOnSphereWithSpinDensityOperator(InputSpace, m + LzBoost, 1, m, 0);
	  else if (SpinIndices == 2) //down-up
 	      TmpOperator = new ParticleOnSphereWithSpinDensityOperator(InputSpace, m + LzBoost, 0, m, 1);
	  else //down-down
 	      TmpOperator = new ParticleOnSphereWithSpinDensityOperator(InputSpace, m + LzBoost, 0, m, 0);

	  VectorOperatorMultiplyOperation Operation(TmpOperator, &InputState, &TmpState);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  char* OutputNameLz = new char [strlen(OutputNamePrefix)+ 16];
	  sprintf (OutputNameLz, "%s.%d.vec", OutputNamePrefix, (m - MinLzValue));
	  TmpState.WriteVector(OutputNameLz);
	}
      return 0;
    }


  cout << " Target Hilbert space dimension = " << OutputSpace->GetHilbertSpaceDimension() << endl;
  cout << " Groundstate Hilbert space dimension = " << InputSpace->GetHilbertSpaceDimension() << endl;
 
  return 0;
}
