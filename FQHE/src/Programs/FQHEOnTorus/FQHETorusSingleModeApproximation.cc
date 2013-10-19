#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Vector/ComplexVector.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"

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
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

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
  OptionManager Manager ("FQHETorusFermionsSingleModeApproximation" , "0.01");
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

  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistics (default value))");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "max-momentum", "maximum momentum for a single particle", 18);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "sma");
  (*SystemGroup) += new SingleIntegerOption  ('y', "y-momentum", "total momentum in the y direction of the ground state (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "kx", "momentum along x direction of the density operator (negative if none)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ky", "momentum along y direction of the density operator (negative if none)", -1);

  (*SystemGroup) += new SingleDoubleOption   ('r', "ratio", 
					      "ratio between lengths along the x and y directions (-1 if has to be taken equal to nbr-particles/4)", 
					      -1);
  (*SystemGroup) += new SingleDoubleOption   ('c', "costheta", "cosine of the angle between the sides of torus (between 0 and 1, 0 for rectangular cell)", 0.0);
  (*MiscGroup) += new SingleStringOption('\n', "ground-state", "name of the file containing the ground state vector upon which rho_k acts");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 
						      500);
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

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int MaxMomentum = ((SingleIntegerOption*) Manager["max-momentum"])->GetInteger();
  int YMomentum = ((SingleIntegerOption*) Manager["y-momentum"])->GetInteger();
  int Kx = ((SingleIntegerOption*) Manager["kx"])->GetInteger();
  int Ky = ((SingleIntegerOption*) Manager["ky"])->GetInteger();
  double XRatio = ((SingleDoubleOption*) Manager["ratio"])->GetDouble();
  double CosTheta = ((SingleDoubleOption*) Manager["costheta"])->GetDouble();
  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  int MomentumModulo = FindGCD(NbrParticles, MaxMomentum);

  int ResultingYMomentum = (YMomentum + Ky) % MaxMomentum; 

  char* OutputNameLz = new char [512];
  ParticleOnTorus* TotalSpace = 0;
  ParticleOnTorus* TargetSpace = 0;
  if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
    {
#ifdef  __64_BITS__
      if ((MaxMomentum + NbrParticles - 1) < 63)
#else
	if ((MaxMomentum + NbrParticles - 1) < 31)	
#endif
	  {
	    TotalSpace = new BosonOnTorusShort(NbrParticles, MaxMomentum, YMomentum);
	    TargetSpace = new BosonOnTorusShort(NbrParticles, MaxMomentum, ResultingYMomentum);
            ((BosonOnTorusShort*)TotalSpace)->SetTargetSpace(TargetSpace);
	  }
	else
	  {
	    TotalSpace = new BosonOnTorus(NbrParticles, MaxMomentum, YMomentum);
	    TargetSpace = new BosonOnTorusShort(NbrParticles, MaxMomentum, ResultingYMomentum);
            ((BosonOnTorus*)TotalSpace)->SetTargetSpace(TargetSpace);
	  }

        sprintf (OutputNameLz, "bosons_sma_n_%d_2s_%d_ky_%d.0.vec", NbrParticles, MaxMomentum, ResultingYMomentum);
    }
  else
    {
      TotalSpace = new FermionOnTorus (NbrParticles, MaxMomentum, YMomentum);
      TargetSpace = new FermionOnTorus (NbrParticles, MaxMomentum, ResultingYMomentum);
      ((FermionOnTorus*)TotalSpace)->SetTargetSpace(TargetSpace);

      sprintf (OutputNameLz, "fermions_sma_n_%d_2s_%d_ky_%d.0.vec", NbrParticles, MaxMomentum, ResultingYMomentum);
    }




  cout << "*******************************************************************"<<endl;
  cout << " Calculating Psi_k = rho_k Psi_0, where rho_k = sum_i exp(ik.R_i) " << endl;
  cout << " Ground state is at Ky= " << YMomentum << " Resulting state is at Ky= " << ResultingYMomentum << endl;
  cout << "*******************************************************************"<<endl;

  cout << " Target Hilbert space dimension = " << TargetSpace->GetHilbertSpaceDimension() << endl;

  cout << " Groundstate Hilbert space dimension = " << TotalSpace->GetHilbertSpaceDimension() << endl;
 
  Architecture.GetArchitecture()->SetDimension(TotalSpace->GetHilbertSpaceDimension());

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
  ComplexVector State(InputState.GetVectorDimension(), true);
  for (int i = 0; i < InputState.GetVectorDimension(); i++)
   {
     Complex Tmp;
     Tmp.Re = InputState[i];
     Tmp.Im = 0.0;
     State[i] = Tmp;
   }
  

  if (State.GetVectorDimension() != TotalSpace->GetHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions "<<State.GetVectorDimension() << " "<< TotalSpace->GetHilbertSpaceDimension() << endl;
      return -1;
    }

  double InvRatio = 1.0 / XRatio;
  double SinTheta = sqrt(1.0 - CosTheta * CosTheta);
  double Lx = sqrt(2.0 * M_PI * (double)MaxMomentum * XRatio/SinTheta);
  double Ly = sqrt(2.0 * M_PI * (double)MaxMomentum * InvRatio/SinTheta);
  double Gx = 2.0 * M_PI / Lx;
  double Gy = 2.0 * M_PI / Ly;
  cout << "-------------------------------------------------"<<endl;
  cout << "-     Geometry: Lx = " << Lx << " , Ly = " << Ly<< "  ; cos angle =   " << CosTheta << " -" << endl;
  cout << "-------------------------------------------------"<<endl;

  Complex Phase, MatEl;
  ComplexVector TmpState(TargetSpace->GetHilbertSpaceDimension(), true);
  int Index;
  double Coefficient;

  for (int m = 0; m < MaxMomentum; ++m)
   {
      MatEl.Re = cos(0.5 * (2.0 * M_PI/(double)MaxMomentum) * Kx * (2.0 * m + Ky));
      MatEl.Im = -sin(0.5 * (2.0 * M_PI/(double)MaxMomentum) * Kx * (2.0 * m + Ky));
      MatEl /= sqrt(NbrParticles);
      cout << "m= " << m << " Mat el " << MatEl << endl;
      for (int i = 0; i < TotalSpace->GetHilbertSpaceDimension(); ++i)
        {
           Index = TotalSpace->AdA(i, (m + Ky)%MaxMomentum, m, Coefficient);
           if ((Index < TargetSpace->GetHilbertSpaceDimension()) && (Coefficient != 0))
             {
                TmpState[Index] += (MatEl * Coefficient * State[i]);
             }
        }
   }
  
  cout << "check the norm: " << endl;
  double TmpNorm = 0.0;
  for (int i = 0; i < TotalSpace->GetHilbertSpaceDimension(); i++)
   {
     Complex Tmp = TmpState[i];
     TmpNorm += (Tmp.Re * Tmp.Re + Tmp.Im * Tmp.Im);
   }
  cout << "Norm " << TmpNorm << endl;
  if (TmpNorm > 1e-10) 
    TmpState /= sqrt(TmpNorm);
  else
    cout << "Warning: Norm " << TmpNorm << endl;
  

  TmpState.WriteVector(OutputNameLz);

   
  return 0;
}
