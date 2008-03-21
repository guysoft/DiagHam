#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include "Options/Options.h"

#include "Tools/FQHEWaveFunction/FQHESphereSymmetrizedSUKToU1WaveFunction.h"
#include "Tools/FQHEWaveFunction/HalperinOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/JainCFFilledLevelOnSphereWaveFunction.h"

#include "Vector/ComplexVector.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


void RandomUV (ComplexVector& uv, int nbrParticles, StdlibRandomNumberGenerator* RandomNumberGenerator);

void FlipCoordinates (ComplexVector& uv, int i, int j);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereSUKToU1MCOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;


  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 10);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 9);
  (*SystemGroup) += new SingleIntegerOption  ('k', "k-value", "k index of the SU(k) symmetry group", 2);
  (*SystemGroup) += new SingleStringOption  ('\n', "exact-state", "name of the file containing the vector obtained using exact diagonalization");
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  
  (*SystemGroup) += new SingleStringOption  ('\n', "use-exact", "file name of an exact state that has to be used as test wave function");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereSUKToU1MCOverlap -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int KValue = Manager.GetInteger("k-value");

  HalperinOnSphereWaveFunction* BaseFunction = new HalperinOnSphereWaveFunction (NbrParticles >> 1, NbrParticles >> 1, 2, 2, 1);
  FQHESphereSymmetrizedSUKToU1WaveFunction* SymmetrizedFunction = new FQHESphereSymmetrizedSUKToU1WaveFunction (NbrParticles, KValue, BaseFunction, true);
  PfaffianOnSphereWaveFunction* ReferenceFunction = new PfaffianOnSphereWaveFunction(NbrParticles);
  JainCFFilledLevelOnSphereWaveFunction* CFFunction = new JainCFFilledLevelOnSphereWaveFunction(NbrParticles, 2, 2);

  StdlibRandomNumberGenerator RandomNumberGenerator(29457);
  Complex Sum = 0.0;
  double RefNorm = 0.0;
  double TestNorm = 0.0;
  for (int i = 0; i < 100000; ++i)
    {
      ComplexVector TmpUV (NbrParticles * 2, true);
      RandomUV (TmpUV, NbrParticles, &RandomNumberGenerator);      
      Complex Tmp1 = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
      Complex Tmp2 = CFFunction->CalculateFromSpinorVariables(TmpUV);
      Sum += (Conj(Tmp2) * Tmp1);
      RefNorm +=  SqrNorm(Tmp2);
      TestNorm +=  SqrNorm(Tmp1);
//       cout << SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV) << endl;;
//       cout << ReferenceFunction->CalculateFromSpinorVariables(TmpUV) << endl;;      
      if ((i % 100) == 0)
	cout << (i / 100) << " " << (SqrNorm(Sum) / (RefNorm * TestNorm)) << endl;
//       cout << SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV) << endl;;
//       cout << ReferenceFunction->CalculateFromSpinorVariables(TmpUV) << endl;;      
    }
  
//   ComplexVector UV (NbrParticles * 2, true);
//   RandomUV (UV, NbrParticles, &RandomNumberGenerator);
//   cout << SymmetrizedFunction->CalculateFromSpinorVariables(UV) << endl;;
      
//   for (int i = 0; i < NbrParticles; ++i)
//     {
//       for (int j = i + 1; j < NbrParticles; ++j)
// 	{
//  	  cout << SymmetrizedFunction->CalculateFromSpinorVariables(UV) << " | " << CFFunction->CalculateFromSpinorVariables(UV) << " | " ;
//  	  FlipCoordinates(UV, i, j);
// 	  cout << SymmetrizedFunction->CalculateFromSpinorVariables(UV) << " | " << CFFunction->CalculateFromSpinorVariables(UV) << endl;
// 	  FlipCoordinates(UV, i, j);

// 	}
//     }
  return 0;
}

void RandomUV (ComplexVector& uv, int nbrParticles, StdlibRandomNumberGenerator* RandomNumberGenerator)
{
  for (int j = 0; j < nbrParticles; ++j)
    {
      double x = acos (1.0 - (2.0 * RandomNumberGenerator->GetRealRandomNumber()));
      double y = 2.0 * M_PI * RandomNumberGenerator->GetRealRandomNumber();
      uv.Re(2 * j) = cos(0.5 * x);
      uv.Im(2 * j) = uv.Re(2 * j) * sin(0.5 * y);
      uv.Re(2 * j) *= cos(0.5 * y);
      uv.Re(2 * j + 1) = sin(0.5 * x);
      uv.Im(2 * j + 1) = - uv.Re(2 * j + 1) * sin(0.5 * y);
      uv.Re(2 * j + 1) *= cos(0.5 * y);      
    }
}

void FlipCoordinates (ComplexVector& uv, int i, int j)
{
  Complex Tmp = uv[2 * i];
  uv.Re(2 * i) = uv.Re(2 * j);
  uv.Re(2 * j) = Tmp.Re;
  uv.Im(2 * i) = uv.Im(2 * j);
  uv.Im(2 * j) = Tmp.Im;
  Tmp = uv[2 * i + 1];
  uv.Re(2 * i + 1) = uv.Re(2 * j + 1);
  uv.Re(2 * j + 1) = Tmp.Re;
  uv.Im(2 * i + 1) = uv.Im(2 * j + 1);
  uv.Im(2 * j + 1) = Tmp.Im;
}
