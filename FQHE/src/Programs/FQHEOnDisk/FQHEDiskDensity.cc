#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnDisk.h"
#include "HilbertSpace/FermionOnDisk.h"
#include "HilbertSpace/FermionOnDiskUnlimited.h"
#include "HilbertSpace/ParticleOnDisk.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Operator/ParticleOnDiskDensityOperator.h"
#include "FunctionBasis/ParticleOnDiskFunctionBasis.h"

#include "Tools/FQHEFiles/FQHEOnDiskFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHEDiskDensity" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PlotOptionGroup = new OptionGroup ("plot options");  
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += SystemGroup;
  Manager += PlotOptionGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (overriding the one found in the vector file name if greater than 0)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "momentum", "single particle momentum (overriding the one found in the vector file name if greater than 0)", 0, true, 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "force-maxmomentum", "force the maximum single particle momentum to a particular value (overriding the one found in the vector file name if greater than 0)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use boson statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermion statistics");
  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the vector file describing the state whose density has to be plotted");
  (*SystemGroup) += new BooleanOption  ('\n', "coefficients-only", "only compute the one body coefficients that are requested to evaluate the density profile", false);

  (*PlotOptionGroup) += new SingleStringOption ('\n', "output", "output file name", "density.dat");
  (*PlotOptionGroup) += new SingleIntegerOption ('\n', "nbr-samples", "number of samples in radial direction", 1000, true, 10);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskDensity -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  if (((SingleStringOption*) Manager["state"])->GetString() == 0)
    {
      cout << "QHEBosonsCorrelation requires a state" << endl;
      return -1;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int Lz = ((SingleIntegerOption*) Manager["momentum"])->GetInteger();
  int ForceMaxMomentum = ((SingleIntegerOption*) Manager["force-maxmomentum"])->GetInteger();
  bool CoefficientOnlyFlag = ((BooleanOption*) Manager["coefficients-only"])->GetBoolean();
  char* OutputName = ((SingleStringOption*) Manager["output"])->GetString();
  int NbrSamples = ((SingleIntegerOption*) Manager["nbr-samples"])->GetInteger();
  bool Statistics = true;

  if (FQHEOnDiskFindSystemInfoFromFileName(((SingleStringOption*) Manager["state"])->GetString(), NbrParticles, ForceMaxMomentum, Lz, Statistics) == false)
    {
      return -1;      
    }
  if ((((BooleanOption*) Manager["boson"])->GetBoolean() == true) || (((BooleanOption*) Manager["fermion"])->GetBoolean() == true))
    {
      if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
        Statistics = false;
      else
        Statistics = true;
    }
  
  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["state"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["state"])->GetString() << endl;
      return -1;      
    }

  ParticleOnDiskFunctionBasis Basis(Lz);
  ParticleOnDisk* Space = 0;
  int TmpMaxMomentum = 0;
  if (Statistics == false)
    {
      TmpMaxMomentum = Lz;
      if ((ForceMaxMomentum >= 0) && (ForceMaxMomentum < TmpMaxMomentum))
	TmpMaxMomentum = ForceMaxMomentum;      
      Space = new BosonOnDisk(NbrParticles, Lz);
    }
  else
    {
      TmpMaxMomentum = (Lz - (((NbrParticles - 1) * (NbrParticles - 2)) / 2));
      if ((ForceMaxMomentum >= 0) && (ForceMaxMomentum < TmpMaxMomentum))
	TmpMaxMomentum = ForceMaxMomentum;
#ifdef __64_BITS__
      if (TmpMaxMomentum < 63)      
#else
      if (TmpMaxMomentum < 31)
#endif
	Space = new FermionOnDisk (NbrParticles, Lz, TmpMaxMomentum);
      else
	Space = new FermionOnDiskUnlimited (NbrParticles, Lz, TmpMaxMomentum);      
    }

  if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
    {
      cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }

  Complex* PrecalculatedValues = new Complex [TmpMaxMomentum + 1];
	  
  for (int i = 0; i <= TmpMaxMomentum; ++i)
    {
      ParticleOnDiskDensityOperator Operator (Space, i);
      PrecalculatedValues[i] = Operator.MatrixElement(State, State);
    }
  
  ofstream File;
  File.precision(14);
  File.open(OutputName, ios::binary | ios::out);
  if (CoefficientOnlyFlag == false)
    {
      Complex Tmp (0.0);
      double RInc = 1.2 * ((double) Lz) / ((double) NbrSamples);
      NbrSamples *= 2;      
      RealVector Value(2, true);
      Value[0] = - 1.2 * ((double) Lz);
      Complex TmpValue;
      for (int i = 0; i <= NbrSamples; ++i)
	{
	  Tmp = 0.0;
	  for (int j = 0; j <= Lz; ++j)
	    {
	      Basis.GetFunctionValue(Value, TmpValue, j);
	      Tmp += PrecalculatedValues[j] * SqrNorm(TmpValue);
	    }
	  File << Value[0] << " " << (Tmp.Re * exp (-0.5 * (Value[0] * Value[0]))) << endl;
	  Value[0] += RInc;
	}
    }
  else
    {
      File << "# density coefficients for " << ((SingleStringOption*) Manager["state"])->GetString() << endl;
      File << "#" << endl << "# m    c_m" << endl;
      for (int i = 0; i <= TmpMaxMomentum; ++i)
        File << i << " " << PrecalculatedValues[i].Re << endl;
    }
  File.close();

  delete[] PrecalculatedValues;
  delete Space;
}
