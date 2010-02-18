#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphere.h"

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <cstring>
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

  // some running options and help
  OptionManager Manager ("FQHESphereBosonsCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('e', "eigenstate", "name of the file containing the eigenstate");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the lz value corresponding to the eigenvector", 0, true, 0);
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 1000);
  (*SystemGroup) += new BooleanOption  ('r', "radians", "set units to radians instead of magnetic lengths", false);
  (*SystemGroup) += new BooleanOption  ('c', "chord", "use chord distance instead of distance on the sphere", false);
  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density insted of density-density correlation", false);
  (*SystemGroup) += new BooleanOption  ('\n', "coefficients-only", "only compute the one or two body coefficients that are requested to evaluate the density-density correlation", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "huge-memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "large-memory", "maximum memory (in kBytes) that can allocated for precalculations when using huge mode", 1);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with rhorho extension");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsCorrelation -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrBosons = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int Lz = Manager.GetInteger("lz-value");
  int NbrPoints = Manager.GetInteger("nbr-points");
  if (Manager.GetString("state") == 0)
    {
      cout << "FQHESphereBosonsCorrelation requires a state" << endl;
      return -1;
    }
  RealVector State;
  char* OutputNameCorr = new char [256 + strlen (Manager.GetString("interaction-name"))];

  BosonOnSphere Space (NbrBosons, Lz, LzMax);
  ParticleOnSphereFunctionBasis Basis(LzMax);

  for (int i = 0; i < Space.GetHilbertSpaceDimension(); ++i)
    {
      Space.PrintState(cout, i) << endl;
    }
  return 0;

  Complex Sum (0.0, 0.0);
  Complex Sum2 (0.0, 0.0);
  Complex TmpValue;
  RealVector Value(2, true);
  double X = 0.0;
  double XInc = M_PI / ((double) NbrPoints);
  Complex* PrecalculatedValues = new Complex [LzMax + 1];
	  
  for (int i = 0; i <= LzMax; ++i)
    {
      Basis.GetFunctionValue(Value, TmpValue, LzMax);
      ParticleOnSphereDensityDensityOperator Operator (&Space, i, LzMax, i, LzMax);
      PrecalculatedValues[i] = Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
    }

  ofstream File;
  File.precision(14);
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = 0;
      if (DensityFlag == false)
	{
	  if (Manager.GetBoolean("coefficients-only"))
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rhorho-c");
	  else
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rhorho");
	}
      else
	{
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rho");
	}
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << Manager.GetString("eigenstate") << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  if (DensityFlag == true)      
    File << "# density correlation coefficients for " << Manager.GetString("eigenstate") << endl;
  else
    File << "# density-density correlation coefficients for " << Manager.GetString("eigenstate") << endl;
  File << "#" << endl << "# (l+S)    n_l" << endl;
  for (int i = 0; i <= LzMax; ++i)
    File << i << " " << PrecalculatedValues[i]<< endl;
  if (CoefficientOnlyFlag == false)
    {
      double Factor1 = (16.0 * M_PI * M_PI) / ((double) (NbrParticles * NbrParticles));
      if (DensityFlag == true)
	Factor1 = 1.0;//4.0 * M_PI;
      double Factor2;
      if (Manager.GetBoolean("radians") == true)
	Factor2 = 1.0;
      else
	Factor2 = sqrt (0.5 * LzMax);
      for (int x = 0; x < NbrPoints; ++x)
	{
	  Value[0] = X;
	  Sum = 0.0;
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      Basis->GetFunctionValue(Value, TmpValue, i);
	      Sum += PrecalculatedValues[i] * (Conj(TmpValue) * TmpValue);
	    }
	  if (ChordFlag == false)
	    File << (X * Factor2) << " " << (Norm(Sum)  * Factor1) << endl;
	  else
	    File << (2.0 * Factor2 * sin (X * 0.5)) << " " << Norm(Sum)  * Factor1 << endl;
	  X += XInc;
	}
    }
  File.close();
 
  delete[] PrecalculatedValues;

  return 0;
}


