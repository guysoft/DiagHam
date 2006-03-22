#include "Vector/RealVector.h"

#include "HilbertSpace/QHEHilbertSpace/FermionOnSphereWithSpin.h"

#include "Operator/QHEOperator/ParticleOnSphereWithSpinDensityDensityOperator.h"
#include "FunctionBasis/QHEFunctionBasis/ParticleOnSphereFunctionBasis.h"

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
  OptionManager Manager ("QHEFermionsWithSpinCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the lz value corresponding to the eigenvector", 0, true, 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "szmax", "twice the total z-component of the spin", 0);
  (*SystemGroup) += new SingleIntegerOption  ('S', "SpinCode", "Code for Spin-channel in bits 0bxxxx, x= u/d = 1/0, (0=dddd, 9=uddu, 10=udud, 15=uuuu) ", 0);
  (*SystemGroup) += new SingleStringOption  ('e', "eigenstate", "name of the file containing the eigenstate");
  (*SystemGroup) += new SingleStringOption  ('i', "interaction-name", "name of the interaction (used for output file name)", "sphere_spin");
  (*SystemGroup) += new SingleStringOption ('a', "add-filename", "add a string with additional informations to the output file name(just before the .dat extension)");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 1000);
  (*SystemGroup) += new BooleanOption  ('r', "radians", "set units to radians instead of magnetic lengths", false);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsWithSpinCorrelation -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int Sz = ((SingleIntegerOption*) Manager["szmax"])->GetInteger();
  int Lz = ((SingleIntegerOption*) Manager["lz-value"])->GetInteger();
  int NbrPoints = ((SingleIntegerOption*) Manager["nbr-points"])->GetInteger();
  unsigned SpinCode = (unsigned) ((SingleIntegerOption*) Manager["SpinCode"])->GetInteger();
  if (((SingleStringOption*) Manager["eigenstate"])->GetString() == 0)
    {
      cout << "QHEFermionsWithSpinCorrelation requires a state" << endl;
      return -1;
    }
  char typeStr[16]="corr-";
  if (SpinCode&8u) sprintf(typeStr,"%su",typeStr); else sprintf(typeStr,"%sd",typeStr);
  if (SpinCode&4u) sprintf(typeStr,"%su",typeStr); else sprintf(typeStr,"%sd",typeStr);
  sprintf(typeStr,"%s-",typeStr);
  if (SpinCode&2u) sprintf(typeStr,"%su",typeStr); else sprintf(typeStr,"%sd",typeStr);
  if (SpinCode&1u) sprintf(typeStr,"%su",typeStr); else sprintf(typeStr,"%sd",typeStr);
  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["eigenstate"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["eigenstate"])->GetString() << endl;
      return -1;      
    }
  char* OutputNameCorr = new char [256 + strlen (((SingleStringOption*) Manager["interaction-name"])->GetString())];
  if (((SingleStringOption*) Manager["add-filename"])->GetString() == 0)
    {
      sprintf(OutputNameCorr,"%s",((SingleStringOption*) Manager["eigenstate"])->GetString());
      OutputNameCorr[strlen(OutputNameCorr)-4]='\0';
      sprintf(OutputNameCorr,"%s.%s.dat",OutputNameCorr,typeStr);
    }
  else
    {
      sprintf (OutputNameCorr, "fermions_%s_n_%d_2s_%d_Sz_%d_lz_%d_%s.%s.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrFermions, LzMax, Sz, Lz, 
	       ((SingleStringOption*) Manager["add-filename"])->GetString(), typeStr);
    }

      ParticleOnSphereWithSpin* Space;
#ifdef __64_BITS__
      if (LzMax <= 31)
        {
          Space = new FermionOnSphereWithSpin(NbrFermions, Lz, LzMax);
        }
      else
        {
	  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	  return -1;
        }
#else
      if (LzMax <= 15)
        {
          Space = new FermionOnSphereWithSpin(NbrFermions, Lz, LzMax, Sz, 0);
        }
      else
        {
	  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	  return -1;
        }
#endif
  ParticleOnSphereFunctionBasis Basis(LzMax);

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
      ParticleOnSphereWithSpinDensityDensityOperator Operator (Space, i, LzMax, i, LzMax,SpinCode);
      PrecalculatedValues[i] = Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
    }
  ofstream File;
  File.precision(14);
  File.open(OutputNameCorr, ios::binary | ios::out);
  double Factor1 = (16.0 * M_PI * M_PI) / ((double) (NbrFermions * NbrFermions));
  double Factor2 = sqrt (0.5 * LzMax );
  if (((BooleanOption*) Manager["radians"])->GetBoolean() == true) Factor2 = 1.0;
  for (int x = 0; x < NbrPoints; ++x)
    {
      Value[0] = X;
      int Pos = 0;
      Sum = 0.0;
      for (int i = 0; i <= LzMax; ++i)
	{
	  Basis.GetFunctionValue(Value, TmpValue, i);
	  Sum += PrecalculatedValues[Pos] * (Conj(TmpValue) * TmpValue);
	  ++Pos;
	}
      File << (X * Factor2) << " " << Norm(Sum)  * Factor1 << endl;
      X += XInc;
    }
  File.close();

  delete[] OutputNameCorr;	  
  delete[] PrecalculatedValues;

  return 0;
}


