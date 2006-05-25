#include "Vector/RealVector.h"

#include "HilbertSpace/QHEHilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/QHEHilbertSpace/FermionOnSphereUnlimited.h"

#include "Operator/QHEOperator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/QHEOperator/ParticleOnSphereDensityOperator.h"
#include "FunctionBasis/QHEFunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/QHEFunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

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
  OptionManager Manager ("QHEFermionsCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the lz value corresponding to the eigenvector", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "landau-level", "index of the Landau level (0 being the LLL)", 0);
  (*SystemGroup) += new SingleStringOption  ('s', "state", "name of the file containing the eigenstate");
  (*SystemGroup) += new SingleStringOption  ('i', "interaction-name", "name of the interaction (used for output file name)", "laplaciandelta");
  (*SystemGroup) += new SingleStringOption ('a', "add-filename", "add a string with additional informations to the output file name(just before the .dat extension)");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 1000);
  (*SystemGroup) += new BooleanOption  ('r', "radians", "set units to radians instead of magnetic lengths", false);
  (*SystemGroup) += new BooleanOption  ('c', "chord", "use chord distance instead of distance on the sphere", false);
  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density insted of density-density correlation", false);
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrize", "use symmetrize combination of the lz and -lz eigenstate (assuming lz -lz symmetry)", false);
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsCorrelation -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int Lz = ((SingleIntegerOption*) Manager["lz-value"])->GetInteger();
  bool InverseLzFlag = false;
  if (Lz < 0)
    {
      InverseLzFlag = true;
      Lz *= -1;
    }
  int LandauLevel = ((SingleIntegerOption*) Manager["landau-level"])->GetInteger();
  int NbrPoints = ((SingleIntegerOption*) Manager["nbr-points"])->GetInteger();
  bool DensityFlag = ((BooleanOption*) Manager["density"])->GetBoolean();
  bool ChordFlag = ((BooleanOption*) Manager["chord"])->GetBoolean();

  if (((SingleStringOption*) Manager["state"])->GetString() == 0)
    {
      cout << "QHEFermionsCorrelation requires a state" << endl;
      return -1;
    }
  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["state"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["state"])->GetString() << endl;
      return -1;      
    }
  char* OutputNameCorr = new char [256 + strlen (((SingleStringOption*) Manager["interaction-name"])->GetString())];
  if (DensityFlag == false)
    if (((SingleStringOption*) Manager["add-filename"])->GetString() == 0)
      {
	sprintf (OutputNameCorr, "fermions_%s_n_%d_2s_%d.rho_rho.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrFermions, LzMax);
      }
    else
      {
	sprintf (OutputNameCorr, "fermions_%s_n_%d_2s_%d_%s.rho_rho.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrFermions, LzMax,
		 ((SingleStringOption*) Manager["add-filename"])->GetString());
      }
  else
    if (((SingleStringOption*) Manager["add-filename"])->GetString() == 0)
      {
	sprintf (OutputNameCorr, "fermions_%s_n_%d_2s_%d.rho.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrFermions, LzMax);
      }
    else
      {
	sprintf (OutputNameCorr, "fermions_%s_n_%d_2s_%d_%s.rho.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), NbrFermions, LzMax,
		 ((SingleStringOption*) Manager["add-filename"])->GetString());
      }

  ParticleOnSphere* Space;
#ifdef __64_BITS__
  if (LzMax <= 63)
    {
      Space = new FermionOnSphere(NbrFermions, Lz, LzMax);
    }
  else
    {
      Space = new FermionOnSphereUnlimited(NbrFermions, Lz, LzMax);
    }
#else
  if (LzMax <= 31)
    {
      Space = new FermionOnSphere(NbrFermions, Lz, LzMax);
    }
  else
    {
      Space = new FermionOnSphereUnlimited(NbrFermions, Lz, LzMax);
    }
#endif
  AbstractFunctionBasis* Basis;
  if (LandauLevel == 0)
    Basis = new ParticleOnSphereFunctionBasis(LzMax);
  else
    Basis = new ParticleOnSphereGenericLLFunctionBasis(LzMax - (2 * LandauLevel), LandauLevel);

  Complex Sum (0.0, 0.0);
  Complex Sum2 (0.0, 0.0);
  Complex TmpValue;
  RealVector Value(2, true);
  double X = 0.0;
  double XInc = M_PI / ((double) NbrPoints);

  Complex* SymmetrizePrecalculatedValues = new Complex [LzMax + 1];
  if ((DensityFlag == true) && (((BooleanOption*) Manager["symmetrize"])->GetBoolean() == true))
    {
      RealVector SymmetricVector;
      SymmetricVector.Copy(State);
      ParticleOnSphere* TargetSpace;
#ifdef __64_BITS__
      if (LzMax <= 63)
	{
	  TargetSpace = new FermionOnSphere(NbrFermions, -Lz, LzMax);
	}
      else
	{
	  TargetSpace = new FermionOnSphereUnlimited(NbrFermions, -Lz, LzMax);
	}
#else
      if (LzMax <= 31)
	{
	  TargetSpace = new FermionOnSphere(NbrFermions, -Lz, LzMax);
	}
      else
	{
	  TargetSpace = new FermionOnSphereUnlimited(NbrFermions, -Lz, LzMax);
	}
#endif
      Space->SetTargetSpace(TargetSpace);
      int TmpIndex;
       for (int i = 0; i <= LzMax; ++i)
 	{
	  TmpIndex = abs(LzMax + i - Lz);
	  if (TmpIndex <= LzMax)
 	    {
 	      ParticleOnSphereDensityOperator Operator (Space, TmpIndex, i);
 	      SymmetrizePrecalculatedValues[i] = Operator.MatrixElement(SymmetricVector, State);
 	    }
	  else
	    SymmetrizePrecalculatedValues[i] = 0.0;
	  cout << i << " " << SymmetrizePrecalculatedValues[i] << endl;
 	}      
    }
  Complex* PrecalculatedValues = new Complex [LzMax + 1];
  if (DensityFlag == false)
    for (int i = 0; i <= LzMax; ++i)
      {
	Basis->GetFunctionValue(Value, TmpValue, LzMax);
	ParticleOnSphereDensityDensityOperator Operator (Space, i, LzMax, i, LzMax);
	PrecalculatedValues[i] = Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
      }
  else
    for (int i = 0; i <= LzMax; ++i)
      {
	ParticleOnSphereDensityOperator Operator (Space, i);
	PrecalculatedValues[i] = Operator.MatrixElement(State, State);
      }
  ofstream File;
  File.precision(14);
  File.open(OutputNameCorr, ios::binary | ios::out);
  double Factor1 = (16.0 * M_PI * M_PI) / ((double) (NbrFermions * NbrFermions));
  if (DensityFlag == true)
    Factor1 = 1.0;//4.0 * M_PI;
  double Factor2;
  if (((BooleanOption*) Manager["radians"])->GetBoolean() == true)
    Factor2 = 1.0;
  else
    Factor2 = sqrt (0.5 * LzMax);
  for (int x = 0; x < NbrPoints; ++x)
    {
      Value[0] = X;
      int Pos = 0;
      Sum = 0.0;
      if (((BooleanOption*) Manager["symmetrize"])->GetBoolean() == true)
	{
	  Complex TmpValue2;
	  int TmpIndex;
	  for (int i = 0; i <= LzMax; ++i)
	    {
	      Basis->GetFunctionValue(Value, TmpValue, LzMax - i);	      
	      Basis->GetFunctionValue(Value, TmpValue2, i);	    
	      Sum += 0.5 * PrecalculatedValues[i] * ((Conj(TmpValue) * TmpValue) + (Conj(TmpValue2) * TmpValue2));
	      TmpIndex = LzMax + i - Lz;
	      if ((TmpIndex >= 0) && (TmpIndex <= LzMax))
		{
		  Basis->GetFunctionValue(Value, TmpValue, TmpIndex);	      
		  Sum -= 0.5 * SymmetrizePrecalculatedValues[i] * (Conj(TmpValue2) * TmpValue);
		}
	      else		
		{
		  TmpIndex *= -1;
		  if ((TmpIndex >= 0) && (TmpIndex <= LzMax))
		    {
		      Basis->GetFunctionValue(Value, TmpValue, TmpIndex);	      
		      Sum -= 0.5 * SymmetrizePrecalculatedValues[i] * (Conj(TmpValue) * TmpValue2);
		    }		  
		}
	    }
	}
      else
	for (int i = 0; i <= LzMax; ++i)
	  {
	    Basis->GetFunctionValue(Value, TmpValue, i);
	    Sum += PrecalculatedValues[Pos] * (Conj(TmpValue) * TmpValue);
	    ++Pos;
	  }
      if (ChordFlag == false)
	File << (X * Factor2) << " " << Norm(Sum)  * Factor1 << endl;
      else
	File << (2.0 * Factor2 * sin (X * 0.5)) << " " << Norm(Sum)  * Factor1 << endl;
      X += XInc;
    }


  delete[] OutputNameCorr;	  
  delete[] PrecalculatedValues;
//  delete Basis;

  return 0;
}


