#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/RandomNumber/FileRandomNumberGenerator.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"

#include "Options/Options.h"

#include "Tools/FQHEWaveFunction/FQHESphereSymmetrizedSUKToU1WaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESU2HalperinPermanentOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESU3HalperinPermanentOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESU3GeneralizedGaffnianOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESU4HalperinPermanentOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/JainCFFilledLevelOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/HundRuleCFStates.h"
#include "Tools/FQHEWaveFunction/SU3HalperinOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/MooreReadOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereTwoQuasiholeWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereTwoQuasielectronWaveFunction.h"

#include "Vector/ComplexVector.h"

#include "GeneralTools/ConfigurationParser.h"

#include "GeneralTools/Endian.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


void RandomUV (ComplexVector& uv, RealVector& positions, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator);

void RandomUVOneCoordinate(ComplexVector& uv, RealVector& positions, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator);

void FlipCoordinates (ComplexVector& uv, int i, int j);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereDensityMC" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MonteCarloGroup = new OptionGroup ("Monte Carlo options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += MonteCarloGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;


  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 10);
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  
  (*SystemGroup) += new BooleanOption ('\n', "test-symmetry", "check the test wave function is symmetric/antisymmetric");  
  (*SystemGroup) += new SingleStringOption  ('\n', "load-permutations", "read all the permutations needed to compute the reference wave function from a file");  
  (*SystemGroup) += new SingleStringOption  ('\n', "save-permutations", "file name where all the permutations needed to compute the reference wave function have to be stored");  
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistic");
  (*SystemGroup) += new BooleanOption ('\n', "quasielectron", "plot quasiparticles instead of quasiholes");  
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output file name", "density.dat"); 

  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "nbr-orbitals", "number of orbitals used to sample the sphere geometry", 100);
  (*MonteCarloGroup) += new SingleIntegerOption  ('s', "nbr-step", "number of steps for the density profil", 100);
  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "nbr-warmup", "number of Monte Carlo iterations that have to be done before evaluating the energy (i.e. warm up sequence)", 10000);
  (*MonteCarloGroup) += new BooleanOption  ('r', "resume", "resume from a previous run");
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "record-step", "number of iteration between two consecutive result recording of energy value (0 if no on-disk recording is needed)", 0);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "record-file", "name of the file where energy recording has to be done", "montecarlo.dat");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "with-timecoherence", "use time coherence between two successive evaluation of the wave function");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "show-details", "show intermediate values used for overlap calculation", false);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "random-file", "name of the file where random number to use are stored (use internal random generator if no file name is provided)");
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "random-seek", "if usage of a random number file is activiated, jump the first random numbers up to the seek position", 0);
  (*MonteCarloGroup) += new BooleanOption  ('\n', "weight-symmetrized" , "use the norm of the symmetrized wave fonction as probalbility density instead of the exact state");
  (*MonteCarloGroup) += new  SingleStringOption ('\n', "record-wavefunctions", "optional file where each wavefunctions will be tabulated and recorded");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereDensityMC -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  bool OverlapFlag = true;
  if (((BooleanOption*) Manager["test-symmetry"])->GetBoolean() == true)
    {
      OverlapFlag = false;
    }
  int NbrWarmUpIter = ((SingleIntegerOption*) Manager["nbr-warmup"])->GetInteger();
  int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  int NbrSteps = ((SingleIntegerOption*) Manager["nbr-step"])->GetInteger();
  int NbrOrbitals = ((SingleIntegerOption*) Manager["nbr-orbitals"])->GetInteger();
  bool ResumeFlag = Manager.GetBoolean("resume");
  bool StatisticFlag = !(((BooleanOption*) Manager["boson"])->GetBoolean());
  bool QuasielectronFlag = (((BooleanOption*) Manager["quasielectron"])->GetBoolean());
  char* RecordWaveFunctions = ((SingleStringOption*) Manager["record-wavefunctions"])->GetString();
  if (RecordWaveFunctions != 0)
    {
      ofstream RecordFile;
      RecordFile.open(RecordWaveFunctions, ios::out | ios::binary);
      RecordFile.close();
    }


  int LzMax = 2 * NbrParticles - 2;
  if (QuasielectronFlag == true)
    LzMax -= 2;

//   AbstractQHEParticle* ExactSpace = new FermionOnSphere (NbrParticles, 0, LzMax);
//   RealVector ExactState;
//   ExactState.ReadVector ("fermions_hardcore_nbody_3_n_6_2s_9_lz_0.0.vec");

  Abstract1DComplexFunctionOnSphere* SymmetrizedFunction = 0;
  if (QuasielectronFlag == true)
    {
      if (((SingleStringOption*) Manager["load-permutations"])->GetString() == 0)
	SymmetrizedFunction = new PfaffianOnSphereTwoQuasielectronWaveFunction(NbrParticles, 0.0, 0.0, M_PI, 0.0, true);
      else
	SymmetrizedFunction = new PfaffianOnSphereTwoQuasielectronWaveFunction(((SingleStringOption*) Manager["load-permutations"])->GetString(),
									       0.0, 0.0, M_PI, 0.0, true);
      if (((SingleStringOption*) Manager["save-permutations"])->GetString() != 0)
	{
	  ((PfaffianOnSphereTwoQuasielectronWaveFunction*) SymmetrizedFunction)->WritePermutations(((SingleStringOption*) Manager["save-permutations"])->GetString());
	  return 0;
	}

    }
  else
    {
      SymmetrizedFunction = new PfaffianOnSphereTwoQuasiholeWaveFunction(NbrParticles, 0.0, 0.0, M_PI, 0.0, true);
    }
  Abstract1DComplexFunctionOnSphere* TestFunction;
  AbstractRandomNumberGenerator* RandomNumber = 0;
  if (((SingleStringOption*) Manager["random-file"])->GetString() != 0)
    {
      if (ResumeFlag == true)
	{
	  ifstream MCState;
	  MCState.open("mcstate.dat", ios::in | ios::binary);
	  int TmpNbrIter;
	  ReadLittleEndian(MCState, TmpNbrIter);
	  unsigned long TmpNumber;
	  ReadLittleEndian(MCState, TmpNumber);
	  MCState.close();
	  RandomNumber = new FileRandomNumberGenerator(((SingleStringOption*) Manager["random-file"])->GetString(), TmpNumber, 
						       ((SingleIntegerOption*) Manager["random-seek"])->GetInteger());	  
	}
      else
	RandomNumber = new FileRandomNumberGenerator(((SingleStringOption*) Manager["random-file"])->GetString(), (NbrWarmUpIter * 4) + (NbrIter * 4) + 2000, 
						     ((SingleIntegerOption*) Manager["random-seek"])->GetInteger());
    }
  else
    {
      //      RandomNumber = new StdlibRandomNumberGenerator (29457);
      RandomNumber = new NumRecRandomGenerator(29457);
    }

  if (OverlapFlag == true)
    {
      RealVector TmpPositions (NbrParticles * 2, true);
      RealVector PreviousTmpPositions (NbrParticles * 2, true);
      ComplexVector TmpUV (NbrParticles * 2, true);
      ComplexVector PreviousTmpUV (NbrParticles * 2, true);
      
      ParticleOnSphereFunctionBasis FunctionBasis(NbrOrbitals - 1);
      Complex* FunctionBasisEvaluation = new Complex [NbrOrbitals];
      double* FunctionBasisDecomposition = new double [NbrOrbitals];
      double* FunctionBasisDecompositionError  = new double [NbrOrbitals];
      for (int k = 0; k < NbrOrbitals; ++k)
	{
	  FunctionBasisDecomposition[k] = 0.0;
	  FunctionBasisDecompositionError[k] = 0.0;
	}

      double TotalProbability = 0.0;
      double TotalProbabilityError = 0.0;
      double ThetaStep = M_PI / NbrSteps;

      double Theta = 0.0;
      double* Density = new double [NbrSteps];
      for (int i = 0; i < NbrSteps; ++i)
	Density[i] = 0.0;

      RandomUV (TmpUV, TmpPositions, NbrParticles, RandomNumber);
      
      int NextCoordinates = 0;
      double PreviousProbabilities = 0.0;
      double CurrentProbabilities = 0.0;
      TotalProbability = 0.0;
      int Acceptance = 0;	  
      
      Complex Tmp = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
      CurrentProbabilities = SqrNorm(Tmp);
      PreviousProbabilities = CurrentProbabilities;
      
      cout << "starting warm-up sequence" << endl;
      for (int i = 1; i < NbrWarmUpIter; ++i)
	{
	  for (int j = 0; j < (NbrParticles * 2); ++j)
	    {
	      PreviousTmpUV[j] = TmpUV[j];
	      PreviousTmpPositions[j] = TmpPositions[j];
	    }
	  RandomUV(TmpUV, TmpPositions, NbrParticles, RandomNumber);
	  Complex TmpMetropolis = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
	  CurrentProbabilities = SqrNorm(TmpMetropolis);
	  if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	    {
	      PreviousProbabilities = CurrentProbabilities;
	      ++Acceptance;
	    }
	  else
	    {
	      for (int j = 0; j < (NbrParticles * 2); ++j)
		{
		  TmpUV[j] = PreviousTmpUV[j];
		  TmpPositions[j] = PreviousTmpPositions[j];
		}
	      CurrentProbabilities = PreviousProbabilities;
	    }
	} 
      cout << "acceptance rate = " <<  ((((double) Acceptance) / ((double) NbrWarmUpIter)) * 100.0) << "%" <<endl;
      
      for (int k = 0; k < NbrOrbitals; ++k)
	FunctionBasis.GetFunctionValue(TmpUV[0], TmpUV[1], FunctionBasisEvaluation[k], k);

      cout << "starting MC sequence" << endl;
      Acceptance = 0;

      for (int i = 0; i < NbrIter; ++i)
	{
	  for (int j = 0; j < (NbrParticles * 2); ++j)
	    {
	      PreviousTmpUV[j] = TmpUV[j];
	      PreviousTmpPositions[j] = TmpPositions[j];
	    }
	  RandomUV(TmpUV, TmpPositions, NbrParticles, RandomNumber);
	  Complex TmpMetropolis = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
	  CurrentProbabilities = SqrNorm(TmpMetropolis);
//  	       if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
//  		 {
//  		   PreviousProbabilities = CurrentProbabilities;
//  		   ++Acceptance;
// 		   for (int k = 0; k < NbrOrbitals; ++k)
// 		     FunctionBasis.GetFunctionValue(TmpUV[0], TmpUV[1], FunctionBasisEvaluation[k], k);
//  		 }
//  	       else
//  		 {
// 		   for (int j = 0; j < (NbrParticles * 2); ++j)
// 		     {
// 		       TmpUV[j] = PreviousTmpUV[j];
// 		       TmpPositions[j] = PreviousTmpPositions[j];
// 		     }
//  		   CurrentProbabilities = PreviousProbabilities;
//  		 }
	  TotalProbability += CurrentProbabilities;
	  TotalProbabilityError += CurrentProbabilities * CurrentProbabilities;
	  for (int k = 0; k < NbrOrbitals; ++k)
	    {
	      double TmpTruc1 = 0.0;
	      Complex TmpTruc2 = 0.0;
	      for (int j = 0; j < NbrParticles; ++j)
		{
		  FunctionBasis.GetFunctionValue(TmpUV[j << 1], TmpUV[(j << 1) + 1], TmpTruc2, k);
		  TmpTruc1 +=  SqrNorm(TmpTruc2);
		}
	      FunctionBasisDecomposition[k] += TmpTruc1 * CurrentProbabilities;
	      FunctionBasisDecompositionError[k] += (TmpTruc1 * CurrentProbabilities * TmpTruc1 * CurrentProbabilities);
	    }
	}
      cout << "acceptance rate = " <<  ((((double) Acceptance) / ((double) NbrIter)) * 100.0) << "%" <<endl;

      TotalProbabilityError /= (double) NbrIter;
      TotalProbabilityError -= (TotalProbability * TotalProbability) / (((double) NbrIter) * ((double) NbrIter));
      TotalProbabilityError = sqrt(TotalProbabilityError) / (TotalProbability / ((double) NbrIter));
      TotalProbabilityError /= sqrt ((double) NbrIter);
      for (int k = 0; k < NbrOrbitals; ++k)
	{
	  FunctionBasisDecompositionError[k] *= ((double) NbrIter);
	  FunctionBasisDecompositionError[k] -= FunctionBasisDecomposition[k] * FunctionBasisDecomposition[k];
	  FunctionBasisDecompositionError[k] = sqrt(FunctionBasisDecompositionError[k] / ((double) NbrIter)) / FunctionBasisDecomposition[k];
	  FunctionBasisDecompositionError[k] += TotalProbabilityError;
	}

      TotalProbability =  1.0 / TotalProbability;
      for (int k = 0; k < NbrOrbitals; ++k)
	{
	  FunctionBasisDecomposition[k] *= TotalProbability;
	  FunctionBasisDecompositionError[k] *= FunctionBasisDecomposition[k];
	}

      ofstream DensityRecordFile;
      DensityRecordFile.precision(14);
      DensityRecordFile.open(((SingleStringOption*) Manager["output"])->GetString(), ios::out);

      RealVector TmpPos(2, true);
      Theta = 0.0;
      double Sum = 0.0; 
      double TmpFactor = 4.0 * M_PI / ((double) NbrOrbitals);
      DensityRecordFile << "# decomposition of the wave on " << NbrOrbitals << " orbitals" << endl; 
      DensityRecordFile << "# orbital_index     component    error" << endl;
      for (int k = 0; k < NbrOrbitals; ++k)
	{
	  DensityRecordFile << "# " << k << " " << FunctionBasisDecomposition[k] << " " <<  FunctionBasisDecompositionError [k] << endl;
	  cout << k << " " << FunctionBasisDecomposition[k] << " " <<  FunctionBasisDecompositionError [k] << endl;
	}
      DensityRecordFile << "#" << endl << "# density wave function " << endl << "# theta density" << endl;
      for (int i = 0; i <= NbrSteps; ++i)
	{
	  TmpPos[0] = Theta;
	  double Tmp = 0.0;
	  Complex Tmp2;
	  for (int k = 0; k < NbrOrbitals; ++k)
	    {
	      FunctionBasis.GetFunctionValue(TmpPos, Tmp2, k);
	      //	      Tmp += SqrNorm(Tmp2);
	      Tmp += FunctionBasisDecomposition[k] * SqrNorm(Tmp2);
	    }
	  Tmp *= TmpFactor;
	  Sum += sin(Theta) * Tmp;
	  DensityRecordFile << Theta << " " << Tmp << endl;
	  Theta += ThetaStep;
	}
      DensityRecordFile.close();
      cout << "Sum = " << (Sum * ThetaStep * 2.0 * M_PI) << endl;
    }
   else
     {
       
       ComplexVector UV (NbrParticles * 2, true);
       RealVector TmpPositions (NbrParticles * 2, true);
       RandomUV (UV, TmpPositions, NbrParticles, RandomNumber);
       cout << SymmetrizedFunction->CalculateFromSpinorVariables(UV) << endl;;
       for (int i = 0; i < NbrParticles; ++i)
	 {
	   for (int j = i + 1; j < NbrParticles; ++j)
	     {
	       cout << i << "<->" << j << " : ";
	       cout << SymmetrizedFunction->CalculateFromSpinorVariables(UV) << " | " ;
	       FlipCoordinates(UV, i, j);
	       cout << SymmetrizedFunction->CalculateFromSpinorVariables(UV) << endl;
	       FlipCoordinates(UV, i, j);
	       
	     }
	 }
     }
  return 0;
}

void RandomUV (ComplexVector& uv, RealVector& positions, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  for (int j = 0; j < nbrParticles; ++j)
    {
      double x = acos (1.0 - (2.0 * randomNumberGenerator->GetRealRandomNumber()));
      double y = 2.0 * M_PI * randomNumberGenerator->GetRealRandomNumber();
      positions[2 * j] = x;
      positions[(2 * j) + 1] = y;
      uv.Re(2 * j) = cos(0.5 * x);
      uv.Im(2 * j) = uv.Re(2 * j) * sin(0.5 * y);
      uv.Re(2 * j) *= cos(0.5 * y);
      uv.Re(2 * j + 1) = sin(0.5 * x);
      uv.Im(2 * j + 1) = - uv.Re(2 * j + 1) * sin(0.5 * y);
      uv.Re(2 * j + 1) *= cos(0.5 * y);      
    }
}

void RandomUVOneCoordinate(ComplexVector& uv, RealVector& positions, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  coordinate *= 2;
  double x = acos (1.0 - (2.0 * randomNumberGenerator->GetRealRandomNumber()));
  double y = 2.0 * M_PI * randomNumberGenerator->GetRealRandomNumber();
  positions[coordinate] = x;
  uv.Re(coordinate) = cos(0.5 * x);
  uv.Im(coordinate) = uv.Re(coordinate) * sin(0.5 * y);
  uv.Re(coordinate) *= cos(0.5 * y);
  ++coordinate;
  positions[coordinate] = y;
  uv.Re(coordinate) = sin(0.5 * x);
  uv.Im(coordinate) = - uv.Re(coordinate) * sin(0.5 * y);
  uv.Re(coordinate) *= cos(0.5 * y);      
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


