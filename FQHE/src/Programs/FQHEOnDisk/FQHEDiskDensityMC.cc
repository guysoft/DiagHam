#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "HilbertSpace/BosonOnSphere.h"
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

#include "Tools/FQHEWaveFunction/PfaffianOnDiskWaveFunction.h"
#include "Tools/FQHEWaveFunction/LaughlinOnDiskWaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHEDiskLaughlinOneQuasiholeWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnDiskTwoQuasiholeWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnDiskTwoQuasielectronWaveFunction.h"

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
using std::flush;
using std::ofstream;


void RandomZ (RealVector& positions, double scale, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator);

void RandomZOneCoordinate(RealVector& positions, double scale, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator);

bool RandomZOneCoordinateWithJump(RealVector& positions, double scale, double jump, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator);

// flip two one-body coordinates
//
// coordinates = reference on the n-body coordinate vector
// i = index of the first one body coordinate
// j = index of the second one body coordinate
void FlipCoordinates (RealVector& coordinates, int i, int j);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHEDiskDensityMC" , "0.01");
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
  (*SystemGroup) += new BooleanOption ('\n', "laughlin", "do the calculation for the Laughlin state instead of the Moore-Read state");  
  (*SystemGroup) += new BooleanOption ('\n', "test-symmetry", "check the test wave function is symmetric/antisymmetric");  
  (*SystemGroup) += new SingleStringOption  ('\n', "load-permutations", "read all the permutations needed to compute the reference wave function from a file");  
  (*SystemGroup) += new SingleStringOption  ('\n', "save-permutations", "file name where all the permutations needed to compute the reference wave function have to be stored");  
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistic");
  (*SystemGroup) += new BooleanOption ('\n', "quasielectron", "plot quasiparticles instead of quasiholes");  
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output file name", "density.dat"); 

  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "nbr-orbitals", "number of orbitals used to sample the disk geometry", 100);
  (*MonteCarloGroup) += new SingleIntegerOption  ('s', "nbr-step", "number of steps for the density profil in each direction", 20);
  (*MonteCarloGroup) += new SingleDoubleOption  ('\n', "grid-length", "size of the sampling grid (in unit of the magnetic length", 20);
  (*MonteCarloGroup) += new SingleDoubleOption  ('\n', "jump", "length of the jump used for the metropolis algorithm", 0.3);
  
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
  (*MonteCarloGroup) += new  SingleStringOption ('\n', "record-wavefunctions", "optional file where each wavefunctions will be tabulated and recorded");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskDensityMC -h" << endl;
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
  bool LaughlinFlag = (((BooleanOption*) Manager["laughlin"])->GetBoolean());
  double StepSize = 0.25;
  char* RecordWaveFunctions = ((SingleStringOption*) Manager["record-wavefunctions"])->GetString();
  if (RecordWaveFunctions != 0)
    {
      ofstream RecordFile;
      RecordFile.open(RecordWaveFunctions, ios::out | ios::binary);
      RecordFile.close();
    }
  double GridScale =  ((SingleDoubleOption*) Manager["grid-length"])->GetDouble();
  double InvGridStep = ((double) NbrSteps) / GridScale;
  double GridStep = 1.0 / InvGridStep;
  double GridShift = -0.5 * GridScale;
  double MCJump = ((SingleDoubleOption*) Manager["jump"])->GetDouble();

  Abstract1DComplexFunction* SymmetrizedFunction = 0;
  if (QuasielectronFlag == true)
    {
//      if (LaughlinFlag == true)
//	SymmetrizedFunction = new FQHEDiskLaughlinOneQuasielectronWaveFunction(NbrParticles, 0.0, 0.0, StatisticFlag);
//      else
// 	{
// 	  if (((SingleStringOption*) Manager["load-permutations"])->GetString() == 0)
// 	    SymmetrizedFunction = new PfaffianOnDiskTwoQuasielectronWaveFunction(NbrParticles, 0.0, 0.0, M_PI, 0.0, StatisticFlag);
// 	  else
// 	    SymmetrizedFunction = new PfaffianOnDiskTwoQuasielectronWaveFunction(((SingleStringOption*) Manager["load-permutations"])->GetString(),
// 										   0.0, 0.0, M_PI, 0.0, StatisticFlag);
// 	  if (((SingleStringOption*) Manager["save-permutations"])->GetString() != 0)
// 	    {
// 	      ((PfaffianOnDiskTwoQuasielectronWaveFunction*) SymmetrizedFunction)->WritePermutations(((SingleStringOption*) Manager["save-permutations"])->GetString());
// 	      return 0;
// 	    }
// 	}
    }
  else
    {
      if (LaughlinFlag == true)
	if (StatisticFlag == true)
	  SymmetrizedFunction = new FQHEDiskLaughlinOneQuasiholeWaveFunction(NbrParticles, 0.0, 3);
	else
	  SymmetrizedFunction = new FQHEDiskLaughlinOneQuasiholeWaveFunction(NbrParticles, 0.0, 2);
      else
	SymmetrizedFunction = new PfaffianOnDiskTwoQuasiholeWaveFunction(NbrParticles, -0.25 * GridScale, 0.25 * GridScale, StatisticFlag);
    }

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
      int TwiceNbrParticles = NbrParticles * 2;
      RealVector TmpZ (TwiceNbrParticles, true);
      double PreviousTmpZRe;
      double PreviousTmpZIm;
      double** FunctionBasisDecomposition = new double* [NbrSteps];
      for (int k = 0; k < NbrSteps; ++k)
	{
	  FunctionBasisDecomposition[k] = new double [NbrSteps];
	  for (int j = 0; j < NbrSteps; ++j)
	    FunctionBasisDecomposition[k][j] = 0.0;
	}
      double* FunctionBasisDecompositionError  = new double [NbrOrbitals];
      double* ExpTable = new double [NbrParticles];
      double PreviousExp = 0.0;

      int* GridLocations = new int [TwiceNbrParticles];

      double TotalProbability = 0.0;
      double TotalProbabilityError = 0.0;
      int CurrentPercent = 0;

      RandomZ (TmpZ, GridScale, NbrParticles, RandomNumber);
      int NextCoordinate = 0;
      double PreviousProbabilities = 0.0;
      double CurrentProbabilities = 0.0;
      TotalProbability = 0.0;
      int Acceptance = 0;	  
      
      Complex Tmp = (*SymmetrizedFunction)(TmpZ);
      CurrentProbabilities = SqrNorm(Tmp);
      for (int k = 0; k < NbrParticles; ++k)
	{
	  ExpTable[k] = exp(-0.5 * ((TmpZ[(k << 1)] * TmpZ[(k << 1)])
				    + (TmpZ[(k << 1) + 1] * TmpZ[(k << 1) + 1])));
	  CurrentProbabilities *= ExpTable[k];
	}
      PreviousProbabilities = CurrentProbabilities;
      
      cout << "starting warm-up sequence" << endl;
      for (int i = 1; i < NbrWarmUpIter; ++i)
	{
	  PreviousTmpZRe = TmpZ[NextCoordinate << 1];
	  PreviousTmpZIm = TmpZ[(NextCoordinate << 1) + 1];
	  if (RandomZOneCoordinateWithJump(TmpZ, GridScale, MCJump, NextCoordinate, RandomNumber))
	    {
	      PreviousExp = ExpTable[NextCoordinate];
	      Complex TmpMetropolis = (*SymmetrizedFunction)(TmpZ);
	      CurrentProbabilities = SqrNorm(TmpMetropolis);
 	      ExpTable[NextCoordinate] = exp(-0.5 * ((TmpZ[(NextCoordinate << 1)] * TmpZ[(NextCoordinate << 1)])
 						     + (TmpZ[(NextCoordinate << 1) + 1] * TmpZ[(NextCoordinate << 1) + 1])));
//	      ExpTable[NextCoordinate] = exp(-0.5 * (((TmpZ[(NextCoordinate << 1)] * TmpZ[(NextCoordinate << 1)]) - (PreviousTmpZRe * PreviousTmpZRe))
//						     + ((TmpZ[(NextCoordinate << 1) + 1] * TmpZ[(NextCoordinate << 1) + 1])) - (PreviousTmpZIm * PreviousTmpZIm)));
//	      cout << ExpTable[NextCoordinate] << endl;
 	      for (int k = 0; k < NbrParticles; ++k)
 		CurrentProbabilities *= ExpTable[k];
//	      CurrentProbabilities *= ExpTable[NextCoordinate];
	      if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
		{
		  PreviousProbabilities = CurrentProbabilities;
		  ++Acceptance;
		}
	      else
		{
		  TmpZ[NextCoordinate << 1] = PreviousTmpZRe;
		  TmpZ[(NextCoordinate << 1) + 1] = PreviousTmpZIm;
		  ExpTable[NextCoordinate] = PreviousExp;
		  CurrentProbabilities = PreviousProbabilities;
		}
	    }
	  else
	    {
	      TmpZ[NextCoordinate << 1] = PreviousTmpZRe;
	      TmpZ[(NextCoordinate << 1) + 1] = PreviousTmpZIm;
	    }
	  NextCoordinate = (int) (RandomNumber->GetRealRandomNumber() * (double) NbrParticles);
	  if (((i * 20) / NbrWarmUpIter) != CurrentPercent)
	    {
	      CurrentPercent = (i * 20) / NbrWarmUpIter;
	      cout << (CurrentPercent * 5) << "% " << flush;
	    }
	} 
      cout << endl << "acceptance rate = " <<  ((((double) Acceptance) / ((double) NbrWarmUpIter)) * 100.0) << "%" <<endl;
      
      CurrentPercent = 0;
      for (int i = 0; i < TwiceNbrParticles; ++i)
	GridLocations[i] = (int) ((TmpZ[i] - GridShift) * InvGridStep);

      cout << "starting MC sequence" << endl;
      Acceptance = 0;

      for (int i = 0; i < NbrIter; ++i)
	{
	  PreviousTmpZRe = TmpZ[NextCoordinate << 1];
	  PreviousTmpZIm = TmpZ[(NextCoordinate << 1) + 1];
	  if (RandomZOneCoordinateWithJump(TmpZ, GridScale, MCJump, NextCoordinate, RandomNumber))
	    {
	      PreviousExp = ExpTable[NextCoordinate];
	      Complex TmpMetropolis = (*SymmetrizedFunction)(TmpZ);
	      CurrentProbabilities = SqrNorm(TmpMetropolis);
// 	      ExpTable[NextCoordinate] = exp(-0.5 * (((TmpZ[(NextCoordinate << 1)] * TmpZ[(NextCoordinate << 1)]) - (PreviousTmpZRe * PreviousTmpZRe))
// 						     + ((TmpZ[(NextCoordinate << 1) + 1] * TmpZ[(NextCoordinate << 1) + 1])) - (PreviousTmpZIm * PreviousTmpZIm)));
	      ExpTable[NextCoordinate] = exp(-0.5 * ((TmpZ[(NextCoordinate << 1)] * TmpZ[(NextCoordinate << 1)])
						     + (TmpZ[(NextCoordinate << 1) + 1] * TmpZ[(NextCoordinate << 1) + 1])));
	      for (int k = 0; k < NbrParticles; ++k)
		CurrentProbabilities *= ExpTable[k];
//	      CurrentProbabilities *= ExpTable[NextCoordinate];
	      if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
		{
		  PreviousProbabilities = CurrentProbabilities;
		  GridLocations[NextCoordinate << 1] = (int) ((TmpZ[NextCoordinate << 1] - GridShift) * InvGridStep);
		  if ((GridLocations[NextCoordinate << 1] < 0) || (GridLocations[NextCoordinate << 1] >= NbrSteps))
		    {
		      cout << "1 : " << GridLocations[NextCoordinate << 1] << " " << TmpZ[NextCoordinate << 1] << " " <<  GridShift << " " << InvGridStep << endl;
		    }
		  GridLocations[(NextCoordinate << 1) + 1] = (int) ((TmpZ[(NextCoordinate << 1) + 1] - GridShift) * InvGridStep);
		  if ((GridLocations[(NextCoordinate << 1) + 1] < 0) || (GridLocations[(NextCoordinate << 1) + 1] >= NbrSteps))
		    {
		      cout << "2 : " << GridLocations[(NextCoordinate << 1) + 1] << " " << TmpZ[(NextCoordinate << 1) + 1] << " " <<  GridShift << " " << InvGridStep << endl;
		    }
		  ++Acceptance;	      
		}
	      else
		{
		  TmpZ[NextCoordinate << 1] = PreviousTmpZRe;
		  TmpZ[(NextCoordinate << 1) + 1] = PreviousTmpZIm;
		  ExpTable[NextCoordinate] = PreviousExp;
		  CurrentProbabilities = PreviousProbabilities;
		}
	    }
	  else
	    {
	      TmpZ[NextCoordinate << 1] = PreviousTmpZRe;
	      TmpZ[(NextCoordinate << 1) + 1] = PreviousTmpZIm;
	    }

	  for (int j = 0; j < TwiceNbrParticles; j += 2)
	    FunctionBasisDecomposition[GridLocations[j]][GridLocations[j + 1]] += 1.0;

	  TotalProbability++;	  
	  TotalProbabilityError++;

	  NextCoordinate = (int) (RandomNumber->GetRealRandomNumber() * (double) NbrParticles);

	  if (((i * 20) / NbrIter) != CurrentPercent)
	    {
	      CurrentPercent = (i * 20) / NbrIter;
	      cout << (CurrentPercent * 5) << "% " << flush;
	    }
	}
      cout << endl << "acceptance rate = " <<  ((((double) Acceptance) / ((double) NbrIter)) * 100.0) << "%" <<endl;

      TotalProbabilityError /= (double) NbrIter;
      TotalProbabilityError -= (TotalProbability * TotalProbability) / (((double) NbrIter) * ((double) NbrIter));
      TotalProbabilityError = sqrt(TotalProbabilityError) / (TotalProbability / ((double) NbrIter));
      TotalProbabilityError /= sqrt ((double) NbrIter);

      TotalProbability =  1.0 / TotalProbability;
      for (int i = 0; i < NbrSteps; ++i)
	for (int j = 0; j < NbrSteps; ++j)	  
	  FunctionBasisDecomposition[i][j] *= TotalProbability;

      ofstream DensityRecordFile;
      DensityRecordFile.precision(14);
      DensityRecordFile.open(((SingleStringOption*) Manager["output"])->GetString(), ios::out);

      DensityRecordFile << "#" << endl << "# density wave function " << endl << "# x y  density density_error" << endl;
      double TmpX = GridShift;
      double Sum = 0.0;
      for (int i = 0; i < NbrSteps; ++i)
	{
	  double TmpY = GridShift;
	  for (int j = 0; j < NbrSteps; ++j)
	    {
	      DensityRecordFile << TmpX << " " << TmpY << " " << FunctionBasisDecomposition[i][j] << endl;
	      Sum += FunctionBasisDecomposition[i][j];
	      TmpY += GridStep;
	    }
	  DensityRecordFile << endl;
	  TmpX += GridStep;
	}
      DensityRecordFile.close();
      cout << "Sum = " << Sum << endl;
    }
   else
     {
       RealVector TmpZ (2 * NbrParticles, true);
       RandomZ (TmpZ, GridScale, NbrParticles, RandomNumber);
       cout << (*SymmetrizedFunction)(TmpZ) << endl;;
       for (int i = 0; i < NbrParticles; ++i)
	 {
	   for (int j = i + 1; j < NbrParticles; ++j)
	     {
	       cout << i << "<->" << j << " : ";
	       cout << (*SymmetrizedFunction)(TmpZ) << " | " ;
	       FlipCoordinates(TmpZ, i, j);
	       cout << (*SymmetrizedFunction)(TmpZ) << endl;
	       FlipCoordinates(TmpZ, i, j);	       
	     }
	 }
     }
  return 0;
}

void RandomZ (RealVector& positions, double scale, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  for (int j = 0; j < nbrParticles; ++j)
    {
      positions[2 * j] = scale * (0.5 - randomNumberGenerator->GetRealRandomNumber());
      positions[(2 * j) + 1] = scale * (0.5 - randomNumberGenerator->GetRealRandomNumber());
    }
}

void RandomZOneCoordinate(RealVector& positions, double scale, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  coordinate *= 2;
  positions[coordinate] = scale * (0.5 - randomNumberGenerator->GetRealRandomNumber());
  ++coordinate;
  positions[coordinate] = scale * (0.5 - randomNumberGenerator->GetRealRandomNumber());
}

bool RandomZOneCoordinateWithJump(RealVector& positions, double scale, double jump, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  coordinate *= 2;
  double OldCoordinate =  positions[coordinate];
  double TmpJump = scale * jump * (0.5 - randomNumberGenerator->GetRealRandomNumber());
  if (((OldCoordinate + TmpJump) > (0.5 * scale)) || ((OldCoordinate + TmpJump) < -(0.5 * scale)))
    return false;
  positions[coordinate] += TmpJump;
  ++coordinate;
  OldCoordinate =  positions[coordinate];
  TmpJump = scale * jump * (0.5 - randomNumberGenerator->GetRealRandomNumber()); 
  if (((OldCoordinate + TmpJump) > (0.5 * scale)) || ((OldCoordinate + TmpJump) < -(0.5 * scale)))
    return false;
  positions[coordinate] += TmpJump;
  return true;
}

// flip two one-body coordinates
//
// coordinates = reference on the n-body coordinate vector
// i = index of the first one body coordinate
// j = index of the second one body coordinate

void FlipCoordinates (RealVector& coordinates, int i, int j)
{
  double Tmp = coordinates[i << 1];  
  coordinates[i << 1] = coordinates[j << 1];
  coordinates[j << 1] = Tmp;
  Tmp = coordinates[(i << 1) + 1];  
  coordinates[(i << 1) + 1] = coordinates[(j << 1) + 1];
  coordinates[(j << 1) + 1] = Tmp;
}


// get locations on the grid from a position vector
//
// coordinates = reference on the n-body coordinate vector
// locations = array of integer location on the grid
// invGridStep = sinvert of the grid cell size
// gridShift = X (or Y) coordinate of the grid upper left corner

void GetGridCoordinates (RealVector& coordinates, int* locations, double invGridStep, double gridShift)
{
  for (int i = 0; i < coordinates.GetVectorDimension(); ++i)
    locations[i] = (int) ((coordinates[i] - gridShift) * invGridStep);
}
