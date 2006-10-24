#include "Vector/RealVector.h"

#include "Tools/FQHEWaveFunction/LaughlinOnSphereWaveFunction.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

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
  OptionManager Manager ("QHEFermionsOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MonteCarloGroup = new OptionGroup ("Monte Carlo options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += MonteCarloGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleStringOption  ('\n', "test-wavefunction", "name of the test wave fuction", "laughlin");
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  

  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "nbr-warmup", "number of Monte Carlo iterations that have to be done before evaluating the energy (i.e. warm up sequence)", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MonteCarloGroup) += new BooleanOption  ('\n', "show-details", "show intermediate values used for energy calculation", false);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionMonteCarloEnergy -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((BooleanOption*) Manager["list-wavefunctions"])->GetBoolean() == true)
    {
      return 0;
    }

  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  int NbrWarmUpIter = ((SingleIntegerOption*) Manager["nbr-warmup"])->GetInteger();

  Abstract1DComplexFunction* WaveFunction = new LaughlinOnSphereWaveFunction(NbrFermions, 3);
  RealVector Location(2 * NbrFermions, true);

  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (29457);

  double Energy = 0.0;
  double EnergyError = 0.0;
  double Normalization = 0.0;
  double NormalizationError = 0.0;
  double Tmp = 0.0;
  double Dist;

  Complex Tmp3;
  double Tmp2;
  double Tmp2bis;
  int NextCoordinates = 0;
  double Radius;
  for (int j = 0; j < NbrFermions; ++j)
    {
      Location[j << 1] = acos (1.0- (2.0 * RandomNumber->GetRealRandomNumber()));
      Location[1 + (j << 1)] = 2.0 * M_PI * RandomNumber->GetRealRandomNumber();
    }
  double PreviousProbabilities = SqrNorm((*WaveFunction)(Location));
  double CurrentProbabilities = PreviousProbabilities;
  double PreviousCoordinates1;
  double PreviousCoordinates2;
  double Theta;
  double Phi;
  long New = 0l;
  long Kept = 0l;
  if (NbrWarmUpIter > 0)
    cout << "starting warm-up sequence" << endl;
  for (int i = 1; i < NbrWarmUpIter; ++i)
    {      
      PreviousCoordinates1 = Location[NextCoordinates << 1];
      PreviousCoordinates2 = Location[1 + (NextCoordinates << 1)];
      Location[NextCoordinates << 1] = acos (1.0- (2.0 * RandomNumber->GetRealRandomNumber()));	  
      Location[1 + (NextCoordinates << 1)] = 2.0 * M_PI * RandomNumber->GetRealRandomNumber();
      CurrentProbabilities = SqrNorm((*WaveFunction)(Location));
//      cout << CurrentProbabilities << " " << PreviousProbabilities << endl;
      if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	{
	  PreviousProbabilities = CurrentProbabilities;
	  ++New;
	}
      else
 	{
 	  Location[NextCoordinates << 1] = PreviousCoordinates1;
 	  Location[1 + (NextCoordinates << 1)] = PreviousCoordinates2;
 	  CurrentProbabilities = PreviousProbabilities;
	  ++Kept;
 	}
      NextCoordinates = (int) (((double) NbrFermions) * RandomNumber->GetRealRandomNumber());
    }
  if (NbrWarmUpIter > 0)
    cout << "warm-up sequence is over" << endl;
  double MeanRejectRate = 0.7;
  for (int i = 1; i < NbrIter; ++i)
    {      
      PreviousCoordinates1 = Location[NextCoordinates << 1];
      PreviousCoordinates2 = Location[1 + (NextCoordinates << 1)];
      Location[NextCoordinates << 1] = acos (1.0- (2.0 * RandomNumber->GetRealRandomNumber()));	  
      Location[1 + (NextCoordinates << 1)] = 2.0 * M_PI * RandomNumber->GetRealRandomNumber();
      CurrentProbabilities = SqrNorm((*WaveFunction)(Location));
//      cout << CurrentProbabilities << " " << PreviousProbabilities << endl;
      double AcceptanceRate = RandomNumber->GetRealRandomNumber();
//      RejectRate = ((6.0 * MeanRejectRate - 3.0) * RejectRate + (4.0 - 6.0 * MeanRejectRate)) * RejectRate;
      if ((CurrentProbabilities > PreviousProbabilities) || ((AcceptanceRate * PreviousProbabilities) < CurrentProbabilities))
	{
	  PreviousProbabilities = CurrentProbabilities;
	  ++New;
	}
      else
 	{
 	  Location[NextCoordinates << 1] = PreviousCoordinates1;
 	  Location[1 + (NextCoordinates << 1)] = PreviousCoordinates2;
 	  CurrentProbabilities = PreviousProbabilities;
	  ++Kept;
 	}
      NextCoordinates = (int) (((double) NbrFermions) * RandomNumber->GetRealRandomNumber());
      if (NextCoordinates == NbrFermions)
	--NextCoordinates;
      
      Tmp = 0.0;
      for (int j = 0; j < NbrFermions; ++j)
	{
	  Theta = Location[j << 1];
	  Phi = Location[1 + (j << 1)];
	  for (int k = j + 1; k < NbrFermions; ++k)
	    {
	      Dist = sqrt ((sin (0.5 * (Theta - Location[k << 1])) * sin (0.5 * (Theta - Location[k << 1])))
			   + (sin (Theta) * sin (Location[k << 1]) * sin (0.5 * (Phi - Location[1 + (k << 1)])) * sin (0.5 * (Phi - Location[1 + (k << 1)]))));	      
	      if (Dist > 1e-6)
		{
		  Dist = 0.5 / (Dist * sqrt(0.5 * ((double) LzMax)));
		  Tmp += Dist;
		}
	    }
	}
      Energy += Tmp;
      EnergyError += Tmp * Tmp;
      if ((i % (((SingleIntegerOption*) Manager["display-step"])->GetInteger())) == 0)
 	{
	  cout << " i = " << i << endl;
	  double TmpEnergy = Energy / ((double) i);
	  double TmpEnergyError = sqrt (((EnergyError  / ((double) i)) - (TmpEnergy * TmpEnergy)) / ((double) i));
	  cout << ((TmpEnergy - ((0.5 * ((double) (NbrFermions * NbrFermions))) / sqrt(0.5 * ((double) LzMax)))) / ((double) NbrFermions)) << " +/- " << (TmpEnergyError / ((double) NbrFermions)) << endl;
 	  cout << "raw : " << Energy << " " << (TmpEnergy * TmpEnergy) << " " << (EnergyError  / ((double) i)) << " " << TmpEnergy << " +/- " << TmpEnergyError << endl;
	  cout << New << " " << Kept << endl;
	  New = 0l;
	  Kept = 0l;
 	  cout << "-----------------------------------------------" << endl;
 	}
    } 
  cout << " final results :" << endl;
  double TmpEnergy = Energy / ((double) NbrIter);
  double TmpEnergyError = sqrt (((EnergyError  / ((double) NbrIter)) - (TmpEnergy * TmpEnergy)) / ((double) NbrIter));
  cout << TmpEnergy << " +/- " << TmpEnergyError << endl;
  cout << ((TmpEnergy - ((0.5 * ((double) (NbrFermions * NbrFermions))) / sqrt(0.5 * ((double) LzMax)))) / ((double) NbrFermions)) << " +/- " << (TmpEnergyError / ((double) NbrFermions)) << endl;
  cout << "-----------------------------------------------" << endl;
 return 0;
}


