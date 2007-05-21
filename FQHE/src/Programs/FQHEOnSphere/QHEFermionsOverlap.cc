#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "Tools/FQHEWaveFunction/QHEWaveFunctionManager.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "MCObservables/RealObservable.h"
#include "MCObservables/ComplexObservable.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>

double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;

Complex OverlapValue(ComplexObservable &ScalarProduct, RealObservable &Norm);
double OverlapError(ComplexObservable &ScalarProduct, RealObservable &Norm);

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("QHEFermionsOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MonteCarloGroup = new OptionGroup ("Monte Carlo options");

  ArchitectureManager Architecture;
  QHEWaveFunctionManager WaveFunctionManager;

  Manager += SystemGroup;
  WaveFunctionManager.AddOptionGroup(&Manager);
  Manager += MonteCarloGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 7);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 12);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "lz", "twice the momentum projection", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "exact-state", "name of the file containing the vector obtained using exact diagonalization");
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  
  (*SystemGroup) += new SingleStringOption  ('\n', "use-exact", "file name of an exact state that has to be used as test wave function");

  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "record-step", "number of iteration between two consecutive result recording of energy value (0 if no on-disk recording is needed)", 0);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "record-file", "name of the file where energy recording has to be done", "montecarlo.dat");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "with-timecoherence", "use time coherence between two successive evaluation of the wave function");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "show-details", "show intermediate values used for overlap calculation", false);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionOverlap -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetBoolean("list-wavefunctions") == true)
    {
      WaveFunctionManager.ShowAvalaibleWaveFunctions(cout);
      return 0;
    }

  int NbrFermions = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int NbrIter = Manager.GetInteger("nbr-iter");
  int Lz = Manager.GetInteger("lz");

  if (Manager.GetString("exact-state") == 0)
    {
      cout << "QHEFermionOverlap requires an exact state" << endl;
      return -1;
    }
  RealVector State;
  if (State.ReadVector (Manager.GetString("exact-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("exact-state") << endl;
      return -1;      
    }
  if (Manager.GetString("use-exact") != 0)
    {
      RealVector TestState;
      if (TestState.ReadVector (Manager.GetString("use-exact")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("use-exact") << endl;
	  return -1;      
	}
      if (State.GetVectorDimension() != TestState.GetVectorDimension())
	{
	  cout << "dimension mismatch" << endl;
	  return -1;      
	}
      cout << "overlap = " << (TestState * State) << endl;
      return 0;
    }

  Abstract1DComplexFunction* TestWaveFunction = WaveFunctionManager.GetWaveFunction();
  
  if (TestWaveFunction == 0)
    {
      cout << "no or unknown analytical wave function" << endl;
      return -1;
    }
  FermionOnSphere Space (NbrFermions, Lz, LzMax);
  ParticleOnSphereFunctionBasis Basis(LzMax);
  
  RealVector Location(2 * NbrFermions, true);

  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (29457);

  
  int RecordStep = Manager.GetInteger("record-step");

  Complex* RecordedOverlap = 0;
  double* RecordedOverlapError = 0;
  if (RecordStep != 0)
    {
      RecordedOverlap = new Complex [(NbrIter / RecordStep) + 1];
      RecordedOverlapError = new double [(NbrIter / RecordStep) + 1];
    }
  int RecordIndex = 0;
  double Factor = 1.0;
  for (int j = 0; j < NbrFermions; ++j)
    {
      Factor *= 4.0 * M_PI;
    }
  Complex Overlap(0.0);
  ComplexObservable ScalarProduct(100);    
  RealObservable NormExact(100);
  Complex TrialValue;
  Complex ValueExact;
  Complex TmpMetropolis;
  int NextCoordinates = 0;
  // initialize with random coordinate positions:
  for (int j = 0; j < NbrFermions; ++j)
    {
      Location[j << 1] = acos (1.0- (2.0 * RandomNumber->GetRealRandomNumber()));
      Location[1 + (j << 1)] = 2.0 * M_PI * RandomNumber->GetRealRandomNumber();
    }
  TrialValue = (*TestWaveFunction)(Location);
  double PreviousSamplingAmplitude = SqrNorm(TrialValue);
  double CurrentSamplingAmplitude = PreviousSamplingAmplitude;

  for (int i = 0; i < NbrIter; ++i)
    {
      // store old position
      double PreviousCoordinates1 = Location[NextCoordinates << 1];
      double PreviousCoordinates2 = Location[1 + (NextCoordinates << 1)];
      // make a random move
      Location[NextCoordinates << 1] = acos (1.0- (2.0 * RandomNumber->GetRealRandomNumber()));	  
      Location[1 + (NextCoordinates << 1)] = 2.0 * M_PI * RandomNumber->GetRealRandomNumber();
      TmpMetropolis = (*TestWaveFunction)(Location);      
      // Complex TmpMetropolis2 = (*TestWaveFunction2)(Location);
      //       cout << TmpMetropolis << " " << (TmpMetropolis2) << endl;
      // accept or reject move according to probability |Psi_new|^2  / |Psi_old|^2
      CurrentSamplingAmplitude = SqrNorm(TmpMetropolis);
      if ((CurrentSamplingAmplitude > PreviousSamplingAmplitude) ||
	  ((RandomNumber->GetRealRandomNumber() * PreviousSamplingAmplitude) < CurrentSamplingAmplitude))
	{
	  PreviousSamplingAmplitude = CurrentSamplingAmplitude;
	  TrialValue = TmpMetropolis;
	}
      else
	{
	  Location[NextCoordinates << 1] = PreviousCoordinates1;
	  Location[1 + (NextCoordinates << 1)] = PreviousCoordinates2;
	  CurrentSamplingAmplitude = PreviousSamplingAmplitude;
	}
      // determine next particle to move
      NextCoordinates = (int) (((double) NbrFermions) * RandomNumber->GetRealRandomNumber());
      if (NextCoordinates == NbrFermions) --NextCoordinates;
      // calculate value of the exact wavefunction
      int TimeCoherence = NextCoordinates;
      if (Manager.GetBoolean("with-timecoherence") == false) TimeCoherence = -1;
      QHEParticleWaveFunctionOperation Operation(&Space, &State, &Location, &Basis, TimeCoherence);
      Operation.ApplyOperation(Architecture.GetArchitecture());      
      ValueExact = Operation.GetScalar() ;
      // note observations:
      double norm = SqrNorm(ValueExact)/CurrentSamplingAmplitude;
      NormExact << norm;
      Overlap = Conj(TrialValue)*ValueExact/CurrentSamplingAmplitude;
      ScalarProduct << Overlap;
      if ((i > 0) && ((RecordStep != 0) && ((i % RecordStep) == 0)))
	{
	  RecordedOverlap[RecordIndex] = OverlapValue(ScalarProduct, NormExact);
	  RecordedOverlapError[RecordIndex] = OverlapError(ScalarProduct, NormExact);
	  ++RecordIndex;
	}
      if ((i > 0) && ((i % (Manager.GetInteger("display-step"))) == 0))
	{
	  cout << " i = " << i << endl;
	  cout << OverlapValue(ScalarProduct, NormExact) << " +/- " << OverlapError(ScalarProduct, NormExact) << endl;
	  cout << "-----------------------------------------------" << endl;
	}
    } 
  cout << " final results :" << endl;
  cout << OverlapValue(ScalarProduct, NormExact) << " +/- " << OverlapError(ScalarProduct, NormExact) << endl;
  cout << "-----------------------------------------------" << endl;

  if (RecordStep != 0)
    {
      RecordedOverlap[RecordIndex] = OverlapValue(ScalarProduct, NormExact);
      RecordedOverlapError[RecordIndex] = OverlapError(ScalarProduct, NormExact);
      ofstream OverlapRecordFile;
      OverlapRecordFile.precision(14);
      OverlapRecordFile.open(Manager.GetString("record-file"), ios::out | ios::binary);
      int NbrRecords = NbrIter / RecordStep;
      OverlapRecordFile << "# Monte Carlo overlap calculation" << endl
		       << "# step overlap.Re overlap.Im overlap2 error error2" << endl;
      for (int i = 0; i < NbrRecords; ++i)
	OverlapRecordFile << i << " " << RecordedOverlap[i].Re << " " << RecordedOverlap[i].Im << " " << SqrNorm(RecordedOverlap[i])
			  << " " << RecordedOverlapError[i] << " " <<  DSQR(RecordedOverlapError[i]) << endl;
      OverlapRecordFile.close();
    }

  return 0;
}


Complex OverlapValue(ComplexObservable &ScalarProduct, RealObservable &NormObs)
{
  return Complex(ScalarProduct.Average()/sqrt(NormObs.Average()));
}
  
double OverlapError(ComplexObservable &ScalarProduct, RealObservable &NormObs)
{
  double norm = NormObs.Average();
  double prod = Norm(ScalarProduct.Average());
  return sqrt( DSQR(ScalarProduct.ErrorEstimate())/norm + DSQR(prod*NormObs.ErrorEstimate()/2.0/norm)/norm);
}
