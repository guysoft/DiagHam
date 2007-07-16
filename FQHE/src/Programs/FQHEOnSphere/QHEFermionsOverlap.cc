#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "Tools/FQHEWaveFunction/QHEWaveFunctionManager.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexFunction.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"
//#include "Tools/FQHEWaveFunction/PairedCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/ParticleOnSphereCollection.h"
#include "Tools/FQHEWaveFunction/WaveFunctionOverlapOptimizer.h"
#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "MCObservables/RealObservable.h"
#include "MCObservables/ComplexObservable.h"
#include "MCObservables/WeightedRealObservable.h"
#include "MCObservables/WeightedComplexObservable.h"
#include "MCObservables/MCHistoryRecord.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/Endian.h"

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

Complex OverlapValue(WeightedComplexObservable &ScalarProduct, WeightedRealObservable &Norm1, WeightedRealObservable &Norm2);
double OverlapError(WeightedComplexObservable &ScalarProduct, WeightedRealObservable &Norm1, WeightedRealObservable &Norm2);

double MinDist (RealVector &Positions)
{
  int NbrFermions = Positions.GetVectorDimension()/2;
  double sp, spmax=-1.0;
  for ( int i=0; i<NbrFermions; ++i)
    for ( int j=i+1; j<NbrFermions; ++j)
      {
	sp=sin(Positions[2*i])*sin(Positions[2*j])*cos(Positions[2*i+1]-Positions[2*j+1]) +
	  cos(Positions[2*i])*cos(Positions[2*j]);
	if (sp>spmax) spmax=sp;
      }
  return acos(spmax);
}

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
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "randomSeed", "number of iteration between two consecutive result displays", -1);
  (*MonteCarloGroup) += new SingleIntegerOption ('H', "history-mode", "use on-file history: (0=off, 1=generate new, 2=read history, 3=optimize with history, 4=continue to generate given history)", 1);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "history-file", "name of the file where overlap recording has to be done", NULL);
  (*MonteCarloGroup) += new BooleanOption ('\n', "varyMR", "vary coefficient of 1/z in pair wavefunction");  
  (*MonteCarloGroup) += new SingleIntegerOption ('d', "sample-density", "spacing of samples to be saved in History-mode", 1);
  (*MonteCarloGroup) += new SingleIntegerOption ('\n', "maxPoints", "number of simultaneous function evaltuation in optimising History mode ", 50);
  (*MonteCarloGroup) += new SingleIntegerOption ('\n', "limitSamples", "maximal number of samples to be used from history-file", 10000000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "record-step", "number of iterations between two consecutive result recording the overlap value (0 if no on-disk recording is needed)", 0);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "record-file", "name of the file where overlap recording has to be done", "montecarlo.dat");
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
  RealVector State;

  if (Manager.GetInteger("history-mode")<2)
    {
      if (Manager.GetString("exact-state") == 0)
	{
	  cout << "QHEFermionOverlap requires an exact state" << endl;
	  return -1;
	}  
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
	  Complex ov = (TestState * State);
	  cout << "overlap = " << ov << endl;
	  cout << "Sqrovrl = " << SqrNorm(ov) << endl;
	  return 0;
	}
    }

  if (Manager.GetInteger("history-mode")==2)
    {
      if (Manager.GetString("exact-state") != 0)
	if (State.ReadVector (Manager.GetString("exact-state")))
	  {
	    cout << "Using exact state to check eventual outliers" << endl;
	  }
	else
	  {
	    cout << "can't open vector file " << Manager.GetString("exact-state") << endl;
	    return -1;      
	  }
      
    }
  
  Abstract1DComplexFunction* TestWaveFunction = WaveFunctionManager.GetWaveFunction();
  
  if (TestWaveFunction == 0)
    {
      cout << "no or unknown analytical wave function" << endl;
      return -1;
    }

  ParticleOnSphereCollection * Particles = new ParticleOnSphereCollection(NbrFermions, Manager.GetInteger("randomSeed"));
  Complex ValueExact;
  Complex TrialValue;
  double CurrentSamplingAmplitude;
  
  int HistoryMode = Manager.GetInteger("history-mode");
  int SampleDensity = Manager.GetInteger("sample-density");
  MCHistoryRecord *History=NULL;
  char *HistoryFileName=NULL;
  if (HistoryMode>0)
    {
      HistoryFileName=Manager.GetString("history-file");      
      if (HistoryMode==1)
	{
	  if (HistoryFileName==NULL)
	    {
	      // default filename: add extension to exact vector
	      HistoryFileName = new char[strlen(Manager.GetString("exact-state"))+6];
	      sprintf(HistoryFileName,"%s.samp",Manager.GetString("exact-state"));
	    }
	  char *tmpC = WaveFunctionManager.GetDescription();
	  History=new MCHistoryRecord(NbrIter, 2*NbrFermions, Manager.GetString("exact-state"), tmpC, HistoryFileName
				      /* could add additional observables here */);
	  delete [] tmpC;
	}
      else if ((HistoryFileName==NULL)&&(HistoryMode>1))
	{
	  cout << "History mode "<<HistoryMode<<" requires a history file!" << endl;
	  return -1;
	}
    }
  
  if (HistoryMode ==2)
    {
      // code to process available samples
      int sampleCount;
      int totalSampleCount=0;
      RealVector Positions(2*NbrFermions);
      WeightedRealObservable NormTrialObs(100);
      WeightedRealObservable NormExactObs(100);
      WeightedComplexObservable OverlapObs(100);
      History = new MCHistoryRecord(HistoryFileName, 2*NbrFermions /* could add additional observables here */);
      double typicalSA=0.0, typicalWF=0.0, typicalTV=0.0;
      int averageTypical=50;
      int initialSkip=10;
      for (int i=0; i<initialSkip; ++i)
	History->GetMonteCarloStep(sampleCount, CurrentSamplingAmplitude, &(Positions[0]), ValueExact);
      for (int i=0; i<averageTypical; ++i)
	{
	  History->GetMonteCarloStep(sampleCount, CurrentSamplingAmplitude, &(Positions[0]), ValueExact);
	  typicalSA+=CurrentSamplingAmplitude;
	  typicalWF+=Norm(ValueExact);
	  typicalTV+=Norm((*TestWaveFunction)(Positions));
	}
      typicalSA/=averageTypical;
      typicalWF/=averageTypical;
      typicalTV/=averageTypical;
      History->RewindHistory();
      cout << "typicalSA= " << typicalSA << ", typicalWF="<<typicalWF<<", typicalTV="<<typicalTV<<endl;
      int i=0;
      Complex rawExact;
      while ( (i++<NbrIter) && (History->GetMonteCarloStep(sampleCount, CurrentSamplingAmplitude, &(Positions[0]), ValueExact)))
	{
	  rawExact=ValueExact;
	  totalSampleCount+=sampleCount;
	  TrialValue = (*TestWaveFunction)(Positions)/typicalTV;	  
	  CurrentSamplingAmplitude /= typicalSA;
	  ValueExact /= typicalWF;
	  if (Norm(ValueExact)>100.0)
	    {
	      cout << i << ": excluding large Psi: " << ValueExact << endl;
	      if (State.GetVectorDimension()>0) // have exact state!
		{
		  cout << "Checking if abnormal result is reproduced!" << endl;
		  FermionOnSphere Space (NbrFermions, Lz, LzMax);
		  ParticleOnSphereFunctionBasis Basis(LzMax,ParticleOnSphereFunctionBasis::LeftHanded);  
		  QHEParticleWaveFunctionOperation Operation(&Space, &State, &Positions, &Basis, /* TimeCoherence */ -1);
		  Operation.ApplyOperation(Architecture.GetArchitecture());      
		  ValueExact = Operation.GetScalar();
		  cout << "Comparing: " << ValueExact << " (new) to "<< rawExact << " (old) (ratio " << ValueExact/rawExact << ")" <<endl;
		  cout << "Minimal distance: " << MinDist(Positions);
		}
	    }
	  else
	    {
	      NormTrialObs.Observe(SqrNorm(TrialValue)/CurrentSamplingAmplitude,(double)sampleCount);
	      NormExactObs.Observe(SqrNorm(ValueExact)/CurrentSamplingAmplitude,(double)sampleCount);
	      OverlapObs.Observe(Conj(TrialValue)*ValueExact/CurrentSamplingAmplitude,(double)sampleCount);
	    }
//              if (sampleCount>10)
// 	    {
// 	      cout << "Total "<<sampleCount<<" samples in these coordinates: " << endl;
// 	      cout << "Psi^2=" << SqrNorm(ValueExact) <<" trial^2="<<CurrentSamplingAmplitude<<endl;
// 	      cout << "Minimal distance: " << MinDist(Positions);
// 	      //cout << Positions<< endl;
// 	    }		
	}
      if (i>NbrIter) cout << "Attention, step number limited by NbrIter!" << endl;
      History->RewindHistory();      
      
      cout << " final results :" << endl;
      cout << "overlap:   " << OverlapValue(OverlapObs, NormTrialObs, NormExactObs) << " +/- "
	   << OverlapError(OverlapObs, NormTrialObs, NormExactObs) << endl;
      cout << "overlap^2: " << SqrNorm(OverlapValue(OverlapObs, NormTrialObs, NormExactObs)) << " +/- " <<
	2.0*Norm(OverlapValue(OverlapObs, NormTrialObs, NormExactObs))*OverlapError(OverlapObs, NormTrialObs, NormExactObs) << endl;
      cout << "Processed a total of "<<totalSampleCount<<" MC samples" << endl;
      delete History;
      return 0;
    }


  if (HistoryMode ==3)
    {
      if (!(WaveFunctionManager.GetWaveFunctionType() & QHEWaveFunctionManager::TrialWaveFunction))
	{
	  cout << "This type of wavefunction cannot be optimized." << endl;
	  return -1;
	}
      // code to optimize paired state with available samples:
      bool varyMR = Manager.GetBoolean("varyMR");
      WaveFunctionOverlapOptimizer *Optimizer =
	new WaveFunctionOverlapOptimizer( (Abstract1DComplexTrialFunction*) TestWaveFunction,
					  HistoryFileName, NbrFermions, /* excludeLastParameter */ !varyMR,
					  Manager.GetInteger("maxPoints"), Manager.GetInteger("limitSamples"));
      RealVector optimalParameters( ((Abstract1DComplexTrialFunction*) TestWaveFunction)->GetNbrParameters());
      Complex optimalOverlap;
      Optimizer->GetMaximumSqrOverlap(optimalParameters, optimalOverlap);
      return 0;
    }
      
  
  FermionOnSphere Space (NbrFermions, Lz, LzMax);
  ParticleOnSphereFunctionBasis Basis(LzMax,ParticleOnSphereFunctionBasis::LeftHanded);  
  AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (29457);
  double PreviousSamplingAmplitude;
  
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

  // MC observables:
  Complex Overlap(0.0);
  ComplexObservable ScalarProduct(100);    
  RealObservable NormExact(100);
  // stuff
  Complex TmpMetropolis;
  int NextCoordinates = 0;
  int Accepted = 0;
  bool NoTimeCoherence=!(Manager.GetBoolean("with-timecoherence"));
  int TimeCoherence;

  // testing
  ofstream Output("test",ios::out | ios::binary);
  char CC='c';
  for(int i=0; i<10; ++i)
    WriteLittleEndian(Output, CC);
  Output.close();

  Output.open("test", ios::in | ios::out | ios::binary );
  Output.seekp(6);
  CC='d';
  WriteLittleEndian(Output, CC);
  Output.close();

  if (HistoryMode == 4) // continuing to work on old History
    {
      if (Manager.GetString("exact-state") == 0)
	{
	  cout << "QHEFermionOverlap requires an exact state" << endl;
	  return -1;
	}  
      if (State.ReadVector (Manager.GetString("exact-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("exact-state") << endl;
	  return -1;      
	}
      History = new MCHistoryRecord (HistoryFileName, 2*NbrFermions, PreviousSamplingAmplitude, &(Particles->GetPositions()[0]), ValueExact, NULL);
      // initialize function values at initial positions: - trial function
      TrialValue = (*TestWaveFunction)(Particles->GetPositions());  

      // - exact function
      TimeCoherence = -1;
      QHEParticleWaveFunctionOperation Operation(&Space, &State, &(Particles->GetPositions()), &Basis, TimeCoherence);
      Operation.ApplyOperation(Architecture.GetArchitecture());      
      if (SqrNorm(ValueExact - Operation.GetScalar()) > 1e-10)
	{
	  cout << "Exact wavefunction from record not exactly reproduced: "
	       << ValueExact << " vs " <<Operation.GetScalar()<<endl;
	  
	  cout << "Exact state may have changed!" << endl;
	}
      
      CurrentSamplingAmplitude = SqrNorm(TrialValue);;
      if ( fabs(PreviousSamplingAmplitude -  CurrentSamplingAmplitude) > 1e-10)
	{
	  cout << "SamplingAmplitude from record not exactly reproduced: "
	       << PreviousSamplingAmplitude << " vs " <<CurrentSamplingAmplitude<<endl;
	  cout << "Sampling function may have changed" << endl;
	}
    }
  else
    {
      // - exact function
      TimeCoherence = -1;
      QHEParticleWaveFunctionOperation Operation(&Space, &State, &(Particles->GetPositions()), &Basis, TimeCoherence);
      Operation.ApplyOperation(Architecture.GetArchitecture());      
      ValueExact = Operation.GetScalar();
      // initialize function values at initial positions: - trial function
      TrialValue = (*TestWaveFunction)(Particles->GetPositions());  
      PreviousSamplingAmplitude = SqrNorm(TrialValue);
      CurrentSamplingAmplitude = PreviousSamplingAmplitude;
    } 
  
  
  
  for (int i = 0; i < NbrIter; ++i)
    {
      // make a random move of particle NextCoordinates
      Particles->Move(NextCoordinates);
      // evaluate new trial wavefunction value
      TmpMetropolis = (*TestWaveFunction)(Particles->GetPositions());
      // accept or reject move according to probability |Psi_new|^2  / |Psi_old|^2
      CurrentSamplingAmplitude = SqrNorm(TmpMetropolis);
      if ((CurrentSamplingAmplitude > PreviousSamplingAmplitude) ||
	  ((RandomNumber->GetRealRandomNumber() * PreviousSamplingAmplitude) < CurrentSamplingAmplitude))
	{
	  PreviousSamplingAmplitude = CurrentSamplingAmplitude;
	  TrialValue = TmpMetropolis;
	  ++Accepted;	  
	  if (SampleDensity==1) // do not recalculate Exact or count steps here if using spread out samples
	    {
	      // recalculate exact function value:
	      TimeCoherence = NextCoordinates;
	      if (NoTimeCoherence) TimeCoherence = -1;      
	      QHEParticleWaveFunctionOperation Operation(&Space, &State, &(Particles->GetPositions()), &Basis, TimeCoherence);
	      Operation.ApplyOperation(Architecture.GetArchitecture());      
	      ValueExact = Operation.GetScalar();
	      if (History)
		History->RecordAcceptedStep( CurrentSamplingAmplitude, Particles->GetPositions(), ValueExact);
	    }
	}
      else
	{
	  Particles->RestoreMove();
	  CurrentSamplingAmplitude = PreviousSamplingAmplitude;
	  if ((History)&&(SampleDensity==1)) // do not count steps here if using spread out samples
	    History->RecordRejectedStep();
	}
      if ((SampleDensity>1)&&( (i % SampleDensity)==0 )) // recalculate and record spread out samples here if SampleDensity>1
	{
	  // recalculate exact function value:
	  TimeCoherence = NextCoordinates;
	  if (NoTimeCoherence) TimeCoherence = -1;      
	  QHEParticleWaveFunctionOperation Operation(&Space, &State, &(Particles->GetPositions()), &Basis, TimeCoherence);
	  Operation.ApplyOperation(Architecture.GetArchitecture());      
	  ValueExact = Operation.GetScalar();
	  if (History) History->RecordAcceptedStep( CurrentSamplingAmplitude, Particles->GetPositions(), ValueExact);
	}
      // determine next particle to move
      NextCoordinates = (int) (((double) NbrFermions) * RandomNumber->GetRealRandomNumber());
      if (NextCoordinates == NbrFermions) --NextCoordinates;
      
      // note observations:
      if ( (i % SampleDensity)==0 )
	{
	  double norm = SqrNorm(ValueExact)/CurrentSamplingAmplitude;
	  NormExact << norm;
	  Overlap = Conj(TrialValue)*ValueExact/CurrentSamplingAmplitude;
	  ScalarProduct << Overlap;
	}
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
  delete History;
  cout << " final results :" << endl;
  cout << "overlap:   " << OverlapValue(ScalarProduct, NormExact) << " +/- "
       << OverlapError(ScalarProduct, NormExact) << endl;
  cout << "overlap^2: " << SqrNorm(OverlapValue(ScalarProduct, NormExact)) << " +/- " <<
    2.0*Norm(OverlapValue(ScalarProduct, NormExact))*OverlapError(ScalarProduct, NormExact) << endl;
  cout << "Fraction of moves accepted: " << (double) Accepted / NbrIter *100.0<< " %"<< endl;
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
  delete Particles;
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

Complex OverlapValue(WeightedComplexObservable &ScalarProduct, WeightedRealObservable &NormObs1, WeightedRealObservable &NormObs2)
//Complex OverlapValue(ComplexObservable &ScalarProduct, RealObservable &NormObs1, RealObservable &NormObs2)
{
  return Complex(ScalarProduct.Average()/sqrt(NormObs1.Average()*NormObs2.Average()));
}

double OverlapError(WeightedComplexObservable &ScalarProduct, WeightedRealObservable &NormObs1, WeightedRealObservable &NormObs2)
//double OverlapError(ComplexObservable &ScalarProduct, RealObservable &NormObs1, RealObservable &NormObs2)
{
  double norm1 = NormObs1.Average();
  double norm2 = NormObs2.Average();
  double norm = sqrt(norm1*norm2);
  double prod = Norm(ScalarProduct.Average());
  return sqrt( DSQR(ScalarProduct.ErrorEstimate()/norm) + DSQR(prod*NormObs1.ErrorEstimate()*NormObs2.Average()/2.0/norm/norm/norm) + DSQR(prod*NormObs2.ErrorEstimate()*NormObs1.Average()/2.0/norm/norm/norm));
}
  
