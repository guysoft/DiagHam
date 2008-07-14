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

#include "Options/Options.h"

#include "Tools/FQHEWaveFunction/PfaffianOnDiskWaveFunction.h"
#include "Tools/FQHEWaveFunction/LaughlinOnDiskWaveFunction.h"

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


void RandomZ (RealVector& positions, double scale, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator);

void RandomZOneCoordinate(RealVector& positions, double scale, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator);



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
  (*SystemGroup) += new SingleIntegerOption  ('l', "max-momentum", "maximum angular monentum a particle can have", 10);
  (*SystemGroup) += new SingleIntegerOption  ('k', "k-value", "k index of the SU(k) symmetry group", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "intra-corr", "power of the intra-color correlations", 3);  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "inter-corr", "power of the inter-color correlations", 2);  
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  
  (*SystemGroup) += new BooleanOption ('\n', "reverse-flux", "use reverse flux attachment for the composite fermions");  
  (*SystemGroup) += new BooleanOption ('\n', "jain-cf", "use composite fermion state instead of the symetrized state");  
  (*SystemGroup) += new SingleStringOption  ('\n', "load-permutations", "read all the permutations needed to compute the reference wave function from a file");  
  (*SystemGroup) += new SingleStringOption  ('\n', "save-permutations", "file name where all the permutations needed to compute the reference wave function have to be stored");  
 
  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "nbr-warmup", "number of Monte Carlo iterations that have to be done before evaluating the energy (i.e. warm up sequence)", 10000);
  (*MonteCarloGroup) += new BooleanOption  ('r', "resume", "resume from a previous run");
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "record-step", "number of iteration between two consecutive result recording of energy value (0 if no on-disk recording is needed)", 0);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "record-file", "name of the file where energy recording has to be done", "montecarlo.dat");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "with-timecoherence", "use time coherence between two successive evaluation of the wave function");
  (*MonteCarloGroup) += new SingleStringOption ('\n', "random-file", "name of the file where random number to use are stored (use internal random generator if no file name is provided)");
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "random-seek", "if usage of a random number file is activiated, jump the first random numbers up to the seek position", 0);
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
  int KValue = Manager.GetInteger("k-value");
  if ((NbrParticles % KValue) != 0)
    {
      cout << "the number of particles has to be a multiple of the number of colors (i.e. k)" << endl;
      return -1;
    }
  int NbrParticlePerColor = NbrParticles / KValue;
  int IntraCorrelation = Manager.GetInteger("intra-corr");
  int InterCorrelation = Manager.GetInteger("inter-corr");
  int NbrWarmUpIter = ((SingleIntegerOption*) Manager["nbr-warmup"])->GetInteger();
  int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  bool InvertFlag = Manager.GetBoolean("reverse-flux");
  bool ResumeFlag = Manager.GetBoolean("resume");
  int LzMax = NbrParticlePerColor * (((KValue - 1) * InterCorrelation) + IntraCorrelation) - IntraCorrelation - KValue + 1;
  bool StatisticFlag = true;
  int NbrMomenta = ((SingleIntegerOption*) Manager["max-momentum"])->GetInteger() + 1;
  int RecordStep = Manager.GetInteger("record-step");
  double InvSqrLengthScale = sqrt(3.0 * ((double) NbrParticles)) / sqrt(2.0 * ((double) ((NbrMomenta - 1))));

  AbstractQHEParticle* ExactSpace = 0;
  RealVector* ExactState = 0;
  AbstractFunctionBasis* ExactBasis = 0;



//  Abstract1DComplexFunction* WaveFunction = new PfaffianOnDiskWaveFunction(NbrParticles);
  Abstract1DComplexFunction* WaveFunction = new LaughlinOnDiskWaveFunction(NbrParticles, 3);
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
      RandomNumber = new StdlibRandomNumberGenerator (29457);
    }


       
  double Normalization = 0.0;
  double ErrorNormalization = 0.0;
  Complex Tmp;
  double Tmp2;
  double* Momentum = new double[NbrMomenta];
  double* ErrorMomentum = new double[NbrMomenta];
  double* Tmp3 = new double[NbrMomenta];
  int NextCoordinates = 0;
  RealVector TmpZ (NbrParticles * 2, true);
  RealVector TmpPositions (NbrParticles * 2, true);
  double* ValueExact = new double[NbrMomenta];
  double PreviousProbabilities = 0.0;
  double CurrentProbabilities = 0.0;
  double TotalProbability = 0.0;
  int Acceptance = 0;
  double AcceptanceRate = 1.0;
  int InitialNbrIter = 0;
  if (ResumeFlag == false)
    {
      for (int j = 0; j < NbrMomenta; ++j)
	{
	  Momentum[j] = 0.0;
	  ErrorMomentum[j] = 0.0;
	  Tmp3[j] = 0.0;	   
	}
      RandomZ (TmpZ, InvSqrLengthScale, NbrParticles, RandomNumber);
      Tmp = (*WaveFunction)(TmpZ);
      PreviousProbabilities = Norm(Tmp);
      CurrentProbabilities = PreviousProbabilities;
      TotalProbability = PreviousProbabilities;
      Acceptance = 0;
      AcceptanceRate = 1.0;
      if (NbrWarmUpIter > 0)
	cout << "starting warm-up sequence" << endl;
      for (int i = 1; i < NbrWarmUpIter; ++i)
	{      
	  double PreviousCoordinates1 = TmpZ[NextCoordinates << 1];
	  double PreviousCoordinates2 = TmpZ[1 + (NextCoordinates << 1)];
	  Complex TmpMetropolis;
	  RandomZOneCoordinate(TmpZ, InvSqrLengthScale, NextCoordinates, RandomNumber);
	  TmpMetropolis = (*WaveFunction)(TmpZ);
	  CurrentProbabilities = SqrNorm(TmpMetropolis);
	  cout << CurrentProbabilities << " " << TmpZ[2 * NextCoordinates] << " " << TmpZ[1 + 2 * NextCoordinates] << " " << endl;
	  if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	    {
	      PreviousProbabilities = CurrentProbabilities;
	      ++Acceptance;
	    }
	  else
	    {
	      TmpZ[NextCoordinates << 1] = PreviousCoordinates1;
	      TmpZ[1 + (NextCoordinates << 1)] = PreviousCoordinates2;
	      CurrentProbabilities = PreviousProbabilities;
	    }
	  NextCoordinates = (int) (((double) NbrParticles) * RandomNumber->GetRealRandomNumber());
	  if ((i % 1000) == 0)
	    {
	      AcceptanceRate = ((double) Acceptance) / ((double) i);
	      cout << Acceptance << " / " << i << " = " << ((100.0 * ((double) Acceptance)) / ((double) i)) << "%" << endl;
	    }
	}
      
      if (NbrWarmUpIter > 0)
	cout << "warm-up sequence is over" << endl;
      Acceptance = 0;
      if ((RecordStep != 0) && (ResumeFlag == false))
	{
	  ofstream MomentumRecordFile;
	  MomentumRecordFile.precision(14);
	  MomentumRecordFile.open(((SingleStringOption*) Manager["record-file"])->GetString(), ios::out | ios::binary);
	  MomentumRecordFile << "# Monte Carlo calculation of the " << endl
			    << "# step overlap.Re overlap.Im error.Re error.Im [(scalar_product.Re scalar_product.Im error_scalar_product.Re error_scalar_product.Im normalization_exact error_normalization_exact) per state] normalization error_normalization" << endl;
	}
    }
  else
    {
      ifstream MCState;
      MCState.open("mcstate.dat", ios::in | ios::binary);
      ReadLittleEndian(MCState, InitialNbrIter);
      unsigned long TmpNumber;
      ReadLittleEndian(MCState, TmpNumber);
      ReadLittleEndian(MCState, Acceptance);
      ReadLittleEndian(MCState, PreviousProbabilities);
      ReadLittleEndian(MCState, CurrentProbabilities);
      ReadLittleEndian(MCState, TotalProbability);
      ReadLittleEndian(MCState, NextCoordinates);
      for (int j = 0; j < (NbrParticles >> 1); ++j)
	{		   	       
	  ReadLittleEndian(MCState, TmpZ[j]);
	}
      ReadLittleEndian(MCState, Normalization);
      ReadLittleEndian(MCState, ErrorNormalization);
      for (int j = 0; j < NbrMomenta; ++j)
	{		   	       
	  ReadLittleEndian(MCState, Momentum[j]);
	  ReadLittleEndian(MCState, ErrorMomentum[j]);
	}
      MCState.close();	     
    }
  for (int i = InitialNbrIter; i < NbrIter; ++i)
    {
      double PreviousCoordinates1 = TmpZ[NextCoordinates << 1];
      double PreviousCoordinates2 = TmpZ[1 + (NextCoordinates << 1)];
      RandomZOneCoordinate(TmpZ, 1.0,NextCoordinates, RandomNumber);
      Complex TmpMetropolis;
      TmpMetropolis = (*WaveFunction)(TmpZ);
      CurrentProbabilities = SqrNorm(TmpMetropolis);
      if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	{
	  PreviousProbabilities = CurrentProbabilities;
	  Tmp = TmpMetropolis;
	  double TmpSqrRadius = ((TmpZ[0] * TmpZ[0]) + (TmpZ[1] * TmpZ[1]));
	  ValueExact[0] = 1.0;//exp (-0.5 * TmpSqrRadius);
	  TmpSqrRadius /= InvSqrLengthScale * InvSqrLengthScale;
	  for (int j = 1; j < NbrMomenta; ++j)
	    ValueExact[j] = ValueExact[j - 1] * TmpSqrRadius / (((double) j));
	  ++Acceptance;
	}
      else
	{
	  TmpZ[NextCoordinates << 1] = PreviousCoordinates1;
	  TmpZ[1 + (NextCoordinates << 1)] = PreviousCoordinates2;
	  CurrentProbabilities = PreviousProbabilities;
	}
      TotalProbability += CurrentProbabilities;
      NextCoordinates = (int) (((double) NbrParticles) * RandomNumber->GetRealRandomNumber());
      if (NextCoordinates == NbrParticles)
	--NextCoordinates;
      Tmp2 = (Tmp.Re * Tmp.Re) + (Tmp.Im * Tmp.Im);
      Tmp2 /= CurrentProbabilities;
      Tmp2 = 1.0;
      for (int j = 0; j < NbrMomenta; ++j)
	{
	  Tmp3[j] = (Tmp2 * ValueExact[j]);
	  Momentum[j] += Tmp3[j];
	  ErrorMomentum[j] += Tmp3[j] * Tmp3[j];
	}
      Normalization += Tmp2;
      ErrorNormalization += Tmp2 * Tmp2;
      if ((i > 0) && ((RecordStep != 0) && ((i % RecordStep) == 0)))
	{
	  ofstream MomentumRecordFile;
	  MomentumRecordFile.precision(14);
	  MomentumRecordFile.open(((SingleStringOption*) Manager["record-file"])->GetString(), ios::out | ios::binary | ios::app);
	  MomentumRecordFile << i ;
	  double Tmp6 = Normalization  / ((double) i);
	  double Tmp7 = sqrt( ((ErrorNormalization / ((double) i))  -  (Tmp6 * Tmp6)) / ((double) i) );	  
	  for (int j = 0; j < NbrMomenta; ++j)
	    {
	      double Tmp4 = Momentum[j] / ((double) i);
	      double Tmp5 = sqrt(((ErrorMomentum[j]/ ((double) i)) - (Tmp4 * Tmp4)) / ((double) i) );
	      Tmp5 /= Tmp4;
	      Tmp5 = fabs(Tmp5);
	      Tmp5 += Tmp7 / Tmp6;
	      Tmp4 /= Tmp6;	  
	      Tmp5 *= Tmp4;
	      MomentumRecordFile << " " << Tmp4 << " " << Tmp5 << " " << Momentum[j] << " " << ErrorMomentum[j];
	    }
	  MomentumRecordFile << " " << Normalization << " " << ErrorNormalization << endl;
	  MomentumRecordFile.close();
	  ofstream MCState;
	  MCState.open("mcstate.dat", ios::out | ios::binary);
	  WriteLittleEndian(MCState, i);
	  unsigned long TmpNumber = RandomNumber->GetNbrGeneratedNumbers();
	  WriteLittleEndian(MCState, TmpNumber);
	  WriteLittleEndian(MCState, Acceptance);
	  WriteLittleEndian(MCState, PreviousProbabilities);
	  WriteLittleEndian(MCState, CurrentProbabilities);
	  WriteLittleEndian(MCState, TotalProbability);
	  WriteLittleEndian(MCState, NextCoordinates);
	  for (int j = 0; j < (NbrParticles >> 1); ++j)
	    {		   	       
	      WriteLittleEndian(MCState, TmpZ[j]);
	      WriteLittleEndian(MCState, TmpZ[j]);
	    }
	  WriteLittleEndian(MCState, Normalization);
	  WriteLittleEndian(MCState, ErrorNormalization);
	  for (int j = 0; j < NbrMomenta; ++j)
	    {		   	       
	      WriteLittleEndian(MCState, Momentum[j]);
	      WriteLittleEndian(MCState, ErrorMomentum[j]);
	    }
	  MCState.close();	     
	}
      if ((i > 0) && ((i % (((SingleIntegerOption*) Manager["display-step"])->GetInteger())) == 0))
	{
	  cout << " i = " << i << endl;
	  double Tmp6 = Normalization  / ((double) i);
	  double Tmp7 = sqrt( ((ErrorNormalization / ((double) i))  -  (Tmp6 * Tmp6)) / ((double) i) );	  
	  for (int j = 0; j < NbrMomenta; ++j)
	    {		   
	      double Tmp4 = Momentum[j] / ((double) i);
	      double Tmp5 = sqrt( ((ErrorMomentum[j] / ((double) i)) - (Tmp4 * Tmp4)) / ((double) i) );
	      Tmp5 /= Tmp4;
	      Tmp5 = fabs(Tmp5);
	      Tmp5 += Tmp7 / Tmp6;
	      Tmp4 /= Tmp6;	  
	      Tmp5 *= Tmp4;
	      cout << "n_" << j << " : " << Tmp4 << " +/- " << Tmp5 << endl;
	    }
	  cout << "acceptance rate = " << (((double) Acceptance) / (((double) i))) << endl;
	  cout << "-----------------------------------------------" << endl;
	}
    } 
  if (((RecordStep != 0) && ((NbrIter % RecordStep) == 0)))
    {
      ofstream MomentumRecordFile;
      MomentumRecordFile.precision(14);
      MomentumRecordFile.open(((SingleStringOption*) Manager["record-file"])->GetString(), ios::out | ios::binary | ios::app);
      MomentumRecordFile << NbrIter ;
      double Tmp6 = Normalization  / ((double) NbrIter);
      double Tmp7 = sqrt( ((ErrorNormalization / ((double) NbrIter))  -  (Tmp6 * Tmp6)) / ((double) NbrIter) );	  
      for (int j = 0; j < NbrMomenta; ++j)
	{
	  double Tmp4 = Momentum[j] / ((double) NbrIter);
	  double Tmp5 = sqrt(((ErrorMomentum[j]/ ((double) NbrIter)) - (Tmp4 * Tmp4)) / ((double) NbrIter) );
	  Tmp5 /= Tmp4;
	  Tmp5 = fabs(Tmp5);
	  Tmp5 += Tmp7 / Tmp6;
	  Tmp4 /= Tmp6;	  
	  Tmp5 *= Tmp4;
	  MomentumRecordFile << " " << Tmp4 << " " << Tmp5 << " " << Momentum[j] << " " << ErrorMomentum[j];
	}
      MomentumRecordFile << " " << Normalization << " " << ErrorNormalization << endl;
      MomentumRecordFile.close();
    }
  cout << " final results :" << endl;
  double Tmp6 = Normalization  / ((double) NbrIter);
  double Tmp7 = sqrt( ((ErrorNormalization / ((double) NbrIter))  -  (Tmp6 * Tmp6)) / ((double) NbrIter) );	  
  double Tmp4Ref = Momentum[NbrMomenta - 1] / ((double) NbrIter);
  double Tmp5Ref = sqrt(((ErrorMomentum[NbrMomenta - 1]/ ((double) NbrIter)) - (Tmp4Ref * Tmp4Ref)) / ((double) NbrIter) );
  Tmp5Ref /= Tmp4Ref;
  Tmp5Ref = fabs(Tmp5Ref);
  Tmp5Ref += Tmp7 / Tmp6;
  Tmp4Ref /= Tmp6;	  
  Tmp5Ref *= Tmp4Ref;
  for (int j = 0; j < NbrMomenta; ++j)
    {
      double Tmp4 = Momentum[j] / ((double) NbrIter);
      double Tmp5 = sqrt(((ErrorMomentum[j]/ ((double) NbrIter)) - (Tmp4 * Tmp4)) / ((double) NbrIter) );
      Tmp5 /= Tmp4;
      Tmp5 = fabs(Tmp5);
      Tmp5 += Tmp7 / Tmp6;
      Tmp4 /= Tmp6;	  
      Tmp5 *= Tmp4;
      double TmpCoef = 1.0;
      for (int m = NbrMomenta - 1; m > j; --m)
	TmpCoef *= 2.0 * (InvSqrLengthScale * InvSqrLengthScale);
      cout << "n_" << j << " : " << Tmp4 << " +/- " << Tmp5 << "   ,   " << (TmpCoef * Tmp4 / Tmp4Ref) << " +/- " << ((Tmp5 / Tmp4Ref) + (Tmp5Ref * Tmp4 / (Tmp4Ref * Tmp4Ref))) << endl;
    }
  cout << "-----------------------------------------------" << endl;
  
  

  return 0;
}

void RandomZ (RealVector& positions, double scale, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  for (int j = 0; j < nbrParticles; ++j)
    {
      double x = scale * sqrt(-2.0 * log(1.0 - randomNumberGenerator->GetRealRandomNumber()));
      double y = 2.0 * M_PI * randomNumberGenerator->GetRealRandomNumber();
      positions[2 * j] = x * cos(y);
      positions[(2 * j) + 1] = x * sin(y);
    }
}

void RandomZOneCoordinate(RealVector& positions, double scale, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  coordinate *= 2;
  double x = scale * sqrt(-2.0 * log(1.0 - randomNumberGenerator->GetRealRandomNumber()));
  double y = 2.0 * M_PI * randomNumberGenerator->GetRealRandomNumber();
  positions[coordinate] = x * cos(y);
  ++coordinate;
  positions[coordinate] = x * sin(y);
}

