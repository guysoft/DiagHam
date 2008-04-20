#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"

#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/RandomNumber/FileRandomNumberGenerator.h"

#include "Options/Options.h"

#include "Tools/FQHEWaveFunction/FQHESphereSymmetrizedSUKToU1WaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESU2HalperinPermanentOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESU3HalperinPermanentOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESU4HalperinPermanentOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/JainCFFilledLevelOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/HundRuleCFStates.h"
#include "Tools/FQHEWaveFunction/SU3HalperinOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/MooreReadOnSphereWaveFunction.h"

#include "Vector/ComplexVector.h"

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
  OptionManager Manager ("FQHESphereSUKToU1MCOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MonteCarloGroup = new OptionGroup ("Monte Carlo options");

   ArchitectureManager Architecture;

 Manager += SystemGroup;
  Manager += MonteCarloGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;


  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 10);
  (*SystemGroup) += new SingleIntegerOption  ('k', "k-value", "k index of the SU(k) symmetry group", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "intra-corr", "power of the intra-color correlations", 3);  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "inter-corr", "power of the inter-color correlations", 2);  
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  
  (*SystemGroup) += new BooleanOption ('\n', "test-symmetry", "check the test wave function is symmetric/antisymmetric");  
  (*SystemGroup) += new BooleanOption ('\n', "reverse-flux", "use reverse flux attachment for the composite fermions");  
  (*SystemGroup) += new BooleanOption ('\n', "jain-cf", "use composite fermion state instead of the symetrized state");  
  (*SystemGroup) += new SingleStringOption  ('\n', "use-exact", "file name of an exact state that has to be used as test wave function");  
 
  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "nbr-warmup", "number of Monte Carlo iterations that have to be done before evaluating the energy (i.e. warm up sequence)", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "record-step", "number of iteration between two consecutive result recording of energy value (0 if no on-disk recording is needed)", 0);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "record-file", "name of the file where energy recording has to be done", "montecarlo.dat");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "with-timecoherence", "use time coherence between two successive evaluation of the wave function");
  (*MonteCarloGroup) += new BooleanOption  ('\n', "show-details", "show intermediate values used for overlap calculation", false);
  (*MonteCarloGroup) += new SingleStringOption ('\n', "random-file", "name of the file where random number to use are stored (use internal random generator if no file name is provided)");
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "random-seek", "if usage of a random number file is activiated, jump the first random numbers up to the seek position", 0);
  (*MonteCarloGroup) += new BooleanOption  ('\n', "weight-symmetrized" , "use the norm of the symmetrized wave fonction as probalbility density instead of the exact state");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereSUKToU1MCOverlap -h" << endl;
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
  bool OverlapFlag = true;
  if (((BooleanOption*) Manager["test-symmetry"])->GetBoolean() == true)
    {
      OverlapFlag = false;
    }
  int NbrWarmUpIter = ((SingleIntegerOption*) Manager["nbr-warmup"])->GetInteger();
  int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
  bool InvertFlag = Manager.GetBoolean("reverse-flux");
  int LzMax = NbrParticlePerColor * (((KValue - 1) * InterCorrelation) + IntraCorrelation) - IntraCorrelation - KValue + 1;
  bool UseExactFlag = false;
  bool StatisticFlag = true;
  bool UseBaseAsWeightFlag = ((BooleanOption*) Manager["weight-symmetrized"])->GetBoolean();

  AbstractQHEParticle* ExactSpace = 0;
  RealVector* ExactState = 0;
  AbstractFunctionBasis* ExactBasis = 0;
  if (((SingleStringOption*) Manager["use-exact"])->GetString() != 0)
    {
      UseExactFlag = true;
      if (StatisticFlag == true)
	ExactSpace = new FermionOnSphere (NbrParticles, 0, LzMax);
      else
	ExactSpace = new BosonOnSphere (NbrParticles, 0, LzMax);
      ExactState = new RealVector;
      if (ExactState->ReadVector (((SingleStringOption*) Manager["use-exact"])->GetString()) == false)
	{
	  cout << "can't open vector file " << ((SingleStringOption*) Manager["use-exact"])->GetString() << endl;
	  return -1;      
	}
      if (ExactSpace->GetHilbertSpaceDimension() != ExactState->GetVectorDimension())
	{
	  cout << "dimension mismatch : hilbert space = " << ExactSpace->GetHilbertSpaceDimension() << ", exact state = " << ExactState->GetVectorDimension() << endl;
	  return -1;
	}
      ExactBasis = new ParticleOnSphereFunctionBasis (LzMax);	
    }

  Abstract1DComplexFunctionOnSphere* BaseFunction = 0;
  switch (KValue)
    {
    case 2:
      {
	BaseFunction = new FQHESU2HalperinPermanentOnSphereWaveFunction (NbrParticlePerColor, NbrParticlePerColor, 
									 IntraCorrelation - 1, IntraCorrelation - 1, InterCorrelation - 1, InvertFlag);
      }
      break;
    case 3:
      {
	BaseFunction = new FQHESU3HalperinPermanentOnSphereWaveFunction(NbrParticlePerColor, NbrParticlePerColor, NbrParticlePerColor, 
									IntraCorrelation - 1, IntraCorrelation - 1, IntraCorrelation - 1, 
									InterCorrelation - 1, InterCorrelation - 1, InterCorrelation - 1, InvertFlag);
      }
      break;
    case 4:
      {
	BaseFunction = new FQHESU4HalperinPermanentOnSphereWaveFunction(NbrParticlePerColor, NbrParticlePerColor, NbrParticlePerColor, 
									NbrParticlePerColor,
									IntraCorrelation - 1, IntraCorrelation - 1, IntraCorrelation - 1,
                                                                        IntraCorrelation - 1, InterCorrelation - 1, InterCorrelation - 1,
									InterCorrelation - 1, InterCorrelation - 1, InterCorrelation - 1,
									InterCorrelation - 1, InvertFlag);
      }
      break;
    default:
      {
	cout << "invalid or unsupported number of colors (i.e. k)" << endl;
	return -1;
      }
    }
  Abstract1DComplexFunctionOnSphere* SymmetrizedFunction = 0;
  if (((BooleanOption*) Manager["jain-cf"])->GetBoolean() == true)
    {
      SymmetrizedFunction = new JainCFFilledLevelOnSphereWaveFunction(NbrParticles, KValue, 2);
    }
  else
    {
      SymmetrizedFunction = new FQHESphereSymmetrizedSUKToU1WaveFunction (NbrParticles, KValue, BaseFunction, true);      
    }
  Abstract1DComplexFunctionOnSphere* TestFunction;
  AbstractRandomNumberGenerator* RandomNumber = 0;
  if (((SingleStringOption*) Manager["random-file"])->GetString() != 0)
    {
      RandomNumber = new FileRandomNumberGenerator(((SingleStringOption*) Manager["random-file"])->GetString(), (NbrWarmUpIter * 4) + (NbrIter * 4) + 2000, 
						   ((SingleIntegerOption*) Manager["random-seek"])->GetInteger());
    }
  else
    {
      RandomNumber = new StdlibRandomNumberGenerator (29457);
    }

  if (UseExactFlag == true)
    {
      TestFunction = 0;
    }
  else
    {
      if (InvertFlag == false)
	{
	  TestFunction = new JainCFFilledLevelOnSphereWaveFunction(NbrParticles, KValue, 2);
	}
      else
	{
	  TestFunction = new HundRuleCFStates (NbrParticles,  - (NbrParticles / 2) + 2, 1);
	}
    }

  StdlibRandomNumberGenerator RandomNumberGenerator(29457);
  
  if (OverlapFlag == true)
    {
       int RecordStep = Manager.GetInteger("record-step");
       
       Complex* RecordedOverlap = 0;
       Complex* RecordedOverlapError = 0;
       if (RecordStep != 0)
	 {
	   RecordedOverlap = new Complex [(NbrIter / RecordStep) + 1];
	   RecordedOverlapError = new Complex [(NbrIter / RecordStep) + 1];
	 }
       int RecordIndex = 0;
       Complex Overlap;
       Complex ErrorOverlap;
       double Normalization = 0.0;
       double ErrorNormalization = 0.0;
       double NormalizationExact = 0.0;
       double ErrorNormalizationExact = 0.0;
       Complex Tmp;
       Complex Tmp3;
       double Tmp2;
       double Tmp2bis;
       int NextCoordinates = 0;
       ComplexVector TmpUV (NbrParticles * 2, true);
       RealVector TmpPositions (NbrParticles * 2, true);
       RandomUV (TmpUV, TmpPositions, NbrParticles, RandomNumber);
       if (UseExactFlag == false)
	 Tmp = TestFunction->CalculateFromSpinorVariables(TmpUV);
       else
	 {
	   if (UseBaseAsWeightFlag == false)
	     {
	       QHEParticleWaveFunctionOperation Operation(ExactSpace, ExactState, &TmpPositions, ExactBasis);
	       Operation.ApplyOperation(Architecture.GetArchitecture());      
	       Tmp = Operation.GetScalar();
	     }
	   else
	     Tmp = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
	 }
       Complex ValueExact;
       if (UseBaseAsWeightFlag == false)
	 ValueExact = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
       else
	 {
	   QHEParticleWaveFunctionOperation Operation(ExactSpace, ExactState, &TmpPositions, ExactBasis);
	   Operation.ApplyOperation(Architecture.GetArchitecture());      
	   ValueExact = Operation.GetScalar();
	 }
       double PreviousProbabilities = Norm(Tmp);
       double CurrentProbabilities = PreviousProbabilities;
       double TotalProbability = PreviousProbabilities;
       int Acceptance = 0;
       double AcceptanceRate = 1.0;
       if (NbrWarmUpIter > 0)
	 cout << "starting warm-up sequence" << endl;
       for (int i = 1; i < NbrWarmUpIter; ++i)
	 {      
	   Complex PreviousCoordinatesU = TmpUV[NextCoordinates << 1];
	   Complex PreviousCoordinatesV = TmpUV[1 + (NextCoordinates << 1)];
	   double PreviousCoordinates1 = TmpPositions[NextCoordinates << 1];
	   double PreviousCoordinates2 = TmpPositions[1 + (NextCoordinates << 1)];
	   Complex TmpMetropolis;
	   RandomUVOneCoordinate(TmpUV, TmpPositions, NextCoordinates, RandomNumber);
	   if (UseExactFlag == false)
	     TmpMetropolis = TestFunction->CalculateFromSpinorVariables(TmpUV);
	   else
	     {
	       if (UseBaseAsWeightFlag == false)
		 {
		   QHEParticleWaveFunctionOperation Operation(ExactSpace, ExactState, &TmpPositions, ExactBasis);
		   Operation.ApplyOperation(Architecture.GetArchitecture());      
		   TmpMetropolis = Operation.GetScalar();
		 }
	       else
		 TmpMetropolis = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
	     }
	   CurrentProbabilities = Norm(TmpMetropolis);
	   if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	     {
	       PreviousProbabilities = CurrentProbabilities;
	       ++Acceptance;
	     }
	   else
	     {
	       TmpUV.Re(NextCoordinates << 1) = PreviousCoordinatesU.Re;
	       TmpUV.Im(NextCoordinates << 1) = PreviousCoordinatesU.Im;
	       TmpUV.Re(1 + (NextCoordinates << 1)) = PreviousCoordinatesV.Re;
	       TmpUV.Im(1 + (NextCoordinates << 1)) = PreviousCoordinatesV.Im;
	       TmpPositions[NextCoordinates << 1] = PreviousCoordinates1;
	       TmpPositions[1 + (NextCoordinates << 1)] = PreviousCoordinates2;
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
      
       for (int i = 0; i < NbrIter; ++i)
	 {
	   Complex PreviousCoordinatesU = TmpUV[NextCoordinates << 1];
	   Complex PreviousCoordinatesV = TmpUV[1 + (NextCoordinates << 1)];
	   double PreviousCoordinates1 = TmpPositions[NextCoordinates << 1];
	   double PreviousCoordinates2 = TmpPositions[1 + (NextCoordinates << 1)];
	   RandomUVOneCoordinate(TmpUV, TmpPositions, NextCoordinates, RandomNumber);
	   Complex TmpMetropolis;
	   if (UseExactFlag == false)
	     TmpMetropolis = TestFunction->CalculateFromSpinorVariables(TmpUV);
	   else
	     {
	       if (UseBaseAsWeightFlag == false)
		 {
		   QHEParticleWaveFunctionOperation Operation(ExactSpace, ExactState, &TmpPositions, ExactBasis);
		   Operation.ApplyOperation(Architecture.GetArchitecture());      
		   TmpMetropolis = Operation.GetScalar();
		 }
	       else
		 TmpMetropolis = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
	     }
	   CurrentProbabilities = Norm(TmpMetropolis);
	   if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	     {
	       PreviousProbabilities = CurrentProbabilities;
	       Tmp = TmpMetropolis;
	       if (UseBaseAsWeightFlag == false)
		 ValueExact = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
	       else
		 {
		   QHEParticleWaveFunctionOperation Operation(ExactSpace, ExactState, &TmpPositions, ExactBasis);
		   Operation.ApplyOperation(Architecture.GetArchitecture());      
		   ValueExact = Operation.GetScalar();
		 }
	       ++Acceptance;
	     }
	   else
	     {
	       TmpUV.Re(NextCoordinates << 1) = PreviousCoordinatesU.Re;
	       TmpUV.Im(NextCoordinates << 1) = PreviousCoordinatesU.Im;
	       TmpUV.Re(1 + (NextCoordinates << 1)) = PreviousCoordinatesV.Re;
	       TmpUV.Im(1 + (NextCoordinates << 1)) = PreviousCoordinatesV.Im;
	       TmpPositions[NextCoordinates << 1] = PreviousCoordinates1;
	       TmpPositions[1 + (NextCoordinates << 1)] = PreviousCoordinates2;
	       CurrentProbabilities = PreviousProbabilities;
	     }
	   TotalProbability += CurrentProbabilities;
	   NextCoordinates = (int) (((double) NbrParticles) * RandomNumber->GetRealRandomNumber());
	   if (NextCoordinates == NbrParticles)
	     --NextCoordinates;
	   Tmp2 = (Tmp.Re * Tmp.Re) + (Tmp.Im * Tmp.Im);
	   Tmp2bis = (ValueExact.Re * ValueExact.Re) + (ValueExact.Im * ValueExact.Im);
	   Tmp3 = (Conj(Tmp) * ValueExact);
	   Tmp2 /= CurrentProbabilities;
	   Tmp3 /= CurrentProbabilities;      
	   Tmp2bis /= CurrentProbabilities;  
	   Overlap += Tmp3;
	   ErrorOverlap.Re += Tmp3.Re * Tmp3.Re;
	   ErrorOverlap.Im += Tmp3.Im * Tmp3.Im;
	   Normalization += Tmp2;
	   ErrorNormalization += Tmp2 * Tmp2;
	   NormalizationExact += Tmp2bis;
	   ErrorNormalizationExact += Tmp2bis * Tmp2bis;
	   if ((i > 0) && ((RecordStep != 0) && ((i % RecordStep) == 0)))
	     {
	       Complex Tmp4 = Overlap / ((double) i);
	       Complex Tmp5 (sqrt( ((ErrorOverlap.Re / ((double) i)) - (Tmp4.Re * Tmp4.Re)) / ((double) i) ),
			     sqrt( ((ErrorOverlap.Im / ((double) i)) - (Tmp4.Im * Tmp4.Im)) / ((double) i) ));
	       double Tmp6 = Normalization  / ((double) i);
	       double Tmp7 = sqrt( ((ErrorNormalization / ((double) i))  -  (Tmp6 * Tmp6)) / ((double) i) );	  
	       double Tmp8 = NormalizationExact  / ((double) i);
	       double Tmp9 = sqrt( ((ErrorNormalizationExact / ((double) i))  -  (Tmp8 * Tmp8)) / ((double) i) );	  
	       Tmp5.Re /= Tmp4.Re;
	       Tmp5.Im /= Tmp4.Im;
	       Tmp5.Re = fabs(Tmp5.Re);
	       Tmp5.Im = fabs(Tmp5.Im);
	       Tmp5.Re += (Tmp7 / Tmp6);
	       Tmp5.Im += (Tmp7 / Tmp6);
	       Tmp5.Re += (Tmp9 / Tmp8);
	       Tmp5.Im += (Tmp9 / Tmp8);
	       Tmp4 /= sqrt(Tmp6 * Tmp8);	  
	       Tmp5.Re *= Tmp4.Re;
	       Tmp5.Im *= Tmp4.Im;
	       RecordedOverlap[RecordIndex] = Tmp4;
	       RecordedOverlapError[RecordIndex] = Tmp5;
	       ++RecordIndex;
	     }
	   if ((i > 0) && ((i % (((SingleIntegerOption*) Manager["display-step"])->GetInteger())) == 0))
	     {
	       cout << " i = " << i << endl;
	       Complex Tmp4 = Overlap / ((double) i);
	       Complex Tmp5 (sqrt( ((ErrorOverlap.Re / ((double) i)) - (Tmp4.Re * Tmp4.Re)) / ((double) i) ),
			     sqrt( ((ErrorOverlap.Im / ((double) i)) - (Tmp4.Im * Tmp4.Im)) / ((double) i) ));
	       double Tmp6 = Normalization  / ((double) i);
	       double Tmp7 = sqrt( ((ErrorNormalization / ((double) i))  -  (Tmp6 * Tmp6)) / ((double) i) );	  
	       double Tmp8 = NormalizationExact  / ((double) i);
	       double Tmp9 = sqrt( ((ErrorNormalizationExact / ((double) i))  -  (Tmp8 * Tmp8)) / ((double) i) );	  
	       
	       if (((BooleanOption*) Manager["show-details"])->GetBoolean() == true)
		 {
		   cout << Tmp4;
		   cout << " +/- " << Tmp5 << endl;
		   cout << Tmp6;
		   cout << " +/- " << Tmp7 << endl;	  
		   cout << Tmp8;
		   cout << " +/- " << Tmp9 << endl;	  
		 }
	       Tmp5.Re /= Tmp4.Re;
	       Tmp5.Im /= Tmp4.Im;
	       Tmp5.Re = fabs(Tmp5.Re);
	       Tmp5.Im = fabs(Tmp5.Im);
	       Tmp5.Re += (Tmp7 / Tmp6);
	       Tmp5.Im += (Tmp7 / Tmp6);
	       Tmp5.Re += (Tmp9 / Tmp8);
	       Tmp5.Im += (Tmp9 / Tmp8);
	       Tmp4 /= sqrt(Tmp6 * Tmp8);	  
	       Tmp5.Re *= Tmp4.Re;
	       Tmp5.Im *= Tmp4.Im;
	       cout << Tmp4 << " +/- " << Tmp5 << endl;
	       cout << "acceptance rate = " << (((double) Acceptance) / (((double) i))) << endl;
	       cout << "-----------------------------------------------" << endl;
	     }
	 } 
       cout << " final results :" << endl;
       Complex Tmp4 = Overlap / ((double) NbrIter);
       Complex Tmp5 (sqrt( ((ErrorOverlap.Re / ((double) NbrIter)) - (Tmp4.Re * Tmp4.Re)) / ((double) NbrIter) ),
		     sqrt( ((ErrorOverlap.Im / ((double) NbrIter)) - (Tmp4.Im * Tmp4.Im)) / ((double) NbrIter) ));
       double Tmp6 = Normalization  / ((double) NbrIter);
       double Tmp7 = sqrt( ((ErrorNormalization / ((double) NbrIter))  -  (Tmp6 * Tmp6)) / ((double) NbrIter) );	  
       double Tmp8 = NormalizationExact  / ((double) NbrIter);
       double Tmp9 = sqrt( ((ErrorNormalizationExact / ((double) NbrIter))  -  (Tmp8 * Tmp8)) / ((double) NbrIter) );	  
       
       Tmp5.Re /= Tmp4.Re;
       Tmp5.Im /= Tmp4.Im;
       Tmp5.Re = fabs(Tmp5.Re);
       Tmp5.Im = fabs(Tmp5.Im);
       Tmp5.Re += (Tmp7 / Tmp6);
       Tmp5.Im += (Tmp7 / Tmp6);
       Tmp5.Re += (Tmp9 / Tmp8);
       Tmp5.Im += (Tmp9 / Tmp8);
       Tmp4 /= sqrt(Tmp6 * Tmp8);	  
       Tmp5.Re *= Tmp4.Re;
       Tmp5.Im *= Tmp4.Im;
       cout << Tmp4 << " +/- " << Tmp5 << endl;
       cout << Norm(Tmp4) << " +/- " << (2.0 * (fabs(Tmp4.Re * Tmp5.Re) + fabs(Tmp4.Im * Tmp5.Im))  / Norm(Tmp4)) << endl;
       cout << "-----------------------------------------------" << endl;
       
       
       if (RecordStep != 0)
	 {
	   RecordedOverlap[RecordIndex] = Tmp4;
	   RecordedOverlapError[RecordIndex] = Tmp5;
	   ofstream OverlapRecordFile;
	   OverlapRecordFile.precision(14);
	   OverlapRecordFile.open(((SingleStringOption*) Manager["record-file"])->GetString(), ios::out | ios::binary);
	   int NbrRecords = NbrIter / RecordStep;
	   OverlapRecordFile << "# Monte Carlo overlap calculation" << endl
			     << "# step overlap.Re overlap.Im overlap2 error.Re error.Im error2" << endl;
	   for (int i = 0; i < NbrRecords; ++i)
	     OverlapRecordFile << i << " " << RecordedOverlap[i].Re << " " << RecordedOverlap[i].Im << " " << SqrNorm(RecordedOverlap[i])
			       << " " << RecordedOverlapError[i].Re << " " <<  RecordedOverlapError[i].Im << " " 
			       << (2 * ((RecordedOverlap[i].Re * RecordedOverlapError[i].Re) + (RecordedOverlap[i].Im * RecordedOverlapError[i].Im))) << endl;
	   OverlapRecordFile.close();
	 }
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
