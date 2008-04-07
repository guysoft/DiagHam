#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"

#include "Options/Options.h"

#include "Tools/FQHEWaveFunction/FQHESphereSymmetrizedSUKToU1WaveFunction.h"
#include "Tools/FQHEWaveFunction/FQHESU3HalperinPermanentOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/HalperinOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/JainCFFilledLevelOnSphereWaveFunction.h"
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


void RandomUV (ComplexVector& uv, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator);

void RandomUVOneCoordinate(ComplexVector& uv, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator);

void FlipCoordinates (ComplexVector& uv, int i, int j);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereSUKToU1MCOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MonteCarloGroup = new OptionGroup ("Monte Carlo options");

  Manager += SystemGroup;
  Manager += MonteCarloGroup;
  Manager += MiscGroup;


  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 10);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 9);
  (*SystemGroup) += new SingleIntegerOption  ('k', "k-value", "k index of the SU(k) symmetry group", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "intra-corr", "power of the intra-color correlations", 3);  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "inter-corr", "power of the inter-color correlations", 2);  
  (*SystemGroup) += new SingleStringOption  ('\n', "exact-state", "name of the file containing the vector obtained using exact diagonalization");
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  
  (*SystemGroup) += new BooleanOption ('\n', "overlap", "list all available test wave fuctions");  
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
  bool OverlapFlag = ((BooleanOption*) Manager["overlap"])->GetBoolean();

  Abstract1DComplexFunctionOnSphere* BaseFunction = 0;
  switch (KValue)
    {
    case 2:
      {
	BaseFunction = new HalperinOnSphereWaveFunction (NbrParticlePerColor, NbrParticlePerColor, 
							 IntraCorrelation - 1, IntraCorrelation - 1, InterCorrelation - 1);
      }
      break;
    case 3:
      {
// 	BaseFunction = new SU3HalperinOnSphereWaveFunction (NbrParticlePerColor, NbrParticlePerColor, NbrParticlePerColor, 
// 							    IntraCorrelation - 1, IntraCorrelation - 1, IntraCorrelation - 1, 
// 							    InterCorrelation - 1, InterCorrelation - 1, InterCorrelation - 1);
	BaseFunction = new FQHESU3HalperinPermanentOnSphereWaveFunction(NbrParticlePerColor, NbrParticlePerColor, NbrParticlePerColor, 
									IntraCorrelation - 1, IntraCorrelation - 1, IntraCorrelation - 1, 
									InterCorrelation - 1, InterCorrelation - 1, InterCorrelation - 1);
      }
      break;
    default:
      {
	cout << "invalid or unsupported number of colors (i.e. k)" << endl;
	return -1;
      }
    }
  FQHESphereSymmetrizedSUKToU1WaveFunction* SymmetrizedFunction = new FQHESphereSymmetrizedSUKToU1WaveFunction (NbrParticles, KValue, BaseFunction, true);
  PfaffianOnSphereWaveFunction* ReferenceFunction = new PfaffianOnSphereWaveFunction(NbrParticles);
  JainCFFilledLevelOnSphereWaveFunction* CFFunction = new JainCFFilledLevelOnSphereWaveFunction(NbrParticles, KValue, 2);

   StdlibRandomNumberGenerator RandomNumberGenerator(29457);

   if (OverlapFlag == true)
     {
       int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();
       AbstractRandomNumberGenerator* RandomNumber = new StdlibRandomNumberGenerator (29457);
  
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
       RandomUV (TmpUV, NbrParticles, RandomNumber);
       Tmp = CFFunction->CalculateFromSpinorVariables(TmpUV);
       double PreviousProbabilities = Norm(Tmp);
       double CurrentProbabilities = PreviousProbabilities;
       double TotalProbability = PreviousProbabilities;
       for (int i = 0; i < NbrIter; ++i)
	 {
	   Complex PreviousCoordinatesU = TmpUV[NextCoordinates << 1];
	   Complex PreviousCoordinatesV = TmpUV[1 + (NextCoordinates << 1)];
	   RandomUVOneCoordinate(TmpUV, NextCoordinates, RandomNumber);
	   Complex TmpMetropolis = CFFunction->CalculateFromSpinorVariables(TmpUV);
	   CurrentProbabilities = Norm(TmpMetropolis);
	   if ((CurrentProbabilities > PreviousProbabilities) || ((RandomNumber->GetRealRandomNumber() * PreviousProbabilities) < CurrentProbabilities))
	     {
	       PreviousProbabilities = CurrentProbabilities;
	       Tmp = TmpMetropolis;
	     }
	   else
	     {
	       TmpUV[NextCoordinates << 1] = PreviousCoordinatesU;
	       TmpUV[1 + (NextCoordinates << 1)] = PreviousCoordinatesV;
	       CurrentProbabilities = PreviousProbabilities;
	     }
	   TotalProbability += CurrentProbabilities;
	   NextCoordinates = (int) (((double) NbrParticles) * RandomNumber->GetRealRandomNumber());
	   if (NextCoordinates == NbrParticles)
	     --NextCoordinates;
	   
	   Complex ValueExact = CFFunction->CalculateFromSpinorVariables(TmpUV);//SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
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
	       cout << "-----------------------------------------------" << endl;
	     }
	 } 
       cout << " final results :" << endl;
       Complex Tmp4 = Overlap / ((double) NbrIter);
       cout << Tmp4;
       Complex Tmp5 (sqrt( ((ErrorOverlap.Re / ((double) NbrIter)) - (Tmp4.Re * Tmp4.Re)) / ((double) NbrIter) ),
		     sqrt( ((ErrorOverlap.Im / ((double) NbrIter)) - (Tmp4.Im * Tmp4.Im)) / ((double) NbrIter) ));
       cout << " +/- " << Tmp5 << endl;
       double Tmp6 = Normalization  / ((double) NbrIter);
       cout << Tmp6;
       double Tmp7 = sqrt( ((ErrorNormalization / ((double) NbrIter))  -  (Tmp6 * Tmp6)) / ((double) NbrIter) );	  
       cout << " +/- " << Tmp7 << endl;	  
       double Tmp8 = NormalizationExact  / ((double) NbrIter);
       cout << Tmp8;
       double Tmp9 = sqrt( ((ErrorNormalizationExact / ((double) NbrIter))  -  (Tmp8 * Tmp8)) / ((double) NbrIter) );	  
       cout << " +/- " << Tmp9 << endl;	  
       
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
//       Complex Sum = 0.0;
//        double RefNorm = 0.0;
//        double TestNorm = 0.0;
//        for (int i = 0; i < 100000; ++i)
// 	 {
// 	   ComplexVector TmpUV (NbrParticles * 2, true);
// 	   RandomUV (TmpUV, NbrParticles, &RandomNumberGenerator);      
// 	   Complex Tmp1 = SymmetrizedFunction->CalculateFromSpinorVariables(TmpUV);
// 	   Complex Tmp2 = CFFunction->CalculateFromSpinorVariables(TmpUV);
// 	   Sum += (Conj(Tmp2) * Tmp1);
// 	   RefNorm +=  SqrNorm(Tmp2);
// 	   TestNorm +=  SqrNorm(Tmp1);
// 	   if ((i % 100) == 0)
// 	     cout << (i / 100) << " " << (SqrNorm(Sum) / (RefNorm * TestNorm)) << endl;
// 	 }
     }
   else
     {
       ComplexVector UV (NbrParticles * 2, true);
       RandomUV (UV, NbrParticles, &RandomNumberGenerator);
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

void RandomUV (ComplexVector& uv, int nbrParticles, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  for (int j = 0; j < nbrParticles; ++j)
    {
      double x = acos (1.0 - (2.0 * randomNumberGenerator->GetRealRandomNumber()));
      double y = 2.0 * M_PI * randomNumberGenerator->GetRealRandomNumber();
      uv.Re(2 * j) = cos(0.5 * x);
      uv.Im(2 * j) = uv.Re(2 * j) * sin(0.5 * y);
      uv.Re(2 * j) *= cos(0.5 * y);
      uv.Re(2 * j + 1) = sin(0.5 * x);
      uv.Im(2 * j + 1) = - uv.Re(2 * j + 1) * sin(0.5 * y);
      uv.Re(2 * j + 1) *= cos(0.5 * y);      
    }
}

void RandomUVOneCoordinate(ComplexVector& uv, int coordinate, AbstractRandomNumberGenerator* randomNumberGenerator)
{
  coordinate *= 2;
  double x = acos (1.0 - (2.0 * randomNumberGenerator->GetRealRandomNumber()));
  double y = 2.0 * M_PI * randomNumberGenerator->GetRealRandomNumber();
  uv.Re(coordinate) = cos(0.5 * x);
  uv.Im(coordinate) = uv.Re(coordinate) * sin(0.5 * y);
  uv.Re(coordinate) *= cos(0.5 * y);
  ++coordinate;
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
