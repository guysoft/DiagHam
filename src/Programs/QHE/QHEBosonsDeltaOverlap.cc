#include "Vector/RealVector.h"

#include "HilbertSpace/QHEHilbertSpace/BosonOnSphere.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereDeltaHamiltonian.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereDeltaModifiedHamiltonian.h"
#include "Hamiltonian/QHEHamiltonian/ParticleOnSphereCoulombDeltaHamiltonian.h"
#include "FunctionBasis/QHEFunctionBasis/ParticleOnSphereFunctionBasis.h"

#include "Tools/QHE/QHEWaveFunction/LaughlinOnSphereWaveFunction.h"
#include "Tools/QHE/QHEWaveFunction/PfaffianOnSphereWaveFunction.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "MainTask/QHEMainTask/QHEOnSphereMainTask.h"

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
  OptionManager Manager ("QHEBosonsDelta" , "0.01");
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
  (*SystemGroup) += new SingleStringOption  ('\n', "exact-state", "name of the file containing the vector obtained using exact diagonalization");
  (*SystemGroup) += new SingleStringOption  ('\n', "test-wavefunction", "name of the test wave fuction", "laughlin");
  (*SystemGroup) += new BooleanOption ('\n', "list-wavefunctions", "list all available test wave fuctions");  

  (*MonteCarloGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*MonteCarloGroup) += new SingleIntegerOption  ('\n', "display-step", "number of iteration between two consecutive result displays", 1000);
  (*MonteCarloGroup) += new BooleanOption  ('\n', "with-timecoherence", "use time coherence between two successive evaluation of the wave function");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsDeltaOverlap -h" << endl;
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

  int NbrBosons = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int NbrIter = ((SingleIntegerOption*) Manager["nbr-iter"])->GetInteger();

  if (((SingleStringOption*) Manager["exact-state"])->GetString() == 0)
    {
      cout << "QHEBosonsDeltaOverlap requires an exact state" << endl;
      return -1;
    }
  RealVector State;
  if (State.ReadVector (((SingleStringOption*) Manager["exact-state"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["exact-state"])->GetString() << endl;
      return -1;      
    }
  BosonOnSphere Space (NbrBosons, 0, LzMax);
  ParticleOnSphereFunctionBasis Basis(LzMax);
//  Abstract1DComplexFunction* WaveFunction = new LaughlinOnSphereWaveFunction(NbrBosons, 2);
  Abstract1DComplexFunction* WaveFunction = new PfaffianOnSphereWaveFunction(NbrBosons);
  RealVector Location(2 * NbrBosons, true);
  srand48(29457);
/*  for (int k = 0; k < 10; ++k)
    {
      for (int i = 0; i < NbrBosons; ++i)
	{
	  Location[i << 1] = M_PI * drand48();
	  Location[(i << 1) + 1] = 2.0 * M_PI * drand48();
	  cout << Location[i << 1] << " " << Location[(i << 1) + 1] << endl;
	}
      //  Location[4] = Location[0];
      //  Location[5] = Location[1];
      ParticleOnSphereFunctionBasis Basis(LzMax);
      QHEParticleWaveFunctionOperation Operation(&Space, &State, &Location, &Basis);
      Architecture.GetArchitecture()->ExecuteOperation(&Operation);      
      Complex ValueExact (Operation.GetScalar());
      //      Complex ValueExact = Space.EvaluateWaveFunction(State, Location, Basis);
      Complex ValueLaughlin = LaughlinWaveFunction(Location, NbrBosons) * 0.36563112422012;
//      cout << ValueExact  << endl; 
      cout << ValueExact  << " " << ValueLaughlin << " " << (Norm(ValueExact) / Norm(ValueLaughlin)) << endl;        
      cout << "-------------------------------------" << endl;
    }
  return 0;*/
  double Factor = 1.0;
  for (int j = 0; j < NbrBosons; ++j)
    {
      Factor *= 4.0 * M_PI; //2.0 * M_PI * M_PI;
    }
  Complex Overlap;
  Complex ErrorOverlap;
  double Normalization = 0.0;
  double ErrorNormalization = 0.0;
  Complex Tmp;
  Complex Tmp3;
  double Tmp2;
  int NextCoordinates = 0;
  for (int j = 0; j < NbrBosons; ++j)
    {
      Location[j << 1] = acos (1.0- (2.0 * drand48()));
      Location[1 + (j << 1)] = 2.0 * M_PI * drand48();
    }

  for (int i = 0; i < NbrIter; ++i)
    {
      /*      for (int j = 0; j < NbrBosons; ++j)
	{
	  Location[NextCoordinates << 1] = acos (1.0- (2.0 * drand48()));
	  Location[1 + (NextCoordinates << 1)] = 2.0 * M_PI * drand48();
	  }*/
      Location[NextCoordinates << 1] = acos (1.0- (2.0 * drand48()));
      Location[1 + (NextCoordinates << 1)] = 2.0 * M_PI * drand48();
      NextCoordinates = (int) (((double) NbrBosons) * drand48());
      if (NextCoordinates == NbrBosons)
	--NextCoordinates;
      Tmp = (*WaveFunction)(Location);
      int TimeCoherence = NextCoordinates;
      if (((BooleanOption*) Manager["with-timecoherence"])->GetBoolean() == false)
	TimeCoherence = -1;
      QHEParticleWaveFunctionOperation Operation(&Space, &State, &Location, &Basis, TimeCoherence);
      Architecture.GetArchitecture()->ExecuteOperation(&Operation);      
      Complex ValueExact (Operation.GetScalar());
//      ValueExact.Re *= -1.0;
      Tmp2 = (Tmp.Re * Tmp.Re) + (Tmp.Im * Tmp.Im);
      Tmp3 = (Conj(Tmp) * ValueExact);
      Overlap += Tmp3;// * Factor;
      ErrorOverlap.Re += Tmp3.Re * Tmp3.Re;//  * Factor * Factor;
      ErrorOverlap.Im += Tmp3.Im * Tmp3.Im;//  * Factor * Factor;
//      ErrorOverlap += Tmp3 * Tmp3  * Factor * Factor;
      Normalization += Tmp2;// * Factor;
      ErrorNormalization += Tmp2 * Tmp2;// *  Factor * Factor;
      if ((i > 0) && ((i % (((SingleIntegerOption*) Manager["display-step"])->GetInteger())) == 0))
	{
	  cout << " i = " << i << endl;
	  Complex Tmp4 = Overlap / ((double) i);
	  cout << (Tmp4 * Factor);
//	  Complex Tmp5 (sqrt((((ErrorOverlap.Re / ((double) (i))) - (Overlap.Re * Overlap.Re)) / ((double) (i))) / Factor),
//			sqrt((((ErrorOverlap.Im / ((double) (i))) - (Overlap.Im * Overlap.Im)) / ((double) (i))) / Factor));
	  Complex Tmp5 (sqrt( ((ErrorOverlap.Re / ((double) i)) - (Tmp4.Re * Tmp4.Re)) / ((double) i) ),
			sqrt( ((ErrorOverlap.Im / ((double) i)) - (Tmp4.Im * Tmp4.Im)) / ((double) i) ));
	  cout << " +/- " << (Tmp5 * Factor) << endl;
	  double Tmp6 = Normalization / ((double) i);
	  cout << Factor * Tmp6;
//	  Tmp6 = sqrt((((ErrorNormalization / ((double) (i))) - (Tmp6 * Tmp6)) / ((double) (i))) / Factor);
	  double Tmp7 = sqrt( ((ErrorNormalization / ((double) i))  -  (Tmp6 * Tmp6)) / ((double) i) );	  
	  cout << " +/- " << (Tmp7  * Factor) << endl;	  
	  Tmp5.Re /= Tmp4.Re;
	  Tmp5.Im /= Tmp4.Im;
	  Tmp5.Re = fabs(Tmp5.Re);
	  Tmp5.Im = fabs(Tmp5.Im);
	  Tmp5.Re += (Tmp7 / Tmp6);
	  Tmp5.Im += (Tmp7 / Tmp6);
	  Tmp4 *= sqrt(Factor / Tmp6);	  
	  Tmp5.Re *= Tmp4.Re;
	  Tmp5.Im *= Tmp4.Im;
	  cout << Tmp4 << " " << Tmp5 << endl;
	  cout << "-----------------------------------------------" << endl;
	}
    } 
  cout << " final results :" << endl;
  Complex Tmp4 = Overlap / ((double) NbrIter);
  cout << (Tmp4 * Factor);
  Complex Tmp5 (sqrt( ((ErrorOverlap.Re / ((double) NbrIter)) - (Tmp4.Re * Tmp4.Re)) / ((double) NbrIter) ),
		sqrt( ((ErrorOverlap.Im / ((double) NbrIter)) - (Tmp4.Im * Tmp4.Im)) / ((double) NbrIter) ));
  cout << " +/- " << (Tmp5 * Factor) << endl;
  double Tmp6 = Normalization / ((double) NbrIter);
  cout << Factor * Tmp6;
  double Tmp7 = sqrt( ((ErrorNormalization / ((double) NbrIter))  -  (Tmp6 * Tmp6)) / ((double) NbrIter) );	  
  cout << " +/- " << (Tmp7  * Factor) << endl;	  
  Tmp5.Re /= Tmp4.Re;
  Tmp5.Im /= Tmp4.Im;
  Tmp5.Re = fabs(Tmp5.Re);
  Tmp5.Im = fabs(Tmp5.Im);
  Tmp5.Re += (Tmp7 / Tmp6);
  Tmp5.Im += (Tmp7 / Tmp6);
  Tmp4 *= sqrt(Factor / Tmp6);	  
  Tmp5.Re *= Tmp4.Re;
  Tmp5.Im *= Tmp4.Im;
  cout << Tmp4 << " " << Tmp5 << endl;
  cout << "-----------------------------------------------" << endl;
 return 0;
}


