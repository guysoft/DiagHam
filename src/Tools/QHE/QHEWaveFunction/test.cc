#include "JainCFOnSphereOrbitals.h"
#include "Vector/RealVector.h"
#include "MathTools/RandomNumber/RanluxRandomNumberGenerator.h"

#include <sys/time.h>
#include <iostream>
using std::cout;
using std::endl;

int main()
{
  unsigned int seed=3;
  AbstractRandomNumberGenerator *R= new RanluxRandomNumberGenerator(seed);
  JainCFOnSphereOrbitals test(/* nbrParticles */ 20, /* nbrLandauLevels */ 3, /*nbrEffectiveFlux*/ -2, /* jastrowPower */ 2);
  RealVector x(2*test.GetNbrParticles());
  for (int i=0; i<test.GetNbrParticles(); i++)
    {
      x[2*i]=M_PI*R->GetRealRandomNumber();
      x[2*i+1]=2*M_PI*R->GetRealRandomNumber();
    }
  cout<<test(x);
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  double Dt2;
  gettimeofday (&(TotalStartingTime), 0);
  int numEval=1000;
  for (int i=0; i< numEval; ++i)
    test(x);
  gettimeofday (&(TotalEndingTime), 0);
  Dt2 = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
  cout << "Number of Orbitals = " << test.GetNbrOrbitals() << endl;
  cout << "time per evaluation = " << Dt2/numEval << endl;
  return 0;
}
