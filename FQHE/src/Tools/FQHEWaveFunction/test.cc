#include "JainCFOnSphereOrbitals.h"
#include "Vector/RealVector.h"
#include "MathTools/RandomNumber/RanluxRandomNumberGenerator.h"
#include "CFOrbitals.h"
#include "particle.h"
#include <sys/time.h>
#include <iostream>
#include <complex>
using std::cout;
using std::endl;
using std::complex;

complex <double> mixedDiv(complex<double> std, Complex DiagHam)
{
  complex <double> toStd(Real(DiagHam),Imag(DiagHam));
  return std/toStd;
}


int main()
{
  unsigned int seed=3;
  AbstractRandomNumberGenerator *R= new RanluxRandomNumberGenerator(seed);
  JainCFOnSphereOrbitals test(/* nbrParticles */ 8, /* nbrLandauLevels */ 5, /*nbrEffectiveFlux*/ 5, /* jastrowPower */ 2);
  CFOrbitals oldOrbitals(/* N1 */ 8, /* N_eff */ 5, /* p */ 1, /* highestLL */ 5 , CFOrbitals::increasing, /* theta_c */ 0.5);
  RealVector x(2*test.GetNbrParticles());
  for (int i=0; i<test.GetNbrParticles(); i++)
    {
      x[2*i]=M_PI*R->GetRealRandomNumber();
      x[2*i+1]=2*M_PI*R->GetRealRandomNumber();
    }
  ComplexMatrix phiNEW(test.GetNbrParticles(),test.GetNbrOrbitals());
  phiNEW = test(x);

  particle * particles = oldOrbitals.getParticles();
  for (int i=0; i<test.GetNbrParticles(); i++)
    particles[i].set_position(x[2*i],x[2*i+1]);
  oldOrbitals.recalculateDistances();
  cout << "oldOrbitals:"<< endl;
  for (int i=0; i<test.GetNbrParticles(); i++)
    cout << particles[i].u << ", " << particles[i].v<<endl;
  oldOrbitals.printDistances();
  oldOrbitals.calcOrbitals();
  complex<double> **phi=oldOrbitals.getOrbitals();
  cout << "New orbitals:" << endl;
  cout<<phiNEW;
  cout << "Old orbitals:" << endl;
  for (int i=0; i<test.GetNbrParticles(); i++)
    {
      for (int j=0; j<oldOrbitals.getNbrOrbitals(); ++j)
	cout << phi[i][j] << " ";
      cout << endl;
    }
  cout << "Comparing ratios:" << endl;
  for (int j=0; j<oldOrbitals.getNbrOrbitals(); ++j)
    {
      cout << "Orbital " <<j<<": ";
      for (int i=0; i<test.GetNbrParticles(); i++)
	cout << mixedDiv(phi[i][j],phiNEW[j][i]) << ", ";
      cout << endl;
    }
      
//   timeval TotalStartingTime;
//   timeval TotalEndingTime;
//   double Dt2;
//   gettimeofday (&(TotalStartingTime), 0);
//   int numEval=1000;
//   for (int i=0; i< numEval; ++i)
//     test(x);
//   gettimeofday (&(TotalEndingTime), 0);
//   Dt2 = (double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
//     ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0);
//   cout << "Number of Orbitals = " << test.GetNbrOrbitals() << endl;
//   cout << "time per evaluation = " << Dt2/numEval << endl;
  delete  R;
  return 0;
}
