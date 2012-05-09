#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereLong.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "GeneralTools/ConfigurationParser.h"

#include "Operator/ParticleOnSphereDensityOperator.h"
#include "FunctionBasis/ParticleOnCylinderFunctionBasis.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"

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
  OptionManager Manager ("FQHECylinderCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('e', "eigenstate", "name of the file containing the eigenstate");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "ky-max", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "total-y", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "landau-level", "index of the Landau level (0 being the LLL)", 0);
  (*SystemGroup) += new SingleDoubleOption  ('r', "ratio", "aspect ratio of the cylinder", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 1000);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with rhorho extension");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsCorrelation -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int KyMax = Manager.GetInteger("ky-max");
  int TotalKy = Manager.GetInteger("total-y");
  int LandauLevel = Manager.GetInteger("landau-level");
  int NbrPoints = Manager.GetInteger("nbr-points");
  double Ratio = Manager.GetDouble("ratio");
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  bool Statistics = true;
  if (Manager.GetString("eigenstate") == 0)
    {
      cout << "FQHECylinderFermionsCorrelation requires a state" << endl;
      return -1;
    }

  //if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("eigenstate"),
//						  NbrParticles, KyMax, TotalKy, Statistics) == false)
//    {
//      cout << "error while retrieving system parameters from file name " << Manager.GetString("eigenstate") << endl;
//      return -1;
//    }

  if (IsFile(Manager.GetString("eigenstate")) == false)
    {
      cout << "can't find vector file " << Manager.GetString("eigenstate") << endl;
      return -1;      
    }

  ParticleOnSphere* Space = 0;
#ifdef __64_BITS__
	  if (TotalKy <= 62)
#else
	  if (TotalKy <= 30)
#endif
  	    Space = new FermionOnSphere(NbrParticles, TotalKy, KyMax);

	  else
#ifdef __128_BIT_LONGLONG__
	    if (TotalKy <= 126)
#else
	      if (TotalKy <= 62)
#endif
 	        Space = new FermionOnSphereLong(NbrParticles, TotalKy, KyMax);
	      else
		Space = new FermionOnSphereUnlimited(NbrParticles, TotalKy, KyMax);


  cout << " Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;

  ofstream File;
  File.precision(14);
  if (Manager.GetString("output-file") != 0)
     File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
   {
    cout << "Enter output file! " << endl;
    exit(1);
   }
  
  ComplexVector State;

  if (State.ReadVector (Manager.GetString("eigenstate")) == false)
   {
     cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
     return -1;      
   }

  cout << "NbrParticles "<<NbrParticles<<" "<<KyMax<<" "<<LandauLevel<<" "<<Ratio<<endl;
  cout << "NbrPoints "<<NbrPoints<<endl;
  
  ParticleOnCylinderFunctionBasis Basis (KyMax, LandauLevel, Ratio);

  double H = sqrt(2.0 * M_PI * (KyMax + 1.0))/sqrt(Ratio);
  double Step = (H/(double)NbrPoints);
  Complex DensityIntegral (0.0,0.0);

  Complex* Occupations= new Complex [KyMax+1];
      
  cout<<"Orbital occupations: ";
  Complex CheckOccupations(0,0);
  for (int i = 0; i <= KyMax; ++i)
   {
     ParticleOnSphereDensityOperator Operator (Space, i);
     Occupations[i] = Operator.MatrixElement(State, State);
     CheckOccupations += Occupations[i];
     cout << Occupations[i] <<" ";
   }
  cout<<endl;
  cout<<"Sum of all occupations: "<<CheckOccupations<<endl;

  double XPosition = -0.5*H;
  while (XPosition <= 0.5*H)
    {
      Complex Density (0.0, 0.0);        
      for (int i = 0; i <= KyMax; ++i)
        {
          ParticleOnSphereDensityOperator Operator (Space, i);
          Density += Conj(Basis.GetFunctionValue(XPosition, 0.0, (double)i-0.5*KyMax)) * Basis.GetFunctionValue(XPosition, 0.0, (double)i-0.5*KyMax) * Occupations[i]/((double)NbrParticles);
	  //cout<<"i= "<<i<<" "<<DensityMatEl[i]<<" ";
        }

      DensityIntegral += Density * Step * Ratio * H;
      File << XPosition <<" "<< Density.Re << " " << Density.Im << endl;
      //cout << XPosition <<" "<< Density << endl;
      XPosition += Step;
    }

  cout << "Integrated density along the cylinder: "<<DensityIntegral <<endl;

  delete[] Occupations;

  File.close();
 
  return 0;
}
