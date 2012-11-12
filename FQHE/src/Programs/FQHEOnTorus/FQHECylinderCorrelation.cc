#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereLong.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "GeneralTools/ConfigurationParser.h"

#include "Operator/ParticleOnSphereDensityOperator.h"
#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "FunctionBasis/ParticleOnCylinderFunctionBasis.h"
#include "Hamiltonian/ParticleOnCylinderStructureFactor.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

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
  (*SystemGroup) += new BooleanOption  ('\n', "rho-rho","evaluate rho-rho correlation", false);
  (*SystemGroup) += new BooleanOption  ('\n', "structure-factor","evaluate the guiding center structure factor", false);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "x-points", "number of points along the cylinder", 50);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
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
  int XPoints = Manager.GetInteger("x-points");
  double Ratio = Manager.GetDouble("ratio");
  bool EvaluateS0Q = Manager.GetBoolean("structure-factor");
  bool EvaluateRhoRho = Manager.GetBoolean("rho-rho");

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

  if (EvaluateS0Q == true)
    {
      cout << "Evaluating structure factor...."<<endl;
      double Length = sqrt(2.0 * M_PI * Ratio * (KyMax + 1));
      double kappaL = 2.0 * M_PI/Length;
      double Height = Length/Ratio;
      double kappaH = 2.0 * M_PI/Height;
      long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;

      int counter = 1;
      for (int s = 0; s <= KyMax; s++)
        for (int t = 0; t <= KyMax; t++)
          {
            cout << "Doing sector " << counter << " out of " << (KyMax+1)*(KyMax+1) << endl;
            counter++;

            int mint = t;
            double QyValue = kappaL * t;
            double QyValueRescaled = kappaL * (t - (KyMax + 1));
            if (QyValue*QyValue > QyValueRescaled*QyValueRescaled)
              {
               mint = t - (KyMax+1);
               QyValue = QyValueRescaled;
              }
   

            int mins = s;
            double QxValue = kappaH * s;
            double QxValueRescaled = kappaH * (s - (KyMax + 1));
            if (QxValue*QxValue > QxValueRescaled*QxValueRescaled)
              {
               mins = s - (KyMax+1);
               QxValue = QxValueRescaled;
              }
   
            double Q = sqrt(QxValue * QxValue + QyValue * QyValue);
       
            AbstractHamiltonian* Hamiltonian = new ParticleOnCylinderStructureFactor (Space, NbrParticles, KyMax, Ratio, QxValue, mint, Architecture.GetArchitecture(), Memory);

           ComplexVector TmpState(Space->GetHilbertSpaceDimension(), true);
           VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
           Operation.ApplyOperation(Architecture.GetArchitecture());
           Complex S0Q = State*TmpState;

           //Add the term \sum_i <exp(iq.R_i)exp(-iq.R_i)>=N_e 
           S0Q += (double)NbrParticles/(double)(KyMax+1);

           if (Q == 0.0)
           {  
            //Subtract the term 1/N_phi \sum_i <exp(iq.R_i)> \sum_j <exp(-iq.R_j)>
            S0Q -= (double)(NbrParticles*NbrParticles)/(double)(KyMax+1);
           }

           delete Hamiltonian;

           cout << "Qx= "<<QxValue<<" , Qy= "<<QyValue<<" "<<S0Q<<endl; 
           File << Q << " " << QxValue << " " << QyValue <<" "<<S0Q.Re<<" "<<S0Q.Im<<endl; 
         }  
	
      return 0;
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

  if (EvaluateRhoRho)
   {
     Complex TmpValue;
     RealVector Value(2, true);
     Complex* PrecalculatedValues = new Complex [KyMax + 1];	  

      cout<<"density-density precalculate ";
      for (int i = 0; i <= KyMax; ++i)
	{
             TmpValue = Basis.GetFunctionValue(-0.5*H, 0.0, -0.5*KyMax);
	    ParticleOnSphereDensityDensityOperator Operator (Space, i, 0, i, 0);
	    PrecalculatedValues[i] = Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
             cout<< i <<" " << PrecalculatedValues[i] << " "<<TmpValue * Conj(TmpValue)<<endl;
	}
      cout<<"done."<<endl;
  
    File << "# pair correlation coefficients " << endl;
    File << "#" << endl << "# (l+S)    n_l" << endl;
      for (int i = 0; i <= KyMax; ++i)
	File << "# " << i << " " << PrecalculatedValues[i].Re << endl;
    double XInc = (H + 2.0)/ ((double) NbrPoints);

        for (int k = 0; k <= NbrPoints; ++k)
	  {
            double X = -0.5*(H + 2.0) + (double)k * XInc;
	    Complex Sum (0.0, 0.0);
	    for (int i = 0; i <= KyMax; ++i)
	      {
	        Complex TmpValue = Basis.GetFunctionValue(X, 0.0, (double)i - 0.5 * KyMax);
	        Sum += PrecalculatedValues[i] * (Conj(TmpValue) * TmpValue);
	      }
            File << X << " " << Norm(Sum) << endl;
	  }
       

   File.close();
 
  delete[] PrecalculatedValues;

  return 0;

  }


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
