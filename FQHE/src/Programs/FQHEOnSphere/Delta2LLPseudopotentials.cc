#include "Options/Options.h"

#include "Tools/FQHESpectrum/PseudoPotentials.h"
#include "Tools/FQHESpectrum/AbstractZDensityProfile.h"

#include "Vector/RealVector.h"
 
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include "GeneralTools/FilenameTools.h"

#include "MathTools/FactorialCoefficient.h"
#include "MathTools/ClebschGordanCoefficients.h"

#define		my_min(X, Y)  ((X) < (Y) ? (X) : (Y))

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <stdlib.h>
#include <math.h>

using std::ofstream;
using std::ios;
using std::cout;
using std::endl;
using std::string;

double CalculateBetaFunction(long x, long y);
double CalculateNormalization(double Q,double l,double m);
double CalculateDeltaInteractionFactor(double Q,double l1,double m1, double l2, double m2, double l3, double m3, double l4, double m4);


int main(int argc, char** argv)
{
  OptionManager Manager ("Delta2LLPseudopotentials" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (i.e. twice the maximum momentum for a single particle on LLL)", 8);
  
 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type Delta2LLPseudopotentials -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFluxQuanta = Manager.GetInteger("nbr-flux");
  int LzMaxUp = NbrFluxQuanta + 2;
  int LzMaxDown = NbrFluxQuanta + 1;
  char* OutputFile = Manager.GetFormattedString("pseudopotential_delta_s_%nbr-flux%.dat");
  
    
  ClebschGordanCoefficients ClebschDownDown (LzMaxDown - 1, LzMaxDown - 1);
  ClebschGordanCoefficients ClebschUpUp (LzMaxUp, LzMaxUp);
  ClebschGordanCoefficients ClebschUpDown (LzMaxUp, LzMaxDown - 1);
  ClebschGordanCoefficients ClebshDownUp (LzMaxDown - 1, LzMaxUp);
  
  // these are the labels of the arrays as they will be in the file.
  string PseudoLabels[9] = {"PseudopotentialsUpUpUpUp","PseudopotentialsUpUpDownDown","PseudopotentialsUpUpUpDown",
			    "PseudopotentialsDownDownUpUp","PseudopotentialsDownDownDownDown","PseudopotentialsDownDownUpDown",
			    "PseudopotentialsUpDownUpUp","PseudopotentialsUpDownDownDown","PseudopotentialsUpDownUpDown"};
  
  // these are the lenghts of the arrays corresponding to the labels above. 			    
  int PseudoLengths[9] = { LzMaxUp + 1, LzMaxUp - 1 , LzMaxUp - 1, LzMaxUp - 1, LzMaxUp - 1, LzMaxUp - 2, LzMaxUp - 1, LzMaxUp - 2, LzMaxUp - 1}; 
  
  // these are the lenghts of the arrays corresponding to the labels above. 			    
  int PseudoMins[9] = { 0, 0, 1, 0, 0, 1, 1, 1, 1}; 
  
  int LLs[9][4] = {{1,1,1,1}, {1,1,0,0}, {1,1,1,0}, {0,0,1,1}, {0,0,0,0}, {0,0,1,0}, {1,0,1,1}, {1,0,0,0}, {1,0,1,0}}; 
  
  
  double *Pseudopotentials;
  Pseudopotentials = new double[LzMaxUp+1];
  ofstream File;
  File.open(OutputFile, ios::binary | ios::out);
  File.precision(14);
  
  double Q = (double)NbrFluxQuanta/2.0;
  
  ClebschGordanCoefficients *LeftCG = &ClebschUpUp;
  ClebschGordanCoefficients *RightCG = &ClebschUpUp;
  for ( int j = 0; j < 9 ; j++ ) 
    {
    File << PseudoLabels[j]<< "=";        
    for ( int L = PseudoLengths[j] - 1 + PseudoMins[j] ; L >= PseudoMins[j] ; L--)
      {
	int idx = PseudoLengths[j] - 1 + PseudoMins[j] - L;	  
	Pseudopotentials[idx] = 0.0;
	int l1 = LLs[j][0], l2 = LLs[j][1], l3 = LLs[j][2], l4 = LLs[j][3];
	cout << l1 << ", " << l2 << ", " << l3 << ", " << l4 << endl;
	if ( l1 == 1 && l2 == 1 ) LeftCG = &ClebschUpUp;
	if ( l1 == 0 && l2 == 0 ) LeftCG = &ClebschDownDown;
	if ( l1 == 1 && l2 == 0 ) LeftCG = &ClebschUpDown;
	if ( l3 == 1 && l4 == 1 ) RightCG = &ClebschUpUp;
	if ( l3 == 0 && l4 == 0 ) RightCG = &ClebschDownDown;
	if ( l3 == 1 && l4 == 0 ) RightCG = &ClebschUpDown;
		
	for ( double m1 = - Q - my_min(l1,l2) ; m1 <= Q + my_min(l1,l2) ; m1+=1.0  ) 
	  {
	    for ( double m2 = - Q - my_min(l3,l4) ; m2 <= Q + my_min(l3,l4) ; m2+=1.0  ) 
	      {
		  //cout << "LeftCG: " << (int)(m1*2) << ", " << -(int)(m1*2) << ", " << L * 2 << ": " << ClebschDownDown.GetCoefficient((int)(m1*2),-(int)(m1*2),L * 2) << endl;
		  //cout << "L: " << L << ", m1: " << m1 << ", m2: " << m2 << ", LCG: " << LeftCG->GetCoefficient((int)(m1*2),-(int)(m1*2),L * 2) << ", RCG: " << RightCG->GetCoefficient((int)(m2*2),-(int)(m2*2),L * 2) << endl;
		  //cout << "Interactionfactor: " << CalculateDeltaInteractionFactor(Q, l1, m1, l2, -m1, l3, m2, l4, -m2) << endl;
		  Pseudopotentials[idx] += LeftCG->GetCoefficient((int)(m1*2),-(int)(m1*2),L * 2)  * RightCG->GetCoefficient((int)(m2*2),-(int)(m2*2),L * 2) 
		  * CalculateDeltaInteractionFactor(Q, l1, m1, l2, -m1, l3, m2, l4, -m2);
	      }
	  }		
      } 
      for (int i = 0; i < PseudoLengths[j] ; ++i)
	  File << " " << Pseudopotentials[i];	
      File << endl ;
  }
        
  File << endl;      
  File.close();

    
  delete[] OutputFile;
  delete[] Pseudopotentials;
    

  return 0;
}


// evaluate a particular interaction factor <l1,m1,l2,m2|\delta|l3,m3,l4.m4>
//
// Q = max angular momentum on LLL.
// l1 = LL of m1, 0 or 1 supported.
// m1 = angular momentum of first operator.
// l2 = LL of m2, 0 or 1 supported.
// m2 = angular momentum of second operator.
// l3 = LL of m3, 0 or 1 supported.
// m3 = angular momentum of third operator.
// l4 = LL of m4, 0 or 1 supported.
// m4 = angular momentum of fourth operator.
//
// return value = interaction factor with delta interaction

double  CalculateDeltaInteractionFactor(double Q,double l1,double m1, double l2, double m2, double l3, double m3, double l4, double m4)
{
   /*cout << "Q: " << Q << ", l1: " << l1 << ", m1:" << m1 ;
   cout << ", l2: " << l2 << ", m2:" << m2 ;
   cout << ", l3: " << l3 << ", m3:" << m3 ;
   cout << ", l4: " << l4 << ", m4:" << m4  << endl;*/
   
  
   double NormCoeff = 1.0;
   double Value = 0.0;
   
   //work out normalisation.
   NormCoeff *= CalculateNormalization(Q, Q+l1, m1)*pow(-1.0,(double)l1);
   NormCoeff *= CalculateNormalization(Q, Q+l2, m2)*pow(-1.0,(double)l2);
   NormCoeff *= CalculateNormalization(Q, Q+l3, m3)*pow(-1.0,(double)l3);
   NormCoeff *= CalculateNormalization(Q, Q+l4, m4)*pow(-1.0,(double)l4);
    
   NormCoeff *= 4.0*M_PI;
   
   for ( int S1 = 0 ; S1 <= l1 ; S1++ ) 
    {
      double S1Coeff = 0.0;
      if ( (l1 + Q - m1 - S1) >= 0 && (Q + m1 + S1) >= 0 ) 
	{
	  FactorialCoefficient* S1Factorial = new FactorialCoefficient(); 
	  S1Factorial->FactorialMultiply((long)(l1+2*Q));
	  S1Factorial->FactorialDivide((long)(l1 + Q - m1 - S1));
	  S1Factorial->FactorialDivide((long)(Q + m1 + S1));
	  S1Coeff = pow(-1.0,(double)S1)*S1Factorial->GetNumericalValue();
	}
	  
      for ( int S2 = 0 ; S2 <= l2 ; S2++ )
        {
	  double S2Coeff = 0.0;
	  if ( (l2 + Q - m2 - S2) >= 0 && (Q + m2 + S2) >= 0 ) 
	    {
	      FactorialCoefficient* S2Factorial = new FactorialCoefficient(); 
	      S2Factorial->FactorialMultiply((long)(l2+2*Q));
	      S2Factorial->FactorialDivide((long)(l2 + Q - m2 - S2));
	      S2Factorial->FactorialDivide((long)(Q + m2 + S2));
	      S2Coeff = pow(-1.0,(double)S2)*S2Factorial->GetNumericalValue();
	    }
          for ( int S3 = 0 ; S3 <= l3 ; S3++ )
            {
	      double S3Coeff = 0.0;
              if ( (l3 + Q - m3 - S3) >= 0 && (Q + m3 + S3) >= 0 ) 
		{
		  FactorialCoefficient* S3Factorial = new FactorialCoefficient(); 
		  S3Factorial->FactorialMultiply((long)(l3+2*Q));
		  S3Factorial->FactorialDivide((long)(l3 + Q - m3 - S3));
		  S3Factorial->FactorialDivide((long)(Q + m3 + S3));
		  S3Coeff = pow(-1.0,(double)S3)*S3Factorial->GetNumericalValue();
		}
		  
              for ( int S4 = 0 ; S4 <= l4 ; S4++ ) 
                {
		    double S4Coeff = 0.0;
                    if ( (l4 + Q - m4 - S4) >= 0 && (Q + m4 + S4) >= 0 ) 
		      {
			FactorialCoefficient* S4Factorial = new FactorialCoefficient(); 
			S4Factorial->FactorialMultiply((long)(l4+2*Q));
			S4Factorial->FactorialDivide((long)(l4 + Q - m4 - S4));
			S4Factorial->FactorialDivide((long)(Q + m4 + S4));
			S4Coeff = pow(-1.0,(double)S4)*S4Factorial->GetNumericalValue();
		      }
                    Value += NormCoeff * S1Coeff * S2Coeff * S3Coeff * S4Coeff * CalculateBetaFunction((long)(2*Q - m1 - m2 + 1 + l1 + l2 + l3 + l4 - S1 - S2 - S3 - S4),(long)(2*Q + m1 + m2 + 1 + S1 + S2 + S3 + S4));
                }  
            }
        }
    }
  return Value;
}

// evaluate normalisation for Q, l, m
//
// Q = max angular momentum on LLL.
// l = Q + LL 
// m = angular momentum
//
// return value = normalisation

double CalculateNormalization(double Q,double l,double m)
{
  double Value = (2.0*(double)l + 1.0)/(4.0*M_PI);
  FactorialCoefficient* MyFactorial = new FactorialCoefficient(); 
  MyFactorial->FactorialMultiply((long)(l-m));
  MyFactorial->FactorialMultiply((long)(l+m));
  MyFactorial->FactorialDivide((long)(l-Q));
  MyFactorial->FactorialDivide((long)(l+Q));
  return sqrt(Value*MyFactorial->GetNumericalValue());
}

// evaluate the beta function B(x,y) = (x-1)!(y-1)!/(x+y-1)!
//
// x = first arg
// y = second arg
//
// return value = value of beta function with args x,y

double CalculateBetaFunction(long x, long y)
{
  FactorialCoefficient* MyFactorial = new FactorialCoefficient(); 
  MyFactorial->FactorialMultiply(x-1);
  MyFactorial->FactorialMultiply(y-1);
  MyFactorial->FactorialDivide(x+y-1);
  return MyFactorial->GetNumericalValue();
}
