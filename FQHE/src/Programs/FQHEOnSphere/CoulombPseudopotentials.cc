#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "MathTools/BinomialCoefficients.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "Vector/RealVector.h"
 
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using std::ofstream;
using std::ios;
using std::cout;
using std::endl;


// evalute pseudopotentials for coulomb interaction in a given Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// layerSeparation = layer separation d in bilayer, or layer thickness d modeled by interaction 1/sqrt(r^2+d^2)
// return value = array that conatins the pseudopotentials
double* EvaluatePseudopotentials(int nbrFlux, int landauLevel, double layerSeparation);

// evalute one body potentials for two impurities located at the poles in a given Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// northPolePotential = potential of the impurity located at the north pole
// southPolePotential = potential of the impurity located at the south pole
// return value = array that conatins the pseudopotentials
double* EvaluateOneBodyPotentials(int nbrFlux, int landauLevel, double northPolePotential, double southPolePotential);


int main(int argc, char** argv)
{
  OptionManager Manager ("CoulombPseudopotentials" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;


  (*SystemGroup) += new SingleIntegerOption  ('l', "landau-level", "index of the Landau level (0 for the lowest Landau level)", 0, true, 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (i.e. twice the maximum momentum for a single particle)", 8);
  (*SystemGroup) += new BooleanOption ('\n', "add-impurities", "add two impurities (one at each pole)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "north-potential", "potential assosciated to the impurity at the north pole", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "south-potential", "potential assosciated to the impurity at the south pole", 0.0);
  (*SystemGroup) += new  BooleanOption ('\n', "relativistic-fermions", "assume relativistic fermions");

  (*SystemGroup) += new  SingleDoubleOption ('d', "layer-separation", "assume finite layer separation / thickness",0.0);

  (*SystemGroup) += new SingleStringOption  ('o', "output", "output file name (default is pseudopotential_coulomb_l_x_2s_y.dat)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type CoulombPseudopotentials -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int LandauLevel = ((SingleIntegerOption*) Manager["landau-level"])->GetInteger();
  int NbrFlux = ((SingleIntegerOption*) Manager["nbr-flux"])->GetInteger();
  double layerSeparation = Manager.GetDouble("layer-separation");
  int MaxMomentum = NbrFlux + (LandauLevel << 1);
  
  char* OutputFile;
  if (((SingleStringOption*) Manager["output"])->GetString() == 0l)
    {
      if (((BooleanOption*) Manager["relativistic-fermions"])->GetBoolean() == true)
	OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_relativistic_l_%landau-level%_2s_%nbr-flux%.dat");
      else
	{
	  if (layerSeparation==0.0)
	    OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%.dat");
	  else
	    OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%_d_%layer-separation%.dat");
	}
    }
  else
    {
      OutputFile = new char [strlen(((SingleStringOption*) Manager["output"])->GetString()) + 1];
      strcpy (OutputFile, ((SingleStringOption*) Manager["output"])->GetString());
    }

  double* Pseudopotentials = EvaluatePseudopotentials(NbrFlux, LandauLevel, layerSeparation);
  if (((BooleanOption*) Manager["relativistic-fermions"])->GetBoolean() == true)
    {
      double* PseudopotentialsNMinus1 = EvaluatePseudopotentials(NbrFlux, LandauLevel - 1, layerSeparation);
      for (int i = 0; i <= MaxMomentum; ++i)
	Pseudopotentials[i] = 0.5 * (Pseudopotentials[i] + PseudopotentialsNMinus1[i]);
      delete[] PseudopotentialsNMinus1;
    }
  double* OneBodyPotentials = 0;
  if (((BooleanOption*) Manager["add-impurities"])->GetBoolean() == true)
    OneBodyPotentials = EvaluateOneBodyPotentials(NbrFlux, LandauLevel, 
						  ((SingleDoubleOption*) Manager["north-potential"])->GetDouble(),
						  ((SingleDoubleOption*) Manager["south-potential"])->GetDouble());
  ofstream File;
  File.open(OutputFile, ios::binary | ios::out);
  File.precision(14);
  File << "# pseudopotentials on the sphere for coulomb interaction ";
  if (((BooleanOption*) Manager["relativistic-fermions"])->GetBoolean() == true)
    File << " for relativistic fermions";
  File << endl
       << "# in the Landau level N=" << LandauLevel << " for 2S=" << NbrFlux << " flux quanta" << endl;
  if (layerSeparation != 0.0)
    File << "# with finite layer separation / thickness d=" << layerSeparation << endl;
  if (OneBodyPotentials != 0)
    {
      File << "# with two impurities (V_north = " << ((SingleDoubleOption*) Manager["north-potential"])->GetDouble() 
	   << " and V_south = " << ((SingleDoubleOption*) Manager["south-potential"])->GetDouble() << ")" << endl;
    }
  File << "#" << endl
       << "# Pseudopotentials = V_0 V_1 ..." << endl << endl
       << "Pseudopotentials =";
  for (int i = 0; i <= MaxMomentum; ++i)
    File << " " << Pseudopotentials[i];
  File << endl;
  if (OneBodyPotentials != 0)
    {
      File << endl << "Onebodypotentials =";
      for (int i = 0; i <= MaxMomentum; ++i)
	File << " " << OneBodyPotentials[i];
     File << endl;
    }
  File.close();
  


  
  delete[] OutputFile;
  return 0;
}

// evalute pseudopotentials for coulomb interaction in a given Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// layerSeparation = layer separation d in bilayer, or layer thickness d modeled by interaction 1/sqrt(r^2+d^2)
// return value = array that conatins the pseudopotentials

double* EvaluatePseudopotentials(int nbrFlux, int landauLevel, double layerSeparation)
{
  cout.precision(14);
  int MaxMomentum = nbrFlux + (landauLevel << 1);
  double* Pseudopotentials = new double [MaxMomentum + 1];
  ClebschGordanCoefficients MainCoefficients(MaxMomentum, MaxMomentum);
  ClebschGordanCoefficients* Coefficients = new ClebschGordanCoefficients[MaxMomentum + 1];
  // new formfactors for finite thickness/separation:
  double *FormFactors = new double[MaxMomentum + 1];
  double dd = layerSeparation*layerSeparation;
  double base = ( sqrt(2.0*nbrFlux + dd) - layerSeparation ) / ( sqrt(2.0*nbrFlux + dd) + layerSeparation );
  FormFactors[0]=sqrt(base);
  for (int l = 1; l <= MaxMomentum; ++l)
    FormFactors[l] = FormFactors[l-1]*base;
  for (int l = 0; l <= MaxMomentum; ++l)
    Coefficients[l] = ClebschGordanCoefficients(MaxMomentum, l << 1);  
  for (int l = 0; l <= MaxMomentum; ++l)
    {
      double TmpPseudopotentials = 0.0;
      for (int m1 = -MaxMomentum; m1 <= MaxMomentum; m1 +=2)
	for (int m2 = -MaxMomentum; m2 <= MaxMomentum; m2 +=2)
	  {
	    double TmpCoef = 0.0;
	    double TmpCoef2;
	    int Min = abs(m1 - m2) >> 1;
	    double Sign = 1.0;
	    for (int j = Min; j <= MaxMomentum; ++j)
	      {
		TmpCoef2 = Coefficients[j].GetCoefficient(m1, m2 - m1, MaxMomentum) * Coefficients[j].GetCoefficient(nbrFlux, 0, MaxMomentum);
		TmpCoef += FormFactors[j] * Sign * TmpCoef2 * TmpCoef2;
		Sign *= -1.0;
	      }
	    TmpPseudopotentials += (MainCoefficients.GetCoefficient(m1, -m1, l << 1) * 
				    MainCoefficients.GetCoefficient(m2, -m2, l << 1) * TmpCoef);
	  }
      Pseudopotentials[MaxMomentum - l] = TmpPseudopotentials / sqrt (0.5 * nbrFlux);
      cout << "V[" << (MaxMomentum - l) << "] = " << Pseudopotentials[MaxMomentum - l] << endl;
    }
  delete[] Coefficients;
  delete[] FormFactors;
  return Pseudopotentials;
}

// evalute one body potentials for two impurities located at the poles in a given Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// northPolePotential = potential of the impurity located at the north pole
// southPolePotential = potential of the impurity located at the south pole
// return value = array that conatins the pseudopotentials

double* EvaluateOneBodyPotentials(int nbrFlux, int landauLevel, double northPolePotential, double southPolePotential)
{
  int MaxMomentum = nbrFlux + (landauLevel << 1);
  double* OneBodyPotentials = new double [MaxMomentum + 1];
  for (int i = 0; i <= MaxMomentum; ++i)
   OneBodyPotentials[i] = 0.0;
  ParticleOnSphereGenericLLFunctionBasis Basis (nbrFlux, landauLevel);
  RealVector Value(2, true);
  Complex TmpValue;
  Value[0] = M_PI;
  Basis.GetFunctionValue(Value, TmpValue, landauLevel);
  OneBodyPotentials[landauLevel] = southPolePotential * SqrNorm(TmpValue);
  Value[0] = 0.0;
  Basis.GetFunctionValue(Value, TmpValue, MaxMomentum - landauLevel);
  OneBodyPotentials[MaxMomentum - landauLevel] = northPolePotential * SqrNorm(TmpValue);  
  return OneBodyPotentials;
}

