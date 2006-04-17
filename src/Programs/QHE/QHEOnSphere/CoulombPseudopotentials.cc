#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "MathTools/BinomialCoefficients.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "FunctionBasis/QHEFunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"


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
// return value = array that conatins the pseudopotentials
double* EvaluatePseudopotentials(int nbrFlux, int landauLevel);


int main(int argc, char** argv)
{
  OptionManager Manager ("CoulombPseudopotentials" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;


  (*SystemGroup) += new SingleIntegerOption  ('l', "landau-level", "index of the Landau level (0 for the lowest Landau level)", 1, true, 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (i.e. twice the maximum momentum for a single particle)", 8);
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
  
  char* OutputFile;
  if (((SingleStringOption*) Manager["output"])->GetString() == 0l)
    {
      OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%.dat");
    }
  else
    {
      OutputFile = new char [strlen(((SingleStringOption*) Manager["output"])->GetString()) + 1];
      strcpy (OutputFile, ((SingleStringOption*) Manager["output"])->GetString());
    }

  cout <<   OutputFile << endl;
  double* Pseudopotentials = EvaluatePseudopotentials(NbrFlux, LandauLevel);

  ofstream File;
  File.open(OutputFile, ios::binary | ios::out);
  File.precision(14);
  File << "// pseudopotentials on the sphere for coulomb interaction " << endl
       << "// in the Landau level N=" << LandauLevel << " for 2S=" << NbrFlux << " flux quanta" << endl
       << "//" << endl
       << "// PseudoPotentials = V_0 V_1 ..." << endl << endl
       << "PseudoPotentials =";
  for (int i = 0; i <= NbrFlux; ++i)
    File << " " << Pseudopotentials[i];
  File << endl;
  File.close();
  


  
  delete[] OutputFile;
  return 0;
}

// evalute pseudopotentials for coulomb interaction in a given Landau level
//
// nbrFlux = number of flux quanta (i.e. twice the maximum momentum for a single particle)
// landauLevel = index of the Landau level (0 for the lowest Landau level)
// return value = array that conatins the pseudopotentials

double* EvaluatePseudopotentials(int nbrFlux, int landauLevel)
{
  double* Pseudopotentials = new double [nbrFlux + 1];
  ParticleOnSphereGenericLLFunctionBasis Basis(nbrFlux, landauLevel);
  int MaxSum = 2 * nbrFlux;
  int** Indices = new int* [MaxSum + 1];
  for (int i = 0; i <= MaxSum; ++i)
    {
      int NbrSum = (MaxSum >> 1) + 1;
      int* TmpIndices = new int[2 * NbrSum];
      Indices[i] = TmpIndices;
      for (int j = 0; j < NbrSum; ++j)
	{
	  (*TmpIndices) = MaxSum - j;
	  ++TmpIndices;
	  (*TmpIndices) = j;
	  ++TmpIndices;
	}
    }
  return Pseudopotentials;
}

