#include "Options/Options.h"

#include "Tools/FQHESpectrum/PseudoPotentials.h"
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

  (*SystemGroup) += new  BooleanOption ('n', "nbody", "add n-body potentials");
  (*SystemGroup) += new  MultipleDoubleOption ('p', "nbody-potentials", "values of n-body potentials to be added (separated by ','",',');

  (*SystemGroup) += new SingleStringOption  ('o', "output", "output file name (default is pseudopotential_coulomb_l_x_2s_y.dat)");
  (*SystemGroup) += new BooleanOption ('\n', "std-output", "use standard output instead of an output file");
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

  double* Pseudopotentials = EvaluatePseudopotentials(NbrFlux, LandauLevel, layerSeparation, true);
  if (((BooleanOption*) Manager["relativistic-fermions"])->GetBoolean() == true)
    {
      double* PseudopotentialsNMinus1 = EvaluatePseudopotentials(NbrFlux, LandauLevel - 1, layerSeparation, true);
      for (int i = 0; i <= MaxMomentum; ++i)
	Pseudopotentials[i] = 0.5 * (Pseudopotentials[i] + PseudopotentialsNMinus1[i]);
      delete[] PseudopotentialsNMinus1;
    }
  double* OneBodyPotentials = 0;
  if (((BooleanOption*) Manager["add-impurities"])->GetBoolean() == true)
    OneBodyPotentials = EvaluateOneBodyPotentials(NbrFlux, LandauLevel, 
						  ((SingleDoubleOption*) Manager["north-potential"])->GetDouble(),
						  ((SingleDoubleOption*) Manager["south-potential"])->GetDouble());

  if (((BooleanOption*) Manager["std-output"])->GetBoolean() == false)
    {
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
      if (Manager.GetBoolean("nbody"))
	{
	  int Length;
	  double *potentials = Manager.GetDoubles("nbody-potentials",Length);
	  if (potentials != 0)
	    {
	      File << endl << "NbrNBody = "<<Length-1<<endl;
	      File << "Weights =";
	      for (int i=0; i<Length; ++i)
		File <<" "<<potentials[i];
	      File <<endl;
	    }
	  else
	    {
	      File << endl << "NbrNBody = 2"<<endl;
	      File << "Weights = 0 0 0";
	    }
	}
      
      File.close();
    }
  else
    {
      cout.precision(14);
      cout << "# pseudopotentials on the sphere for coulomb interaction ";
      if (((BooleanOption*) Manager["relativistic-fermions"])->GetBoolean() == true)
	cout << " for relativistic fermions";
      cout << endl
	   << "# in the Landau level N=" << LandauLevel << " for 2S=" << NbrFlux << " flux quanta" << endl;
      if (layerSeparation != 0.0)
	cout << "# with finite layer separation / thickness d=" << layerSeparation << endl;
      if (OneBodyPotentials != 0)
	{
	  cout << "# with two impurities (V_north = " << ((SingleDoubleOption*) Manager["north-potential"])->GetDouble() 
	       << " and V_south = " << ((SingleDoubleOption*) Manager["south-potential"])->GetDouble() << ")" << endl;
	}
      cout << "#" << endl
	   << "# Pseudopotentials = V_0 V_1 ..." << endl << endl
	   << "Pseudopotentials =";
      for (int i = 0; i <= MaxMomentum; ++i)
	cout << " " << Pseudopotentials[i];
      cout << endl;
      if (OneBodyPotentials != 0)
	{
	  cout << endl << "Onebodypotentials =";
	  for (int i = 0; i <= MaxMomentum; ++i)
	    cout << " " << OneBodyPotentials[i];
	  cout << endl;
	}
    }  


  
  delete[] OutputFile;
  return 0;
}


