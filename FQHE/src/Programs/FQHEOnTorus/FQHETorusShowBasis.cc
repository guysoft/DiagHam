#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Vector/ComplexVector.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETorusShowBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('q', "nbr-flux", "number of flux quanta", 8);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new SingleStringOption ('\n', "state", "name of an optional vector state whose component values can be displayed behind each corresponding n-body state");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hide-component", "hide state components (and thus the corresponding n-body state) whose absolute value is lower than a given error (0 if all components have to be shown", 0.0);
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_sphere_n_nbrparticles_q_nbrfluxquanta_z_totallz.basis");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusShowBasis -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int NbrFluxQuanta = ((SingleIntegerOption*) Manager["nbr-flux"])->GetInteger(); 
  int MomentumModulo = FindGCD(NbrParticles, NbrFluxQuanta);

  for (int x = 0; x < MomentumModulo; ++x)
    for (int y = 0; y < MomentumModulo; ++y)
      {
	if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	  {
	    BosonOnTorusWithMagneticTranslations Space (NbrParticles, NbrFluxQuanta, x, y);
	    cout << " (k_x = " << x << ", k_y = " << y << ") : " << endl;
	    for (int i = 0; i <  Space.GetHilbertSpaceDimension(); ++i)
	      Space.PrintState(cout, i) << endl;
	    cout << endl;
	  }
	else
	  {
	    FermionOnTorusWithMagneticTranslations Space (NbrParticles, NbrFluxQuanta, x, y);
	    cout << " (k_x = " << x << ", k_y = " << y << ") : " << endl;
	    if (((SingleStringOption*) Manager["state"])->GetString() == 0)
	      {
		for (int i = 0; i <  Space.GetHilbertSpaceDimension(); ++i)
		  Space.PrintState(cout, i) << endl;;
		cout << endl;
	      }
	    else
	      {
		int NbrHiddenComponents = 0;
		double WeightHiddenComponents = 0.0;
		double Normalization = 0.0;
		ComplexVector State;
		if (State.ReadVector(((SingleStringOption*) Manager["state"])->GetString()) == false)
		  {
		    cout << "error while reading " << ((SingleStringOption*) Manager["state"])->GetString() << endl;
		    return -1;
		  }
		if (Space.GetHilbertSpaceDimension() != State.GetVectorDimension())
		  {
		    cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space.GetHilbertSpaceDimension() << ")" << endl;
		    return -1;
		  }
		if (((SingleDoubleOption*) Manager["hide-component"])->GetDouble() > 0.0)
		  {
		    double Error = ((SingleDoubleOption*) Manager["hide-component"])->GetDouble();
		    for (int i = 0; i < Space.GetHilbertSpaceDimension(); ++i)
		      if (Norm(State[i]) > Error)
			Space.PrintState(cout, i) << " : "  << State[i] << endl;;
		  }
		else
		  for (int i = 0; i < Space.GetHilbertSpaceDimension(); ++i)
		    Space.PrintState(cout, i) << " : "  << State[i] << endl;;
	      }
	  }
      }
}

