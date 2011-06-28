#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"

#include "HilbertSpace/FermionOnTorusWithSpinNew.h"
#include "HilbertSpace/FermionOnTorusWithSpin.h"

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
  cout.precision(14);

  OptionManager Manager ("FQHETorusShowBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta", 8);
  (*SystemGroup) += new SingleIntegerOption ('x', "kx-momentum", "the total momentum along the x axis (only available when using magnetic translations, should be lower than GDC(p, l))", 0);
  (*SystemGroup) += new SingleIntegerOption ('y', "ky-momentum", "the total momentum along the y axis", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "consider particles with SU(2) spin");
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (only useful in su(2) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "do not consider magnetic translation (only trivial translations along one axis)");
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
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrFluxQuanta = Manager.GetInteger("nbr-flux"); 
  int MomentumModulo = FindGCD(NbrParticles, NbrFluxQuanta);

  if (Manager.GetBoolean("su2-spin") == false)
    {
      if (Manager.GetBoolean("no-translation") == false)
	{
	  int Kx = Manager.GetInteger("kx-momentum") % MomentumModulo;
	  int Ky = Manager.GetInteger("ky-momentum") % NbrFluxQuanta;
	  if (Manager.GetBoolean("boson") == true)
	    {
	      BosonOnTorusWithMagneticTranslations Space (NbrParticles, NbrFluxQuanta, Kx, Ky);
	      for (int i = 0; i <  Space.GetHilbertSpaceDimension(); ++i)
		Space.PrintState(cout, i) << endl;
	      cout << endl;
	    }
	  else
	    {
	      FermionOnTorusWithMagneticTranslations Space (NbrParticles, NbrFluxQuanta, Kx, Ky);
	      if (Manager.GetString("state") == 0)
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
		  if (State.ReadVector(Manager.GetString("state")) == false)
		    {
		      cout << "error while reading " << Manager.GetString("state") << endl;
		      return -1;
		    }
		  if (Space.GetHilbertSpaceDimension() != State.GetVectorDimension())
		    {
		      cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space.GetHilbertSpaceDimension() << ")" << endl;
		      return -1;
		    }
		  if (Manager.GetDouble("hide-component") > 0.0)
		    {
		      double Error = Manager.GetDouble("hide-component");
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
      else
	{
	  int MinKy = 0;
	  int MaxKy = NbrFluxQuanta;
	  if (Manager.GetInteger("ky-momentum") >= 0)
	    {
	      MinKy = Manager.GetInteger("ky-momentum") % NbrFluxQuanta;
	      MaxKy = MinKy + 1;
	    }
	  for (int y = MinKy; y < MaxKy; ++y)
	    {
	      if (Manager.GetBoolean("boson") == true)
		{
		  BosonOnTorusShort Space (NbrParticles, NbrFluxQuanta, y);
		  cout << " (k_y = " << y << ") : " << endl;
		  if (Manager.GetString("state") == 0)
		    {
		      for (int i = 0; i <  Space.GetHilbertSpaceDimension(); ++i)
			Space.PrintState(cout, i) << endl;
		      cout << endl;
		    }
		  else
		    {
		      int NbrHiddenComponents = 0;
		      double WeightHiddenComponents = 0.0;
		      double Normalization = 0.0;
		      RealVector State;
		      if (State.ReadVector(Manager.GetString("state")) == false)
			{
			  cout << "error while reading " << Manager.GetString("state") << endl;
			  return -1;
			}
		      if (Space.GetHilbertSpaceDimension() != State.GetVectorDimension())
			{
			  cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space.GetHilbertSpaceDimension() << ")" << endl;
			  return -1;
			}
		      if (Manager.GetDouble("hide-component") > 0.0)
			{
			  double Error = Manager.GetDouble("hide-component");
			  for (int i = 0; i < Space.GetHilbertSpaceDimension(); ++i)
			    if (Norm(State[i]) > Error)
			      Space.PrintState(cout, i) << " : "  << State[i] << endl;;
			}
		      else
			for (int i = 0; i < Space.GetHilbertSpaceDimension(); ++i)
			  Space.PrintState(cout, i) << " : "  << State[i] << endl;;
		    }
		}
	      else
		{
		  FermionOnTorus Space (NbrParticles, NbrFluxQuanta, y);
		  cout << " (k_y = " << y << ") : " << endl;
		  if (Manager.GetString("state") == 0)
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
		      if (State.ReadVector(Manager.GetString("state")) == false)
			{
			  cout << "error while reading " << Manager.GetString("state") << endl;
			  return -1;
			}
		      if (Space.GetHilbertSpaceDimension() != State.GetVectorDimension())
			{
			  cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space.GetHilbertSpaceDimension() << ")" << endl;
			  return -1;
			}
		      if (Manager.GetDouble("hide-component") > 0.0)
			{
			  double Error = Manager.GetDouble("hide-component");
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
    }
  else
    {
      if (Manager.GetBoolean("boson") == true)
	{
	  cout << "bosons with spin on torus are not implemented" << endl;
	  return -1;
	}
      else
	{
	  if (Manager.GetBoolean("no-translation") == false)
	    {
	      for (int x = 0; x < MomentumModulo; ++x)
		for (int y = 0; y < MomentumModulo; ++y)
		  {
		  }
	    }
	  else
	    {
	      int MinKy = 0;
	      int MaxKy = NbrFluxQuanta;
	      if (Manager.GetInteger("ky-momentum") >= 0)
		{
		  MinKy = Manager.GetInteger("ky-momentum") % NbrFluxQuanta;
		  MaxKy = MinKy + 1;
		}
	      for (int y = MinKy; y < MaxKy; ++y)
		{
		  //		  FermionOnTorusWithSpinNew Space (NbrParticles, Manager.GetInteger("total-sz"), NbrFluxQuanta, y);
		  FermionOnTorusWithSpin Space (NbrParticles, NbrFluxQuanta, Manager.GetInteger("total-sz"), y);
		  cout << " (k_y = " << y << ") : " << endl;
		  if (Manager.GetString("state") == 0)
		    {
		      for (int i = 0; i <  Space.GetHilbertSpaceDimension(); ++i)
			Space.PrintState(cout, i) << endl;
		      cout << endl;
		    }
		  else
		    {
		      int NbrHiddenComponents = 0;
		      double WeightHiddenComponents = 0.0;
		      double Normalization = 0.0;
		      RealVector State;
		      if (State.ReadVector(Manager.GetString("state")) == false)
			{
			  cout << "error while reading " << Manager.GetString("state") << endl;
			  return -1;
			}
		      if (Space.GetHilbertSpaceDimension() != State.GetVectorDimension())
			{
			  cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space.GetHilbertSpaceDimension() << ")" << endl;
			  return -1;
			}
		      if (Manager.GetDouble("hide-component") > 0.0)
			{
			  double Error = Manager.GetDouble("hide-component");
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
    }
}

