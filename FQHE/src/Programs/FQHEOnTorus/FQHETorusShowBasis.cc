#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"
#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"

#include "HilbertSpace/FermionOnTorusWithSpinNew.h"
#include "HilbertSpace/FermionOnTorusWithSpin.h"
#include "HilbertSpace/BosonOnTorusWithSpin.h"
#include "HilbertSpace/BosonOnTorusWithSpinAndMagneticTranslations.h"
#include "HilbertSpace/FermionOnTorusWithSpinAndMagneticTranslations.h"

#include "HilbertSpace/BosonOnTorusWithSU3Spin.h"
#include "HilbertSpace/BosonOnTorusWithSU3SpinAndMagneticTranslations.h"
#include "HilbertSpace/BosonOnTorusWithSU4Spin.h"
#include "HilbertSpace/BosonOnTorusWithSU4SpinAndMagneticTranslations.h"

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
  (*SystemGroup) += new BooleanOption  ('\n', "2-ll", "consider particles with 2 Landau levels");
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (only useful in su(2) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "su3-spin", "consider particles with SU(3) spin");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-tz", "twice the quantum number of the system associated to the Tz generator (only useful in su(3) mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "total-y", "three time the quantum number of the system associated to the Y generator (only useful in su(3) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "su4-spin", "consider particles with SU(4) spin");
  (*SystemGroup) += new SingleIntegerOption  ('i', "total-isosz", "twice the z component of the total isospin of the system (only usefull in su(4) mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('e', "total-entanglement", "twice the projection of the total spin-isopsin entanglement of the system (only usefull in su(4) mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n1", "number of type 1 particles (only useful in su(3) mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n2", "number of type 2 particles (only useful in su(3) mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n3", "number of type 3 particles (only useful in su(3) mode)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-n4", "number of type 4 particles (only useful in su(4) mode)", 0);
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

  if (Manager.GetBoolean("su3-spin") == true)
    {
      int Ky = Manager.GetInteger("ky-momentum") % NbrFluxQuanta;
      int TotalTz = Manager.GetInteger("total-tz");
      int TotalY = Manager.GetInteger("total-y");
      if ((Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n2") + Manager.GetInteger("nbr-n3")) == NbrParticles)
	{
	  TotalTz = (Manager.GetInteger("nbr-n1") - Manager.GetInteger("nbr-n2"));
	  TotalY = (Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n2") - (2 * Manager.GetInteger("nbr-n3")));
	}
      if (Manager.GetBoolean("no-translation") == true)
	{
	  ParticleOnSphereWithSU3Spin* Space;
	  if (Manager.GetBoolean("boson") == true)
	    {
	      Space = new BosonOnTorusWithSU3Spin(NbrParticles, TotalTz, TotalY, NbrFluxQuanta, Ky);
	    }
	  else
	    {
	      cout << "fermions with SU(3) spin is not implemented" << endl;
	      return 0;
	    }
	  if (Manager.GetString("state") == 0)
	    {
	      for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
		Space->PrintState(cout, i) << endl;;
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
	      if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
		{
		  cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
		  return -1;
		}
	      if (Manager.GetDouble("hide-component") > 0.0)
		{
		  double Error = Manager.GetDouble("hide-component");
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    if (fabs(State[i]) > Error)
		      Space->PrintState(cout, i) << " : "  << State[i] << endl;;
		}
	      else
		for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		  Space->PrintState(cout, i) << " : "  << State[i] << endl;
	    }
	  delete Space;
	  return 0;
	}	  
      else
	{
	  int Kx = Manager.GetInteger("kx-momentum") % MomentumModulo;
	  ParticleOnTorusWithSU3SpinAndMagneticTranslations* Space;
	  if (Manager.GetBoolean("boson") == true)
	    {
	      Space = new BosonOnTorusWithSU3SpinAndMagneticTranslations(NbrParticles, TotalTz, TotalY, NbrFluxQuanta, Kx, Ky);
	    }
	  else
	    {
	      cout << "fermions with SU(3) spin and mangnetic translations is not implemented" << endl;
	      return 0;
	    }
	  if (Manager.GetString("state") == 0)
	    {
	      for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
		Space->PrintState(cout, i) << endl;;
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
	      if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
		{
		  cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
		  return -1;
		}
	      if (Manager.GetDouble("hide-component") > 0.0)
		{
		  double Error = Manager.GetDouble("hide-component");
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    if (Norm(State[i]) > Error)
		      Space->PrintState(cout, i) << " : "  << State[i] << endl;;
		}
	      else
		for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		  Space->PrintState(cout, i) << " : "  << State[i] << endl;
	    }
	  delete Space;
	  return 0;
	}	  
    }

  if (Manager.GetBoolean("su4-spin") == true)
    {
      int TotalSz = Manager.GetInteger("total-sz");
      int TotalIz = Manager.GetInteger("total-isosz");
      int TotalPz = Manager.GetInteger("total-entanglement");
      if ((Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n2") + Manager.GetInteger("nbr-n3") + Manager.GetInteger("nbr-n4")) == NbrParticles)
	{
	  TotalSz = (Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n2")) - (Manager.GetInteger("nbr-n3") + Manager.GetInteger("nbr-n4"));
	  TotalIz = (Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n3")) - (Manager.GetInteger("nbr-n2") + Manager.GetInteger("nbr-n4"));
	  TotalPz = (Manager.GetInteger("nbr-n1") + Manager.GetInteger("nbr-n4")) - (Manager.GetInteger("nbr-n2") + Manager.GetInteger("nbr-n3"));
	}
      int Ky = Manager.GetInteger("ky-momentum") % NbrFluxQuanta;
      if (Manager.GetBoolean("no-translation") == true)
	{
	  ParticleOnSphereWithSU3Spin* Space;
	  if (Manager.GetBoolean("boson") == true)
	    {
	      Space = new BosonOnTorusWithSU3Spin(NbrParticles, TotalSz, TotalIz, TotalPz, NbrFluxQuanta, Ky);
	    }
	  else
	    {
	      cout << "fermions with SU(4) spin is not implemented" << endl;
	      return 0;
	    }
	  if (Manager.GetString("state") == 0)
	    {
	      for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
		Space->PrintState(cout, i) << endl;;
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
	      if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
		{
		  cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
		  return -1;
		}
	      if (Manager.GetDouble("hide-component") > 0.0)
		{
		  double Error = Manager.GetDouble("hide-component");
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    if (fabs(State[i]) > Error)
		      Space->PrintState(cout, i) << " : "  << State[i] << endl;;
		}
	      else
		for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		  Space->PrintState(cout, i) << " : "  << State[i] << endl;
	    }
	  delete Space;
	  return 0;
	}
      else
	{
	  int Kx = Manager.GetInteger("kx-momentum") % MomentumModulo;
	  ParticleOnTorusWithSU4SpinAndMagneticTranslations* Space;
	  if (Manager.GetBoolean("boson") == true)
	    {
	      Space = new BosonOnTorusWithSU4SpinAndMagneticTranslations(NbrParticles, TotalSz, TotalIz, TotalPz, NbrFluxQuanta, Kx, Ky);
	    }
	  else
	    {
	      cout << "fermions with SU(4) spin and mangnetic translations is not implemented" << endl;
	      return 0;
	    }
	  if (Manager.GetString("state") == 0)
	    {
	      for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
		Space->PrintState(cout, i) << endl;;
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
	      if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
		{
		  cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
		  return -1;
		}
	      if (Manager.GetDouble("hide-component") > 0.0)
		{
		  double Error = Manager.GetDouble("hide-component");
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    if (Norm(State[i]) > Error)
		      Space->PrintState(cout, i) << " : "  << State[i] << endl;;
		}
	      else
		for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		  Space->PrintState(cout, i) << " : "  << State[i] << endl;
	    }
	  delete Space;
	  return 0;
	}	  
    }

  if (Manager.GetBoolean("su2-spin") == false)
    {
      if(Manager.GetBoolean("2-ll") == false)
	{
	  if (Manager.GetBoolean("no-translation") == false)
	    {
	      int Kx = Manager.GetInteger("kx-momentum") % MomentumModulo;
	      int Ky = Manager.GetInteger("ky-momentum") % NbrFluxQuanta;
	      ParticleOnTorusWithMagneticTranslations* Space;
	      if (Manager.GetBoolean("boson") == true)
		{
		  Space = new BosonOnTorusWithMagneticTranslationsShort  (NbrParticles, NbrFluxQuanta, Kx, Ky);
		}
	      else
		{
		  Space = new FermionOnTorusWithMagneticTranslations (NbrParticles, NbrFluxQuanta, Kx, Ky);
		}
	      if (Manager.GetString("state") == 0)
		{
		  for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
		    Space->PrintState(cout, i) << endl;;
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
		  if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
		    {
		      cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
		      return -1;
		    }
		  if (Manager.GetDouble("hide-component") > 0.0)
		    {
		      double Error = Manager.GetDouble("hide-component");
		      for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			if (Norm(State[i]) > Error)
			  Space->PrintState(cout, i) << " : "  << State[i] << endl;;
		    }
		  else
		    for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		      Space->PrintState(cout, i) << " : "  << State[i] << endl;;
		}
	      delete Space;
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
      else
	{
	  if (Manager.GetBoolean("no-translation") == false)
	    {
	      cout << "particles on torus with with 2 Landau levels are not implemented" << endl;
	      return -1;
	      // 	      for (int x = 0; x < MomentumModulo; ++x)
	      // 		for (int y = 0; y < MomentumModulo; ++y)
	      // 		  {
	      // 		  }
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
		  ParticleOnSphereWithSpin* Space;
		  if (Manager.GetBoolean("boson") == true)
		    {
		      Space = new BosonOnTorusWithSpin  (NbrParticles, NbrFluxQuanta, y);
		    }
		  else
		    {
		      //Space = new FermionOnTorusWithSpin  (NbrParticles, NbrFluxQuanta,  y);
		      
		      cout <<"fermions on torus with 2 Landau levels are not implemented"<<endl;
		    }
		  cout << " (k_y = " << y << ") : " << endl;
		  if (Manager.GetString("state") == 0)
		    {
		      for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
			Space->PrintState(cout, i) << endl;
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
		      if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
			{
			  cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
			  return -1;
			}
		      if (Manager.GetDouble("hide-component") > 0.0)
			{
			  double Error = Manager.GetDouble("hide-component");
			  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			    if (Norm(State[i]) > Error)
			      Space->PrintState(cout, i) << " : "  << State[i] << endl;;
			}
		      else
			for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			  Space->PrintState(cout, i) << " : "  << State[i] << endl;;
		    }
		  delete Space;
		  
		}
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("no-translation") == false)
	{
	  int Kx = Manager.GetInteger("kx-momentum") % MomentumModulo;
	  int Ky = Manager.GetInteger("ky-momentum") % NbrFluxQuanta;
	  ParticleOnTorusWithSpinAndMagneticTranslations* Space;
	  if (Manager.GetBoolean("boson") == true)
	    {
	      Space = new BosonOnTorusWithSpinAndMagneticTranslations(NbrParticles, Manager.GetInteger("total-sz"), NbrFluxQuanta, Kx, Ky);
	    }
	  else
	    {
	      Space = new FermionOnTorusWithSpinAndMagneticTranslations(NbrParticles, Manager.GetInteger("total-sz"), NbrFluxQuanta, Kx, Ky);
	    }
	  if (Manager.GetString("state") == 0)
	    {
	      for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
		Space->PrintState(cout, i) << endl;;
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
	      if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
		{
		  cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
		  return -1;
		}
	      if (Manager.GetDouble("hide-component") > 0.0)
		{
		  double Error = Manager.GetDouble("hide-component");
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    if (Norm(State[i]) > Error)
		      Space->PrintState(cout, i) << " : "  << State[i] << endl;;
		}
	      else
		for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		  Space->PrintState(cout, i) << " : "  << State[i] << endl;
	    }
	  delete Space;
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
	      ParticleOnSphereWithSpin* Space;
	      if (Manager.GetBoolean("boson") == true)
		{
		  Space = new BosonOnTorusWithSpin  (NbrParticles, NbrFluxQuanta, Manager.GetInteger("total-sz"), y);
		}
	      else
		{
		  Space = new FermionOnTorusWithSpin  (NbrParticles, NbrFluxQuanta, Manager.GetInteger("total-sz"), y);
		}
	      cout << " (k_y = " << y << ") : " << endl;
	      if (Manager.GetString("state") == 0)
		{
		  for (int i = 0; i <  Space->GetHilbertSpaceDimension(); ++i)
		    Space->PrintState(cout, i) << endl;
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
		  if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
		    {
			  cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
			  return -1;
		    }
		  if (Manager.GetDouble("hide-component") > 0.0)
		    {
		      double Error = Manager.GetDouble("hide-component");
		      for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
			if (Norm(State[i]) > Error)
			  Space->PrintState(cout, i) << " : "  << State[i] << endl;;
		    }
		  else
		    for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		      Space->PrintState(cout, i) << " : "  << State[i] << endl;;
		}
	      delete Space;
	    }
	}
    }
}

