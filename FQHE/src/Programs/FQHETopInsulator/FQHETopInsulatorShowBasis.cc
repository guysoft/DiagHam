#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"

#include "Vector/ComplexVector.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <fstream>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHETopInsulatorShowBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "kx", "total momentum along the x direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ky", "total momentum along the y direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-subbands", "number of subbands", 1);
  (*SystemGroup) += new SingleStringOption ('\n', "state", "name of an optional vector state whose component values can be displayed behind each corresponding n-body state");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hide-component", "hide state components (and thus the corresponding n-body state) whose absolute value is lower than a given error (0 if all components have to be shown", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "no-autodetect", "do not autdetect system parameter from state file name");
  //  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  //  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_topinsulator_nbrsubbands_n_nbrparticles_x_nbrsitex_y_nbrsitey.dim");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETopInsulatorShowBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSiteX = Manager.GetInteger("nbr-sitex"); 
  int NbrSiteY = Manager.GetInteger("nbr-sitey"); 
  int TotalKx = Manager.GetInteger("kx"); 
  int TotalKy = Manager.GetInteger("ky");
  if ((Manager.GetString("state") != 0) && (Manager.GetBoolean("no-autodetect") == false))
    {
      double Mass = 0.0;
      bool Statistics = false;
      if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(Manager.GetString("state"),
							      NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy, Mass, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("state") << endl;
	  return -1;
	}	  
    }
 
  AbstractQHEParticle* Space;
  if (Manager.GetInteger("nbr-subbands") == 1)
    {
      Space = new FermionOnSquareLatticeMomentumSpace (NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
    }
  else
    {
      if (Manager.GetInteger("nbr-subbands") == 2)
	{
	  Space = new FermionOnSquareLatticeWithSpinMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalKx, TotalKy);
	}
      else
	{
	  return 0;
	}
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
	    {
	      if (Norm(State[i]) > Error)
		Space->PrintState(cout, i) << " : "  << State[i] << endl;
	      else
		{
		  WeightHiddenComponents += SqrNorm(State[i]);
		  NbrHiddenComponents++;
		}
	      Normalization += SqrNorm(State[i]);
	    }
	  cout << NbrHiddenComponents << " hidden components (square normalization error = " << WeightHiddenComponents << " / " << Normalization << ")" << endl;
	}
      else
	for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	  Space->PrintState(cout, i) << " : "  << State[i] << endl;;
    }
  delete Space;
  return 0;
}
