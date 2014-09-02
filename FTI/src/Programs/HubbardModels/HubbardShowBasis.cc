#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Options/Options.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("HubbardShowBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sites", "number of flux quanta", 20);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "gutzwiller", "use the Gutzwiller projection");
  (*SystemGroup) += new SingleStringOption ('\n', "get-index", "find the index of a given n-body state");
  (*SystemGroup) += new BooleanOption  ('\n', "xperiodic-boundary", "use periodic boundary conditions in the x direction");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "x-momentum", "momentum along the x direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "x-periodicity", "periodicity in the number of site index that implements the periodic boundary condition in the x direction", 4);
  (*SystemGroup) += new BooleanOption  ('\n', "2dperiodic-boundaries", "use periodic boundary conditions in the x and y directions");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "y-momentum", "set the momentum along the y direction", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "y-periodicity", "periodicity in the number of site index that implements the periodic boundary condition in the y direction", 2);
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use the Sz <-> -Sz symmetry");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-parity", "select the  Sz <-> -Sz parity (can be 1 or -1", 1);
  
  (*SystemGroup) += new BooleanOption  ('\n', "add-index", "add index of the Hilbert space vectors");
  
  (*SystemGroup) += new SingleStringOption ('\n', "state", "name of an optional vector state whose component values can be displayed behind each corresponding n-body state");
  (*SystemGroup) += new BooleanOption  ('c', "complex-vector" , "vector state is complex instead of real");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hide-component", "hide state components (and thus the corresponding n-body state) whose absolute value is lower than a given error (0 if all components have to be shown", 0.0);
//   (*SystemGroup) += new SingleStringOption ('\n', "get-index", "find the index of a given n-body state");
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_hubbard_suN_n_nbrparticles_q_nbrfluxquanta_z_totallz.basis");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereShowBasis -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSites = Manager.GetInteger("nbr-sites"); 
  bool AddIndex = Manager.GetBoolean("add-index");
  bool ComplexFlag = Manager.GetBoolean("complex-vector");
  
  int NbrOrbitals;
  
  ParticleOnSphere* Space = 0;
  if (Manager.GetBoolean("boson") == true)
    {
      cout << "bosonic Hubbard model not implemented" << endl;
      return -1;
    }
  else
    {
      if (Manager.GetBoolean("xperiodic-boundary") == false)
	{
	  if (Manager.GetBoolean("2dperiodic-boundaries") == false)
	    {
	      if (Manager.GetBoolean("gutzwiller") == false)
		{
		  Space = new FermionOnLatticeWithSpinRealSpace(NbrParticles, NbrSites);
		}
	      else
		{
		  Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpace(NbrParticles, NbrSites);
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("gutzwiller") == false)
		{
		  Space = new FermionOnLatticeWithSpinRealSpaceAnd2DTranslation(NbrParticles, NbrSites, Manager.GetInteger("x-momentum"), Manager.GetInteger("x-periodicity"),
										Manager.GetInteger("y-momentum"), Manager.GetInteger("y-periodicity"));
		}
	      else
		{
		  Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd2DTranslation(NbrParticles, NbrSites, 
												       Manager.GetInteger("x-momentum"), Manager.GetInteger("x-periodicity"),
												       Manager.GetInteger("y-momentum"), Manager.GetInteger("y-periodicity"));
		}
	    }
	}
      else
	{
	  if (Manager.GetBoolean("szsymmetrized-basis") == false)
	    {
	      if (Manager.GetBoolean("gutzwiller") == false)
		{
		  Space = new FermionOnLatticeWithSpinRealSpaceAnd1DTranslation(NbrParticles, NbrSites, Manager.GetInteger("x-momentum"), Manager.GetInteger("x-periodicity"));
		}
	      else
		{
		  Space = new FermionOnLatticeWithSpinAndGutzwillerProjectionRealSpaceAnd1DTranslation(NbrParticles, NbrSites, Manager.GetInteger("x-momentum"), Manager.GetInteger("x-periodicity"));
		}
	    }
	  else
	    {
	      bool MinusParitySector = true;
	      if (Manager.GetInteger("sz-parity") == 1)
		MinusParitySector = false;
	      if (Manager.GetBoolean("gutzwiller") == false)
		Space = new FermionOnLatticeWithSpinSzSymmetryRealSpaceAnd1DTranslation (NbrParticles, NbrSites, Manager.GetInteger("x-momentum"), Manager.GetInteger("x-periodicity"), MinusParitySector);
	      else
		Space = new FermionOnLatticeWithSpinSzSymmetryAndGutzwillerProjectionRealSpaceAnd1DTranslation (NbrParticles, NbrSites, Manager.GetInteger("x-momentum"), 
														Manager.GetInteger("x-periodicity"),MinusParitySector);
	    }
	}
    }
  
  if (Manager.GetString("get-index") != 0)
    {
      long TmpIndex = Space->FindStateIndex(Manager.GetString("get-index"));
      if (TmpIndex == Space->GetHilbertSpaceDimension())
	{
	  cout << "state " << Manager.GetString("get-index") << " not found" << endl;
	}
      else
	{
	  cout << TmpIndex << " : ";
	  Space->PrintState(cout, TmpIndex) << endl;	   
	}
      return 0;
    }


  ofstream File;
  char* OutputFileName = 0;
  if (Manager.GetBoolean("save-disk") == true)
    {
      if (Manager.GetString("output-file") == 0)
	{
	  OutputFileName = new char[512];
	  if (Manager.GetBoolean("boson") == true)
	  {
	   cout << "bosonic model not implemented" << endl;
	   return -1;
	  }
	  else
	    sprintf (OutputFileName, "fermions_hubbard_n_%d_x_%d.basis", NbrParticles, NbrSites);
	}
      else
	{
	  OutputFileName = new char[strlen(Manager.GetString("output-file")) + 1];
	  strcpy (OutputFileName, Manager.GetString("output-file"));
	}
      File.open(OutputFileName, ios::binary | ios::out);
      if (!File.is_open())
	{
	  cout << "Cannot create file " << OutputFileName << endl;
	  return -1;
	}
      File.precision(14);
    }
  else
    cout.precision(14);
  if (Manager.GetString("state") == 0)
    {
      if (Manager.GetBoolean("save-disk") == true)
	{
	  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	    {
	      
		  if (AddIndex == true) 
		    File << i << " ";
		  Space->PrintState(File, i);
		  File<<endl;
		
	    }
	}
      else
	{
	  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	    {
	      
		  if (AddIndex == true) 
		    cout << i <<" ";
		  Space->PrintState(cout, i);
		  cout<<endl;
		
	    }
	}
    }
  else
   if (ComplexFlag == false)
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
	      if (Manager.GetBoolean("save-disk") == true)
		{
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    {
		      if (fabs(State[i]) > Error)
			{
			    if (AddIndex == true) 
			      File << i << " ";	
			    Space->PrintState(File, i) << " : " << State[i];
			    File<<endl;
			}
		      else		     
			{
			  WeightHiddenComponents += State[i] * State[i];
			  ++NbrHiddenComponents; 
			}
		      Normalization += State[i] * State[i];
		    }
		}
	      else
		{
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    {
		      if (fabs(State[i]) > Error)
			{
			  
			  if (AddIndex == true) 
			    cout << i <<" ";
			  Space->PrintState(cout, i) << " : " << State[i];
			  cout<<endl;
			  
			  
			}
		      else		     
			{
			  WeightHiddenComponents += State[i] * State[i];
			  ++NbrHiddenComponents; 
			}
		      Normalization += State[i] * State[i];
		    }
		  cout << NbrHiddenComponents << " hidden components (square normalization error = " << WeightHiddenComponents << " / " << Normalization << ")" << endl;
		}
	    }
    cout << NbrHiddenComponents << " hidden components (square normalization error = " << WeightHiddenComponents << " / " << Normalization << ")" << endl;
    }

  else //complex vector
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
	      if (Manager.GetBoolean("save-disk") == true)
		{
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    {
		      if (Norm(State[i]) > Error)
			{
			  if (AddIndex == true) 
			    File << i << " ";	
			  Space->PrintState(File, i) << " : " << State.Re(i)<<" "<<State.Im(i);
			  File<<endl;
			}
		      else		     
			{
			  WeightHiddenComponents += State.Re(i) * State.Re(i) + State.Im(i) * State.Im(i);
			  ++NbrHiddenComponents; 
			}
		      Normalization += State.Re(i) * State.Re(i) + State.Im(i) * State.Im(i);
		    }
		}
	      else
		{
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    {
		      if (Norm(State[i]) > Error)
			{
			  if (AddIndex == true) 
			    cout << i <<" ";
			  Space->PrintState(cout, i) << " : " << State.Re(i)<<" "<<State.Im(i);
			  cout<<endl;
			}
		      else		     
			{
			  WeightHiddenComponents += State.Re(i) * State.Re(i) + State.Im(i) * State.Im(i);
			  ++NbrHiddenComponents; 
			}
		      Normalization += State.Re(i) * State.Re(i) + State.Im(i) * State.Im(i);
		    }
		  cout << NbrHiddenComponents << " hidden components (square normalization error = " << WeightHiddenComponents << " / " << Normalization << ")" << endl;
		}
	    }
	  else
	    {    
	      if (Manager.GetBoolean("save-disk") == true)
		{		  
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    {
		      if (AddIndex == true) 
			File << i << " ";	
		      Space->PrintState(File, i) << " : " << " " << State.Re(i) << " " << State.Im(i) << endl;		  
		    }
		}
	      else
		{
		  
		  for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		    {
		      if (AddIndex == true) 
			cout << i << " ";	
		      Space->PrintState(cout, i) << " : " << " " << State.Re(i) << " " << State.Im(i) << endl;
		    }
		}	      	      
	    }
    }

 
  if (OutputFileName != 0)
    {
      File.close();
      delete[] OutputFileName;
    }
}
