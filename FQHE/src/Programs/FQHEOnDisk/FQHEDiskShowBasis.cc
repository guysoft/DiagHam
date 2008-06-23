#include "HilbertSpace/BosonOnDisk.h"
#include "HilbertSpace/FermionOnDisk.h"
#include "HilbertSpace/FermionOnDiskUnlimited.h"

#include "Vector/RealVector.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

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
  OptionManager Manager ("FQHEDiskShowBasis" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "total lz value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new SingleStringOption ('\n', "state", "name of an optional vector state whose component values can be displayed behind each corresponding n-body state");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "hide-component", "hide state components (and thus the corresponding n-body state) whose absolute value is lower than a given error (0 if all components have to be shown", 0.0);
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_disk_suN_n_nbrparticles_q_nbrfluxquanta_z_totallz.basis");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEDiskShowBasis -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int TotalLz = ((SingleIntegerOption*) Manager["lz-value"])->GetInteger();
    
  ParticleOnDisk* Space = 0;
  if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
    {
      Space = new BosonOnDisk(NbrParticles, TotalLz);
    }
  else
    {
#ifdef __64_BITS__
      if ((TotalLz - (((NbrParticles - 1) * (NbrParticles - 2)) / 2)) < 63)      
#else
      if ((TotalLz - (((NbrParticles - 1) * (NbrParticles - 2)) / 2)) < 31)
#endif
	Space = new FermionOnDisk(NbrParticles, TotalLz);
      else
	Space = new FermionOnDiskUnlimited(NbrParticles, TotalLz);
    }


  ofstream File;
  char* OutputFileName = 0;
  if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
    {
      if (((SingleStringOption*) Manager["output-file"])->GetString() == 0)
	{
	  OutputFileName = new char[512];
	  if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
	    sprintf (OutputFileName, "bosons_disk_n_%d_lz_%d.basis", NbrParticles, TotalLz);
	  else
	    sprintf (OutputFileName, "fermions_disk_n_%d_lz_%d.basis", NbrParticles, TotalLz);
	}
      else
	{
	  OutputFileName = new char[strlen(((SingleStringOption*) Manager["output-file"])->GetString()) + 1];
	  strcpy (OutputFileName, ((SingleStringOption*) Manager["output-file"])->GetString());
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
  if (((SingleStringOption*) Manager["state"])->GetString() == 0)
    if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
      for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	Space->PrintState(File, i) << endl;
    else
      for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	Space->PrintState(cout, i) << endl;
      
   else
     {
       int NbrHiddenComponents = 0;
       double WeightHiddenComponents = 0.0;
       double Normalization = 0.0;
       RealVector State;
       if (State.ReadVector(((SingleStringOption*) Manager["state"])->GetString()) == false)
	 {
	   cout << "error while reading " << ((SingleStringOption*) Manager["state"])->GetString() << endl;
	   return -1;
	 }
       if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
	 {
	   cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
	   return -1;
	 }
       if (((SingleDoubleOption*) Manager["hide-component"])->GetDouble() > 0.0)
	 {
	   double Error = ((SingleDoubleOption*) Manager["hide-component"])->GetDouble();
	   if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
	     {
	       for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
		 {
		   if (fabs(State[i]) > Error)
		     Space->PrintState(File, i) << " : " << State[i] << endl;	   
		   else		     
		     {
		       WeightHiddenComponents += State[i] * State[i];
		       ++NbrHiddenComponents; 
		     }
		   Normalization += State[i] * State[i];
		 }
	     }
	   else
	     for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	       {
		 if (fabs(State[i]) > Error)
		   Space->PrintState(cout, i) << " : " << State[i] << endl;	   
		 else		     
		   {
		     WeightHiddenComponents += State[i] * State[i];
		     ++NbrHiddenComponents; 
		   }
		 Normalization += State[i] * State[i];
	       }
	   cout << NbrHiddenComponents << " hidden components (square normalization error = " << WeightHiddenComponents << " / " << Normalization << ")" << endl;
	 }
       else
	 {
	   if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
	     for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	       Space->PrintState(File, i) << " : " << State[i] << endl;	   
	   else
	     for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	       Space->PrintState(cout, i) << " : " << State[i] << endl;	   
	 }
     }

  if (OutputFileName != 0)
    {
      File.close();
      delete[] OutputFileName;
    }
}

