#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"

#include "MathTools/ClebschGordanCoefficients.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

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
  OptionManager Manager ("GetBosonsDimension" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta", 20);
  (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the total lz value", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "use fermionic statistic instead of bosonic statistic");
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "consider particles with SU(2) spin");
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system (only usefull in su(2) mode)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "su4-spin", "consider particles with SU(4) spin");
  (*SystemGroup) += new SingleIntegerOption  ('i', "total-isosz", "twice the z component of the total isospin of the system (only usefull in su(4) mode)", 0);
  (*OutputGroup) += new BooleanOption  ('\n', "save-disk", "save output on disk");
  (*OutputGroup) += new SingleStringOption ('\n', "output-file", "use this file name instead of statistics_sphere_n_nbrparticles_q_nbrfluxquanta_z_totallz.basis");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereShowBasis -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int NbrFluxQuanta = ((SingleIntegerOption*) Manager["nbr-flux"])->GetInteger(); 
  int TotalLz = ((SingleIntegerOption*) Manager["lz-value"])->GetInteger();
  if (((NbrParticles * NbrFluxQuanta) & 1) != (TotalLz & 1)) 
    {
      cout << "incompatible values for the number of particles, the number of flux quanta and twice the total lz value (nbr-particles * nbr-flux and lz-value should have the same parity)" << endl;
      return -1;
    }

  ParticleOnSphere* Space;
  if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
    {
      Space = new BosonOnSphere(NbrParticles, TotalLz, NbrFluxQuanta);
    }
  else
    {
      if ((((BooleanOption*) Manager["su4-spin"])->GetBoolean() == false) && (((BooleanOption*) Manager["su2-spin"])->GetBoolean() == false))
	{
#ifdef __64_BITS__
      if (NbrFluxQuanta <= 63)
        {
          Space = new FermionOnSphere(NbrParticles, TotalLz, NbrFluxQuanta);
        }
      else
        {
          Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, NbrFluxQuanta);
        }
#else
      if (NbrFluxQuanta <= 31)
        {
          Space = new FermionOnSphere(NbrParticles, TotalLz, NbrFluxQuanta);
        }
      else
        {
          Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, NbrFluxQuanta);
        }
#endif
	}
//       else
// 	if (((BooleanOption*) Manager["su4-spin"])->GetBoolean() == false)
// 	  {
// 	    int SzTotal = ((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
// 	    Space = new FermionOnSphereWithSU4Spin(NbrParticles, TotalLz, NbrFluxQuanta, SzTotal);
// 	  }
// 	else
// 	  {
// 	    int SzTotal = ((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
// 	    int IsoSzTotal = ((SingleIntegerOption*) Manager["total-isosz"])->GetInteger();
// 	    Space = new FermionOnSphereWithSU4Spin(NbrParticles, TotalLz, NbrFluxQuanta, SzTotal, IsoSzTotal);
// 	  }
    }
  

  ClebschGordanCoefficients Clebsch (NbrFluxQuanta, NbrFluxQuanta);
  for (int m1 = -NbrFluxQuanta; m1 <= NbrFluxQuanta; m1 += 2)
    for (int m2 =  -NbrFluxQuanta; m2 < m1; m2 += 2)
      for (int i = (abs(m1 + m2) >> 1); i <= NbrFluxQuanta; ++i)
	cout << m1 << " " << m2 << " " << i << " " << Clebsch.GetCoefficient(m1, m2, i << 1) << endl;

   if (((BooleanOption*) Manager["save-disk"])->GetBoolean() == true)
     {
//       char* OutputFileName = 0;
//       if (((SingleStringOption*) Manager["output-file"])->GetString() == 0)
// 	{
// 	  OutputFileName = new char[256];
// 	  if (((BooleanOption*) Manager["boson"])->GetBoolean() == true)
// 	    sprintf (OutputFileName, "bosons_sphere_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
// 	  else
// 	    sprintf (OutputFileName, "fermions_sphere_n_%d_2s_%d.dim", NbrParticles, NbrFluxQuanta);
// 	}
//       else
// 	{
// 	  OutputFileName = new char[strlen(((SingleStringOption*) Manager["output-file"])->GetString()) + 1];
// 	  strcpy (OutputFileName, ((SingleStringOption*) Manager["output-file"])->GetString());
// 	}
//       delete[] OutputFileName;
     }
   else
     {
       for (int i = 0; i < Space->GetHilbertSpaceDimension(); ++i)
	 {
	   Space->PrintState(cout, i) << endl;
	 }
     }

}

