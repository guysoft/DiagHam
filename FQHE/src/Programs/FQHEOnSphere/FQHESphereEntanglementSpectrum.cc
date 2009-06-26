#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereEntanglementSpectrum" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "density-matrix", "file containing the reduced density matrix");  
  (*SystemGroup) += new SingleIntegerOption ('n', "nbr-particles", "compute the Jack polynomial from the min-index-th component (require an initial state)", 0l);
  (*SystemGroup) += new SingleIntegerOption ('l', "nbr-orbitals", "compute the Jack polynomial from the max-index-th component (require an initial state, 0 if it has computed up to the end)", 0l);
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output name for the entanglement spectrum (default name replace density-matrix full.ent extension with la_x_na_y.entspec)");
  (*SystemGroup) += new SingleDoubleOption  ('e', "eigenvalue-error", "lowest acceptable reduced density matrix eignvalue", 1e-14);  
  (*SystemGroup) += new BooleanOption ('\n', "show-minmaxlza", "show minimum an maximum Lz value that can be reached");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereEntanglementSpectrum -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrOrbitalsInPartition = Manager.GetInteger("nbr-orbitals");
  int NbrParticlesInPartition = Manager.GetInteger("nbr-particles");
  double Error = Manager.GetDouble("eigenvalue-error");
  int NbrParticles = 0;
  int NbrFluxQuanta = 0;
  int TotalLz = 0;
  bool Statistics = true;
  
  if (Manager.GetString("density-matrix") == 0)
    {
      cout << "a reduced density matrix has to be provided, see man page for option syntax or type FQHESphereEntanglementSpectrum -h" << endl;
      return -1;
    }

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("density-matrix"), NbrParticles, NbrFluxQuanta, TotalLz, Statistics) == false)
    {
      cout << "can't retrieve system informations from the reduced density matrix file name" << endl;
      return -1;
    } 

  MultiColumnASCIIFile DensityMatrix;
  if (DensityMatrix.Parse(Manager.GetString("density-matrix")) == false)
    {
      DensityMatrix.DumpErrors(cout);
      return -1;
    }
  if (DensityMatrix.GetNbrColumns() < 4)
    {
      cout << "wrong number of columns in " << Manager.GetString("density-matrix") << endl;
      return -1;
    }

  int* LaValues = DensityMatrix.GetAsIntegerArray(0);
  int* NaValues = DensityMatrix.GetAsIntegerArray(1);
  int* LzValues = DensityMatrix.GetAsIntegerArray(2);
  double* Coefficients = DensityMatrix.GetAsDoubleArray(3);
  double MinLza = 1e300;
  double MaxLza = -MinLza; 
  int Index = 0l;
  int MaxIndex = DensityMatrix.GetNbrLines();
  while ((Index < MaxIndex) && (LaValues[Index] != NbrOrbitalsInPartition))
    ++Index;


  if (Index < MaxIndex)
    {
      char* OutputFileName = Manager.GetString("output");
      if (OutputFileName == 0)
	{
	  char* TmpExtension = new char[256];
	  sprintf(TmpExtension, "la_%d_na_%d.entspec", NbrOrbitalsInPartition, NbrParticlesInPartition);
	  if (strcasestr(Manager.GetString("density-matrix"), "bz2") == 0)
	    {
	      OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.ent", TmpExtension);
	    }
	  else
	    {
	      OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.ent.bz2", TmpExtension);
	    }
	}
      ofstream File;
      File.open(OutputFileName, ios::out);
      File.precision(14);
      File << "# la na lz shifted_lz lambda -log(lambda)" << endl;
      if (NbrParticlesInPartition == 0)
	{
	  while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition))
	    {
	      double Tmp = Coefficients[Index];
	      if (Tmp > Error)
		{
		  double TmpLza = (-0.5 * (LzValues[Index] + ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NaValues[Index])));
		  File << NbrOrbitalsInPartition << " " << NaValues[Index] << " " << LzValues[Index] << " " <<  TmpLza<< " " << Tmp << " " << (-log(Tmp)) << endl;
		  if (TmpLza < MinLza)
		    MinLza = TmpLza;
		  if (TmpLza > MaxLza)
		    MaxLza = TmpLza;
		}
	      ++Index;
	    }
	}
      else
	{
	  while ((Index < MaxIndex) && (NaValues[Index] != NbrParticlesInPartition))
	    ++Index;
	  if (Index < MaxIndex)
	    {
	      double Shift = ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NbrParticlesInPartition);
	      while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition) && (NaValues[Index] == NbrParticlesInPartition))
		{
		  double Tmp = Coefficients[Index];
		  if (Tmp > Error)
		    {
		      double TmpLza = (-0.5 * (LzValues[Index] + Shift));
		      File << NbrOrbitalsInPartition << " " << NbrParticlesInPartition << " " << LzValues[Index] << " " << TmpLza << " " << Tmp << " " << (-log(Tmp)) << endl;
		      if (TmpLza < MinLza)
			MinLza = TmpLza;
		      if (TmpLza > MaxLza)
			MaxLza = TmpLza;
		    }
		  ++Index;
		}	      
	    }
	  else
	    {
	      cout << "error, no entanglement spectrum can be computed from current data (invalid number of particles)" << endl;	      
	      return -1;
	    }
	}
      File.close();
    }
  else
    {
      cout << "error, no entanglement spectrum can be computed from current data (invalid number of orbitals)" << endl;
      return -1;
    }

  if (Manager.GetBoolean("show-minmaxlza"))
    {
      cout << "min Lza = " << MinLza << endl;
      cout << "max Lza = " << MaxLza << endl;
    }

  return 0;
}

