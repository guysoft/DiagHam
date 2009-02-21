#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


// get the root parition from a file
// 
// rootFileName = name of the file that contains the root description
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the maximum Lz value
// referenceState = array where the root partition description will be stored
// return value = true if no error occured
bool GetRootPartition (char* rootFileName, int& nbrParticles, int& lzMax, int*& referenceState);


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereUnnormalizeState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "file that describes states to fuse");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state of the output state");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "name of the unnormalized vector that will be generated");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the vector into a text file");
  (*OutputGroup) += new BooleanOption ('\n', "normalize", "normalize the state instead of unnormalizing");  
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereUnnormalizeState -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int LzMax = 0;
  int TotalLz = 0;
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  bool SymmetrizedBasis = ((BooleanOption*) Manager["symmetrized-basis"])->GetBoolean();
  
  bool Statistics = true;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						   NbrParticles, LzMax, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from " << Manager.GetString("input-state") << endl;
      return -1;
    }
  
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName = new char [256];
      if (Statistics == true)
	{
	  if (Manager.GetBoolean("normalize"))	
	    sprintf (OutputFileName, "fermions_normalized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);
	  else
	    sprintf (OutputFileName, "fermions_unnormalized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);
	}
      else
	{
	  if (Manager.GetBoolean("normalize"))	
	    sprintf (OutputFileName, "bosons_normalized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);
	  else
	    sprintf (OutputFileName, "bosons_unnormalized_n_%d_2s_%d_lz_%d.0.vec", NbrParticles, LzMax, TotalLz);
	}
    }

  RealVector OutputState;
  if (OutputState.ReadVector (Manager.GetString("input-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
      return -1;      
    }
  
  ParticleOnSphere* OutputBasis = 0;
  if (Statistics == true)
    {
      if (HaldaneBasisFlag == false)
	OutputBasis = new FermionOnSphere(NbrParticles, TotalLz, LzMax);
      else
	{
	  int* ReferenceState = 0;
	  if (GetRootPartition(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceState) == false)
	    return -1;
	  if (Manager.GetString("load-hilbert") != 0)
	    OutputBasis = new FermionOnSphereHaldaneBasis(Manager.GetString("load-hilbert"));	  
	  else
	    OutputBasis = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, LzMax, ReferenceState);
	}
    }
  else
    {
      if (Manager.GetBoolean("huge-basis") == true)
	{
	  if (Manager.GetString("load-hilbert") == 0)
	    {
	      cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
	      return -1;
	    }
	  OutputBasis = new BosonOnSphereHaldaneHugeBasisShort (Manager.GetString("load-hilbert"), Manager.GetInteger("memory"));
	}
      else
	{
	  if (HaldaneBasisFlag == false)
	    OutputBasis = new BosonOnSphereShort(NbrParticles, TotalLz, LzMax);
	  else
	    {
	      int* ReferenceState = 0;
	      if (GetRootPartition(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceState) == false)
		return -1;
	      if (Manager.GetString("load-hilbert") != 0)
		OutputBasis = new BosonOnSphereHaldaneBasisShort(Manager.GetString("load-hilbert"));	  
	      else
		OutputBasis = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, LzMax, ReferenceState);	  
	    }
	}
    }

  if (Manager.GetBoolean("normalize"))
    OutputBasis->ConvertFromUnnormalizedMonomial(OutputState);
  else
    OutputBasis->ConvertToUnnormalizedMonomial(OutputState);
  
  if (OutputTxtFileName != 0)
    {
      ofstream File;
      File.open(OutputTxtFileName, ios::binary | ios::out);
      File.precision(14);
      for (long i = 0; i < OutputBasis->GetLargeHilbertSpaceDimension(); ++i)
	{
	  File << OutputState[i] << " ";
	  OutputBasis->PrintStateMonomial(File, i) << endl;
	}
      File.close();
    }
  if (OutputFileName != 0)
    {
      if (OutputState.WriteVector(OutputFileName) == false)
	{
	  cout << "error while writing output state " << OutputFileName << endl;
	  return -1;
	}	  
    }

  return 0;
}


// get the root parition from a file
// 
// rootFileName = name of the file that contains the root description
// nbrParticles = reference on the number of particles
// lzMax = reference on twice the maximum Lz value
// referenceState = array where the root partition description will be stored
// return value = true if no error occured

bool GetRootPartition (char* rootFileName, int& nbrParticles, int& lzMax, int*& referenceState)
{
  ConfigurationParser ReferenceStateDefinition;
  if (ReferenceStateDefinition.Parse(rootFileName) == false)
    {
      ReferenceStateDefinition.DumpErrors(cout) << endl;
      return false;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", nbrParticles) == false) || (nbrParticles <= 0))
    {
      cout << "NbrParticles is not defined or as a wrong value" << endl;
      return false;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", lzMax) == false) || (lzMax < 0))
    {
      cout << "LzMax is not defined or as a wrong value" << endl;
      return false;
    }
  int MaxNbrLz;
  if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', referenceState, MaxNbrLz) == false)
    {
      cout << "error while parsing ReferenceState in " << rootFileName << endl;
      return false;     
    }
  if (MaxNbrLz != (lzMax + 1))
    {
      cout << "wrong LzMax value in ReferenceState" << endl;
      return false;     
    }
  return true;
}
