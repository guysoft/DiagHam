#include "Vector/RealVector.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

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
  OptionManager Manager ("FQHESphereJackGenerator" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleDoubleOption  ('a', "alpha", "alpha coefficient of the Jack polynomial", -2.0);
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new SingleStringOption  ('\n', "initial-state", "use an optional state where some of the components have already been computed, improving computation time");
  (*SystemGroup) += new BooleanOption  ('\n', "fermion", "compute the slater decomposition of the Jack polynomial times Vandermonde");
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the Jack polynomial decomposition into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the Jack polynomial decomposition into a text file");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereJackGenerator -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int NbrParticles = 0; 
  int NbrFluxQuanta = 0; 
  bool SymmetrizedBasis = ((BooleanOption*) Manager["symmetrized-basis"])->GetBoolean();
  double Alpha = ((SingleDoubleOption*) Manager["alpha"])->GetDouble();
  int TotalLz = 0;
  char* OutputFileName = ((SingleStringOption*) Manager["bin-output"])->GetString();
  char* OutputTxtFileName = ((SingleStringOption*) Manager["txt-output"])->GetString();

  if ((OutputTxtFileName == 0) && (OutputFileName == 0))
    {
      cout << "error, an output file (binary or text) has to be provided" << endl;
      return 0;
    }

  int* ReferenceState = 0;
  if (((SingleStringOption*) Manager["reference-file"])->GetString() == 0)
    {
      cout << "error, a reference file is needed" << endl;
      return 0;
    }
  ConfigurationParser ReferenceStateDefinition;
  if (ReferenceStateDefinition.Parse(((SingleStringOption*) Manager["reference-file"])->GetString()) == false)
    {
      ReferenceStateDefinition.DumpErrors(cout) << endl;
      return 0;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
    {
      cout << "NbrParticles is not defined or as a wrong value" << endl;
      return 0;
    }
  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", NbrFluxQuanta) == false) || (NbrFluxQuanta < 0))
    {
      cout << "LzMax is not defined or as a wrong value" << endl;
      return 0;
    }
  int MaxNbrLz;
  if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
    {
      cout << "error while parsing ReferenceState in " << ((SingleStringOption*) Manager["reference-file"])->GetString() << endl;
      return 0;     
    }
  if (MaxNbrLz != (NbrFluxQuanta + 1))
    {
      cout << "wrong LzMax value in ReferenceState" << endl;
      return 0;     
    }

  if (Manager.GetBoolean("fermion") == false)
    {
      BosonOnSphereHaldaneBasisShort* InitialSpace;
      if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
	InitialSpace = new BosonOnSphereHaldaneBasisShort(((SingleStringOption*) Manager["load-hilbert"])->GetString());
      else
	{
	  InitialSpace = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);	  
	  if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
	    {
	      InitialSpace->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
	      return 0;
	    }
	}
      RealVector OutputState;
      if (Manager.GetString("initial-state") == 0)
	OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
      else
	if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
	  {
	    cout << "can't open " << Manager.GetString("initial-state") << endl;
	    return -1;
	  }
      if (SymmetrizedBasis == false)    
	InitialSpace->GenerateJackPolynomial(OutputState, Alpha);
      else
	InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, Alpha);
      
      if (OutputTxtFileName != 0)
	{
	  ofstream File;
	  File.open(OutputTxtFileName, ios::binary | ios::out);
	  File.precision(14);
	  for (long i = 0; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      File << OutputState[i] << " ";
	      InitialSpace->PrintStateMonomial(File, i) << endl;
	    }
	  File.close();
	}
      if (OutputFileName != 0)
	{
	  OutputState.WriteVector(OutputFileName);
	}
    }
  else
    {
      FermionOnSphereHaldaneBasis* InitialSpace;
      if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
	InitialSpace = new FermionOnSphereHaldaneBasis(((SingleStringOption*) Manager["load-hilbert"])->GetString());
      else
	{
	  InitialSpace = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);	  
	  if (((SingleStringOption*) Manager["save-hilbert"])->GetString() != 0)
	    {
	      InitialSpace->WriteHilbertSpace(((SingleStringOption*) Manager["save-hilbert"])->GetString());
	      return 0;
	    }
	}
      RealVector OutputState;
      if (Manager.GetString("initial-state") == 0)
	OutputState = RealVector(InitialSpace->GetLargeHilbertSpaceDimension(), true);
      else
	if (OutputState.ReadVector(Manager.GetString("initial-state")) == false)
	  {
	    cout << "can't open " << Manager.GetString("initial-state") << endl;
	    return -1;
	  }
      if (SymmetrizedBasis == false)    
	InitialSpace->GenerateJackPolynomial(OutputState, Alpha);
      else
	InitialSpace->GenerateSymmetrizedJackPolynomial(OutputState, Alpha);
      
      if (OutputTxtFileName != 0)
	{
	  ofstream File;
	  File.open(OutputTxtFileName, ios::binary | ios::out);
	  File.precision(14);
	  for (long i = 0; i < InitialSpace->GetLargeHilbertSpaceDimension(); ++i)
	    {
	      File << OutputState[i] << " ";
	      InitialSpace->PrintStateMonomial(File, i) << endl;
	    }
	  File.close();
	}
      if (OutputFileName != 0)
	{
	  OutputState.WriteVector(OutputFileName);
	}
    }    
  return 0;
}

