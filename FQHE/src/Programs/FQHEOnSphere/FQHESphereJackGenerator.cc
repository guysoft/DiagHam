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
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleDoubleOption  ('a', "alpha", "alpha coefficient of the Jack polynomial", -2.0);
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the Jack polynomial decomposition into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the Jack polynomial decomposition into a text file");
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
  double Alpha = ((BooleanOption*) Manager["alpha"])->GetBoolean();
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
  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", NbrFluxQuanta) == false) || (NbrFluxQuanta <= 0))
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

  RealVector InputState;
  if (InputState.ReadVector (((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;      
    }

  BosonOnSphereHaldaneBasisShort InitialSpace (NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);	  

  RealVector OutputState = InitialSpace.GenerateJackPolynomial(Alpha);
//  InitialSpace.ConvertToUnnormalizedMonomial(InputState);
  OutputState.WriteVector(((SingleStringOption*) Manager["output-file"])->GetString());

  if (OutputTxtFileName != 0)
    {
      ofstream File;
      File.open(OutputTxtFileName, ios::binary | ios::out);
      for (int i = 0; i < InitialSpace.GetHilbertSpaceDimension(); ++i)
	{
	  File << OutputState[i] << " ";
	  InitialSpace.PrintStateMonomial(File, i) << endl;
	}
      File.close();
    }
  if (OutputFileName != 0)
    {
      OutputState.WriteVector(OutputFileName);
    }
  return 0;
}

