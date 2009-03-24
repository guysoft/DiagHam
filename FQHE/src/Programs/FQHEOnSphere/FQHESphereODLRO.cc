#include "Vector/RealVector.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

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
  OptionManager Manager ("FQHESphereODLRO" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "vector file that describes the smaller system");
  (*SystemGroup) += new SingleStringOption  ('o', "output-state", "vector file that describes the smaller system");
  (*SystemGroup) += new BooleanOption  ('\n', "input-haldane", "use the squeezed basis or the input state");
  (*SystemGroup) += new BooleanOption  ('\n', "output-haldane", "use the squeezed basis or the output state");  
  (*SystemGroup) += new SingleStringOption  ('\n', "input-reference", "use a file as the definition of the reference state of the input state");
  (*SystemGroup) += new SingleStringOption  ('\n', "output-reference", "use a file as the definition of the reference state of the output state");
  (*SystemGroup) += new BooleanOption  ('\n', "inputhuge-basis", "use huge Hilbert space support for the input state");
  (*SystemGroup) += new BooleanOption  ('\n', "outputhuge-basis", "use huge Hilbert space support for the output state");
  (*SystemGroup) += new SingleStringOption  ('\n', "pattern", "pattern that has to be shared between the two n-body states");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "shift-pattern", "shift the pattern away from the pole from a given number of orbitals", 0);

  (*PrecalculationGroup) += new SingleStringOption  ('\n', "inputload-hilbert", "load Hilbert space description from the input file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "outputload-hilbert", "load Hilbert space description from the output file",0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereODLRO -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int OutputNbrParticles = 0;
  int OutputLzMax = 0;
  int OutputTotalLz = 0;
  
  int InputNbrParticles = 0;
  int InputLzMax = 0;
  int InputTotalLz = 0;
  bool Statistics = true;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						   InputNbrParticles, InputLzMax, InputTotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from input state name " << Manager.GetString("input-state") << endl;
      return -1;
    }
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("output-state"),
						   OutputNbrParticles, OutputLzMax, OutputTotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from input state name " << Manager.GetString("input-state") << endl;
      return -1;
    }

  int PatternNbrParticles = 0;
  int PatternLzMax = 0;
  int* Pattern = 0;
  if (FQHEGetRootPartition(Manager.GetString("pattern"), PatternNbrParticles, PatternLzMax, Pattern) == false)
    return -1;

  RealVector InputState;
  if (InputState.ReadVector (Manager.GetString("input-state")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("input-state") << endl;
      return -1;      
    }

  ParticleOnSphere* InputBasis = 0;
  if (Statistics == false)
    {
      if (Manager.GetBoolean("inputhuge-basis") == true)
	{
	  if (Manager.GetString("inputload-hilbert") == 0)
	    {
	      cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
	      return -1;
	    }
	  InputBasis = new BosonOnSphereHaldaneHugeBasisShort (Manager.GetString("inputload-hilbert"), Manager.GetInteger("memory"));
	}
      else
	{
	  if (Manager.GetBoolean("input-haldane") == true)
	    {
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(Manager.GetString("input-reference"), InputNbrParticles, InputLzMax, ReferenceState) == false)
		return -1;
	      if (Manager.GetString("inputload-hilbert") != 0)
		InputBasis = new BosonOnSphereHaldaneBasisShort(Manager.GetString("inputload-hilbert"));	  
	      else
		InputBasis = new BosonOnSphereHaldaneBasisShort(InputNbrParticles, InputTotalLz, InputLzMax, ReferenceState);	  
	    }
	  else
	    {
	      InputBasis = new BosonOnSphereShort(InputNbrParticles, InputTotalLz, InputLzMax);
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("input-haldane") == true)
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("input-reference"), InputNbrParticles, InputLzMax, ReferenceState) == false)
	    return -1;
	  if (Manager.GetString("inputload-hilbert") != 0)
	    InputBasis = new FermionOnSphereHaldaneBasis(Manager.GetString("inputload-hilbert"));	  
	  else
	    InputBasis = new FermionOnSphereHaldaneBasis(InputNbrParticles, InputTotalLz, InputLzMax, ReferenceState);
	}
      else
	{
	  InputBasis = new FermionOnSphere(InputNbrParticles, InputTotalLz, InputLzMax);
	}
    }

  ParticleOnSphere* OutputBasis = 0;
  if (Statistics == false)
    {
      if (Manager.GetBoolean("outputhuge-basis") == true)
	{
	  if (Manager.GetString("outputload-hilbert") == 0)
	    {
	      cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
	      return -1;
	    }
	  OutputBasis = new BosonOnSphereHaldaneHugeBasisShort (Manager.GetString("outputload-hilbert"), Manager.GetInteger("memory"));
	}
      else
	{
	  if (Manager.GetBoolean("output-haldane") == true)
	    {
	      int* ReferenceState = 0;
	      if (FQHEGetRootPartition(Manager.GetString("output-reference"), OutputNbrParticles, OutputLzMax, ReferenceState) == false)
		return -1;
	      if (Manager.GetString("outputload-hilbert") != 0)
		OutputBasis = new BosonOnSphereHaldaneBasisShort(Manager.GetString("outputload-hilbert"));	  
	      else
		OutputBasis = new BosonOnSphereHaldaneBasisShort(OutputNbrParticles, OutputTotalLz, OutputLzMax, ReferenceState);	  
	    }
	  else
	    {
		OutputBasis = new BosonOnSphereShort(OutputNbrParticles, OutputTotalLz, OutputLzMax);	  
	    }
	}
    }
  else
    {
      if (Manager.GetBoolean("output-haldane") == true)
	{
	  int* ReferenceState = 0;
	  if (FQHEGetRootPartition(Manager.GetString("output-reference"), OutputNbrParticles, OutputLzMax, ReferenceState) == false)
	    return -1;
	  if (Manager.GetString("outputload-hilbert") != 0)
	    OutputBasis = new FermionOnSphereHaldaneBasis(Manager.GetString("outputload-hilbert"));	  
	  else
	    OutputBasis = new FermionOnSphereHaldaneBasis(OutputNbrParticles, OutputTotalLz, OutputLzMax, ReferenceState);
	}
      else
	{
	  OutputBasis = new FermionOnSphere(OutputNbrParticles, OutputTotalLz, OutputLzMax);
	}
    }

  RealVector TruncatedState = InputBasis->TruncateStateWithPatternConstraint(InputState, OutputBasis, Pattern, PatternLzMax + 1, Manager.GetInteger("shift-pattern"));
  
//   for (long i = 0; i < OutputBasis->GetLargeHilbertSpaceDimension(); ++i)
//     {
//       cout << TruncatedState[i] << " ";
//       OutputBasis->PrintStateMonomial(cout, i) << endl;
//     }
  if (TruncatedState.Norm() > 1e-6)
    {
      OutputBasis->ConvertFromUnnormalizedMonomial(TruncatedState);
      TruncatedState /= TruncatedState.Norm();
    }

  RealVector OutputState;
  if (OutputState.ReadVector(Manager.GetString("output-state")) == false)
    {
      cout << "can't open " << Manager.GetString("output-state") << endl;
      return -1;
    }

  cout.precision(14); 
  cout << "ODLRO=" << (OutputState * TruncatedState) << endl;



  return 0;
}


