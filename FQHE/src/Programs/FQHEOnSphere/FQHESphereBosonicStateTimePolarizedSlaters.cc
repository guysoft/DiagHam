#include "Vector/RealVector.h"

#include "HilbertSpace/ParticleOnSphere.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/BosonOnSphereTwoLandauLevels.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinHaldaneBasis.h"

#include "Options/Options.h"
#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"


#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <stdio.h>
#include <string.h>

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using std::ifstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  
  OptionManager Manager ("FQHESphere2LLBosonicStateTimePolarizedSlaters" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  
	ArchitectureManager Architecture;
	
  Manager += SystemGroup;
  Manager += MiscGroup;
	Architecture.AddOptionGroup(&Manager);
  Manager += OutputGroup;
	
  (*SystemGroup) += new SingleStringOption ('\0', "state", "name of the vector file in the 2LL");
  (*SystemGroup) += new BooleanOption  ('\n', "2-ll", "consider particles within two Landau levels");    
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the Jack polynomial decomposition into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the Jack polynomial decomposition into a text file");
  // (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  // (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  // (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "lzsymmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "szsymmetrized-basis", "use Sz <-> -Sz symmetrized version of the basis (only valid if total-sz=0, override auto-detection from file name)");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-szparity", "select the  Sz <-> -Sz symmetric sector with negative parity");
  (*SystemGroup) += new BooleanOption  ('\n', "minus-lzparity", "select the  Lz <-> -Lz symmetric sector with negative parity");  
  (*OutputGroup) += new BooleanOption ('\n', "normalize", "the output vector will be normalize on the factory",false);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
	
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsTimesFermions -h" << endl;
      return -1;
    }
	
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
 
  bool LL2 = Manager.GetBoolean("2-ll"); 
 
  RealVector InitialState;
  if (InitialState.ReadVector(Manager.GetString("state")) == false)
    {
	cout << "error while reading " << Manager.GetString("state") << endl;
	return -1;
    }
 
  int NbrParticles = 0;
  int LzMax = 0;
  int TotalLz = 0;
  int TotalSz = 0;
  bool FermionFlag = false;
  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("state"), NbrParticles, LzMax, TotalLz, FermionFlag) == false)
   {
     return -1;
   }
 
   ParticleOnSphere * InitialSpace;
   if ( LL2 == true)
     InitialSpace = new BosonOnSphereTwoLandauLevels (NbrParticles,TotalLz,LzMax+ 2,LzMax);
   else
     InitialSpace = new BosonOnSphereShort (NbrParticles,TotalLz,LzMax);
	
  if (InitialSpace->GetHilbertSpaceDimension() != InitialState.GetVectorDimension())
    {
      cout << "dimension mismatch between the state (" << InitialState.GetVectorDimension() << ") and the Hilbert space (" << InitialSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  int HalfNbrParticles = NbrParticles >> 1;
  FermionOnSphere * SlaterSpace = new FermionOnSphere (HalfNbrParticles,0,HalfNbrParticles-1);
	
  FermionOnSphereWithSpin * FinalSpace;
  if ( Manager.GetBoolean("haldane") == false)
    {
       FinalSpace = new FermionOnSphereWithSpin (NbrParticles,TotalLz, LzMax + HalfNbrParticles - 1,0);
    }
  else if (Manager.GetString("reference-file") != 0 )
    {
      int** ReferenceStates = 0;
      int NbrReferenceStates;		    
      bool TexturelessFlag;
      if (FQHEGetRootPartitionSU2(Manager.GetString("reference-file"), NbrParticles, LzMax, ReferenceStates, NbrReferenceStates, TexturelessFlag) == false)
	{
	  cout << "error while parsing " << Manager.GetString("reference-file") << endl;	      
	  return 0;
	}      		  		         
      if (TexturelessFlag == false ) 
	{	
	  if ( Manager.GetBoolean("lzsymmetrized-basis") && Manager.GetBoolean("szsymmetrized-basis") )
	    {		      
	      bool LzMinusParity = ((BooleanOption*) Manager["minus-lzparity"])->GetBoolean();      
	      bool SzMinusParity = ((BooleanOption*) Manager["minus-szparity"])->GetBoolean();	  	  
	      FinalSpace = new FermionOnSphereWithSpinHaldaneLzSzSymmetry(NbrParticles, LzMax, SzMinusParity, LzMinusParity, ReferenceStates, NbrReferenceStates);
	    }
	  else
	    {
	      FinalSpace = new FermionOnSphereWithSpinHaldaneBasis(NbrParticles, TotalLz, LzMax, TotalSz, ReferenceStates, NbrReferenceStates);
	    }	    
	}
      else
	{			  			    
	  int **TexturelessReferenceState = new int*[NbrReferenceStates];
	  for ( int j = 0 ; j < NbrReferenceStates ; j++ ) 
	    {
	      TexturelessReferenceState[j] = new int[LzMax+1];
	      for ( int i = 0 ; i < (LzMax + 1) ; i++ )
		{
		  if ( ReferenceStates[j][i] == 3 ) 
		    {
			TexturelessReferenceState[j][i] = 2;
		    }
		  else if ( (ReferenceStates[j][i] == 1) || (ReferenceStates[j][i] == 2) ) 
		    {
			TexturelessReferenceState[j][i] = 1;
		    }
		  else
		    {
			TexturelessReferenceState[j][i] = 0;
		    }
		}
	    }	
	  if ( Manager.GetBoolean("lzsymmetrized-basis") && Manager.GetBoolean("szsymmetrized-basis") )
	    {		      
	      bool LzMinusParity = ((BooleanOption*) Manager["minus-lzparity"])->GetBoolean();      
	      bool SzMinusParity = ((BooleanOption*) Manager["minus-szparity"])->GetBoolean();
	      FinalSpace = new FermionOnSphereWithSpinHaldaneLzSzSymmetry(NbrParticles, LzMax, SzMinusParity, LzMinusParity, TexturelessReferenceState, NbrReferenceStates, true);
	    }
	  else
	    {
	      FinalSpace = new FermionOnSphereWithSpinHaldaneBasis(NbrParticles, TotalLz, LzMax, TotalSz, TexturelessReferenceState, NbrReferenceStates, true);
	    }
	}
    }    
  else
    {
      cout << "error, no reference file." << endl;
      return 0;
    }
	
  RealVector OutputVector(FinalSpace->GetHilbertSpaceDimension(),true);
	
  FQHESphereBosonicStateTimesPolarizedSlaterProjectionOperation Operation(InitialSpace, SlaterSpace, FinalSpace, &InitialState, &OutputVector,LL2);
	
  Operation.ApplyOperation(Architecture.GetArchitecture());

  if(Manager.GetBoolean("normalize"))
    FinalSpace->ConvertFromUnnormalizedMonomial(OutputVector,0l,true);
	
  OutputVector.WriteVector(Manager.GetString("bin-output"));	
}

