#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <cstring> 


using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHETorusBosonsFilterOccupation" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  //  ArchitectureManager Architecture;

  Manager += SystemGroup;
  //  Architecture.AddOptionGroup(&Manager);
//  Manager += PrecalculationGroup;
  Manager += MiscGroup;
//  Manager += ToolsGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "vector file that corresponds to the input state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-occ", "maximum occupation per obital", 1);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusBosonsFilterOccupation -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles;
  int NbrFluxQuanta;
  int KxMomentum;
  int KyMomentum;
  double Ratio = 1.0;
  bool Statistics = true;
  if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						  NbrParticles, NbrFluxQuanta, KxMomentum, KyMomentum, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
      return -1;
    }

  cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << NbrFluxQuanta << " Ky=" << KyMomentum << " ";
  if (Statistics == true)
    {
      cout << "FQHETorusBosonsFilterOccupation is only relevant for bosonic states" << endl;
    }
  BosonOnTorusWithMagneticTranslationsShort* SpaceWithTranslations = 0;
  SpaceWithTranslations = new BosonOnTorusWithMagneticTranslationsShort(NbrParticles, NbrFluxQuanta, KxMomentum, KyMomentum);
  ComplexVector InputStateWithTranslations;
  if (InputStateWithTranslations.ReadVector(Manager.GetString("input-state")) == false)
    {
      cout << "error while reading " << Manager.GetString("input-state") << endl;
      return -1;
    }
  if (InputStateWithTranslations.GetVectorDimension() != SpaceWithTranslations->GetHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions " << InputStateWithTranslations.GetVectorDimension() 
	   << " " << SpaceWithTranslations->GetHilbertSpaceDimension() << endl;
      return -1;
    }
  return 0;
}
