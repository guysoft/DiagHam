#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/QuantumDot/Spectra/PeriodicSpectra.h"


#include <iostream>
#include <stdlib.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ostream;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{  
  cout.precision(14);
  OptionManager Manager ("PeriodicWaveFunctionValue" , "0.01");
  OptionGroup* PositionGroup = new OptionGroup ("Position options");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* FileGroup = new OptionGroup ("file options");

  Manager += PositionGroup;
  Manager += HilbertSpaceGroup;
  Manager += FileGroup;
  Manager += MiscGroup;

  (*PositionGroup) += new SingleDoubleOption ('x', "x-position", "position in X direction in the length of the big box unit", 0.5);
  (*PositionGroup) += new SingleDoubleOption ('y', "y-position", "position in Y direction in the length of the big box unit", 0.5);
  (*PositionGroup) += new SingleDoubleOption ('z', "z-position", "position in Z direction in the length of the big box unit", 0.5);

  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statex", "number of states in x direction", 31);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowx", "lower impulsion in x direction", -15);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statey", "number of states in y direction", 31);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowy", "lower impulsion in y direction", -15);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statez", "number of states in z direction", 31);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -15);

  (*FileGroup) += new SingleStringOption ('f', "input", "the name of the input file", "eigenvector.0");

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsLaplacianDelta -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  double PositionX = ((SingleDoubleOption*) Manager["x-position"])->GetDouble();
  double PositionY = ((SingleDoubleOption*) Manager["y-position"])->GetDouble();
  double PositionZ = ((SingleDoubleOption*) Manager["z-position"])->GetDouble();

  int NbrStateX = ((SingleIntegerOption*) Manager["nbr-statex"])->GetInteger();
  int LowImpulsionX = ((SingleIntegerOption*) Manager["lowx"])->GetInteger();
  int NbrStateY = ((SingleIntegerOption*) Manager["nbr-statey"])->GetInteger();
  int LowImpulsionY = ((SingleIntegerOption*) Manager["lowy"])->GetInteger();
  int NbrStateZ = ((SingleIntegerOption*) Manager["nbr-statez"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();
  
  char* FileName = ((SingleStringOption*) Manager["input"])->GetString();

  Periodic3DOneParticle* Space = new Periodic3DOneParticle(NbrStateX, LowImpulsionX, NbrStateY, LowImpulsionY, NbrStateZ, LowImpulsionZ);

  PeriodicSpectra Spectra (Space, FileName);

  double Real = 0.0, Imaginary = 0.0;
  
  Spectra.WaveFunctionValue (PositionX, 1.0, PositionY, 1.0, PositionZ, 1.0, Real, Imaginary);
  
  cout << "The probability at considered the point: " << Real << " " << Imaginary << endl;

  return 1;
}
