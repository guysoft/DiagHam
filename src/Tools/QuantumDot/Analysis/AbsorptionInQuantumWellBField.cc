#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Tools/QuantumDot/Spectra/Spectra.h"
#include "Tools/QuantumDot/Spectra/DOSSpectra.h"

#include <iostream>
#include <fstream>
#ifdef __SSTREAM_STYLE__
#include <sstream>
#else
#include <strstream>
#endif
#include <string>
#include <unistd.h>


using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);  
  OptionManager Manager ("AbsorptionInQuantumWellBField" , "0.01");
  OptionGroup* AbsorptionInQuantumWellBFieldGroup = new OptionGroup ("AbsorptionInQuantumWellBField");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  
  Manager += AbsorptionInQuantumWellBFieldGroup;
  Manager += MiscGroup;

  (*AbsorptionInQuantumWellBFieldGroup) += new SingleStringOption('\n', "input", "name of the input file (directory mode only)", "eigenvalues");
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleStringOption('\n', "prefix", "prefix of the directories containing spectrum (directory mode only)", "");
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleIntegerOption('\n', "begin", "number of the first directory (directory mode only)", 0);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleIntegerOption('\n', "end", "number of the last directory (directory mode only)", 0);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleIntegerOption('n', "nbr-state", "number of states", 10);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleDoubleOption('\n', "min", "lower limit of the spectrum (in eV unit)", 0.0);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleDoubleOption('\n', "max", "upper limit of the spectrum (in eV unit)", 0.0);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleDoubleOption('g', "gamma", "full width at half maximum of each Lorentzian peak (in eV unit)", 0.01);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleDoubleOption('\n', "step", "length of each discretized step (in eV unit) in the spectrum", 2e-4);
  (*AbsorptionInQuantumWellBFieldGroup) += new SingleStringOption('\n', "output", "name of the output file", "AbsorptionInQuantumWellBField.txt");

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");
  (*MiscGroup) += new BooleanOption ('v', "verbose", "verbose mode", false);

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type AsorptionInQuantumWellBFiled -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrState = ((SingleIntegerOption*) Manager["nbr-state"])->GetInteger();
  double Min = ((SingleDoubleOption*) Manager["min"])->GetDouble();
  double Max = ((SingleDoubleOption*) Manager["max"])->GetDouble();
  double Gamma = ((SingleDoubleOption*) Manager["gamma"])->GetDouble();
  double Step = ((SingleDoubleOption*) Manager["step"])->GetDouble();
  char* OutputFile = ((SingleStringOption*) Manager["output"])->GetString();

  bool VerboseFlag = ((BooleanOption*) Manager["verbose"])->GetBoolean();


  char* InputFile = ((SingleStringOption*) Manager["input"])->GetString();
  char* Prefix = ((SingleStringOption*) Manager["prefix"])->GetString();
  int Begin = ((SingleIntegerOption*) Manager["begin"])->GetInteger();
  int End = ((SingleIntegerOption*) Manager["end"])->GetInteger();
  int Number = End - Begin + 1;
  char* Prefixbis = new char [100]; char* InputFilebis = new char [100];
  AddString (Prefixbis, Prefix, 0, "");  strcpy(InputFilebis, "/"); strcat(InputFilebis, InputFile);
  char** Files = new char* [Number]; int* State = new int[Number];
  for (int i = Begin; i <= End; ++i)
    {
      State[i - Begin] = NbrState;
      Files[i - Begin] = new char[200];
      if (i < 10)
	AddString(Files[i - Begin], Prefixbis, i, InputFilebis);
      else
	AddString(Files[i - Begin], Prefix, i, InputFilebis);
      if (VerboseFlag)
	cout << Files[i - Begin] << endl;
    }
  DOSSpectra AbsorptionInQuantumWellBField(Number, Files, State, Gamma, Min, Max, Step);
  AbsorptionInQuantumWellBField.WriteSpectra(OutputFile);

  return 1;
}
