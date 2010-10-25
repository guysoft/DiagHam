#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "Options/Options.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("VectorPhaseMultiply" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-vector", "name of the file containing the input vector");
  (*SystemGroup) += new SingleStringOption  ('o', "output-vector", "name of the file where the output vector has to be stored");
  (*SystemGroup) += new SingleDoubleOption  ('p', "phase", "phase (in radian unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "pi-unit", "the phase is expressed in pi unit instead of radian unit"); 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type VectorPhaseMultiply -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("input-vector") == 0)
    {
      cout << "no input file" << endl;
    }
  if (Manager.GetString("output-vector") == 0)
    {
      cout << "no output file" << endl;
    }

  Complex Phase;
  if (Manager.GetBoolean("pi-unit"))  
    {
      Phase.Re = cos (M_PI * Manager.GetDouble("phase"));
      Phase.Im = sin (M_PI * Manager.GetDouble("phase"));
    }
  else
    {
      Phase.Re = cos (Manager.GetDouble("phase"));
      Phase.Im = sin (Manager.GetDouble("phase"));
    }

  RealVector InputVector1;
  if (InputVector1.ReadVectorTest(Manager.GetString("input-vector")))  
    {
      if (InputVector1.ReadVector(Manager.GetString("input-vector")) == false)
	{
	  return -1;
	}
      ComplexVector OutputVector (InputVector1, true);
      OutputVector *= Phase;
      OutputVector.WriteVector(Manager.GetString("output-vector"));
    }
  else
    {
      ComplexVector OutputVector;
      if (OutputVector.ReadVector(Manager.GetString("input-vector")) == false)
	{
	  return -1;
	}
      OutputVector *= Phase;
      OutputVector.WriteVector(Manager.GetString("output-vector"));
    }
  return 0;
}
