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
  (*SystemGroup) += new BooleanOption  ('c', "complex", "run for a complex vector");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type VectorPhaseMultiply -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  

  if (Manager.GetBoolean("complex"))
    {
      ComplexVector State;
      if (State.ReadVector (Manager.GetString("vector")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("vector") << endl;
	  return -1;      
	}
      if (Manager.GetBoolean("normalize-phase"))
	{
	  double MaxNorm=0.0;
	  int MaxIndex=0;
	  for (int i = 0; i < State.GetVectorDimension();++i)
	    {
	      if (Norm(State[i])>MaxNorm)
		{
		  MaxNorm=Norm(State[i]);
		  MaxIndex=i;
		}
	    }
	  cout << "MaxNorm="<<MaxNorm<<", correcting phase with "<<Polar(1.0,-Arg(State[MaxIndex]))<<endl;
	  State *= Polar(1.0,-Arg(State[MaxIndex]));
	  cout << endl<<State;
	}
      if (Manager.GetBoolean("careful-normalize"))
	{
	  cout << "Attention, careful-normalize is not implemented yet for complex numbers" << endl;
	}
      State /= State.Norm();

      if (State.WriteVector (Manager.GetString("vector")) == false)
	{
	  cout << "can't overwrite vector file " << Manager.GetString("vector") << endl;
	  return -1;
	}
    }
  else
    {
      RealVector State;
      if (State.ReadVector (Manager.GetString("vector")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("vector") << endl;
	  return -1;      
	}
      
      if (Manager.GetBoolean("careful-normalize") == false)
	{
	  State /= State.Norm();
	}
      else
	{
	  long double Norm = 0.0;
	  for (long i = 0l; i < State.GetLargeVectorDimension();++i)
	    Norm += ((long double) State[i]) * ((long double) State[i]);
	  Norm = sqrtl(Norm);
	  State /= (double) Norm;      
	}
      if (State.WriteVector (Manager.GetString("vector")) == false)
	{
	  cout << "can't overwrite vector file " << Manager.GetString("vector") << endl;
	  return -1;
	}
    }
  return 0;
}
