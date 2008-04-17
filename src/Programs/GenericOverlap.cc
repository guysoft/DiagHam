#include "Vector/RealVector.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

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
  OptionManager Manager ("GenericOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('1', "state1", "name of the file containing the 1st vector obtained using exact diagonalization");
  (*SystemGroup) += new SingleStringOption  ('2', "state2", "name of the file containing the 2nd vector obtained using exact diagonalization");
  (*SystemGroup) += new BooleanOption  ('\n', "discard-sign", "compute sum_i |v1_i * v2_i| instead of sum_i v1_i * v2_i");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericOverlap -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (((SingleStringOption*) Manager["state1"])->GetString() == 0)
    {
      cout << "GenericOverlap requires an exact state" << endl;
      return -1;
    }
  RealVector State1;
  if (State1.ReadVector (((SingleStringOption*) Manager["state1"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["state1"])->GetString() << endl;
      return -1;      
    }
  if (((SingleStringOption*) Manager["state2"])->GetString() == 0)
    {
      cout << "GenericOverlap requires a 2nd state to compare to." << endl;
      return -1;
    }
  RealVector State2;
  if (State2.ReadVector (((SingleStringOption*) Manager["state2"])->GetString()) == false)
    {
      cout << "can't open vector file " << ((SingleStringOption*) Manager["state2"])->GetString() << endl;
      return -1;      
    }

  if (State1.GetVectorDimension() != State2.GetVectorDimension() )
    {
      cout << "Dimension of Hilbert spaces in input files does not coincide" << endl;
      return -2;
    }

  double sp=0.0;
  if (((BooleanOption*) Manager["discard-sign"])->GetBoolean() == true)
    for (int i=0; i<State1.GetVectorDimension(); ++i)
      sp+=fabs(State1[i]*State2[i]);
  else
    for (int i=0; i<State1.GetVectorDimension(); ++i)
      sp+= State1[i]*State2[i];


  cout << "The overlap is: |<1|2>|^2 = " << sp*sp << endl;
  
}
