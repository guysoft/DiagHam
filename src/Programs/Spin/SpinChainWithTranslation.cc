#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/BlockDiagonalMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Vector/RealVector.h"
#include "HilbertSpace/SpinHilbertSpace/Spin1ChainWithTranslations.h"
#include "GeneralTools/List.h"
#include "GeneralTools/ListIterator.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithm.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "Architecture/SMPArchitecture.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);
 
  int NbrSpinPerChain = 3;
  if (argc >= 2)
    NbrSpinPerChain = atoi (argv[1]);

  for (int i = 0; i < NbrSpinPerChain; ++i)
    {
      Spin1ChainWithTranslations Chain(NbrSpinPerChain, i, 0, 10000000, 10000000);
      cout << "--------------------------------------------" << endl;
      cout << "momentum = " << i << endl;
      cout << "dimension = " << Chain.GetHilbertSpaceDimension() << endl;
//      for (int i = 0; i < Chain.GetHilbertSpaceDimension(); ++i)
//	Chain.PrintState(cout, i) << endl;
    }
  return 0;
}

