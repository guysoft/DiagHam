#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/FermionOnSphereMPSWrapper.h"
#include "HilbertSpace/FermionOnCylinderMPSWrapper.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "Tools/FQHEMPS/FQHEMPSMatrixManager.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Hamiltonian/TensorProductSparseMatrixHamiltonian.h"
#include "Hamiltonian/TensorProductSparseMatrixSelectedBlockHamiltonian.h"

#include "LanczosAlgorithm/BasicArnoldiAlgorithm.h"
#include "LanczosAlgorithm/BasicArnoldiAlgorithmWithDiskStorage.h"

#include "Matrix/SparseRealMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "MainTask/FQHEMPSEMatrixMainTask.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Options/Options.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14); 
  
  OptionManager Manager ("FQHETorusMPSOverlap" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  FQHEMPSMatrixManager MPSMatrixManager (false, true);

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  OptionGroup* ArnoldiGroup  = new OptionGroup ("Arnoldi options");
  Architecture.AddOptionGroup(&Manager);
  Manager += ArnoldiGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-fluxquanta", "set the total number of flux quanta and deduce the number of particles", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n', "topologicalsector", "set the topological sector of state", 0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");  



  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusMPSOverlap -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFluxQuanta = 0;
  if (Manager.GetInteger("nbr-fluxquanta") <= 0)
    {
      cout << "invalid number of flux quanta" << endl;
      return -1;
    }
   NbrFluxQuanta = Manager.GetInteger("nbr-fluxquanta");

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 

  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  int NbrBMatrices = MPSMatrix->GetNbrMatrices();
  cout << "handling " << NbrBMatrices << " B matrices" << endl;
  int NbrParticles = MPSMatrix->GetMatrixNaturalNbrParticles(NbrFluxQuanta, true);
  cout << "Nbr of particles = " << NbrParticles << " " << ", Nbr of flux quanta=" << NbrFluxQuanta << endl;
  SparseRealMatrix* BMatrices = MPSMatrix->GetMatrices();

  cout << "B matrix size = " << BMatrices[0].GetNbrRow() << "x" << BMatrices[0].GetNbrColumn() << endl;
  
   RealVector CoefficientVector (NbrBMatrices,true);
 
   for (int i = 0; i <NbrBMatrices; i++)
   {
      CoefficientVector[i] = ((double) rand() / (RAND_MAX) - 0.5);
   }
   CoefficientVector.Normalize();

double lamba = 0.0;
double Newlamba = 100.0;

 cout <<"Start algorithm"<<endl;
 while ( fabs(lamba - Newlamba) > 1e-4  ) 
{
  SparseRealMatrix TmpMatrix =  BMatrices[0] * CoefficientVector[0];

  for (int i = 1; i<NbrBMatrices; i++)
  {
      TmpMatrix = TmpMatrix + (BMatrices[i] * CoefficientVector[i]);
  }
 

  SparseRealMatrix TmpMatrix2;
  TmpMatrix2.Copy(TmpMatrix);
  for (int i = 1; i< NbrParticles - 1; i++)
  {
      TmpMatrix2.Multiply(TmpMatrix);
  }

  for (int i = 0; i < NbrBMatrices; i++)
 {
   SparseRealMatrix TmpMatrix3 = Multiply(BMatrices[i],TmpMatrix2);
   CoefficientVector[i] = TmpMatrix3.Tr();
 }
  lamba = Newlamba ;
  Newlamba = CoefficientVector.Norm();
  CoefficientVector.Normalize();

  cout <<Newlamba << " "<<lamba<<endl;
}


  SparseRealMatrix TmpMatrix =  BMatrices[0] * CoefficientVector[0];

  for (int i = 1; i<NbrBMatrices; i++)
  {
      TmpMatrix = TmpMatrix + (BMatrices[i] * CoefficientVector[i]);
  }
 

  SparseRealMatrix TmpMatrix2;
  TmpMatrix2.Copy(TmpMatrix);
  for (int i = 1; i < NbrParticles ; i++)
  {
     TmpMatrix2.Multiply(TmpMatrix);
  }

 double FinalResult = TmpMatrix2.Tr();

 cout <<" FinalResult = "<<FinalResult <<endl;
 unsigned long *  PhysicalIndice = MPSMatrix->GetPhysicalIndices();

 for (int i = 0; i < NbrBMatrices;  i++)
   cout << PhysicalIndice[i] <<" " <<CoefficientVector[i]<<endl;

  return 0;
}

