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
  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-fusedblock", "number of blocks to be considered simultaneously", 1);
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
  int SizeBlock =  Manager.GetInteger("nbr-fusedblock");

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 

  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  int NbrParticles = MPSMatrix->GetMatrixNaturalNbrParticles(NbrFluxQuanta, true);
  int NbrBlock =  NbrParticles/SizeBlock;
  if ( (NbrParticles  % NbrBlock) != 0)
    {
      cout << "invalid number of flux quanta" << endl;
      return -1;
    }

  int NbrBMatrices = MPSMatrix->GetNbrMatrices();

  int NbrOrbitals = MPSMatrix->GetNbrOrbitals();
  int NbrStatesPerOrbital = MPSMatrix->GetMaximumOccupation() + 1;
  cout << " MPSMatrix->GetMaximumOccupation()  = "<< MPSMatrix->GetMaximumOccupation()<<endl;
  cout <<"NbrOrbitals = "<< NbrOrbitals<<endl;
  int DimensionPhysicalHilbertSpace = 1;

  for (int i = 1; i < NbrOrbitals * SizeBlock; i++)
     DimensionPhysicalHilbertSpace *= NbrStatesPerOrbital;

  cout << "handling " << NbrBMatrices << " B matrices" << endl;
  cout <<" Nbr of particles = " << NbrParticles << " " << ", Nbr of flux quanta=" << NbrFluxQuanta << endl;
  cout <<"Size block = "<< SizeBlock<< " , Nbr blocks = " << NbrBlock<<endl;
  cout <<"Dimension of the physical Hilbert space = " << DimensionPhysicalHilbertSpace <<endl;

  SparseRealMatrix* BMatrices = MPSMatrix->GetMatrices();

  cout << "B matrix size = " << BMatrices[0].GetNbrRow() << "x" << BMatrices[0].GetNbrColumn() << endl;
  
  SparseRealMatrix* FusedBMatrices = new SparseRealMatrix [DimensionPhysicalHilbertSpace];

  int TmpI;
  for(int i =0 ; i <DimensionPhysicalHilbertSpace; i++)
  { 
     TmpI = i,			
     FusedBMatrices[i].Copy(BMatrices[TmpI %  NbrStatesPerOrbital]);
     TmpI /=NbrStatesPerOrbital;
     for(int p = 1; p < SizeBlock - 1; p++)
	{
           FusedBMatrices[i].Multiply(BMatrices[TmpI %  NbrStatesPerOrbital]);
           TmpI /=NbrStatesPerOrbital;
	}
  }

   RealVector CoefficientVector (DimensionPhysicalHilbertSpace,true);
 
   for (int i = 0; i <DimensionPhysicalHilbertSpace; i++)
   {
      CoefficientVector[i] = ((double) rand() / (RAND_MAX) - 0.5);
   }
   CoefficientVector.Normalize();

double lamba = 0.0;
double Newlamba = 100.0;

 cout <<"Start algorithm"<<endl;

 while ( fabs(lamba - Newlamba) > 1e-8  ) 
{
  SparseRealMatrix TmpMatrix =  FusedBMatrices[0] * CoefficientVector[0];

  for (int i = 1; i<DimensionPhysicalHilbertSpace; i++)
  {
      TmpMatrix = TmpMatrix + (FusedBMatrices[i] * CoefficientVector[i]);
  }
 

  SparseRealMatrix TmpMatrix2;
  TmpMatrix2.Copy(TmpMatrix);
  for (int i = 1; i< NbrBlock - 1; i++)
  {
      TmpMatrix2.Multiply(TmpMatrix);
  }

  for (int i = 0; i < DimensionPhysicalHilbertSpace; i++)
 {
   SparseRealMatrix TmpMatrix3 = Multiply(FusedBMatrices[i],TmpMatrix2);
   CoefficientVector[i] = TmpMatrix3.Tr();
 }
  lamba = Newlamba ;
  Newlamba = CoefficientVector.Norm();
  CoefficientVector.Normalize();

  cout <<Newlamba << " "<<lamba<<endl;
}


  SparseRealMatrix TmpMatrix =  FusedBMatrices[0] * CoefficientVector[0];

  for (int i = 1; i<DimensionPhysicalHilbertSpace; i++)
  {
      TmpMatrix = TmpMatrix + (FusedBMatrices[i] * CoefficientVector[i]);
  }
 

  SparseRealMatrix TmpMatrix2;
  TmpMatrix2.Copy(TmpMatrix);
  for (int i = 1; i < NbrBlock ; i++)
  {
     TmpMatrix2.Multiply(TmpMatrix);
  }

 double FinalResult = TmpMatrix2.Tr();

 cout <<" FinalResult = "<<FinalResult <<endl;

 for (int i = 0; i < DimensionPhysicalHilbertSpace;  i++)
   cout << i <<" " <<CoefficientVector[i]<<endl;

  return 0;
}

