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

#include "Tools/FQHEMPS/FQHEMPODensityOperator.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Hamiltonian/TensorProductSparseMatrixHamiltonian.h"
#include "Hamiltonian/TensorProductSparseMatrixSelectedBlockHamiltonian.h"
#include "Hamiltonian/TripleTensorProductSparseMatrixHamiltonian.h"

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
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

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
  
  OptionManager Manager ("FQHECylinderMPSEdgeSpectrum" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  
  ArchitectureManager Architecture;
   FQHEMPSMatrixManager MPSMatrixManager (false, false);

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  OptionGroup* ArnoldiGroup  = new OptionGroup ("Arnoldi options");
  Architecture.AddOptionGroup(&Manager);
  Manager += ArnoldiGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new BooleanOption ('\n', "infinite-cylinder", "evaluate the entanglement spectrum on the infinite cylinder");  
  (*SystemGroup) += new SingleStringOption  ('\n', "left-eigenstate", "file containing the transfer matrix left eigenstate");
  (*SystemGroup) += new SingleStringOption  ('\n', "right-eigenstate", "file containing the transfer matrix right eigenstate");
  (*SystemGroup) += new BooleanOption  ('\n',"use-padding","use-padding");
  (*SystemGroup) += new SingleStringOption ('\n', "left-interaction", "file describing the confining potential of the left hand side of the cylinder (optional)");
  (*SystemGroup) += new SingleStringOption ('\n', "right-interaction", "file describing the confining potential of the right hand side of the cylinder");
  (*SystemGroup) += new SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");  


  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHECylinderMPSEdgeSpectrum -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int LeftNbrFluxQuanta = -1;
  double* LeftOneBodyPotential = 0;
  int RightNbrFluxQuanta = 0;
  double* RightOneBodyPotential = 0;
  if (Manager.GetString("right-interaction") == 0)
    {
      cout << "error, the right interaction description should be provided" << endl;
      return -1;
    }
  if (FQHESphereGetOneBodyPotentials(Manager.GetString("right-interaction"), RightNbrFluxQuanta, RightOneBodyPotential) == false)
    {
      return -1;
    }
  if (Manager.GetString("left-interaction") != 0)
    {
      if (FQHESphereGetOneBodyPotentials(Manager.GetString("left-interaction"), LeftNbrFluxQuanta, LeftOneBodyPotential) == false)
	{
	  return -1;
	}
    }
  int NbrFluxQuanta = (LeftNbrFluxQuanta + RightNbrFluxQuanta) + 1;

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }
  int NbrBMatrices = MPSMatrix->GetNbrMatrices();
  SparseRealMatrix* BMatrices = MPSMatrix->GetMatrices();
  double Normalisation = 1.0;

  if (Manager.GetBoolean("infinite-cylinder") == false)
    {
      double* TotalOneBodyPotentials = new double[NbrFluxQuanta + 1];
      if (Manager.GetString("left-interaction") != 0)
	{
	  for (int i = 0; i <= LeftNbrFluxQuanta; ++i)
	    {
	      TotalOneBodyPotentials[i] = LeftOneBodyPotential[i];
	    }
	}
      for (int i = 0; i <= RightNbrFluxQuanta; ++i)
	{
	  TotalOneBodyPotentials[i + LeftNbrFluxQuanta + 1] = RightOneBodyPotential[i];
	}

      cout << "Number of flux quanta = " << NbrFluxQuanta << endl;
      int MPSRowIndex = 0;
      int MPSColumnIndex = 0;
      MPSMatrix->GetMatrixBoundaryIndices(MPSRowIndex, MPSColumnIndex, Manager.GetBoolean("use-padding"));
      
      RealVector TmpVectorEMatrix (BMatrices[0].GetNbrRow() * BMatrices[0].GetNbrRow(), true);
      RealVector TmpVectorEMatrix2 (BMatrices[0].GetNbrRow() * BMatrices[0].GetNbrRow(), true);
      double* Coefficients = new double[NbrBMatrices];
      for(int i= 0; i < NbrBMatrices; i++)
	Coefficients[i] = 1.0;
      TmpVectorEMatrix[MPSRowIndex + (BMatrices[0].GetNbrRow() * MPSRowIndex)] = 1.0;

//       SparseRealMatrix* TmpMPOMatrices = new SparseRealMatrix [NbrBMatrices];
//       for (int i = 0; i < NbrBMatrices; ++i)
// 	{
// 	  TmpMPOMatrices[i] = SparseRealMatrix(1, 1);
// 	  TmpMPOMatrices[i].SetMatrixElement(0, 0, 1.0);
// 	}

      TensorProductSparseMatrixHamiltonian* TransferMatrix = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, BMatrices, BMatrices, Coefficients,
												      Architecture.GetArchitecture()); 
//       TripleTensorProductSparseMatrixHamiltonian* TransferMatrix = new TripleTensorProductSparseMatrixHamiltonian(NbrBMatrices, BMatrices, TmpMPOMatrices, BMatrices, Coefficients,
// 														  Architecture.GetArchitecture()); 
      for (int i = 0; i <= NbrFluxQuanta; ++i)
	{
	  TransferMatrix->LowLevelMultiply(TmpVectorEMatrix, TmpVectorEMatrix2);
//	  cout << TmpVectorEMatrix2 << "---------" << endl;
	  RealVector TmpVectorEMatrix3 = TmpVectorEMatrix;
	  TmpVectorEMatrix = TmpVectorEMatrix2;
	  TmpVectorEMatrix2 = TmpVectorEMatrix3;
	}      
      Normalisation = 1.0 / TmpVectorEMatrix[MPSColumnIndex + (BMatrices[0].GetNbrColumn() * MPSColumnIndex)];

      cout << "Normalisation = " << Normalisation << endl;
      
//       for (int i = 0; i < NbrBMatrices; ++i)
//  	{
//  	  BMatrices[i] = SparseRealMatrix(1, 1);
//  	  BMatrices[i].SetMatrixElement(0, 0, 1.0);
//  	}
//       MPSRowIndex = 0;
//       MPSColumnIndex = 0;

      FQHEMPODensityOperator TmpMPO (MPSMatrix->GetMaximumOccupation(), 1.0);
      int MPORowIndex = 0;
      int MPOColumnIndex = 0;
      TmpMPO.GetMatrixBoundaryIndices(MPORowIndex, MPOColumnIndex);
      SparseRealMatrix* TmpMPOMatrices = TmpMPO.GetMatrices();
      cout << "MPSRowIndex=" << MPSRowIndex << " MPORowIndex=" << MPORowIndex << " MPSColumnIndex=" << MPSColumnIndex << " MPOColumnIndex=" << MPOColumnIndex << endl;
      TmpVectorEMatrix = RealVector(BMatrices[0].GetNbrRow() * BMatrices[0].GetNbrRow() * TmpMPO.GetBondDimension(), true);
      TmpVectorEMatrix2 = RealVector(BMatrices[0].GetNbrRow() * BMatrices[0].GetNbrRow() * TmpMPO.GetBondDimension(), true);
      TmpVectorEMatrix[MPSColumnIndex + (((TmpMPOMatrices[0].GetNbrColumn() * MPSColumnIndex) + MPOColumnIndex) * BMatrices[0].GetNbrColumn())] = 1.0;
//       for (int i = 0; i < NbrBMatrices; ++i)
// 	{
// 	  cout << TmpMPOMatrices[i] << endl;
// 	}
      for (int i = 0; i <= NbrFluxQuanta; ++i)
	{
	  TmpMPO.SetPrefactor(TotalOneBodyPotentials[i]);
	  TripleTensorProductSparseMatrixHamiltonian* MPOTransferMatrix = new TripleTensorProductSparseMatrixHamiltonian(NbrBMatrices, BMatrices, TmpMPOMatrices, BMatrices, Coefficients,
															 Architecture.GetArchitecture()); 
	  MPOTransferMatrix->LowLevelMultiply(TmpVectorEMatrix, TmpVectorEMatrix2);
 	  RealVector TmpVectorEMatrix3 = TmpVectorEMatrix;
 	  TmpVectorEMatrix = TmpVectorEMatrix2;
 	  TmpVectorEMatrix2 = TmpVectorEMatrix3;	  
 	}
      double TmpNbrParticles = TmpVectorEMatrix[MPSRowIndex + (((TmpMPOMatrices[0].GetNbrRow() * MPSRowIndex) + MPORowIndex) * BMatrices[0].GetNbrRow())];
      cout << (TmpNbrParticles * Normalisation) << endl;
    }

//   cout << "B matrix size = " << BMatrices[0].GetNbrRow() << "x" << BMatrices[0].GetNbrColumn() << endl;
//   unsigned long * ArrayPhysicalIndice = MPSMatrix->GetPhysicalIndices();
  
  
  return 0;
}

