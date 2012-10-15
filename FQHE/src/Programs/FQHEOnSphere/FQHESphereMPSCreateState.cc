#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/BosonOnDiskWithSU2Spin.h"

#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/LongRational.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHEMPSCreateStateOperation.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Matrix/SparseComplexMatrix.h"

#include "Options/Options.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


void CreateLaughlinBMatrices (int laughlinIndex, SparseComplexMatrix* bMatrices,int pLevel, bool cylinderFlag, double kappa);
Complex CreateLaughlinAMatrixElement (int laughlinIndex, unsigned long* partition1, unsigned long* partition2, int p1Level, int p2Level, int nValue, FactorialCoefficient& coef);

void CreateMooreReadBMatrices (int mRIndex, SparseComplexMatrix* bMatrices,int pLevel, bool cylinderFlag, double kappa);

LongRational ComputeDescendantScalarProduct (long* partition, int partitionLength, int position, LongRational& centralCharge12, LongRational& weight);

LongRational ComputeDescendantMatrixElement (long* partition, int partitionLength, int descendantPosition, int position, 
					     LongRational& centralCharge12, LongRational& weight1, LongRational& weight2, 
					     LongRational& weight);

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereMPSCreateState" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

//   (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
//   (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta", 20);
//   (*SystemGroup) += new SingleIntegerOption  ('z', "lz-value", "twice the total lz value", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "file that describes the root configuration");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "p-truncation", "truncation level", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "laughlin-index", "index of the Laughlin state to generate", 3);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption  ('\n', "full-basis", "express the final vector in the full Haldane basis for the given root partition");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "precalculation-blocksize", " indicates the size of the block (i.e. number of B matrices) for precalculations", 1);
  (*OutputGroup) += new SingleStringOption ('o', "bin-output", "output the MPS state into a binary file");
  (*OutputGroup) += new SingleStringOption ('t', "txt-output", "output the MPS state into a text file");
  (*OutputGroup) += new BooleanOption ('n', "normalize-sphere", "express the MPS in the normalized sphere basis");
  (*OutputGroup) += new BooleanOption ('c', "normalize-cylinder", "express the MPS in the normalized cylinder basis");
  (*OutputGroup) += new SingleDoubleOption  ('r', "aspect-ratio", "aspect ratio of the cylinder", 1);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereMPSCreateState -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int NbrFluxQuanta = 0;
  int TotalLz = 0;
  int* ReferenceState = 0;
  char* OutputFileName = Manager.GetString("bin-output");
  char* OutputTxtFileName = Manager.GetString("txt-output");
  if ((OutputTxtFileName == 0) && (OutputFileName == 0))
    {
      cout << "error, an output file (binary or text) has to be provided" << endl;
      return 0;
    }

  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
    return -1;

  bool CylinderFlag = Manager.GetBoolean("normalize-cylinder");
  double AspectRatio = Manager.GetDouble("aspect-ratio");
  double kappa = 0.0;
  if (CylinderFlag)
    {
       kappa = (2.0 * M_PI)/sqrt(2.0 * M_PI * (NbrFluxQuanta + 1) * AspectRatio);
       cout<<"Cylinder geometry, kappa= "<<kappa<<endl;
    }


//  ParticleOnSphere* Space = 0;
  FermionOnSpherePTruncated* Space = 0;
  if (Manager.GetBoolean("boson") == true)
    {
      cout << "bosons are not yet implemented" << endl;
      return 0;
    }
  else
    {
#ifdef __64_BITS__
	  if (NbrFluxQuanta <= 62)
#else
	  if (NbrFluxQuanta <= 30)
#endif
	    {
	      Space = new FermionOnSpherePTruncated(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
	    }
	  else
	    {
#ifdef __128_BIT_LONGLONG__
	      if (NbrFluxQuanta <= 126)
#else
		if (NbrFluxQuanta <= 62)
#endif
		  {
		    Space = 0;//new FermionOnSpherePTruncatedLong(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
		  }
		else
		  {
#ifdef __128_BIT_LONGLONG__
		    cout << "cannot generate an Hilbert space when nbr-flux > 126" << endl;
#else
		    cout << "cannot generate an Hilbert space when nbr-flux > 62" << endl;
#endif
		    return 0;
		  }
	    }

   }

  cout << "Hilbert space dimension : " << Space->GetLargeHilbertSpaceDimension() << endl;

  int LaughlinIndex = Manager.GetInteger("laughlin-index");
  int NbrBMatrices = 2;

  
//   CreateMooreReadBMatrices (LaughlinIndex, 0, Manager.GetInteger("p-truncation"), CylinderFlag, kappa);
//   return 0;

  SparseComplexMatrix* SparseBMatrices = new SparseComplexMatrix[NbrBMatrices];
  SparseComplexMatrix* SparseConjugateBMatrices = new SparseComplexMatrix[NbrBMatrices];
  CreateLaughlinBMatrices (LaughlinIndex, SparseBMatrices, Manager.GetInteger("p-truncation"), CylinderFlag, kappa);

  return 0;


  for (int i = 0; i < NbrBMatrices; ++i)
    {
      SparseConjugateBMatrices[i] = SparseBMatrices[i].HermitianTranspose();
    }
  cout << "B matrix size = " << SparseBMatrices[0].GetNbrRow() << "x" << SparseBMatrices[0].GetNbrColumn() << endl;

  SparseComplexMatrix** SparseTensorProductBMatrices = new SparseComplexMatrix*[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      SparseTensorProductBMatrices[i] = new SparseComplexMatrix[NbrBMatrices];
      for (int j = 0; j < NbrBMatrices; ++j)
	SparseTensorProductBMatrices[i][j] = TensorProduct(SparseBMatrices[i], SparseConjugateBMatrices[j]);
    }
  
  RealVector State (Space->GetHilbertSpaceDimension(), true);


  FQHEMPSCreateStateOperation Operation(Space, SparseBMatrices, &State, Manager.GetInteger("p-truncation") + ((LaughlinIndex - 1) / 2),
					Manager.GetInteger("precalculation-blocksize"));
  Operation.ApplyOperation(Architecture.GetArchitecture());
  //    Space->CreateStateFromMPSDescription(SparseBMatrices, State, Manager.GetInteger("p-truncation") + ((LaughlinIndex - 1) / 2));

  if (Architecture.GetArchitecture()->CanWriteOnDisk() == true)
    {
      if (Manager.GetBoolean("normalize-sphere"))
	Space->ConvertFromUnnormalizedMonomial(State);
      if (CylinderFlag)
	State /= State.Norm();
      
      if (Manager.GetBoolean("full-basis") == true)  
	{
	  FermionOnSphereHaldaneBasis* SpaceHaldane = 0;
	  if (Manager.GetBoolean("boson") == true)
	    {
	      cout << "bosons are not yet implemented" << endl;
	      return 0;
	    }
	  else
	    {
#ifdef __64_BITS__
	      if (NbrFluxQuanta <= 62)
#else
		if (NbrFluxQuanta <= 30)
#endif
		  {
		    SpaceHaldane = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState);
		  }
		else
		  {
#ifdef __128_BIT_LONGLONG__
		    if (NbrFluxQuanta <= 126)
#else
		      if (NbrFluxQuanta <= 62)
#endif
			{
			  SpaceHaldane = 0;//new FermionOnSpherePTruncatedLong(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
			}
		      else
			{
#ifdef __128_BIT_LONGLONG__
			  cout << "cannot generate an Hilbert space when nbr-flux > 126" << endl;
#else
			  cout << "cannot generate an Hilbert space when nbr-flux > 62" << endl;
#endif
			  return 0;
			}
		  }
	      
	    }
	  
	  RealVector NewState;
	  NewState = Space->ConvertToHaldaneBasis(State, *SpaceHaldane);
	  
	  if (OutputTxtFileName != 0)
	    {
	      ofstream File;
	      File.open(OutputTxtFileName, ios::binary | ios::out);
	      File.precision(14);	
	      for (long i = 0; i < SpaceHaldane->GetLargeHilbertSpaceDimension(); ++i)
		{
		  NewState.PrintComponent(File, i) << " ";
		  Space->PrintStateMonomial(File, i) << endl;
		}
	      File.close();
	    }
	  if (OutputFileName != 0)
	    {
	      NewState.WriteVector(OutputFileName);
	    }
	}
      else 
	{ 
	  if (OutputTxtFileName != 0)
	    {
	      ofstream File;
	      File.open(OutputTxtFileName, ios::binary | ios::out);
	      File.precision(14);	
	      for (long i = 0; i < Space->GetLargeHilbertSpaceDimension(); ++i)
		{
		  State.PrintComponent(File, i) << " ";
		  Space->PrintStateMonomial(File, i) << endl;
		}
	      File.close();
	    }
	  if (OutputFileName != 0)
	    {
	      State.WriteVector(OutputFileName);
	    }
	}
    }
  
  return 0;
}

void CreateLaughlinBMatrices (int laughlinIndex, SparseComplexMatrix* bMatrices, int pLevel, bool cylinderFlag, double kappa)
{
  int NbrBMatrices = 2;
  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [pLevel + 1];
  ComplexMatrix* BMatrices = new ComplexMatrix[NbrBMatrices];
  for (int i = 0; i <= pLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort(i, i, pLevel + 1);
    }
  int* StartingIndexPerPLevel = new int [pLevel + 1];
  int* NbrIndicesPerPLevel = new int [pLevel + 1];
  StartingIndexPerPLevel[0] = 0;
  int NbrNValue = ((2 * pLevel) + laughlinIndex);
  int NValueShift = NbrNValue - 1;
  NbrIndicesPerPLevel[0] = U1BosonBasis[0]->GetHilbertSpaceDimension() * NbrNValue;
  for (int i = 1; i <= pLevel; ++i)
    {
      StartingIndexPerPLevel[i] = StartingIndexPerPLevel[i - 1] + NbrIndicesPerPLevel[i - 1];
      NbrIndicesPerPLevel[i] = U1BosonBasis[i]->GetHilbertSpaceDimension()  * NbrNValue;
    }
  int MatrixSize = NbrIndicesPerPLevel[pLevel] + StartingIndexPerPLevel[pLevel];

  BMatrices[0] = ComplexMatrix(MatrixSize, MatrixSize, true);
  for (int i = 0; i <= pLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace = U1BosonBasis[i];
      for (int j = 1; j < NbrNValue; ++j)
	{
	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
	    {
              int N1 = (j - NValueShift/2);
	      Complex Tmp (1.0, 0.0);
              if (cylinderFlag)
                Tmp *= exp(-kappa*kappa*(i + (N1 - 1) * (N1 - 1)/(4.0 * laughlinIndex)+ (N1 * N1)/(4.0 * laughlinIndex)));
	      BMatrices[0].SetMatrixElement(StartingIndexPerPLevel[i] + ((k * NbrNValue) + j - 1), StartingIndexPerPLevel[i] + ((k * NbrNValue) + j), Tmp);
	    }
	}
    }
   
  BMatrices[1] = ComplexMatrix(MatrixSize, MatrixSize, true);
  unsigned long* Partition1 = new unsigned long [pLevel + 2];
  unsigned long* Partition2 = new unsigned long [pLevel + 2];
  FactorialCoefficient Coef;

  for (int i = 0; i <= pLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace1 = U1BosonBasis[i];
      //int MaxN1 = (2 * i) + laughlinIndex;
      for (int j = 0; j <= pLevel; ++j)
	{
	  BosonOnDiskShort* TmpSpace2 = U1BosonBasis[j];
	  //int MaxN2 = (2 * 2) + laughlinIndex;
	  int N1 = (2 * (j - i) + laughlinIndex - 1 + NValueShift) / 2;
	  int N2 = (2 * (j - i) - laughlinIndex + 1 + NValueShift) / 2;
	  for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
	    {
	      TmpSpace1->GetOccupationNumber(k1, Partition1);
	      for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
		{
		  TmpSpace2->GetOccupationNumber(k2, Partition2);
		  Complex Tmp = CreateLaughlinAMatrixElement(laughlinIndex, Partition1, Partition2, i, j, - (N1 + N2 - NValueShift) / 2, Coef);
                  if (cylinderFlag)
                    Tmp *= exp(-kappa*kappa*(0.5 * i + 0.5 * j + pow(N1 - NValueShift/2,2.0)/(4.0 * laughlinIndex) + pow(N2 - NValueShift/2,2.0)/(4.0 * laughlinIndex)));
		  BMatrices[1].SetMatrixElement(StartingIndexPerPLevel[i] + ((k1 * NbrNValue) + N1), StartingIndexPerPLevel[j] + ((k2 * NbrNValue) + N2), Tmp);
//		  cout << i << " " << j << " | " << k1 << " " << k2 << " | " << N1 << " " << N2 << " " << (-(N1 + N2 - NValueShift) / 2) << " : " << Tmp << endl;
		}
	    }
	}
    }

  delete[] Partition1;
  delete[] Partition2;

//   HermitianMatrix TmpDensityMatrix(BMatrices[0].GetNbrRow(), true);
//   HermitianMatrix TmpDensityMatrix2(BMatrices[0].GetNbrRow(), true);
//   int OrbitalCut = 25;
//   TmpDensityMatrix.SetToIdentity();
//   for (int i = 0; i < 50; ++i)
//     {
//       cout << "iteration " << i << endl;
//       TmpDensityMatrix2.ClearMatrix();
//       TmpDensityMatrix2 += TmpDensityMatrix.Conjugate(BMatrices[0]);
//       TmpDensityMatrix2 += TmpDensityMatrix.Conjugate(BMatrices[1]);
//       TmpDensityMatrix.Copy(TmpDensityMatrix2);      
//     }
//   RealDiagonalMatrix TmpDiag;
//   ComplexMatrix TmpBasis (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow(), true);
//   TmpBasis.SetToIdentity();
// #ifdef __LAPACK__
//   TmpDensityMatrix2.LapackDiagonalize(TmpDiag, TmpBasis);
// #else
//   TmpDensityMatrix2.Diagonalize(TmpDiag, TmpBasis);
// #endif
//   int Count = TmpDiag.GetNbrRow();
//   for (int n = 0; n < TmpDiag.GetNbrRow(); ++n)
//     {
//       if (fabs(TmpDiag(n, n)) > 1e-10)
// 	++Count;
//       cout << (TmpDiag(n, n)  / TmpDensityMatrix.Tr()) << " ";
//     }
//   cout << endl;

//   ComplexMatrix TmpBasis2 (BMatrices[0].GetNbrRow(), Count, true);
//   Count = 0;
//   for (int n = 0; n < TmpDiag.GetNbrRow(); ++n)
//     {
//       if (fabs(TmpDiag(n, n)) > 1e-10)
// 	{
// 	  TmpBasis2[Count].Copy(TmpBasis[n]);
// 	  ++Count;
// 	}
//     }


//   TmpDensityMatrix.SetToIdentity();
//   for (int i = 0; i < OrbitalCut; ++i)
//     {
//       cout << "iteration " << i << endl;
//       TmpDensityMatrix2.ClearMatrix();
//       TmpDensityMatrix2 += TmpDensityMatrix.InvConjugate(BMatrices[0]);
//       TmpDensityMatrix2 += TmpDensityMatrix.InvConjugate(BMatrices[1]);
//       TmpDensityMatrix.Copy(TmpDensityMatrix2);      
//     }
//   HermitianMatrix TmpDensityMatrix3 = TmpDensityMatrix.Conjugate(TmpBasis2);
  


  for (int i = 0; i < NbrBMatrices; ++i)
    {
      bMatrices[i] = BMatrices[i];
    }
  delete[] BMatrices;
}


Complex CreateLaughlinAMatrixElement (int laughlinIndex, unsigned long* partition1, unsigned long* partition2, int p1Level, int p2Level, int nValue, FactorialCoefficient& coef)
{
  Complex Tmp = 1.0;
  if (nValue != (p1Level - p2Level))
    {
      Tmp = 0.0;
      return Tmp;
    }
  int PMax = p1Level;
  if (p2Level > p1Level)
    PMax = p2Level;
  for (int i = 1; i <= PMax; ++i)
    {
      Complex Tmp2 = 0.0;
//      cout << partition1[i] << " " << partition2[i] << " " << i << " " << PMax << endl;
      for (int j = 0; j <= partition1[i]; ++j)
	{
	  int k = partition2[i] + j - partition1[i];
	  if ((k >= 0) && (k <= partition2[i]))
	    {
	      int Sum = k + j;
	      coef.SetToOne();
	      coef.PartialFactorialMultiply(partition1[i] - j + 1, partition1[i]);
	      coef.PartialFactorialMultiply(partition2[i] - k + 1, partition2[i]);
	      coef.FactorialDivide(j);
	      coef.FactorialDivide(k);
	      coef.FactorialDivide(j);
	      coef.FactorialDivide(k);
	      coef.PowerNMultiply(laughlinIndex, Sum);
	      coef.PowerNDivide(i, Sum);
//	      cout << "Sum=" << Sum << endl;
	      switch  (Sum & 0x3)
		{
		case 0:
		  Tmp2.Re += sqrt(coef.GetNumericalValue());
		  break;
		case 1:
		  Tmp2.Im += sqrt(coef.GetNumericalValue());
		  break;
		case 2:
		  Tmp2.Re -= sqrt(coef.GetNumericalValue());
		  break;
		case 3:
		  Tmp2.Im -= sqrt(coef.GetNumericalValue());
		  break;
		}
	    }
	}
//      cout << Tmp2 << endl;
      Tmp *= Tmp2;
    }
  return Tmp;
}

void CreateMooreReadBMatrices (int mRIndex, SparseComplexMatrix* bMatrices,int pLevel, bool cylinderFlag, double kappa)
{

//   LongRational CentralCharge12 (-3l, 5l);

  LongRational CentralCharge12 (1l, 2l);
  CentralCharge12 /= 12l;
  LongRational WeightIdentity (0l, 1l);
  LongRational WeightPsi (1l, 2l);
  long* Partition = new long[2 * (pLevel + 1)];
  unsigned long* TmpPartition = new unsigned long [pLevel + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [pLevel + 1];
  RealSymmetricMatrix* ScalarProductIdentity = new RealSymmetricMatrix[pLevel + 1];
  RealSymmetricMatrix* ScalarProductPsi = new RealSymmetricMatrix[pLevel + 1];
  RealMatrix** MatrixPsi = new RealMatrix*[pLevel + 1];
  RealMatrix* OrthogonalBasisIdentity = new RealMatrix[pLevel + 1];
  RealMatrix* OrthogonalBasisPsi = new RealMatrix[pLevel + 1];

  for (int i = 0; i <= pLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort (i, i, pLevel + 1);
      MatrixPsi[i] = new RealMatrix[pLevel + 1];
    }
  
  for (int i = 0; i <= pLevel; ++i)
    {
      cout << "Level = " <<  i << endl;
      ScalarProductIdentity[i] = RealSymmetricMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      ScalarProductPsi[i] = RealSymmetricMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	for (int m = n; m < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++m)
	  {
	    int PartitionLength = 0;
	    U1BosonBasis[i]->GetOccupationNumber(n, TmpPartition);	    
	    for (int k = 1; k <= i; ++k)
	      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		{
		  ++PartitionLength;
		}
	    int Position = PartitionLength;
	    PartitionLength = 0;
	    for (int k = 1; k <= i; ++k)
	      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		{
		   Partition[Position - PartitionLength - 1] = (long) k;
		  ++PartitionLength;
		}
	    U1BosonBasis[i]->GetOccupationNumber(m, TmpPartition);	    
	    for (int k = 1; k <= i; ++k)
	      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		{
		  Partition[PartitionLength] = -(long) k;
		  ++PartitionLength;		  
		}
// 	    cout << "< ";
// 	    for (int k = 0; k < Position; ++k)
// 	      cout << Partition[k] << " ";
//  	    cout << "| ";
//  	    for (int k = Position; k < PartitionLength; ++k)
//  	      cout << Partition[k] << " ";
//  	    cout << "> = ";	    
	    LongRational Tmp = ComputeDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightIdentity);
//	    cout << Tmp << endl;
	    ScalarProductIdentity[i].SetMatrixElement(m, n, Tmp.GetNumericalValue());
	    Tmp = ComputeDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightPsi);
	    ScalarProductPsi[i].SetMatrixElement(m, n, Tmp.GetNumericalValue());
	  }
      
      RealSymmetricMatrix TmpMatrix;
      TmpMatrix.Copy(ScalarProductIdentity[i]);
      RealMatrix TmpBasis(U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension());
      TmpBasis.SetToIdentity();
      RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
      TmpMatrix.LapackDiagonalize(TmpDiag, TmpBasis);
#else
      TmpMatrix.Diagonalize(TmpDiag, TmpBasis);
#endif
//       for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
// 	cout << TmpDiag(n, n) << " ";
      double Error = 0.0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) > Error)
	  Error = fabs(TmpDiag(n, n));
      Error *= 1e-14;
      if (Error < 1e-14)
	Error = 1e-14;
      int Count  = 0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) < Error)
	  ++Count;
      cout << "nbr of null vectors identity sector = " << Count << endl;
      if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	{
	  OrthogonalBasisIdentity[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisIdentity[i][Count].Copy(TmpBasis[n]);
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisIdentity[i] = RealMatrix();
	}

      TmpMatrix.Copy(ScalarProductPsi[i]);
      TmpBasis.SetToIdentity();
#ifdef __LAPACK__
      TmpMatrix.LapackDiagonalize(TmpDiag, TmpBasis);
#else
      TmpMatrix.Diagonalize(TmpDiag, TmpBasis);
#endif
      Error = 0.0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) > Error)
	  Error = fabs(TmpDiag(n, n));
      Error *= 1e-14;
      if (Error < 1e-14)
	Error = 1e-14;
      Count  = 0;
      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	if (fabs(TmpDiag(n, n)) < Error)
	  ++Count;
      cout << "nbr of null vectors Psi sector = " << Count << endl;
      if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	{
	  OrthogonalBasisPsi[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisPsi[i][Count].Copy(TmpBasis[n]);
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisPsi[i] = RealMatrix();
	}
      cout << "---------------------------------" << endl;
    }

  cout << "computing Psi matrix elements" << endl;
  LongRational Weight (WeightPsi);
  for (int i = 0; i <= pLevel; ++i)
    {
      for (int j = 0; j <= pLevel; ++j)
	{
	  MatrixPsi[i][j] = RealMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  cout << "Levels = " <<  i << " " << j << endl;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    for (int m = n; m < U1BosonBasis[j]->GetHilbertSpaceDimension(); ++m)
	  {
	    int PartitionLength = 0;
	    U1BosonBasis[i]->GetOccupationNumber(n, TmpPartition);	    
	    for (int k = 1; k <= i; ++k)
	      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		{
		  ++PartitionLength;
		}
	    int Position = PartitionLength;
	    PartitionLength = 0;
	    for (int k = 1; k <= i; ++k)
	      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		{
		   Partition[Position - PartitionLength - 1] = (long) k;
		  ++PartitionLength;
		}
	    U1BosonBasis[j]->GetOccupationNumber(m, TmpPartition);	    
	    for (int k = 1; k <= j; ++k)
	      for (unsigned long  l = 0ul; l < TmpPartition[k]; ++l)
		{
		  Partition[PartitionLength] = -(long) k;
		  ++PartitionLength;		  
		}
// 	    cout << "< ";
// 	    for (int k = 0; k < Position; ++k)
// 	      cout << Partition[k] << " ";
//  	    cout << "| ";
//  	    for (int k = Position; k < PartitionLength; ++k)
//  	      cout << Partition[k] << " ";
//  	    cout << "> = ";	    
	    LongRational Tmp = ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightIdentity, WeightPsi, Weight);
//	    cout << Tmp << endl;
	    MatrixPsi[i][j].SetMatrixElement(m, n, Tmp.GetNumericalValue());
	  }
	}
    }


  int NbrBMatrices = 2;
  ComplexMatrix* BMatrices = new ComplexMatrix[NbrBMatrices];

  int* StartingIndexPerPLevel = new int [pLevel + 1];
  int* NbrIndicesPerPLevel = new int [pLevel + 1];
  StartingIndexPerPLevel[0] = 0;
  int NbrNValue = (((2 * pLevel) + mRIndex));
  int NValueShift = NbrNValue - 1;
  NbrIndicesPerPLevel[0] = (U1BosonBasis[0]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentity[0].GetNbrColumn() + OrthogonalBasisPsi[0].GetNbrColumn())) * NbrNValue;
  for (int i = 1; i <= pLevel; ++i)
    {
      int Tmp = 0;
      for (int j = 0; j <= i; ++j)
	{
	  cout << i << " " << j << endl;
	  cout << OrthogonalBasisIdentity[j - i].GetNbrColumn() << endl;
	  cout << OrthogonalBasisPsi[j - i].GetNbrColumn() << endl;
	  Tmp += U1BosonBasis[j]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentity[j - i].GetNbrColumn() + OrthogonalBasisPsi[j - i].GetNbrColumn());
	}
      StartingIndexPerPLevel[i] = StartingIndexPerPLevel[i - 1] + NbrIndicesPerPLevel[i - 1];
      NbrIndicesPerPLevel[i] =  Tmp * NbrNValue;
    }
  int MatrixSize = NbrIndicesPerPLevel[pLevel] + StartingIndexPerPLevel[pLevel];
  
  BMatrices[0] = ComplexMatrix(MatrixSize, MatrixSize, true);
  for (int i = 0; i <= pLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[p - i];
	  for (int j = 1; j < NbrNValue; ++j)
	    {
	      for (int k = 0; k < TmpSpaceCharged->GetHilbertSpaceDimension(); ++k)
		{
		  Complex Tmp (0.0, 0.0);		   
		  BMatrices[0].SetMatrixElement(StartingIndexPerPLevel[i] + ((k * NbrNValue) + j - 1), StartingIndexPerPLevel[i] + ((k * NbrNValue) + j), Tmp);
		}
	    }
	}
    }
  
  BMatrices[1] = ComplexMatrix(MatrixSize, MatrixSize, true);
}


LongRational ComputeDescendantScalarProduct (long* partition, int partitionLength, int position, LongRational& centralCharge12, LongRational& weight)
{
  if (partitionLength == 0)
    {
      return 1l;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
  if ((position == partitionLength) || (position == 0))
    return 0l;
  if (partitionLength == 2)
    {
      if (partition[0] != -partition[1])
	{
	  return 0l;
	}
      else
	{
	  LongRational Tmp1 (centralCharge12);
	  Tmp1 *= partition[0] * (partition[0] * partition[0] - 1l);
	  LongRational Tmp2 (weight);
	  Tmp2 *= 2l * partition[0];
	  Tmp1 += Tmp2;
	  return Tmp1;
	}
    }
  LongRational Tmp(0l);
  if ((partition[position - 1] + partition[position]) == 0)
    {
      long TmpLength = 0l;
      long Store = partition[position - 1];
      for (int i = position + 1; i < partitionLength; ++i)
	TmpLength += partition[i];
      for (int i = position + 1; i < partitionLength; ++i)
	partition[i - 2] = partition[i];
      Tmp += ((((Store * (Store * Store - 1l)) * centralCharge12)
	       + (2l * Store) * (weight - TmpLength)) * 
	      ComputeDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, weight));
      for (int i = partitionLength - 1; i > position; --i)
	partition[i] = partition[i - 2];
      partition[position - 1] = Store;
      partition[position] = -Store;
    }
  else
    {
      long Store1 = partition[position - 1];
      long Store2 = partition[position];
      partition[position - 1] += partition[position];
      for (int i = position + 1; i < partitionLength; ++i)
	partition[i - 1] = partition[i];
      if ((Store1 + Store2) > 0)
	{
	  Tmp += ((Store1 - Store2) 
		  * ComputeDescendantScalarProduct(partition, partitionLength - 1, position, centralCharge12, weight));
	}
      else
	{
	  Tmp += ((Store1 - Store2) 
		  * ComputeDescendantScalarProduct(partition, partitionLength - 1, position - 1, centralCharge12, weight));
	}
      for (int i = partitionLength - 1; i > position; --i)
	partition[i] = partition[i - 1];
      partition[position] = Store2;
      partition[position - 1] = Store1;
    }

  long Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  Tmp += ComputeDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, weight);
  Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  return Tmp;
}


LongRational ComputeDescendantMatrixElement (long* partition, int partitionLength, int descendantPosition, int position, 
					     LongRational& centralCharge12, LongRational& weight1, LongRational& weight2, 
					     LongRational& weight)
{
  if (partitionLength == 0)
    {
      return 1l;
    }
  if ((descendantPosition < partitionLength) && (partition[partitionLength - 1] > 0))
    {
      return 0l;
    }
  while ((position > 0) && (partition[position - 1] < 0l))
    --position;
//   for (int i = 0; i < partitionLength; ++i)
//     cout << partition[i] << " ";
//   cout << "| " << partitionLength << " " << descendantPosition << " " << position << endl;
  if (descendantPosition == partitionLength) 
    {
      LongRational Tmp(0l);
      LongRational TmpSum = weight1;
      TmpSum -= weight2;
      LongRational Tmp2;
      for (int i = 0; i < partitionLength ; ++i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  Tmp2 += TmpSum;
	  Tmp += Tmp2;
	}
      return Tmp;
    }
  if (position == 0)
    {
      LongRational Tmp(0l);
      LongRational TmpSum = weight1;
      TmpSum -= weight2;
      LongRational Tmp2;
      for (int i = 0; i < partitionLength ; ++i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  Tmp2 += TmpSum;
	  Tmp += Tmp2;
	}
      Tmp *= 1l - (2l * (partitionLength & 1l));
      return Tmp;
    }
  if (descendantPosition < position)
    {
      LongRational Tmp(0l);
      if ((partition[position - 1] + partition[position]) == 0)
	{
	  long TmpLength = 0l;
	  long Store = partition[position - 1];
	  for (int i = position + 1; i < partitionLength; ++i)
	    TmpLength += partition[i];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 2] = partition[i];
	  Tmp += ((((Store * (Store * Store - 1l)) * centralCharge12)
		   + (2l * Store) * (weight2 - TmpLength)) * 
		  ComputeDescendantMatrixElement(partition, partitionLength - 2, descendantPosition, position - 1, centralCharge12, weight1, weight2, weight));
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 2];
	  partition[position - 1] = Store;
	  partition[position] = -Store;
	}
      else
	{
	  long Store1 = partition[position - 1];
	  long Store2 = partition[position];
	  partition[position - 1] += partition[position];
	  for (int i = position + 1; i < partitionLength; ++i)
	    partition[i - 1] = partition[i];
	  if ((Store1 + Store2) > 0)
	    {
	      Tmp += ((Store1 - Store2) 
		      * ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, position, centralCharge12, weight1, weight2, weight));
	    }
	  else
	    {
	      Tmp += ((Store1 - Store2) 
		      * ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition, position - 1, centralCharge12, weight1, weight2, weight));
	    }
	  for (int i = partitionLength - 1; i > position; --i)
	    partition[i] = partition[i - 1];
	  partition[position] = Store2;
	  partition[position - 1] = Store1;
	}
      
      long Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      Tmp += ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition, position + 1, centralCharge12, weight1, weight2, weight);
      Store1 = partition[position - 1];
      partition[position - 1] = partition[position];
      partition[position] = Store1;
      return Tmp;  
    }

  LongRational Tmp1 = ComputeDescendantMatrixElement(partition, partitionLength, descendantPosition - 1, position, 
						     centralCharge12, weight1, weight2, weight);
  LongRational Tmp2 = weight;
  Tmp2 *= partition[position - 1];
  Tmp2 += weight1;
  Tmp2 -= weight2;
  long TmpLength = 0l;
  for (int i = position; i < partitionLength; ++i)
    TmpLength -= partition[i];
  for (int i = 0; i < (position - 1); ++i)
    TmpLength += partition[i];
  Tmp2 += TmpLength;
  long Store = partition[position - 1];
  for (int i = position; i < partitionLength; ++i)
    partition[i - 1] = partition[i];
  Tmp2 *= ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition - 1, position - 1, 
					 centralCharge12, weight1, weight2, weight);
  for (int i = partitionLength - 1; i > position; --i)
    partition[i] = partition[i - 1];
  partition[position - 1] = Store;
  Tmp1 += Tmp2;
  return Tmp1;
}
