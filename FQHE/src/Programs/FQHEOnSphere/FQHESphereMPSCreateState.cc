#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnSpherePTruncated.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/BosonOnDiskWithSU2Spin.h"

#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/LongRational.h"
#include "MathTools/IntegerAlgebraTools.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/FQHEMPSCreateStateOperation.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Matrix/SparseComplexMatrix.h"
#include "Matrix/SparseRealMatrix.h"

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


// create the B matrices for the laughlin state
//
// laughlinIndex = power of the Laughlin part 
// bMatrices = array where the B matrices will be stored
// nbrBMatrices = number of B matrices to compute (max occupation + 1)
// pLevel = |P| level truncation
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
void CreateLaughlinBMatrices (int laughlinIndex, SparseRealMatrix* bMatrices, int nbrBMatrices, int pLevel, bool cylinderFlag, double kappa);

// create the matrix element of the B matrix U(1) part
//
// chargeNumerator = numerator of the charge (in sqrt(q) unit)
// chargeDenominator = denominator of the charge (in sqrt(q) unit)
// partition1 = U(1) partition associated to the left state
// p1Level = length of partition1
// partition2 = U(1) partition associated to the left state
// p1Level = length of partition2
// coef = reference on a temporary factorial coefficient
// return value = matrix element
double CreateLaughlinAMatrixElement (int chargeNumerator, int chargeDenominator, unsigned long* partition1, unsigned long* partition2, int p1Level, int p2Level, FactorialCoefficient& coef);

// create the B matrices for the (k,r) clustered with k=2
//
// laughlinIndex = power of the Laughlin part 
// rIndex = r value
// bMatrices = array where the B matrices will be stored
// pLevel = |P| level truncation
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
void CreateK2ClusteredStateBMatrices (int laughlinIndex, int rIndex, SparseRealMatrix* bMatrices,int pLevel, bool cylinderFlag, double kappa);

// 
LongRational ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, LongRational& centralCharge12, LongRational& weight);

LongRational ComputeDescendantMatrixElement (long* partition, int partitionLength, int descendantPosition, int position, 
					     LongRational& centralCharge12, LongRational& weight1, LongRational& weight2, 
					     LongRational& weight);

// compute the linearized index of the B matrix for the k=2 case
//
// return value = linearized index
int GetK2MatrixIndex(int charge, int chargedPartitionIndex, int nbrCharges, int chargeSectorDimension, int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int globalIndexShift);

// create the B matrices for the k=3 Read-Rezayi state
//
// laughlinIndex = power of the Laughlin part 
// bMatrices = array where the B matrices will be stored
// pLevel = |P| level truncation
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
void CreateReadRezayiK3StateBMatrices (int laughlinIndex, SparseRealMatrix* bMatrices,int pLevel, bool cylinderFlag, double kappa);


// compute the linearized index of the B matrix for the k=3 Read-Rezayi state
//
// return value = linearized index
int GetReadRezayiK3MatrixIndex(int charge, int chargedPartitionIndex, int nbrCharges, int chargeSectorDimension, int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int nbrPsiDescendant, int globalIndexShift);


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
  (*SystemGroup) += new BooleanOption  ('\n', "k-2", "consider the (k=2,r) series of clustered states");
 (*SystemGroup) += new BooleanOption  ('\n', "rr-3", "consider the k= 3 Read-Rezayi state");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "r-index", "r index of the (k,r) clustered state", 2);
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


  ParticleOnSphere* Space = 0;
  if (Manager.GetBoolean("boson") == true)
    {
      Space = new BosonOnSpherePTruncated(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
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
  int RIndex = Manager.GetInteger("r-index");
  int NbrBMatrices = 2;
  int MatrixElement = Manager.GetInteger("p-truncation") + ((LaughlinIndex - 1) / 2);
  if (Manager.GetBoolean("k-2") == true)
    {
      if ((Manager.GetInteger("r-index") & 1) == 0)
	MatrixElement = Manager.GetInteger("p-truncation") + (Manager.GetInteger("r-index") / 2);
      else
	MatrixElement = 2 * Manager.GetInteger("p-truncation") + Manager.GetInteger("r-index") - 1;
    }
  if (Manager.GetBoolean("rr-3") == true)
    {
      MatrixElement = 3 * (Manager.GetInteger("p-truncation") + 1);
    }
  if (Manager.GetBoolean("boson"))
    {
      NbrBMatrices = NbrParticles + 1;
    }

//   CreateK2ClusteredStateBMatrices (2, RIndex, 0, Manager.GetInteger("p-truncation"), CylinderFlag, kappa);
//   return 0;
  
  SparseRealMatrix* SparseBMatrices = new SparseRealMatrix[NbrBMatrices];
  if (Manager.GetBoolean("k-2") == false)
    {
      if (Manager.GetBoolean("rr-3") == false)
	CreateLaughlinBMatrices (LaughlinIndex, SparseBMatrices, NbrBMatrices, Manager.GetInteger("p-truncation"), CylinderFlag, kappa);
      else
	CreateReadRezayiK3StateBMatrices (2, SparseBMatrices, Manager.GetInteger("p-truncation"), CylinderFlag, kappa);
    }
  else
    {
      CreateK2ClusteredStateBMatrices (2, RIndex, SparseBMatrices, Manager.GetInteger("p-truncation"), CylinderFlag, kappa);
    }

  cout << "B matrix size = " << SparseBMatrices[0].GetNbrRow() << "x" << SparseBMatrices[0].GetNbrColumn() << endl;

  RealVector State (Space->GetHilbertSpaceDimension(), true);


  FQHEMPSCreateStateOperation Operation(Space, SparseBMatrices, &State, MatrixElement,
					Manager.GetInteger("precalculation-blocksize"));
  Operation.ApplyOperation(Architecture.GetArchitecture());

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
	  NewState = ((FermionOnSpherePTruncated*) Space)->ConvertToHaldaneBasis(State, *SpaceHaldane);
	  
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

// creates the B matrices for the laughlin state
//
// laughlinIndex = power of the Laughlin part 
// bMatrices = array where the B matrices will be stored
// nbrBMatrices = number of B matrices to compute (max occupation + 1)
// pLevel = |P| level truncation
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

void CreateLaughlinBMatrices (int laughlinIndex, SparseRealMatrix* bMatrices, int nbrBMatrices, int pLevel, bool cylinderFlag, double kappa)
{
  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [pLevel + 1];
  RealMatrix* BMatrices = new RealMatrix[nbrBMatrices];
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

  BMatrices[0] = RealMatrix(MatrixSize, MatrixSize, true);
  for (int i = 0; i <= pLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace = U1BosonBasis[i];
      for (int j = 1; j < NbrNValue; ++j)
	{
	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
	    {
              int N1 = (j - NValueShift/2);
	      double Tmp = 1.0;
              if (cylinderFlag)
                Tmp *= exp(-kappa*kappa*(i + (N1 - 1) * (N1 - 1)/(4.0 * laughlinIndex)+ (N1 * N1)/(4.0 * laughlinIndex)));
	      BMatrices[0].SetMatrixElement(StartingIndexPerPLevel[i] + ((k * NbrNValue) + j - 1), StartingIndexPerPLevel[i] + ((k * NbrNValue) + j), Tmp);
	    }
	}
    }
   
  for (int m = 1; m < nbrBMatrices; ++m)
    BMatrices[m] = RealMatrix(MatrixSize, MatrixSize, true);

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
	  int N2 = (2 * (j - i) - laughlinIndex + 1 + NValueShift) / 2;
	  int N1 = N2 + (laughlinIndex - 1);
	  for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
	    {
	      TmpSpace1->GetOccupationNumber(k1, Partition1);
	      for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
		{
		  TmpSpace2->GetOccupationNumber(k2, Partition2);
		  for (int m = 1; m < nbrBMatrices; ++m)
		    {
		      double Tmp = CreateLaughlinAMatrixElement(laughlinIndex * m * m, 1, Partition1, Partition2, i, j, Coef);
		      if (cylinderFlag)
			Tmp *= exp(-kappa*kappa*(0.5 * i + 0.5 * j + pow(N1 - NValueShift/2,2.0)/(4.0 * laughlinIndex) + pow(N2 - NValueShift/2,2.0)/(4.0 * laughlinIndex)));
		      BMatrices[m].SetMatrixElement(StartingIndexPerPLevel[i] + ((k1 * NbrNValue) + N1), StartingIndexPerPLevel[j] + ((k2 * NbrNValue) + N2), Tmp);
		    }
		}
	    }
	}
    }

  delete[] Partition1;
  delete[] Partition2;

  for (int i = 0; i < nbrBMatrices; ++i)
    {
      cout << i << " " << nbrBMatrices << endl;
      bMatrices[i] = BMatrices[i];
    }
  delete[] BMatrices;
}

double CreateLaughlinAMatrixElement (int chargeNumerator, int chargeDenominator, unsigned long* partition1, unsigned long* partition2, 
				     int p1Level, int p2Level, FactorialCoefficient& coef)
{
  double Tmp = 1.0;
  int PMax = p1Level;
  if (p2Level > p1Level)
    PMax = p2Level;
  for (int i = 1; i <= PMax; ++i)
    {
      double Tmp2 = 0.0;
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
	      coef.PowerNMultiply(chargeNumerator, Sum);
	      coef.PowerNDivide(chargeDenominator, Sum);
	      coef.PowerNDivide(i, Sum);
	      Tmp2 += ((double) (1 - ((j & 1) << 1))) * sqrt(coef.GetNumericalValue());
	    }
	}
      Tmp *= Tmp2;
    }
  return Tmp;
}

// creates the B matrices for the (k,r) clustered with k=2
//
// laughlinIndex = power of the Laughlin part 
// rIndex = r value
// bMatrices = array where the B matrices will be stored
// pLevel = |P| level truncation
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio
// return value = matrix element

void CreateK2ClusteredStateBMatrices (int laughlinIndex, int rIndex, SparseRealMatrix* bMatrices, int pLevel, bool cylinderFlag, double kappa)
{
  LongRational CentralCharge12 ((rIndex + 2l) - (2l * (rIndex - 1l) * (rIndex - 1l)), rIndex + 2l);
  cout << "central charge = " << CentralCharge12 << endl;
  CentralCharge12 /= 12l;
  LongRational WeightIdentity (0l, 1l);
  LongRational WeightPsi (rIndex, 4l);
  long* Partition = new long[2 * (pLevel + 1)];
  unsigned long* TmpPartition = new unsigned long [pLevel + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [pLevel + 1];
  RealSymmetricMatrix* ScalarProductIdentity = new RealSymmetricMatrix[pLevel + 1];
  RealSymmetricMatrix* ScalarProductPsi = new RealSymmetricMatrix[pLevel + 1];
  RealMatrix** MatrixPsi01 = new RealMatrix*[pLevel + 1];
  RealMatrix** MatrixPsi10 = new RealMatrix*[pLevel + 1];
  RealMatrix* OrthogonalBasisIdentityLeft = new RealMatrix[pLevel + 1];
  RealMatrix* OrthogonalBasisPsiLeft = new RealMatrix[pLevel + 1];
  RealMatrix* OrthogonalBasisIdentityRight = new RealMatrix[pLevel + 1];
  RealMatrix* OrthogonalBasisPsiRight = new RealMatrix[pLevel + 1];

  for (int i = 0; i <= pLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort (i, i, pLevel + 1);
      MatrixPsi01[i] = new RealMatrix[pLevel + 1];
      MatrixPsi10[i] = new RealMatrix[pLevel + 1];
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
	    LongRational Tmp = ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightIdentity);
	    ScalarProductIdentity[i].SetMatrixElement(m, n, Tmp.GetNumericalValue());
	    Tmp = ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightPsi);
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
	  OrthogonalBasisIdentityLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  OrthogonalBasisIdentityRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisIdentityLeft[i][Count].Copy(TmpBasis[n]);
		OrthogonalBasisIdentityRight[i][Count].Copy(TmpBasis[n]);
		if (TmpDiag(n, n) > 0)
		  {
		    OrthogonalBasisIdentityLeft[i][Count] /=  sqrt(TmpDiag(n, n));
		    OrthogonalBasisIdentityRight[i][Count] /=  sqrt(TmpDiag(n, n));
		  }
		else
		  {
		    OrthogonalBasisIdentityLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
		    OrthogonalBasisIdentityRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		  }
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisIdentityLeft[i] = RealMatrix();
	  OrthogonalBasisIdentityRight[i] = RealMatrix();
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
	  OrthogonalBasisPsiLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	  OrthogonalBasisPsiRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisPsiLeft[i][Count].Copy(TmpBasis[n]);
		OrthogonalBasisPsiRight[i][Count].Copy(TmpBasis[n]);
		if (TmpDiag(n, n) > 0)
		  {
		    OrthogonalBasisPsiLeft[i][Count] /=  sqrt(TmpDiag(n, n));
		    OrthogonalBasisPsiRight[i][Count] /=  sqrt(TmpDiag(n, n));
		  }
		else
		  {
		    OrthogonalBasisPsiLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
		    OrthogonalBasisPsiRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		  }
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisPsiLeft[i] = RealMatrix();
	  OrthogonalBasisPsiRight[i] = RealMatrix();
	}
      cout << "---------------------------------" << endl;
    }

  cout << "computing Psi matrix elements" << endl;
  LongRational Weight (WeightPsi);
  for (int i = 0; i <= pLevel; ++i)
    {
      for (int j = 0; j <= pLevel; ++j)
	{
	  MatrixPsi01[i][j] = RealMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  MatrixPsi10[i][j] = RealMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  cout << "Levels = " <<  i << " " << j << endl;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    for (int m = 0; m < U1BosonBasis[j]->GetHilbertSpaceDimension(); ++m)
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
		LongRational Tmp = ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightIdentity, WeightPsi, Weight);
		MatrixPsi01[i][j].SetMatrixElement(n, m, Tmp.GetNumericalValue());
		Tmp = ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightPsi, WeightIdentity, Weight);
		MatrixPsi10[i][j].SetMatrixElement(n, m, Tmp.GetNumericalValue());
	      }
	}
    }


  int NbrBMatrices = 2;
  RealMatrix* BMatrices = new RealMatrix[NbrBMatrices];

  int** StartingIndexPerPLevel = new int* [pLevel + 1];
  int* TotalStartingIndexPerPLevel = new int [pLevel + 1];
  int* NbrIndicesPerPLevel = new int [pLevel + 1];
  TotalStartingIndexPerPLevel[0] = 0;
  StartingIndexPerPLevel[0] = new int [1];
  StartingIndexPerPLevel[0][0] = 0;

  
  int NbrNValue;
  int NValueShift;
  int QValue;
  int QValueDenominator;
  if ((rIndex & 1) == 0)
    {
      QValue = 1 + (rIndex / 2);
      NbrNValue = ((2 * pLevel) + QValue) + rIndex / 2 + 1;
      NValueShift = 2 * pLevel - 1;
      QValueDenominator = 1;
    }
  else
    {
      QValue = 2 + rIndex;
      NbrNValue = ((4 * pLevel) + QValue) + rIndex + 1;
      NValueShift = 4 * pLevel - 2;
      QValueDenominator = 2;
    }

     
  NbrIndicesPerPLevel[0] = (U1BosonBasis[0]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[0].GetNbrColumn() + OrthogonalBasisPsiLeft[0].GetNbrColumn())) * NbrNValue;
  for (int i = 1; i <= pLevel; ++i)
    {
      TotalStartingIndexPerPLevel[i] = TotalStartingIndexPerPLevel[i - 1] + NbrIndicesPerPLevel[i - 1];
      StartingIndexPerPLevel[i] = new int [i + 1];      
      StartingIndexPerPLevel[i][0] = TotalStartingIndexPerPLevel[i];
      int Tmp = 0;
      int Tmp2;
      for (int j = 0; j < i; ++j)
	{
	  Tmp2 = U1BosonBasis[j]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[i - j].GetNbrColumn() + OrthogonalBasisPsiLeft[i - j].GetNbrColumn()) * NbrNValue;
	  StartingIndexPerPLevel[i][j + 1] = Tmp2 + StartingIndexPerPLevel[i][j];
	  Tmp += Tmp2;
	}
      Tmp += U1BosonBasis[i]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[0].GetNbrColumn() + OrthogonalBasisPsiLeft[0].GetNbrColumn()) * NbrNValue;
      NbrIndicesPerPLevel[i] =  Tmp;
    }
  
  int MatrixSize = NbrIndicesPerPLevel[pLevel] + TotalStartingIndexPerPLevel[pLevel];
  cout << "B matrix size = " << MatrixSize << endl;
  BMatrices[0] = RealMatrix(MatrixSize, MatrixSize, true);
  for (int i = 0; i <= pLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentityRight = OrthogonalBasisIdentityRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiRight = OrthogonalBasisPsiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductIdentity = ScalarProductIdentity[i - p];
	  RealSymmetricMatrix& TmpScalarProductPsi = ScalarProductPsi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      if ((rIndex & 1) == 0)
		{
		  for (int j = 1; j < NbrNValue; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductIdentity(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisIdentityRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisIdentityLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      BMatrices[0].SetMatrixElement(GetK2MatrixIndex(j - 1, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							    GetK2MatrixIndex(j, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			    }
			}
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductPsi(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsiRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisPsiLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      BMatrices[0].SetMatrixElement(GetK2MatrixIndex(j - 1, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							    GetK2MatrixIndex(j, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			    }
			}
		    }
		}
	      else
		{
		  for (int j = 2; j < NbrNValue; ++j)
		    {
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductIdentity(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisIdentityRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisIdentityLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      BMatrices[0].SetMatrixElement(GetK2MatrixIndex(j - 2, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							    GetK2MatrixIndex(j, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			    }
			}
		      for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
			{
			  for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			    {
			      double Tmp = 0.0;
			      for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
				{
				  double Tmp1 = 0.0;			      
				  for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				    {
				      Tmp1 += TmpScalarProductPsi(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsiRight(NeutralIndex4, NeutralIndex2);				  
				    }
				  Tmp += TmpOrthogonalBasisPsiLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
				}
			      BMatrices[0].SetMatrixElement(GetK2MatrixIndex(j - 2, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							    GetK2MatrixIndex(j, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			    }
			}
		    }
		}	      
	    }
	}
    }
  
  FactorialCoefficient Coef;
  unsigned long* Partition1 = new unsigned long [pLevel + 2];
  unsigned long* Partition2 = new unsigned long [pLevel + 2];
  BMatrices[1] = RealMatrix(MatrixSize, MatrixSize, true);
  for (int i = 0; i <= pLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	  for (int j = 0; j <= pLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		  RealMatrix& TmpMatrixPsi01 = MatrixPsi01[i - p][j - q];
		  RealMatrix& TmpMatrixPsi10 = MatrixPsi10[i - p][j - q];	
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  if ((rIndex & 1) == 0)
			    {
			      N2 = (2 * (j - i) + rIndex + 1 + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (4 * (j - i) + 2 * rIndex + 2 + NValueShift) / 2;
			      N1 = N2 + QValue - 2;
			    }			  
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentity1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsi2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisIdentity1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  BMatrices[1].SetMatrixElement(GetK2MatrixIndex(N1, ChargedIndex1, NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								GetK2MatrixIndex(N2, ChargedIndex2, NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}

			    }
			  if ((rIndex & 1) == 0)
			    {
			      N2 = (2 * (j - i) + 1 + NValueShift) / 2;
			      N1 = N2 + QValue - 1;
			    }
			  else
			    {
			      N2 = (4 * (j - i) + 2 + NValueShift) / 2;
			      N1 = N2 + QValue - 2;
			    }
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentity2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi10(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisIdentity2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisPsi1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  BMatrices[1].SetMatrixElement(GetK2MatrixIndex(N1, ChargedIndex1, NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								GetK2MatrixIndex(N2, ChargedIndex2, NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}
			    }
			}
		    }
		}	      
	    }
	}
    }
  
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      bMatrices[i] = BMatrices[i];
    }
  delete[] BMatrices;
}


inline int GetK2MatrixIndex(int charge, int chargedPartitionIndex, int nbrCharges, int chargeSectorDimension, int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int globalIndexShift)
{
  return (((((nbrIdentityDescendant * fieldIndex) + neutralPartitionIndex) * chargeSectorDimension + chargedPartitionIndex) * nbrCharges + charge) + globalIndexShift);
}

LongRational ComputeVirasoroDescendantScalarProduct (long* partition, int partitionLength, int position, LongRational& centralCharge12, LongRational& weight)
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
	      ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 2, position - 1, centralCharge12, weight));
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
		  * ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, position, centralCharge12, weight));
	}
      else
	{
	  Tmp += ((Store1 - Store2) 
		  * ComputeVirasoroDescendantScalarProduct(partition, partitionLength - 1, position - 1, centralCharge12, weight));
	}
      for (int i = partitionLength - 1; i > position; --i)
	partition[i] = partition[i - 1];
      partition[position] = Store2;
      partition[position - 1] = Store1;
    }

  long Store1 = partition[position - 1];
  partition[position - 1] = partition[position];
  partition[position] = Store1;
  Tmp += ComputeVirasoroDescendantScalarProduct(partition, partitionLength, position + 1, centralCharge12, weight);
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
  if (descendantPosition == partitionLength) 
    {
      LongRational Tmp(1l);
      LongRational TmpSum = weight1;
      TmpSum -= weight2;
      LongRational Tmp2;
      for (int i = partitionLength - 1; i >= 0; --i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  long Tmp3 = 0l;
	  for (int j = 0; j < i; ++j)
	    Tmp3 += partition[j];
	  Tmp2 += TmpSum;
	  Tmp2 += Tmp3;	  
	  Tmp *= Tmp2;
	}
      return Tmp;
    }
  if (position == 0)
    {
      LongRational Tmp(1l);
      LongRational TmpSum = weight1;
      TmpSum -= weight2;
      LongRational Tmp2;
      for (int i = 0; i < partitionLength ; ++i)
	{
	  Tmp2 = weight;
	  Tmp2 *= partition[i];
	  long Tmp3 = 0l;
	  for (int j = i + 1; j < partitionLength; ++j)
	    Tmp3 -= partition[j];
	  Tmp2 += TmpSum;
	  Tmp2 -= Tmp3;
	  Tmp *= Tmp2;
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
    TmpLength += partition[i];
  for (int i = 0; i < (position - 1); ++i)
    TmpLength += partition[i];
  Tmp2 += TmpLength;
  long Store = partition[position - 1];
  for (int i = position; i < partitionLength; ++i)
    partition[i - 1] = partition[i];
  Tmp2 *= ComputeDescendantMatrixElement(partition, partitionLength - 1, descendantPosition - 1, position - 1, 
					 centralCharge12, weight1, weight2, weight);
  for (int i = partitionLength - 1; i >= position; --i)
    partition[i] = partition[i - 1];
  partition[position - 1] = Store;
  Tmp1 += Tmp2;
  return Tmp1;
}

// create the B matrices for the k=3 Read-Rezayi state
//
// laughlinIndex = power of the Laughlin part 
// bMatrices = array where the B matrices will be stored
// pLevel = |P| level truncation
// cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
// kappa = cylinder aspect ratio

void CreateReadRezayiK3StateBMatrices (int laughlinIndex, SparseRealMatrix* bMatrices,int pLevel, bool cylinderFlag, double kappa)
{
  LongRational CentralCharge (4l, 5l);
  cout << "central charge = " << CentralCharge << endl;
  LongRational CentralCharge12(CentralCharge);
  CentralCharge12 /= 12l;
  LongRational WeightIdentity (0l, 1l);
  LongRational WeightPsi (2l, 3l);
  LongRational WeightW (3l);
  long* Partition = new long[2 * (pLevel + 1)];
  unsigned long* TmpPartition = new unsigned long [pLevel + 2];

  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [pLevel + 1];
  RealSymmetricMatrix* ScalarProductIdentity = new RealSymmetricMatrix[pLevel + 1];
  RealSymmetricMatrix* ScalarProductW = new RealSymmetricMatrix[pLevel + 1];
  RealSymmetricMatrix* ScalarProductPsi = new RealSymmetricMatrix[pLevel + 1];
  RealMatrix** MatrixPsi01 = new RealMatrix*[pLevel + 1];
  RealMatrix** MatrixPsi10 = new RealMatrix*[pLevel + 1];
  RealMatrix** MatrixPsi11 = new RealMatrix*[pLevel + 1];
  RealMatrix** MatrixPsi12 = new RealMatrix*[pLevel + 1];
  RealMatrix** MatrixPsi21 = new RealMatrix*[pLevel + 1];
  RealMatrix* OrthogonalBasisIdentityLeft = new RealMatrix[pLevel + 1];
  RealMatrix* OrthogonalBasisWLeft = new RealMatrix[pLevel + 1];
  RealMatrix* OrthogonalBasisPsiLeft = new RealMatrix[pLevel + 1];
  RealMatrix* OrthogonalBasisIdentityRight = new RealMatrix[pLevel + 1];
  RealMatrix* OrthogonalBasisWRight = new RealMatrix[pLevel + 1];
  RealMatrix* OrthogonalBasisPsiRight = new RealMatrix[pLevel + 1];

  for (int i = 0; i <= pLevel; ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort (i, i, pLevel + 1);
      MatrixPsi01[i] = new RealMatrix[pLevel + 1];
      MatrixPsi10[i] = new RealMatrix[pLevel + 1];
      MatrixPsi11[i] = new RealMatrix[pLevel + 1];
      MatrixPsi21[i] = new RealMatrix[pLevel + 1];
      MatrixPsi12[i] = new RealMatrix[pLevel + 1];
    }
  
  for (int i = 0; i <= pLevel; ++i)
    {
      cout << "Level = " <<  i << endl;
      ScalarProductIdentity[i] = RealSymmetricMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      if ((3 + i) <= pLevel)
	ScalarProductW[3 + i] = RealSymmetricMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(), true);
      if (i < 3)
	ScalarProductW[i] = RealSymmetricMatrix();
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
	    LongRational Tmp = ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightIdentity);
	    ScalarProductIdentity[i].SetMatrixElement(m, n, Tmp.GetNumericalValue());
	    if ((3 + i) <= pLevel)
	      {
		Tmp = ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightW);
		ScalarProductW[3 + i].SetMatrixElement(m, n, Tmp.GetNumericalValue());
	      }
	    Tmp = ComputeVirasoroDescendantScalarProduct (Partition, PartitionLength, Position, CentralCharge12, WeightPsi);
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
	  OrthogonalBasisIdentityLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  OrthogonalBasisIdentityRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count, true);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisIdentityLeft[i][Count].Copy(TmpBasis[n]);
		OrthogonalBasisIdentityRight[i][Count].Copy(TmpBasis[n]);
		if (TmpDiag(n, n) > 0)
		  {
		    OrthogonalBasisIdentityLeft[i][Count] /=  sqrt(TmpDiag(n, n));
		    OrthogonalBasisIdentityRight[i][Count] /=  sqrt(TmpDiag(n, n));
		  }
		else
		  {
		    OrthogonalBasisIdentityLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
		    OrthogonalBasisIdentityRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		  }
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisIdentityLeft[i] = RealMatrix();
	  OrthogonalBasisIdentityRight[i] = RealMatrix();
	}
      if ((3 + i) <= pLevel)
	{	  
	  TmpMatrix.Copy(ScalarProductW[3 + i]);
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
	  cout << "nbr of null vectors W_-3 sector = " << Count << endl;
	  if (Count < U1BosonBasis[i]->GetHilbertSpaceDimension())
	    {
	      OrthogonalBasisWLeft[3 + i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	      OrthogonalBasisWRight[3 + i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	      Count = 0;
	      for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
		if (fabs(TmpDiag(n, n)) > Error)
		  {
		    OrthogonalBasisWLeft[3 + i][Count].Copy(TmpBasis[n]);
		    OrthogonalBasisWRight[3 + i][Count].Copy(TmpBasis[n]);
		    if (TmpDiag(n, n) > 0)
		      {
			OrthogonalBasisWLeft[3 + i][Count] /=  sqrt(TmpDiag(n, n));
			OrthogonalBasisWRight[3 + i][Count] /=  sqrt(TmpDiag(n, n));
		      }
		    else
		      {
			OrthogonalBasisWLeft[3 + i][Count] /=  sqrt(-TmpDiag(n, n));
			OrthogonalBasisWRight[3 + i][Count] /=  -sqrt(-TmpDiag(n, n));
		      }
		    ++Count;
	      }
	    }
	  else
	    {
	      OrthogonalBasisWLeft[3 + i] = RealMatrix();
	      OrthogonalBasisWRight[3 + i] = RealMatrix();
	    }
	}
      if (i < 3)
	{
	  OrthogonalBasisWLeft[i] = RealMatrix();
	  OrthogonalBasisWRight[i] = RealMatrix();
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
	  OrthogonalBasisPsiLeft[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	  OrthogonalBasisPsiRight[i] = RealMatrix (U1BosonBasis[i]->GetHilbertSpaceDimension(), U1BosonBasis[i]->GetHilbertSpaceDimension() - Count);
	  Count = 0;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    if (fabs(TmpDiag(n, n)) > Error)
	      {
		OrthogonalBasisPsiLeft[i][Count].Copy(TmpBasis[n]);
		OrthogonalBasisPsiRight[i][Count].Copy(TmpBasis[n]);
		if (TmpDiag(n, n) > 0)
		  {
		    OrthogonalBasisPsiLeft[i][Count] /=  sqrt(TmpDiag(n, n));
		    OrthogonalBasisPsiRight[i][Count] /=  sqrt(TmpDiag(n, n));
		  }
		else
		  {
		    OrthogonalBasisPsiLeft[i][Count] /=  sqrt(-TmpDiag(n, n));
		    OrthogonalBasisPsiRight[i][Count] /=  -sqrt(-TmpDiag(n, n));
		  }
		++Count;
	      }
	}
      else
	{
	  OrthogonalBasisPsiLeft[i] = RealMatrix();
	  OrthogonalBasisPsiRight[i] = RealMatrix();
	}
      cout << "---------------------------------" << endl;
    }

  cout << "computing Psi matrix elements" << endl;
  LongRational Weight (WeightPsi);
  for (int i = 0; i <= pLevel; ++i)
    {
      for (int j = 0; j <= pLevel; ++j)
	{
	  MatrixPsi01[i][j] = RealMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  MatrixPsi10[i][j] = RealMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  MatrixPsi11[i][j] = RealMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	  if ((3 + j) <= pLevel)
	    {	  
	      MatrixPsi12[i][3 + j] = RealMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	    }
	  if (j < 3)
	    {
	      MatrixPsi12[i][j] = RealMatrix();
	    }
	  if ((3 + i) <= pLevel)
	    {	  
	      MatrixPsi21[i + 3][j] = RealMatrix(U1BosonBasis[i]->GetHilbertSpaceDimension(),  U1BosonBasis[j]->GetHilbertSpaceDimension(), true);
	    }
	  if (i < 3)
	    {
	      MatrixPsi21[i][j] = RealMatrix();
	    }
	  cout << "Levels = " <<  i << " " << j << endl;
	  for (int n = 0; n < U1BosonBasis[i]->GetHilbertSpaceDimension(); ++n)
	    for (int m = 0; m < U1BosonBasis[j]->GetHilbertSpaceDimension(); ++m)
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
		LongRational Tmp = ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightIdentity, WeightPsi, Weight);
		MatrixPsi01[i][j].SetMatrixElement(n, m, Tmp.GetNumericalValue());
		Tmp = ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightPsi, WeightIdentity, Weight);
		MatrixPsi10[i][j].SetMatrixElement(n, m, Tmp.GetNumericalValue());
		Tmp = ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightPsi, WeightPsi, Weight);
		MatrixPsi11[i][j].SetMatrixElement(n, m, Tmp.GetNumericalValue());
		if ((3 + j) <= pLevel)
		  {	  
		    Tmp = ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightPsi, WeightW, Weight);
		    MatrixPsi12[i][3 + j].SetMatrixElement(n, m, Tmp.GetNumericalValue() * -sqrt(26.0) / 9.0);
		  }
		if ((3 + i) <= pLevel)
		  {	  
		    Tmp = ComputeDescendantMatrixElement (Partition, PartitionLength, Position, Position, CentralCharge12, WeightW, WeightPsi, Weight);
		    MatrixPsi21[3 + i][j].SetMatrixElement(n, m, Tmp.GetNumericalValue() * sqrt(26.0) / 9.0);
		  }
	      }
	}
    }

  int NbrBMatrices = 2;
  SparseRealMatrix* BMatrices = new SparseRealMatrix[NbrBMatrices];

  int** StartingIndexPerPLevel = new int* [pLevel + 1];
  int* TotalStartingIndexPerPLevel = new int [pLevel + 1];
  int* NbrIndicesPerPLevel = new int [pLevel + 1];
  TotalStartingIndexPerPLevel[0] = 0;
  StartingIndexPerPLevel[0] = new int [1];
  StartingIndexPerPLevel[0][0] = 0;

  
  int NbrNValue;
  int NValueShift;
  int QValue;
  int QValueDenominator;
  QValue = 5;
  QValueDenominator = 3;
  NbrNValue = QValueDenominator * (2 * pLevel) + 4 + 2 + 1;
  NValueShift = QValueDenominator * pLevel;

     
  NbrIndicesPerPLevel[0] = (U1BosonBasis[0]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[0].GetNbrColumn() 
									   + OrthogonalBasisPsiLeft[0].GetNbrColumn() 
									   + OrthogonalBasisPsiLeft[0].GetNbrColumn() 
									   + OrthogonalBasisWLeft[0].GetNbrColumn())) * NbrNValue;
  for (int i = 1; i <= pLevel; ++i)
    {
      TotalStartingIndexPerPLevel[i] = TotalStartingIndexPerPLevel[i - 1] + NbrIndicesPerPLevel[i - 1];
      StartingIndexPerPLevel[i] = new int [i + 1];      
      StartingIndexPerPLevel[i][0] = TotalStartingIndexPerPLevel[i];
      int Tmp = 0;
      int Tmp2;
      for (int j = 0; j < i; ++j)
	{
	  Tmp2 = U1BosonBasis[j]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[i - j].GetNbrColumn() 
								+ OrthogonalBasisPsiLeft[i - j].GetNbrColumn() 
								+ OrthogonalBasisPsiLeft[i - j].GetNbrColumn() 
								+ OrthogonalBasisWLeft[i - j].GetNbrColumn()) * NbrNValue;
	  StartingIndexPerPLevel[i][j + 1] = Tmp2 + StartingIndexPerPLevel[i][j];
	  Tmp += Tmp2;
	}
      Tmp += U1BosonBasis[i]->GetHilbertSpaceDimension() * (OrthogonalBasisIdentityLeft[0].GetNbrColumn() 
							    + OrthogonalBasisPsiLeft[0].GetNbrColumn()
							    + OrthogonalBasisPsiLeft[0].GetNbrColumn() 
							    + OrthogonalBasisWLeft[0].GetNbrColumn()) * NbrNValue;
      NbrIndicesPerPLevel[i] =  Tmp;
    }
  
  int MatrixSize = NbrIndicesPerPLevel[pLevel] + TotalStartingIndexPerPLevel[pLevel];
  cout << "B matrix size = " << MatrixSize << endl;
  BMatrices[0] = SparseRealMatrix(MatrixSize, MatrixSize);
//  BMatrices[0] = RealMatrix(MatrixSize, MatrixSize, true);
  for (int i = 0; i <= pLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentityRight = OrthogonalBasisIdentityRight[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiRight = OrthogonalBasisPsiRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductIdentity = ScalarProductIdentity[i - p];
	  RealSymmetricMatrix& TmpScalarProductPsi = ScalarProductPsi[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      for (int j = QValueDenominator; j < NbrNValue; ++j)
		{
		  // identity 
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentityLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  double Tmp = 0.0;
			  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
			    {
			      double Tmp1 = 0.0;			      
			      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				{
				  Tmp1 += TmpScalarProductIdentity(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisIdentityRight(NeutralIndex4, NeutralIndex2);				  
				}
				  Tmp += TmpOrthogonalBasisIdentityLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
			    }
			  BMatrices[0].SetMatrixElement(GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							GetReadRezayiK3MatrixIndex(j, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			}
		    }
		  // Psi_{+1} and Psi_{-1}
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsiLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  double Tmp = 0.0;
			  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
			    {
			      double Tmp1 = 0.0;			      
			      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				{
				  Tmp1 += TmpScalarProductPsi(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsiRight(NeutralIndex4, NeutralIndex2);				  
				}
			      Tmp += TmpOrthogonalBasisPsiLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
			    }
			  BMatrices[0].SetMatrixElement(GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							GetReadRezayiK3MatrixIndex(j, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			  BMatrices[0].SetMatrixElement(GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 3, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							GetReadRezayiK3MatrixIndex(j, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 3, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			}
		    }
		}
	    }
	}
    }

  for (int i = 0; i <= pLevel; ++i)
    {
      for (int p = 0; p <= (i - 3); ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral = U1BosonBasis[i - p - 3];
	  RealMatrix& TmpOrthogonalBasisIdentityLeft = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsiLeft = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisWLeft = OrthogonalBasisWLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisWRight = OrthogonalBasisWRight[i - p];
	  RealSymmetricMatrix& TmpScalarProductW = ScalarProductW[i - p];
	  for (int ChargedIndex = 0; ChargedIndex < TmpSpaceCharged->GetHilbertSpaceDimension(); ++ChargedIndex)
	    {	      
	      for (int j = QValueDenominator; j < NbrNValue; ++j)
		{
		  // W
		  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisWLeft.GetNbrColumn(); ++NeutralIndex1)
		    {
		      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisWLeft.GetNbrColumn(); ++NeutralIndex2)
			{
			  double Tmp = 0.0;
			  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex3)
			    {
			      double Tmp1 = 0.0;			      
			      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral->GetHilbertSpaceDimension(); ++NeutralIndex4)
				{
				  Tmp1 += TmpScalarProductW(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisWRight(NeutralIndex4, NeutralIndex2);
				}
			      Tmp += TmpOrthogonalBasisWLeft(NeutralIndex3, NeutralIndex1) * Tmp1;
			    }
			  BMatrices[0].SetMatrixElement(GetReadRezayiK3MatrixIndex(j - QValueDenominator, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 5, NeutralIndex1, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
							GetReadRezayiK3MatrixIndex(j, ChargedIndex, NbrNValue, TmpSpaceCharged->GetHilbertSpaceDimension(), 5, NeutralIndex2, TmpOrthogonalBasisIdentityLeft.GetNbrColumn(), TmpOrthogonalBasisPsiLeft.GetNbrColumn(), StartingIndexPerPLevel[i][p]), Tmp);
			}
		    }
		}
	    }
	}
    }

  FactorialCoefficient Coef;
  unsigned long* Partition1 = new unsigned long [pLevel + 2];
  unsigned long* Partition2 = new unsigned long [pLevel + 2];
  BMatrices[1] = SparseRealMatrix(MatrixSize, MatrixSize);
//  BMatrices[1] = RealMatrix(MatrixSize, MatrixSize, true);
  for (int i = 0; i <= pLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	  for (int j = 0; j <= pLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		  RealMatrix& TmpMatrixPsi01 = MatrixPsi01[i - p][j - q];
		  RealMatrix& TmpMatrixPsi10 = MatrixPsi10[i - p][j - q];	
		  RealMatrix& TmpMatrixPsi11 = MatrixPsi11[i - p][j - q];
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  N2 = QValueDenominator * (j - i) + 4 + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisIdentity1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi01(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsi2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisIdentity1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  BMatrices[1].SetMatrixElement(GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 0, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), TmpOrthogonalBasisPsi1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 3, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), TmpOrthogonalBasisPsi2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}

			    }
			  N2 = QValueDenominator * (j - i) + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisIdentity2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi10(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisIdentity2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisPsi1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  BMatrices[1].SetMatrixElement(GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), TmpOrthogonalBasisPsi1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 0, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), TmpOrthogonalBasisPsi2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}
			    }

			  N2 = QValueDenominator * (j - i) + 2 + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi11(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsi2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisPsi1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  BMatrices[1].SetMatrixElement(GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 3, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), TmpOrthogonalBasisPsi1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 1, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), TmpOrthogonalBasisPsi2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}
			    }
			}
		    }
		}	      
	    }
	}
    }
  
  for (int i = 0; i <= pLevel; ++i)
    {
      for (int p = 0; p <= (i - 3); ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p - 3];
	  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisW1 = OrthogonalBasisWLeft[i - p];
	  for (int j = 0; j <= pLevel; ++j)
	    {
	      for (int q = 0; q <= j; ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q];
		  RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		  RealMatrix& TmpOrthogonalBasisW2 = OrthogonalBasisWRight[j - q];
		  RealMatrix& TmpMatrixPsi21 = MatrixPsi21[i - p][j - q];	
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  N2 = QValueDenominator * (j - i) + 4 + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisW1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisPsi2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi21(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisPsi2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisW1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  BMatrices[1].SetMatrixElement(GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 5, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), TmpOrthogonalBasisPsi1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 3, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), TmpOrthogonalBasisPsi2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}

			    }
			}
		    }
		}	      
	    }
	}
    }
  
  for (int i = 0; i <= pLevel; ++i)
    {
      for (int p = 0; p <= i; ++p)
	{
	  BosonOnDiskShort* TmpSpaceCharged1 = U1BosonBasis[p];
	  BosonOnDiskShort* TmpSpaceNeutral1 = U1BosonBasis[i - p];
	  RealMatrix& TmpOrthogonalBasisIdentity1 = OrthogonalBasisIdentityLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisPsi1 = OrthogonalBasisPsiLeft[i - p];
	  RealMatrix& TmpOrthogonalBasisW1 = OrthogonalBasisWLeft[i - p];
	  for (int j = 0; j <= pLevel; ++j)
	    {
	      for (int q = 0; q <= (j - 3); ++q)
		{
		  BosonOnDiskShort* TmpSpaceCharged2 = U1BosonBasis[q];
		  BosonOnDiskShort* TmpSpaceNeutral2 = U1BosonBasis[j - q - 3];
		  RealMatrix& TmpOrthogonalBasisIdentity2 = OrthogonalBasisIdentityRight[j - q];
		  RealMatrix& TmpOrthogonalBasisPsi2 = OrthogonalBasisPsiRight[j - q];
		  RealMatrix& TmpOrthogonalBasisW2 = OrthogonalBasisWRight[j - q];
		  RealMatrix& TmpMatrixPsi12 = MatrixPsi12[i - p][j - q];
	  
		  for (int ChargedIndex1 = 0; ChargedIndex1 < TmpSpaceCharged1->GetHilbertSpaceDimension(); ++ChargedIndex1)
		    {	      
		      TmpSpaceCharged1->GetOccupationNumber(ChargedIndex1, Partition1);
		      for (int ChargedIndex2 = 0; ChargedIndex2 < TmpSpaceCharged2->GetHilbertSpaceDimension(); ++ChargedIndex2)
			{	      
			  TmpSpaceCharged2->GetOccupationNumber(ChargedIndex2, Partition2);
			  int N2;
			  int N1;
			  N2 = QValueDenominator * (j - i) + NValueShift;
			  N1 = N2 + QValue - QValueDenominator;
			  for (int NeutralIndex1 = 0; NeutralIndex1 < TmpOrthogonalBasisPsi1.GetNbrColumn(); ++NeutralIndex1)
			    {
			      for (int NeutralIndex2 = 0; NeutralIndex2 < TmpOrthogonalBasisW2.GetNbrColumn(); ++NeutralIndex2)
				{
				  double Tmp = 0.0;
				  for (int NeutralIndex3 = 0; NeutralIndex3 < TmpSpaceNeutral1->GetHilbertSpaceDimension(); ++NeutralIndex3)
				    {
				      double Tmp1 = 0.0;			      
				      for (int NeutralIndex4 = 0; NeutralIndex4 < TmpSpaceNeutral2->GetHilbertSpaceDimension(); ++NeutralIndex4)
					{
					  Tmp1 += TmpMatrixPsi12(NeutralIndex3, NeutralIndex4) * TmpOrthogonalBasisW2(NeutralIndex4, NeutralIndex2);				  
					}
				      Tmp += TmpOrthogonalBasisPsi1(NeutralIndex3, NeutralIndex1) * Tmp1;
				    }
				  Tmp *= CreateLaughlinAMatrixElement(QValue, QValueDenominator, Partition1, Partition2, i, j, Coef);
				  BMatrices[1].SetMatrixElement(GetReadRezayiK3MatrixIndex(N1, ChargedIndex1, NbrNValue, TmpSpaceCharged1->GetHilbertSpaceDimension(), 1, NeutralIndex1, TmpOrthogonalBasisIdentity1.GetNbrColumn(), TmpOrthogonalBasisPsi1.GetNbrColumn(), StartingIndexPerPLevel[i][p]), 
								GetReadRezayiK3MatrixIndex(N2, ChargedIndex2, NbrNValue, TmpSpaceCharged2->GetHilbertSpaceDimension(), 5, NeutralIndex2, TmpOrthogonalBasisIdentity2.GetNbrColumn(), TmpOrthogonalBasisPsi2.GetNbrColumn(), StartingIndexPerPLevel[j][q]), Tmp);
				}
			    }
			}
		    }
		}	      
	    }
	}
    }

   for (int i = 0; i < NbrBMatrices; ++i)
     {
       bMatrices[i] = BMatrices[i];
//       cout << "B" << i << " = " << endl;
//       bMatrices[i].PrintNonZero(cout) << endl;
     }
//   delete[] BMatrices;
}

// compute the linearized index of the B matrix for the k=3 Read-Rezayi state
//
// fieldIndex = field index (0 for the identity, 1 for psi_{+1}, 3 for psi_{-1}, 5 for W_-3)
// return value = linearized index

int GetReadRezayiK3MatrixIndex(int charge, int chargedPartitionIndex, int nbrCharges, int chargeSectorDimension, int fieldIndex, int neutralPartitionIndex, int nbrIdentityDescendant, int nbrPsiDescendant, int globalIndexShift)
{
  return ((((((nbrIdentityDescendant * (fieldIndex & 1)) + (nbrPsiDescendant * (fieldIndex >> 1))) + neutralPartitionIndex) * chargeSectorDimension + chargedPartitionIndex) * nbrCharges + charge) + globalIndexShift);
}
