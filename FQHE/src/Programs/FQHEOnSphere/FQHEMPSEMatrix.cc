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

#include "MainTask/FQHEMPSEMatrixMainTask.h"

#include "GeneralTools/ArrayTools.h"

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
  
  OptionManager Manager ("FQHEMPSEMatrix" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  FQHEMPSMatrixManager MPSMatrixManager;

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  OptionGroup* ArnoldiGroup  = new OptionGroup ("Arnoldi options");
  Architecture.AddOptionGroup(&Manager);
  Manager += ArnoldiGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new BooleanOption  ('\n', "diagonal-block", "consider only the block diagonal in P and Q");
  (*SystemGroup) += new BooleanOption  ('\n', "fixed-qsector", "consider only the block diagonal in P and Q");
  
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "memory", "amount of memory that can used for precalculations (in Mb)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "ematrix-memory", "amount of memory that can used for precalculations of the E matrix (in Mb)", 500);
  (*OutputGroup) += new SingleStringOption  ('o', "output-file", "output file name");
  (*OutputGroup) += new BooleanOption  ('\n', "show-bmatrices", "show the B matrices");
  (*OutputGroup) += new BooleanOption  ('\n', "show-ematrix", "show the transfer matrix");
  (*ArnoldiGroup) += new SingleIntegerOption  ('\n', "full-diag", 
					       "maximum Hilbert space dimension for which full diagonalization is applied", 1000);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "disk", "enable disk storage for the Arnoldi algorithm", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "resume", "resume from disk datas", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Arnoldi iteration", false); 
  (*ArnoldiGroup) += new BooleanOption  ('\n', "power-method", "use the power method instead of the Arnoldi algorithm. A single eigenvalue is computed (the one with the largest real part)", false); 
  (*ArnoldiGroup) += new  SingleIntegerOption ('\n', "arnoldi-memory", "amount of memory when using the Arnoldi algorithm (in Mb)", 500); 
  (*ArnoldiGroup) += new BooleanOption  ('\n', "implicitly-restarted", "use the implicitly restarted Arnoldi algorithm", false); 
  (*ArnoldiGroup) += new  SingleIntegerOption ('\n', "nbr-excited", "number of eigenvalues to compute above the groundstate", 0);
  (*ArnoldiGroup) += new  SingleIntegerOption ('\n', "nbr-eigenstates", "number of eigenvalues to compute (if set to zero, this number will be deduced from the state and nbr-excited)", 0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEMPSEMatrix -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0; 
  int NbrFluxQuanta = 1;
  int TotalLz = 0;

  bool CylinderFlag = Manager.GetBoolean("normalize-cylinder");

  int LandauLevel = 0;

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  int NbrBMatrices = 2;
  
  int MPSRowIndex = 0;
  int MPSColumnIndex = 0;

  SparseRealMatrix* SparseBMatrices = MPSMatrix->GetMatrices();
  cout << "B matrix size = " << SparseBMatrices[0].GetNbrRow() << "x" << SparseBMatrices[0].GetNbrColumn() << endl;
  int NbrEigenstates = MPSMatrix->GetTransferMatrixLargestEigenvalueDegeneracy();
  NbrEigenstates += Manager.GetInteger("nbr-excited");
  if (Manager.GetInteger("nbr-eigenstates") > 0)
    {
      NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
    }

  char* StateName = new char [256];
  strcpy(StateName, MPSMatrix->GetName());

  double EnergyShift = 0.0;
  if (Manager.GetBoolean("power-method") == true)
    {
      EnergyShift = 10.0;
      NbrEigenstates = 1;
    }

  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
    }
  else
    {
      OutputFileName  = new char [512];
      if (CylinderFlag == true)
	{
	  if (Manager.GetBoolean("diagonal-block"))
	    {
	      sprintf(OutputFileName, "ematrix_diagblock_cylinder_%s_perimeter_%f_plevel_%ld.dat", StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"));
	    }
	  else
	    {
	      sprintf(OutputFileName, "ematrix_cylinder_%s_perimeter_%f_plevel_%ld.dat", StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"));
	    }
	}
      else
	{
	  if (Manager.GetBoolean("diagonal-block"))
	    {
	      sprintf(OutputFileName, "ematrix_diagblock_%s_plevel_%ld.dat", StateName,
		      Manager.GetInteger("p-truncation"));
	    }
	  else
	    {
	      sprintf(OutputFileName, "ematrix_%s_plevel_%ld.dat", StateName,
		      Manager.GetInteger("p-truncation"));
	    }
	}
   }


  double Error = 1e-13;
  
  SparseRealMatrix* SparseTransposeBMatrices = new SparseRealMatrix[NbrBMatrices];
  double* Coefficients = new double[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      Coefficients[i] = 1.0;
      SparseTransposeBMatrices[i] = SparseBMatrices[i].Transpose();
    }
  
  int MinQValue = 0;
  int MaxQValue = 0;
  MPSMatrix->GetChargeIndexRange(MinQValue, MaxQValue);
  int TmpBMatrixDimension = SparseBMatrices[0].GetNbrRow();
  
  TensorProductSparseMatrixHamiltonian* ETransposeHamiltonian = 0;
  if (Manager.GetBoolean("diagonal-block"))
    {
      if (Manager.GetBoolean("show-bmatrices"))
	{
	  cout << "----------------------------------" << endl;
	  for (int i = 0; i < NbrBMatrices; ++i)
	    {
	      cout << "B^[" << i << "]" << endl;
	      SparseBMatrices[i].PrintNonZero(cout);
	      cout << "----------------------------------" << endl;
	    }
	}
      long EffectiveDimension = 0l;
      for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	{
	  for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
	    {
	      long Tmp = MPSMatrix->GetBondIndexRange(PLevel, QValue);
	      EffectiveDimension += Tmp * Tmp;
	    }
	}
      
      cout << "computing effective E matrix indices " << endl;
      long** BlockIndexProductTable = new long* [TmpBMatrixDimension];
      int* BlockIndexProductTableNbrElements = new int [TmpBMatrixDimension];
      int* BlockIndexProductTableShift = new int [TmpBMatrixDimension];
      for (long i = 0; i < TmpBMatrixDimension; ++i)
	{
	  BlockIndexProductTableNbrElements[i] = 0;
	  BlockIndexProductTableShift[i] = -1;	  
	}
      
      long* EffectiveBlockIndices = new long [EffectiveDimension];
      EffectiveDimension = 0l;
      for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	{
	  long Tmp = MPSMatrix->GetBondIndexRange(PLevel, MaxQValue);
	  for (int i = 0; i < Tmp; ++i)
	    {
	      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
		{
		  long Tmp2 = ((long) MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(i, PLevel, QValue));
		  BlockIndexProductTableNbrElements[Tmp2] = Tmp;
		  BlockIndexProductTableShift[Tmp2] = EffectiveDimension;
		  BlockIndexProductTable[Tmp2] = new long[Tmp];
		  long* TmpBlockIndexProductTable = BlockIndexProductTable[Tmp2];
		  Tmp2 *= TmpBMatrixDimension;
		  for (int j = 0; j < Tmp; ++j)
		    {
		      TmpBlockIndexProductTable[j] = Tmp2 + MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(j, PLevel, QValue);
		      EffectiveBlockIndices[EffectiveDimension] = TmpBlockIndexProductTable[j];		      
		      ++EffectiveDimension;
		    }
		}
	    }
	}
      for (long i = 1l; i < EffectiveDimension; ++i)
	if (EffectiveBlockIndices[i] < EffectiveBlockIndices[i - 1l])
	  cout << "error, unsorted indices" << endl;
      cout << "E matrix effective dimension = " << EffectiveDimension << "( vs " << (SparseBMatrices[0].GetNbrRow() * SparseBMatrices[0].GetNbrRow()) << ")" << endl;
      Architecture.GetArchitecture()->SetDimension(EffectiveDimension);
      long Memory = Manager.GetInteger("ematrix-memory") << 20;
       if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
     ETransposeHamiltonian = new TensorProductSparseMatrixSelectedBlockHamiltonian(NbrBMatrices, SparseBMatrices, SparseBMatrices, Coefficients, 
										   EffectiveDimension, EffectiveBlockIndices, 
										   BlockIndexProductTable, BlockIndexProductTableNbrElements, BlockIndexProductTableShift, 
										   Architecture.GetArchitecture(), Manager.GetInteger("ematrix-memory") << 20);    
    }
  else
    {
      Architecture.GetArchitecture()->SetDimension(((long) TmpBMatrixDimension) * ((long) TmpBMatrixDimension));
      if (Manager.GetBoolean("fixed-qsector"))
	{
	  int GroupSize = 1;
	  int TmpNumerator;
	  MPSMatrix->GetFillingFactor(TmpNumerator, GroupSize);
	  int NbrGroupBMatrices = 1 << GroupSize;
	  SparseRealMatrix* TmpSparseGroupBMatrices = new SparseRealMatrix[NbrGroupBMatrices];
	  double* GroupCoefficients = new double[NbrGroupBMatrices];
	  for (int i = 0; i < NbrGroupBMatrices; ++i)
	    {
	      GroupCoefficients[i] = 1.0;
	      TmpSparseGroupBMatrices[i].Copy(SparseBMatrices[i & 1]);
	      for (int j = 1; j < GroupSize; ++j)
		TmpSparseGroupBMatrices[i].Multiply(SparseBMatrices[(i >> j) & 1]);
	    }
	  if (Manager.GetBoolean("show-bmatrices"))
	    {
	      cout << "----------------------------------" << endl;
	      for (int i = 0; i < NbrGroupBMatrices; ++i)
		{
		  cout << "B^[" << i << "]" << endl;
		  TmpSparseGroupBMatrices[i].PrintNonZero(cout);
		  cout << "----------------------------------" << endl;
		}
	    }	  
	  int MinQ;
	  int MaxQ;
	  MPSMatrix->GetChargeIndexRange(MinQ, MaxQ);
	  MinQ = 0;
	  int GroupBMatrixDimension = 0l;
	  for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	    {
	      for (int QValue = MinQ; QValue <= MaxQ; QValue += GroupSize)
		GroupBMatrixDimension += MPSMatrix->GetBondIndexRange(PLevel, QValue);
	    }
	  cout << "group B matrix size = " <<  GroupBMatrixDimension << "x" << GroupBMatrixDimension << endl;
	  bool* GlobalIndices = new bool [TmpBMatrixDimension];
	  for (int i = 0; i < TmpBMatrixDimension; ++i)
	    GlobalIndices[i] = false;
	  GroupBMatrixDimension = 0;
	  for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	    {
	      for (int QValue = MinQ; QValue <= MaxQ; QValue += GroupSize)
		{
		  int MaxLocalIndex = MPSMatrix->GetBondIndexRange(PLevel, QValue);
		  for (int i = 0; i < MaxLocalIndex; ++i)
		    {
		      GlobalIndices[MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(i, PLevel, QValue)] = true;
		      ++GroupBMatrixDimension;
		    }
		}
	    }
	  SparseRealMatrix* TmpSparseGroupBMatrices2 = new SparseRealMatrix[NbrGroupBMatrices];
	  int NbrGroupBMatrices2 = 0;
	  for (int i = 0; i < NbrGroupBMatrices; ++i)
	    {
	      TmpSparseGroupBMatrices2[i] = TmpSparseGroupBMatrices[i].ExtractMatrix(GroupBMatrixDimension, GroupBMatrixDimension, GlobalIndices, GlobalIndices);
	      if (TmpSparseGroupBMatrices2[i].GetNbrRow() > 0)
		++NbrGroupBMatrices2;
	    }
	  delete[] GlobalIndices;
	  SparseRealMatrix* TmpSparseGroupBMatrices3 = new SparseRealMatrix[NbrGroupBMatrices2];
	  NbrGroupBMatrices2 = 0;
	  for (int i = 0; i < NbrGroupBMatrices; ++i)
	    {
	      if (TmpSparseGroupBMatrices2[i].GetNbrRow() > 0)
		{
		  TmpSparseGroupBMatrices3[NbrGroupBMatrices2] = TmpSparseGroupBMatrices2[i];
		  ++NbrGroupBMatrices2;
		}
	    }
	  if (Manager.GetBoolean("show-bmatrices"))
	    {
	      cout << "----------------------------------" << endl;
	      for (int i = 0; i < NbrGroupBMatrices2; ++i)
		{
		  cout << "B^[" << i << "]" << endl;
		  TmpSparseGroupBMatrices3[i].PrintNonZero(cout);
		  cout << "----------------------------------" << endl;
		}
	    }	  
	  ETransposeHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrGroupBMatrices2, TmpSparseGroupBMatrices3, TmpSparseGroupBMatrices3, GroupCoefficients); 
	}
      else
	{
	  ETransposeHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, SparseBMatrices, SparseBMatrices, Coefficients); 
	  if (Manager.GetBoolean("show-bmatrices"))
	    {
	      cout << "----------------------------------" << endl;
	      for (int i = 0; i < NbrBMatrices; ++i)
		{
		  cout << "B^[" << i << "]" << endl;
		  SparseBMatrices[i].PrintNonZero(cout);
		  cout << "----------------------------------" << endl;
		}
	    }	  
	}
    }
  if (Manager.GetBoolean("power-method") == true)
    ETransposeHamiltonian->ShiftHamiltonian(EnergyShift);
  
  if (Manager.GetBoolean("show-ematrix"))
    {
      ComplexMatrix EMatrix (ETransposeHamiltonian->GetHilbertSpaceDimension(), ETransposeHamiltonian->GetHilbertSpaceDimension());
      ETransposeHamiltonian->GetHamiltonian(EMatrix);
      cout << EMatrix << endl;
    }

  FQHEMPSEMatrixMainTask TaskLeft(&Manager, ETransposeHamiltonian, NbrEigenstates, false, true, 1e-10, EnergyShift, OutputFileName);
  MainTaskOperation TaskOperationLeft (&TaskLeft);
  TaskOperationLeft.ApplyOperation(Architecture.GetArchitecture());
  return 0;
}

