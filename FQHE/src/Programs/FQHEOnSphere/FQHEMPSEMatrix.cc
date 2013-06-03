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
#include "GeneralTools/FilenameTools.h"

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
  (*SystemGroup) += new BooleanOption  ('\n', "right-eigenstates", "compute the right eigenstates");
  (*SystemGroup) += new BooleanOption  ('\n', "left-eigenstates", "compute the left eigenstates");
  
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "memory", "amount of memory that can used for precalculations (in Mb)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "ematrix-memory", "amount of memory that can used for precalculations of the E matrix (in Mb)", 500);
  (*OutputGroup) += new SingleStringOption  ('o', "output-file", "output file name");
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

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  
  int MPSRowIndex = 0;
  int MPSColumnIndex = 0;

  SparseRealMatrix* SparseBMatrices = MPSMatrix->GetMatrices();
  int NbrBMatrices = MPSMatrix->GetNbrMatrices();
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

  char* PrefixOutputFileName = 0;
  char* OutputFileName = 0;
  if (Manager.GetString("output-file") != 0)
    {
      OutputFileName = new char [strlen(Manager.GetString("output-file")) + 1];
      strcpy (OutputFileName, Manager.GetString("output-file"));
      if (GetExtensionFromFileName(OutputFileName) != 0)
	{
	  int TmpLength = strlen(OutputFileName) - strlen(GetExtensionFromFileName(OutputFileName));
	  PrefixOutputFileName = new char[TmpLength + 1];
	  strncpy(PrefixOutputFileName, OutputFileName, TmpLength);
	  PrefixOutputFileName[TmpLength] = '\0';
	}
    }
  else
    {
      PrefixOutputFileName  = new char [1024];
      if (CylinderFlag == true)
	{
	  if (Manager.GetBoolean("diagonal-block"))
	    {
	      sprintf(PrefixOutputFileName, "ematrix_diagblock_cylinder_%s_perimeter_%f_plevel_%ld", StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"));
	    }
	  else
	    {
	      sprintf(PrefixOutputFileName, "ematrix_cylinder_%s_perimeter_%f_plevel_%ld", StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"));
	    }
	}
      else
	{
	  if (Manager.GetBoolean("diagonal-block"))
	    {
	      sprintf(PrefixOutputFileName, "ematrix_diagblock_%s_plevel_%ld", StateName,
		      Manager.GetInteger("p-truncation"));
	    }
	  else
	    {
	      sprintf(PrefixOutputFileName, "ematrix_%s_plevel_%ld", StateName,
		      Manager.GetInteger("p-truncation"));
	    }
	}
      OutputFileName = new char[strlen(PrefixOutputFileName) + 64];
      sprintf(OutputFileName, "%s.dat", PrefixOutputFileName);
   }


  double Error = 1e-13;
  
  SparseRealMatrix* SparseTransposeBMatrices = new SparseRealMatrix[NbrBMatrices];
  double* Coefficients = new double[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    {
      Coefficients[i] = 1.0;
      SparseTransposeBMatrices[i] = SparseBMatrices[i].Transpose();
    }
  
  int TmpBMatrixDimension = SparseBMatrices[0].GetNbrRow();
  
  TensorProductSparseMatrixHamiltonian* ETransposeHamiltonian = 0;
  TensorProductSparseMatrixHamiltonian* EHamiltonian = 0;
  if (Manager.GetBoolean("diagonal-block"))
    {
      long EffectiveDimension = 0l;
      for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	{
	  int MinQValue = 0;
	  int MaxQValue = 0;
	  MPSMatrix->GetChargeIndexRange(PLevel, MinQValue, MaxQValue);
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
	  int MinQValue = 0;
	  int MaxQValue = 0;
	  MPSMatrix->GetChargeIndexRange(PLevel, MinQValue, MaxQValue);
	  long Tmp = MPSMatrix->GetBondIndexRange(PLevel, MaxQValue);
	  for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
	    {
	      for (int i = 0; i < Tmp; ++i)
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
      if ((Manager.GetBoolean("right-eigenstates") == true) || (Manager.GetBoolean("left-eigenstates") == false))
	ETransposeHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, SparseBMatrices, SparseBMatrices, Coefficients); 
      if (Manager.GetBoolean("left-eigenstates") == true)
	{
	  SparseRealMatrix* ConjugateSparseBMatrices = new SparseRealMatrix[NbrBMatrices];
	  for (int i = 0; i < NbrBMatrices; ++i)
	    ConjugateSparseBMatrices[i] = SparseBMatrices[i].Transpose();
	  RealMatrix Test1(SparseBMatrices[1]);
	  RealMatrix Test2(ConjugateSparseBMatrices[1]);

	  EHamiltonian = new TensorProductSparseMatrixHamiltonian(NbrBMatrices, ConjugateSparseBMatrices, ConjugateSparseBMatrices, Coefficients); 
	}
     }
  if (Manager.GetBoolean("power-method") == true)
    {
      if (ETransposeHamiltonian != 0)
	ETransposeHamiltonian->ShiftHamiltonian(EnergyShift);
      if (EHamiltonian != 0)
	EHamiltonian->ShiftHamiltonian(EnergyShift);
    }
  
  if (Manager.GetBoolean("show-ematrix"))
    {
      if (ETransposeHamiltonian != 0)
	{
	  ComplexMatrix EMatrix (ETransposeHamiltonian->GetHilbertSpaceDimension(), ETransposeHamiltonian->GetHilbertSpaceDimension());
	  ETransposeHamiltonian->GetHamiltonian(EMatrix);
	  cout << EMatrix << endl;
	}
      if (EHamiltonian != 0)
	{
	  ComplexMatrix EMatrix (EHamiltonian->GetHilbertSpaceDimension(), EHamiltonian->GetHilbertSpaceDimension());
	  EHamiltonian->GetHamiltonian(EMatrix);
	  cout << EMatrix << endl;
	}      
    }

  if ((Manager.GetBoolean("right-eigenstates") == false) && (Manager.GetBoolean("left-eigenstates") == false))
    {
      FQHEMPSEMatrixMainTask TaskLeft(&Manager, ETransposeHamiltonian, NbrEigenstates, false, true, 1e-10, EnergyShift, OutputFileName);
      MainTaskOperation TaskOperationLeft (&TaskLeft);
      TaskOperationLeft.ApplyOperation(Architecture.GetArchitecture());
    }
  if (Manager.GetBoolean("right-eigenstates") == true)
    {
      cout << "computing right eigenstates" << endl;
      char* EigenvectorFileName = new char [strlen(PrefixOutputFileName) + 128];
      sprintf(EigenvectorFileName, "%s_right", PrefixOutputFileName);
      FQHEMPSEMatrixMainTask TaskLeft(&Manager, ETransposeHamiltonian, NbrEigenstates, true, true, 1e-10, EnergyShift, OutputFileName, 0, EigenvectorFileName);
      MainTaskOperation TaskOperationLeft (&TaskLeft);
      TaskOperationLeft.ApplyOperation(Architecture.GetArchitecture());
    }
  if (Manager.GetBoolean("left-eigenstates") == true)
    {
      cout << "computing left eigenstates" << endl;
      char* EigenvectorFileName = new char [strlen(PrefixOutputFileName) + 128];
      sprintf(EigenvectorFileName, "%s_left", PrefixOutputFileName);
      FQHEMPSEMatrixMainTask TaskLeft(&Manager, EHamiltonian, NbrEigenstates, true, false, 1e-10, EnergyShift, OutputFileName, 0, EigenvectorFileName);
      MainTaskOperation TaskOperationLeft (&TaskLeft);
      TaskOperationLeft.ApplyOperation(Architecture.GetArchitecture());
    }
  return 0;
}

