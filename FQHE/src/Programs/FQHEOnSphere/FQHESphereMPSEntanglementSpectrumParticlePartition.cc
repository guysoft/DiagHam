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
  
  OptionManager Manager ("FQHESphereMPSEntanglementSpectrumParticlePartition" , "0.01");
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

  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "file that describes the root configuration");
  (*SystemGroup) += new BooleanOption  ('\n', "use-padding", "root partitions use the extra zero padding");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "na", "number of particles in subsystem A", 0);
  (*SystemGroup) += new BooleanOption ('\n', "all-na", "print all charge sectors");
  (*SystemGroup) += new BooleanOption ('\n', "infinite-cylinder", "evaluate the entanglement spectrum on the infinite cylinder");
  (*SystemGroup) += new BooleanOption ('\n', "use-singlestate", "use a single real eigenstate of the E matrix  when evaluating the infinite entanglement spectrum");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "memory", "amount of memory that can used for precalculations (in Mb)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "ematrix-memory", "amount of memory that can used for precalculations of the E matrix (in Mb)", 500);
  (*OutputGroup) += new SingleStringOption  ('o', "output-file", "output file name");
  (*OutputGroup) += new BooleanOption ('\n', "suppress-output", "minimize the amount of output information");

  (*OutputGroup) += new BooleanOption ('n', "normalize-sphere", "express the MPS in the normalized sphere basis");
  (*ArnoldiGroup) += new SingleIntegerOption  ('\n', "full-diag", 
					       "maximum Hilbert space dimension for which full diagonalization is applied", 1000);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "disk", "enable disk storage for the Arnoldi algorithm", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "resume", "resume from disk datas", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Arnoldi iteration", false); 
  (*ArnoldiGroup) += new  SingleIntegerOption ('\n', "arnoldi-memory", "amount of memory when using the Arnoldi algorithm (in Mb)", 500); 
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");  

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereMPSEntanglementSpectrumParticlePartition -h" << endl;
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
  if ((Manager.GetBoolean("infinite-cylinder") == false) && (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false))
    return -1;

  int EntCut = Manager.GetInteger("la");
  int Na = Manager.GetInteger("na");

  bool MinimizeOutput = Manager.GetBoolean("suppress-output");

  bool CylinderFlag = Manager.GetBoolean("normalize-cylinder");

  int LandauLevel = 0;

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  int NbrBMatrices = MPSMatrix->GetNbrMatrices();
  SparseRealMatrix* BMatrices = MPSMatrix->GetMatrices();
  SparseRealMatrix* ConjugateBMatrices = new SparseRealMatrix[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    ConjugateBMatrices[i] = BMatrices[i].Transpose();

  cout << "B matrix size = " << BMatrices[0].GetNbrRow() << "x" << BMatrices[0].GetNbrColumn() << endl;
  
  int MPSRowIndex = 0;
  int MPSColumnIndex = 0;
  int NbrEigenstates = 0;
  int MinQ;
  int MaxQ;
  MPSMatrix->GetChargeIndexRange(0, MinQ, MaxQ);
  MPSMatrix->GetMatrixBoundaryIndices(MPSRowIndex, MPSColumnIndex, Manager.GetBoolean("use-padding"));
  NbrEigenstates = MPSMatrix->GetTransferMatrixLargestEigenvalueDegeneracy();

  ofstream File;
  File.precision(14);
  
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = new char [512];
      char* StateName = new char [256];
      strcpy(StateName, MPSMatrix->GetName());
      if (CylinderFlag == true)
	{
	  if (Manager.GetBoolean("infinite-cylinder"))
	    {
	      sprintf(TmpFileName, "fermions_infinite_cylinder_%s_perimeter_%f_plevel_%ld_n_0_2s_0_lz_0.0.full.parent", StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"));
	    }
	  else
	    {
	      sprintf(TmpFileName, "fermions_cylinder_%s_plevel_%ld_n_%d_2s_%d_lz_%d.0.full.parent", StateName,
		      Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz);
	    }
	}
      else
	{
	  if (Manager.GetBoolean("infinite-cylinder"))
	    {
	      sprintf(TmpFileName, "fermions_infinite_%s_plevel_%ld_n_0_2s_0_lz_0.0.full.parent", StateName,
		      Manager.GetInteger("p-truncation"));
	    }
	  else
	    {
	      sprintf(TmpFileName, "fermions_%s_plevel_%ld_n_%d_2s_%d_lz_%d.0.full.parent", StateName,
		      Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz);
	    }
	}
      File.open(TmpFileName, ios::binary | ios::out);     
   }


  if (Manager.GetBoolean("infinite-cylinder"))
    {
      return 0;
    }

  int QValue = 3;
  int PLevel = Manager.GetInteger("p-truncation");


  long MaxTmpMatrixElements = (((long) BMatrices[0].GetNbrRow()) * 
				((long) BMatrices[0].GetNbrRow() / 1l));
  cout << "Requested memory for sparse matrix multiplications = " << ((MaxTmpMatrixElements * (2l * sizeof(double) + sizeof(int))) >> 20) << "Mb" << endl;
  double* TmpMatrixElements = new double [MaxTmpMatrixElements];
  int* TmpColumnIndices = new int [MaxTmpMatrixElements];
  double* TmpElements = new double [BMatrices[0].GetNbrRow()];

  SparseRealMatrix FullOverlapMatrixA (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
  FullOverlapMatrixA.SetMatrixElement(MPSRowIndex, MPSRowIndex, 1.0);  
  int MaxNbrFluxQuantaA = PLevel + QValue * (Na - 1);
  if (MaxNbrFluxQuantaA > NbrFluxQuanta)
    {
      MaxNbrFluxQuantaA = NbrFluxQuanta;
    }
  for (int i = 0; i <= MaxNbrFluxQuantaA; ++i)
    {
      SparseRealMatrix TmpMatrix2;
      SparseRealMatrix TmpMatrix3;
      TmpMatrix2 = Conjugate(ConjugateBMatrices[0], FullOverlapMatrixA, BMatrices[0], 
			     TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix3 = Conjugate(ConjugateBMatrices[1], FullOverlapMatrixA, BMatrices[1], 
			     TmpMatrixElements, TmpColumnIndices, TmpElements); 
      FullOverlapMatrixA = TmpMatrix2 + TmpMatrix3;
    }

  SparseRealMatrix FullOverlapMatrixB (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
  FullOverlapMatrixB.SetMatrixElement(MPSColumnIndex, MPSColumnIndex, 1.0);  
  int MaxNbrFluxQuantaB = PLevel + QValue * (Na - 1);
  if (MaxNbrFluxQuantaB > NbrFluxQuanta)
    {
      MaxNbrFluxQuantaB = NbrFluxQuanta;
    }
  for (int i = 0; i <= MaxNbrFluxQuantaB; ++i)
    {
      SparseRealMatrix TmpMatrix2;
      SparseRealMatrix TmpMatrix3;
      TmpMatrix2 = Conjugate(ConjugateBMatrices[0], FullOverlapMatrixB, BMatrices[0], 
			     TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix3 = Conjugate(ConjugateBMatrices[1], FullOverlapMatrixB, BMatrices[1], 
			     TmpMatrixElements, TmpColumnIndices, TmpElements); 
      FullOverlapMatrixB = TmpMatrix2 + TmpMatrix3;
    }


  for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
    {
      int NbrFluxQuantaA = CurrentPLevel + QValue * (Na - 1);
      int NbrFluxQuantaB = NbrFluxQuanta + CurrentPLevel - QValue * Na;
      if ((NbrFluxQuantaA >= 0) && (NbrFluxQuantaA <= NbrFluxQuanta) && 
	  (NbrFluxQuantaB >= 0) && (NbrFluxQuantaB <= NbrFluxQuanta))
	{
	  int QValueA = 0;
	  int QValueB = 0;
	  SparseRealMatrix OverlapMatrixA = MPSMatrix->ExtractBlock (FullOverlapMatrixA, CurrentPLevel, QValueA, CurrentPLevel, QValueA);
	  SparseRealMatrix OverlapMatrixB = MPSMatrix->ExtractBlock (FullOverlapMatrixB, CurrentPLevel, QValueB, CurrentPLevel, QValueB);
	  
	}
      
    }
 
  return 0;
}

