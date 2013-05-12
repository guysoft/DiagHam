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

   int Na = Manager.GetInteger("na");

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
      File << "#  N    Lz    lambda" << endl;
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

  cout << "MPSRowIndex = " << MPSRowIndex << endl;
  cout << "MPSColumnIndex = " << MPSColumnIndex << endl;
  
//   SparseRealMatrix TmpNormalizationMatrix (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
//   TmpNormalizationMatrix.SetMatrixElement(MPSRowIndex, MPSRowIndex, 1.0);  
//   for (int i = 0; i <= NbrFluxQuanta; ++i)
//     {
//       SparseRealMatrix TmpMatrix2;
//       SparseRealMatrix TmpMatrix3;
//       TmpMatrix2 = Conjugate(ConjugateBMatrices[0], TmpNormalizationMatrix, BMatrices[0], 
// 			     TmpMatrixElements, TmpColumnIndices, TmpElements); 
//       TmpMatrix3 = Conjugate(ConjugateBMatrices[1], TmpNormalizationMatrix, BMatrices[1], 
// 			     TmpMatrixElements, TmpColumnIndices, TmpElements); 
//       TmpNormalizationMatrix  = TmpMatrix2 + TmpMatrix3;
//     }
//   double Normalization;
//   TmpNormalizationMatrix.GetMatrixElement(MPSColumnIndex, MPSColumnIndex, Normalization);
//   cout << "Normalization = " << Normalization << endl;

  SparseRealMatrix FullLeftOverlapMatrix (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
  FullLeftOverlapMatrix.SetMatrixElement(MPSRowIndex, MPSRowIndex, 1.0);  
  int MaxNbrFluxQuantaA = PLevel + QValue * (Na - 1) + ((QValue - 1) / 2);
  if (MaxNbrFluxQuantaA > NbrFluxQuanta)
    {
      MaxNbrFluxQuantaA = NbrFluxQuanta;
    }
  for (int i = 0; i <= MaxNbrFluxQuantaA; ++i)
    {
      SparseRealMatrix TmpMatrix2;
      SparseRealMatrix TmpMatrix3;
      TmpMatrix2 = Conjugate(ConjugateBMatrices[0], FullLeftOverlapMatrix, BMatrices[0], 
			     TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix3 = Conjugate(ConjugateBMatrices[1], FullLeftOverlapMatrix, BMatrices[1], 
			     TmpMatrixElements, TmpColumnIndices, TmpElements); 
      FullLeftOverlapMatrix = TmpMatrix2 + TmpMatrix3;
    }

  SparseRealMatrix FullRightOverlapMatrix (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
  FullRightOverlapMatrix.SetMatrixElement(MPSColumnIndex, MPSColumnIndex, 1.0);  
  int MaxNbrFluxQuantaB = PLevel + QValue * (NbrParticles - Na - 1)+ ((QValue - 1) / 2);
  if (MaxNbrFluxQuantaB > NbrFluxQuanta)
    {
      MaxNbrFluxQuantaB = NbrFluxQuanta;
    }
  for (int i = 0; i <= MaxNbrFluxQuantaB; ++i)
    {
      SparseRealMatrix TmpMatrix2;
      SparseRealMatrix TmpMatrix3;
      TmpMatrix2 = Conjugate(BMatrices[0], FullRightOverlapMatrix, ConjugateBMatrices[0], 
			     TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix3 = Conjugate(BMatrices[1], FullRightOverlapMatrix, ConjugateBMatrices[1], 
			     TmpMatrixElements, TmpColumnIndices, TmpElements); 
      FullRightOverlapMatrix = TmpMatrix2 + TmpMatrix3;
    }
  
  SparseRealMatrix NormalizationMatrix (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
  NormalizationMatrix.SetToIdentity();
  int TmpPower = 2 * PLevel;
  for (int i = 0; i < TmpPower; ++i)
    {
      NormalizationMatrix.Multiply(BMatrices[0]);
    }
  double Error = 1e-13;
  double LeftEigenvalueError = 0.0;
  double RightEigenvalueError = 0.0;
  LeftEigenvalueError = Error;
  RightEigenvalueError = Error;
  double TotalTraceThoA = 0;
  int MinQValue;
  int MaxQValue;
  MPSMatrix->GetChargeIndexRange(0, MinQValue, MaxQValue);
  for (int CurrentPLevel = 1; CurrentPLevel <= PLevel; ++CurrentPLevel)
    {
      int TmpMinQValue;
      int TmpMaxQValue;
      MPSMatrix->GetChargeIndexRange(CurrentPLevel, TmpMinQValue, TmpMaxQValue);
    }
  double** EntanglementSpectrum = new double*[PLevel + 1];
  int* EntanglementSpectrumDimension = new int[PLevel + 1];
  for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
    EntanglementSpectrumDimension[CurrentPLevel] = 0;


//   FullLeftOverlapMatrix.PrintNonZero(cout) << endl;
//   FullRightOverlapMatrix.PrintNonZero(cout) << endl;
//  cout << NormalizationMatrix << endl;
  for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
    {
      cout << "computing level " << CurrentPLevel << endl;
      int LeftQValue = PLevel + PLevel + (QValue - 1) / 2; // charge Q=+1
      int RightQValue = PLevel - PLevel + (QValue - 1) / 2; // charge Q=+1
      SparseRealMatrix RightOverlapMatrix = MPSMatrix->ExtractBlock(FullRightOverlapMatrix, CurrentPLevel, RightQValue, CurrentPLevel, RightQValue);
      SparseRealMatrix LocalNormalizationMatrix = MPSMatrix->ExtractBlock(NormalizationMatrix, CurrentPLevel, RightQValue, CurrentPLevel, LeftQValue);
      for (int i = 0; i < LocalNormalizationMatrix.GetNbrRow(); ++i)
	{
	  double Tmp;
	  LocalNormalizationMatrix.GetMatrixElement(i, i, Tmp);
	  LocalNormalizationMatrix.SetMatrixElement(i, i, 1.0 / Tmp);
	}
      SparseRealMatrix LeftOverlapMatrix =MPSMatrix->ExtractBlock(FullLeftOverlapMatrix, CurrentPLevel, LeftQValue, CurrentPLevel, LeftQValue);
      cout << LocalNormalizationMatrix << endl;
      cout << LeftOverlapMatrix << endl;
      cout << RightOverlapMatrix << endl;
      if ((LeftOverlapMatrix.GetNbrRow() > 0) && (RightOverlapMatrix.GetNbrRow() > 0))
	{
	  cout << "scalar product matrix for the left part : " << endl;
	  RealSymmetricMatrix SymLeftOverlapMatrix (LeftOverlapMatrix);
	  RealMatrix TmpLeftBasis(SymLeftOverlapMatrix.GetNbrRow(), SymLeftOverlapMatrix.GetNbrRow());
	  TmpLeftBasis.SetToIdentity();
	  RealDiagonalMatrix TmpLeftDiag;
#ifdef __LAPACK__
	  SymLeftOverlapMatrix.LapackDiagonalize(TmpLeftDiag, TmpLeftBasis);
#else
	  SymLeftOverlapMatrix.Diagonalize(TmpLeftDiag, TmpLeftBasis);
#endif
	  double LocalLeftEigenvalueError = 0.0;
	  for (int i = 0; i < TmpLeftDiag.GetNbrColumn(); ++i)
	    if (TmpLeftDiag(i, i) > LocalLeftEigenvalueError)
	      LocalLeftEigenvalueError = TmpLeftDiag(i, i);
	  LocalLeftEigenvalueError *= Error;
	  int NbrZeroLeftEigenvalues = 0;
	  for (int i = 0; i < TmpLeftDiag.GetNbrRow(); ++i)
	    {
	      if (TmpLeftDiag(i, i) < LocalLeftEigenvalueError)
		{
		  ++NbrZeroLeftEigenvalues;	    
		}
	    }
	  cout << "nbr non zero eigenvalues = " << (TmpLeftDiag.GetNbrRow() - NbrZeroLeftEigenvalues) << " (full dim = " << TmpLeftDiag.GetNbrRow() << ")" << endl;
	  
	  RealSymmetricMatrix SymRightOverlapMatrix (RightOverlapMatrix);
	  RealMatrix TmpRightBasis(SymRightOverlapMatrix.GetNbrRow(), SymRightOverlapMatrix.GetNbrRow());
	  RealDiagonalMatrix TmpRightDiag;
	  TmpRightBasis.SetToIdentity();
#ifdef __LAPACK__
	  SymRightOverlapMatrix.LapackDiagonalize(TmpRightDiag, TmpRightBasis);
#else
	  SymRightOverlapMatrix.Diagonalize(TmpRightDiag, TmpRightBasis);
#endif
	  int NbrZeroRightEigenvalues = 0;
	  cout << "scalar product matrix for the right part : " << endl;
	  double LocalRightEigenvalueError = 0.0;
	  for (int i = 0; i < TmpRightDiag.GetNbrColumn(); ++i)
	    if (TmpRightDiag(i, i) > LocalRightEigenvalueError)
	      LocalRightEigenvalueError = TmpRightDiag(i, i);
	  LocalRightEigenvalueError *= Error;
	  for (int i = 0; i < TmpRightDiag.GetNbrRow(); ++i)
	    {
	      if (TmpRightDiag(i, i) < LocalRightEigenvalueError)
		{
		  ++NbrZeroRightEigenvalues;	    
		}
	    }
	  cout << "nbr non zero eigenvalues = " << (TmpRightDiag.GetNbrRow() - NbrZeroRightEigenvalues) << " (full dim = " << TmpRightDiag.GetNbrRow() << ")"  << endl;
	  if ((NbrZeroLeftEigenvalues < SymLeftOverlapMatrix.GetNbrRow()) && (NbrZeroRightEigenvalues < SymRightOverlapMatrix.GetNbrRow()))
	    {
	      RealMatrix TruncatedLeftBasis (TmpLeftDiag.GetNbrRow(), TmpLeftDiag.GetNbrRow() -  NbrZeroLeftEigenvalues, true);
	      NbrZeroLeftEigenvalues = 0;
	      for (int i = 0; i < TmpLeftBasis.GetNbrColumn(); ++i)
		{
		  if (TmpLeftDiag(i, i) > LocalLeftEigenvalueError)
//		  if (fabs(TmpLeftDiag(i, i)) > LeftEigenvalueError)
		    {
		      TruncatedLeftBasis[NbrZeroLeftEigenvalues].Copy(TmpLeftBasis[i]);
		      TruncatedLeftBasis[NbrZeroLeftEigenvalues] *= sqrt(TmpLeftDiag(i, i));
		      ++NbrZeroLeftEigenvalues;
		    }
		}
	      
	      RealMatrix TruncatedRightBasis (TmpRightDiag.GetNbrRow(), TmpRightDiag.GetNbrRow() -  NbrZeroRightEigenvalues, true);
	      NbrZeroRightEigenvalues = 0;
	      for (int i = 0; i < TmpRightBasis.GetNbrColumn(); ++i)
		{
		  cout << TmpRightDiag(i, i) << endl;
//		  if (fabs(TmpRightDiag(i, i)) > RightEigenvalueError)
		  if (TmpRightDiag(i, i) > LocalRightEigenvalueError)
		    {
		      TruncatedRightBasis[NbrZeroRightEigenvalues].Copy(TmpRightBasis[i]);
		      TruncatedRightBasis[NbrZeroRightEigenvalues] *= sqrt(TmpRightDiag(i, i));
		      ++NbrZeroRightEigenvalues;
		    }
		}
	      
	      
	      RealDiagonalMatrix TmpRhoADiag;
	      RealMatrix TmpLocalNormalizationMatrix (LocalNormalizationMatrix);
// 	      for (int i = 0; i < TmpLocalNormalizationMatrix.GetNbrRow(); ++i)
// 		{
// 		  TmpLocalNormalizationMatrix[i][i] = 1.0 / TmpLocalNormalizationMatrix[i][i];
// 		}
	      if (false)
		{
		  RealMatrix TranposedTruncatedRightBasis = TruncatedRightBasis.DuplicateAndTranspose();
		  TruncatedRightBasis.Multiply(TranposedTruncatedRightBasis);		      
		  RealMatrix TranposedTruncatedLeftBasis = TruncatedLeftBasis.DuplicateAndTranspose();
		  TranposedTruncatedLeftBasis.Multiply(TruncatedRightBasis);
		  TranposedTruncatedLeftBasis.Multiply(TruncatedLeftBasis);
		  
		  RealSymmetricMatrix ReducedDensityMatrix ((Matrix&) TranposedTruncatedLeftBasis);
	      
		  if (ReducedDensityMatrix.IsDiagonal() == true)
		    {
		      TmpRhoADiag = RealDiagonalMatrix(TranposedTruncatedLeftBasis);
		    }
		  else
		    {
#ifdef __LAPACK__
		      ReducedDensityMatrix.LapackDiagonalize(TmpRhoADiag);
#else
		      ReducedDensityMatrix.Diagonalize(TmpRhoADiag);
#endif
		    }
		}
	      else
		{
		  cout << TruncatedLeftBasis  << endl;
		  cout <<  TruncatedRightBasis << endl;
		  RealMatrix TmpEntanglementMatrix = TruncatedLeftBasis.DuplicateAndTranspose();
		  TmpEntanglementMatrix.Multiply(TmpLocalNormalizationMatrix);
		  TmpEntanglementMatrix.Multiply(TruncatedRightBasis);
		  cout <<  TmpEntanglementMatrix << endl;
#ifdef __LAPACK__
		  double* TmpValues = TmpEntanglementMatrix.SingularValueDecomposition();
		  int TmpDimension = TmpEntanglementMatrix.GetNbrColumn();
		  if (TmpDimension > TmpEntanglementMatrix.GetNbrRow())
		    {
		      TmpDimension = TmpEntanglementMatrix.GetNbrRow();
		    }
		  for (int i = 0; i < TmpDimension; ++i)
		    {
		      TmpValues[i] *= TmpValues[i];
		    }
		  TmpRhoADiag = RealDiagonalMatrix(TmpValues, TmpDimension);
#endif		  
		}
	      int NbrNonZeroEigenvalues = 0;
	      double Sum = 0.0;
	      EntanglementSpectrum[CurrentPLevel] = new double[TmpRhoADiag.GetNbrColumn()];
	      for (int i = 0; i < TmpRhoADiag.GetNbrColumn(); ++i)
		{
		  if (TmpRhoADiag[i] > 0.0)
		    {
		      cout << "xi = " << TmpRhoADiag[i] << endl;
		      EntanglementSpectrum[CurrentPLevel][EntanglementSpectrumDimension[CurrentPLevel]] = TmpRhoADiag[i];
		      ++EntanglementSpectrumDimension[CurrentPLevel];
		      Sum += TmpRhoADiag(i, i);
		      ++NbrNonZeroEigenvalues;
		    }
		}
	      SortArrayDownOrdering<double>(EntanglementSpectrum[CurrentPLevel], EntanglementSpectrumDimension[CurrentPLevel]);
	      cout << "P=" << PLevel << " " << " Q=" << QValue << " NbrStates=" << NbrNonZeroEigenvalues << " Tr(rho_A)=" << Sum << endl;
	      TotalTraceThoA += Sum;
	    }
	}
    }
	
  for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
    {
      for (int i = 0; i < EntanglementSpectrumDimension[CurrentPLevel]; ++i)
	File << Na << " " << CurrentPLevel << " " 
	     <<  (EntanglementSpectrum[CurrentPLevel][i] / EntanglementSpectrum[0][0])  
	     <<  " " << (-log(EntanglementSpectrum[CurrentPLevel][i] / TotalTraceThoA)) << endl;
//	     <<  (EntanglementSpectrum[CurrentPLevel][i] / EntanglementSpectrum[0][0] * (50.0 / 3.0))  
    }
  File.close();

  cout << "Tr(rho_A)=" << TotalTraceThoA << endl;
  
  return 0;
}

