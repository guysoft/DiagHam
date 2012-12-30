#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/FermionOnSphereMPSWrapper.h"
#include "HilbertSpace/FermionOnCylinderMPSWrapper.h"

#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "MathTools/ClebschGordanCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEMPS/FQHEMPSMatrixManager.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Hamiltonian/TensorProductSparseMatrixHamiltonian.h"

#include "LanczosAlgorithm/BasicArnoldiAlgorithm.h"

#include "Matrix/SparseRealMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "GeneralTools/ArrayTools.h"

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




void MPSDiagonalizeEMatrix(AbstractHamiltonian* hamiltonian, int nbrEigenstates, Complex*& eigenvalues, ComplexVector*& eigenstates, 
			   AbstractArchitecture* architecture, int fullDIagonalizationLimit, double error, bool leftFlag);


int main(int argc, char** argv)
{
  cout.precision(14); 
  
  OptionManager Manager ("FQHESphereMPSEntanglementSpectrum" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  FQHEMPSMatrixManager MPSMatrixManager;

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "file that describes the root configuration");
  (*SystemGroup) += new BooleanOption  ('\n', "use-padding", "root partitions use the extra zero padding");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "la", "number of orbitals in subsystem A", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "na", "number of particles in subsystem A", 0);
  (*SystemGroup) += new BooleanOption ('\n', "infinite-cylinder", "evaluate the entnaglement spectrum on the infinite cylinder");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "memory", "amount of memory that can used for precalculations (in Mb)", 500);
  (*OutputGroup) += new SingleStringOption  ('o', "output-file", "output file name");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereMPSEntanglementSpectrum -h" << endl;
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
  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
    return -1;

  int EntCut = Manager.GetInteger("la");
  int Na = Manager.GetInteger("na");

  bool CylinderFlag = Manager.GetBoolean("normalize-cylinder");
  double AspectRatio = Manager.GetDouble("aspect-ratio");
  double kappa = 0.0;
  if (CylinderFlag)
    {
       kappa = (2.0 * M_PI)/sqrt(2.0 * M_PI * (NbrFluxQuanta + 1) * AspectRatio);
       cout << "Cylinder geometry, kappa= "<<kappa<<endl;
    }

  int LandauLevel = 0;

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  int NbrBMatrices = 2;
  SparseRealMatrix* BMatrices = MPSMatrix->GetMatrices();
  SparseRealMatrix* ConjugateBMatrices = new SparseRealMatrix[NbrBMatrices];
  for (int i = 0; i < NbrBMatrices; ++i)
    ConjugateBMatrices[i] = BMatrices[i].Transpose();

  cout << "B matrix size = " << BMatrices[0].GetNbrRow() << "x" << BMatrices[0].GetNbrColumn() << endl;
  
  int MPSRowIndex = 0;
  int MPSColumnIndex = 0;
  int NbrEigenstates = 0;
  if (Manager.GetBoolean("use-padding") == true)
    {
      if (Manager.GetBoolean("k-2") == true)
	{
	  if ((Manager.GetInteger("r-index") & 1) == 0)
	    MPSRowIndex = Manager.GetInteger("p-truncation") + (Manager.GetInteger("r-index") / 2);
	  else
	    MPSRowIndex = 2 * Manager.GetInteger("p-truncation") + Manager.GetInteger("r-index") - 1;
	  NbrEigenstates = Manager.GetInteger("r-index") + 2;
	}
      else
	{
	  if (Manager.GetBoolean("rr-3") == true)
	    {
	      MPSRowIndex = 3 * (Manager.GetInteger("p-truncation") + 1);
	      NbrEigenstates = 5;
	    }
	  else
	    {
	      MPSRowIndex = Manager.GetInteger("p-truncation") + ((Manager.GetInteger("laughlin-index") - 1) / 2);
	      NbrEigenstates = Manager.GetInteger("laughlin-index");
	    }
	}
      MPSColumnIndex = MPSRowIndex;
    }
  else
    {
      if (Manager.GetBoolean("k-2") == true)
	{
	  if ((Manager.GetInteger("r-index") & 1) == 0)
	    {
	      MPSRowIndex = Manager.GetInteger("p-truncation") + Manager.GetInteger("r-index");
	      MPSColumnIndex = Manager.GetInteger("p-truncation");
	    }
	  else
	    {
	      MPSRowIndex = 2 * (Manager.GetInteger("p-truncation") + Manager.GetInteger("r-index"));
	      MPSColumnIndex = 2 * Manager.GetInteger("p-truncation");
	    }
	  NbrEigenstates = Manager.GetInteger("r-index") + 2;
	}
      else
	{
	  if (Manager.GetBoolean("rr-3") == true)
	    {
	      MPSRowIndex = 3 * (Manager.GetInteger("p-truncation") + 2);
	      MPSColumnIndex = 3 * Manager.GetInteger("p-truncation");
	      NbrEigenstates = 5;
	    }
	  else
	    {
	      MPSRowIndex = Manager.GetInteger("p-truncation") + (Manager.GetInteger("laughlin-index") - 1);
	      MPSColumnIndex = Manager.GetInteger("p-truncation");
	      NbrEigenstates = Manager.GetInteger("laughlin-index");
	    }
	}
    }

  ofstream File;
  File.precision(14);
  
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = new char [512];
      char* StateName = new char [256];
      if (Manager.GetBoolean("k-2") == true)
	{
	  sprintf (StateName, "clustered_k_2_r_%ld", Manager.GetInteger("r-index"));
	}
      else
	{
	  if (Manager.GetBoolean("rr-3") == true)
	    {
	      sprintf (StateName, "readrezayi3");
	    }
	  else
	    {
	      sprintf (StateName, "laughlin%ld", Manager.GetInteger("laughlin-index"));
	    }
	}      
      if (CylinderFlag == true)
	{
	  if (Manager.GetBoolean("infinite-cylinder"))
	    {
	      sprintf(TmpFileName, "fermions_infinite_cylinder_%s_plevel_%ld_n_0_2s_0_lz_0.0.full.ent", StateName,
		      Manager.GetInteger("p-truncation"));
	    }
	  else
	    {
	      sprintf(TmpFileName, "fermions_cylinder_%s_plevel_%ld_n_%d_2s_%d_lz_%d.0.full.ent", StateName,
		      Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz);
	    }
	}
      else
	{
	  if (Manager.GetBoolean("infinite-cylinder"))
	    {
	      sprintf(TmpFileName, "fermions_infinite_%s_plevel_%ld_n_0_2s_0_lz_0.0.full.ent", StateName,
		      Manager.GetInteger("p-truncation"));
	    }
	  else
	    {
	      sprintf(TmpFileName, "fermions_%s_plevel_%ld_n_%d_2s_%d_lz_%d.0.full.ent", StateName,
		      Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz);
	    }
	}
      File.open(TmpFileName, ios::binary | ios::out);     
   }


  if (Manager.GetBoolean("infinite-cylinder"))
    {
      int NbrBMatrices = 2;
      double Error = 1e-10;
      int FullDiagonalizationLimit = 1000;

      SparseRealMatrix* SparseBMatrices = MPSMatrix->GetMatrices();
      SparseRealMatrix* SparseTransposeBMatrices = new SparseRealMatrix[NbrBMatrices];
      double* Coefficients = new double[NbrBMatrices];
      for (int i = 0; i < NbrBMatrices; ++i)
	{
	  Coefficients[i] = 1.0;
	  SparseTransposeBMatrices[i] = SparseBMatrices[i].Transpose();
	}
      
      
      TensorProductSparseMatrixHamiltonian ETransposeHamiltonian(NbrBMatrices, SparseBMatrices, SparseBMatrices, Coefficients);
      ComplexVector* LeftEigenstates = 0;
      Complex* LeftEigenvalues = 0;
      cout << "computing left eigenstates : " << endl;
      MPSDiagonalizeEMatrix(&ETransposeHamiltonian, NbrEigenstates, LeftEigenvalues, LeftEigenstates, Architecture.GetArchitecture(), FullDiagonalizationLimit, Error, true);


      TensorProductSparseMatrixHamiltonian EHamiltonian(NbrBMatrices, SparseTransposeBMatrices, SparseTransposeBMatrices, Coefficients);
      ComplexVector* RightEigenstates = 0;
      Complex* RightEigenvalues = 0;
      cout << "computing right eigenstates : " << endl;
      MPSDiagonalizeEMatrix(&EHamiltonian, NbrEigenstates, RightEigenvalues, RightEigenstates, Architecture.GetArchitecture(), FullDiagonalizationLimit, Error, false);

      cout << "eigenvalues : " << endl;
      for (int i = 0; i < NbrEigenstates; ++i)
	cout << LeftEigenvalues[i] << " " << RightEigenvalues[i] 
	     << "   moduli : " << Norm(LeftEigenvalues[i]) << " " << Norm(RightEigenvalues[i]) <<endl;
      cout << endl;
      
      cout << "checking scalar products between left and right E eigenstates:" << endl;
      for (int i = 0; i < NbrEigenstates; ++i)
	for (int j = 0; j < NbrEigenstates; ++j)
	  {
	    Complex Test = 0.0;
	    for (int k = 0; k < LeftEigenstates[i].GetVectorDimension(); ++k)
	      Test += LeftEigenstates[i][k] * RightEigenstates[j][k];
	    cout << "< " << i << " | " << j << " > = " << EuclidianScalarProduct(LeftEigenstates[i], RightEigenstates[j]) << " " << Test << endl;
	  }
      
      int MinQValue = 0;
      int MaxQValue = 0;
      MPSMatrix->GetChargeIndexRange(MinQValue, MaxQValue);
      char* OutputFileName = Manager.GetString("output-file");
      if (OutputFileName == 0)
	{
	  char* TmpExtension = new char[512];
	  sprintf(TmpExtension, "fermions_infinite_cylinder_n_0_2s_0_lz.0.0.full.ent");
	}
      File << "# la na lz shifted_lz lambda -log(lambda)" << endl;
      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
	{
	  for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	    {
	      int IndexRange = MPSMatrix->GetBondIndexRange(PLevel, QValue);
	      if (IndexRange >= 0)
		{
		  int TmpBMatrixDimension = SparseBMatrices[0].GetNbrRow();
		  ComplexMatrix LeftMDaggerM (IndexRange, IndexRange);		    
		  Complex* TmpFactors = new Complex [NbrEigenstates];
		  for (int i = 0;i < NbrEigenstates; ++i)
		    TmpFactors[i] = RightEigenstates[i][MPSRowIndex * TmpBMatrixDimension + MPSRowIndex] / EuclidianScalarProduct(LeftEigenstates[i], RightEigenstates[i]);
		  for (int i = 0; i < IndexRange; ++i)
		    for (int j = 0; j < IndexRange; ++j)
		      {
			LeftMDaggerM[j][i] = 0.0;
			for (int k = 0; k < NbrEigenstates; ++k)
			  LeftMDaggerM[j][i] += (LeftEigenstates[k][MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(j, PLevel, QValue) * TmpBMatrixDimension 
								    + MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(i, PLevel, QValue)] * TmpFactors[k]);
		      }
		  
		  
		  ComplexMatrix RightMDaggerM (IndexRange, IndexRange);
		  for (int i = 0;i < NbrEigenstates; ++i)
		    TmpFactors[i] = LeftEigenstates[i][MPSColumnIndex * TmpBMatrixDimension + MPSColumnIndex] / EuclidianScalarProduct(LeftEigenstates[i], RightEigenstates[i]);
		  for (int i = 0; i < IndexRange; ++i)
		    for (int j = 0; j < IndexRange; ++j)
		      {
			RightMDaggerM[j][i] = 0.0;
			for (int k = 0; k < NbrEigenstates; ++k)
			  RightMDaggerM[j][i] += (RightEigenstates[k][MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(j, PLevel, QValue) * TmpBMatrixDimension 
								      + MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(i, PLevel, QValue)] * TmpFactors[k]);
		      }
		  
		  
		  RealSymmetricMatrix SymLeftMDaggerM (LeftMDaggerM);
		  RealMatrix TmpLeftBasis(SymLeftMDaggerM.GetNbrRow(), SymLeftMDaggerM.GetNbrRow());
		  TmpLeftBasis.SetToIdentity();
		  RealDiagonalMatrix TmpLeftDiag;
#ifdef __LAPACK__
		  SymLeftMDaggerM.LapackDiagonalize(TmpLeftDiag, TmpLeftBasis);
#else
		  SymLeftMDaggerM.Diagonalize(TmpLeftDiag, TmpLeftBasis);
#endif
		  int NbrZeroLeftEigenvalues = 0;
		  for (int i = 0; i < TmpLeftDiag.GetNbrRow(); ++i)
		    {
		      if (fabs(TmpLeftDiag(i, i)) < Error)
			{
			  ++NbrZeroLeftEigenvalues;	    
			}
		    }

		  RealSymmetricMatrix SymRightMDaggerM (RightMDaggerM);
		  RealMatrix TmpRightBasis(SymRightMDaggerM.GetNbrRow(), SymRightMDaggerM.GetNbrRow());
		  RealDiagonalMatrix TmpRightDiag;
		  TmpRightBasis.SetToIdentity();
#ifdef __LAPACK__
		  SymRightMDaggerM.LapackDiagonalize(TmpRightDiag, TmpRightBasis);
#else
		  SymRightMDaggerM.Diagonalize(TmpRightDiag, TmpRightBasis);
#endif
		  int NbrZeroRightEigenvalues = 0;
		  for (int i = 0; i < TmpRightDiag.GetNbrRow(); ++i)
		    {
		      if (fabs(TmpRightDiag(i, i)) < Error)
			{
			  ++NbrZeroRightEigenvalues;	    
			}
		    }


		  
		  if ((NbrZeroLeftEigenvalues < SymLeftMDaggerM.GetNbrRow()) && (NbrZeroRightEigenvalues < SymRightMDaggerM.GetNbrRow()))
		    {
		      RealMatrix TruncatedLeftBasis (TmpLeftDiag.GetNbrRow(),TmpLeftDiag.GetNbrRow() -  NbrZeroLeftEigenvalues, true);
		      NbrZeroLeftEigenvalues = 0;
		      for (int i = 0; i < TmpLeftBasis.GetNbrColumn(); ++i)
			{
			  if (fabs(TmpLeftDiag(i, i)) > Error)
			    {
			      TruncatedLeftBasis[NbrZeroLeftEigenvalues].Copy(TmpLeftBasis[i]);
			      TruncatedLeftBasis[NbrZeroLeftEigenvalues] *= sqrt(TmpLeftDiag(i, i));
			      ++NbrZeroLeftEigenvalues;
			    }
			}
		      
		      RealMatrix TruncatedRightBasis (TmpRightDiag.GetNbrRow(),TmpRightDiag.GetNbrRow() -  NbrZeroRightEigenvalues, true);
		      NbrZeroRightEigenvalues = 0;
		      for (int i = 0; i < TmpRightBasis.GetNbrColumn(); ++i)
			{
			  if (fabs(TmpRightDiag(i, i)) > Error)
			    {
			      TruncatedRightBasis[NbrZeroRightEigenvalues].Copy(TmpRightBasis[i]);
			      TruncatedRightBasis[NbrZeroRightEigenvalues] *= sqrt(TmpRightDiag(i, i));
			      ++NbrZeroRightEigenvalues;
			    }
			}
		      
		      RealMatrix TranposedTruncatedRightBasis = TruncatedRightBasis.DuplicateAndTranspose();
		      TruncatedRightBasis.Multiply(TranposedTruncatedRightBasis);		      
		      RealMatrix TranposedTruncatedLeftBasis = TruncatedLeftBasis.DuplicateAndTranspose();
		      TranposedTruncatedLeftBasis.Multiply(TruncatedRightBasis);
		      TranposedTruncatedLeftBasis.Multiply(TruncatedLeftBasis);
		      
		      RealSymmetricMatrix ReducedDensityMatrix (TruncatedLeftBasis.GetNbrColumn(), true);
		      
		      if (ReducedDensityMatrix.IsDiagonal() == true)
			{
			  cout << "the reduced density matrix is diagonal" << endl;
			}
		      else
			{
			  cout << "the reduced density matrix is not diagonal" << endl;
			}
		      int NbrNonZeroEigenvalues = 0;
		      for (int i = 0; i < TranposedTruncatedLeftBasis.GetNbrColumn(); ++i)
			{
			  if (fabs(TranposedTruncatedLeftBasis(i, i)) > Error)
			    {
			      File << "0 " << QValue << " " << PLevel << " " << PLevel << " " << TranposedTruncatedLeftBasis(i, i) << endl;
			      ++NbrNonZeroEigenvalues;
			    }
			}
		      cout << PLevel << " " << QValue << " " << NbrNonZeroEigenvalues << endl;
		    }
		}
	    }
	}
      File.close();
      return 0;
    }

  int MatDim = BMatrices[0].GetNbrRow();
  int LambdaMax = Manager.GetInteger("p-truncation");
  int LaughlinIndex = Manager.GetInteger("laughlin-index");

  cout << "B matrix size = " << MatDim << "x" << MatDim << endl;

  SparseRealMatrix* FinalBMatrices = new SparseRealMatrix[NbrBMatrices];
  for (int j = 0; j < NbrBMatrices; ++j)
    FinalBMatrices[j] = SparseRealMatrix (MatDim, MatDim);

  for(int i = 0; i < MatDim; i++)
   {
    double Tmp;
    for (int j = 0; j < NbrBMatrices; ++j)
      {
	BMatrices[j].GetMatrixElement(i, MPSColumnIndex, Tmp);
	if (Tmp != 0.0)
	  FinalBMatrices[j].SetMatrixElement(i, MPSColumnIndex, Tmp);
      }
   }

  cout<<"Done preparing B matrices and the vectors at 0 and Nphi orbital"<<endl;

  double CutOff = 1e-14;
  double TmpTrace;

  cout<<"Proceed to calculate overlap matrix (full space dimension, and it will be stored) "<<endl;

  long MaxTmpMatrixElements = (((long) BMatrices[0].GetNbrRow()) * 
				((long) BMatrices[0].GetNbrRow() / 1l));
  if (MaxTmpMatrixElements  > (1l << 28))
    MaxTmpMatrixElements  = 1l << 28; 
  cout << "Requested memory for sparse matrix multiplications = " << ((MaxTmpMatrixElements * (2l * sizeof(double) + sizeof(int))) >> 20) << "Mb" << endl;
  double* TmpMatrixElements = new double [MaxTmpMatrixElements];
  int* TmpColumnIndices = new int [MaxTmpMatrixElements];
  double * TmpElements = new double [BMatrices[0].GetNbrRow()];

  double* NormalizationCoefficients = new double[NbrFluxQuanta + 1];
  BinomialCoefficients Binomial(NbrFluxQuanta);
  for (int i = 0; i <= NbrFluxQuanta; ++i)
    {      
      NormalizationCoefficients[i] = ((double) (NbrFluxQuanta + 1)) / Binomial.GetNumericalCoefficient(NbrFluxQuanta, i);
    }

  SparseRealMatrix OverlapMatrix (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
  OverlapMatrix.SetMatrixElement(MPSRowIndex, MPSRowIndex, 1.0);
  SparseRealMatrix TmpMatrix1;
  SparseRealMatrix TmpMatrix2;
  SparseRealMatrix TmpMatrix3;
  for (int i = 0; i < EntCut; i++)
    {
      TmpMatrix1 = Multiply(ConjugateBMatrices[0], OverlapMatrix, 
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix2 = Multiply(TmpMatrix1, BMatrices[0],  
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix1 = Multiply(ConjugateBMatrices[1], OverlapMatrix, 
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix3 = Multiply(TmpMatrix1, BMatrices[1],  
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 

      if (CylinderFlag == false)
	TmpMatrix3 *= NormalizationCoefficients[i];
      
      OverlapMatrix = TmpMatrix2 + TmpMatrix3;
      
      TmpTrace = OverlapMatrix.Tr();
      if (fabs(TmpTrace) > CutOff)
        OverlapMatrix /= TmpTrace;
      else
        {
          cout << "Warning: trying to normalize with 0! " << endl;
          exit(1);
        } 
      
    }


  cout<<"Compute density matrix in nonorthogonal basis (full space dimension, will be stored)"<<endl;

  SparseRealMatrix RhoA;

  TmpMatrix2 = Multiply(FinalBMatrices[0], FinalBMatrices[0].Transpose(), 
 			TmpMatrixElements, TmpColumnIndices, TmpElements); 

  if (FinalBMatrices[1].ComputeNbrNonZeroMatrixElements() != 0)
   {
     TmpMatrix3 = Multiply(FinalBMatrices[1], FinalBMatrices[1].Transpose(), 
			   TmpMatrixElements, TmpColumnIndices, TmpElements); 
      if (CylinderFlag == false)
	TmpMatrix3 *= NormalizationCoefficients[NbrFluxQuanta];
     RhoA = TmpMatrix2 + TmpMatrix3;
   }
  else
    {
      RhoA = TmpMatrix2;
    }

  for (int i = (NbrFluxQuanta - 1); i >= EntCut; i--)
    {
      TmpMatrix1 = Multiply(BMatrices[0], RhoA, 
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix2 = Multiply(TmpMatrix1, ConjugateBMatrices[0],  
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix1 = Multiply(BMatrices[1], RhoA, 
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix3 = Multiply(TmpMatrix1, ConjugateBMatrices[1],  
			    TmpMatrixElements, TmpColumnIndices, TmpElements);   
      if (CylinderFlag == false)
	TmpMatrix3 *= NormalizationCoefficients[i];
      RhoA = TmpMatrix2 + TmpMatrix3;

      TmpTrace = RhoA.Tr();
      if (fabs(TmpTrace) > CutOff)
        RhoA /= TmpTrace;
      else
        {
          cout << "Warning: trying to normalize with 0! " << endl;
          exit(1);
        } 
    }

  //Free up some space that is no longer needed (this needs to be done in a cleaner way)

  delete[] BMatrices; 
  delete[] ConjugateBMatrices; 
  delete[] FinalBMatrices;
  delete[] TmpElements;

  cout<<"Proceed to calculate ES per momentum P sector (all N sectors)"<<endl;

  double TraceRho = 0.0;
  double* RhoEigenvalues = new double [MatDim];
  int* RhoNSector = new int [MatDim];
  int* RhoPSector = new int [MatDim];

  int MinNValue;
  int MaxNValue;
  MPSMatrix->GetChargeIndexRange(MinNValue, MaxNValue);
  int NbrNValue = MaxNValue - MinNValue + 1;

  for (int i =0 ; i < MatDim; i++)
   {
     RhoEigenvalues[i] = 0.0; 
     RhoNSector[i] = 0; 
     RhoPSector[i] = 0;
   }

  long NbrEigenvalues = 0l;
  OverlapMatrix.PrintNonZero(cout) << endl;
  for (int NSector = 0; NSector < NbrNValue; NSector++)
    {
      for (int MomentumSector = 0; MomentumSector <= LambdaMax; MomentumSector++)
	{
	  SparseRealMatrix TmpOverlapBlock = MPSMatrix->ExtractBlock(OverlapMatrix, MomentumSector, NSector, MomentumSector, NSector);
	  SparseRealMatrix RhoABlock = MPSMatrix->ExtractBlock(RhoA, MomentumSector, NSector, MomentumSector, NSector);

	  if ((TmpOverlapBlock.ComputeNbrNonZeroMatrixElements() != 0) && (RhoABlock.ComputeNbrNonZeroMatrixElements()))
	    {
 	      cout << "NSector=" << NSector << "  MomentumSector=" << MomentumSector << " : "<< endl;
	      // 	      TmpOverlapBlock.PrintNonZero(cout) << endl;
	      int TmpSectorDim = TmpOverlapBlock.GetNbrRow();
	      
	      RealSymmetricMatrix HRep (TmpOverlapBlock);
	      RealDiagonalMatrix TmpDiag (TmpSectorDim);
	      RealMatrix Q(TmpSectorDim, TmpSectorDim);
	      HRep.LapackDiagonalize(TmpDiag, Q);
	      
	      int NbrNonZeroVectors = 0;
	      for (int j = 0; j < TmpSectorDim; ++j)
		{
		  if (fabs(TmpDiag[j]) > CutOff)
		    {
		      ++NbrNonZeroVectors;
		    }
		}
	      if (NbrNonZeroVectors > 0)
		{
		  RealMatrix BlockBasisLeftMatrix (TmpSectorDim, NbrNonZeroVectors);
		  RealMatrix BlockBasisRightMatrix (TmpSectorDim, NbrNonZeroVectors);
		  NbrNonZeroVectors = 0;
		  for (int j = 0; j < TmpSectorDim; ++j)
		    {
		      if (fabs(TmpDiag[j]) > CutOff)
			{  			
			  BlockBasisLeftMatrix[NbrNonZeroVectors].Copy(Q[j]);
			  if (TmpDiag[j] > 0.0)
			    BlockBasisLeftMatrix[NbrNonZeroVectors] /= sqrt(fabs(TmpDiag[j]));
			  else
			    BlockBasisLeftMatrix[NbrNonZeroVectors] /= -sqrt(fabs(TmpDiag[j]));
			  BlockBasisRightMatrix[NbrNonZeroVectors].Copy(Q[j]);
			  BlockBasisRightMatrix[NbrNonZeroVectors] /= sqrt(fabs(TmpDiag[j]));
			  ++NbrNonZeroVectors;
			}  
		    }
	      
		  
		  //cout<<"-------------------- Start computing rho in the new basis --------------------"<<endl;

		  TmpElements = new double [TmpOverlapBlock.GetNbrRow()];
		  SparseRealMatrix TmpMat = Conjugate(TmpOverlapBlock, RhoABlock, TmpOverlapBlock.Transpose(),
						      TmpMatrixElements, TmpColumnIndices, TmpElements);
		  delete[] TmpElements;
		  RealSymmetricMatrix HRepRho = TmpMat.Conjugate(BlockBasisLeftMatrix, BlockBasisRightMatrix);

		  RealDiagonalMatrix TmpDiagRho (NbrNonZeroVectors);
		  HRepRho.LapackDiagonalize(TmpDiagRho);
		  
		  cout<<"------------sector P = "<<MomentumSector<<" N = "<<((EntCut - (NSector - (2 * LambdaMax + LaughlinIndex - 1)/2)))/LaughlinIndex << "---------------"<<endl;
		  
		  
		  double Sum = 0.0;
		  for (int j = 0; j < NbrNonZeroVectors; ++j)
		    {
		      TraceRho += TmpDiagRho[j];
		      RhoEigenvalues[NbrEigenvalues]=TmpDiagRho[j];
		      RhoNSector[NbrEigenvalues]=NSector;
		      RhoPSector[NbrEigenvalues]=MomentumSector; 
		      NbrEigenvalues++;
		    }
		}
	    }
	}
    }
  
  cout<<"Trace rho = "<<TraceRho<<endl;

  File << "# l_a    N    Lz    lambda" << endl;
  for (int i=0; i<MatDim; ++i)
    if (((fabs(RhoEigenvalues[i]) > CutOff)) && ((((EntCut - (RhoNSector[i] - (2 * LambdaMax + LaughlinIndex - 1)/2)))/LaughlinIndex) == Na))
      {
        cout<< "Na= " << Na << " P= " << RhoPSector[i] << " " << RhoEigenvalues[i]/TraceRho << endl;  
        File << EntCut << " " << Na << " " << RhoPSector[i] << " " << (RhoEigenvalues[i] / TraceRho) << endl;
      }

  delete[] TmpMatrixElements;
  delete[] TmpColumnIndices;
  delete[] RhoEigenvalues;
  delete[] RhoPSector;
  delete[] RhoNSector;

  File.close();
 
  return 0;
}


void MPSDiagonalizeEMatrix(AbstractHamiltonian* hamiltonian, int nbrEigenstates, Complex*& eigenvalues, ComplexVector*& eigenstates, 
			   AbstractArchitecture* architecture, int fullDiagonalizationLimit, double error, bool leftFlag)
{
  eigenstates = new ComplexVector[nbrEigenstates];
  eigenvalues = new Complex[nbrEigenstates];
  double* TmpRealPart = new double [nbrEigenstates];
  int* TmpIndices = new int[nbrEigenstates];
  Complex* TmpEigenvalues = 0;
  ComplexVector* TmpEigenstates = 0;

  if (hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension() < fullDiagonalizationLimit)
    {	  
      RealMatrix HRepresentation (hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension(), 
				  hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension());
      hamiltonian->GetHamiltonian(HRepresentation);
      ComplexDiagonalMatrix TmpDiag (HRepresentation.GetNbrRow(), true);  
      ComplexMatrix TmpEigenstateMatrix (HRepresentation.GetNbrRow(), HRepresentation.GetNbrRow());  
      HRepresentation.LapackDiagonalize(TmpDiag, TmpEigenstateMatrix, true);
      TmpDiag.SortMatrixDownOrder(TmpEigenstateMatrix, true);
      
      TmpEigenvalues = new Complex [nbrEigenstates];
      TmpEigenstates = new ComplexVector [nbrEigenstates];
      for (int i = 0; i < nbrEigenstates; ++i)
	{
	  TmpEigenvalues[i] = TmpDiag[i];
	  TmpEigenstates[i] = TmpEigenstateMatrix[i];
	}
    }
  else
    {	 
      BasicArnoldiAlgorithm Arnoldi(architecture, nbrEigenstates, 3000, true, false, false);
      Arnoldi.SetHamiltonian(hamiltonian);
      Arnoldi.InitializeLanczosAlgorithm();
      Arnoldi.RunLanczosAlgorithm(nbrEigenstates);
      while (Arnoldi.TestConvergence() == false)
	{
	  Arnoldi.RunLanczosAlgorithm(1);
	}
      TmpEigenstates = (ComplexVector*) Arnoldi.GetEigenstates(nbrEigenstates);  
      Arnoldi.GetEigenvalues(TmpEigenvalues, nbrEigenstates);
     }

  for (int i = 0; i < nbrEigenstates; ++i)
    {
      TmpRealPart[i] = TmpEigenvalues[i].Re;
      TmpIndices[i] = i;
    }
  SortArrayDownOrdering<int>(TmpRealPart, TmpIndices, nbrEigenstates);
  for (int i = 1; i < nbrEigenstates; ++i)
    {
      if ((fabs(TmpEigenvalues[TmpIndices[i - 1]].Re - TmpEigenvalues[TmpIndices[i]].Re) < error) && 
	  (TmpEigenvalues[TmpIndices[i - 1]].Im < TmpEigenvalues[TmpIndices[i]].Im))
	{
	  int Tmp = TmpIndices[i - 1];
	  TmpIndices[i - 1] = TmpIndices[i];
	  TmpIndices[i] = Tmp;
	}
    }
  for (int i = 0; i < nbrEigenstates; ++i)
    {
      eigenstates[i] = TmpEigenstates[TmpIndices[i]];
      eigenvalues[i] = TmpEigenvalues[TmpIndices[i]];
    }
  for (int i = 0; i < nbrEigenstates; ++i)
    {
      ComplexVector TestE (eigenstates[i].GetVectorDimension());
      hamiltonian->Multiply(eigenstates[i], TestE);
      cout << (TestE * eigenstates[i]) << " " << eigenvalues[i] << endl;
    }
  delete[] TmpRealPart;
  delete[] TmpIndices;
  delete[] TmpEigenstates;
  delete[] TmpEigenvalues;
}
