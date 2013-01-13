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




void MPSDiagonalizeEMatrix(OptionManager* manager, AbstractHamiltonian* hamiltonian, int nbrEigenstates, Complex*& eigenvalues, ComplexVector*& eigenstates, 
			   AbstractArchitecture* architecture, double error, bool leftFlag);


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
  OptionGroup* ArnoldiGroup  = new OptionGroup ("Arnoldi options");
  Architecture.AddOptionGroup(&Manager);
  Manager += ArnoldiGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "file that describes the root configuration");
  (*SystemGroup) += new BooleanOption  ('\n', "use-padding", "root partitions use the extra zero padding");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "la", "number of orbitals in subsystem A", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "na", "number of particles in subsystem A", 0);
  (*SystemGroup) += new BooleanOption ('\n', "infinite-cylinder", "evaluate the entnaglement spectrum on the infinite cylinder");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "memory", "amount of memory that can used for precalculations (in Mb)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "ematrix-memory", "amount of memory that can used for precalculations of the E matrix (in Mb)", 500);
  (*OutputGroup) += new SingleStringOption  ('o', "output-file", "output file name");
  (*ArnoldiGroup) += new SingleIntegerOption  ('\n', "full-diag", 
					       "maximum Hilbert space dimension for which full diagonalization is applied", 1000);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "disk", "enable disk storage for the Arnoldi algorithm", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "resume", "resume from disk datas", false);
  (*ArnoldiGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Arnoldi iteration", false); 
  (*ArnoldiGroup) += new  SingleIntegerOption ('\n', "arnoldi-memory", "amount of memory when using the Arnoldi algorithm (in Mb)", 500); 
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
	      sprintf(TmpFileName, "fermions_infinite_cylinder_%s_perimeter_%f_plevel_%ld_n_0_2s_0_lz_0.0.full.ent", StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"));
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
      double Error = 1e-13;

      SparseRealMatrix* SparseBMatrices = MPSMatrix->GetMatrices();
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
      //      SortArrayUpOrdering(EffectiveBlockIndices, EffectiveDimension);
      cout << "E matrix effective dimension = " << EffectiveDimension << "( vs " << (SparseBMatrices[0].GetNbrRow() * SparseBMatrices[0].GetNbrRow()) << ")" << endl;
      
      TensorProductSparseMatrixSelectedBlockHamiltonian* ETransposeHamiltonian = new TensorProductSparseMatrixSelectedBlockHamiltonian(NbrBMatrices, SparseBMatrices, SparseBMatrices, Coefficients, 
																       EffectiveDimension, EffectiveBlockIndices, 
																       BlockIndexProductTable, BlockIndexProductTableNbrElements, BlockIndexProductTableShift, 
																       Architecture.GetArchitecture(), Manager.GetInteger("ematrix-memory") << 20);
      ComplexVector* LeftEigenstates = 0;
      Complex* LeftEigenvalues = 0;
      cout << "computing left eigenstates : " << endl;
      MPSDiagonalizeEMatrix(&Manager, ETransposeHamiltonian, NbrEigenstates, LeftEigenvalues, LeftEigenstates, Architecture.GetArchitecture(), 1e-10, true);
      delete ETransposeHamiltonian;

      TensorProductSparseMatrixSelectedBlockHamiltonian* EHamiltonian = new TensorProductSparseMatrixSelectedBlockHamiltonian(NbrBMatrices, SparseTransposeBMatrices, SparseTransposeBMatrices, Coefficients, 
															       EffectiveDimension, EffectiveBlockIndices, 
															       BlockIndexProductTable, BlockIndexProductTableNbrElements, BlockIndexProductTableShift, 
															       Architecture.GetArchitecture(), Manager.GetInteger("ematrix-memory") << 20);
      ComplexVector* RightEigenstates = 0;
      Complex* RightEigenvalues = 0;
      cout << "computing right eigenstates : " << endl;
      MPSDiagonalizeEMatrix(&Manager, EHamiltonian, NbrEigenstates, RightEigenvalues, RightEigenstates, Architecture.GetArchitecture(), 1e-10, false);
      delete EHamiltonian;

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

      File << "# la na lz shifted_lz lambda -log(lambda)" << endl;

      
      Complex* TmpLeftFactors = new Complex [NbrEigenstates];
      Complex* TmpRightFactors = new Complex [NbrEigenstates];
      int ReducedBoundaryIndex = SearchInArray<long>((((long) MPSRowIndex) * TmpBMatrixDimension) + MPSRowIndex, 
						    EffectiveBlockIndices, EffectiveDimension);
      for (int i = 0; i < NbrEigenstates; ++i)
	TmpLeftFactors[i] = RightEigenstates[i][ReducedBoundaryIndex] / EuclidianScalarProduct(LeftEigenstates[i], RightEigenstates[i]);
      ReducedBoundaryIndex = SearchInArray<long>((((long) MPSColumnIndex) * TmpBMatrixDimension) + MPSColumnIndex, 
						EffectiveBlockIndices, EffectiveDimension);
      for (int i = 0; i < NbrEigenstates; ++i)
	TmpRightFactors[i] = LeftEigenstates[i][ReducedBoundaryIndex] / EuclidianScalarProduct(LeftEigenstates[i], RightEigenstates[i]);

      double LeftEigenvalueError = 0.0;
      double RightEigenvalueError = 0.0;

//       for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
// 	{
// 	  for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
// 	    {
// 	      int IndexRange = MPSMatrix->GetBondIndexRange(PLevel, QValue);
// 	      if (IndexRange >= 0)
// 		{
// 		  ComplexMatrix LeftMDaggerM (IndexRange, IndexRange);		    
// 		  for (int i = 0; i < IndexRange; ++i)
// 		    for (int j = 0; j < IndexRange; ++j)
// 		      {
// 			LeftMDaggerM[j][i] = 0.0;
// 			for (int k = 0; k < NbrEigenstates; ++k)
// 			  {
// 			    LeftMDaggerM[j][i] += (LeftEigenstates[k][MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(j, PLevel, QValue) * TmpBMatrixDimension 
// 								      + MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(i, PLevel, QValue)] * TmpLeftFactors[k]);
// 			  }
// 		      }
		  
		  
// 		  ComplexMatrix RightMDaggerM (IndexRange, IndexRange);
// 		  for (int i = 0; i < IndexRange; ++i)
// 		    for (int j = 0; j < IndexRange; ++j)
// 		      {
// 			RightMDaggerM[j][i] = 0.0;
// 			for (int k = 0; k < NbrEigenstates; ++k)
// 			  {
// 			    RightMDaggerM[j][i] += (RightEigenstates[k][MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(j, PLevel, QValue) * TmpBMatrixDimension 
// 									+ MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(i, PLevel, QValue)] * TmpRightFactors[k]);
// 			  }
// 		      }
		  
		  
// 		  RealSymmetricMatrix SymLeftMDaggerM (LeftMDaggerM);
// 		  RealDiagonalMatrix TmpLeftDiag;
// #ifdef __LAPACK__
// 		  SymLeftMDaggerM.LapackDiagonalize(TmpLeftDiag);
// #else
// 		  SymLeftMDaggerM.Diagonalize(TmpLeftDiag,);
// #endif
// 		  for (int i = 0; i < TmpLeftDiag.GetNbrRow(); ++i)
// 		    {
// 		      if (fabs(TmpLeftDiag(i, i)) > LeftEigenvalueError)
// 			{
// 			  LeftEigenvalueError = fabs(TmpLeftDiag(i, i));	    
// 			}
// 		    }
// 		  RealSymmetricMatrix SymRightMDaggerM (RightMDaggerM);
// 		  RealDiagonalMatrix TmpRightDiag;
// #ifdef __LAPACK__
// 		  SymRightMDaggerM.LapackDiagonalize(TmpRightDiag);
// #else
// 		  SymRightMDaggerM.Diagonalize(TmpRightDiag);
// #endif
// 		  for (int i = 0; i < TmpRightDiag.GetNbrRow(); ++i)
// 		    {
// 		      if (fabs(TmpRightDiag(i, i)) > RightEigenvalueError)
// 			{
// 			  RightEigenvalueError = fabs(TmpRightDiag(i, i));	    
// 			}
// 		    }
// 		}
// 	    }
// 	}

      LeftEigenvalueError = Error;
      RightEigenvalueError = Error;

      double TotalTraceThoA = 0;

      double*** EntanglementSpectrum = new double**[MaxQValue - MinQValue + 1];
      int** EntanglementSpectrumDimension = new int*[MaxQValue - MinQValue + 1];
      
      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
	{
	  EntanglementSpectrum[QValue - MinQValue]  = new double* [Manager.GetInteger("p-truncation") + 1];
	  EntanglementSpectrumDimension[QValue - MinQValue]  = new int [Manager.GetInteger("p-truncation") + 1];
	  for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	    {
	      int IndexRange = MPSMatrix->GetBondIndexRange(PLevel, QValue);
	      if (IndexRange >= 0)
		{
		  EntanglementSpectrumDimension[QValue - MinQValue][PLevel] = 0;
		  EntanglementSpectrum[QValue - MinQValue][PLevel] = new double[IndexRange];
		}	      
	    }
	}


      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
	{
	  for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	    {
	      int IndexRange = MPSMatrix->GetBondIndexRange(PLevel, QValue);
	      if (IndexRange >= 0)
		{
		  ComplexMatrix LeftMDaggerM (IndexRange, IndexRange);		    
		  for (int i = 0; i < IndexRange; ++i)
		    for (int j = 0; j < IndexRange; ++j)
		      {
			LeftMDaggerM[j][i] = 0.0;
			for (int k = 0; k < NbrEigenstates; ++k)
			  {
			    int TmpIndex = SearchInArray<long>((((long) MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(j, PLevel, QValue)) * TmpBMatrixDimension) 
							       + MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(i, PLevel, QValue), 
							       EffectiveBlockIndices, EffectiveDimension);
			    LeftMDaggerM[j][i] += (LeftEigenstates[k][TmpIndex] * TmpLeftFactors[k]);
			  }
		      }
		  
		  
		  ComplexMatrix RightMDaggerM (IndexRange, IndexRange);
		  for (int i = 0; i < IndexRange; ++i)
		    for (int j = 0; j < IndexRange; ++j)
		      {
			RightMDaggerM[j][i] = 0.0;
			for (int k = 0; k < NbrEigenstates; ++k)
			  {
			    int TmpIndex = SearchInArray<long>((((long) MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(j, PLevel, QValue)) * TmpBMatrixDimension) 
							       + MPSMatrix->GetBondIndexWithFixedChargeAndPLevel(i, PLevel, QValue), 
							       EffectiveBlockIndices, EffectiveDimension);
			    RightMDaggerM[j][i] += (RightEigenstates[k][TmpIndex] * TmpRightFactors[k]);
			  }
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
		      if (fabs(TmpLeftDiag(i, i)) < LeftEigenvalueError)
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
		      if (fabs(TmpRightDiag(i, i)) < RightEigenvalueError)
			{
			  ++NbrZeroRightEigenvalues;	    
			}
		    }
		 
		  if ((NbrZeroLeftEigenvalues < SymLeftMDaggerM.GetNbrRow()) && (NbrZeroRightEigenvalues < SymRightMDaggerM.GetNbrRow()))
		    {
		      RealMatrix TruncatedLeftBasis (TmpLeftDiag.GetNbrRow(), TmpLeftDiag.GetNbrRow() -  NbrZeroLeftEigenvalues, true);
		      NbrZeroLeftEigenvalues = 0;
		      for (int i = 0; i < TmpLeftBasis.GetNbrColumn(); ++i)
			{
			  if (fabs(TmpLeftDiag(i, i)) > LeftEigenvalueError)
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
			  if (fabs(TmpRightDiag(i, i)) > RightEigenvalueError)
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
		      
		      RealSymmetricMatrix ReducedDensityMatrix ((Matrix&) TranposedTruncatedLeftBasis);

		      RealDiagonalMatrix TmpRhoADiag;
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
		      int NbrNonZeroEigenvalues = 0;
		      double Sum = 0.0;
		      for (int i = 0; i < TmpRhoADiag.GetNbrColumn(); ++i)
			{
			  if (TmpRhoADiag[i] > 0.0)
			    {
			      EntanglementSpectrum[QValue - MinQValue][PLevel][EntanglementSpectrumDimension[QValue - MinQValue][PLevel]] = TmpRhoADiag[i];
			      ++EntanglementSpectrumDimension[QValue - MinQValue][PLevel];
			      Sum += TmpRhoADiag(i, i);
			      ++NbrNonZeroEigenvalues;
			    }
			}
		      cout << "P=" << PLevel << " " << " Q=" << QValue << " NbrStates=" << NbrNonZeroEigenvalues << " Tr(rho_A)=" << Sum << endl;
		      TotalTraceThoA += Sum;
		    }
		}
	    }
	}

      for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
	{
	  for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	    {
	      for (int i = 0; i < EntanglementSpectrumDimension[QValue - MinQValue][PLevel]; ++i)
		File << "0 " << QValue << " " << PLevel << " " << PLevel << " " 
		     <<  (EntanglementSpectrum[QValue - MinQValue][PLevel][i] / TotalTraceThoA)  
		     <<  " " << (-log(EntanglementSpectrum[QValue - MinQValue][PLevel][i] / TotalTraceThoA)) <<endl;
	    }
	}

      cout << "Tr(rho_A)=" << TotalTraceThoA << endl;
      delete[] TmpLeftFactors;
      delete[] TmpRightFactors;
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
  int* RhoQSector = new int [MatDim];
  int* RhoPSector = new int [MatDim];

  for (int i =0 ; i < MatDim; i++)
   {
     RhoEigenvalues[i] = 0.0; 
     RhoQSector[i] = 0; 
     RhoPSector[i] = 0;
   }

  long NbrEigenvalues = 0l;
  OverlapMatrix.PrintNonZero(cout) << endl;
  int MinQValue = 0;
  int MaxQValue = 0;
  MPSMatrix->GetChargeIndexRange(MinQValue, MaxQValue);
  cout << "MinQValue = " << MinQValue << " MaxQValue= " << MaxQValue << endl;
  for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
    {
      for (int PLevel = 0; PLevel <= Manager.GetInteger("p-truncation"); ++PLevel)
	{
	  SparseRealMatrix TmpOverlapBlock = MPSMatrix->ExtractBlock(OverlapMatrix, PLevel, QValue, PLevel, QValue);
	  SparseRealMatrix RhoABlock = MPSMatrix->ExtractBlock(RhoA, PLevel, QValue, PLevel, QValue);

	  if ((TmpOverlapBlock.ComputeNbrNonZeroMatrixElements() != 0) && (RhoABlock.ComputeNbrNonZeroMatrixElements()))
	    {
 	      //cout << "QValue=" << QValue << "  PLevel=" << PLevel << " : "<< endl;
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
                      if (TmpDiag[j] < 0)
                        cout << "******* Negative j = " << j << " " << TmpDiag[j] << " ******* "<< endl;
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
		  
		  cout<<"------------sector P = "<<PLevel<<" Q = "<< QValue << "---------------"<<endl;
		  
		  
		  double Sum = 0.0;
		  for (int j = 0; j < NbrNonZeroVectors; ++j)
		    {
		      TraceRho += TmpDiagRho[j];
		      RhoEigenvalues[NbrEigenvalues] = TmpDiagRho[j];
		      RhoQSector[NbrEigenvalues] = QValue;
		      RhoPSector[NbrEigenvalues] = PLevel; 
		      NbrEigenvalues++;
		    }
		}
	    }
	}
    }
  
  cout<<"Trace rho = "<<TraceRho<<endl;


  int p = 0;
  int q = 0;
  if (Manager.GetBoolean("k-2") == true) 
    {
       p = 2;
       q = 2 + Manager.GetInteger("r-index");
    }
  else
   {
    if (Manager.GetBoolean("rr-3") == true)
      {
        p = 3;
        q = 5;
      }
     else
      {
        p = 1;
        q = Manager.GetInteger("laughlin-index");
      }  
   } 
  cout << "Filling factor nu = p/q = "<<p<<"/"<<q<<endl;
  int gcd = FindGCD(p, q);
  if (gcd > 1)
    {
      p /= gcd;
      q /= gcd;
      cout << "Filling factor nu* = p/q = "<<p<<"/"<<q<<endl;
    }

  int TmpNa;
  File << "# l_a    Na    Lz    lambda" << endl;
  for (int i = 0; i < MatDim; ++i)
   { 
    if (((fabs(RhoEigenvalues[i]) > CutOff))) 
      {
        TmpNa = RhoQSector[i];
        if ((p == 2) && (q == 4)) //Moore-Read
          TmpNa -= (MaxQValue - 1)/2;
        else
          TmpNa -= MaxQValue/2;

        if (Manager.GetBoolean("use-padding") == false)
          TmpNa -= p;

        TmpNa =  (p * EntCut - TmpNa)/q; 

        if (TmpNa == Na)
          {
            cout<< "Na= " << TmpNa << " P= " << RhoPSector[i] << " " << RhoEigenvalues[i]/TraceRho << endl;  
            File << EntCut << " " << TmpNa << " " << RhoPSector[i] << " " << (RhoEigenvalues[i] / TraceRho) << endl;
          }
        //else
        // cout << "Na = " << TmpNa << " " << RhoEigenvalues[i]/TraceRho<<endl;
      }
   }
  delete[] TmpMatrixElements;
  delete[] TmpColumnIndices;
  delete[] RhoEigenvalues;
  delete[] RhoPSector;
  delete[] RhoQSector;

  File.close();
 
  return 0;
}


void MPSDiagonalizeEMatrix(OptionManager* manager, AbstractHamiltonian* hamiltonian, int nbrEigenstates, Complex*& eigenvalues, ComplexVector*& eigenstates, 
			   AbstractArchitecture* architecture, double error, bool leftFlag)
{
  eigenstates = new ComplexVector[nbrEigenstates];
  eigenvalues = new Complex[nbrEigenstates];
  double* TmpRealPart = new double [nbrEigenstates];
  int* TmpIndices = new int[nbrEigenstates];
  Complex* TmpEigenvalues = 0;
  ComplexVector* TmpEigenstates = 0;

  if (hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension() < manager->GetInteger("full-diag"))
    {	  
      RealMatrix HRepresentation (hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension(), 
				  hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension());
      hamiltonian->GetHamiltonian(HRepresentation);
//      cout << HRepresentation << endl;
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
      BasicArnoldiAlgorithm* Arnoldi = 0;
      if (manager->GetBoolean("disk"))
	{
	  long TmpMemory = (((long) manager->GetInteger("arnoldi-memory")) << 17) / hamiltonian->GetHilbertSpace()->GetHilbertSpaceDimension();
	  if (TmpMemory == 0)
	    TmpMemory = 1;
	  Arnoldi = new BasicArnoldiAlgorithmWithDiskStorage (architecture, nbrEigenstates, 3000, true, false, manager->GetBoolean("resume"), TmpMemory, false);
	}
      else
	{
	  Arnoldi = new BasicArnoldiAlgorithm (architecture, nbrEigenstates, 3000, true, false, false);
	}
      Arnoldi->SetHamiltonian(hamiltonian);
      Arnoldi->InitializeLanczosAlgorithm();
      Arnoldi->RunLanczosAlgorithm(nbrEigenstates);
      bool ShowTimeFlag = manager->GetBoolean("show-itertime");
      while (Arnoldi->TestConvergence() == false)
	{
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalStartingTime), 0);
	    }
	  Arnoldi->RunLanczosAlgorithm(1);
	  if (ShowTimeFlag == true)
	    {
	      gettimeofday (&(TotalEndingTime), 0);
	      double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				    ((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	      cout << "iteration done in " << Dt << "s" << endl;
	    }
	}
      TmpEigenstates = (ComplexVector*) Arnoldi->GetEigenstates(nbrEigenstates);  
      Arnoldi->GetEigenvalues(TmpEigenvalues, nbrEigenstates);
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
