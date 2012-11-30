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
  if (Manager.GetBoolean("use-padding") == true)
    {
      if (Manager.GetBoolean("k-2") == true)
	{
	  if ((Manager.GetInteger("r-index") & 1) == 0)
	    MPSRowIndex = Manager.GetInteger("p-truncation") + (Manager.GetInteger("r-index") / 2);
	  else
	    MPSRowIndex = 2 * Manager.GetInteger("p-truncation") + Manager.GetInteger("r-index") - 1;
	}
      else
	{
	  if (Manager.GetBoolean("rr-3") == true)
	    {
	      MPSRowIndex = 3 * (Manager.GetInteger("p-truncation") + 1);
	    }
	  else
	    {
	      MPSRowIndex = Manager.GetInteger("p-truncation") + ((Manager.GetInteger("laughlin-index") - 1) / 2);
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
	}
      else
	{
	  if (Manager.GetBoolean("rr-3") == true)
	    {
	      MPSRowIndex = 3 * (Manager.GetInteger("p-truncation") + 2);
	      MPSColumnIndex = 3 * Manager.GetInteger("p-truncation");
	    }
	  else
	    {
	      MPSRowIndex = Manager.GetInteger("p-truncation") + (Manager.GetInteger("laughlin-index") - 1);
	      MPSColumnIndex = Manager.GetInteger("p-truncation");
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
	  sprintf(TmpFileName, "fermions_cylinder_%s_plevel_%ld_n_%d_2s_%d_lz_%d.0.full.ent", StateName,
		  Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz);
	}
      else
	{
	  sprintf(TmpFileName, "fermions_%s_plevel_%ld_n_%d_2s_%d_lz_%d.0.full.ent", StateName,
		  Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz);
	}
      File.open(TmpFileName, ios::binary | ios::out);     
   }


  if (Manager.GetBoolean("infinite-cylinder"))
    {
      int NbrBMatrices = 2;
      SparseRealMatrix** SparseTensorProductBMatrices = new SparseRealMatrix*[NbrBMatrices];
      for (int i = 0; i < NbrBMatrices; ++i)
	{
	  SparseTensorProductBMatrices[i] = new SparseRealMatrix[NbrBMatrices];
	  for (int j = 0; j < NbrBMatrices; ++j)
	    {
	      SparseTensorProductBMatrices[i][j] = TensorProduct(BMatrices[i], BMatrices[j]);
	    }
	}
      SparseRealMatrix NormalizedB0B0B1B1 = SparseRealMatrixLinearCombination(1.0, SparseTensorProductBMatrices[0][0], 1.0, SparseTensorProductBMatrices[1][1]);      
      RealMatrix NormalizedB0B0B1B1Full (NormalizedB0B0B1B1);

      ComplexDiagonalMatrix TmpLeftDiag (NormalizedB0B0B1B1Full.GetNbrRow(), true);  
      ComplexMatrix LeftEigenstates (NormalizedB0B0B1B1Full.GetNbrRow(), NormalizedB0B0B1B1Full.GetNbrRow(), true);  
      NormalizedB0B0B1B1Full.LapackDiagonalize(TmpLeftDiag, LeftEigenstates, true);
      int LargestLeftEigenvalueIndex = -1;
      int LargestLeftEigenvalueIndex2 = -1;
      int LargestLeftEigenvalueIndex3 = -1;
      int NbrNonZeroEigenvalues = 0;
      double Error = 1e-10;
      for (int i = 0; i < TmpLeftDiag.GetNbrRow(); ++i)
	if (Norm(TmpLeftDiag[i]) > Error)
	  ++NbrNonZeroEigenvalues;
      int* NonZeroEigenvalueIndices = new int [NbrNonZeroEigenvalues];
      double* NonZeroEigenvalueNorms = new double [NbrNonZeroEigenvalues];
      NbrNonZeroEigenvalues = 0;
      for (int i = 0; i < TmpLeftDiag.GetNbrRow(); ++i)
	if (Norm(TmpLeftDiag[i]) > Error)
	  {
	    NonZeroEigenvalueIndices[NbrNonZeroEigenvalues] = i;
	    NonZeroEigenvalueNorms[NbrNonZeroEigenvalues] = Norm(TmpLeftDiag[i]);
	    ++NbrNonZeroEigenvalues;
	  }
      SortArrayDownOrdering<int>(NonZeroEigenvalueNorms, NonZeroEigenvalueIndices, NbrNonZeroEigenvalues);
//       for (int i = 0; i < NbrNonZeroEigenvalues; ++i)
// 	cout << NonZeroEigenvalueIndices[i] << " " << NonZeroEigenvalueNorms[i] << " " << TmpLeftDiag[NonZeroEigenvalueIndices[i]] << endl;
      cout << "selected eigenvalues " <<  NonZeroEigenvalueIndices[0] << " " << NonZeroEigenvalueIndices[1] << " " << NonZeroEigenvalueIndices[2] << endl;
      int RealLeftEigenvalueIndex = NonZeroEigenvalueIndices[0];
      int ComplexLeftEigenvalueIndex = NonZeroEigenvalueIndices[1];
      int ComplexConjLeftEigenvalueIndex = NonZeroEigenvalueIndices[2];
      if (fabs(TmpLeftDiag[NonZeroEigenvalueIndices[0]].Im) < Error)
	{
	  RealLeftEigenvalueIndex = NonZeroEigenvalueIndices[0];
	  ComplexLeftEigenvalueIndex  = NonZeroEigenvalueIndices[1];
	  ComplexConjLeftEigenvalueIndex  = NonZeroEigenvalueIndices[2];
	}
      else
	{
	  if (fabs(TmpLeftDiag[NonZeroEigenvalueIndices[1]].Im) < Error)
	    {
	      RealLeftEigenvalueIndex = NonZeroEigenvalueIndices[1];
	      ComplexLeftEigenvalueIndex  = NonZeroEigenvalueIndices[0];
	      ComplexConjLeftEigenvalueIndex  = NonZeroEigenvalueIndices[2];
	    }
	  else
	    {
	      if (fabs(TmpLeftDiag[NonZeroEigenvalueIndices[2]].Im) < Error)
		{
		  RealLeftEigenvalueIndex = NonZeroEigenvalueIndices[2];
		  ComplexLeftEigenvalueIndex  = NonZeroEigenvalueIndices[0];
		  ComplexConjLeftEigenvalueIndex  = NonZeroEigenvalueIndices[1];
		}
	    }	  
	}
      cout << "real eigenvalue = " << RealLeftEigenvalueIndex << endl;
      int TmpBMatrixDimension = BMatrices[0].GetNbrRow();
      ComplexMatrix MDaggerM (TmpBMatrixDimension, TmpBMatrixDimension);
      for (int i = 0; i < TmpBMatrixDimension; ++i)
	for (int j = 0; j < TmpBMatrixDimension; ++j)
	  MDaggerM[j][i] = (LeftEigenstates[RealLeftEigenvalueIndex][j * TmpBMatrixDimension + i] + LeftEigenstates[ComplexLeftEigenvalueIndex][j * TmpBMatrixDimension + i] + LeftEigenstates[ComplexConjLeftEigenvalueIndex][j * TmpBMatrixDimension + i]) / sqrt(3.0);
      if (MDaggerM.IsReal(Error) == true)
	cout << "M^+M is real" << endl;
      if (MDaggerM.IsSymmetric(Error) == true)
	cout << "M^+M is symmetric" << endl;
      if (MDaggerM.IsHermitian(Error) == true)
	cout << "M^+M is hermitian" << endl;

      RealSymmetricMatrix SymMDaggerM (MDaggerM);
      RealMatrix TmpBasis(SymMDaggerM.GetNbrRow(), SymMDaggerM.GetNbrRow());
      TmpBasis.SetToIdentity();
      RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
      SymMDaggerM.LapackDiagonalize(TmpDiag, TmpBasis);
#else
      SymMDaggerM.Diagonalize(TmpDiag, TmpBasis);
#endif
      int NbrZeroEigenvalues = 0;
      for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
	if (fabs(TmpDiag(i, i)) < Error)
	  ++NbrZeroEigenvalues;
      cout << "nbr of zero eigenvalues = " << NbrZeroEigenvalues <<  " / " << TmpDiag.GetNbrRow() << endl;

      RealMatrix LeftM = RealMatrix (TmpDiag.GetNbrRow(),TmpDiag.GetNbrRow() -  NbrZeroEigenvalues);
      RealMatrix RightMInv = RealMatrix (TmpDiag.GetNbrRow(),TmpDiag.GetNbrRow() -  NbrZeroEigenvalues);
      int Count = 0;
      for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
	if (fabs(TmpDiag(i, i)) > Error)
	  {
	    LeftM[Count].Copy(TmpBasis[i]);
	    RightMInv[Count].Copy(TmpBasis[i]);
	    if (TmpDiag(i, i) > 0)
	      {
		LeftM[Count] *=  sqrt(TmpDiag(i, i));
		RightMInv[Count] /= sqrt(TmpDiag(i, i));		
	      }
	    else
	      {
		LeftM[Count] *=  sqrt(-TmpDiag(i, i));
		RightMInv[Count] /= sqrt(-TmpDiag(i, i));		
	      }
	    ++Count;
	  }
      LeftM.Transpose();
//      cout << RightMInv << endl;
//       RealMatrix TmpMatrix = LeftM * RightMInv;
//       cout << TmpMatrix << endl;

      SparseRealMatrix SparseLeftM(LeftM, Error);
      SparseRealMatrix SparseRightMInv(RightMInv, Error);
//       SparseLeftM.PrintNonZero(cout) << endl;
//       SparseRightMInv.PrintNonZero(cout) << endl;
      SparseRealMatrix* BTildaMatrices = new SparseRealMatrix[NbrBMatrices];
      for (int i = 0; i < NbrBMatrices; ++i)
	BTildaMatrices[i] = Conjugate(SparseLeftM, BMatrices[i], SparseRightMInv);
      
      SparseRealMatrix** SparseTensorProductBTildaMatrices = new SparseRealMatrix*[NbrBMatrices];
      for (int i = 0; i < NbrBMatrices; ++i)
	{
	  SparseTensorProductBTildaMatrices[i] = new SparseRealMatrix[NbrBMatrices];
	  for (int j = 0; j < NbrBMatrices; ++j)
	    {
	      SparseTensorProductBTildaMatrices[i][j] = TensorProduct(BTildaMatrices[i], BTildaMatrices[j]);
	    }
	}
      SparseRealMatrix NormalizedB0B0B1B12 = SparseRealMatrixLinearCombination(1.0, SparseTensorProductBTildaMatrices[0][0], 1.0, SparseTensorProductBTildaMatrices[1][1]);      
      RealMatrix NormalizedB0B0B1B1Full2 (NormalizedB0B0B1B12);
      LeftEigenstates.SetToIdentity();
      NormalizedB0B0B1B1Full2.LapackDiagonalize(TmpLeftDiag, LeftEigenstates, true);
      LargestLeftEigenvalueIndex = -1;
      LargestLeftEigenvalueIndex2 = -1;
      LargestLeftEigenvalueIndex3 = -1;
      NbrNonZeroEigenvalues = 0;
      for (int i = 0; i < TmpLeftDiag.GetNbrRow(); ++i)
	if (Norm(TmpLeftDiag[i]) > Error)
	  ++NbrNonZeroEigenvalues;
      int* NonZeroEigenvalueIndices2 = new int [NbrNonZeroEigenvalues];
      double* NonZeroEigenvalueNorms2 = new double [NbrNonZeroEigenvalues];
      NbrNonZeroEigenvalues = 0;
      for (int i = 0; i < TmpLeftDiag.GetNbrRow(); ++i)
	if (Norm(TmpLeftDiag[i]) > Error)
	  {
	    NonZeroEigenvalueIndices2[NbrNonZeroEigenvalues] = i;
	    NonZeroEigenvalueNorms2[NbrNonZeroEigenvalues] = Norm(TmpLeftDiag[i]);
	    ++NbrNonZeroEigenvalues;
	  }
      SortArrayDownOrdering<int>(NonZeroEigenvalueNorms2, NonZeroEigenvalueIndices2, NbrNonZeroEigenvalues);
//       for (int i = 0; i < NbrNonZeroEigenvalues; ++i)
// 	cout << NonZeroEigenvalueIndices2[i] << " " << NonZeroEigenvalueNorms2[i] << " " << TmpLeftDiag[NonZeroEigenvalueIndices2[i]] << endl;
      cout << "selected eigenvalues " <<  NonZeroEigenvalueIndices2[0] << " " << NonZeroEigenvalueIndices2[1] << " " << NonZeroEigenvalueIndices2[2] << endl;
      RealLeftEigenvalueIndex = NonZeroEigenvalueIndices2[0];
      ComplexLeftEigenvalueIndex = NonZeroEigenvalueIndices2[1];
      ComplexConjLeftEigenvalueIndex = NonZeroEigenvalueIndices2[2];
      if (fabs(TmpLeftDiag[NonZeroEigenvalueIndices2[0]].Im) < Error)
	{
	  RealLeftEigenvalueIndex = NonZeroEigenvalueIndices2[0];
	  ComplexLeftEigenvalueIndex  = NonZeroEigenvalueIndices2[1];
	  ComplexConjLeftEigenvalueIndex  = NonZeroEigenvalueIndices2[2];
	}
      else
	{
	  if (fabs(TmpLeftDiag[NonZeroEigenvalueIndices2[1]].Im) < Error)
	    {
	      RealLeftEigenvalueIndex = NonZeroEigenvalueIndices2[1];
	      ComplexLeftEigenvalueIndex  = NonZeroEigenvalueIndices2[0];
	      ComplexConjLeftEigenvalueIndex  = NonZeroEigenvalueIndices2[2];
	    }
	  else
	    {
	      if (fabs(TmpLeftDiag[NonZeroEigenvalueIndices2[2]].Im) < Error)
		{
		  RealLeftEigenvalueIndex = NonZeroEigenvalueIndices2[2];
		  ComplexLeftEigenvalueIndex  = NonZeroEigenvalueIndices2[0];
		  ComplexConjLeftEigenvalueIndex  = NonZeroEigenvalueIndices2[1];
		}
	    }	  
	}
      cout << "real eigenvalue = " << RealLeftEigenvalueIndex << endl;

      TmpBMatrixDimension = BTildaMatrices[0].GetNbrRow();
      ComplexMatrix RhoA (TmpBMatrixDimension, TmpBMatrixDimension);
      for (int i = 0; i < TmpBMatrixDimension; ++i)
	for (int j = 0; j < TmpBMatrixDimension; ++j)
	  RhoA[j][i] = (LeftEigenstates[RealLeftEigenvalueIndex][j * TmpBMatrixDimension + i] + LeftEigenstates[ComplexLeftEigenvalueIndex][j * TmpBMatrixDimension + i] + LeftEigenstates[ComplexConjLeftEigenvalueIndex][j * TmpBMatrixDimension + i]) / sqrt(3.0);
      if (RhoA.IsReal(Error) == true)
	cout << "rho_A is real" << endl;
      if (RhoA.IsSymmetric(Error) == true)
	cout << "rho_A is symmetric" << endl;
      if (RhoA.IsHermitian(Error) == true)
	cout << "rho_A is hermitian" << endl;
      if (RhoA.IsDiagonal(Error) == true)
	{
	  cout << "rho_A is diagonal" << endl;
	  int Count = 0;
	  for (int i = 0; i < TmpBMatrixDimension; ++i)
	    {
	      double Tmp;
	      RhoA.GetMatrixElement(i, i, Tmp);
	      cout << i << " " << Tmp << endl;
	      if (fabs(Tmp) > Error)
		{
		  cout << Count << " : " << Tmp << endl;
		  ++Count;
		}
	    }
	}

//      cout << RhoA << endl;

//       cout << TmpMatrix << endl;
//       for (int i = 0; i < TmpBMatrixDimension; ++i)
// 	for (int j = 0; j < TmpBMatrixDimension; ++j)
// 	  if (Norm(MDaggerM[i][j]) > Error)
// 	    cout << i << " " << j << " = " << MDaggerM[i][j] << endl;
// 	  cout << endl << endl;

//       for (int i = 0; i < TmpBMatrixDimension; ++i)
// 	for (int j = 0; j < TmpBMatrixDimension; ++j)
// 	  MDaggerM[j][i] = (Phase(2.0 *M_PI/3.0) * LeftEigenstates[LargestLeftEigenvalueIndex][j * TmpBMatrixDimension + i] + Phase(4.0 *M_PI/3.0) * LeftEigenstates[LargestLeftEigenvalueIndex2][j * TmpBMatrixDimension + i] + LeftEigenstates[LargestLeftEigenvalueIndex3][j * TmpBMatrixDimension + i]) / sqrt(3.0);
//       if (MDaggerM.IsDiagonal(1e-12) == true)
// 	{
// 	  for (int i = 0; i < TmpBMatrixDimension; ++i)
// 	    cout << MDaggerM[i][i] << " ";
// 	  cout << endl << endl;
// 	}
//       for (int i = 0; i < TmpBMatrixDimension; ++i)
// 	for (int j = 0; j < TmpBMatrixDimension; ++j)
// 	  MDaggerM[j][i] = (LeftEigenstates[LargestLeftEigenvalueIndex][j * TmpBMatrixDimension + i] + Phase(2.0 *M_PI/3.0) * LeftEigenstates[LargestLeftEigenvalueIndex2][j * TmpBMatrixDimension + i] + Phase(4.0 *M_PI/3.0) * LeftEigenstates[LargestLeftEigenvalueIndex3][j * TmpBMatrixDimension + i]) / sqrt(3.0);
// //       if (MDaggerM.IsDiagonal(1e-10) == true)
// // 	{
//       if (MDaggerM.IsReal(1e-10) == true)
// 	cout << "M^+M is real" << endl;
//       if (MDaggerM.IsSymmetric(1e-10) == true)
// 	cout << "M^+M is symmetric" << endl;
//       for (int i = 0; i < TmpBMatrixDimension; ++i)
// 	for (int j = 0; j < TmpBMatrixDimension; ++j)
// 	  if (Norm(MDaggerM[i][j]) > 1e-10)
// 	    cout << i << " " << j << " = " << MDaggerM[i][j] << endl;
// 	  cout << endl << endl;
//	}
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
        cout<< "P= " << RhoPSector[i] << " N= " << RhoNSector[i] << " " << RhoEigenvalues[i]/TraceRho << endl;  
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
