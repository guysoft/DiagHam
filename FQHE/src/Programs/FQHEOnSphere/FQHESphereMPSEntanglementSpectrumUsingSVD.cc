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


/*
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
      SparseRealMatrix NormalizedB0B0B1B1V2 = Multiply (NormalizedB0B0B1B1, NormalizedB0B0B1B1);
      SparseRealMatrix NormalizedB0B0B1B1V3 = Multiply (NormalizedB0B0B1B1V2, NormalizedB0B0B1B1);
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
// 	cout << i << " : " << TmpLeftDiag[NonZeroEigenvalueIndices[i]] << " " << TmpLeftDiag[NonZeroEigenvalueIndices[i]] << " " 
// 	     << NonZeroEigenvalueNorms[i] << endl;
   
      cout << "selected left eigenvalues " <<  TmpLeftDiag[NonZeroEigenvalueIndices[0]] << " (" << NonZeroEigenvalueIndices[0] << "), " 
	   << " " << TmpLeftDiag[NonZeroEigenvalueIndices[1]] << " (" << NonZeroEigenvalueIndices[1] << "), " 
	   << " " << TmpLeftDiag[NonZeroEigenvalueIndices[2]] << " (" << NonZeroEigenvalueIndices[2] << ")" << endl;
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
      cout << "real left eigenvalue = " << RealLeftEigenvalueIndex << endl;

      ComplexDiagonalMatrix TmpRightDiag (NormalizedB0B0B1B1Full.GetNbrRow(), true);  
      ComplexMatrix RightEigenstates (NormalizedB0B0B1B1Full.GetNbrRow(), NormalizedB0B0B1B1Full.GetNbrRow(), true);  
      NormalizedB0B0B1B1Full.LapackDiagonalize(TmpRightDiag, RightEigenstates, false);
      int LargestRightEigenvalueIndex = -1;
      int LargestRightEigenvalueIndex2 = -1;
      int LargestRightEigenvalueIndex3 = -1;
      int NbrNonZeroRightEigenvalues = 0;
      for (int i = 0; i < TmpRightDiag.GetNbrRow(); ++i)
	{
	  if (Norm(TmpRightDiag[i]) > Error)
	    ++NbrNonZeroRightEigenvalues;
	}
      int* NonZeroRightEigenvalueIndices = new int [NbrNonZeroRightEigenvalues];
      double* NonZeroRightEigenvalueNorms = new double [NbrNonZeroRightEigenvalues];
      NbrNonZeroRightEigenvalues = 0;
      for (int i = 0; i < TmpRightDiag.GetNbrRow(); ++i)
	if (Norm(TmpRightDiag[i]) > Error)
	  {
	    NonZeroRightEigenvalueIndices[NbrNonZeroRightEigenvalues] = i;
	    NonZeroRightEigenvalueNorms[NbrNonZeroRightEigenvalues] = Norm(TmpRightDiag[i]);
	    ++NbrNonZeroRightEigenvalues;
	  }
      SortArrayDownOrdering<int>(NonZeroRightEigenvalueNorms, NonZeroRightEigenvalueIndices, NbrNonZeroRightEigenvalues);
      cout << "selected right eigenvalues " <<  NonZeroRightEigenvalueIndices[0] << " " << NonZeroRightEigenvalueIndices[1] << " " << NonZeroRightEigenvalueIndices[2] << endl;
      int RealRightEigenvalueIndex = NonZeroRightEigenvalueIndices[0];
      int ComplexRightEigenvalueIndex = NonZeroRightEigenvalueIndices[1];
      int ComplexConjRightEigenvalueIndex = NonZeroRightEigenvalueIndices[2];
      if (fabs(TmpRightDiag[NonZeroRightEigenvalueIndices[0]].Im) < Error)
	{
	  RealRightEigenvalueIndex = NonZeroRightEigenvalueIndices[0];
	  ComplexRightEigenvalueIndex  = NonZeroRightEigenvalueIndices[1];
	  ComplexConjRightEigenvalueIndex  = NonZeroRightEigenvalueIndices[2];
	}
      else
	{
	  if (fabs(TmpRightDiag[NonZeroRightEigenvalueIndices[1]].Im) < Error)
	    {
	      RealRightEigenvalueIndex = NonZeroRightEigenvalueIndices[1];
	      ComplexRightEigenvalueIndex  = NonZeroRightEigenvalueIndices[0];
	      ComplexConjRightEigenvalueIndex  = NonZeroRightEigenvalueIndices[2];
	    }
	  else
	    {
	      if (fabs(TmpRightDiag[NonZeroRightEigenvalueIndices[2]].Im) < Error)
		{
		  RealRightEigenvalueIndex = NonZeroRightEigenvalueIndices[2];
		  ComplexRightEigenvalueIndex  = NonZeroRightEigenvalueIndices[0];
		  ComplexConjRightEigenvalueIndex  = NonZeroRightEigenvalueIndices[1];
		}
	    }	  
	}
      if ((TmpRightDiag[ComplexConjRightEigenvalueIndex].Im * TmpLeftDiag[ComplexConjLeftEigenvalueIndex].Im) < 0.0)
	{
	  int Tmp = ComplexConjRightEigenvalueIndex;
	  ComplexConjRightEigenvalueIndex = ComplexRightEigenvalueIndex;
	  ComplexRightEigenvalueIndex = Tmp;
	}
      cout << "real right eigenvalue = " << RealRightEigenvalueIndex << endl;
      
      cout << (LeftEigenstates[RealLeftEigenvalueIndex] * RightEigenstates[RealRightEigenvalueIndex]) << endl;
      cout << (LeftEigenstates[ComplexLeftEigenvalueIndex] * RightEigenstates[ComplexRightEigenvalueIndex]) << endl;
      cout << (LeftEigenstates[ComplexConjLeftEigenvalueIndex] * RightEigenstates[ComplexRightEigenvalueIndex]) << endl;
      cout << (LeftEigenstates[ComplexLeftEigenvalueIndex] * RightEigenstates[ComplexConjRightEigenvalueIndex]) << endl;
      cout << (LeftEigenstates[ComplexConjLeftEigenvalueIndex] * RightEigenstates[ComplexConjRightEigenvalueIndex]) << endl;

      cout << (LeftEigenstates[RealLeftEigenvalueIndex] * LeftEigenstates[ComplexLeftEigenvalueIndex]) << endl;
      cout << (LeftEigenstates[RealLeftEigenvalueIndex] * LeftEigenstates[ComplexConjLeftEigenvalueIndex]) << endl;


      int TmpBMatrixDimension = BMatrices[0].GetNbrRow();
      ComplexMatrix LeftMDaggerM (TmpBMatrixDimension, TmpBMatrixDimension);
      cout << "test " << endl;
//      ++MPSRowIndex;
//      ++MPSColumnIndex;
//       ++MPSRowIndex;
//       ++MPSColumnIndex;
//       ++MPSRowIndex;
//       ++MPSColumnIndex;
//       cout << RightEigenstates[RealRightEigenvalueIndex][MPSRowIndex * TmpBMatrixDimension + MPSColumnIndex] << endl;
//       cout << RightEigenstates[ComplexRightEigenvalueIndex][MPSRowIndex * TmpBMatrixDimension + MPSColumnIndex] << endl;
//       cout << RightEigenstates[ComplexConjRightEigenvalueIndex][MPSRowIndex * TmpBMatrixDimension + MPSColumnIndex] << endl;

      cout << (RightEigenstates[RealRightEigenvalueIndex][MPSRowIndex * TmpBMatrixDimension + MPSColumnIndex] 
	        /  (LeftEigenstates[RealRightEigenvalueIndex] * RightEigenstates[RealRightEigenvalueIndex])) << endl;
      cout << (RightEigenstates[ComplexRightEigenvalueIndex][MPSRowIndex * TmpBMatrixDimension + MPSColumnIndex] 
	       /  (LeftEigenstates[ComplexLeftEigenvalueIndex] * RightEigenstates[ComplexRightEigenvalueIndex])) << endl;
      cout << (RightEigenstates[ComplexConjRightEigenvalueIndex][MPSRowIndex * TmpBMatrixDimension + MPSColumnIndex] 
	       /  (LeftEigenstates[ComplexConjLeftEigenvalueIndex] * RightEigenstates[ComplexConjRightEigenvalueIndex])) << endl;
      for (int i = 0; i < TmpBMatrixDimension; ++i)
	for (int j = 0; j < TmpBMatrixDimension; ++j)
	  {
 	    LeftMDaggerM[j][i] = (LeftEigenstates[RealLeftEigenvalueIndex][j * TmpBMatrixDimension + i] * RightEigenstates[RealRightEigenvalueIndex][MPSRowIndex * TmpBMatrixDimension + MPSColumnIndex] 
 				  /  (LeftEigenstates[RealRightEigenvalueIndex] * RightEigenstates[RealRightEigenvalueIndex]));
  	    LeftMDaggerM[j][i] += (LeftEigenstates[ComplexLeftEigenvalueIndex][j * TmpBMatrixDimension + i] * RightEigenstates[ComplexRightEigenvalueIndex][MPSRowIndex * TmpBMatrixDimension + MPSColumnIndex] 
  				  /  (LeftEigenstates[ComplexLeftEigenvalueIndex] * RightEigenstates[ComplexRightEigenvalueIndex]));
  	    LeftMDaggerM[j][i] += (LeftEigenstates[ComplexConjLeftEigenvalueIndex][j * TmpBMatrixDimension + i] * RightEigenstates[ComplexConjRightEigenvalueIndex][MPSRowIndex * TmpBMatrixDimension + MPSColumnIndex] 
  				  /  (LeftEigenstates[ComplexConjLeftEigenvalueIndex] * RightEigenstates[ComplexConjRightEigenvalueIndex]));
	  }

//       SparseRealMatrix SparseLeftMDaggerM (LeftMDaggerM);
//       for (int i = 0;i < 7; ++i)
// 	{
// 	  SparseRealMatrix TmpOverlapBlock = MPSMatrix->ExtractBlock(SparseLeftMDaggerM, 2, i, 2, 2);      
// 	  cout << TmpOverlapBlock << endl << endl;
// 	}

      if (LeftMDaggerM.IsReal(Error) == true)
	cout << "left M^+M is real" << endl;
      if (LeftMDaggerM.IsSymmetric(Error) == true)
	cout << "left M^+M is symmetric" << endl;
      if (LeftMDaggerM.IsHermitian(Error) == true)
	cout << "left M^+M is hermitian" << endl;

      ComplexMatrix RightMDaggerM (TmpBMatrixDimension, TmpBMatrixDimension);
      for (int i = 0; i < TmpBMatrixDimension; ++i)
	for (int j = 0; j < TmpBMatrixDimension; ++j)
	  {
	    RightMDaggerM[j][i] = (RightEigenstates[RealRightEigenvalueIndex][j * TmpBMatrixDimension + i] * LeftEigenstates[RealLeftEigenvalueIndex][MPSRowIndex * TmpBMatrixDimension + MPSColumnIndex] 
				  /  (LeftEigenstates[RealLeftEigenvalueIndex] * RightEigenstates[RealRightEigenvalueIndex]));
	    RightMDaggerM[j][i] += (RightEigenstates[ComplexRightEigenvalueIndex][j * TmpBMatrixDimension + i] * LeftEigenstates[ComplexLeftEigenvalueIndex][MPSRowIndex * TmpBMatrixDimension + MPSColumnIndex] 
				  /  (LeftEigenstates[ComplexLeftEigenvalueIndex] * RightEigenstates[ComplexRightEigenvalueIndex]));
	    RightMDaggerM[j][i] += (RightEigenstates[ComplexConjRightEigenvalueIndex][j * TmpBMatrixDimension + i] * LeftEigenstates[ComplexConjLeftEigenvalueIndex][MPSRowIndex * TmpBMatrixDimension + MPSColumnIndex] 
				  /  (LeftEigenstates[ComplexConjLeftEigenvalueIndex] * RightEigenstates[ComplexConjRightEigenvalueIndex]));
	  }


//	  RightMDaggerM[j][i] = RightEigenstates[RealRightEigenvalueIndex][j * TmpBMatrixDimension + i];// + LeftEigenstates[ComplexLeftEigenvalueIndex][j * TmpBMatrixDimension + i] + LeftEigenstates[ComplexConjLeftEigenvalueIndex][j * TmpBMatrixDimension + i]) / sqrt(3.0);
      if (RightMDaggerM.IsReal(Error) == true)
	cout << "right M^+M is real" << endl;
      if (RightMDaggerM.IsSymmetric(Error) == true)
	cout << "right M^+M is symmetric" << endl;
      if (RightMDaggerM.IsHermitian(Error) == true)
	cout << "right M^+M is hermitian" << endl;

      RealSymmetricMatrix SymMDaggerM (LeftMDaggerM);
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
	{
	  if (fabs(TmpDiag(i, i)) < Error)
	    {
	      ++NbrZeroEigenvalues;	    
	    }
	  else
	    {
	      cout << "ev="  << TmpDiag(i, i) << endl;
	    }
	}
      cout << "nbr of zero eigenvalues = " << NbrZeroEigenvalues <<  " / " << TmpDiag.GetNbrRow() << endl;

//       ComplexMatrix LeftM = ComplexMatrix (TmpDiag.GetNbrRow(),TmpDiag.GetNbrRow() -  NbrZeroEigenvalues, true);
//       ComplexMatrix RightMInv = ComplexMatrix (TmpDiag.GetNbrRow(),TmpDiag.GetNbrRow() -  NbrZeroEigenvalues, true);
//       int Count = 0;
//       for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
// 	if (fabs(TmpDiag(i, i)) > Error)
// 	  {
// 	    ComplexVector& TmpLeftVector = LeftM[Count];
// 	    ComplexVector& TmpRightVector = RightMInv[Count];
// 	    RealVector& TmpVector = TmpBasis[i];
// 	    if (TmpDiag(i, i) > 0)
// 	      {
// 		double Factor = sqrt(TmpDiag(i, i));
// 		double InvFactor = 1.0 / sqrt(TmpDiag(i, i));
// 		for (int j = 0; j < TmpDiag.GetNbrRow(); ++j)
// 		  {
// 		    TmpLeftVector[j].Re = TmpVector[j] * Factor;
// 		    TmpRightVector[j].Re = TmpVector[j] * InvFactor;
// 		  }
// 	      }
// 	    else
// 	      {
// 		double Factor = sqrt(-TmpDiag(i, i));
// 		double InvFactor = 1.0 / sqrt(-TmpDiag(i, i));
// 		for (int j = 0; j < TmpDiag.GetNbrRow(); ++j)
// 		  {
// 		    TmpLeftVector[j].Im = TmpVector[j] * Factor;
// 		    TmpRightVector[j].Im = TmpVector[j] * InvFactor;
// 		  }
// 	      }
// 	    ++Count;
// 	  }
//       LeftM.HermitianTranspose();
// //      cout << RightMInv << endl;
//       ComplexMatrix TmpMatrix = LeftM * RightMInv;
// //      cout << TmpMatrix << endl;

//       SparseComplexMatrix SparseLeftM(LeftM, Error);
//       SparseComplexMatrix SparseRightMInv(RightMInv, Error);
// //       SparseLeftM.PrintNonZero(cout) << endl;
// //       SparseRightMInv.PrintNonZero(cout) << endl;
//       SparseComplexMatrix* BTildaMatrices = new SparseComplexMatrix[NbrBMatrices];
//       for (int i = 0; i < NbrBMatrices; ++i)
// 	BTildaMatrices[i] = Conjugate(SparseLeftM, BMatrices[i], SparseRightMInv);
      
//       SparseComplexMatrix** SparseTensorProductBTildaMatrices = new SparseComplexMatrix*[NbrBMatrices];
//       for (int i = 0; i < NbrBMatrices; ++i)
// 	{
// 	  SparseTensorProductBTildaMatrices[i] = new SparseComplexMatrix[NbrBMatrices];
// 	  for (int j = 0; j < NbrBMatrices; ++j)
// 	    {
// 	      SparseTensorProductBTildaMatrices[i][j] = TensorProductWithConjugation(BTildaMatrices[i], BTildaMatrices[j]);
// 	    }
// 	}
//       SparseComplexMatrix NormalizedB0B0B1B12 = SparseComplexMatrixLinearCombination(1.0, SparseTensorProductBTildaMatrices[0][0], 1.0, SparseTensorProductBTildaMatrices[1][1]);      
//       ComplexMatrix NormalizedB0B0B1B1Full2 (NormalizedB0B0B1B12);
//       ComplexMatrix ComplexLeftEigenstates (NormalizedB0B0B1B12.GetNbrRow(), NormalizedB0B0B1B12.GetNbrRow());
//       ComplexLeftEigenstates.SetToIdentity();
//       NormalizedB0B0B1B1Full2.LapackDiagonalize(TmpLeftDiag, ComplexLeftEigenstates);//, true);
//       LargestLeftEigenvalueIndex = -1;
//       LargestLeftEigenvalueIndex2 = -1;
//       LargestLeftEigenvalueIndex3 = -1;
//       NbrNonZeroEigenvalues = 0;
//       for (int i = 0; i < TmpLeftDiag.GetNbrRow(); ++i)
// 	if (Norm(TmpLeftDiag[i]) > Error)
// 	  ++NbrNonZeroEigenvalues;
//       int* NonZeroEigenvalueIndices2 = new int [NbrNonZeroEigenvalues];
//       double* NonZeroEigenvalueNorms2 = new double [NbrNonZeroEigenvalues];
//       NbrNonZeroEigenvalues = 0;
//       for (int i = 0; i < TmpLeftDiag.GetNbrRow(); ++i)
// 	if (Norm(TmpLeftDiag[i]) > Error)
// 	  {
// 	    NonZeroEigenvalueIndices2[NbrNonZeroEigenvalues] = i;
// 	    NonZeroEigenvalueNorms2[NbrNonZeroEigenvalues] = Norm(TmpLeftDiag[i]);
// 	    ++NbrNonZeroEigenvalues;
// 	  }
//       SortArrayDownOrdering<int>(NonZeroEigenvalueNorms2, NonZeroEigenvalueIndices2, NbrNonZeroEigenvalues);
//       for (int i = 0; i < NbrNonZeroEigenvalues; ++i)
// 	cout << NonZeroEigenvalueIndices2[i] << " " << NonZeroEigenvalueNorms2[i] << " " << TmpLeftDiag[NonZeroEigenvalueIndices2[i]] << endl;
//       cout << "selected eigenvalues " <<  NonZeroEigenvalueIndices2[0] << " " << NonZeroEigenvalueIndices2[1] << " " << NonZeroEigenvalueIndices2[2] << endl;
//       RealLeftEigenvalueIndex = NonZeroEigenvalueIndices2[0];
//       ComplexLeftEigenvalueIndex = NonZeroEigenvalueIndices2[1];
//       ComplexConjLeftEigenvalueIndex = NonZeroEigenvalueIndices2[2];
//       if (fabs(TmpLeftDiag[NonZeroEigenvalueIndices2[0]].Im) < Error)
// 	{
// 	  RealLeftEigenvalueIndex = NonZeroEigenvalueIndices2[0];
// 	  ComplexLeftEigenvalueIndex  = NonZeroEigenvalueIndices2[1];
// 	  ComplexConjLeftEigenvalueIndex  = NonZeroEigenvalueIndices2[2];
// 	}
//       else
// 	{
// 	  if (fabs(TmpLeftDiag[NonZeroEigenvalueIndices2[1]].Im) < Error)
// 	    {
// 	      RealLeftEigenvalueIndex = NonZeroEigenvalueIndices2[1];
// 	      ComplexLeftEigenvalueIndex  = NonZeroEigenvalueIndices2[0];
// 	      ComplexConjLeftEigenvalueIndex  = NonZeroEigenvalueIndices2[2];
// 	    }
// 	  else
// 	    {
// 	      if (fabs(TmpLeftDiag[NonZeroEigenvalueIndices2[2]].Im) < Error)
// 		{
// 		  RealLeftEigenvalueIndex = NonZeroEigenvalueIndices2[2];
// 		  ComplexLeftEigenvalueIndex  = NonZeroEigenvalueIndices2[0];
// 		  ComplexConjLeftEigenvalueIndex  = NonZeroEigenvalueIndices2[1];
// 		}
// 	    }	  
// 	}
//       cout << "real eigenvalue = " << RealLeftEigenvalueIndex << endl;

//       TmpBMatrixDimension = BTildaMatrices[0].GetNbrRow();
//       ComplexMatrix RhoA (TmpBMatrixDimension, TmpBMatrixDimension);
//       for (int i = 0; i < TmpBMatrixDimension; ++i)
// 	for (int j = 0; j < TmpBMatrixDimension; ++j)
// 	  RhoA[j][i] = (ComplexLeftEigenstates[RealLeftEigenvalueIndex][j * TmpBMatrixDimension + i] + ComplexLeftEigenstates[ComplexLeftEigenvalueIndex][j * TmpBMatrixDimension + i] + ComplexLeftEigenstates[ComplexConjLeftEigenvalueIndex][j * TmpBMatrixDimension + i]) / sqrt(3.0);
//       if (RhoA.IsReal(Error) == true)
// 	cout << "rho_A is real" << endl;
//       if (RhoA.IsSymmetric(Error) == true)
// 	cout << "rho_A is symmetric" << endl;
//       if (RhoA.IsHermitian(Error) == true)
// 	cout << "rho_A is hermitian" << endl;
//       if (RhoA.IsDiagonal(Error) == true)
// 	{
// 	  cout << "rho_A is diagonal" << endl;
// 	  int Count = 0;
// 	  for (int i = 0; i < TmpBMatrixDimension; ++i)
// 	    {
// 	      double Tmp;
// 	      RhoA.GetMatrixElement(i, i, Tmp);
// 	      cout << i << " " << Tmp << endl;
// 	      if (fabs(Tmp) > Error)
// 		{
// 		  cout << Count << " : " << Tmp << endl;
// 		  ++Count;
// 		}
// 	    }
// 	}

// //      cout << RhoA << endl;
//       for (int i = 0; i < TmpBMatrixDimension; ++i)
// 	for (int j = 0; j < TmpBMatrixDimension; ++j)
// 	  if (Norm(RhoA[i][j]) > Error)
// 	    cout << i << " " << j << " = " << RhoA[i][j] << endl;
// 	  cout << endl << endl;

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
*/

  //cout<<"B0"<<endl;
  //cout<<BMatrices[0]<<endl;
  //cout<<"B1"<<endl;
  //cout<<BMatrices[1]<<endl;

  int MatDim = BMatrices[0].GetNbrRow();
  int LambdaMax = Manager.GetInteger("p-truncation");
  int LaughlinIndex = Manager.GetInteger("laughlin-index");

  cout << "B matrix size = " << MatDim << "x" << MatDim << endl;

  double CutOff = 1e-14;

  double* NormalizationCoefficients = new double[NbrFluxQuanta + 1];
  BinomialCoefficients Binomial(NbrFluxQuanta);
  for (int i = 0; i <= NbrFluxQuanta; ++i)
    {      
      NormalizationCoefficients[i] = ((double) (NbrFluxQuanta + 1)) / Binomial.GetNumericalCoefficient(NbrFluxQuanta, i);
    }


  RealMatrix* TmpBUD;
  RealMatrix* UMatrix;
  RealMatrix* VMatrix;
  RealMatrix* NewDMatrix;
  RealMatrix* NewDVMatrix;
  RealMatrix* NewUDMatrix;
  double* SingularValues; 
  int MaxDim;

  RealMatrix DenseB0 (BMatrices[0]);
  RealMatrix DenseB1 (BMatrices[1]);


  cout<<"----------Start SVD from right end-------------"<<endl;

  TmpBUD = new RealMatrix (MatDim, 2, true);
  double Tmp;
  for (int i = 0; i < MatDim; i++)
    {
      BMatrices[0].GetMatrixElement(i, MPSColumnIndex, Tmp);
      TmpBUD->SetMatrixElement(i, 0, Tmp);
      BMatrices[1].GetMatrixElement(i, MPSColumnIndex, Tmp);
      TmpBUD->SetMatrixElement(i, 1, Tmp);
    }
  //cout<<"TmpBUD"<<endl;
  //cout<<*TmpBUD<<endl;

/*
  RealMatrix TestSVD (10, 2, true);
  TestSVD.SetMatrixElement(0, 0, 0.76966541249324);
  TestSVD.SetMatrixElement(3, 1, 0.87730576909835);

  cout<<"TestSVD "<<TestSVD<<endl;
  RealMatrix* TestU = new RealMatrix(10,10,true);

  RealMatrix* TestV = new RealMatrix(2,2,true);
  double* SingTest = TestSVD.SingularValueDecomposition(*TestU, *TestV, false);

  RealMatrix* TestD = new RealMatrix (10, 2, true);
  for (int i = 0; i < 2; i++)
    TestD->SetMatrixElement(i,i,SingTest[i]);  

  cout<<"D "<<*TestD<<endl;

  RealMatrix TestUDV (10, 2, true);
  TestUDV.Copy(*TestU);
  
  TestUDV.Multiply(*TestD);
  TestUDV.Multiply(*TestV);
  cout<<"Test"<<endl;
  cout<<TestUDV<<endl;
  exit(1);
*/

  //cout<<"U"<<endl;
  //cout<<(*TestU)<<endl;
  //cout<<"V"<<endl;
  //cout<<(*TestV)<<endl;
  

  UMatrix = new RealMatrix (MatDim, MatDim, true);
  VMatrix = new RealMatrix (2, 2, true);

  SingularValues = TmpBUD->SingularValueDecomposition(*UMatrix, *VMatrix, false);

  int SingularDimension = 2;
  if (MatDim < 2)
    SingularDimension = MatDim;
 
  //cout<<"Singular values: ";
  //for (int i = 0; i < SingularDimension; i++)
  //  cout<<SingularValues[i]<<" ";
  //cout<<endl;

  int NbrNonZeroSingularValues = 0;
  for (int i = 0; i < SingularDimension; i++)
   if (SingularValues[i] > CutOff)
     NbrNonZeroSingularValues++;

  //cout<<"non zero = "<<NbrNonZeroSingularValues<<endl;

  //cout<<"UMatrix "<<endl;
  //cout<<*UMatrix<<endl;

  //cout<<"VMatrix "<<endl;
  //cout<<*VMatrix<<endl;

  UMatrix->Resize(MatDim, NbrNonZeroSingularValues);
  VMatrix->Resize(NbrNonZeroSingularValues, 2);
  
  NewDMatrix = new RealMatrix (NbrNonZeroSingularValues, NbrNonZeroSingularValues);
  for (int i = 0; i < NbrNonZeroSingularValues; ++i)
    NewDMatrix->SetMatrixElement(i, i, SingularValues[i]);

  NewUDMatrix = new RealMatrix (MatDim, NbrNonZeroSingularValues, true);
  NewUDMatrix->Copy(UMatrix->Multiply(*NewDMatrix));

  //cout<<"NewUD"<<*NewUDMatrix<<endl;

  delete UMatrix;
  delete VMatrix;
  delete NewDMatrix;
  delete TmpBUD;
  delete[] SingularValues;


  int OldNbrNonZeroSingularValues = NbrNonZeroSingularValues;

  for (int TmpCut = NbrFluxQuanta; TmpCut > EntCut; TmpCut--)
    {
       RealMatrix TmpB0 (MatDim, MatDim);
       TmpB0.Copy(DenseB0);
       TmpB0.Multiply(*NewUDMatrix);
       TmpB0.Resize(MatDim, OldNbrNonZeroSingularValues);

       RealMatrix TmpB1 (MatDim, MatDim);
       TmpB1.Copy(DenseB1);
       TmpB1.Multiply(*NewUDMatrix);
       TmpB1.Resize(MatDim, OldNbrNonZeroSingularValues);

       TmpBUD = new RealMatrix (MatDim, 2 * OldNbrNonZeroSingularValues, true);

       for (int i = 0; i < MatDim; i++)
         for (int j = 0; j < OldNbrNonZeroSingularValues; ++j)
           {
             TmpB0.GetMatrixElement(i, j, Tmp);
             TmpBUD->SetMatrixElement(i, j, Tmp);
             TmpB1.GetMatrixElement(i, j, Tmp);
             TmpBUD->SetMatrixElement(i, j + OldNbrNonZeroSingularValues, Tmp);
           }


       delete NewUDMatrix;

       //cout<<"TmpBUD"<<endl;
       //cout<<*TmpBUD<<endl;


       UMatrix = new RealMatrix (MatDim, MatDim, true);
       VMatrix = new RealMatrix (2 * OldNbrNonZeroSingularValues, 2 * OldNbrNonZeroSingularValues, true);

       SingularValues = TmpBUD->SingularValueDecomposition(*UMatrix, *VMatrix, false);
 
       SingularDimension = 2 * OldNbrNonZeroSingularValues;
       if (MatDim < 2 * OldNbrNonZeroSingularValues)
          SingularDimension = MatDim;

       //cout<<"Singular values: ";
       //for (int i = 0; i < SingularDimension; i++)
       //  cout<<SingularValues[i]<<" ";
       //cout<<endl;

       int NbrNonZeroSingularValues = 0;
       for (int i = 0; i < SingularDimension; i++)
         if (SingularValues[i] > CutOff)
           NbrNonZeroSingularValues++; 

       //cout<<"non zero = "<<NbrNonZeroSingularValues<<endl;

/*

       cout<<"UMatrix "<<endl;
       cout<<*UMatrix<<endl;

       cout<<"VMatrix "<<endl;
       cout<<*VMatrix<<endl;

        RealMatrix DTest (MatDim, 2 * OldNbrNonZeroSingularValues, true);
        for (int i = 0; i < SingularDimension; ++i)
          DTest.SetMatrixElement(i, i, SingularValues[i]);

       cout<<"DTest "<<endl;
       cout<<DTest<<endl;


        RealMatrix SingTest (MatDim, 2 * OldNbrNonZeroSingularValues, true);
        SingTest.Copy(*UMatrix);
        SingTest.Multiply(DTest);
        SingTest.Multiply(*VMatrix);
        cout<<"Test SVD"<<endl;
        cout<<SingTest<<endl;

        exit(1);

*/
        UMatrix->Resize(MatDim, NbrNonZeroSingularValues);
        VMatrix->Resize(NbrNonZeroSingularValues, 2 * OldNbrNonZeroSingularValues);
       
        NewDMatrix = new RealMatrix (NbrNonZeroSingularValues, NbrNonZeroSingularValues, true);
        for (int i = 0; i < NbrNonZeroSingularValues; ++i)
          NewDMatrix->SetMatrixElement(i, i, SingularValues[i]);
    
        NewUDMatrix = new RealMatrix (MatDim, 2 * OldNbrNonZeroSingularValues, true);
        NewUDMatrix->Copy(UMatrix->Multiply(*NewDMatrix));
        NewUDMatrix->Resize(MatDim, NbrNonZeroSingularValues);

        //cout<<"NewUDMatrix "<<*NewUDMatrix<<endl;
 
        delete UMatrix;
        delete VMatrix;
        delete TmpBUD;
        delete NewDMatrix;
        delete[] SingularValues;

        OldNbrNonZeroSingularValues = NbrNonZeroSingularValues;
       
    }

  
  RealMatrix NewUDOld (*NewUDMatrix);
  //cout<<"Carry away NewUD "<<endl;
  //cout<<NewUDOld<<endl;

  
  cout<<"----------Start SVD from left end-------------"<<endl;



  TmpBUD = new RealMatrix (2, MatDim, true);
  for (int i = 0; i < MatDim; i++)
    {
      BMatrices[0].GetMatrixElement(MPSRowIndex, i, Tmp);
      TmpBUD->SetMatrixElement(0, i, Tmp);
      BMatrices[1].GetMatrixElement(MPSRowIndex, i, Tmp);
      TmpBUD->SetMatrixElement(1, i, Tmp);
    }
  //cout<<"TmpBUD"<<endl;
  //cout<<*TmpBUD<<endl;


  UMatrix = new RealMatrix (2, 2, true);
  VMatrix = new RealMatrix (MatDim, MatDim, true);

  SingularValues = TmpBUD->SingularValueDecomposition(*UMatrix, *VMatrix, false);

  SingularDimension = 2;
  if (MatDim < 2)
    SingularDimension = MatDim;
 
  //cout<<"Singular values: ";
  //for (int i = 0; i < SingularDimension; i++)
  //  cout<<SingularValues[i]<<" ";
  //cout<<endl;

  NbrNonZeroSingularValues = 0;
  for (int i = 0; i < SingularDimension; i++)
   if (SingularValues[i] > CutOff)
     NbrNonZeroSingularValues++;

  //cout<<"non zero = "<<NbrNonZeroSingularValues<<endl;

  //cout<<"UMatrix "<<endl;
  //cout<<*UMatrix<<endl;

  //cout<<"VMatrix "<<endl;
  //cout<<*VMatrix<<endl;


  UMatrix->Resize(2, NbrNonZeroSingularValues);
  VMatrix->Resize(NbrNonZeroSingularValues, MatDim);
  
  NewDMatrix = new RealMatrix (NbrNonZeroSingularValues, NbrNonZeroSingularValues, true);
  for (int i = 0; i < NbrNonZeroSingularValues; ++i)
    NewDMatrix->SetMatrixElement(i, i, SingularValues[i]);


  NewDVMatrix = new RealMatrix (NbrNonZeroSingularValues, MatDim, true);
  NewDVMatrix->Copy(NewDMatrix->Multiply(*VMatrix));

  //cout<<"NewDV"<<*NewDVMatrix<<endl;


  delete UMatrix;
  delete VMatrix;
  delete NewDMatrix;
  delete TmpBUD;
  delete[] SingularValues;

  OldNbrNonZeroSingularValues = NbrNonZeroSingularValues;

  for (int TmpCut = 2; TmpCut <= EntCut; TmpCut++)
    {
       RealMatrix TmpB0 (MatDim, MatDim);
       TmpB0.Copy(*NewDVMatrix);
       TmpB0.Multiply(DenseB0);
       TmpB0.Resize(OldNbrNonZeroSingularValues,MatDim);

       RealMatrix TmpB1 (MatDim, MatDim);
       TmpB1.Copy(*NewDVMatrix);
       TmpB1.Multiply(DenseB1);
       TmpB1.Resize(OldNbrNonZeroSingularValues, MatDim);

       TmpBUD = new RealMatrix (2 * OldNbrNonZeroSingularValues, MatDim, true);

       for (int i = 0; i < OldNbrNonZeroSingularValues; i++)
         for (int j = 0; j < MatDim; ++j)
           {
             TmpB0.GetMatrixElement(i, j, Tmp);
             TmpBUD->SetMatrixElement(i, j, Tmp);
             TmpB1.GetMatrixElement(i, j, Tmp);
             TmpBUD->SetMatrixElement(i + OldNbrNonZeroSingularValues, j, Tmp);
           }


       delete NewDVMatrix;

       //cout<<"TmpB0"<<endl;
       //cout<<TmpB0<<endl;
       //cout<<"TmpB1"<<endl;
       //cout<<TmpB1<<endl;

       //cout<<"TmpBUD"<<endl;
       //cout<<*TmpBUD<<endl;


       UMatrix = new RealMatrix (2 * OldNbrNonZeroSingularValues, 2 * OldNbrNonZeroSingularValues, true);
       VMatrix = new RealMatrix (MatDim, MatDim, true);


       //cout<<"TmpBUD "<<TmpBUD->GetNbrRow()<<" "<<TmpBUD->GetNbrColumn()<<" ; "<<(2 * OldNbrNonZeroSingularValues)<<" "<<MatDim<<endl;

       SingularValues = TmpBUD->SingularValueDecomposition(*UMatrix, *VMatrix, false);
 
       SingularDimension = 2 * OldNbrNonZeroSingularValues;
       if (MatDim < 2 * OldNbrNonZeroSingularValues)
          SingularDimension = MatDim;

       //cout<<"Singular values: ";
       //for (int i = 0; i < SingularDimension; i++)
       //  cout<<SingularValues[i]<<" ";
       //cout<<endl;

       int NbrNonZeroSingularValues = 0;
       for (int i = 0; i < SingularDimension; i++)
         if (SingularValues[i] > CutOff)
           NbrNonZeroSingularValues++; 

       //cout<<"non zero = "<<NbrNonZeroSingularValues<<endl;

       UMatrix->Resize(2 * OldNbrNonZeroSingularValues, NbrNonZeroSingularValues);
       VMatrix->Resize(NbrNonZeroSingularValues, MatDim);

       //cout<<"UMatrix "<<endl;
       //cout<<*UMatrix<<endl;

       //cout<<"VMatrix "<<endl;
       //cout<<*VMatrix<<endl;
       

       NewDMatrix = new RealMatrix (NbrNonZeroSingularValues, NbrNonZeroSingularValues, true);
       for (int i = 0; i < NbrNonZeroSingularValues; ++i)
          NewDMatrix->SetMatrixElement(i, i, SingularValues[i]);

       //cout<<"DMatrix "<<endl;
       //cout<<*NewDMatrix<<endl;

       NewDVMatrix = new RealMatrix (NbrNonZeroSingularValues, MatDim, true);
       NewDVMatrix->Copy(NewDMatrix->Multiply(*VMatrix));

       //cout<<"NewDV"<<*NewDVMatrix<<endl;

        delete UMatrix;
        delete VMatrix;
        delete TmpBUD;
        delete NewDMatrix;
        delete[] SingularValues;

        OldNbrNonZeroSingularValues = NbrNonZeroSingularValues;    
 
    }

  cout<<"Completed..........................."<<endl;

  RealMatrix EntMatrixOnSite (MatDim, MatDim, true);
  EntMatrixOnSite.Copy(*NewDVMatrix);
  EntMatrixOnSite.Multiply(NewUDOld);
  EntMatrixOnSite.Resize(NewDVMatrix->GetNbrRow(), NewUDOld.GetNbrColumn());
  SingularValues = EntMatrixOnSite.SingularValueDecomposition();

  SingularDimension = NewDVMatrix->GetNbrRow();
  if (NewUDOld.GetNbrColumn() < NewDVMatrix->GetNbrRow())
    SingularDimension = NewUDOld.GetNbrColumn();

  double Trace = 0.0;
  for (int i = 0; i < SingularDimension; i++)
    {
      SingularValues[i] *= SingularValues[i];
      Trace += SingularValues[i];
    }  

  cout<<"Entanglement levels: "<<endl;
  for (int i = 0; i < SingularDimension; i++)
    {
      if (fabs(SingularValues[i]) > CutOff)
        cout<<SingularValues[i]/Trace<<endl;
    } 
  cout<<endl; 

  delete NewDVMatrix;
  delete[] SingularValues;  

/*
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

*/
  File.close();
 
  return 0;
}
