#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/FermionOnSphereMPSWrapper.h"
#include "HilbertSpace/FermionOnCylinderMPSWrapper.h"

#include "MathTools/ClebschGordanCoefficients.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"
#include "MathTools/FactorialCoefficient.h"

#include "Tools/FQHEMPS/FQHEMPSMatrixManager.h"
#include "Tools/FQHEMPS/AbstractFQHEMPSMatrix.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Matrix/SparseRealMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

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

double OverlapEig(RealVector& EigV, SparseRealMatrix& TmpM, int TmpSectorDim, int i);

int main(int argc, char** argv)
{
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
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "memory", "amount of memory that can used for precalculations (in Mb)", 500);
  (*OutputGroup) += new SingleStringOption  ('o', "output-file", "output file name");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereMPSCorrelation -h" << endl;
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
      sprintf(TmpFileName, "fermions_laughlin%ld_plevel_%ld_n_%d_2s_%d_lz_%d.0.ent", Manager.GetInteger("laughlin-index"),
		  Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz);
      File.open(TmpFileName, ios::binary | ios::out);     
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
  cout << "Requested memory for sparse matrix multiplications = " << ((MaxTmpMatrixElements * (2l * sizeof(double) + sizeof(int))) >> 20) << "Mb" << endl;
  double* TmpMatrixElements = new double [MaxTmpMatrixElements];
  int* TmpColumnIndices = new int [MaxTmpMatrixElements];
  double * TmpElements = new double [BMatrices[0].GetNbrRow()];

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

      OverlapMatrix = TmpMatrix2 + TmpMatrix3;
    }


  cout<<"Compute density matrix in nonorthogonal basis (full space dimension, will be stored)"<<endl;

  SparseRealMatrix rhoM;

  TmpMatrix2 = Multiply(FinalBMatrices[0], FinalBMatrices[0].Transpose(), 
 			TmpMatrixElements, TmpColumnIndices, TmpElements); 

  if (FinalBMatrices[1].ComputeNbrNonZeroMatrixElements() != 0)
   {
     TmpMatrix3 = Multiply(FinalBMatrices[1], FinalBMatrices[1].Transpose(), 
			   TmpMatrixElements, TmpColumnIndices, TmpElements); 
     rhoM = TmpMatrix2 + TmpMatrix3;
   }
  else
    {
      rhoM = TmpMatrix2;
    }

  for (int i = (NbrFluxQuanta - 1); i >= EntCut; i--)
    {
      TmpMatrix1 = Multiply(BMatrices[0], rhoM, 
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix2 = Multiply(TmpMatrix1, ConjugateBMatrices[0],  
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix1 = Multiply(BMatrices[1], rhoM, 
			    TmpMatrixElements, TmpColumnIndices, TmpElements); 
      TmpMatrix3 = Multiply(TmpMatrix1, ConjugateBMatrices[1],  
			    TmpMatrixElements, TmpColumnIndices, TmpElements);   
      rhoM = TmpMatrix2 + TmpMatrix3;
    }

  //Free up some space that is no longer needed (this needs to be done in a cleaner way)

  delete[] BMatrices; 
  delete[] ConjugateBMatrices; 
  delete[] FinalBMatrices;
  delete[] TmpMatrixElements;
  delete[] TmpColumnIndices;
  delete[] TmpElements;

  cout<<"Proceed to calculate ES per momentum P sector (all N sectors)"<<endl;

  double TraceRho = 0.0;
  double* RhoEigenvalues = new double [MatDim];
  int* RhoNSector = new int [MatDim];
  int* RhoPSector = new int [MatDim];

  int NbrNValue = ((2 * LambdaMax) + LaughlinIndex);

  for (int i =0 ; i < MatDim; i++)
   {
     RhoEigenvalues[i] = 0.0; 
     RhoNSector[i] = 0; 
     RhoPSector[i]=0;
   }

  int eigenvaluecounter = 0;
  for (int NSector = 0; NSector < NbrNValue; NSector++)
    {
      for (int MomentumSector = 0; MomentumSector <= LambdaMax; MomentumSector++)
	{
	  SparseRealMatrix TmpMBlock = MPSMatrix->ExtractBlock(OverlapMatrix, MomentumSector, NSector, MomentumSector, NSector);
	  
	  int TmpSectorDim = TmpMBlock.GetNbrRow();
	  
	  RealSymmetricMatrix HRep (TmpMBlock);
	  RealDiagonalMatrix TmpDiag (TmpSectorDim);
	  RealMatrix Q(TmpSectorDim, TmpSectorDim);
	  HRep.LapackDiagonalize(TmpDiag, Q);
      
	  int NonZeroEig = 0;
	  for (int j = TmpSectorDim - 1; j >= 0; --j)
	    {
	      if (fabs(TmpDiag[j]) > CutOff)
		{  
		  NonZeroEig++;
		  Q[j] /= sqrt(fabs(TmpDiag[j]));
		}  
	    }
	  
	  if (NonZeroEig > 0)
	    {
	      
	      //cout<<"-------------------- Start computing rho in the new basis --------------------"<<endl;
	      
	      RealMatrix NewrhoM(NonZeroEig, NonZeroEig, true);
	      
	      SparseRealMatrix rhoMBlock = MPSMatrix->ExtractBlock(rhoM, MomentumSector, NSector, MomentumSector, NSector);
	      //cout<<"rhoMBlock "<<rhoMBlock<<endl;
	      
	      
	      int rowIndex, columnIndex;
	      
	      int offset = TmpSectorDim - NonZeroEig;
	      for (int i = TmpSectorDim - 1; i >= offset; i--)
		for (int j = TmpSectorDim - 1; j >= offset; j--)
		  {
		    double TmpEl = 0.0;
		    for (int ii = 0; ii < TmpSectorDim; ++ii)
		      for (int jj = 0; jj < TmpSectorDim; ++jj)
			{
			  double Tmp;
			  rhoMBlock.GetMatrixElement(ii, jj, Tmp);
			  TmpEl += Tmp * OverlapEig(Q[j], TmpMBlock, TmpSectorDim, jj) * OverlapEig(Q[i], TmpMBlock, TmpSectorDim, ii);
			}
		    
		    NewrhoM.SetMatrixElement(i-offset, j-offset, TmpEl); 
		  }
	      
	      RealSymmetricMatrix HRepRho ((Matrix&) NewrhoM);
	      
	      RealDiagonalMatrix TmpDiagRho (NonZeroEig);
	      RealMatrix QRho(NonZeroEig, NonZeroEig);
	      HRepRho.LapackDiagonalize(TmpDiagRho, QRho);
	      
	      cout<<"------------sector P = "<<MomentumSector<<" N = "<<((EntCut - (NSector - (2 * LambdaMax + LaughlinIndex - 1)/2)))/LaughlinIndex << "---------------"<<endl;
	      
	      cout.precision(14); 
	      
	      double Sum = 0.0;
	      for (int j = 0; j < NonZeroEig; ++j)
		{
		  //cout<<"Eigenvalue "<<TmpDiagRho[j]<<endl;
		  TraceRho += TmpDiagRho[j];
		  RhoEigenvalues[eigenvaluecounter]=TmpDiagRho[j];
		  RhoNSector[eigenvaluecounter]=NSector;
		  RhoPSector[eigenvaluecounter]=MomentumSector; 
		  eigenvaluecounter++;
		}
	    }
	}
    }
  
  cout<<"Trace rho = "<<TraceRho<<endl;

  for (int i=0; i<MatDim; ++i)
    if ((fabs(RhoEigenvalues[i]) > CutOff) && ((((EntCut - (RhoNSector[i] - (2 * LambdaMax + LaughlinIndex - 1)/2)))/LaughlinIndex) == Na))
      {
        cout<<"P= "<<RhoPSector[i]<<" N= "<<RhoNSector[i]<<" "<<RhoEigenvalues[i]/TraceRho<<endl;  
        File<<RhoPSector[i]<<" "<<RhoNSector[i]<<" "<<RhoEigenvalues[i]/TraceRho<<endl;
      }

 delete[] RhoEigenvalues;
 delete[] RhoPSector;
 delete[] RhoNSector;

  File.close();
 
  return 0;
}

double OverlapEig(RealVector& EigV, SparseRealMatrix& TmpM, int TmpSectorDim, int i)
{
  double Tmp;
  double ov = 0.0;
  for(int ii = 0; ii < TmpSectorDim; ii++)
   {
      TmpM.GetMatrixElement(ii, i, Tmp); 
      ov += EigV[ii] * Tmp;
   } 
 return ov;
}

