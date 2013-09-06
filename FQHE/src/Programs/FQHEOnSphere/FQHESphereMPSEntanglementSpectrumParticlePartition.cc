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
#include "GeneralTools/ConfigurationParser.h"

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


// compute the entanglement spectrum from the overlap matrices of the left and right parts
// 
// mPSMatrix = pointer tothe MPS matrix
// leftOverlapMatrix= reference on the left overlap matrix
// rightOverlapMatrix= reference on the right overlap matrix
// leftPSector = P sector that has to be selected for the left part
// leftQSector = Q sector that has to be selected for the left part
// rightPSector = P sector that has to be selected for the right part
// rightQSector = Q sector that has to be selected for the right part
// eigenvalueError = relative error on the eigenvalues below which an eigenvalue is considered to be equal to zero
RealDiagonalMatrix FQHEMPSEvaluatePartialEntanglementSpectrum(AbstractFQHEMPSMatrix* mPSMatrix, SparseRealMatrix& leftOverlapMatrix, SparseRealMatrix& rightOverlapMatrix, 
						       int leftPSector, int leftQSector, int rightPSector, int rightQSector, double eigenvalueError);

// compute the entanglement spectrum from the overlap matrices of the left and right parts
// 
// mPSMatrix = pointer tothe MPS matrix
// leftOverlapMatrix= reference on the left overlap matrix
// rightOverlapMatrix= reference on the right overlap matrix
// normalizationMatrix = optional normalization matrix (for PES and RES)
// leftPSector = P sector that has to be selected for the left part
// leftCFTSector = CFT sector that has to be selected for the left part
// leftQSector = Q sector that has to be selected for the left part
// rightPSector = P sector that has to be selected for the right part
// rightCFTSector = CFT sector that has to be selected for the right part
// rightQSector = Q sector that has to be selected for the right part
// eigenvalueError = relative error on the eigenvalues below which an eigenvalue is considered to be equal to zero
RealDiagonalMatrix FQHEMPSEvaluatePartialEntanglementSpectrum(AbstractFQHEMPSMatrix* mPSMatrix, SparseRealMatrix& leftOverlapMatrix, SparseRealMatrix& rightOverlapMatrix, 
							      SparseRealMatrix& normalizationMatrix,
							      int leftPSector, int leftCFTSector, int leftQSector, int rightPSector, int rightCFTSector, int rightQSector, 
							      double eigenvalueError);


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
  (*SystemGroup) += new BooleanOption  ('\n', "realspace-cut", "use real space partition instead of particle partition");
  (*SystemGroup) += new SingleStringOption  ('\n', "realspace-partition", "geometrical weights that define the real spac partition");
  (*SystemGroup) += new BooleanOption ('\n', "use-singlestate", "use a single real eigenstate of the E matrix  when evaluating the infinite entanglement spectrum");
  (*SystemGroup) += new BooleanOption ('\n', "orbital-es", "compute the orbital entanglement spectrum");
  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-orbitals", "number of orbitals for the A part (i negative, use (N_phi+1) / 2)", -1);
  (*SystemGroup) += new SingleIntegerOption ('\n', "nbr-fluxquanta", "set the total number of flux quanta and deduce the root partition instead of using the reference-file", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "left-eigenstate", "file containing the transfer matrix left eigenstate");
  (*SystemGroup) += new SingleStringOption  ('\n', "right-eigenstate", "file containing the transfer matrix right eigenstate");
  (*SystemGroup) += new BooleanOption  ('\n', "diagonal-block", "transfer matrix eigenstates are computed only from the block diagonal in P, CFT sector and Q (override autodetect from eigenvector file names)");
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
  if ((Manager.GetBoolean("infinite-cylinder") == false) && (Manager.GetInteger("nbr-fluxquanta") <= 0) && 
      (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false))
    {
      return -1;
    }
  if (Manager.GetInteger("nbr-fluxquanta") > 0)
    {
      NbrFluxQuanta = Manager.GetInteger("nbr-fluxquanta");
    }

  int Na = Manager.GetInteger("na");

  bool CylinderFlag = Manager.GetBoolean("normalize-cylinder");

  int LandauLevel = 0;

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta, Architecture.GetArchitecture()); 
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
  NbrParticles = MPSMatrix->GetMatrixNaturalNbrParticles(NbrFluxQuanta, Manager.GetBoolean("use-padding"));
  NbrEigenstates = MPSMatrix->GetTransferMatrixLargestEigenvalueDegeneracy();
  int PLevel = MPSMatrix->GetTruncationLevel();
  int NbrCFTSectors = MPSMatrix->GetNbrCFTSectors();
  ofstream File;
  File.precision(14);
  ofstream File2;
  File2.precision(14);
  
  char* Extension = new char[8];  
  if (Manager.GetBoolean("orbital-es"))
    {
      sprintf (Extension, "ent");
    }
  else
    {
      sprintf (Extension, "parent");
    }
  if (Manager.GetString("output-file") != 0)
    {
      File.open(Manager.GetString("output-file"), ios::binary | ios::out);
    }
  else
    {
      char* TmpFileName = new char [512];
      char* TmpFileName2 = new char [512];
      char* StateName = new char [256];
      strcpy(StateName, MPSMatrix->GetName());
      if (CylinderFlag == true)
	{
	  if (Manager.GetBoolean("infinite-cylinder"))
	    {
	      sprintf(TmpFileName, "fermions_infinite_cylinder_%s_perimeter_%f_plevel_%ld_n_0_2s_0_lz_0.0.full.%s", StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"), Extension);
	      sprintf(TmpFileName2, "fermions_infinite_cylinder_%s_perimeter_%f_plevel_%ld_n_0_2s_0_lz_0.0.%s", StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"), Extension);
	    }
	  else
	    {
	      sprintf(TmpFileName, "fermions_cylinder_%s_perimeter_%f_plevel_%ld_n_%d_2s_%d_lz_%d.0.full.%s", StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz, Extension);
	      sprintf(TmpFileName2, "fermions_cylinder_%s_perimeter_%f_plevel_%ld_n_%d_2s_%d_lz_%d.0.%s", StateName,
		      MPSMatrixManager.GetCylinderPerimeter(NbrFluxQuanta), Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz, Extension);
	    }
	}
      else
	{
	  if (Manager.GetBoolean("infinite-cylinder"))
	    {
	      sprintf(TmpFileName, "fermions_infinite_%s_plevel_%ld_n_0_2s_0_lz_0.0.full.%s", StateName,
		      Manager.GetInteger("p-truncation"), Extension);
	      sprintf(TmpFileName2, "fermions_infinite_%s_plevel_%ld_n_0_2s_0_lz_0.0.%s", StateName,
		      Manager.GetInteger("p-truncation"), Extension);
	    }
	  else
	    {
	      sprintf(TmpFileName, "fermions_%s_plevel_%ld_n_%d_2s_%d_lz_%d.0.full.%s", StateName,
		      Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz, Extension);
	      sprintf(TmpFileName2, "fermions_%s_plevel_%ld_n_%d_2s_%d_lz_%d.0.%s", StateName,
		      Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz, Extension);
	    }
	}
      File.open(TmpFileName, ios::binary | ios::out);     
      File2.open(TmpFileName2, ios::binary | ios::out);     
   }

  if ((Manager.GetBoolean("orbital-es")) || (Manager.GetBoolean("realspace-cut")))
    {
      File << "# x    Q    P    lambda    -ln(lambda)" << endl;
    }
  else
    {
      File << "#  N    Lz    lambda" << endl;
    }

  if (Manager.GetBoolean("infinite-cylinder"))
    {
      if (Manager.GetString("left-eigenstate") == 0)
	{
	  cout << "The transfer matrix left eigenstate has to be provided in infinite-cylinder mode" << endl;
	  return 0;
	}
      if (Manager.GetString("right-eigenstate") == 0)
	{
	  cout << "The transfer matrix right eigenstate has to be provided in infinite-cylinder mode" << endl;
	  return 0;
	}
      ComplexVector LeftEigenstate;
      ComplexVector RightEigenstate;
      if (LeftEigenstate.ReadVector(Manager.GetString("left-eigenstate")) == false)
	{
	  cout << "can't read " << Manager.GetString("left-eigenstate") << endl;
	  return 0;
	}      
      if (RightEigenstate.ReadVector(Manager.GetString("right-eigenstate")) == false)
	{
	  cout << "can't read " << Manager.GetString("right-eigenstate") << endl;
	  return 0;
	}
      Complex Factor = EuclidianScalarProduct(LeftEigenstate, RightEigenstate);
      double InvFactor = 1.0 / Factor.Re;
      int TmpDimension = BMatrices[0].GetNbrRow();
      RealMatrix TmpFullLeftOverlapMatrix(TmpDimension, TmpDimension, true);
      RealMatrix TmpFullRightOverlapMatrix(TmpDimension, TmpDimension, true);
      SparseRealMatrix NormalizationMatrix;
      if ((Manager.GetBoolean("diagonal-block") == false) && 
	  ((strstr(Manager.GetString("left-eigenstate"), "_diagblock_") == 0) || (strstr(Manager.GetString("right-eigenstate"), "_diagblock_") == 0)))
	{
	  if (LeftEigenstate.GetVectorDimension() != (TmpDimension * TmpDimension))
	    {
	      cout << "error, left eigenstate does not have the expected dimension" << endl;
	      return 0;
	    }
	  if (RightEigenstate.GetVectorDimension() != (TmpDimension * TmpDimension))
	    {
	      cout << "error, right eigenstate does not have the expected dimension" << endl;
	      return 0;
	    }
	  for (int i = 0; i < TmpDimension; ++i)
	    for (int j = 0; j < TmpDimension; ++j)
	      {
		TmpFullLeftOverlapMatrix.SetMatrixElement(i, j, InvFactor * LeftEigenstate[i * TmpDimension + j].Re); 
		TmpFullRightOverlapMatrix.SetMatrixElement(i, j, InvFactor * RightEigenstate[i * TmpDimension + j].Re); 
	      }
	}
      else
	{
	  int TmpBlockDimension = 0;
	  for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
	    {
	      for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
		{
		  int MinQValue = 0;
		  int MaxQValue = 0;
		  MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
		  for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
		    {
		      int TmpLocalDimension = MPSMatrix->GetBondIndexRange(CurrentPLevel, QValue, CurrentCFTSector);
		      TmpBlockDimension += TmpLocalDimension * TmpLocalDimension;
		    }
		}
	    }
	  if (LeftEigenstate.GetVectorDimension() != TmpBlockDimension)
	    {
	      cout << "error, left eigenstate does not have the expected dimension" << endl;
	      return 0;
	    }
	  if (RightEigenstate.GetVectorDimension() != TmpBlockDimension)
	    {
	      cout << "error, right eigenstate does not have the expected dimension" << endl;
	      return 0;
	    }

	  int TmpBlockPosition = 0;
	  for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
	    {
	      for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
		{
		  int MinQValue = 0;
		  int MaxQValue = 0;
		  MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
		  for (int QValue = MinQValue; QValue <= MaxQValue; ++QValue)
		    {
		      int TmpLocalDimension = MPSMatrix->GetBondIndexRange(CurrentPLevel, QValue, CurrentCFTSector);
		      for (int i = 0; i < TmpLocalDimension; ++i)
			{
			  int Tmp1 = MPSMatrix->GetBondIndexWithFixedChargePLevelCFTSector(i, CurrentPLevel, QValue, CurrentCFTSector);
			  for (int j = 0; j < TmpLocalDimension; ++j)
			    {
			      int Tmp2 = MPSMatrix->GetBondIndexWithFixedChargePLevelCFTSector(j, CurrentPLevel, QValue, CurrentCFTSector);
			      TmpFullLeftOverlapMatrix.SetMatrixElement(Tmp1, Tmp2, InvFactor * LeftEigenstate[TmpBlockPosition].Re); 
			      TmpFullRightOverlapMatrix.SetMatrixElement(Tmp1, Tmp2, InvFactor * RightEigenstate[TmpBlockPosition].Re); 
			      ++TmpBlockPosition;
			    }
			}
		    }
		}
	    }
	}
      if (TmpFullLeftOverlapMatrix.Tr() < 0.0)
	TmpFullLeftOverlapMatrix *= -1.0;
      if (TmpFullRightOverlapMatrix.Tr() < 0.0)
	TmpFullRightOverlapMatrix *= -1.0;
      if (TmpFullLeftOverlapMatrix.IsSymmetric())
	{
	  cout << "FullLeftOverlapMatrix is symmetric" << endl;
	}
      if (TmpFullRightOverlapMatrix.IsSymmetric())
	{
	  cout << "FullRightOverlapMatrix is symmetric" << endl;
	}
      SparseRealMatrix FullLeftOverlapMatrix (TmpFullLeftOverlapMatrix);
      SparseRealMatrix FullRightOverlapMatrix (TmpFullRightOverlapMatrix);
      double Error = 1e-13;
      double LeftEigenvalueError = 0.0;
      double RightEigenvalueError = 0.0;
      LeftEigenvalueError = Error;
      RightEigenvalueError = Error;
      double TotalTraceThoA = 0;
      if (Manager.GetBoolean("orbital-es"))
	{
	  double**** EntanglementSpectrum = new double***[PLevel + 1];
	  int*** EntanglementSpectrumDimension = new int**[PLevel + 1];
	  for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
	    {
	      EntanglementSpectrumDimension[CurrentPLevel] = new int*[NbrCFTSectors];
	      EntanglementSpectrum[CurrentPLevel] = new double**[NbrCFTSectors];	      	      
	      for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
		{
		  int LocalMinQValue;
		  int LocalMaxQValue;
		  MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, LocalMinQValue, LocalMaxQValue);
		  if (LocalMinQValue <=  LocalMaxQValue)
		    {
		      EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector] = new int[LocalMaxQValue - LocalMinQValue + 1];
		      EntanglementSpectrum[CurrentPLevel][CurrentCFTSector] = new double*[LocalMaxQValue - LocalMinQValue + 1];	      
		    }
		  else
		    {
		      EntanglementSpectrum[CurrentPLevel][CurrentCFTSector] = 0;
		    }
		  for (int LocalQValue =  LocalMinQValue; LocalQValue <= LocalMaxQValue; ++LocalQValue)
		    {
		      cout << "computing sector P=" << CurrentPLevel<< " CFT=" << CurrentCFTSector << " Q=" << LocalQValue << endl;
		      RealDiagonalMatrix TmpRhoADiag = FQHEMPSEvaluatePartialEntanglementSpectrum(MPSMatrix, FullLeftOverlapMatrix, FullRightOverlapMatrix, 
												  NormalizationMatrix,
												  CurrentPLevel, CurrentCFTSector, LocalQValue, 
												  CurrentPLevel, CurrentCFTSector, LocalQValue, Error);
		      EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = 0;
		      if (TmpRhoADiag.GetNbrRow() > 0)
			{
			  for (int i = 0 ; i < TmpRhoADiag.GetNbrRow(); ++i)
			    {
			      if (TmpRhoADiag[i] > 0.0)
				EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]++;
			    }
			  if (EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] > 0)
			    {
			      EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = new double[EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]];
			      EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = 0;
			      for (int i = 0 ; i < TmpRhoADiag.GetNbrRow(); ++i)
				{
				  if (TmpRhoADiag[i] > 0.0)
				    {
				      EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]] = TmpRhoADiag[i];
				      EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]++;
				      TotalTraceThoA += TmpRhoADiag[i];
				    }
				}
			      SortArrayDownOrdering<double>(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue],
							    EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]);
			    }
			}
		    }
		}
	    }
	  
	  
	  int GlobalMinQValue = 100000;
	  int GlobalMaxQValue = 0;
	  for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
	    {
	      for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
		{
		  int MinQValue = 0;
		  int MaxQValue = 0;
		  MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
		  if (MinQValue < GlobalMinQValue)
		    GlobalMinQValue = MinQValue;
		  if (MaxQValue > GlobalMaxQValue)
		    GlobalMaxQValue = MaxQValue;		  
		}
	    }
	  double EntanglementEntropy = 0.0;
	  for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
	    {
	      for (int LocalQValue =  GlobalMinQValue; LocalQValue <= GlobalMaxQValue; ++LocalQValue)
		{
		  for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
		    {
		      int LocalMinQValue;
		      int LocalMaxQValue;
		      MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, LocalMinQValue, LocalMaxQValue);
		      if ((LocalQValue >= LocalMinQValue) && (LocalQValue <= LocalMaxQValue))
			{
			  for (int i = 0; i < EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]; ++i)
			    {
			      File << CurrentCFTSector  << " " << LocalQValue << " " 
				   << CurrentPLevel << " "
				   <<  (EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)  
				   <<  " " << (-log(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)) << endl;
			      EntanglementEntropy -= (log(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)
						      * EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA);
			    }
			}
		    }
		}
	    }
	  cout << "S_A=" << EntanglementEntropy << endl;
	  File2 << "inf" << " " << EntanglementEntropy << " 1" << endl;
	}
      File.close();
      File2.close();
      
      cout << "Tr(rho_A)=" << TotalTraceThoA << endl;
      return 0;
    }
  
  
  int QValue = 3;
  
  cout << "MPSRowIndex = " << MPSRowIndex << endl;
  cout << "MPSColumnIndex = " << MPSColumnIndex << endl;
  
  SparseRealMatrix FullLeftOverlapMatrix (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
  SparseRealMatrix FullRightOverlapMatrix (BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());

  long MaxTmpMatrixElements = (((long) BMatrices[0].GetNbrRow()) * 
				((long) BMatrices[0].GetNbrRow() / 1l));
  cout << "Requested memory for sparse matrix multiplications = " << ((MaxTmpMatrixElements * (2l * sizeof(double) + sizeof(int))) >> 20) << "Mb" << endl;
  double* TmpMatrixElements = new double [MaxTmpMatrixElements];
  int* TmpColumnIndices = new int [MaxTmpMatrixElements];
  double* TmpElements = new double [BMatrices[0].GetNbrRow()];

  int MaxNbrFluxQuantaA = PLevel + QValue * (Na - 1) + ((QValue - 1) / 2);
  if (MaxNbrFluxQuantaA > NbrFluxQuanta)
    {
      MaxNbrFluxQuantaA = NbrFluxQuanta;
    }
  int MaxNbrFluxQuantaB = PLevel + QValue * (NbrParticles - Na - 1)+ ((QValue - 1) / 2);
  if (MaxNbrFluxQuantaB > NbrFluxQuanta)
    {
      MaxNbrFluxQuantaB = NbrFluxQuanta;
    }

  double* WeightAOrbitals = 0;
  double* WeightBOrbitals = 0;
  if (Manager.GetBoolean("realspace-cut") == true)
    {
      ConfigurationParser RealSpaceWeights;
      if (RealSpaceWeights.Parse(Manager.GetString("realspace-partition")) == false)
	{
	  RealSpaceWeights.DumpErrors(cout) << endl;
	  return -1;
	}
      double* TmpSquareWeights = 0;
      int TmpNbrOrbitals = 0;
      if (RealSpaceWeights.GetAsDoubleArray("OrbitalSquareWeights", ' ', TmpSquareWeights, TmpNbrOrbitals) == false)
	{
	  cout << "OrbitalSquareWeights is not defined or as a wrong value" << endl;
	  return -1;
	}
      if (TmpNbrOrbitals > (NbrFluxQuanta + 1))
	{
	  cout << "error, the number of weights (" << TmpNbrOrbitals << ") cannot exceed the number of orbitals (" << (NbrFluxQuanta + 1) << ")" << endl;
	  return -1;
	}
      int NbrAOrbitals = (NbrFluxQuanta + 1 + TmpNbrOrbitals) / 2;
      WeightAOrbitals = new double [NbrAOrbitals];
      for (int i = 0; i < (NbrAOrbitals - TmpNbrOrbitals); ++i)
	WeightAOrbitals[i] = 1.0;
      for (int i = NbrAOrbitals - TmpNbrOrbitals; i < NbrAOrbitals; ++i)
	{
	  WeightAOrbitals[i] = sqrt(TmpSquareWeights[i - NbrAOrbitals + TmpNbrOrbitals]);
	}
      int NbrBOrbitals = (NbrFluxQuanta + 1 + TmpNbrOrbitals) / 2;
      WeightBOrbitals = new double [NbrBOrbitals];
      for (int i = 0; i < TmpNbrOrbitals; ++i)
	{
	  WeightBOrbitals[i] = sqrt(1.0 - TmpSquareWeights[i]);
	}
      for (int i = TmpNbrOrbitals; i < NbrBOrbitals; ++i)
	{
	  WeightBOrbitals[i] = 1.0;
	}     
      MaxNbrFluxQuantaA = NbrAOrbitals - 1;
      MaxNbrFluxQuantaB = NbrBOrbitals - 1;
    }


  if (Manager.GetBoolean("orbital-es") == true)
    {
      if ((Manager.GetInteger("nbr-orbitals") < 0) || (Manager.GetInteger("nbr-orbitals") > (NbrFluxQuanta + 1)))
	{
	  MaxNbrFluxQuantaA = ((NbrFluxQuanta + 1) / 2) - 1;
	}
      else
	{
	  MaxNbrFluxQuantaA = Manager.GetInteger("nbr-orbitals") - 1;
	}
      MaxNbrFluxQuantaB = NbrFluxQuanta - MaxNbrFluxQuantaA - 1;
    }

  cout << "keeping " << (MaxNbrFluxQuantaA + 1) << " orbitals in A" << endl; 
  FullLeftOverlapMatrix.SetMatrixElement(MPSRowIndex, MPSRowIndex, 1.0);  
  for (int i = 0; i <= MaxNbrFluxQuantaA; ++i)
    {
      SparseRealMatrix TmpMatrix2;
      SparseRealMatrix TmpMatrix3;
      TmpMatrix2 = Conjugate(&(ConjugateBMatrices[0]), &FullLeftOverlapMatrix, &(BMatrices[0]), 
			     TmpMatrixElements, TmpColumnIndices, MaxTmpMatrixElements, Architecture.GetArchitecture()); 
      TmpMatrix3 = Conjugate(&(ConjugateBMatrices[1]), &FullLeftOverlapMatrix, &(BMatrices[1]), 
			     TmpMatrixElements, TmpColumnIndices, MaxTmpMatrixElements, Architecture.GetArchitecture()); 
      if (WeightAOrbitals != 0)
 	{
	  cout << "WeightAOrbitals[" << i << "]=" << WeightAOrbitals[i] << endl;
 	  TmpMatrix3 *= WeightAOrbitals[i] * WeightAOrbitals[i];
 	}
      FullLeftOverlapMatrix = TmpMatrix2 + TmpMatrix3;
    }

  FullRightOverlapMatrix.SetMatrixElement(MPSColumnIndex, MPSColumnIndex, 1.0);  
  cout << "keeping " << (MaxNbrFluxQuantaB + 1) << " orbitals in B" << endl; 
  for (int i = 0; i <= MaxNbrFluxQuantaB; ++i)
    {
      SparseRealMatrix TmpMatrix2;
      SparseRealMatrix TmpMatrix3;
      TmpMatrix2 = Conjugate(&(BMatrices[0]), &FullRightOverlapMatrix, &(ConjugateBMatrices[0]), 
 			     TmpMatrixElements, TmpColumnIndices, MaxTmpMatrixElements, Architecture.GetArchitecture()); 
      TmpMatrix3 = Conjugate(&(BMatrices[1]), &FullRightOverlapMatrix, &(ConjugateBMatrices[1]), 
 			     TmpMatrixElements, TmpColumnIndices, MaxTmpMatrixElements, Architecture.GetArchitecture()); 
      if (WeightBOrbitals != 0)
 	{
	  cout << "WeightBOrbitals[" << i << "]=" << WeightBOrbitals[i] << endl;
 	  TmpMatrix3 *= WeightBOrbitals[MaxNbrFluxQuantaB - i] * WeightBOrbitals[MaxNbrFluxQuantaB - i];
 	}
      FullRightOverlapMatrix = TmpMatrix2 + TmpMatrix3;
    }

  
  SparseRealMatrix NormalizationMatrix(BMatrices[0].GetNbrRow(), BMatrices[0].GetNbrRow());
  NormalizationMatrix.SetToIdentity();
  int TmpPower = MaxNbrFluxQuantaA + 1 + MaxNbrFluxQuantaB + 1 - NbrFluxQuanta - 1;
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


  if (Manager.GetBoolean("orbital-es"))
    {
      double**** EntanglementSpectrum = new double***[PLevel + 1];
      int*** EntanglementSpectrumDimension = new int**[PLevel + 1];
      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
	{
	  EntanglementSpectrumDimension[CurrentPLevel] = new int*[NbrCFTSectors];
	  EntanglementSpectrum[CurrentPLevel] = new double**[NbrCFTSectors];	      	      
	  for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
	    {
	      int LocalMinQValue;
	      int LocalMaxQValue;
	      MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, LocalMinQValue, LocalMaxQValue);
	      if (LocalMinQValue <=  LocalMaxQValue)
		{
		  EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector] = new int[LocalMaxQValue - LocalMinQValue + 1];
		  EntanglementSpectrum[CurrentPLevel][CurrentCFTSector] = new double*[LocalMaxQValue - LocalMinQValue + 1];	      
		}
	      else
		{
		  EntanglementSpectrum[CurrentPLevel][CurrentCFTSector] = 0;
		}
	      for (int LocalQValue =  LocalMinQValue; LocalQValue <= LocalMaxQValue; ++LocalQValue)
		{
		  cout << "computing sector P=" << CurrentPLevel<< " CFT=" << CurrentCFTSector << " Q=" << LocalQValue << endl;
		  RealDiagonalMatrix TmpRhoADiag = FQHEMPSEvaluatePartialEntanglementSpectrum(MPSMatrix, FullLeftOverlapMatrix, FullRightOverlapMatrix, NormalizationMatrix,
											      CurrentPLevel, CurrentCFTSector, LocalQValue, 
											      CurrentPLevel, CurrentCFTSector, LocalQValue, Error);
		  EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = 0;
		  if (TmpRhoADiag.GetNbrRow() > 0)
		    {
		      for (int i = 0 ; i < TmpRhoADiag.GetNbrRow(); ++i)
			{
			  if (TmpRhoADiag[i] > 0.0)
			    EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]++;
			}
		      if (EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] > 0)
			{
			  EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = new double[EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]];
			  EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = 0;
			  for (int i = 0 ; i < TmpRhoADiag.GetNbrRow(); ++i)
			    {
			      if (TmpRhoADiag[i] > 0.0)
				{
				  EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]] = TmpRhoADiag[i];
				  EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]++;
				  TotalTraceThoA += TmpRhoADiag[i];
				}
			    }
			  SortArrayDownOrdering<double>(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue],
							EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]);
			}
		    }
		}
	    }
	}
      
      
      int GlobalMinQValue = 100000;
      int GlobalMaxQValue = 0;
      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
	{
	  for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
	    {
	      int MinQValue = 0;
	      int MaxQValue = 0;
	      MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
	      if (MinQValue < GlobalMinQValue)
		GlobalMinQValue = MinQValue;
	      if (MaxQValue > GlobalMaxQValue)
		GlobalMaxQValue = MaxQValue;		  
	    }
	}
      double EntanglementEntropy = 0.0;
      for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
	{
	  for (int LocalQValue =  GlobalMinQValue; LocalQValue <= GlobalMaxQValue; ++LocalQValue)
	    {
	      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
		{
		  int LocalMinQValue;
		  int LocalMaxQValue;
		  MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, LocalMinQValue, LocalMaxQValue);
		  if ((LocalQValue >= LocalMinQValue) && (LocalQValue <= LocalMaxQValue))
		    {
		      for (int i = 0; i < EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]; ++i)
			{
			  File << CurrentCFTSector  << " " << LocalQValue << " " 
			       << CurrentPLevel << " "
			       <<  (EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)  
			       <<  " " << (-log(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)) << endl;
			  EntanglementEntropy -= (log(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)
						  * EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA);
			}
		    }
		}
	    }
	}
      cout << "S_A=" << EntanglementEntropy << endl;
      File2 << (MaxNbrFluxQuantaA + 1) << " " << EntanglementEntropy << " 1" << endl;
      File.close();
      File2.close();
      cout << "Tr(rho_A)=" << TotalTraceThoA << endl;
      return 0;
    }
  if (Manager.GetBoolean("realspace-cut") == false)
    {
      double** EntanglementSpectrum = new double*[PLevel + 1];
      int* EntanglementSpectrumDimension = new int[PLevel + 1];
      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
	EntanglementSpectrumDimension[CurrentPLevel] = 0;
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
	  SparseRealMatrix LeftOverlapMatrix = MPSMatrix->ExtractBlock(FullLeftOverlapMatrix, CurrentPLevel, LeftQValue, CurrentPLevel, LeftQValue);
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

      double EntanglementEntropy = 0.0;
      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
	{
	  for (int i = 0; i < EntanglementSpectrumDimension[CurrentPLevel]; ++i)
	    {
	      File << Na << " " << CurrentPLevel << " " 
		   <<  (EntanglementSpectrum[CurrentPLevel][i] / EntanglementSpectrum[0][0])  
		   <<  " " << (-log(EntanglementSpectrum[CurrentPLevel][i] / TotalTraceThoA)) << endl;
	      //	     <<  (EntanglementSpectrum[CurrentPLevel][i] / EntanglementSpectrum[0][0] * (50.0 / 3.0))  
	      EntanglementEntropy -= (log(EntanglementSpectrum[CurrentPLevel][i] / TotalTraceThoA)
				      * EntanglementSpectrum[CurrentPLevel][i] / TotalTraceThoA);
	    }
	}
      cout << "S_A=" << EntanglementEntropy << endl;
      File2 << Na << " " << EntanglementEntropy << " 1" << endl;
      File.close();
      File2.close();
      cout << "Tr(rho_A)=" << TotalTraceThoA << endl;
    }
  else
    {
      double**** EntanglementSpectrum = new double***[PLevel + 1];
      int*** EntanglementSpectrumDimension = new int**[PLevel + 1];
      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
	{
	  EntanglementSpectrumDimension[CurrentPLevel] = new int*[NbrCFTSectors];
	  EntanglementSpectrum[CurrentPLevel] = new double**[NbrCFTSectors];	      	      
	  for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
	    {
	      int LocalMinQValue;
	      int LocalMaxQValue;
	      MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, LocalMinQValue, LocalMaxQValue);
	      if (LocalMinQValue <=  LocalMaxQValue)
		{
		  EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector] = new int[LocalMaxQValue - LocalMinQValue + 1];
		  EntanglementSpectrum[CurrentPLevel][CurrentCFTSector] = new double*[LocalMaxQValue - LocalMinQValue + 1];	      
		}
	      else
		{
		  EntanglementSpectrum[CurrentPLevel][CurrentCFTSector] = 0;
		}
	      for (int LocalQValue =  LocalMinQValue; LocalQValue <= LocalMaxQValue; ++LocalQValue)
		{
		  cout << "computing sector P=" << CurrentPLevel<< " CFT=" << CurrentCFTSector << " Q=" << LocalQValue << endl;
		  int RightLocalQValue = LocalQValue + (NbrFluxQuanta + 1) - (MaxNbrFluxQuantaB + 1) - (MaxNbrFluxQuantaA + 1);
		  EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = 0;
		  if ((RightLocalQValue <= LocalMaxQValue) && (RightLocalQValue  >= LocalMinQValue))
		    {
		      RealDiagonalMatrix TmpRhoADiag = FQHEMPSEvaluatePartialEntanglementSpectrum(MPSMatrix, FullLeftOverlapMatrix, FullRightOverlapMatrix, NormalizationMatrix,
												  CurrentPLevel, CurrentCFTSector, LocalQValue, 
												  CurrentPLevel, CurrentCFTSector, RightLocalQValue, Error);
		      if (TmpRhoADiag.GetNbrRow() > 0)
			{
			  for (int i = 0 ; i < TmpRhoADiag.GetNbrRow(); ++i)
			    {
			      if (TmpRhoADiag[i] > 0.0)
				EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]++;
			    }
			  if (EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] > 0)
			    {
			      EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = new double[EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]];
			      EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue] = 0;
			      for (int i = 0 ; i < TmpRhoADiag.GetNbrRow(); ++i)
				{
				  if (TmpRhoADiag[i] > 0.0)
				    {
				      EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]] = TmpRhoADiag[i];
				      EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]++;
				      TotalTraceThoA += TmpRhoADiag[i];
				    }
				}
			      SortArrayDownOrdering<double>(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue],
							    EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]);
			    }
			}
		    }
		}
	    }
	}
      
      
      int GlobalMinQValue = 100000;
      int GlobalMaxQValue = 0;
      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
	{
	  for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
	    {
	      int MinQValue = 0;
	      int MaxQValue = 0;
	      MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, MinQValue, MaxQValue);
	      if (MinQValue < GlobalMinQValue)
		GlobalMinQValue = MinQValue;
	      if (MaxQValue > GlobalMaxQValue)
		GlobalMaxQValue = MaxQValue;		  
	    }
	}
      double EntanglementEntropy = 0.0;
      for (int CurrentCFTSector = 0; CurrentCFTSector < NbrCFTSectors; ++CurrentCFTSector)
	{
	  for (int LocalQValue =  GlobalMinQValue; LocalQValue <= GlobalMaxQValue; ++LocalQValue)
	    {
	      for (int CurrentPLevel = 0; CurrentPLevel <= PLevel; ++CurrentPLevel)
		{
		  int LocalMinQValue;
		  int LocalMaxQValue;
		  MPSMatrix->GetChargeIndexRange(CurrentPLevel, CurrentCFTSector, LocalMinQValue, LocalMaxQValue);
		  if ((LocalQValue >= LocalMinQValue) && (LocalQValue <= LocalMaxQValue))
		    {
		      for (int i = 0; i < EntanglementSpectrumDimension[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue]; ++i)
			{
			  File << CurrentCFTSector  << " " << LocalQValue << " " 
			       << CurrentPLevel << " "
			       <<  (EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)  
			       <<  " " << (-log(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)) << endl;
			  EntanglementEntropy -= (log(EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA)
						  * EntanglementSpectrum[CurrentPLevel][CurrentCFTSector][LocalQValue - LocalMinQValue][i] / TotalTraceThoA);
			}
		    }
		}
	    }
	}
      cout << "S_A=" << EntanglementEntropy << endl;
      File2 << (MaxNbrFluxQuantaA + 1) << " " << EntanglementEntropy << " 1" << endl;
      File.close();
      File2.close();
      cout << "Tr(rho_A)=" << TotalTraceThoA << endl;
    }
  
  return 0;
}


// compute the entanglement spectrum from the overlap matrices of the left and right parts
// 
// mPSMatrix = pointer tothe MPS matrix
// leftOverlapMatrix= reference on the left overlap matrix
// rightOverlapMatrix= reference on the right overlap matrix
// leftPSector = P sector that has to be selected for the left part
// leftQSector = Q sector that has to be selected for the left part
// rightPSector = P sector that has to be selected for the right part
// rightQSector = Q sector that has to be selected for the right part
// eigenvalueError = relative error on the eigenvalues below which an eigenvalue is considered to be equal to zero

RealDiagonalMatrix FQHEMPSEvaluatePartialEntanglementSpectrum(AbstractFQHEMPSMatrix* mPSMatrix, SparseRealMatrix& leftOverlapMatrix, 
							      SparseRealMatrix& rightOverlapMatrix, int leftPSector, int leftQSector, 
							      int rightPSector, int rightQSector, double eigenvalueError)
{
  SparseRealMatrix RightOverlapMatrix = mPSMatrix->ExtractBlock(rightOverlapMatrix, rightPSector, rightQSector,
								rightPSector, rightQSector);
  SparseRealMatrix LeftOverlapMatrix = mPSMatrix->ExtractBlock(leftOverlapMatrix, leftPSector, leftQSector, leftPSector, leftQSector);
  RealDiagonalMatrix TmpRhoADiag;
  if ((LeftOverlapMatrix.GetNbrRow() > 0) && (RightOverlapMatrix.GetNbrRow() > 0) && 
      (LeftOverlapMatrix.ComputeNbrNonZeroMatrixElements() > 0l) && (RightOverlapMatrix.ComputeNbrNonZeroMatrixElements() > 0l))
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
      LocalLeftEigenvalueError *= eigenvalueError;
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
      LocalRightEigenvalueError *= eigenvalueError;
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
	      if (TmpRightDiag(i, i) > LocalRightEigenvalueError)
		{
		  TruncatedRightBasis[NbrZeroRightEigenvalues].Copy(TmpRightBasis[i]);
		  TruncatedRightBasis[NbrZeroRightEigenvalues] *= sqrt(TmpRightDiag(i, i));
		  ++NbrZeroRightEigenvalues;
		}
	    }
	  
	  RealMatrix TmpEntanglementMatrix = TruncatedLeftBasis.DuplicateAndTranspose();
	  TmpEntanglementMatrix.Multiply(TruncatedRightBasis);
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
    }
  return TmpRhoADiag;
}

// compute the entanglement spectrum from the overlap matrices of the left and right parts
// 
// mPSMatrix = pointer tothe MPS matrix
// leftOverlapMatrix= reference on the left overlap matrix
// rightOverlapMatrix= reference on the right overlap matrix
// normalizationMatrix = optional normalization matrix (for PES and RES)
// leftPSector = P sector that has to be selected for the left part
// leftCFTSector = CFT sector that has to be selected for the left part
// leftQSector = Q sector that has to be selected for the left part
// rightPSector = P sector that has to be selected for the right part
// rightCFTSector = CFT sector that has to be selected for the right part
// rightQSector = Q sector that has to be selected for the right part
// eigenvalueError = relative error on the eigenvalues below which an eigenvalue is considered to be equal to zero

RealDiagonalMatrix FQHEMPSEvaluatePartialEntanglementSpectrum(AbstractFQHEMPSMatrix* mPSMatrix, SparseRealMatrix& leftOverlapMatrix, SparseRealMatrix& rightOverlapMatrix, 
							      SparseRealMatrix& normalizationMatrix,
							      int leftPSector, int leftCFTSector, int leftQSector, 
							      int rightPSector, int rightCFTSector, int rightQSector, double eigenvalueError)
{
  SparseRealMatrix RightOverlapMatrix = mPSMatrix->ExtractBlock(rightOverlapMatrix, rightPSector, rightCFTSector, rightQSector,
								rightPSector, rightCFTSector, rightQSector);
  SparseRealMatrix LeftOverlapMatrix = mPSMatrix->ExtractBlock(leftOverlapMatrix, leftPSector, leftCFTSector, leftQSector, 
							       leftPSector, leftCFTSector, leftQSector);
  RealDiagonalMatrix TmpRhoADiag;
  if ((LeftOverlapMatrix.GetNbrRow() > 0) && (RightOverlapMatrix.GetNbrRow() > 0) && 
      (LeftOverlapMatrix.ComputeNbrNonZeroMatrixElements() > 0l) && (RightOverlapMatrix.ComputeNbrNonZeroMatrixElements() > 0l))
    {
      SparseRealMatrix NormalizationMatrix;
      if (normalizationMatrix.GetNbrRow() > 0)
	{
	  NormalizationMatrix = mPSMatrix->ExtractBlock(normalizationMatrix, rightPSector, rightCFTSector, rightQSector, leftPSector, leftCFTSector, leftQSector);
	  for (int i = 0; i < NormalizationMatrix.GetNbrRow(); ++i)
	    {
	      double Tmp;
	      NormalizationMatrix.GetMatrixElement(i, i, Tmp);
	      NormalizationMatrix.SetMatrixElement(i, i, 1.0 / Tmp);
	    }
	}
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
      LocalLeftEigenvalueError *= eigenvalueError;
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
      LocalRightEigenvalueError *= eigenvalueError;
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
	      if (TmpRightDiag(i, i) > LocalRightEigenvalueError)
		{
		  TruncatedRightBasis[NbrZeroRightEigenvalues].Copy(TmpRightBasis[i]);
		  TruncatedRightBasis[NbrZeroRightEigenvalues] *= sqrt(TmpRightDiag(i, i));
		  ++NbrZeroRightEigenvalues;
		}
	    }
	  
	  RealMatrix TmpEntanglementMatrix = TruncatedLeftBasis.DuplicateAndTranspose();
	  if (NormalizationMatrix.GetNbrRow() > 0)
	    {
	      RealMatrix TruncatedNormalizationMatrix(NormalizationMatrix);
	      TmpEntanglementMatrix.Multiply(TruncatedNormalizationMatrix);
	    }
	  TmpEntanglementMatrix.Multiply(TruncatedRightBasis);
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
    }
  return TmpRhoADiag;
}
