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

void CompactPrintMatrix(Matrix& TmpMatrix, int MatDim);
void ConvertSparseToDenseRealMatrix(SparseRealMatrix& Mat1, int MatDim, RealMatrix& Mat2, int Dim1, int Dim2);
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
       cout<<"Cylinder geometry, kappa= "<<kappa<<endl;
    }

  int LandauLevel = 0;

  AbstractFQHEMPSMatrix* MPSMatrix = MPSMatrixManager.GetMPSMatrices(NbrFluxQuanta); 
  if (Manager.GetBoolean("only-export"))
    {
      return 0;
    }

  SparseRealMatrix* SparseBMatrices = MPSMatrix->GetMatrices();

  cout << "B matrix size = " << SparseBMatrices[0].GetNbrRow() << "x" << SparseBMatrices[0].GetNbrColumn() << endl;
  
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


  int MatDim = SparseBMatrices[0].GetNbrRow();
  int LambdaMax = Manager.GetInteger("p-truncation");
  int LaughlinIndex = Manager.GetInteger("laughlin-index");
  int Nmax =  LambdaMax + (LaughlinIndex - 1)/2;
  int Nmin =  -LambdaMax - (LaughlinIndex - 1)/2;

  cout << "B matrix size = " << MatDim << "x" << MatDim << endl;

  SparseRealMatrix SparseFinalB0 (MatDim, MatDim);
  SparseRealMatrix SparseFinalB1 (MatDim, MatDim);
  SparseRealMatrix SparseFirstB0 (MatDim, MatDim);
  SparseRealMatrix SparseFirstB1 (MatDim, MatDim);

  for(int i = 0; i < MatDim; i++)
   {
    double Tmp;
    SparseBMatrices[0].GetMatrixElement(i, (Nmax-Nmin)/2, Tmp);
    SparseFinalB0.SetMatrixElement(i, (Nmax-Nmin)/2, Tmp);
    SparseBMatrices[1].GetMatrixElement(i, (Nmax-Nmin)/2, Tmp);
    SparseFinalB1.SetMatrixElement(i, (Nmax-Nmin)/2, Tmp);
    SparseBMatrices[0].GetMatrixElement((Nmax-Nmin)/2, i, Tmp);
    SparseFirstB0.SetMatrixElement((Nmax-Nmin)/2, i, Tmp);
    SparseBMatrices[1].GetMatrixElement((Nmax-Nmin)/2, i, Tmp);
    SparseFirstB1.SetMatrixElement((Nmax-Nmin)/2, i, Tmp);
   }



  //cout<<"  ********************* SparseFirstB0 ********************* " <<endl;
  //CompactPrintMatrix(SparseFirstB0, MatDim);

  //cout<<"  ********************* SparseFirstB1 ********************* " <<endl;
  //CompactPrintMatrix(SparseFirstB1, MatDim);

  //cout<<"  ********************* SparseFinalB0 ********************* " <<endl;
  //CompactPrintMatrix(SparseFinalB0, MatDim);

  //cout<<"  ********************* SparseFinalB1 ********************* " <<endl;
  //CompactPrintMatrix(SparseFinalB1, MatDim);


  cout<<"Done preparing B matrices and the vectors at 0 and Nphi orbital"<<endl;

  double CutOff = 1e-14;

  cout<<"Proceed to calculate overlap matrix (full space dimension, and it will be stored) "<<endl;

  SparseRealMatrix TmpM(MatDim, MatDim);

  SparseRealMatrix TmpM0(MatDim, MatDim);
  TmpM0 = SparseFirstB0.Transpose();
  TmpM0.Multiply(SparseFirstB0);
  TmpM = TmpM + TmpM0;
  SparseRealMatrix TmpM1(MatDim, MatDim);
  if (SparseFirstB1.ComputeNbrNonZeroMatrixElements() != 0)
    TmpM1 = SparseFirstB1.Transpose();
  TmpM1.Multiply(SparseFirstB1);
  TmpM = TmpM + TmpM1;

  for (int i = 1; i < EntCut; i++)
    {
      TmpM0 = SparseBMatrices[0].Transpose();
      TmpM0.Multiply(TmpM);
      TmpM0.Multiply(SparseBMatrices[0]);

      TmpM1 = SparseBMatrices[1].Transpose();
      TmpM1.Multiply(TmpM);
      TmpM1.Multiply(SparseBMatrices[1]);

      TmpM = (TmpM0 + TmpM1);
    }

  //cout<<"-------------------- Overlap matrix computed ------------------ "<<endl;
  //CompactPrintMatrix(TmpM, MatDim);


  cout<<"Compute density matrix in nonorthogonal basis (full space dimension, will be stored)"<<endl;

  SparseRealMatrix rhoM  (MatDim, MatDim);
  SparseRealMatrix rhoM0 (MatDim, MatDim);
  SparseRealMatrix rhoM1 (MatDim, MatDim);

  rhoM0.Copy(SparseFinalB0);
  rhoM0.Multiply(SparseFinalB0.Transpose());

  if (SparseFirstB1.ComputeNbrNonZeroMatrixElements() != 0)
   {
     rhoM1.Copy(SparseFinalB1);
     rhoM1.Multiply(SparseFinalB1.Transpose());
   }

  rhoM.Copy(rhoM0);
  rhoM = rhoM + rhoM1;

  for (int i = (NbrFluxQuanta - 1); i >= EntCut; i--)
    {
      rhoM0.Copy(SparseBMatrices[0]);
      rhoM0.Multiply(rhoM);
      rhoM0.Multiply(SparseBMatrices[0].Transpose());

      rhoM1.Copy(SparseBMatrices[1]);
      rhoM1.Multiply(rhoM);
      rhoM1.Multiply(SparseBMatrices[1].Transpose());

      rhoM = (rhoM0 + rhoM1);
    }

  //cout<<"-------------------- rho in nonorthogonal basis computed ------------------ "<<endl;
  //CompactPrintMatrix(rhoM, MatDim); 
  
  //Free up some space that is no longer needed (this needs to be done in a cleaner way)

  delete[] SparseBMatrices; 
  TmpM0.ResizeAndClean(MatDim, MatDim);
  TmpM1.ResizeAndClean(MatDim, MatDim);
  rhoM0.ResizeAndClean(MatDim, MatDim);
  rhoM1.ResizeAndClean(MatDim, MatDim);
  SparseFirstB0.ResizeAndClean(MatDim, MatDim);
  SparseFirstB1.ResizeAndClean(MatDim, MatDim);
  SparseFinalB0.ResizeAndClean(MatDim, MatDim);
  SparseFinalB0.ResizeAndClean(MatDim, MatDim);

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
    for (int MomentumSector = 0; MomentumSector <= LambdaMax; MomentumSector++)
    {

      //cout<<"##########################################################################"<<endl;

      SparseRealMatrix TmpMBlock = MPSMatrix->ExtractBlock(TmpM, MomentumSector, NSector, MomentumSector, NSector);
 
      int TmpSectorDim = TmpMBlock.GetNbrRow();

      RealMatrix TmpMDense (TmpSectorDim, TmpSectorDim, true);
      ConvertSparseToDenseRealMatrix(TmpMBlock, TmpSectorDim, TmpMDense, 0, TmpSectorDim);

      RealSymmetricMatrix HRep ((Matrix&)TmpMDense);
      RealDiagonalMatrix TmpDiag (TmpSectorDim);
      RealMatrix Q(TmpSectorDim, TmpSectorDim);
      HRep.LapackDiagonalize(TmpDiag, Q);

      int NonZeroEig = 0;
      for (int j = TmpSectorDim - 1; j >= 0; --j)
       {
         //cout<<j<<" "<<TmpDiag[j]<<" ; ";
         if (fabs(TmpDiag[j])>CutOff)
          {  
            NonZeroEig++;
            Q[j] /= sqrt(fabs(TmpDiag[j]));
       	    //for (int i = 0; i < TmpSectorDim; i++)
            //  if (Norm(Q[j][i]) > CutOff) cout<<"i= "<<i<<" "<<Q[j][i]<<"; ";
            //cout<<endl;
          }  
       }
      //cout<<endl;
      //cout<<"Nbr of nonzero vectors= "<<NonZeroEig<<" out of "<<TmpSectorDim<<endl;

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

      //cout<<"------------------ Done with new rho -----------------------"<<endl;
      //CompactPrintMatrix(NewrhoM, NonZeroEig);

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

      } //NonZeroEig
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


void CompactPrintMatrix(Matrix& TmpMatrix, int MatDim)
{
  for (int i = 0; i < MatDim; i++)
    for (int j = 0; j < MatDim; j++)
      {
         double Tmp;
         TmpMatrix.GetMatrixElement(i, j, Tmp);
         if (fabs(Tmp) != 0.0)
           cout<<"i= "<<i<<" j= "<<j<<" "<<Tmp<<endl;
      }
}

void ConvertSparseToDenseRealMatrix(SparseRealMatrix& Mat1, int MatDim, RealMatrix& Mat2, int Dim1, int Dim2)
{
  Mat2.ClearMatrix();
  for (int i = Dim1; i < Dim2; i++)
    for (int j = Dim1; j < Dim2; j++)
      {
         double Tmp;
         Mat1.GetMatrixElement(i, j, Tmp);
         Mat2.SetMatrixElement(i - Dim1, j - Dim1, Tmp);
      }
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

/* Old working version with complex matrices

#include "HilbertSpace/FermionOnSpherePTruncated.h"
#include "HilbertSpace/FermionOnSpherePTruncatedLong.h"
#include "HilbertSpace/BosonOnDiskShort.h"
#include "HilbertSpace/FermionOnSphereMPSWrapper.h"
#include "HilbertSpace/FermionOnCylinderMPSWrapper.h"

#include "MathTools/ClebschGordanCoefficients.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"
#include "MathTools/FactorialCoefficient.h"

#include "Operator/ParticleOnSphereDensityOperator.h"
#include "Operator/ParticleOnSphereDensityDensityOperator.h"

#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"
#include "FunctionBasis/ParticleOnCylinderFunctionBasis.h"

#include "Vector/Vector.h"
#include "Vector/ComplexVector.h"
#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Matrix/SparseComplexMatrix.h"

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


void CreateLaughlinBMatrices (int laughlinIndex, ComplexMatrix* bMatrices, BosonOnDiskShort** u1BosonBasis, int pLevel, bool cylinderFlag, double kappa);
Complex CreateLaughlinAMatrixElement (int laughlinIndex, unsigned long* partition1, unsigned long* partition2, int p1Level, int p2Level, int nValue, FactorialCoefficient& coef);
Complex OverlapEig(ComplexVector& EigV, SparseComplexMatrix& TmpM, int TmpSectorDim, int StartIndex, int i);
void CompactPrintMatrix(Matrix& TmpMatrix, int MatDim);
void ConvertSparseToDenseComplexMatrix(SparseComplexMatrix& Mat1, int MatDim, ComplexMatrix& Mat2, int Dim1, int Dim2);

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereMPSCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "file that describes the root configuration");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "p-truncation", "truncation level", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "laughlin-index", "index of the Laughlin state to generate", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "la", "number of orbitals in subsystem A", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "na", "number of particles in subsystem A", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics");
  (*SystemGroup) += new BooleanOption ('\n', "cylinder", "evaluate density on the cylinder");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "aspect-ratio", "aspect ratio of the cylinder", 1);
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

  cout<<"NbrParticles= "<<NbrParticles<<" NbrFluxQuanta= "<<NbrFluxQuanta<<endl;
 
  bool CylinderFlag = Manager.GetBoolean("cylinder");
  double AspectRatio = Manager.GetDouble("aspect-ratio");
  double kappa = 0.0;
  if (CylinderFlag)
    {
       kappa = (2.0 * M_PI)/sqrt(2.0 * M_PI * (NbrFluxQuanta + 1) * AspectRatio);
       cout<<"Cylinder geometry, kappa= "<<kappa<<endl;
    }


  int EntCut = Manager.GetInteger("la");
  int Na = Manager.GetInteger("na");

  int LandauLevel = 0;

  int LaughlinIndex = Manager.GetInteger("laughlin-index");


  ofstream File;
  File.precision(14);
  char* TmpFileName = new char [512];
  sprintf(TmpFileName, "fermions_laughlin%ld_plevel_%ld_n_%d_2s_%d_lz_%d_la_%d_na_%d.ent", Manager.GetInteger("laughlin-index"), Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz, EntCut, Na);
  File.open(TmpFileName, ios::binary | ios::out);     

  int NbrBMatrices = 2;
  ComplexMatrix* BMatrices = new ComplexMatrix[NbrBMatrices];
  SparseComplexMatrix* SparseBMatrices = new SparseComplexMatrix[NbrBMatrices];
  BosonOnDiskShort** U1BosonBasis = new BosonOnDiskShort* [Manager.GetInteger("p-truncation") + 1];
  for (int i = 0; i <= Manager.GetInteger("p-truncation"); ++i)
    {
      U1BosonBasis[i] = new BosonOnDiskShort(i, i, Manager.GetInteger("p-truncation") + 1);
    }
  CreateLaughlinBMatrices (LaughlinIndex, BMatrices, U1BosonBasis, Manager.GetInteger("p-truncation"), CylinderFlag, kappa);

  for (int i = 0; i < NbrBMatrices; ++i)
    {
      SparseBMatrices[i] = BMatrices[i];
    }

  delete[] BMatrices;


  int MatDim = SparseBMatrices[0].GetNbrRow();
  int LambdaMax = Manager.GetInteger("p-truncation");
 
  int Nmax =  LambdaMax + (LaughlinIndex - 1)/2;
  int Nmin =  -LambdaMax - (LaughlinIndex - 1)/2;

  cout << "B matrix size = " << MatDim << "x" << MatDim << endl;

//  cout<<"  ********************* B0 ********************* " <<endl;
//  CompactPrintMatrix(SparseBMatrices[0], MatDim);

//  cout<<"  ********************* B1 ********************* " <<endl;
//  CompactPrintMatrix(SparseBMatrices[1], MatDim);

  ComplexMatrix FinalB0(MatDim, MatDim, true);
  ComplexMatrix FinalB1(MatDim, MatDim, true);
  ComplexMatrix FirstB0(MatDim, MatDim, true);
  ComplexMatrix FirstB1(MatDim, MatDim, true);


  for(int i = 0; i < MatDim; i++)
   {
    Complex Tmp;
    SparseBMatrices[0].GetMatrixElement(i, (Nmax-Nmin)/2, Tmp);
    FinalB0.SetMatrixElement(i, (Nmax-Nmin)/2, Tmp);
    SparseBMatrices[1].GetMatrixElement(i, (Nmax-Nmin)/2, Tmp);
    FinalB1.SetMatrixElement(i, (Nmax-Nmin)/2, Tmp);
    SparseBMatrices[0].GetMatrixElement((Nmax-Nmin)/2, i, Tmp);
    FirstB0.SetMatrixElement((Nmax-Nmin)/2, i, Tmp);
    SparseBMatrices[1].GetMatrixElement((Nmax-Nmin)/2, i, Tmp);
    FirstB1.SetMatrixElement((Nmax-Nmin)/2, i, Tmp);
   }

  SparseComplexMatrix SparseFinalB0 = FinalB0;
  SparseComplexMatrix SparseFinalB1 = FinalB1;
  SparseComplexMatrix SparseFirstB0 = FirstB0;
  SparseComplexMatrix SparseFirstB1 = FirstB1;


//  cout<<"  ********************* SparseFirstB0 ********************* " <<endl;
//  CompactPrintMatrix(SparseFirstB0, MatDim);

//  cout<<"  ********************* SparseFirstB1 ********************* " <<endl;
//  CompactPrintMatrix(SparseFirstB1, MatDim);

//  cout<<"  ********************* SparseFinalB0 ********************* " <<endl;
//  CompactPrintMatrix(SparseFinalB0, MatDim);

//  cout<<"  ********************* SparseFinalB1 ********************* " <<endl;
//  CompactPrintMatrix(SparseFinalB1, MatDim);



  cout<<"Done preparing B matrices and the vectors at 0 and Nphi orbital"<<endl;

  double CutOff = 1e-14;

  SparseComplexMatrix TmpM(MatDim, MatDim, 0, true);

  SparseComplexMatrix TmpM0(MatDim, MatDim, 0, true);
  TmpM0 = SparseFirstB0.HermitianTranspose();
  TmpM0.Multiply(SparseFirstB0);
  TmpM = TmpM + TmpM0;
  SparseComplexMatrix TmpM1(MatDim, MatDim, 0, true);
  if (SparseFirstB1.ComputeNbrNonZeroMatrixElements() != 0)
    TmpM1 = SparseFirstB1.HermitianTranspose();
  TmpM1.Multiply(SparseFirstB1);
  TmpM = TmpM + TmpM1;

  for (int i = 1; i < EntCut; i++)
    {
      TmpM0 = SparseBMatrices[0].HermitianTranspose();
      TmpM0.Multiply(TmpM);
      TmpM0.Multiply(SparseBMatrices[0]);

      TmpM1 = SparseBMatrices[1].HermitianTranspose();
      TmpM1.Multiply(TmpM);
      TmpM1.Multiply(SparseBMatrices[1]);

      TmpM = (TmpM0 + TmpM1);
    }

  cout<<"-------------------- Overlap matrix computed ------------------ "<<endl;
  CompactPrintMatrix(TmpM, MatDim);

  cout<<"Compute rho in nonorthogonal basis"<<endl;

  SparseComplexMatrix rhoM(MatDim, MatDim, 0, true);
  SparseComplexMatrix rhoM0 (MatDim, MatDim, 0, true);
  rhoM0 = SparseFinalB0;
  rhoM0.Multiply(SparseFinalB0.HermitianTranspose());
  rhoM = rhoM + rhoM0;
  SparseComplexMatrix rhoM1 (MatDim, MatDim, 0, true);
  if (SparseFirstB1.ComputeNbrNonZeroMatrixElements() != 0)
   {
     rhoM1 = SparseFinalB1;
     rhoM1.Multiply(SparseFinalB1.HermitianTranspose());
   }
  rhoM = rhoM + rhoM1;

  for (int i = (NbrFluxQuanta - 1); i >= EntCut; i--)
    {
      rhoM0.Copy(SparseBMatrices[0]);
      rhoM0.Multiply(rhoM);
      rhoM0.Multiply(SparseBMatrices[0].HermitianTranspose());

      rhoM1.Copy(SparseBMatrices[1]);
      rhoM1.Multiply(rhoM);
      rhoM1.Multiply(SparseBMatrices[1].HermitianTranspose());

      rhoM = (rhoM0 + rhoM1);
    }

  cout<<"-------------------- rho in nonorthogonal basis computed ------------------ "<<endl;
  //CompactPrintMatrix(rhoM, MatDim); 
  
  int* RowIndices = new int[rhoM.GetNbrMatrixElements()]; 
  rhoM.GetRowIndices(RowIndices);
      

   delete[] SparseBMatrices; 
   TmpM0.ResizeAndClean(MatDim, MatDim);
   TmpM1.ResizeAndClean(MatDim, MatDim);
   rhoM0.ResizeAndClean(MatDim, MatDim);
   rhoM1.ResizeAndClean(MatDim, MatDim);
   SparseFirstB0.ResizeAndClean(MatDim, MatDim);
   SparseFirstB1.ResizeAndClean(MatDim, MatDim);
   SparseFinalB0.ResizeAndClean(MatDim, MatDim);
   SparseFinalB0.ResizeAndClean(MatDim, MatDim);


  int* StartingIndexPerPLevel = new int [LambdaMax + 1];
  int* NbrIndicesPerPLevel = new int [LambdaMax + 1];
  StartingIndexPerPLevel[0] = 0;
  int NbrNValue = ((2 * LambdaMax) + LaughlinIndex);
  int NValueShift = NbrNValue - 1;
  NbrIndicesPerPLevel[0] = U1BosonBasis[0]->GetHilbertSpaceDimension() * NbrNValue;
  for (int i = 1; i <= LambdaMax; ++i)
    {
      StartingIndexPerPLevel[i] = StartingIndexPerPLevel[i - 1] + NbrIndicesPerPLevel[i - 1];
      NbrIndicesPerPLevel[i] = U1BosonBasis[i]->GetHilbertSpaceDimension()  * NbrNValue;
    }



//     cout<<"++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;
//      unsigned long* TmpPartition = new unsigned long [LambdaMax + 2];
//      for (int p = 0; p <= LambdaMax; ++p)
//      {
//      BosonOnDiskShort* TmpSpace = U1BosonBasis[p];
//      for (int n=0; n<NbrNValue; ++n)
//      for (int i = 0; i < U1BosonBasis[p]->GetHilbertSpaceDimension(); ++i)
//       {
//         TmpSpace->GetOccupationNumber(i, TmpPartition);
//         cout<<"|P|= "<<p<<" P="; 
//         for(int x=0; x <(LambdaMax+2); ++x)
//            cout<<TmpPartition[x]<<" ";
//         cout<<"N = "<<(n-NValueShift/2)<<" ; ";    
//         cout<<endl;  
//       }
//      }
//      delete[] TmpPartition;
//     cout<<"++++++++++++++++++++++++++++++++++++++++++++++++"<<endl;




  double TraceRho = 0.0;
  double* RhoEigenvalues = new double [MatDim];
  int* RhoNSector = new int [MatDim];
  int* RhoPSector = new int [MatDim];


  for (int i =0 ; i < MatDim; i++)
   {
     RhoEigenvalues[i] = 0.0; 
     RhoNSector[i] = 0; 
     RhoPSector[i]=0;
   }

  int eigenvaluecounter = 0;
  for (int NSector = 0; NSector < NbrNValue; NSector++)
    for (int MomentumSector = 0; MomentumSector <= LambdaMax; MomentumSector++)
    {

      int StartIndex = StartingIndexPerPLevel[MomentumSector] + NSector * U1BosonBasis[MomentumSector]->GetHilbertSpaceDimension();
      int EndIndex = StartIndex + U1BosonBasis[MomentumSector]->GetHilbertSpaceDimension();

      int TmpSectorDim = U1BosonBasis[MomentumSector]->GetHilbertSpaceDimension();


//      cout<<"testing indices"<<endl;
//
//      unsigned long* Partition = new unsigned long [LambdaMax + 2];
//      BosonOnDiskShort* TmpSpace = U1BosonBasis[MomentumSector];
//      for (int i = 0; i < TmpSectorDim; ++i)
//       {
//         TmpSpace->GetOccupationNumber(i, Partition);
//         cout<<"|P|= "<<MomentumSector<<" P="; 
//         for(int x=0; x <(LambdaMax+2); ++x)
//            cout<<Partition[x]<<" ";
//         cout<<"N = "<<NSector<<" ; ";    
//         cout<<endl;  
//       }
//      delete[] Partition;


      ComplexMatrix TmpMDense (TmpSectorDim, TmpSectorDim, true);

      ConvertSparseToDenseComplexMatrix(TmpM, MatDim, TmpMDense, StartIndex, EndIndex);

      HermitianMatrix HRep = TmpMDense;

      //cout<<TmpMDense<<endl;

      RealDiagonalMatrix TmpDiag (TmpSectorDim);
      ComplexMatrix Q(TmpSectorDim, TmpSectorDim);
      HRep.LapackDiagonalize(TmpDiag, Q);


      //cout<<"TmpM "<<endl;
      int NonZeroEig = 0;
      for (int j = TmpSectorDim - 1; j >= 0; --j)
       {
         //cout<<j<<" "<<TmpDiag[j]<<" ; ";
         if (fabs(TmpDiag[j])>CutOff)
          {  
            NonZeroEig++;
            Q[j] /= sqrt(fabs(TmpDiag[j]));
       	    //for (int i = 0; i < TmpSectorDim; i++)
            //  if (Norm(Q[j][i]) > CutOff) cout<<"i= "<<i<<" "<<Q[j][i]<<"; ";
            //cout<<endl;
          }  
       }
      //cout<<endl;
      cout<<"Nbr of nonzero vectors= "<<NonZeroEig<<" out of "<<TmpSectorDim<<endl;

      if (NonZeroEig > 0)
      {

      //cout<<"-------------------- Start computing rho in the new basis --------------------"<<endl;

      ComplexMatrix NewrhoM(NonZeroEig, NonZeroEig, true);

      int rowIndex, columnIndex;

      int offset = TmpSectorDim - NonZeroEig;
      for (int i = TmpSectorDim - 1; i >= offset; i--)
        for (int j = TmpSectorDim - 1; j >= offset; j--)
          {
            Complex TmpEl(0.0, 0.0);
            for (int ii = 0; ii < rhoM.GetNbrMatrixElements(); ++ii)
              {
               columnIndex = rhoM.GetColumnIndex(ii);
               rowIndex = RowIndices[ii];
               TmpEl += Conj(rhoM.GetMatrixElement(ii)) * OverlapEig(Q[j], TmpM, TmpSectorDim, StartIndex, columnIndex) * Conj(OverlapEig(Q[i], TmpM, TmpSectorDim, StartIndex, rowIndex));
            }

           NewrhoM.SetMatrixElement(i-offset, j-offset, TmpEl); 
         }

      //cout<<"------------------ Done with new rho -----------------------"<<endl;
      //CompactPrintMatrix(NewrhoM, NonZeroEig);

      HermitianMatrix HRepRho = NewrhoM ;///NewrhoM.Tr();

      RealDiagonalMatrix TmpDiagRho (NonZeroEig);
      ComplexMatrix QRho(NonZeroEig, NonZeroEig);
      HRepRho.LapackDiagonalize(TmpDiagRho, QRho);

      cout<<"--------------------------sector P = "<<MomentumSector<<" N = "<<((EntCut - (NSector - (2 * LambdaMax + LaughlinIndex - 1)/2)))/LaughlinIndex << "--------------------------"<<endl;

      cout.precision(14); 

     double Sum = 0.0;
     for (int j = 0; j < NonZeroEig; ++j)
      {
        
//        if (fabs(TmpDiagRho[j]) > CutOff)
//         {
//           bool found = false;
//        for(int i = 0; i < NonZeroEig; ++i)
//          if (Norm(QRho[j][i]) > sqrt(CutOff)) 
//           {
//             for (int k = 0; k < MatDim; ++k)
//               if (Norm(Q[i + offset][k]) > sqrt(CutOff))
//                 if (NArray[k] == Na)
//                  {
//                    cout<<PArray[k]<<" "<<NArray[k]<<"     ";
//                    File<<PArray[k]<<" "<<NArray[k]<<"     ";
//                    k = MatDim;
//                    i = NonZeroEig;
//                    found = true;
//                  } 
//           }
//       if (found)
//        {
//         
         //cout<<TmpDiagRho[j]<<endl;
         TraceRho += TmpDiagRho[j];
         RhoEigenvalues[eigenvaluecounter]=TmpDiagRho[j];
         RhoNSector[eigenvaluecounter]=NSector;
         RhoPSector[eigenvaluecounter]=MomentumSector; 
         eigenvaluecounter++;
         
         //File<<TmpDiagRho[j]<<endl;
        //}
      //}
   }

    } //NonZeroEig

     // cout<<"--------------------------completed sector P= "<<MomentumSector<<"---------------------"<<endl;
    }
 
  cout<<"Trace rho = "<<TraceRho<<endl;

  for (int i=0; i<MatDim; ++i)
    if ((fabs(RhoEigenvalues[i]) > CutOff) && ((((EntCut - (RhoNSector[i] - (2 * LambdaMax + LaughlinIndex - 1)/2)))/LaughlinIndex) == Na))
      {
        cout<<"P= "<<RhoPSector[i]<<" N= "<<RhoNSector[i]<<" "<<RhoEigenvalues[i]/TraceRho<<endl;  
        File<<RhoPSector[i]<<" "<<RhoNSector[i]<<" "<<RhoEigenvalues[i]/TraceRho<<endl;
      }

 delete[] StartingIndexPerPLevel;
 delete[] NbrIndicesPerPLevel;
 delete[] RhoEigenvalues;
 delete[] RhoPSector;
 delete[] RhoNSector;


//  cout<<"Testing TmpM"<<endl;
//  for(int i =0; i<MatDim; ++i)
//    for(int j =0; j<MatDim; ++j)
//      {
//         Complex Tmp;
//         TmpM.GetMatrixElement(i,j,Tmp);
//         if (Norm(Tmp) > CutOff)
//           {
//             if ((PArray[i] != PArray[j]) || (NArray[i] != NArray[j]))
//                cout<<i<<" "<<PArray[i]<<" "<<NArray[i]<<"; "<<j<<PArray[j]<<" "<<NArray[j]<<" "<<Tmp<<endl;
//           }
//      }

//  cout<<"PArray, NArray "; for(int i=0; i<MatDim; ++i) cout<<i<<" "<<PArray[i]<<" "<<NArray[i]<<endl;

//  cout<<"Analyze TmpM "<<endl;

//  for (int j = (MatDim - NonZeroEig); j < MatDim; ++j)
//   {
//     if (fabs(TmpDiag[j])>CutOff)
//      {  
//       cout<<TmpDiag[j]<<" Eigenvector: ";
//       for (int i = 0; i < MatDim; ++i)
//         if (Norm(Q[j][i]) > sqrt(CutOff)) cout<<Q[j][i]<<" "<<PArray[i]<<" "<<NArray[i]<<" ; ";
//       cout<<endl;
//      }  
//   }
//  cout<<endl;


  File.close();

  return 0;
}

void CreateLaughlinBMatrices (int laughlinIndex, ComplexMatrix* bMatrices, BosonOnDiskShort** u1BosonBasis, int pLevel, bool cylinderFlag, double kappa)
{
  int* StartingIndexPerPLevel = new int [pLevel + 1];
  int* NbrIndicesPerPLevel = new int [pLevel + 1];
  StartingIndexPerPLevel[0] = 0;
  int NbrNValue = ((2 * pLevel) + laughlinIndex);
  int NValueShift = NbrNValue - 1;
  NbrIndicesPerPLevel[0] = u1BosonBasis[0]->GetHilbertSpaceDimension() * NbrNValue;
  for (int i = 1; i <= pLevel; ++i)
    {
      StartingIndexPerPLevel[i] = StartingIndexPerPLevel[i - 1] + NbrIndicesPerPLevel[i - 1];
      NbrIndicesPerPLevel[i] = u1BosonBasis[i]->GetHilbertSpaceDimension()  * NbrNValue;
    }
  int MatrixSize = NbrIndicesPerPLevel[pLevel] + StartingIndexPerPLevel[pLevel];

  bMatrices[0] = ComplexMatrix(MatrixSize, MatrixSize, true);

  for (int i = 0; i <= pLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace = u1BosonBasis[i];
      for (int j = 1; j < NbrNValue; ++j)
	{
	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
	    {
              int N1 = (j - NValueShift/2);
	      Complex Tmp (1.0, 0.0);
              if (cylinderFlag)
                Tmp *= exp(-kappa*kappa*(i + (N1 - 1) * (N1 - 1)/(4.0 * laughlinIndex)+ (N1 * N1)/(4.0 * laughlinIndex)));
	      bMatrices[0].SetMatrixElement(StartingIndexPerPLevel[i] + k + (j - 1) * TmpSpace->GetHilbertSpaceDimension(), StartingIndexPerPLevel[i] + k  + j * TmpSpace->GetHilbertSpaceDimension(), Tmp);
	    }
	}
    }

  bMatrices[1] = ComplexMatrix(MatrixSize, MatrixSize, true);
  unsigned long* Partition1 = new unsigned long [pLevel + 2];
  unsigned long* Partition2 = new unsigned long [pLevel + 2];
  FactorialCoefficient Coef;


//  cout<<"testing indices"<<endl;
//  for (int i = 0; i <= pLevel; ++i)
//    {
//      BosonOnDiskShort* TmpSpace = u1BosonBasis[i];
//      for (int j = 0; j < NbrNValue; ++j)
//	{
//	  for (int k = 0; k < TmpSpace->GetHilbertSpaceDimension(); ++k)
//	    {
//              int N1 = (j - NValueShift/2);
//	      TmpSpace->GetOccupationNumber(k, Partition1);
//
//              cout<<"|P|= "<<i<<" P1="; 
//              for(int x=0; x <(pLevel+2); ++x)
//                cout<<Partition1[x]<<" ";
//              cout<<"N = "<<N1<<" ; ";    
//              cout<<endl;  
//	    }
//	}
//    }


  for (int i = 0; i <= pLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace1 = u1BosonBasis[i];
      for (int j = 0; j <= pLevel; ++j)
	{
	  BosonOnDiskShort* TmpSpace2 = u1BosonBasis[j];
	  int N2 = (2 * (j - i) - laughlinIndex + 1 + NValueShift) / 2;
	  int N1 = N2 + (laughlinIndex - 1);
	  for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
	    {
	      TmpSpace1->GetOccupationNumber(k1, Partition1);
	      for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
		{
		  TmpSpace2->GetOccupationNumber(k2, Partition2);
 	          Complex Tmp = CreateLaughlinAMatrixElement(laughlinIndex, Partition1, Partition2, i, j, -(N1 + N2 - NValueShift)/2, Coef);
		  if (cylinderFlag)
		     Tmp *= exp(-kappa*kappa*(0.5 * i + 0.5 * j + pow(N1 - NValueShift/2,2.0)/(4.0 * laughlinIndex) + pow(N2 - NValueShift/2,2.0)/(4.0 * laughlinIndex)));
		  bMatrices[1].SetMatrixElement(StartingIndexPerPLevel[i] + k1 + N1 * TmpSpace1->GetHilbertSpaceDimension(), StartingIndexPerPLevel[j] + k2 + N2 * TmpSpace2->GetHilbertSpaceDimension(), Tmp);


//                  cout<<(StartingIndexPerPLevel[i] + (k1 + N1 * TmpSpace1->GetHilbertSpaceDimension()))<<" "<<(StartingIndexPerPLevel[j] + (k2 + N2 * TmpSpace2->GetHilbertSpaceDimension()))<<" : ";
//                  cout<<"|P1|= "<<i<<" P1="; 
//                  for(int x=0; x <(pLevel+2); ++x)
//                    cout<<Partition1[x]<<" ";
//                  cout<<"N1 = "<<(N1-NValueShift/2)<<" ; ";    
//                  cout<<"|P2|= "<<j<<" P2="; 
//                  for(int x=0; x <(pLevel+2); ++x)
//                    cout<<Partition2[x]<<" ";
//                  cout<<"N2 = "<<(N2-NValueShift/2)<<" ; ";
//                  cout<<Tmp<<endl;    

		}
	    }
	}
    }

  delete[] Partition1;
  delete[] Partition2;
}

Complex CreateLaughlinAMatrixElement (int laughlinIndex, unsigned long* partition1, unsigned long* partition2, int p1Level, int p2Level, int nValue, FactorialCoefficient& coef)
{
  Complex Tmp = 1.0;
  if (nValue != (p1Level - p2Level))
    {
      Tmp = 0.0;
      return Tmp;
    }
  int PMax = p1Level;
  if (p2Level > p1Level)
    PMax = p2Level;
  for (int i = 1; i <= PMax; ++i)
    {
      Complex Tmp2 = 0.0;
//      cout << partition1[i] << " " << partition2[i] << " " << i << " " << PMax << endl;
      for (int j = 0; j <= partition1[i]; ++j)
	{
	  int k = partition2[i] + j - partition1[i];
	  if ((k >= 0) && (k <= partition2[i]))
	    {
	      int Sum = k + j;
	      coef.SetToOne();
	      coef.PartialFactorialMultiply(partition1[i] - j + 1, partition1[i]);
	      coef.PartialFactorialMultiply(partition2[i] - k + 1, partition2[i]);
	      coef.FactorialDivide(j);
	      coef.FactorialDivide(k);
	      coef.FactorialDivide(j);
	      coef.FactorialDivide(k);
	      coef.PowerNMultiply(laughlinIndex, Sum);
	      coef.PowerNDivide(i, Sum);
//	      cout << "Sum=" << Sum << endl;
	      switch  (Sum & 0x3)
		{
		case 0:
		  Tmp2.Re += sqrt(coef.GetNumericalValue());
		  break;
		case 1:
		  Tmp2.Im += sqrt(coef.GetNumericalValue());
		  break;
		case 2:
		  Tmp2.Re -= sqrt(coef.GetNumericalValue());
		  break;
		case 3:
		  Tmp2.Im -= sqrt(coef.GetNumericalValue());
		  break;
		}
	    }
	}
//      cout << Tmp2 << endl;
      Tmp *= Tmp2;
    }
  return Tmp;
}

Complex OverlapEig(ComplexVector& EigV, SparseComplexMatrix& TmpM, int TmpSectorDim, int StartIndex, int i)
{
  Complex Tmp;
  Complex ov(0.0, 0.0);
  for(int ii = 0; ii < TmpSectorDim; ii++)
   {
      TmpM.GetMatrixElement(ii + StartIndex, i, Tmp); 
      ov += EigV[ii] * Tmp;
   } 
 return ov;
}

void CompactPrintMatrix(Matrix& TmpMatrix, int MatDim)
{
  for (int i = 0; i < MatDim; i++)
    for (int j = 0; j < MatDim; j++)
      {
         Complex Tmp;
         TmpMatrix.GetMatrixElement(i, j, Tmp);
         if (Norm(Tmp) != 0.0)
           cout<<"i= "<<i<<" j= "<<j<<" "<<Tmp<<endl;
      }
}

void ConvertSparseToDenseComplexMatrix(SparseComplexMatrix& Mat1, int MatDim, ComplexMatrix& Mat2, int Dim1, int Dim2)
{
  Mat2.ClearMatrix();
  for (int i = Dim1; i < Dim2; i++)
    for (int j = Dim1; j < Dim2; j++)
      {
         Complex Tmp;
         Mat1.GetMatrixElement(i, j, Tmp);
         Mat2.SetMatrixElement(i - Dim1, j - Dim1, Tmp);
      }
}

*/
