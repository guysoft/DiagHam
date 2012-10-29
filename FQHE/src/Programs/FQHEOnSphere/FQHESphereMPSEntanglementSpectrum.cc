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
Complex OverlapEig(ComplexVector& EigV, ComplexMatrix& TmpM, int MatDim, int i);


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



  int LandauLevel = 0;
  AbstractFunctionBasis* Basis;
  if (CylinderFlag == false)
   {
     if (LandauLevel == 0)
       Basis = new ParticleOnSphereFunctionBasis(NbrFluxQuanta);
     else
       Basis = new ParticleOnSphereGenericLLFunctionBasis(NbrFluxQuanta - (2 * LandauLevel), LandauLevel);
   }
  else
   {
       Basis = new ParticleOnCylinderFunctionBasis(NbrFluxQuanta, LandauLevel, AspectRatio);
   }

  int LaughlinIndex = Manager.GetInteger("laughlin-index");
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



  int MatDim = SparseBMatrices[0].GetNbrRow();
  int LambdaMax = Manager.GetInteger("p-truncation");
 
  int Nmax =  LambdaMax + (LaughlinIndex - 1)/2;
  int Nmin =  -LambdaMax - (LaughlinIndex - 1)/2;

  cout << "B matrix size = " << SparseBMatrices[0].GetNbrRow() << "x" << SparseBMatrices[0].GetNbrColumn() << endl;

  cout<<"  ********************* B0 ********************* " <<endl;
  for (int i = 0; i < MatDim; i++)
    for (int j = 0; j < MatDim; j++)
      {
         Complex Tmp;
         BMatrices[0].GetMatrixElement(i,j,Tmp);
         if (Norm(Tmp) != 0.0)
           cout<<"i= "<<i<<" j= "<<j<<" "<<Tmp<<endl;
      }


  cout<<"  ********************* B1 ********************* " <<endl;
  for (int i = 0; i < MatDim; i++)
    for (int j = 0; j < MatDim; j++)
      {
         Complex Tmp;
         BMatrices[1].GetMatrixElement(i,j,Tmp);
         if (Norm(Tmp) != 0.0)
           cout<<"i= "<<i<<" j= "<<j<<" "<<Tmp<<endl;
      }




  ComplexMatrix FinalB0, FinalB1, FirstB0, FirstB1;

  FinalB0 = ComplexMatrix(MatDim, MatDim, true);
  FinalB1 = ComplexMatrix(MatDim, MatDim, true);
  FirstB0 = ComplexMatrix(MatDim, MatDim, true);
  FirstB1 = ComplexMatrix(MatDim, MatDim, true);

  for(int i = 0; i < MatDim; i++)
   {
    Complex Tmp;
    BMatrices[0].GetMatrixElement(i, (Nmax-Nmin)/2, Tmp);
    FinalB0.SetMatrixElement(i, (Nmax-Nmin)/2, Tmp);
    BMatrices[1].GetMatrixElement(i, (Nmax-Nmin)/2, Tmp);
    FinalB1.SetMatrixElement(i, (Nmax-Nmin)/2, Tmp);
    BMatrices[0].GetMatrixElement((Nmax-Nmin)/2, i, Tmp);
    FirstB0.SetMatrixElement((Nmax-Nmin)/2, i, Tmp);
    BMatrices[1].GetMatrixElement((Nmax-Nmin)/2, i, Tmp);
    FirstB1.SetMatrixElement((Nmax-Nmin)/2, i, Tmp);
   }

  //SparseComplexMatrix SparseFinalB0 = FinalB0;
  //SparseComplexMatrix SparseFinalB1 = FinalB1;
  //SparseComplexMatrix SparseFirstB0 = FirstB0;
  //SparseComplexMatrix SparseFirstB1 = FirstB1;


//  cout<<"  ********************* SparseFirstB0 ********************* " <<endl;
//  cout << FirstB0;

//  cout<<"  ********************* SparseFirstB1 ********************* " <<endl;
//  cout << FirstB1;

//  cout<<"  ********************* SparseFinalB0 ********************* " <<endl;
//  cout << FinalB0;

//  cout<<"  ********************* SparseFinalB1 ********************* " <<endl;
//  cout << FinalB1;

  int EntCut = 5;
  double CutOff = 1e-12;

  ComplexMatrix TmpM(MatDim, MatDim, true);

  ComplexMatrix TmpM0 = FirstB0.GetAdjoint();
  TmpM0.Multiply(FirstB0);
  TmpM = TmpM + TmpM0;
  ComplexMatrix TmpM1 = FirstB1.GetAdjoint();
  TmpM1.Multiply(FirstB1);
  TmpM = TmpM + TmpM1;

  for (int i = 1; i < EntCut; i++)
    {
      TmpM0 = BMatrices[0].GetAdjoint();
      TmpM0.Multiply(TmpM);
      TmpM0.Multiply(BMatrices[0]);

      TmpM1 = BMatrices[1].GetAdjoint();
      TmpM1.Multiply(TmpM);
      TmpM1.Multiply(BMatrices[1]);

      TmpM = (TmpM0 + TmpM1);
    }

  cout<<"-------------------- TmpM ------------------ "<<endl;
  for (int i = 0; i < MatDim; i++)
    for (int j = 0; j < MatDim; j++)
      {
         Complex Tmp;
         TmpM.GetMatrixElement(i,j,Tmp);
         if (Norm(Tmp) != 0.0)
           cout<<"i= "<<i<<" j= "<<j<<" "<<Tmp<<endl;
      }


  ComplexMatrix rhoM(MatDim, MatDim, true);
  ComplexMatrix rhoM0 = FinalB0;
  rhoM0.Multiply(FinalB0.GetAdjoint());
  rhoM = rhoM + rhoM0;
  ComplexMatrix rhoM1 = FinalB1;
  rhoM1.Multiply(FinalB1.GetAdjoint());
  rhoM = rhoM + rhoM1;

  for (int i = (NbrFluxQuanta - 1); i >= EntCut; i--)
    {
      rhoM0.Copy(BMatrices[0]);
      rhoM0.Multiply(rhoM);
      rhoM0.Multiply(BMatrices[0].GetAdjoint());

      rhoM1.Copy(BMatrices[1]);
      rhoM1.Multiply(rhoM);
      rhoM1.Multiply(BMatrices[1].GetAdjoint());

      rhoM = (rhoM0 + rhoM1);
    }

  cout<<"-------------------- rhoM ------------------ "<<endl;
 
  for (int i = 0; i < MatDim; i++)
    for (int j = 0; j < MatDim; j++)
      {
         Complex Tmp;
         rhoM.GetMatrixElement(i,j,Tmp);
         if (Norm(Tmp) != 0.0)
           cout<<"i= "<<i<<" j= "<<j<<" "<<Tmp<<endl;
      }


  delete[] BMatrices;


  HermitianMatrix HRep = TmpM;

  RealDiagonalMatrix TmpDiag (MatDim);
  ComplexMatrix Q(MatDim, MatDim);
  HRep.LapackDiagonalize(TmpDiag, Q);


//  cout<<"-------------------- Q ------------------ "<<endl;
//  for (int i = 0; i < 3; i++)
//      {
//         cout<<"i= "<<i<<" Q= "<<Q[i]<<endl;
//      }


  cout<<"TmpM "<<endl;
  int NonZeroEig = 0;
  for (int j = 0; j < MatDim; ++j)
   {
     cout<<TmpDiag[j]<<" Eigenvector: ";
     if (fabs(TmpDiag[j])>CutOff)
      {  
       NonZeroEig++;
       Q[j] /= sqrt(fabs(TmpDiag[j]));
       for (int i = 0; i < MatDim; i++)
         if (Norm(Q[j][i]) != 0) cout<<"i= "<<i<<" "<<Q[j][i]<<"; ";
      }  
     cout<<endl;
   }
  cout<<endl;
  cout<<"Nbr of nonzero= "<<NonZeroEig<<endl;

  cout<<"Test overlap"<<endl;

  for (int i = MatDim - 1; i >= (MatDim - NonZeroEig); --i)
    {
       Complex ov(0.0,0.0);
       for (int k=0;k<MatDim;k++)
        for (int l=0;l<MatDim;l++)
          {
            Complex Tmp;
            TmpM.GetMatrixElement(k,l,Tmp);
            ov+= Q[i][l] * Conj(Q[i][k]) * Conj(Tmp);
          }
       cout<<i<<" "<<ov<<" ; ";  
    }

  ComplexMatrix NewrhoM(NonZeroEig, NonZeroEig, true);

  int offset = MatDim - NonZeroEig;
  for (int i = MatDim - 1; i >= (MatDim - NonZeroEig); --i)
    for (int j = MatDim - 1; j >= (MatDim - NonZeroEig); --j)
      {
        Complex TmpEl(0.0, 0.0);
        for (int ii = 0; ii < MatDim; ++ii)
          for (int jj = 0; jj < MatDim; ++jj)
            {
               Complex Tmp;
               rhoM.GetMatrixElement(ii, jj, Tmp);
               TmpEl += Conj(Tmp) * OverlapEig(Q[j], TmpM, MatDim, jj) * Conj(OverlapEig(Q[i], TmpM, MatDim, ii));
            }
        NewrhoM.SetMatrixElement(i - offset, j - offset, TmpEl); 
      }

  cout<<"&&&&&&&&&&&&&&& NewrhoM &&&&&&&&&&&&&&&&&"<<endl;
  cout<<NewrhoM<<endl;

  HermitianMatrix HRepRho = NewrhoM;

  RealDiagonalMatrix TmpDiagRho (NonZeroEig);
  ComplexMatrix QRho(NonZeroEig, NonZeroEig);
  HRepRho.LapackDiagonalize(TmpDiagRho, QRho);


  TmpDiagRho /= TmpDiagRho.Tr();

  cout<<"Entanglement spectrum: "<<endl;
  for (int j = 0; j < NonZeroEig; ++j)
   {
     cout<<TmpDiagRho[j]<<" ";
   }
  cout<<endl;

/*
  int TmpIndex = Manager.GetInteger("p-truncation") + ((LaughlinIndex - 1) / 2);
  TmpIndex = TmpIndex *  SparseBMatrices[0].GetNbrRow() + TmpIndex;  

  FermionOnSphereMPSWrapper* SpaceWrapper = 0;
  if (CylinderFlag == false)
    {
      SpaceWrapper = new FermionOnSphereMPSWrapper  (NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState, TmpIndex, TmpIndex, SparseBMatrices);
    }
  else
    {
      SpaceWrapper = new FermionOnCylinderMPSWrapper  (NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState, TmpIndex, TmpIndex, SparseBMatrices, Manager.GetInteger("memory") << 20);
    }
  RealVector DummyState (1);
  DummyState[0] = 1.0;


  Complex TmpValue;
  RealVector Value(2, true);
  Complex* PrecalculatedValues = new Complex [NbrFluxQuanta + 1];	  
  if (DensityFlag == false)
    {
      for (int i = 0; i <= NbrFluxQuanta; ++i)
	{
	  Basis->GetFunctionValue(Value, TmpValue, NbrFluxQuanta);
	  ParticleOnSphereDensityDensityOperator Operator (SpaceWrapper, i, NbrFluxQuanta, i, NbrFluxQuanta);
	  PrecalculatedValues[i] = Operator.MatrixElement(DummyState, DummyState);// * TmpValue * Conj(TmpValue);
	}
    }
  else
    {
      cout<<"density precalculate ";
      Complex CheckSum (0.0,0.0);
      for (int i = 0; i <= NbrFluxQuanta; ++i)
	{
	  ParticleOnSphereDensityOperator Operator (SpaceWrapper, i);
	  PrecalculatedValues[i] = Operator.MatrixElement(DummyState, DummyState);
          CheckSum += PrecalculatedValues[i];
          cout<<i<<" "<<PrecalculatedValues[i]<<endl;
	}
      cout<<"done. Checksum="<<CheckSum<<endl;
    }
  
  ofstream File;
  File.precision(14);

  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = new char [512];
      if (DensityFlag == true)      
	{
	  sprintf(TmpFileName, "fermions_laughlin%ld_plevel_%ld_n_%d_2s_%d_lz_%d.0.rho", Manager.GetInteger("laughlin-index"),
		  Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz);
	}
      else
 	{
	  sprintf(TmpFileName, "fermions_laughlin%ld_plevel_%ld_n_%d_2s_%d_lz_%d.0.rhorho", Manager.GetInteger("laughlin-index"),
		  Manager.GetInteger("p-truncation"), NbrParticles, NbrFluxQuanta, TotalLz);
	}
     File.open(TmpFileName, ios::binary | ios::out);     
   }
  if (DensityFlag == true)      
    File << "# density  coefficients  "  << endl;
  else
    File << "# pair correlation coefficients " << endl;
  File << "#" << endl << "# (l+S)    n_l" << endl;
  if (CoefficientOnlyFlag == false)
    {
      for (int i = 0; i <= NbrFluxQuanta; ++i)
	File << "# " << i << " " << PrecalculatedValues[i].Re<< " " << PrecalculatedValues[i].Im << endl;
    }
  else
    {
      for (int i = 0; i <= NbrFluxQuanta; ++i)
	File << i << " " << PrecalculatedValues[i].Re<< " " << PrecalculatedValues[i].Im << endl;
    }


  if (CylinderFlag == false)
  {
    Complex Sum (0.0, 0.0);
    Complex Sum2 (0.0, 0.0);
    double X = 0.0;
    double XInc = M_PI / ((double) NbrPoints);
    if (CoefficientOnlyFlag == false)
      {
        double Factor1 = (16.0 * M_PI * M_PI) / ((double) (NbrParticles * NbrParticles));
        if (DensityFlag == true)
	  Factor1 = 1.0;//4.0 * M_PI;
        double Factor2;
        if (Manager.GetBoolean("radians") == true)
	  Factor2 = 1.0;
        else
	  Factor2 = sqrt (0.5 * NbrFluxQuanta);
        for (int x = 0; x < NbrPoints; ++x)
	  {
	    Value[0] = X;
	    Sum = 0.0;
	    for (int i = 0; i <= NbrFluxQuanta; ++i)
	      {
	        Basis->GetFunctionValue(Value, TmpValue, i);
	        Sum += PrecalculatedValues[i] * (Conj(TmpValue) * TmpValue);
	      }
	    if (ChordFlag == false)
	      File << (X * Factor2) << " " << (Norm(Sum)  * Factor1) << endl;
	    else
	      File << (2.0 * Factor2 * sin (X * 0.5)) << " " << Norm(Sum)  * Factor1 << endl;
	    X += XInc;
	  }
       }
     } 
 else //cylinder
  {
    double H = sqrt(2.0 * M_PI * (NbrFluxQuanta + 1.0))/sqrt(AspectRatio);
    cout<<"Cylinder H= "<<H<<endl;
    double X = -0.5 * H;
    double XInc = H / ((double) NbrPoints);

    if (CoefficientOnlyFlag == false)
      {
        for (int x = 0; x < NbrPoints; ++x)
	  {
	    Complex Sum (0.0, 0.0);
	    for (int i = 0; i <= NbrFluxQuanta; ++i)
	      {
	        Complex TmpValue = ((ParticleOnCylinderFunctionBasis*)Basis)->GetFunctionValue(X, 0.0, (double)i-0.5*NbrFluxQuanta);
	        Sum += PrecalculatedValues[i] * (Conj(TmpValue) * TmpValue);
	      }
            File << X << " " << Norm(Sum) << endl;
	    X += XInc;
	  }
       }

   }


   File.close();
 
  delete[] PrecalculatedValues;
*/

//   cout << "correlation = ";
//   for (int m1 = 0; m1 < NbrFluxQuanta; ++m1)
//     for (int m2 = m1 + 1; m2 <= NbrFluxQuanta; ++m2)
//       for (int n1 = 0; n1 < NbrFluxQuanta; ++n1)
// 	{
// 	  int n2 = m1 + m2 - n1;
// 	  if ((n2 > n1) && (n2 <= NbrFluxQuanta))
// 	    {
// 	      cout << m1 << "," << m2 << ";" << n1 <<  "," << n2 << " = ";   
// 	      if (Space != 0)
// 		{
// 		  ParticleOnSphereDensityDensityOperator Operator (Space, m1, m2, n1, n2);
// 		  Complex TmpDensityDensity = Operator.MatrixElement(State, State);
// 		  cout << TmpDensityDensity.Re << " ";
// 		  ParticleOnSphereDensityDensityOperator Operator2 (&SpaceWrapper, m1, m2, n1, n2);
// 		  Complex TmpDensityDensity2 = Operator2.MatrixElement(DummyState, DummyState);
// 		  cout << TmpDensityDensity2.Re << " ";
// 		  if (fabs(TmpDensityDensity.Re - TmpDensityDensity2.Re) > 1e-10)
// 		    cout << " error";
// 		  cout << endl;
// 		}
// 	      else
// 		{
// 		  ParticleOnSphereDensityDensityOperator Operator2 (&SpaceWrapper, m1, m2, n1, n2);
// 		  Complex TmpDensityDensity = Operator2.MatrixElement(DummyState, DummyState);
// 		  cout << TmpDensityDensity.Re << " ";
// 		  cout << endl;
// 		}
// 	    }
// 	}

//   for (int i = 0; i <= NbrFluxQuanta; ++i)
//     {
//       cout<< "n(" << i << ") = ";
//       if (Space != 0)
// 	{
// 	  ParticleOnSphereDensityOperator Operator (Space, i);
// 	  Complex TmpDensity = Operator.MatrixElement(State, State);
// 	  cout << TmpDensity.Re << " ";
// 	}
//       ParticleOnSphereDensityOperator Operator2 (&SpaceWrapper, i);
//       Complex TmpDensity2 = Operator2.MatrixElement(DummyState, DummyState);
//       cout << TmpDensity2.Re << " ";

//       cout << endl;
//     }

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
	      bMatrices[0].SetMatrixElement(StartingIndexPerPLevel[i] + ((k * NbrNValue) + j - 1), StartingIndexPerPLevel[i] + ((k * NbrNValue) + j), Tmp);
	    }
	}
    }

  bMatrices[1] = ComplexMatrix(MatrixSize, MatrixSize, true);
  unsigned long* Partition1 = new unsigned long [pLevel + 2];
  unsigned long* Partition2 = new unsigned long [pLevel + 2];
  FactorialCoefficient Coef;

  for (int i = 0; i <= pLevel; ++i)
    {
      BosonOnDiskShort* TmpSpace1 = u1BosonBasis[i];
      //int MaxN1 = (2 * i) + laughlinIndex;
      for (int j = 0; j <= pLevel; ++j)
	{
	  BosonOnDiskShort* TmpSpace2 = u1BosonBasis[j];
	  //int MaxN2 = (2 * 2) + laughlinIndex;
	  int N1 = (2 * (j - i) + laughlinIndex - 1 + NValueShift) / 2;
	  int N2 = (2 * (j - i) - laughlinIndex + 1 + NValueShift) / 2;
	  for (int k1 = 0; k1 < TmpSpace1->GetHilbertSpaceDimension(); ++k1)
	    {
	      TmpSpace1->GetOccupationNumber(k1, Partition1);
	      for (int k2 = 0; k2 < TmpSpace2->GetHilbertSpaceDimension(); ++k2)
		{
		  TmpSpace2->GetOccupationNumber(k2, Partition2);
		  Complex Tmp = CreateLaughlinAMatrixElement(laughlinIndex, Partition1, Partition2, i, j, - (N1 + N2 - NValueShift) / 2, Coef);
                  if (cylinderFlag)
                    Tmp *= exp(-kappa*kappa*(0.5 * i + 0.5 * j + pow(N1 - NValueShift/2,2.0)/(4.0 * laughlinIndex) + pow(N2 - NValueShift/2,2.0)/(4.0 * laughlinIndex)));
		  bMatrices[1].SetMatrixElement(StartingIndexPerPLevel[i] + ((k1 * NbrNValue) + N1), StartingIndexPerPLevel[j] + ((k2 * NbrNValue) + N2), Tmp);
//		  cout << i << " " << j << " | " << k1 << " " << k2 << " | " << N1 << " " << N2 << " " << (-(N1 + N2 - NValueShift) / 2) << " : " << Tmp << endl;
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

Complex OverlapEig(ComplexVector& EigV, ComplexMatrix& TmpM, int MatDim, int i)
{
  Complex Tmp;
  Complex ov(0.0, 0.0);
  for(int ii = 0; ii < MatDim; ii++)
   {
      TmpM.GetMatrixElement(ii, i, Tmp); 
      ov += EigV[ii] * Tmp;
   } 
 return ov;
}
