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


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereMPSCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  FQHEMPSMatrixManager MPSMatrixManager;

  MPSMatrixManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* OutputGroup = Manager.GetOptionGroup("output options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "file that describes the root configuration");
  (*SystemGroup) += new BooleanOption  ('\n', "use-padding", "root partitions use the extra zero padding");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 1000);
  (*SystemGroup) += new BooleanOption  ('r', "radians", "set units to radians instead of magnetic lengths", false);
  (*SystemGroup) += new BooleanOption  ('c', "chord", "use chord distance instead of distance on the sphere", false);
  (*SystemGroup) += new BooleanOption ('\n', "cylinder", "evaluate density on the cylinder");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "aspect-ratio", "aspect ratio of the cylinder", 1);
  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density insted of density-density correlation", false);
  (*SystemGroup) += new BooleanOption  ('\n', "coefficients-only", "only compute the one or two body coefficients that are requested to evaluate the density-density correlation", false);
  (*SystemGroup) += new SingleStringOption  ('\n', "state", "provide an external state for comparison purposes");
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

  int NbrPoints = Manager.GetInteger("nbr-points");
  bool DensityFlag = Manager.GetBoolean("density");
  bool ChordFlag = Manager.GetBoolean("chord");
  bool CoefficientOnlyFlag = Manager.GetBoolean("coefficients-only");

  int NbrParticles = 0; 
  int NbrFluxQuanta = 0;
  int TotalLz = 0;
  int* ReferenceState = 0;
  if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles, NbrFluxQuanta, ReferenceState) == false)
    return -1;

  bool CylinderFlag = Manager.GetBoolean("normalize-cylinder");
  double AspectRatio = Manager.GetDouble("aspect-ratio");
  double kappa = 0.0;
  if (CylinderFlag)
    {
       kappa = (2.0 * M_PI)/sqrt(2.0 * M_PI * (NbrFluxQuanta + 1) * AspectRatio);
       cout<<"Cylinder geometry, kappa= "<<kappa<<endl;
    }


  RealVector State;
  FermionOnSpherePTruncated* Space = 0;
 
  if (Manager.GetString("state") != 0)
    {
      if (State.ReadVector(Manager.GetString("state")) == false)
	{
	  cout << "can't read " << Manager.GetString("state") << endl;
	  return 0;
	}
      if (Manager.GetBoolean("boson") == true)
	{
	  cout << "bosons are not yet implemented" << endl;
	  return 0;
	}
      else
	{
#ifdef __64_BITS__
	  if (NbrFluxQuanta <= 62)
#else
	    if (NbrFluxQuanta <= 30)
#endif
	      {
		Space = new FermionOnSpherePTruncated(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
	      }
	    else
	      {
#ifdef __128_BIT_LONGLONG__
		if (NbrFluxQuanta <= 126)
#else
		  if (NbrFluxQuanta <= 62)
#endif
		    {
		      Space = 0;//new FermionOnSpherePTruncatedLong(NbrParticles, TotalLz, NbrFluxQuanta, Manager.GetInteger("p-truncation"), ReferenceState);
		    }
		  else
		    {
#ifdef __128_BIT_LONGLONG__
		      cout << "cannot generate an Hilbert space when nbr-flux > 126" << endl;
#else
		      cout << "cannot generate an Hilbert space when nbr-flux > 62" << endl;
#endif
		      return 0;
		    }
	      }
	}      
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

  FermionOnSphereMPSWrapper* SpaceWrapper = 0;
  if (CylinderFlag == false)
    {
      SpaceWrapper = new FermionOnSphereMPSWrapper  (NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState, MPSRowIndex, MPSColumnIndex, SparseBMatrices);
    }
  else
    {
      SpaceWrapper = new FermionOnCylinderMPSWrapper  (NbrParticles, TotalLz, NbrFluxQuanta, ReferenceState, MPSRowIndex, MPSColumnIndex, SparseBMatrices, Manager.GetInteger("memory") << 20);
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

