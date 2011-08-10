#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneLargeBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "GeneralTools/ConfigurationParser.h"
#include "MathTools/ClebschGordanCoefficients.h"

#include "Operator/ParticleOnSphereDensityDensityOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereFermionsCorrelation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('e', "eigenstate", "name of the file containing the eigenstate");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "landau-level", "index of the Landau level (0 being the LLL)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
  (*SystemGroup) += new BooleanOption  ('\n', "large-basis", "use large Hilbert space support (i.e. handle non-squeezed Hilbert space larger than 2^31 without hard-drive storage)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "file-size", "maximum file size (in MBytes) when using huge mode", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new SingleIntegerOption  ('n', "nbr-points", "number of point to evaluate", 1000);
  (*SystemGroup) += new BooleanOption  ('r', "radians", "set units to radians instead of magnetic lengths", false);
  (*SystemGroup) += new BooleanOption  ('c', "chord", "use chord distance instead of distance on the sphere", false);
  (*SystemGroup) += new BooleanOption  ('\n', "density", "plot density insted of density-density correlation", false);
  (*SystemGroup) += new BooleanOption  ('\n', "structure-factor", "evaluate LLL structure factor instead of (density-)density (use with radians option)", false);
  (*SystemGroup) += new BooleanOption  ('\n', "coefficients-only", "only compute the one or two body coefficients that are requested to evaluate the density-density correlation", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "huge-memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "large-memory", "maximum memory (in kBytes) that can allocated for precalculations when using huge mode", 1);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with rhorho extension");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsCorrelation -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = Manager.GetInteger("total-lz");
  int LandauLevel = Manager.GetInteger("landau-level");
  int NbrPoints = Manager.GetInteger("nbr-points");
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  bool DensityFlag = Manager.GetBoolean("density");
  bool StructureFactorFlag = Manager.GetBoolean("structure-factor");
  bool ChordFlag = Manager.GetBoolean("chord");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  bool SymmetrizedBasis = Manager.GetBoolean("symmetrized-basis");
  bool CoefficientOnlyFlag = Manager.GetBoolean("coefficients-only");
  bool Statistics = true;
  if (Manager.GetString("eigenstate") == 0)
    {
      cout << "FQHESphereFermionsCorrelation requires a state" << endl;
      return -1;
    }

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("eigenstate"),
						  NbrParticles, LzMax, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("eigenstate") << endl;
      return -1;
    }

  if (IsFile(Manager.GetString("eigenstate")) == false)
    {
      cout << "can't find vector file " << Manager.GetString("eigenstate") << endl;
      return -1;      
    }

  ParticleOnSphere* Space = 0;
  if (HaldaneBasisFlag == false)
    {
#ifdef __64_BITS__
      if (LzMax <= 62)
#else
	if (LzMax <= 30)
#endif
	  if ((SymmetrizedBasis == false) || (TotalLz != 0))
	    {
	      Space = new FermionOnSphere(NbrParticles, TotalLz, LzMax, MemorySpace);
	    }
	  else
	    {
	      if (Manager.GetString("load-hilbert") != 0)
		Space = new FermionOnSphereSymmetricBasis(Manager.GetString("load-hilbert"), MemorySpace);
	      else
		Space = new FermionOnSphereSymmetricBasis(NbrParticles, LzMax, MemorySpace);
	      if (Manager.GetString("save-hilbert") != 0)
		{
		  ((FermionOnSphereSymmetricBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		  return 0;
		}
	    }
	else
#ifdef __128_BIT_LONGLONG__
	  if (LzMax <= 126)
#else
	    if (LzMax <= 62)
#endif
	      {
		if ((SymmetrizedBasis == false) || (TotalLz != 0))
		  Space = new FermionOnSphereLong(NbrParticles, TotalLz, LzMax, MemorySpace);
		else
		  {
		    if (Manager.GetString("load-hilbert") != 0)
		      Space = new FermionOnSphereSymmetricBasisLong(Manager.GetString("load-hilbert"), MemorySpace);
		    else
		      Space = new FermionOnSphereSymmetricBasisLong(NbrParticles, LzMax, MemorySpace);
		    if (Manager.GetString("save-hilbert") != 0)
		      {
			((FermionOnSphereSymmetricBasisLong*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			return 0;
		      }
		  }
	      }
	    else
	      Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, LzMax, MemorySpace);
    }
  else
    {
      int* ReferenceState = 0;
      if (Manager.GetString("reference-file") == 0)
	{
	  return -1;
	}
      else
	if (FQHEGetRootPartition(Manager.GetString("reference-file"), NbrParticles,LzMax, ReferenceState) == false)
	  return -1;
      if (SymmetrizedBasis == false)
	{
	  if (Manager.GetBoolean("large-basis") == true)
	    {
	      if (Manager.GetString("load-hilbert") != 0)
		Space = new FermionOnSphereHaldaneLargeBasis(Manager.GetString("load-hilbert"), Manager.GetInteger("large-memory") << 10);
	      else
		{
		  Space = new FermionOnSphereHaldaneLargeBasis(NbrParticles, TotalLz, LzMax, ReferenceState, Manager.GetInteger("large-memory") << 10);	  
		  if (Manager.GetString("save-hilbert") != 0)
		    {
		      ((FermionOnSphereHaldaneLargeBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		      return 0;
		    }
		}
	    }
	  else
	    {
	      if (Manager.GetBoolean("huge-basis") == true)
		{
		  if (Manager.GetString("save-hilbert") != 0)
		    {
		      Space = new FermionOnSphereHaldaneHugeBasis (NbrParticles, TotalLz, LzMax, Manager.GetInteger("file-size"), ReferenceState, ((unsigned long) Manager.GetInteger("huge-memory")) << 20, false);
		      ((FermionOnSphereHaldaneHugeBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		      return 0;
		    }
		  if (Manager.GetString("load-hilbert") == 0)
		    {
		      cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
		      return -1;
		    }
		  Space = new FermionOnSphereHaldaneHugeBasis (Manager.GetString("load-hilbert"), Manager.GetInteger("huge-memory"));
		}
	      else
		{
#ifdef __64_BITS__
		  if (LzMax <= 62)
#else
		    if (LzMax <= 30)
#endif
		      {
			if (Manager.GetString("load-hilbert") != 0)
			  Space = new FermionOnSphereHaldaneBasis(Manager.GetString("load-hilbert"), MemorySpace);
			else
			  Space = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
			if (Manager.GetString("save-hilbert") != 0)
			  {
			    ((FermionOnSphereHaldaneBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			    return 0;
			  }
		      }
		    else
#ifdef __128_BIT_LONGLONG__
		      if (LzMax <= 126)
#else
			if (LzMax <= 62)
#endif
			  {
			    if (Manager.GetString("load-hilbert") != 0)
			      Space = new FermionOnSphereHaldaneBasisLong(Manager.GetString("load-hilbert"), MemorySpace);
			    else
			      Space = new FermionOnSphereHaldaneBasisLong(NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
			    if (Manager.GetString("save-hilbert") != 0)
			      {
				((FermionOnSphereHaldaneBasisLong*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
				return 0;
			      }
			  }	       
		}
	    }
	}
      else
	{
#ifdef __64_BITS__
	  if (LzMax <= 62)
#else
	    if (LzMax <= 30)
#endif
	      {
		if (Manager.GetString("load-hilbert") != 0)
		  Space = new FermionOnSphereHaldaneSymmetricBasis(Manager.GetString("load-hilbert"), MemorySpace);
		else
		  Space = new FermionOnSphereHaldaneSymmetricBasis(NbrParticles, LzMax, ReferenceState, MemorySpace);
		if (Manager.GetString("save-hilbert") != 0)
		  {
		    ((FermionOnSphereHaldaneSymmetricBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			 return 0;
		  }
	      }
	    else
#ifdef __128_BIT_LONGLONG__
	      if (LzMax <= 126)
#else
		if (LzMax <= 62)
#endif
		  {
		    if (Manager.GetString("load-hilbert") != 0)
		      Space = new FermionOnSphereHaldaneSymmetricBasisLong(Manager.GetString("load-hilbert"), MemorySpace);
		    else
		      Space = new FermionOnSphereHaldaneSymmetricBasisLong(NbrParticles, LzMax, ReferenceState, MemorySpace);
		    if (Manager.GetString("save-hilbert") != 0)
		      {
			((FermionOnSphereHaldaneSymmetricBasisLong*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			return 0;
		      }
		  }
	}
    }

  cout << Space->GetHilbertSpaceDimension() << endl;

  if (StructureFactorFlag == true)
   {
     cout<<"Structure factor is evaluated for LLL; " << endl;
     cout<<"use normalization convention from He, Simon & Halperin." << endl;
     ofstream File;
     File.precision(14);
     if (Manager.GetString("output-file") != 0)
       File.open(Manager.GetString("output-file"), ios::binary | ios::out);
     else
      {
        cout << "Enter output file! " << endl;
        exit(1);
      }
  
     RealVector State;
     if (State.ReadVectorTest(Manager.GetString("eigenstate")) == true)
      {
        if (State.ReadVector (Manager.GetString("eigenstate")) == false)
	  {
	    cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
	    return -1;      
	  }

        Complex* DensityMatEl;
        DensityMatEl = new Complex[LzMax + 1];

        for (int i = 0; i <= LzMax; ++i)
         {
           ParticleOnSphereDensityOperator Operator (Space, i);
           DensityMatEl[i] = Operator.MatrixElement(State, State);
         }

       double S = 0.5 * (double)LzMax;
       ClebschGordanCoefficients CoeffLLS(LzMax, LzMax);  
       for (int L = 0; L <= LzMax; ++L)
        {
          double Factor = pow(-1.0, 3.0 * LzMax + (L << 1)) * pow(LzMax + 1.0, 2.0) /(4.0 * M_PI * (2.0 * L + 1.0));
          Factor *=  pow(CoeffLLS.GetCoefficient(-LzMax, LzMax, L << 1), 2.0);

          Complex SumIJ(0.0,0.0);
          for (int i = 0; i <= LzMax; ++i)
           {
            double FactorI = Factor * CoeffLLS.GetCoefficient(((i << 1) - LzMax), -((i << 1) - LzMax), L << 1);
            for (int j = 0; j <= LzMax; ++j)
	     {        
               double FactorJ = FactorI * CoeffLLS.GetCoefficient(((j << 1) - LzMax), -((j << 1) - LzMax), L << 1) * pow(-1.0, -(i - S)-(j-S));
               ParticleOnSphereDensityDensityOperator Operator (Space, i, j, i, j);
               SumIJ -=  FactorJ * Operator.MatrixElement(State, State);
               if (i == j)
                 SumIJ += FactorJ * DensityMatEl[i]; 
                
	     }
           }     
         cout << L <<" "<< (4.0 * M_PI/(double)NbrParticles) * SumIJ.Re << " " << (4.0 * M_PI/(double)NbrParticles) * SumIJ.Im <<endl;
         File << L <<" "<< (4.0 * M_PI/(double)NbrParticles) * SumIJ.Re << " " << (4.0 * M_PI/(double)NbrParticles) * SumIJ.Im <<endl;
        }
       delete[] DensityMatEl;
       return 0;
     } 
   }


  AbstractFunctionBasis* Basis;
  if (LandauLevel == 0)
    Basis = new ParticleOnSphereFunctionBasis(LzMax);
  else
    Basis = new ParticleOnSphereGenericLLFunctionBasis(LzMax - (2 * LandauLevel), LandauLevel);

  Complex Sum (0.0, 0.0);
  Complex Sum2 (0.0, 0.0);
  Complex TmpValue;
  RealVector Value(2, true);
  double X = 0.0;
  double XInc = M_PI / ((double) NbrPoints);

  Complex* PrecalculatedValues = new Complex [LzMax + 1];
  RealVector State;
  if (State.ReadVectorTest(Manager.GetString("eigenstate")) == true)
    {
      if (State.ReadVector (Manager.GetString("eigenstate")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
	  return -1;      
	}
      if (DensityFlag == false)
	for (int i = 0; i <= LzMax; ++i)
	  {
	    Basis->GetFunctionValue(Value, TmpValue, LzMax);
	    ParticleOnSphereDensityDensityOperator Operator (Space, i, LzMax, i, LzMax);
	    PrecalculatedValues[i] = Operator.MatrixElement(State, State) * TmpValue * Conj(TmpValue);
	  }
      else
	for (int i = 0; i <= LzMax; ++i)
	  {
	    ParticleOnSphereDensityOperator Operator (Space, i);
	    PrecalculatedValues[i] = Operator.MatrixElement(State, State);
	  }
    }
  else
    {
      ComplexVector ComplexState;
      if (ComplexState.ReadVector (Manager.GetString("eigenstate")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
	  return -1;      
	}
      if (DensityFlag == false)
	for (int i = 0; i <= LzMax; ++i)
	  {
	    Basis->GetFunctionValue(Value, TmpValue, LzMax);
	    ParticleOnSphereDensityDensityOperator Operator (Space, i, LzMax, i, LzMax);
	    PrecalculatedValues[i] = Operator.MatrixElement(ComplexState, ComplexState) * TmpValue * Conj(TmpValue);
	  }
      else
	for (int i = 0; i <= LzMax; ++i)
	  {
	    ParticleOnSphereDensityOperator Operator (Space, i);
	    PrecalculatedValues[i] = Operator.MatrixElement(ComplexState, ComplexState);
	  }
    }

  ofstream File;
  File.precision(14);
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = 0;
      if (DensityFlag == false)
	{
	  if (Manager.GetBoolean("coefficients-only"))
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rhorho-c");
	  else
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rhorho");
	}
      else
	{
	    TmpFileName = ReplaceExtensionToFileName(Manager.GetString("eigenstate"), "vec", "rho");
	}
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << Manager.GetString("eigenstate") << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  if (DensityFlag == true)      
    File << "# density correlation coefficients for " << Manager.GetString("eigenstate") << endl;
  else
    File << "# density-density correlation coefficients for " << Manager.GetString("eigenstate") << endl;
  File << "#" << endl << "# (l+S)    n_l" << endl;
  if (CoefficientOnlyFlag == false)
    {
      for (int i = 0; i <= LzMax; ++i)
	File << "# " << i << " " << PrecalculatedValues[i].Re<< endl;
    }
  else
    {
      for (int i = 0; i <= LzMax; ++i)
	File << i << " " << PrecalculatedValues[i].Re<< endl;
    }
  if (CoefficientOnlyFlag == false)
    {
      double Factor1 = (16.0 * M_PI * M_PI) / ((double) (NbrParticles * NbrParticles));
      if (DensityFlag == true)
	Factor1 = 1.0;//4.0 * M_PI;
      double Factor2;
      if (Manager.GetBoolean("radians") == true)
	Factor2 = 1.0;
      else
	Factor2 = sqrt (0.5 * LzMax);
      for (int x = 0; x < NbrPoints; ++x)
	{
	  Value[0] = X;
	  Sum = 0.0;
	  for (int i = 0; i <= LzMax; ++i)
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
  File.close();
 
  delete[] PrecalculatedValues;

  return 0;
}


