#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"

#include "Operator/ParticleOnTorusKxOperator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/IntegerAlgebraTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Tools/FQHEFiles/FQHEOnTorusFileTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>
#include <cstring> 

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHETorusComputeKx" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption('i', "input-state", "name of the file containing the state whose Kx momentum has to be computed");
  (*SystemGroup) += new SingleStringOption('\n', "degenerated-states", "name of the file containing a list of states (override input-state)");
  (*SystemGroup) += new SingleDoubleOption   ('r', "ratio", "ratio between lengths along the x and y directions", 1);
  (*SystemGroup) += new BooleanOption ('\n',  "compute-eigenstate", "compute the eigenstates of th Kx operator in the given basis");
  (*SystemGroup) += new SingleStringOption ('\n',  "interaction-name", "name that should be inserted in the output file names", "dummy");
  (*SystemGroup) += new BooleanOption ('\n',  "no-convertion", "do not convert the final vectors to the (Kx,Ky) n-body basis");
  (*SystemGroup) += new BooleanOption ('\n',  "export-transformation", "export the transformation matrix in a ascii file (one per momentum sector)");
  (*SystemGroup) += new BooleanOption ('\n',  "export-bintransformation", "export the transformation matrix in a binary file (one per momentum sector)");
  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusComputeKx -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int MaxMomentum = 0;
  int YMomentum = -1;
  bool Statistics = true;

  
  RealVector* InputStates = 0;
  int NbrInputStates = 0;
  if (Manager.GetString("degenerated-states") == 0)
    {
      if (Manager.GetString("input-state") == 0)
	{
	  cout << "error, either input-state or degenerated-states has to be provided" << endl;
	  return -1;
	}
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						      NbrParticles, MaxMomentum, YMomentum, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("ground-state") << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Ky=" << YMomentum << " ";
      NbrInputStates = 1;
      InputStates = new RealVector [NbrInputStates];
      if (InputStates[0].ReadVector(Manager.GetString("input-state")) == false)
	{
	  cout << "error while reading " << Manager.GetString("input-state") << endl;
	  return -1;
	}
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerated-states")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	} 
      NbrInputStates = DegeneratedFile.GetNbrLines();
      if (NbrInputStates < 1)
	{
	  cout << "no state found in " << Manager.GetString("degenerated-states") << endl;
	}
      InputStates = new RealVector [NbrInputStates];
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(DegeneratedFile(0, 0),
						      NbrParticles, MaxMomentum, YMomentum, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, 0) << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Ky=" << YMomentum << " ";
      if (InputStates[0].ReadVector(DegeneratedFile(0, 0)) == false)
	{
	  cout << "error while reading " << DegeneratedFile(0, 0) << endl;
	  return -1;
	}
      for (int i = 1; i < NbrInputStates; ++i)
	{
	  int TmpNbrParticles = 0;
	  int TmpMaxMomentum = 0;
	  int TmpYMomentum = -1;
	  bool TmpStatistics = true;
	  if (FQHEOnTorusFindSystemInfoFromVectorFileName(DegeneratedFile(0, i),
							  TmpNbrParticles, TmpMaxMomentum, TmpYMomentum, TmpStatistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, i) << endl;
	      return -1;
	    }
	  if ((TmpNbrParticles != NbrParticles) || (TmpMaxMomentum != MaxMomentum) || 
	      (TmpYMomentum != YMomentum) || (Statistics != TmpStatistics))
	    {
	      cout << "error, " << DegeneratedFile(0, i) << " has different system parameters than " << DegeneratedFile(0, 0) 
		   << ", N=" << TmpNbrParticles << "(" << NbrParticles << "), N_phi=" << TmpMaxMomentum << "(" << MaxMomentum 
		   << "), Ky=" << TmpYMomentum << "(" << YMomentum << ")" << endl;
	    }
	  if (InputStates[i].ReadVector(DegeneratedFile(0, i)) == false)
	    {
	      cout << "error while reading " << DegeneratedFile(0, i) << endl;
	      return -1;
	    }
	  if (InputStates[i].GetVectorDimension() != InputStates[0].GetVectorDimension())
	    {
	      cout << "error, " << DegeneratedFile(0, i) << " has different dimension than " << DegeneratedFile(0, 0) << endl;	      
	    }
	}
    }
  int MomentumModulo = FindGCD(NbrParticles, MaxMomentum);
  ParticleOnTorus* TotalSpace = 0;
  ParticleOnTorusWithMagneticTranslations** TargetSpaces = 0;
  char* OutputNamePrefix = new char [256 + strlen(Manager.GetString("interaction-name"))];

  if (Statistics == false)
    {
#ifdef  __64_BITS__
      if ((MaxMomentum + NbrParticles - 1) < 63)
#else
	if ((MaxMomentum + NbrParticles - 1) < 31)	
#endif
	  {
	    TotalSpace = new BosonOnTorusShort(NbrParticles, MaxMomentum, YMomentum);
	  }
	else
	  {
	    TotalSpace = new BosonOnTorus(NbrParticles, MaxMomentum, YMomentum);
	  }
      sprintf (OutputNamePrefix, "bosons_torus_%s_n_%d_2s_%d", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum);
      if (Manager.GetBoolean("compute-eigenstate") == true)
	{
	  TargetSpaces = new ParticleOnTorusWithMagneticTranslations*[MomentumModulo];
	  for (int i = 0; i < MomentumModulo; ++i)
	    TargetSpaces[i] = new BosonOnTorusWithMagneticTranslationsShort(NbrParticles, MaxMomentum, i, YMomentum);
	}
    }
  else
    {
      TotalSpace = new FermionOnTorus (NbrParticles, MaxMomentum, YMomentum);
      sprintf (OutputNamePrefix, "fermions_torus_%s_n_%d_2s_%d", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum);
      if (Manager.GetBoolean("compute-eigenstate") == true)
	{
	  TargetSpaces = new ParticleOnTorusWithMagneticTranslations*[MomentumModulo];
	  for (int i = 0; i < MomentumModulo; ++i)
	    TargetSpaces[i] = new FermionOnTorusWithMagneticTranslations(NbrParticles, MaxMomentum, i, YMomentum);
	}
    }
  if (InputStates[0].GetVectorDimension() != TotalSpace->GetHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions " << InputStates[0].GetVectorDimension() 
	   << " "<< TotalSpace->GetHilbertSpaceDimension() << endl;
      return -1;
    }

  Architecture.GetArchitecture()->SetDimension(TotalSpace->GetHilbertSpaceDimension());

  RealMatrix KxRep (NbrInputStates, NbrInputStates, true);
  ParticleOnTorusKxOperator KxOperator(TotalSpace);

  RealVector TmpVector(TotalSpace->GetHilbertSpaceDimension());
  for (int i = 0; i < NbrInputStates; ++i)
    {
      KxOperator.Multiply(InputStates[i], TmpVector);
      for (int j = 0; j < NbrInputStates; ++j)
	{
	  KxRep.SetMatrixElement(j, i, (InputStates[j] * TmpVector));
	}
    }

  char* OutputName = new char [256 + strlen(OutputNamePrefix)];
  sprintf (OutputName, "%s_ky_%d.dat", OutputNamePrefix, YMomentum);
  ofstream File;
  File.open(OutputName, ios::binary | ios::out);
  File.precision(14);
  File << "# eigenvalue Norm Arg (GCD(N,N_phi)*Arg/2pi) Kx round(Kx)" << endl;
  if (NbrInputStates == 1)
    {
      double TmpValue;
      KxRep.GetMatrixElement(0, 0, TmpValue);
      File << TmpValue << " " << fabs(TmpValue) << " ";
      if (TmpValue < 0.0)
	{
	  File << M_PI << " " << (MomentumModulo / 2.0) << " " << (MomentumModulo / 2.0) << endl;
	}
      else
	{
	  File << 0.0 << " " << 0.0 <<  " " << 0 << endl;
	}
      if (Manager.GetBoolean("compute-eigenstate") == true)
	{
	  char* VectorOutputName = new char [256 + strlen(OutputNamePrefix)];
	  ComplexVector TmpVector (InputStates[0], true);
	  ComplexVector TmpVector2;
	  if (TmpValue < 0.0)
	    {
	      sprintf (VectorOutputName, "%s_kx_%d_ky_%d.0.vec", OutputNamePrefix, (MomentumModulo / 2), YMomentum);
	      if (Manager.GetBoolean("no-convertion") == false)
		TmpVector2 = TargetSpaces[MomentumModulo / 2]->ConvertToKxKyBasis(TmpVector, TotalSpace);
	      else
		TmpVector2 = TmpVector;

	    }
	  else
	    {
	      sprintf (VectorOutputName, "%s_kx_%d_ky_%d.0.vec", OutputNamePrefix, 0, YMomentum);
	      if (Manager.GetBoolean("no-convertion") == false)
		TmpVector2 = TargetSpaces[0]->ConvertToKxKyBasis(TmpVector, TotalSpace);
	      else
		TmpVector2 = TmpVector;	      
	    }
	  if (TmpVector2.WriteVector(VectorOutputName) == false)
	    {
	      cout << "error, can't write vector " << VectorOutputName << endl;
	    }
	  delete[] VectorOutputName;	  
	}
    }
  else
    {
      if (Manager.GetBoolean("compute-eigenstate") == false)
	{
	  ComplexDiagonalMatrix Eigenvalues(NbrInputStates, true);
	  KxRep.LapackDiagonalize(Eigenvalues);
	  double Factor = ((double) MomentumModulo) / (2.0 * M_PI);
	  for (int i = 0; i < NbrInputStates; ++i)
	    {
	      Complex TmpValue = Eigenvalues[i];
	      File << TmpValue << " " << Norm(TmpValue) << " " << Arg(TmpValue);
	      double TmpValue2 = Arg(TmpValue);
	      if (TmpValue2 < 0.0)
		TmpValue2 += 2.0 * M_PI;
	      TmpValue2 *= Factor;
	      File << " " << TmpValue2 << " " << round(TmpValue2) << endl;
	    }      
	}
      else
	{
	  ComplexDiagonalMatrix Eigenvalues(NbrInputStates, true);
	  ComplexMatrix Eigenstates(NbrInputStates, NbrInputStates);
	  KxRep.LapackDiagonalize(Eigenvalues, Eigenstates);
	  double Factor = ((double) MomentumModulo) / (2.0 * M_PI);
	  int* NbrStatePerKxSector = new int [MomentumModulo];
	  int* KxValues = new int [NbrInputStates] ;
	  for (int i = 0; i < MomentumModulo; ++i)
	    NbrStatePerKxSector[i] = 0;
	  for (int i = 0; i < NbrInputStates; ++i)
	    {
	      Complex TmpValue = Eigenvalues[i];
	      File << TmpValue << " " << Norm(TmpValue) << " " << Arg(TmpValue);
	      double TmpValue2 = Arg(TmpValue);
	      if (TmpValue2 < 0.0)
		TmpValue2 += 2.0 * M_PI;
	      TmpValue2 *= Factor;
	      int TmpKx = round(TmpValue2);
	      File << " " << TmpValue2 << " " << TmpKx << endl;
	      KxValues[i] = TmpKx;
	      NbrStatePerKxSector[TmpKx]++;
	    }
	  ComplexVector** TmpVectors = new ComplexVector*[MomentumModulo];
	  for (int i = 0; i < MomentumModulo; ++i)
	    {
	      TmpVectors[i]  = new ComplexVector[NbrStatePerKxSector[i]];
	      NbrStatePerKxSector[i] = 0;
	    }
	  ComplexVector TmpVector (TotalSpace->GetHilbertSpaceDimension());
 	  for (int i = 0; i < NbrInputStates; ++i)
	    {
	      Complex TmpValue = Eigenvalues[i];
	      for (int j = 0; j < TotalSpace->GetHilbertSpaceDimension(); ++j)
		{
		  Complex Tmp = 0.0;
		  for (int k = 0; k < NbrInputStates; ++k)
		    {
		      Tmp += Conj(Eigenstates[i][k]) * InputStates[k][j];
		    }
		  TmpVector[j] = Tmp;		  
		}
	      if (Manager.GetBoolean("no-convertion") == false)
		TmpVectors[KxValues[i]][NbrStatePerKxSector[KxValues[i]]] = (TargetSpaces[KxValues[i]]->ConvertToKxKyBasis(TmpVector, TotalSpace));
	      else
		TmpVectors[KxValues[i]][NbrStatePerKxSector[KxValues[i]]] = ComplexVector(TmpVector, true);
	      NbrStatePerKxSector[KxValues[i]]++;	      
	    }	  
	  
	  ComplexMatrix* TransformationMatrices = new ComplexMatrix [MomentumModulo]; 
	  for (int i = 0; i < MomentumModulo; ++i)
	    {
	      ComplexMatrix TmpMatrix (TmpVectors[i], NbrStatePerKxSector[i]);
	      TmpMatrix.OrthoNormalizeColumns(TransformationMatrices[i]);
	      for (int j = 0; j < NbrStatePerKxSector[i]; ++j)
		{
		  char* VectorOutputName = new char [256 + strlen(OutputNamePrefix)];
		  sprintf (VectorOutputName, "%s_kx_%d_ky_%d.%d.vec", OutputNamePrefix, i, YMomentum, j);
		  if (TmpMatrix[j].WriteVector(VectorOutputName) == false)
		    {
		      cout << "error, can't write vector " << VectorOutputName << endl;
		    }
		  delete[] VectorOutputName;	  	      
		}
//	      delete[] TmpVectors[i];
	    }	  
	  if ((Manager.GetBoolean("export-transformation") == true) || 
	      (Manager.GetBoolean("export-bintransformation") == true))
	    {
	      for (int i = 0; i < MomentumModulo; ++i)
		{
		  if (NbrStatePerKxSector[i] > 0)
		    {
		      ComplexMatrix TmpMatrix(NbrInputStates, NbrStatePerKxSector[i]);
		      int Tmp = 0;
		      for (int j = 0; j < NbrInputStates; ++j)
			{
			  if (KxValues[j] == i)
			    {
			      TmpMatrix[Tmp] = Eigenstates[j];
			      ++Tmp;
			    }
			}
		      TmpMatrix.ComplexConjugate();
		      TmpMatrix.Multiply(TransformationMatrices[i]);
		      char* MatrixOutputName =  new char [256 + strlen(OutputNamePrefix)];
		      if (Manager.GetBoolean("export-transformation") == true)
			{
			  sprintf (MatrixOutputName, "%s_kx_%d_ky_%d.mat.txt", OutputNamePrefix, i, YMomentum);
			  TmpMatrix.WriteAsciiMatrix(MatrixOutputName);
			}
		      else
			{
			  sprintf (MatrixOutputName, "%s_kx_%d_ky_%d.mat", OutputNamePrefix, i, YMomentum);
			  TmpMatrix.WriteMatrix(MatrixOutputName);
			}
		      delete[] MatrixOutputName;
		    }
		}
	    }
	  delete[] TmpVectors;
	}
    }
  File.close();
  return 0;

}
