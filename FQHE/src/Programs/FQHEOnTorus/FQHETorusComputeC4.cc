#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/BosonOnTorus.h"
#include "HilbertSpace/BosonOnTorusShort.h"
#include "HilbertSpace/BosonOnTorusWithMagneticTranslationsShort.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"

#include "Operator/ParticleOnTorusC4Operator.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

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
  OptionManager Manager ("FQHETorusComputeC4" , "0.01");
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
  (*SystemGroup) += new BooleanOption ('\n',  "compute-eigenstate", "compute the eigenstates of the C4 (or C2) operator in the given basis");
  (*SystemGroup) += new SingleStringOption ('\n',  "interaction-name", "name that should be inserted in the output file names", "dummy");
  (*SystemGroup) += new BooleanOption ('\n',  "c2-only", "only check the C2 symmetry");
  (*SystemGroup) += new BooleanOption ('\n',  "export-transformation", "export the transformation matrix in a ascii file (one per momentum sector)");
  (*SystemGroup) += new BooleanOption ('\n',  "export-bintransformation", "export the transformation matrix in a binary file (one per momentum sector)");
  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETorusComputeC4 -h" << endl;
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
  bool C2OnlyFlag = Manager.GetBoolean("c2-only");
  int XMomentum = -1;
  ComplexVector* InputStates = 0;
  int NbrInputStates = 0;
  char* OutputNamePrefix = new char [256 + strlen(Manager.GetString("interaction-name"))];
  if (Manager.GetString("degenerated-states") == 0)
    {
      if (Manager.GetString("input-state") == 0)
	{
	  cout << "error, either input-state or degenerated-states has to be provided" << endl;
	  return -1;
	}
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(Manager.GetString("input-state"),
						      NbrParticles, MaxMomentum, XMomentum, YMomentum, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << Manager.GetString("ground-state") << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Kx=" << XMomentum << " " << " Ky=" << YMomentum << endl;
      if ((XMomentum != YMomentum) || (!((XMomentum == 0) || (((NbrParticles & 1) == 0) && (XMomentum == (NbrParticles / 2))))))
	{
	  cout << "C4 symmetry can only be computed in the (0,0) or (pi,pi) sectors" << endl;
	  return -1;
	}
      NbrInputStates = 1;
      InputStates = new ComplexVector [NbrInputStates];
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
      InputStates = new ComplexVector [NbrInputStates];
      if (FQHEOnTorusFindSystemInfoFromVectorFileName(DegeneratedFile(0, 0),
						      NbrParticles, MaxMomentum, XMomentum, YMomentum, Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, 0) << endl;
	  return -1;
	}
      cout << "Nbr particles=" << NbrParticles << ", Nbr flux quanta=" << MaxMomentum << " Kx=" << XMomentum << " " << " Ky=" << YMomentum << endl;
      if ((XMomentum != YMomentum) || (!((XMomentum == 0) || (((NbrParticles & 1) == 0) && (XMomentum == (NbrParticles / 2))))))
	{
	  cout << "C4 symmetry can only be computed in the (0,0) or (pi,pi) sectors" << endl;
	  return -1;
	}
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
	  int TmpXMomentum = -1;
	  bool TmpStatistics = true;
	  if (FQHEOnTorusFindSystemInfoFromVectorFileName(DegeneratedFile(0, i),
							  TmpNbrParticles, TmpMaxMomentum, TmpXMomentum, TmpYMomentum, TmpStatistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << DegeneratedFile(0, i) << endl;
	      return -1;
	    }
	  if ((TmpNbrParticles != NbrParticles) || (TmpMaxMomentum != MaxMomentum) || 
	      (TmpYMomentum != YMomentum) || (TmpXMomentum != XMomentum) || (Statistics != TmpStatistics))
	    {
	      cout << "error, " << DegeneratedFile(0, i) << " has different system parameters than " << DegeneratedFile(0, 0) 
		   << ", N=" << TmpNbrParticles << "(" << NbrParticles << "), N_phi=" << TmpMaxMomentum << "(" << MaxMomentum 
		   << "), Kx=" << TmpXMomentum << "(" << XMomentum << "), Ky=" << TmpYMomentum << "(" << YMomentum << ")" << endl;
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
  bool PiPiSectorFlag = false;
  if (XMomentum != 0)
    PiPiSectorFlag = true;
  ParticleOnTorusWithMagneticTranslations* Space = 0;
  if (Statistics == false)
    {
      Space = new BosonOnTorusWithMagneticTranslationsShort(NbrParticles, MaxMomentum, XMomentum, YMomentum);
      sprintf (OutputNamePrefix, "bosons_torus_%s_n_%d_2s_%d", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum);
    }
  else
    {
      Space = new FermionOnTorusWithMagneticTranslations(NbrParticles, MaxMomentum, XMomentum, YMomentum);
      sprintf (OutputNamePrefix, "fermions_torus_%s_n_%d_2s_%d", Manager.GetString("interaction-name"), NbrParticles, MaxMomentum);
    }
  if (InputStates[0].GetVectorDimension() != Space->GetHilbertSpaceDimension())
    {
      cout << "error: vector and Hilbert-space have unequal dimensions " << InputStates[0].GetVectorDimension() 
	   << " "<< Space->GetHilbertSpaceDimension() << endl;
      return -1;
    }

  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

  ComplexMatrix C4Rep (NbrInputStates, NbrInputStates, true);
  ParticleOnTorusC4Operator C4Operator(Space, PiPiSectorFlag, C2OnlyFlag);


  ComplexVector TmpVector(Space->GetHilbertSpaceDimension());
  for (int i = 0; i < NbrInputStates; ++i)
    {
      VectorOperatorMultiplyOperation Operation(&C4Operator, &(InputStates[i]), &TmpVector);
      Operation.ApplyOperation(Architecture.GetArchitecture());
      for (int j = 0; j < NbrInputStates; ++j)
	{
	  C4Rep.SetMatrixElement(j, i, (TmpVector * InputStates[j]));
	}
    }

  char* OutputName = new char [256 + strlen(OutputNamePrefix)];
  sprintf (OutputName, "%s_c4_kx_%d_ky_%d.dat", OutputNamePrefix, XMomentum, YMomentum);
  ofstream File;
  File.open(OutputName, ios::binary | ios::out);
  File.precision(14);
  if (C2OnlyFlag == false)
    {
      File << "# eigenvalue Norm Arg (2 * Arg/pi) round(C4)" << endl;
    }
  else
    {
      File << "# eigenvalue Norm Arg (2 * Arg/pi) round(C2)" << endl;
    }

  if (Manager.GetBoolean("compute-eigenstate") == false)
    {
      ComplexDiagonalMatrix Eigenvalues(NbrInputStates, true);
      C4Rep.LapackDiagonalize(Eigenvalues);
      double Factor = 2.0 /  M_PI;
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  Complex TmpValue = Eigenvalues[i];
	  File << TmpValue << " " << Norm(TmpValue) << " " << Arg(TmpValue);
	  double TmpValue2 = Arg(TmpValue);
	  if (TmpValue2 < 0.0)
	    TmpValue2 += 2.0 * M_PI;
	  TmpValue2 *= Factor;
	  File << " " << TmpValue2 << " " << round(TmpValue2) << endl;
	  if (NbrInputStates == 1)
	    {
	      if (C2OnlyFlag == false)
		{
		  cout << "C4 = " << round(TmpValue2) << " (" << TmpValue2 << ")" << endl;
		}
	      else
		{
		  cout << "C2 = " << round(TmpValue2) << " (" << TmpValue2 << ")" << endl;
		}
	    }
	}      
    }
  else
    {
      ComplexDiagonalMatrix Eigenvalues(NbrInputStates, true);
      ComplexMatrix Eigenstates(NbrInputStates, NbrInputStates);
      C4Rep.LapackDiagonalize(Eigenvalues, Eigenstates);
      double Factor = 2.0 /  M_PI;
      int* NbrStatePerC4Sector = new int [4];
      int* C4Values = new int [NbrInputStates] ;
      for (int i = 0; i < 4; ++i)
	NbrStatePerC4Sector[i] = 0;
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  Complex TmpValue = Eigenvalues[i];
	  File << TmpValue << " " << Norm(TmpValue) << " " << Arg(TmpValue);
	  double TmpValue2 = Arg(TmpValue);
	  if (TmpValue2 < 0.0)
	    TmpValue2 += 2.0 * M_PI;
	  TmpValue2 *= Factor;
	  int TmpC4 = round(TmpValue2);
	  File << " " << TmpValue2 << " " << TmpC4 << endl;
	  C4Values[i] = TmpC4;
	  NbrStatePerC4Sector[TmpC4]++;
	}
      ComplexVector** TmpVectors = new ComplexVector*[4];
      for (int i = 0; i < 4; ++i)
	{
	  TmpVectors[i]  = new ComplexVector[NbrStatePerC4Sector[i]];
	  NbrStatePerC4Sector[i] = 0;
	}
      ComplexVector TmpVector (Space->GetHilbertSpaceDimension());
      for (int i = 0; i < NbrInputStates; ++i)
	{
	  Complex TmpValue = Eigenvalues[i];
	  for (int j = 0; j < Space->GetHilbertSpaceDimension(); ++j)
	    {
	      Complex Tmp = 0.0;
	      for (int k = 0; k < NbrInputStates; ++k)
		{
		  Tmp += Conj(Eigenstates[i][k]) * InputStates[k][j];
		}
	      TmpVector[j] = Tmp;		  
	    }
	  NbrStatePerC4Sector[C4Values[i]]++;	      
	}	  
      
      ComplexMatrix* TransformationMatrices = new ComplexMatrix [4]; 
      for (int i = 0; i < 4; ++i)
	{
	  ComplexMatrix TmpMatrix (TmpVectors[i], NbrStatePerC4Sector[i]);
	  TmpMatrix.OrthoNormalizeColumns(TransformationMatrices[i]);
	  for (int j = 0; j < NbrStatePerC4Sector[i]; ++j)
	    {
	      char* VectorOutputName = new char [256 + strlen(OutputNamePrefix)];
	      sprintf (VectorOutputName, "%s_c4_%d_kx_%d_ky_%d.%d.vec", OutputNamePrefix, i, XMomentum, YMomentum, j);
	      if (TmpMatrix[j].WriteVector(VectorOutputName) == false)
		{
		  cout << "error, can't write vector " << VectorOutputName << endl;
		}
	      delete[] VectorOutputName;	  	      
	    }
	}	  
      if ((Manager.GetBoolean("export-transformation") == true) || 
	  (Manager.GetBoolean("export-bintransformation") == true))
	{
	  for (int i = 0; i < 4; ++i)
	    {
	      if (NbrStatePerC4Sector[i] > 0)
		{
		  ComplexMatrix TmpMatrix(NbrInputStates, NbrStatePerC4Sector[i]);
		  int Tmp = 0;
		  for (int j = 0; j < NbrInputStates; ++j)
		    {
		      if (C4Values[j] == i)
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
		      sprintf (MatrixOutputName, "%s_c4_%d_kx_%d_ky_%d.mat.txt", OutputNamePrefix, i, XMomentum, YMomentum);
		      TmpMatrix.WriteAsciiMatrix(MatrixOutputName);
		    }
		  else
		    {
		      sprintf (MatrixOutputName, "%s_c4_%d_kx_%d_ky_%d.mat", OutputNamePrefix, i, XMomentum, YMomentum);
		      TmpMatrix.WriteMatrix(MatrixOutputName);
		    }
		  delete[] MatrixOutputName;
		}
	    }
	}
      delete[] TmpVectors;
    }
  File.close();
  return 0;

}
