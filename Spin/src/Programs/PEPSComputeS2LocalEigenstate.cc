#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Operator/SpinS2Operator.h"

#include "HilbertSpace/PEPSLocalPhysicalAndVirtualSpin.h"


#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/GenericRealMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Tools/SpinFiles/SpinFileTools.h"

#include "Options/Options.h"

#include "GeneralTools/Endian.h"
#include "GeneralTools/StringTools.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;
using std::ios;




int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("PEPSComputeS2LocalEigenstate" , "0.01");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PEPSComputeS2LocalEigenstate -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int PhysicalSpinValue = 2;
  int Sz = 0;
  int NbrRepresentation = 1;
  int * TableOfRepresentations = new int [NbrRepresentation];
  TableOfRepresentations[0]=1;
//  TableOfRepresentations[0]=1;
  
  
  PEPSLocalPhysicalAndVirtualSpin Space (PhysicalSpinValue,NbrRepresentation,TableOfRepresentations,Sz);
  Space.PrintConversionTable();
  for(int i =0; i < Space.GetHilbertSpaceDimension();i++)
    Space.PrintState(cout,i);

  int NbrStates = Space.GetHilbertSpaceDimension();
  HermitianMatrix S2Matrix(NbrStates, true);
  SpinS2Operator TmpOperator(&Space, 5);
  TmpOperator.GetOperator(S2Matrix) ;

  RealDiagonalMatrix TmpS2Eigenvalues(NbrStates);
  ComplexMatrix TmpBasis (NbrStates, NbrStates);
  TmpBasis.SetToIdentity();
  S2Matrix.LapackDiagonalize(TmpS2Eigenvalues, TmpBasis);
  int NbrS0State = 0;
  for (int i = 0; i < NbrStates; ++i)
    {
      double TmpS2 = TmpS2Eigenvalues[i];
      cout << "<S^2>=" << TmpS2 << " <S>=" << (0.5 * (sqrt((4.0 * TmpS2) + 1.0) - 1.0)) << endl;
      cout << "round(<2S>)=" <<  round(sqrt((4.0 * TmpS2) + 1.0) - 1.0) << endl; 
      if ( fabs(TmpS2) < 1e-9  ) 
	NbrS0State++;
    }
  
  
  ComplexMatrix C2Matrix(NbrS0State,NbrS0State,true);
  ComplexVector TmpVector(Space.GetHilbertSpaceDimension(),true);
  ComplexVector TmpVector1(Space.GetHilbertSpaceDimension(),true);
  for(int i =0; i <NbrS0State; i++)
    {
      Space.ApplyC4Rotation(TmpBasis[i], TmpVector);
      Space.ApplyC4Rotation(TmpVector, TmpVector1); 
      for(int j=0; j<NbrS0State; j++)
	{
	  C2Matrix.SetMatrixElement(i,j, TmpVector1*   TmpBasis[j]);
	}
      TmpVector.ClearVector();
      TmpVector1.ClearVector();
    }
  cout <<C2Matrix<<endl;

  ComplexDiagonalMatrix TmpC2Eigenvalues(NbrS0State);
  ComplexMatrix TmpBasisC2 (NbrS0State,NbrS0State);
  TmpBasisC2.SetToIdentity();
  C2Matrix.LapackDiagonalize(TmpC2Eigenvalues, TmpBasisC2);
  int Nbr1Eigenvalue=0;
  for(int i =0; i <NbrS0State; i++)
    {
      if ( Norm(TmpC2Eigenvalues[i] - 1.0) <1e-8 )
	Nbr1Eigenvalue++;
      cout <<TmpC2Eigenvalues[i]<<" ";
    }

  ComplexMatrix VectorWith1C2 (Space.GetHilbertSpaceDimension(),NbrS0State,true);
  
  cout << TmpBasis<<endl;
  cout <<TmpBasisC2<<endl;
  
  for(int i = 0; i<NbrS0State;i++)
    {
      for(int j = 0; j<NbrS0State;j++)
	{
	  VectorWith1C2[i].AddLinearCombination(TmpBasisC2.GetMatrixElement(j,i), TmpBasis[j]);
	}
    }

  cout <<VectorWith1C2<<endl;

  for(int i =0; i <NbrS0State; i++)
    {
      TmpVector.ClearVector();
      TmpVector1.ClearVector(); 
      Space.ApplyC4Rotation(VectorWith1C2[i], TmpVector);
      Space.ApplyC4Rotation(TmpVector, TmpVector1); 
      cout << TmpVector1 * VectorWith1C2[i]<<endl;
      
    }





  ComplexMatrix C4Matrix(NbrS0State,NbrS0State,true );
  for(int i =0; i <NbrS0State; i++)
    {
      Space.ApplyC4Rotation( TmpBasis[i], TmpVector);
      for(int j=0; j<NbrS0State; j++)
	{
	  C4Matrix.SetMatrixElement(i,j, TmpVector*  TmpBasis[j]);
	}
      TmpVector.ClearVector();
    }
  cout <<C4Matrix<<endl;
  
  ComplexDiagonalMatrix TmpC4Eigenvalues(NbrS0State);
  ComplexMatrix TmpBasisC4 (NbrS0State,NbrS0State);
  TmpBasisC4.SetToIdentity();
  C4Matrix.LapackDiagonalize(TmpC4Eigenvalues, TmpBasisC4);
  
  for(int i =0; i <NbrS0State; i++)
    {
      cout <<TmpC4Eigenvalues[i]<<" ";
    }
  cout <<endl;
  return 0;
}

