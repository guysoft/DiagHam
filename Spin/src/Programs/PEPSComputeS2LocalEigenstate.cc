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


  int PhysicalSpinValue = 1;
  int Sz = 0;
  int NbrRepresentation = 2;
  int * TableOfRepresentations = new int [NbrRepresentation];
  TableOfRepresentations[1]=0;
  TableOfRepresentations[0]=1;
  
  
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
      if ( fabs(TmpS2) < 1e-7  ) 
	NbrS0State++;
    }
  
  ComplexMatrix C4Matrix(NbrS0State,NbrS0State,true);
  ComplexVector Tmp(Space.GetHilbertSpaceDimension(),true);

  for(int i =0; i <NbrS0State; i++)
    {
      Space.ApplyC4Rotation(TmpBasis[i], Tmp);
      for(int j=0; j<NbrS0State; j++)
	{
	  C4Matrix.SetMatrixElement(i,j, Tmp*   TmpBasis[j]);
	}
      Tmp.ClearVector();
    }
  
  ComplexDiagonalMatrix TmpC4Eigenvalues(NbrS0State);
  ComplexMatrix TmpBasisC4 (NbrS0State,NbrS0State);
  TmpBasisC4.SetToIdentity();
  C4Matrix.LapackDiagonalize(TmpC4Eigenvalues, TmpBasisC4);
  
  for(int i =0; i <NbrS0State; i++)
    {
      cout <<TmpC4Eigenvalues[i]<<" ";
    }
  cout <<endl;

  cout <<"END C4 eigenvalue"<<endl;
  
  
  HermitianMatrix C2Matrix(NbrS0State,true);
  ComplexVector TmpVector(Space.GetHilbertSpaceDimension(),true);
  ComplexVector TmpVector1(Space.GetHilbertSpaceDimension(),true);
  for(int i =0; i <NbrS0State; i++)
    {
      Space.ApplyC4Rotation(TmpBasis[i], TmpVector);
      Space.ApplyC4Rotation(TmpVector, TmpVector1); 
      for(int j=i; j<NbrS0State; j++)
	{
	  C2Matrix.SetMatrixElement(i,j, TmpVector1*   TmpBasis[j]);
	}
      TmpVector.ClearVector();
      TmpVector1.ClearVector();
    }
  cout <<C2Matrix<<endl;

  RealDiagonalMatrix TmpC2Eigenvalues(NbrS0State);
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
  cout <<endl;

    
  ComplexMatrix VectorWith1C2 (Space.GetHilbertSpaceDimension(), Nbr1Eigenvalue,true);
  ComplexMatrix VectorWithMinus1C2 (Space.GetHilbertSpaceDimension(), NbrS0State-Nbr1Eigenvalue,true);
  //  cout << TmpBasis<<endl;
  //  cout <<TmpBasisC2<<endl;
  int IndexPlus1=0;
  int IndexMinus1=0;
  cout <<" Nbr C2 = 1 "<<  Nbr1Eigenvalue<<endl;
  for(int i = 0; i< NbrS0State;i++)
    {
      if ( Norm(TmpC2Eigenvalues[i] - 1.0) <1e-8 )
	{
	  for(int j = 0; j<NbrS0State;j++)
	    {
	      VectorWith1C2[IndexPlus1].AddLinearCombination(TmpBasisC2.GetMatrixElement(j,i), TmpBasis[j]);
	    }
	  IndexPlus1++;
	}
      else
	{
	  for(int j = 0; j<NbrS0State;j++)
	    {
	      VectorWithMinus1C2[IndexMinus1].AddLinearCombination(TmpBasisC2.GetMatrixElement(j,i), TmpBasis[j]);
	    }
	  IndexMinus1++;
	}
    }
  cout <<" IndexPlus1 "<<  IndexPlus1<<endl;
  cout <<" IndexMinusPlus1 "<<  IndexMinus1<<endl;

  cout <<  VectorWith1C2<<endl;
  if ( Nbr1Eigenvalue > 0 ) 
    {
      cout <<"Check C2"<<endl;

      for(int i =0; i < Nbr1Eigenvalue; i++)
	{
	  Space.ApplyC4Rotation(VectorWith1C2[i], TmpVector);
	  Space.ApplyC4Rotation(TmpVector, TmpVector1); 
	  cout << TmpVector1 *VectorWith1C2[i]<<endl;
	  TmpVector.ClearVector();
	  TmpVector1.ClearVector();
	}

      HermitianMatrix C4Matrix1C2(Nbr1Eigenvalue,true );
      for(int i =0; i < Nbr1Eigenvalue; i++)
	{
	  Space.ApplyC4Rotation(VectorWith1C2[i], TmpVector);
	  for(int j=i; j < Nbr1Eigenvalue ; j++)
	    {
	      C4Matrix1C2.SetMatrixElement(i,j, TmpVector*  VectorWith1C2[j]);
	    }
	  TmpVector.ClearVector();
	}   
      cout << C4Matrix1C2<<endl;

      RealDiagonalMatrix TmpC4Eigenvalues(Nbr1Eigenvalue);
      ComplexMatrix TmpBasisC4 (Nbr1Eigenvalue, Nbr1Eigenvalue);
      TmpBasisC4.SetToIdentity();
      C4Matrix1C2.LapackDiagonalize(TmpC4Eigenvalues, TmpBasisC4);
      
      int Nbr1EigenvalueC4=0;
      for(int i =0; i < Nbr1Eigenvalue; i++)
	{
	  if ( Norm(TmpC4Eigenvalues[i] - 1.0) <1e-8 )
	    Nbr1EigenvalueC4++;
	  cout <<TmpC4Eigenvalues[i]<<" ";

	}
      cout <<endl;
      
      int NbrMinus1C41C2 =  Nbr1Eigenvalue -Nbr1EigenvalueC4;
      
      ComplexMatrix VectorWith1C21C4 (Space.GetHilbertSpaceDimension(),Nbr1EigenvalueC4,true); 
      ComplexMatrix VectorWith1C2Minus1C4 (Space.GetHilbertSpaceDimension(),  NbrMinus1C41C2 ,true);
      
      int IndexPlus1C4=0;
      int IndexMinus1C4=0;
      for(int i = 0; i<  Nbr1Eigenvalue;i++)
	{
	  if ( Norm(TmpC4Eigenvalues[i] - 1.0) <1e-8 )
	    {
	      for(int j = 0; j< Nbr1Eigenvalue;j++)
		{
		  VectorWith1C21C4[IndexPlus1C4].AddLinearCombination(TmpBasisC4.GetMatrixElement(j,i), VectorWith1C2[j]);
		}
	      IndexPlus1C4++;
	    }
	  else
	    {
	      for(int j = 0; j< Nbr1Eigenvalue;j++)
		{
		  VectorWith1C2Minus1C4[IndexMinus1C4].AddLinearCombination(TmpBasisC4.GetMatrixElement(j,i), VectorWith1C2[j]);
		}
	      IndexMinus1C4++;
	    }
	}
      
      
      
      if (Nbr1EigenvalueC4 > 0 ) 
	{  
	  HermitianMatrix HReflexionMatrix(Nbr1EigenvalueC4,true );
	  for(int i =0; i <  Nbr1EigenvalueC4; i++)
	    {
	      Space.ApplyHReflexion( VectorWith1C21C4[i], TmpVector);
	      for(int j=i; j < Nbr1EigenvalueC4 ; j++)
		{
		  HReflexionMatrix.SetMatrixElement(i,j, TmpVector*   VectorWith1C21C4[j]);
		}
	      TmpVector.ClearVector();
	    }
	  
	  RealDiagonalMatrix TmpHREigenvalues(Nbr1EigenvalueC4);
	  ComplexMatrix TmpBasisHR (Nbr1EigenvalueC4,Nbr1EigenvalueC4);
	  TmpBasisHR.SetToIdentity();
	  HReflexionMatrix.LapackDiagonalize(TmpHREigenvalues, TmpBasisHR);
	  

	  int Nbr1EigenvalueHR = 0;
	  
	  for(int i =0; i < Nbr1EigenvalueC4; i++)
	    {
	      if ( Norm(TmpHREigenvalues[i] - 1.0) <1e-8 )
		Nbr1EigenvalueHR++;
	    }
	  
	  int NbrMinus1EigenvalueHR = Nbr1EigenvalueC4 -  Nbr1EigenvalueHR;


	  ComplexMatrix VectorInA1 (Space.GetHilbertSpaceDimension(), Nbr1EigenvalueHR,true);
	  ComplexMatrix VectorInA2 (Space.GetHilbertSpaceDimension(), NbrMinus1EigenvalueHR,true);
	  int IndexA1 = 0;
	  int IndexA2 = 0;
	  for(int i = 0; i<  Nbr1EigenvalueC4;i++)
	    {
	      if ( Norm(TmpHREigenvalues[i] - 1.0) <1e-8 )
		{
		  for(int j = 0; j<  Nbr1EigenvalueC4; j++)
		    {
		      VectorInA1[IndexA1].AddLinearCombination(TmpBasisHR.GetMatrixElement(j,i), VectorWith1C21C4[j]);
		    }
		  IndexA1++;
		}
	      else
		{
		  for(int j = 0; j<   Nbr1EigenvalueC4 ; j++)
		    {
		      VectorInA2[IndexA2].AddLinearCombination(TmpBasisHR.GetMatrixElement(j,i),  VectorWith1C21C4[j]);
		    }
		  IndexA2++;
		}
	    }
	  

	      
	  if ( Nbr1EigenvalueHR > 0 ) 
	    {

	      cout <<"C2 =1 C4 = 1 Sigma H = 1  A1 representation" <<  Nbr1EigenvalueHR<<endl;
	      
	      cout <<"Check SigmaV" <<endl;
	      for(int i = 0; i<  Nbr1EigenvalueHR;i++)
		{
		  Space.ApplyVReflexion( VectorInA1[i], TmpVector);
		  cout << TmpVector *  VectorInA1[i]<<endl;
		  TmpVector.ClearVector();		  
		}
	      
	      cout <<"Check Sigma DX" <<endl;
	      for(int i = 0; i<  Nbr1EigenvalueHR;i++)
		{
		  Space.ApplyDXReflexion( VectorInA1[i], TmpVector);
		  cout << TmpVector *  VectorInA1[i]<<endl;
		  TmpVector.ClearVector();		  
		}
	      cout <<"Check Sigma Minus DX" <<endl;
	      for(int i = 0; i<  Nbr1EigenvalueHR;i++)
		{
		  Space. ApplyDMinusXReflexion( VectorInA1[i], TmpVector);
		  cout << TmpVector *  VectorInA1[i]<<endl;
		  TmpVector.ClearVector(); 	  
		}	      
	    }
	
	  if (  NbrMinus1EigenvalueHR > 0 ) 
	    {

	      cout <<"C2 =1 C4 = 1 Sigma H = -1  A2 representation" 	<<  NbrMinus1EigenvalueHR<<endl;
	      
	      cout <<"Check SigmaV" <<endl;
	      for(int i = 0; i<   NbrMinus1EigenvalueHR;i++)
		{
		  Space.ApplyVReflexion( VectorInA2[i], TmpVector);
		  cout << TmpVector *  VectorInA2[i]<<endl;
		  TmpVector.ClearVector();		  
		}
	      
	      cout <<"Check Sigma DX" <<endl;
	      for(int i = 0; i<   NbrMinus1EigenvalueHR;i++)
		{
		  Space.ApplyDXReflexion( VectorInA2[i], TmpVector);
		  cout << TmpVector *  VectorInA2[i]<<endl;
		  TmpVector.ClearVector();		  
		}
	      cout <<"Check Sigma Minus DX" <<endl;
	      for(int i = 0; i<   NbrMinus1EigenvalueHR;i++)
		{
		  Space. ApplyDMinusXReflexion( VectorInA2[i], TmpVector);
		  cout << TmpVector *  VectorInA2[i]<<endl;
		  TmpVector.ClearVector();		  
		}
	    }	  
	}
          
      if (NbrMinus1C41C2 > 0 ) 
	{  
	 
	  HermitianMatrix HReflexionMatrix(NbrMinus1C41C2,true );
	  for(int i =0; i <  NbrMinus1C41C2; i++)
	    {
	      Space.ApplyHReflexion(VectorWith1C2Minus1C4[i], TmpVector);
	      for(int j=i; j <  NbrMinus1C41C2 ; j++)
		{
		  HReflexionMatrix.SetMatrixElement(i,j, TmpVector* VectorWith1C2Minus1C4[j]);
		}
	      TmpVector.ClearVector();
	    }
	  
	  cout << HReflexionMatrix<<endl;
	  
	  RealDiagonalMatrix TmpHREigenvalues(NbrMinus1C41C2);
	  ComplexMatrix TmpBasisHR (NbrMinus1C41C2,NbrMinus1C41C2);
	  TmpBasisHR.SetToIdentity();
	  HReflexionMatrix.LapackDiagonalize(TmpHREigenvalues, TmpBasisHR);

	  int Nbr1EigenvalueHR = 0;
	  
	  for(int i =0; i < NbrMinus1C41C2; i++)
	    {
	      if ( Norm(TmpHREigenvalues[i] - 1.0) <1e-8 )
		Nbr1EigenvalueHR++;
	    }
	  
	  int NbrMinus1EigenvalueHR = NbrMinus1C41C2 -  Nbr1EigenvalueHR;
	  


	  ComplexMatrix VectorInB1 (Space.GetHilbertSpaceDimension(), Nbr1EigenvalueHR,true);
	  ComplexMatrix VectorInB2 (Space.GetHilbertSpaceDimension(), NbrMinus1EigenvalueHR,true);
	  int IndexB1 = 0;
	  int IndexB2 = 0;
	  for(int i = 0; i<  Nbr1EigenvalueC4;i++)
	    {
	      if ( Norm(TmpHREigenvalues[i] - 1.0) <1e-8 )
		{
		  for(int j = 0; j<  Nbr1EigenvalueC4; j++)
		    {
		      VectorInB1[IndexB1].AddLinearCombination(TmpBasisHR.GetMatrixElement(j,i),  VectorWith1C2Minus1C4[j]);
		    }
		  IndexB1++;
		}
	      else
		{
		  for(int j = 0; j<   Nbr1EigenvalueC4 ; j++)
		    {
		      VectorInB2[IndexB2].AddLinearCombination(TmpBasisHR.GetMatrixElement(j,i),  VectorWith1C2Minus1C4[j]);
		    }
		  IndexB2++;
		}
	    }
	  




	  if ( Nbr1EigenvalueHR > 0 ) 
	    {
	      cout <<"C2 =1 C4 = -1 Sigma H = 1  B1 representation" <<  Nbr1EigenvalueHR<<endl;
	      
	      cout <<"Check SigmaV" <<endl;
	      for(int i = 0; i<  Nbr1EigenvalueHR;i++)
		{
		  Space.ApplyVReflexion( VectorInB1[i], TmpVector);
		  cout << TmpVector *  VectorInB1[i]<<endl;
		  TmpVector.ClearVector();		  
		}
	      
	      cout <<"Check Sigma DX" <<endl;
	      for(int i = 0; i<  Nbr1EigenvalueHR;i++)
		{
		  Space.ApplyDXReflexion( VectorInB1[i], TmpVector);
		  cout << TmpVector *  VectorInB1[i]<<endl;
		  TmpVector.ClearVector();		  
		}
	      cout <<"Check Sigma Minus DX" <<endl;
	      for(int i = 0; i<  Nbr1EigenvalueHR;i++)
		{
		  Space. ApplyDMinusXReflexion( VectorInB1[i], TmpVector);
		  cout << TmpVector *  VectorInB1[i]<<endl;
		  TmpVector.ClearVector();		  
		}	      
	    }
	  if (  NbrMinus1EigenvalueHR > 0 ) 
	    {
	      cout <<"C2 =1 C4 = -1 Sigma H = -1  B2 representation" 	<<  NbrMinus1EigenvalueHR<<endl;
	      
	      cout <<"Check SigmaV" <<endl;
	      for(int i = 0; i<   NbrMinus1EigenvalueHR;i++)
		{
		  Space.ApplyVReflexion( VectorInB2[i], TmpVector);
		  cout << TmpVector *  VectorInB2[i]<<endl;
		  TmpVector.ClearVector();		  
		}
	      
	      cout <<"Check Sigma DX" <<endl;
	      for(int i = 0; i<   NbrMinus1EigenvalueHR;i++)
		{
		  Space.ApplyDXReflexion( VectorInB2[i], TmpVector);
		  cout << TmpVector *  VectorInB2[i]<<endl;
		  TmpVector.ClearVector();		  
		}
	      cout <<"Check Sigma Minus DX" <<endl;
	      for(int i = 0; i<   NbrMinus1EigenvalueHR;i++)
		{
		  Space. ApplyDMinusXReflexion( VectorInB2[i], TmpVector);
		  cout << TmpVector *  VectorInB2[i]<<endl;
		  TmpVector.ClearVector();		  
		}
	    }
	  
	}
    }
  
  return 0;
}

