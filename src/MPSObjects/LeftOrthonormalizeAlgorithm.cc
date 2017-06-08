#include <iostream>
#include <sys/time.h>

#include "Matrix/ComplexMatrix.h"
#include "LeftOrthonormalizeAlgorithm.h"
#include "CompletelyPositiveMap.h"
#include "LanczosAlgorithm/ComplexArnoldiCPMapsAlgorithm.h"

using std::cout;
using std::endl;

LeftOrthonormalizeAlgorithm::LeftOrthonormalizeAlgorithm (ComplexMatrix * mps, unsigned int physicalDimension, double accuracy, ComplexMatrix & initialGuess)
{
  this->Mps = mps;
  this->PhysicalDimension = physicalDimension;
  this->Accuracy = accuracy;
  this->CenterMatrix = initialGuess;
}


LeftOrthonormalizeAlgorithm::~LeftOrthonormalizeAlgorithm()
{
}

void LeftOrthonormalizeAlgorithm::RunAlgorithm()
{
  ComplexMatrix TmpUnitary (this->CenterMatrix.GetNbrRow(),this->CenterMatrix.GetNbrRow(),true);
  ComplexMatrix TmpUpperTriangular (this->CenterMatrix.GetNbrRow(),this->CenterMatrix.GetNbrRow(),true);
  this->CenterMatrix.QRDecompositionFromLapack (TmpUnitary,  TmpUpperTriangular);
  TmpUpperTriangular/= TmpUpperTriangular.FrobeniusNorm();
  
  for(int i = 0; i < PhysicalDimension; i++)
    this->MpsInLeftForm[i] = TmpUpperTriangular * this->Mps[i];
  
  this->OldCenterMatrix = TmpUpperTriangular;

  ComplexMatrix ACenter (this->MpsInLeftForm[0].GetNbrRow() * PhysicalDimension, this->MpsInLeftForm[0].GetNbrColumn());
  
  for(int i = 0; i < ACenter.GetNbrRow() ; i++)
    {
      for(int j = 0; j <  ACenter.GetNbrColumn() ; j++)
	{
	  ACenter.SetMatrixElement(i,j, this->MpsInLeftForm[i/this->MpsInLeftForm[i].GetNbrRow()].GetMatrixElement(i%this->MpsInLeftForm[i].GetNbrRow(),j));
	}
    }
  
  ComplexMatrix TmpUnitary2 (this->MpsInLeftForm[0].GetNbrRow() * PhysicalDimension, this->MpsInLeftForm[0].GetNbrColumn());
  
  ACenter.QRDecompositionFromLapack (TmpUnitary2, this->CenterMatrix);
  
  for(int i = 0; i < ACenter.GetNbrRow() ; i++)
    {
      for(int j = 0; j <  ACenter.GetNbrColumn() ; j++)
	{
	  this->MpsInLeftForm[i/this->MpsInLeftForm[i].GetNbrRow()].SetMatrixElement(i%this->MpsInLeftForm[i].GetNbrRow(),j , TmpUnitary2.GetMatrixElement(i,j));
	}
    }

  this->Eigenvalue = this->CenterMatrix.FrobeniusNorm();
  this->CenterMatrix/=  this->Eigenvalue;
  double ActualAccuracy = (this->CenterMatrix - this->OldCenterMatrix).FrobeniusNorm();
  while( ActualAccuracy > this->Accuracy )
    {
      CompletelyPositiveMap Map (this->PhysicalDimension, this->Mps, this->MpsInLeftForm);
      ComplexArnoldiCPMapsAlgorithm Arnoldi (&Map,  ActualAccuracy*0.1, 1, 100, true);  
      Arnoldi.InitializeLanczosAlgorithm(this->CenterMatrix);
      Arnoldi.RunLanczosAlgorithm(1);
      while (Arnoldi.TestConvergence() == false)
	{
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  gettimeofday (&(TotalStartingTime), 0);
	  Arnoldi.RunLanczosAlgorithm(1);
	  gettimeofday (&(TotalEndingTime), 0);
	  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	  cout << "iteration done in " << Dt << "s" << endl;
	}
      
      this->CenterMatrix = Arnoldi.GetGroundState();
      
      ComplexMatrix TmpUnitary (this->CenterMatrix.GetNbrRow(),this->CenterMatrix.GetNbrRow(),true);
      ComplexMatrix TmpUpperTriangular (this->CenterMatrix.GetNbrRow(),this->CenterMatrix.GetNbrRow(),true);
      this->CenterMatrix.QRDecompositionFromLapack (TmpUnitary,  TmpUpperTriangular);
      TmpUpperTriangular/= TmpUpperTriangular.FrobeniusNorm();
      
      for(int i = 0; i < PhysicalDimension; i++)
	this->MpsInLeftForm[i] = TmpUpperTriangular * this->Mps[i];
      
      this->OldCenterMatrix = TmpUpperTriangular;
      
  for(int i = 0; i < ACenter.GetNbrRow() ; i++)
    {
      for(int j = 0; j <  ACenter.GetNbrColumn() ; j++)
	{
	  ACenter.SetMatrixElement(i,j, this->MpsInLeftForm[i/this->MpsInLeftForm[i].GetNbrRow()].GetMatrixElement(i%this->MpsInLeftForm[i].GetNbrRow(),j));
	}
    }
      
  ComplexMatrix TmpUnitary2 ( ACenter.GetNbrRow(),  ACenter.GetNbrColumn());
  
  ACenter.QRDecompositionFromLapack (TmpUnitary2, this->CenterMatrix);
  
  for(int i = 0; i < ACenter.GetNbrRow() ; i++)
    {
      for(int j = 0; j <  ACenter.GetNbrColumn() ; j++)
	{
	  this->MpsInLeftForm[i/this->MpsInLeftForm[i].GetNbrRow()].SetMatrixElement(i%this->MpsInLeftForm[i].GetNbrRow(),j , TmpUnitary2.GetMatrixElement(i,j));
	}
    }
  
      
      this->Eigenvalue = this->CenterMatrix.FrobeniusNorm();
      this->CenterMatrix/=  this->Eigenvalue;
      ActualAccuracy = (this->CenterMatrix - this->OldCenterMatrix).FrobeniusNorm();
    }
}
