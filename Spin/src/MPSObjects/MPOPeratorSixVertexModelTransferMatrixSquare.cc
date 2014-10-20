#include "MPOPeratorSixVertexModelTransferMatrixSquare.h"

using std::cout;
using std::endl;



MPOPeratorSixVertexModelTransferMatrixSquare::MPOPeratorSixVertexModelTransferMatrixSquare(int nbrSites)
{
  this->PhysicalDimension = 2;
  this->MPOBondDimension = 4;
  this->InitializeTensorsElements();
  this->NbrSites = nbrSites;

  
  this->LeftVector = new double[4];
  this->RightVector = new double[4];
  this->RightVector[0] = 0;  this->RightVector[1] = 1.0;  this->RightVector[2] = -1.0; this->RightVector[3] = 0;
  this->LeftVector[0] = 0;  this->LeftVector[1] = - 1.0;  this->LeftVector[2] = 1.0; this->LeftVector[3] = 0;
}



MPOPeratorSixVertexModelTransferMatrixSquare::~MPOPeratorSixVertexModelTransferMatrixSquare()
{
  delete LeftVector;
  delete RightVector;
  delete [] ElementsValues;
  delete [] IndexValues;    
}


void MPOPeratorSixVertexModelTransferMatrixSquare::InitializeTensorsElements()
{
  double T[this->PhysicalDimension][this->PhysicalDimension][this->PhysicalDimension][this->PhysicalDimension];
  
  for (int j = 0; j < this->PhysicalDimension ; j++)
    {
      for (int k=0; k < this->PhysicalDimension; k++)
	{
	  for (int l=0; l <  this->PhysicalDimension; l++)
	    {
	      for (int m = 0;  m <  this->PhysicalDimension; m++)
		T[j][k][l][m] = 0;
	    }
	}
    }
  
  T[1][1][1][1] = 0.5;
  T[0][0][0][0] = 0.5;
  T[1][0][1][0] = -0.5;
  T[0][1][0][1] = -0.5;
  T[1][1][0][0] = 1;
  T[0][0][1][1] = 1;
  
  this->NbrNonZeroElements = 0;
  for (int j = 0; j <  this->MPOBondDimension; j++)
    {
      for (int k = 0; k < this->PhysicalDimension; k++)
	{
	  for (int l = 0; l< this->MPOBondDimension ; l++)
	    {
	      for (int m=0; m< this->PhysicalDimension; m++)
		{
		  double Tmp = 0.0;
		  for (int p=0; p < this->PhysicalDimension; p++)
		    {
		      Tmp += T[j% this->PhysicalDimension][k][l% this->PhysicalDimension][p] * T[j/ this->PhysicalDimension][p][l/ this->PhysicalDimension][m];
		    }
		  if (Tmp != 0.0)
		    {
		      this->NbrNonZeroElements++;
		    }
		}
	    }
	}
    }
  cout <<"Nbr Non zero Elements = "<<  this->NbrNonZeroElements<<endl;
  this->ElementsValues = new double [this->NbrNonZeroElements];
  this->IndexValues = new unsigned int[this->NbrNonZeroElements];

  this->NbrNonZeroElements=0;  
  for (int j = 0; j<  this->MPOBondDimension; j++)
    {
      for (int k = 0; k < this->PhysicalDimension; k++)
	{
	  for (int l = 0; l< this->MPOBondDimension ; l++)
	    {
	      for (int m = 0;  m < this->PhysicalDimension; m++)
		{
		  double Tmp = 0.0;
		  for (int p=0; p < this->PhysicalDimension; p++)
		    {
		      Tmp += T[j% this->PhysicalDimension][k][l% this->PhysicalDimension][p] * T[j/ this->PhysicalDimension][p][l/ this->PhysicalDimension][m];
		    }
		  if (Tmp != 0.0)
		    {
		      this->ElementsValues[this->NbrNonZeroElements]=Tmp;
		      this->IndexValues[this->NbrNonZeroElements] = GetTensorIndexFromAllIndices(m, k, l,j);
		      this->NbrNonZeroElements++;
		    }
		}
	    }
	}
    }
   
   cout <<"Nbr Non zero Elements = "<<  this->NbrNonZeroElements<<endl;
   
   /*   for (int j = 0; j < this->PhysicalDimension ; j++)
     {
       for (int k=0; k < this->PhysicalDimension; k++)
	 {
	   for (int l=0; l <  this->PhysicalDimension; l++)
	     {
	       delete T[j][k][l];
	    }
	   	       delete T[j][k];
	}
       delete T[j];
       } */
}


void MPOPeratorSixVertexModelTransferMatrixSquare::PrintTensorElements()
{
  cout <<"#Tensor index indexDown indexUp indexLeft indexRight Check Index Values" <<endl;
  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      int MPOIndiceDown = this->GetIndiceDownFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceLeft = this->GetIndiceLeftFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceUp = this->GetIndiceUpFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceRight = this->GetIndiceRightFromTensorIndex(this->IndexValues[i]);
      
      int Tmp = GetTensorIndexFromAllIndices( MPOIndiceDown,  MPOIndiceUp,  MPOIndiceLeft,  MPOIndiceRight);
      cout << this->IndexValues[i] <<" "<<MPOIndiceDown<< " "<< MPOIndiceUp<< " "<< MPOIndiceLeft<< " "<<MPOIndiceRight<<" " <<Tmp<<" "<< this->ElementsValues[i]<<endl;
    }
}


void MPOPeratorSixVertexModelTransferMatrixSquare::ComputeL(Tensor3<double> & L)
{
  cout <<" I am here"<<endl;;
  if (this->Site->GetSitePosition() == 0)
    {
      int BondDimensionRight = this->Site->GetBondDimensionRight();
      int BondDimensionLeft = this->Site->GetBondDimensionLeft();
      RealMatrix * M = this->Site->GetM();
      cout << M[0] << M[1];
      for (int i = 0; i < this->NbrNonZeroElements; i++)
	{
	  int MPOIndiceDown = this->GetIndiceDownFromTensorIndex(this->IndexValues[i]);
	  int MPOIndiceLeft = this->GetIndiceLeftFromTensorIndex(this->IndexValues[i]);
	  int MPOIndiceUp = this->GetIndiceUpFromTensorIndex(this->IndexValues[i]);
	  int MPOIndiceRight = this->GetIndiceRightFromTensorIndex(this->IndexValues[i]);
	  
	  for (int RightA = 0;RightA < this->Site->GetBondDimensionRight() ; RightA++ )
	    {
	      for (int RightC = 0;RightC < this->Site->GetBondDimensionRight() ; RightC++ )
		{
		  L(RightA, MPOIndiceRight,RightC) +=  M[MPOIndiceUp](0,RightC) * this->ElementsValues[i]*this->LeftVector[MPOIndiceLeft] * M[MPOIndiceDown](0,RightA);
		}
	    }
	}
      cout <<" I have finished"<<endl;
    }
  else
    AbstractMPOperatorOBC::ComputeL(L);
}


void MPOPeratorSixVertexModelTransferMatrixSquare::ComputeR(Tensor3<double> & R)
{
  if (this->Site->GetSitePosition() == this->NbrSites - 1)
    {
      RealMatrix * M = this->Site->GetM();
      int BondDimensionRight = this->Site->GetBondDimensionRight();
      int BondDimensionLeft = this->Site->GetBondDimensionLeft();
      for (int i = 0; i < this->NbrNonZeroElements; i++)
	{
	  int MPOIndiceDown = this->GetIndiceDownFromTensorIndex(this->IndexValues[i]);
	  int MPOIndiceLeft = this->GetIndiceLeftFromTensorIndex(this->IndexValues[i]);
	  int MPOIndiceUp = this->GetIndiceUpFromTensorIndex(this->IndexValues[i]);
	  int MPOIndiceRight = this->GetIndiceRightFromTensorIndex(this->IndexValues[i]);
	  
	  for (int LeftA = 0; LeftA < this->Site->GetBondDimensionLeft() ; LeftA++ )
	    {
	      for (int LeftC = 0;LeftC < this->Site->GetBondDimensionLeft() ; LeftC++ )
		{
		  R(LeftA, MPOIndiceLeft,LeftC) +=  M[MPOIndiceUp](LeftC,0) * this->ElementsValues[i]*this->RightVector[MPOIndiceRight] * M[MPOIndiceDown](LeftA,0);
		}
	    }
	}
    }
  else
    AbstractMPOperatorOBC::ComputeR(R);
}

