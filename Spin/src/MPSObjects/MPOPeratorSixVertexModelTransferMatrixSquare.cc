#include "MPOPeratorSixVertexModelTransferMatrixSquare.h"

using std::cout;
using std::endl;



MPOPeratorSixVertexModelTransferMatrixSquare::MPOPeratorSixVertexModelTransferMatrixSquare()
{
  this->LeftVector = 0;
  this->RightVector = 0; 
  this->ElementsValues = 0;
  this->IndexValues = 0;
}



MPOPeratorSixVertexModelTransferMatrixSquare::MPOPeratorSixVertexModelTransferMatrixSquare(int nbrSites)
{
  this->PhysicalDimension = 2;
  this->MPOBondDimension = 4; 
  this->InitializeTensorsElements();
  this->NbrSites = nbrSites;

  
  this->LeftVector = new double[this->MPOBondDimension];
  this->RightVector = new double[this->MPOBondDimension];
  this->RightVector[0] = 0;  this->RightVector[1] = 1.0;  this->RightVector[2] = -1.0; this->RightVector[3] = 0;
  this->LeftVector[0] = 0;  this->LeftVector[1] = -1.0;  this->LeftVector[2] = 1.0; this->LeftVector[3] = 0;
}



MPOPeratorSixVertexModelTransferMatrixSquare::~MPOPeratorSixVertexModelTransferMatrixSquare()
{
  delete [] LeftVector;
  delete [] RightVector;
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
}


void MPOPeratorSixVertexModelTransferMatrixSquare::PrintTensorElements()
{
  cout <<"#Tensor index indexDown indexUp indexLeft indexRight Check Index Values" <<endl;
  unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);      
      int Tmp = GetTensorIndexFromAllIndices( MPOIndiceDown,  MPOIndiceUp,  MPOIndiceLeft,  MPOIndiceRight);
      cout << this->IndexValues[i] <<" "<<MPOIndiceDown<< " "<< MPOIndiceUp<< " "<< MPOIndiceLeft<< " "<<MPOIndiceRight<<" " <<Tmp<<" "<< this->ElementsValues[i]<<endl;
    }
}


void MPOPeratorSixVertexModelTransferMatrixSquare::ComputeL(Tensor3<double> & L)
{
  if (this->Site->GetSitePosition() == 0)
    {
      int BondDimensionRight = this->Site->GetBondDimensionRight();
      RealMatrix * M = this->Site->GetM();
      unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
      for (int i = 0; i < this->NbrNonZeroElements; i++)
	{
          this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);      
		  
	      for (int RightC = 0;RightC < this->Site->GetBondDimensionRight() ; RightC++ )
		{
	  for (int RightA = 0;RightA < this->Site->GetBondDimensionRight() ; RightA++ )
	    {

		  L(RightA, MPOIndiceRight,RightC) +=  M[MPOIndiceUp](0,RightC) * this->ElementsValues[i] * this->LeftVector[MPOIndiceLeft] * M[MPOIndiceDown](0,RightA);
		}
	    }
	}
    }
  else
    AbstractMPOperatorOBC::ComputeL(L);
}


void MPOPeratorSixVertexModelTransferMatrixSquare::ComputeR(Tensor3<double> & R)
{
  if (this->Site->GetSitePosition() == this->NbrSites - 1)
    {
      RealMatrix * M = this->Site->GetM();

      int BondDimensionLeft = this->Site->GetBondDimensionLeft();
      unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
      for (int i = 0; i < this->NbrNonZeroElements; i++)
	{
           this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);      
		  
	      for (int LeftC = 0;LeftC < this->Site->GetBondDimensionLeft() ; LeftC++ )
		{
	  for (int LeftA = 0; LeftA < this->Site->GetBondDimensionLeft() ; LeftA++ )
	    {
		  R(LeftA, MPOIndiceLeft,LeftC) +=  M[MPOIndiceUp](LeftC,0) * this->ElementsValues[i]*this->RightVector[MPOIndiceRight] * M[MPOIndiceDown](LeftA,0);
		}
	    }
	}
    }
  else
    AbstractMPOperatorOBC::ComputeR(R);
}

RealVector& MPOPeratorSixVertexModelTransferMatrixSquare::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent)
{
 if(this->IDMRGFlag)
  return this->LowLevelMultiplyTwoSites(vSource,vDestination, firstComponent, nbrComponent); 
else 
 return this->LowLevelMultiplyOneSite(vSource,vDestination, firstComponent, nbrComponent);

}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& MPOPeratorSixVertexModelTransferMatrixSquare::LowLevelMultiplyOneSite(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent)
{
  vDestination.ClearVector();
  if (this->Site->GetSitePosition() == 0)
  {
    int BondDimensionRight = this->Site->GetBondDimensionRight(); 
    Tensor3<double> & RightR = this->Site->GetNextR();
    Tensor3<double> * B = new Tensor3<double>[this->PhysicalDimension];

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
        B[i] = Tensor3<double>(this->MPOBondDimension,BondDimensionRight,1,true);
  for (int RightB = 0; RightB < this->MPOBondDimension ; RightB++)
    {
	  for(int RightC = 0; RightC < BondDimensionRight; RightC++)
	    {
	      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		{
		  B[i](RightB,RightC,0) +=  RightR(RightA,RightB,RightC) * vSource[(long int)BondDimensionRight*i + RightA]; 
		}
	}
     }
   }

      unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
 for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);      
      for (int NewRight = 0;  NewRight < BondDimensionRight;  NewRight++)
	{
            vDestination[(long int)BondDimensionRight*MPOIndiceUp + NewRight] +=  B[MPOIndiceDown](MPOIndiceRight,NewRight,0) *  this->ElementsValues[i] * this->LeftVector[MPOIndiceLeft];
        }
   }

  delete [] B;
  return  vDestination;
 
  }

   if (this->Site->GetSitePosition() == this->NbrSites - 1)
   {
    int BondDimensionLeft = this->Site->GetBondDimensionLeft(); 
    Tensor3<double> & LeftL = this->Site->GetPreviousL();
    Tensor3<double> * B = new Tensor3<double>[this->PhysicalDimension];

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      B[i] = Tensor3<double>(this->MPOBondDimension,BondDimensionLeft,1,true);

for(int LeftC = 0; LeftC < BondDimensionLeft; LeftC++)
	    {
   for (int LeftB = 0; LeftB < this->MPOBondDimension ; LeftB++)
    {
	  
	      for(int LeftA  = 0;  LeftA < BondDimensionLeft;  LeftA++)
		{
		  B[i](LeftB,LeftC,0) +=  LeftL(LeftA,LeftB,LeftC) * vSource[(long int)BondDimensionLeft*i + LeftA];
		}
	}
     }
    }

 unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
 for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);      

	  for (int NewLeft = 0;  NewLeft < BondDimensionLeft;  NewLeft++)
	    {
              vDestination[(long int)BondDimensionLeft*MPOIndiceUp + NewLeft] +=  B[MPOIndiceDown](MPOIndiceLeft,NewLeft,0) *  this->ElementsValues[i] * this->RightVector[MPOIndiceRight];
            }
    }

    delete [] B;
    return vDestination;
   }

  return AbstractMPOperatorOBC::LowLevelMultiply(vSource,vDestination,firstComponent,nbrComponent);

}
 





// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& MPOPeratorSixVertexModelTransferMatrixSquare::LowLevelMultiplyTwoSites(RealVector& vSource, RealVector& vDestination,  int firstComponent, int nbrComponent)
{

 vDestination.ClearVector();
  if ((this->SiteLeft->GetSitePosition() == 0)&&(this->SiteRight->GetSitePosition() ==  this->NbrSites - 1))
  {

  int SquarePhysicalDimension = this->PhysicalDimension * this->PhysicalDimension;

  Tensor3<double> * B =  new  Tensor3<double>  [SquarePhysicalDimension];
  for (int i = 0; i < SquarePhysicalDimension; i++)
    {
      B[i] = Tensor3<double>(this->MPOBondDimension,1,1,true);
    }

//  int LinearizedPhysicalIndice = PhysicalIndiceLeft +  this->PhysicalDimension *  PhysicalIndiceRight;

    unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
    for(int  PhysicalIndiceRight = 0;  PhysicalIndiceRight <this->PhysicalDimension ; PhysicalIndiceRight++)
{
    for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight); 
      B[MPOIndiceUp +  this->PhysicalDimension *  PhysicalIndiceRight](MPOIndiceRight,0,0) +=  this->ElementsValues[i] * this->LeftVector[MPOIndiceLeft] * vSource[ (long int) MPOIndiceDown +  this->PhysicalDimension *  PhysicalIndiceRight];
   }

}
 for(int PhysicalIndiceLeft= 0;   PhysicalIndiceLeft <this->PhysicalDimension ; PhysicalIndiceLeft++)
{
    for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);
      vDestination[(long) PhysicalIndiceLeft + this->PhysicalDimension *MPOIndiceUp] += this->ElementsValues[i] * this->RightVector[MPOIndiceRight] *  B[PhysicalIndiceLeft +  this->PhysicalDimension *  MPOIndiceDown](MPOIndiceLeft,0,0);
    }
}

 delete [] B;
 return vDestination;

 }

 return AbstractMPOperatorOBC::LowLevelMultiplyTwoSites(vSource,vDestination,firstComponent,nbrComponent);
}
