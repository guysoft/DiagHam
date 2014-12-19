#include "ComplexMPOperatorOBC.h"
#include "Tensor/Tensor3.h"
#include "Matrix/RealSymmetricMatrix.h"
#include <iostream>
#include <sys/time.h>

using std::cout;
using std::endl;

ComplexMPOperatorOBC::ComplexMPOperatorOBC()
{
  this->ElementsValues = 0;
  this->LeftVector = 0;
  this->RightVector = 0;
}

ComplexMPOperatorOBC::~ComplexMPOperatorOBC()
{
}

void ComplexMPOperatorOBC::ComputeLCore(Tensor3<Complex> & L)
{
  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  
  Tensor3<Complex> & LeftL = ((ComplexMPSSite *) this->Site)->GetPreviousL();
  
  Tensor3<Complex> * B =  new  Tensor3<Complex> [this->PhysicalDimension];
  ComplexMatrix * M = ((ComplexMPSSite *) this->Site)->GetM();

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
     B[i] = Tensor3<Complex>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
     for(int LeftC = 0; LeftC < BondDimensionLeft; LeftC++)
     {
     for(int LeftB = 0; LeftB < this->MPOBondDimension; LeftB++)
     {
      for (int RightA = 0; RightA < BondDimensionRight; RightA++)
	{
              Complex & Tmp =  B[i](RightA,LeftB,LeftC);
	      for(int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
		{
		   Tmp +=  LeftL(LeftA,LeftB,LeftC)* M[i].GetMatrixElement(LeftA,RightA);
		}
	    }
           }
	}
    }
  
  Tensor3<Complex> * A = new Tensor3<Complex>[this->PhysicalDimension];

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
  A[i] = Tensor3<Complex>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
    }

  unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
     this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);
     for (int LeftC = 0;  LeftC < BondDimensionLeft;  LeftC++)
     {
      for (int RightA = 0; RightA < BondDimensionRight; RightA++)
	{
	      A[MPOIndiceUp](RightA, MPOIndiceRight,LeftC) +=  B[MPOIndiceDown](RightA,MPOIndiceLeft,LeftC) * this->ElementsValues[i];
	}
	}
    }

  delete [] B;
  

  for (int RightC = 0; RightC <  BondDimensionRight;  RightC++)
  {
    for (int  RightB = 0;  RightB < this->MPOBondDimension;  RightB++)
	   { 
             for (int RightA = 0; RightA < BondDimensionRight; RightA++)
             {
              Complex & Tmp = L(RightA,RightB,RightC);
              for (int i = 0; i < this->PhysicalDimension; i++)
               {
	      for (int LeftC = 0;  LeftC < BondDimensionLeft; LeftC++)
		{ 
		   Tmp += Conj(M[i].GetMatrixElement(LeftC,RightC)) * A[i](RightA,RightB,LeftC);
		 }
		}
	    }
	}
    } 

  delete [] A;
}


void ComplexMPOperatorOBC::ComputeRCore(Tensor3<Complex> & R)
{
  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  Tensor3<Complex> & RightR = ((ComplexMPSSite *) this->Site)->GetNextR();
  
  Tensor3<Complex> * B =  new  Tensor3<Complex>  [this->PhysicalDimension];
  ComplexMatrix * M = ((ComplexMPSSite *) this->Site)->GetM();

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
     B[i] = Tensor3<Complex>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
  for (int RightB = 0; RightB < this->MPOBondDimension ; RightB++)
    {
      for (int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
	{
	  for(int RightC = 0; RightC < BondDimensionRight; RightC++)
	    {
              Complex & Tmp =  B[i](LeftA,RightB,RightC);
	      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		{
		   Tmp +=  RightR(RightA,RightB,RightC) * M[i].GetMatrixElement(LeftA,RightA);
		}
	    }
	}
    }
   }

  Tensor3<Complex> * A = new Tensor3<Complex> [this->PhysicalDimension];
  
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      A[i] = Tensor3<Complex>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
    }
  
  unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);
      
      for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	{
	  for (int RightC = 0;  RightC < BondDimensionRight;  RightC++)
	    {
	      A[MPOIndiceUp](LeftA, MPOIndiceLeft,RightC) +=  B[MPOIndiceDown](LeftA,MPOIndiceRight,RightC) * this->ElementsValues[i];
	    }
	}
    } 

  delete [] B;
  
  for (int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
    {
      for (int LeftC = 0; LeftC <  BondDimensionLeft;  LeftC++)
	{
	  for (int  LeftB = 0;  LeftB < this->MPOBondDimension;  LeftB++)
	    { 
              Complex & Tmp =  R(LeftA,LeftB,LeftC);
                 for (int i = 0; i < this->PhysicalDimension; i++)
		 {
	      for (int RightC = 0;  RightC < BondDimensionRight; RightC++)
		{ 
	      Tmp += Conj(M[i].GetMatrixElement(LeftC,RightC)) * A[i](LeftA,LeftB,RightC);
		    }
		}
	    }
	}
    }

  delete [] A;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ComplexMPOperatorOBC::LowLevelMultiplyCore(ComplexVector& vSource, ComplexVector& vDestination, 
				       int firstComponent, int nbrComponent)
{
//     cout <<"ComplexVector& ComplexMPOperatorOBC::LowLevelMultiplyCore(ComplexVector& vSource, ComplexVector& vDestination,    int firstComponent, int nbrComponent)"<<endl;
/*  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&TotalStartingTime, 0);*/

  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  Tensor3<Complex> & RightR = ((ComplexMPSSite *) this->Site)->GetNextR();
  Tensor3<Complex> & LeftL = ((ComplexMPSSite *) this->Site)->GetPreviousL();
 
  Tensor3<Complex> * B =  new  Tensor3<Complex>  [this->PhysicalDimension];

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      B[i] = Tensor3<Complex>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
      for(int RightC = 0; RightC < BondDimensionRight; RightC++)
	    {
          for (int RightB = 0; RightB < this->MPOBondDimension ; RightB++)
            {

      for (int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
	{
              Complex & Tmp =  B[i](LeftA,RightB,RightC);
	      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		{
		  Tmp +=  RightR(RightA,RightB,RightC) * vSource[(long int)BondDimensionRight*(BondDimensionLeft*i+ LeftA) + RightA];
		}
	    }
	}
    }
   }

/* gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout <<"First Part " << Dt << "s" << endl;
 gettimeofday (&TotalStartingTime, 0);
*/
  Tensor3<Complex> * A =  new  Tensor3<Complex>  [this->PhysicalDimension];
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      A[i] = Tensor3<Complex>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
    }
 
  unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
 for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);

     
      for (int RightC = 0;  RightC < BondDimensionRight;  RightC++)
	   {
       for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	 {
	      A[MPOIndiceUp](LeftA, MPOIndiceLeft,RightC) +=  B[MPOIndiceDown](LeftA,MPOIndiceRight,RightC) * this->ElementsValues[i];
	    }
	}
   } 

  delete [] B;
  int LastComponent = firstComponent + nbrComponent;

/*  gettimeofday (&TotalEndingTime, 0);
  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout <<"Second Part " << Dt << "s" << endl;
 gettimeofday (&TotalStartingTime, 0);*/

  for(int Index =  firstComponent; Index < LastComponent ;Index++)
  {
     Complex & Tmp = vDestination[Index];
  for (int RightB = 0;  RightB < this->MPOBondDimension;  RightB++)
	{
       for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	{
          Tmp +=  LeftL(LeftA,RightB, (Index/BondDimensionRight)%BondDimensionLeft) * A[Index/(BondDimensionRight*BondDimensionLeft)](LeftA, RightB,Index % BondDimensionRight);
        }
    }
  }

/*gettimeofday (&TotalEndingTime, 0);
Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
 cout <<"Third Part " << Dt << "s" << endl;*/
 delete [] A;
 return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ComplexMPOperatorOBC::LowLevelMultiplyTwoSitesCore(ComplexVector& vSource, ComplexVector& vDestination,  int firstComponent, int nbrComponent)
{

//cout <<"ComplexVector& ComplexMPOperatorOBC::LowLevelMultiplyTwoSitesCore(ComplexVector& vSource, ComplexVector& vDestination,  int firstComponent, int nbrComponent)"<<endl;
//cout <<"Vsource" <<vSource<<endl;

/*  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&TotalStartingTime, 0)
;
*/

  int BondDimensionRight = this->SiteRight->GetBondDimensionRight();
  int BondDimensionLeft = this->SiteLeft->GetBondDimensionLeft();

  Tensor3<Complex> & RightR = ((ComplexMPSSite *) this->SiteRight)->GetNextR();
  Tensor3<Complex> & LeftL = ((ComplexMPSSite *) this->SiteLeft)->GetPreviousL();


//  LeftL.PrintTensor();
//  RightR.PrintTensor();

  int SquarePhysicalDimension = this->PhysicalDimension * this->PhysicalDimension;

  Tensor3<Complex> * B =  new  Tensor3<Complex>  [SquarePhysicalDimension];
  for (int i = 0; i < SquarePhysicalDimension; i++)
    {
      B[i] = Tensor3<Complex>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
    }

//  int LinearizedPhysicalIndice = PhysicalIndiceLeft +  this->PhysicalDimension *  PhysicalIndiceRight;
  for (int LinearizedPhysicalIndice= 0 ; LinearizedPhysicalIndice < SquarePhysicalDimension; LinearizedPhysicalIndice++)
  {
  for (int LeftB = 0; LeftB < this->MPOBondDimension ; LeftB++)
    {
	  for(int LeftC = 0; LeftC < BondDimensionLeft; LeftC++)
	    {
	      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		{

           for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	    {
		  B[LinearizedPhysicalIndice](RightA,LeftB,LeftC) +=  LeftL(LeftA,LeftB,LeftC) * vSource[LinearizedPhysicalIndice + SquarePhysicalDimension*(LeftA + RightA*BondDimensionLeft)];
	}
	    }
	}
    }
}

  Tensor3<Complex> * A =  new  Tensor3<Complex>  [SquarePhysicalDimension];
  for (int i = 0; i < SquarePhysicalDimension ; i++)
    {
      A[i] = Tensor3<Complex>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
    }
 
    unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceMiddle,MPOIndiceRight;

 for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceMiddle);
    
 for(int PhysicalIndiceRight= 0; PhysicalIndiceRight <this->PhysicalDimension ; PhysicalIndiceRight++)
     {
      for (int LeftC = 0;  LeftC < BondDimensionLeft;  LeftC++)
	{
	  for (int RightA = 0;  RightA < BondDimensionRight;  RightA++)
	    {
              A[MPOIndiceUp +  this->PhysicalDimension *  PhysicalIndiceRight](LeftC, MPOIndiceMiddle,RightA) += B[MPOIndiceDown +  this->PhysicalDimension *  PhysicalIndiceRight](RightA,MPOIndiceLeft,LeftC) * this->ElementsValues[i];
	    }
	}
    }
   }

  delete [] B;

  Tensor3<Complex> * C =  new  Tensor3<Complex>  [SquarePhysicalDimension];
  for (int i = 0; i < SquarePhysicalDimension ; i++)
    {
      C[i] = Tensor3<Complex>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
    }

  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
     this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceMiddle,  MPOIndiceRight);
     for(int PhysicalIndiceLeft= 0; PhysicalIndiceLeft <this->PhysicalDimension ; PhysicalIndiceLeft++)
     {
      for (int LeftC = 0;  LeftC < BondDimensionLeft;  LeftC++)
	{
	  for (int RightA = 0;  RightA < BondDimensionRight;  RightA++)
	    {
              C[PhysicalIndiceLeft +  this->PhysicalDimension *   MPOIndiceUp](LeftC, MPOIndiceRight,RightA) += A[PhysicalIndiceLeft +  this->PhysicalDimension *  MPOIndiceDown](LeftC,MPOIndiceMiddle,RightA) * this->ElementsValues[i];
	    }
	}
    }
   } 

 delete [] A;

 // Index = LinearizedPhysicalIndice + SquarePhysicalDimension*(LeftA + RightA*BondDimensionLeft))

  int LastComponent = firstComponent + nbrComponent;

/*  gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout <<"First Part " << Dt << "s" << endl;
  gettimeofday (&TotalStartingTime, 0);*/
 
 for(int Index =  firstComponent; Index < LastComponent ;Index++)
  {
       for (int RightA = 0;  RightA < BondDimensionRight;  RightA++)
	{
      for (int RightB = 0;  RightB < this->MPOBondDimension;  RightB++)
	{
            vDestination[Index] +=  RightR(RightA,RightB,Index/(SquarePhysicalDimension*BondDimensionLeft)) *  C[Index%SquarePhysicalDimension](Index/SquarePhysicalDimension%BondDimensionLeft , RightB,RightA);
        }
    }
 }
/*
gettimeofday (&TotalEndingTime, 0);
Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
 cout <<"Second Part " << Dt << "s" << endl;*/
 delete [] C;
 // cout <<" vDestination = "<<  vDestination<<endl;
 return vDestination;
}

// store Hamiltonian into an hermitian matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix
HermitianMatrix& ComplexMPOperatorOBC::GetTwoSitesHamiltonian (HermitianMatrix & M)
{
  ComplexVector TmpV1 (this->GetTwoSitesHilbertSpaceDimension(), true);
  ComplexVector TmpV2 (this->GetTwoSitesHilbertSpaceDimension(), true);
  for (int i = 0; i < this->GetTwoSitesHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->AbstractMPOperatorOBC::LowLevelMultiplyTwoSites(TmpV1, TmpV2);
      if (this->LeftHamiltonianVectorMultiplicationFlag == false)
	{
	  for (int j = i; j < this->GetTwoSitesHilbertSpaceDimension(); j++)
	    {
	      M.SetMatrixElement(i, j, TmpV2[j]);
	    }
	}
      else
	{
	  for (int j = i; j < this->GetTwoSitesHilbertSpaceDimension(); j++)
	    {
	      M.SetMatrixElement(j, i, TmpV2[j]);
	    }
	}
      TmpV1[i] = 0.0;
    }
  return M;  
}



void ComplexMPOperatorOBC::ComputeL(Tensor3<Complex> & L)
{
  if (this->Site->GetSitePosition() == 0)
    {
      int BondDimensionRight = this->Site->GetBondDimensionRight();
      ComplexMatrix * M =  ((ComplexMPSSite * )this->Site)->GetM();
      unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
     
      for (int i = 0; i < this->NbrNonZeroElements; i++)
	{
          this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);      
		  
	      for (int RightC = 0;RightC < this->Site->GetBondDimensionRight() ; RightC++ )
		{
	  for (int RightA = 0;RightA < this->Site->GetBondDimensionRight() ; RightA++ )
	    {
		  L(RightA, MPOIndiceRight,RightC) +=  Conj(M[MPOIndiceUp].GetMatrixElement(0,RightC)) * this->ElementsValues[i] * this->LeftVector[MPOIndiceLeft] * M[MPOIndiceDown].GetMatrixElement(0,RightA);
	    }
	    }
	}
    }
  else
    this->ComputeLCore(L);
//      L.PrintTensor();
}


void ComplexMPOperatorOBC::ComputeR(Tensor3<Complex> & R)
{
  if (this->Site->GetSitePosition() == this->NbrSites - 1)
    {
      ComplexMatrix * M = ((ComplexMPSSite *)this->Site)->GetM();

      int BondDimensionLeft = this->Site->GetBondDimensionLeft();
      unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
      for (int i = 0; i < this->NbrNonZeroElements; i++)
	{
           this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);      
		  
	      for (int LeftC = 0;LeftC < this->Site->GetBondDimensionLeft() ; LeftC++ )
		{
	  for (int LeftA = 0; LeftA < this->Site->GetBondDimensionLeft() ; LeftA++ )
	    {
		  R(LeftA, MPOIndiceLeft,LeftC) +=  Conj(M[MPOIndiceUp].GetMatrixElement(LeftC,0)) * this->ElementsValues[i]*this->RightVector[MPOIndiceRight] * M[MPOIndiceDown].GetMatrixElement(LeftA,0);
		}
	    }
	}
    }
  else
    this->ComputeRCore(R);
// R.PrintTensor();
}



// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ComplexMPOperatorOBC::LowLevelMultiplyOneSite(ComplexVector& vSource, ComplexVector& vDestination, 
				       int firstComponent, int nbrComponent)
{
//cout <<"ComplexVector& ComplexMPOperatorOBC::LowLevelMultiplyOneSite(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)"<<endl;
  vDestination.ClearVector();
  if (this->Site->GetSitePosition() == 0)
  {
    int BondDimensionRight = this->Site->GetBondDimensionRight(); 
    Tensor3<Complex> & RightR = ((ComplexMPSSite * )this->Site)->GetNextR();
    Tensor3<Complex> * B = new Tensor3<Complex>[this->PhysicalDimension];

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
        B[i] = Tensor3<Complex>(this->MPOBondDimension,BondDimensionRight,1,true);
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
    Tensor3<Complex> & LeftL = ((ComplexMPSSite * )this->Site)->GetPreviousL();
    Tensor3<Complex> * B = new Tensor3<Complex>[this->PhysicalDimension];

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      B[i] = Tensor3<Complex>(this->MPOBondDimension,BondDimensionLeft,1,true);

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

  return this->LowLevelMultiplyCore(vSource,vDestination,firstComponent,nbrComponent);

}
 
// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ComplexMPOperatorOBC::LowLevelMultiplyTwoSites(ComplexVector& vSource, ComplexVector& vDestination,  int firstComponent, int nbrComponent)
{
//  cout<<"using ComplexVector& ComplexMPOperatorOBC::LowLevelMultiplyTwoSites(ComplexVector& vSource, ComplexVector& vDestination,  int firstComponent, int nbrComponent)"<<endl;
  vDestination.ClearVector();
  if ((this->SiteLeft->GetSitePosition() == 0)&&(this->SiteRight->GetSitePosition() ==  this->NbrSites - 1))
  {

  int SquarePhysicalDimension = this->PhysicalDimension * this->PhysicalDimension;

  Tensor3<Complex> * B =  new  Tensor3<Complex>  [SquarePhysicalDimension];
  for (int i = 0; i < SquarePhysicalDimension; i++)
    {
      B[i] = Tensor3<Complex>(this->MPOBondDimension,1,1,true);
    }

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
      vDestination[(long) PhysicalIndiceLeft + this->PhysicalDimension * MPOIndiceUp] += this->ElementsValues[i] * this->RightVector[MPOIndiceRight] *  B[PhysicalIndiceLeft +  this->PhysicalDimension *  MPOIndiceDown](MPOIndiceLeft,0,0);
    }
  } 

 delete [] B;
 return vDestination;
 }

 return this->LowLevelMultiplyTwoSitesCore(vSource,vDestination,firstComponent,nbrComponent);
}


void ComplexMPOperatorOBC::PrintTensorElements()
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
