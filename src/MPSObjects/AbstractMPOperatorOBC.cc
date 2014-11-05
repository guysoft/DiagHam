#include "AbstractMPOperatorOBC.h"
#include "Tensor/Tensor3.h"
#include "Matrix/RealSymmetricMatrix.h"
#include <iostream>
#include <sys/time.h>

using std::cout;
using std::endl;

AbstractMPOperatorOBC::AbstractMPOperatorOBC()
{
  this->NbrNonZeroElements = 0;
  this->ElementsValues = 0;
  this->IndexValues = 0;
  this->PhysicalDimension = 0;
  this->MPOBondDimension = 0;
}

AbstractMPOperatorOBC::~AbstractMPOperatorOBC()
{
}


// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void AbstractMPOperatorOBC::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->HilbertSpace = hilbertSpace;
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* AbstractMPOperatorOBC::GetHilbertSpace ()
{
  return this->HilbertSpace;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractMPOperatorOBC::GetHilbertSpaceDimension ()
{
  return  this->Site->GetBondDimensionRight()* this->Site->GetBondDimensionLeft()*this->PhysicalDimension;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractMPOperatorOBC::GetTwoSitesHilbertSpaceDimension ()
{
  return  this->SiteRight->GetBondDimensionRight()* this->SiteLeft->GetBondDimensionLeft()*this->PhysicalDimension *this->PhysicalDimension;
}
 
// shift Hamiltonian from a given energy
//
// shift = shift value

void AbstractMPOperatorOBC::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}

// set site to be acted on
//
// site = pointer to the siteto use 

void AbstractMPOperatorOBC::SetSite (MPSSite* site)
{
 this->Site = site;
}


// set site to be acted on
//
// site = pointer to the siteto use 

void AbstractMPOperatorOBC::SetSiteLeftAndRight (MPSSite* siteLeft,MPSSite* siteRight)
{
 this->SiteLeft = siteLeft;
 this->SiteRight = siteRight;
}


void AbstractMPOperatorOBC::ComputeL(Tensor3<double> & L)
{
  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  
  Tensor3<double> & LeftL = this->Site->GetPreviousL();
  
  Tensor3<double> * B =  new  Tensor3<double> [this->PhysicalDimension];
  RealMatrix * M = this->Site->GetM();

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
     B[i] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
     for(int LeftC = 0; LeftC < BondDimensionLeft; LeftC++)
     {
     for(int LeftB = 0; LeftB < this->MPOBondDimension; LeftB++)
     {
      for (int RightA = 0; RightA < BondDimensionRight; RightA++)
	{
              double & Tmp =  B[i](RightA,LeftB,LeftC);
	      for(int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
		{
		   Tmp +=  LeftL(LeftA,LeftB,LeftC)* M[i](LeftA,RightA);
		}
	    }
           }
	}
    }
  
  Tensor3<double> * A = new Tensor3<double>[this->PhysicalDimension];

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
  A[i] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
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
              double & Tmp = L(RightA,RightB,RightC);
              for (int i = 0; i < this->PhysicalDimension; i++)
               {
	      for (int LeftC = 0;  LeftC < BondDimensionLeft; LeftC++)
		{ 
		   Tmp += M[i](LeftC,RightC) * A[i](RightA,RightB,LeftC);
		   }
		}
	    }
	}
    } 

  delete [] A;
}


void AbstractMPOperatorOBC::ComputeR(Tensor3<double> & R)
{
  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  Tensor3<double> & RightR = this->Site->GetNextR();
  
  Tensor3<double> * B =  new  Tensor3<double>  [this->PhysicalDimension];
  RealMatrix * M = this->Site->GetM();

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
     B[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
  for (int RightB = 0; RightB < this->MPOBondDimension ; RightB++)
    {
      for (int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
	{
	  for(int RightC = 0; RightC < BondDimensionRight; RightC++)
	    {
              double & Tmp =  B[i](LeftA,RightB,RightC);
	      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		{
		   Tmp +=  RightR(RightA,RightB,RightC) * M[i](LeftA,RightA);
		}
	    }
	}
    }
   }

  Tensor3<double> * A = new Tensor3<double> [this->PhysicalDimension];
  
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      A[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
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
              double & Tmp =  R(LeftA,LeftB,LeftC);
 for (int i = 0; i < this->PhysicalDimension; i++)
{
	      for (int RightC = 0;  RightC < BondDimensionRight; RightC++)
		{ 
	      Tmp += M[i](LeftC,RightC) * A[i](LeftA,LeftB,RightC);
		    }
		}
	    }
	}
    }

  delete [] A;
}




// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored
RealVector& AbstractMPOperatorOBC::LowLevelMultiply(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored
RealVector& AbstractMPOperatorOBC::LowLevelMultiplyTwoSites(RealVector& vSource, RealVector& vDestination)
{
  return this->LowLevelMultiplyTwoSites(vSource, vDestination, 0, this->GetTwoSitesHilbertSpaceDimension());
}



// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractMPOperatorOBC::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent)
{
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&TotalStartingTime, 0);
  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  Tensor3<double> & RightR = this->Site->GetNextR();
  Tensor3<double> & LeftL = this->Site->GetPreviousL();
 
  Tensor3<double> * B =  new  Tensor3<double>  [this->PhysicalDimension];

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      B[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight);
      for(int RightC = 0; RightC < BondDimensionRight; RightC++)
	    {
          for (int RightB = 0; RightB < this->MPOBondDimension ; RightB++)
            {

      for (int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
	{
              double & Tmp =  B[i](LeftA,RightB,RightC);
	      Tmp = 0;	
	      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		{
		  Tmp +=  RightR(RightA,RightB,RightC) * vSource[(long int)BondDimensionRight*(BondDimensionLeft*i+ LeftA) + RightA];
		}
	    }
	}
    }
   }

 gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout <<"First Part " << Dt << "s" << endl;
 gettimeofday (&TotalStartingTime, 0);

  Tensor3<double> * A =  new  Tensor3<double>  [this->PhysicalDimension];
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      A[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
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

  gettimeofday (&TotalEndingTime, 0);
  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout <<"Second Part " << Dt << "s" << endl;
 gettimeofday (&TotalStartingTime, 0);

// new = (long int)this->BondDimensionRight*(this->BondDimensionLeft*i+ LeftA) + RightA

// old = (long int) this->PhysicalDimension*( this->BondDimensionRight*j +k )+i

  for(int Index =  firstComponent; Index < LastComponent ;Index++)
  {
     double & Tmp = vDestination[Index];
  for (int RightB = 0;  RightB < this->MPOBondDimension;  RightB++)
	{
       for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	{
          Tmp +=  LeftL(LeftA,RightB, (Index/BondDimensionRight)%BondDimensionLeft) * A[Index/(BondDimensionRight*BondDimensionLeft)](LeftA, RightB,Index % BondDimensionRight);
        }
    }
  }

gettimeofday (&TotalEndingTime, 0);
Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
 cout <<"Third Part " << Dt << "s" << endl;
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

RealVector& AbstractMPOperatorOBC::LowLevelMultiplyTwoSites(RealVector& vSource, RealVector& vDestination,  int firstComponent, int nbrComponent)
{
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&TotalStartingTime, 0)
;
  int BondDimensionRight = this->SiteRight->GetBondDimensionRight();
  int BondDimensionLeft = this->SiteLeft->GetBondDimensionLeft();

  Tensor3<double> & RightR = this->SiteRight->GetNextR();
  Tensor3<double> & LeftL = this->SiteLeft->GetPreviousL();

  int SquarePhysicalDimension = this->PhysicalDimension * this->PhysicalDimension;

  Tensor3<double> * B =  new  Tensor3<double>  [SquarePhysicalDimension];
  for (int i = 0; i < SquarePhysicalDimension; i++)
    {
      B[i] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
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

  Tensor3<double> * A =  new  Tensor3<double>  [SquarePhysicalDimension];
  for (int i = 0; i < SquarePhysicalDimension ; i++)
    {
      A[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
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

  Tensor3<double> * C =  new  Tensor3<double>  [SquarePhysicalDimension];
  for (int i = 0; i < SquarePhysicalDimension ; i++)
    {
      C[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
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

  gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout <<"First Part " << Dt << "s" << endl;
  gettimeofday (&TotalStartingTime, 0);
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

gettimeofday (&TotalEndingTime, 0);
Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
 cout <<"Second Part " << Dt << "s" << endl;
 delete [] C;
 return vDestination;
}

// store Hamiltonian into an hermitian matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix

RealSymmetricMatrix& AbstractMPOperatorOBC::GetTwoSitesHamiltonian (RealSymmetricMatrix & M)
{
  RealVector TmpV1 (this->GetTwoSitesHilbertSpaceDimension(), true);
  RealVector TmpV2 (this->GetTwoSitesHilbertSpaceDimension(), true);
  for (int i = 0; i < this->GetTwoSitesHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->LowLevelMultiplyTwoSites(TmpV1, TmpV2);
      if (this->LeftHamiltonianVectorMultiplicationFlag == false)
	{
	  for (int j = i; j < this->GetTwoSitesHilbertSpaceDimension(); j++)
	    {
	      M.SetMatrixElement(i, j, TmpV2[j]);
	    }
	}
      else
	{
	  for (int j = i; j < this->GetHilbertSpaceDimension(); j++)
	    {
	      M.SetMatrixElement(j, i, TmpV2[j]);
	    }
	}
      TmpV1[i] = 0.0;
    }
  return M;  
}
