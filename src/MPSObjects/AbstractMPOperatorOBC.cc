#include "AbstractMPOperatorOBC.h"
#include "Tensor/Tensor3.h"
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
 
// shift Hamiltonian from a given energy
//
// shift = shift value

void AbstractMPOperatorOBC::ShiftHamiltonian (double shift)
{
  cout << "it is working"<<endl;
  this->HamiltonianShift = shift;
}

// set site to be acted on
//
// site = pointer to the siteto use 

void AbstractMPOperatorOBC::SetSite (MPSSite* site)
{
 this->Site = site;
}

void AbstractMPOperatorOBC::ComputeL(Tensor3<double> & L)
{
//  cout <<"void AbstractMPOperatorOBC::ComputeL(Tensor3<double> & L)"<<endl;
  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  
  Tensor3<double> & LeftL = this->Site->GetPreviousL();
  
  Tensor3<double> * B =  new  Tensor3<double> [this->PhysicalDimension];
  RealMatrix * M = this->Site->GetM();
   for (int i = 0; i < this->PhysicalDimension; i++)
    {
//      cout <<"M[i] = "<< M[i]<<endl;
      B[i] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
    }
  
  for (int i = 0; i < this->PhysicalDimension; i++)
    {

      for (int RightA = 0; RightA < BondDimensionRight; RightA++)
	{
	  for(int LeftB = 0; LeftB < this->MPOBondDimension; LeftB++)
    {
	  for(int LeftC = 0; LeftC < BondDimensionLeft; LeftC++)
	    {
	      for(int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
		{
		  B[i](RightA,LeftB,LeftC) +=  LeftL(LeftA,LeftB,LeftC)* M[i](LeftA,RightA);
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
  
  
  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      int MPOIndiceDown = this->GetIndiceDownFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceLeft = this->GetIndiceLeftFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceUp = this->GetIndiceUpFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceRight = this->GetIndiceRightFromTensorIndex(this->IndexValues[i]);
      
      for (int RightA = 0; RightA < BondDimensionRight; RightA++)
	{
	  for (int LeftC = 0;  LeftC < BondDimensionLeft;  LeftC++)
	    {
	      A[MPOIndiceUp](RightA, MPOIndiceRight,LeftC) +=  B[MPOIndiceDown](RightA,MPOIndiceLeft,LeftC) * this->ElementsValues[i];
	    }
	}
    }

  delete [] B;
  
  for (int RightA = 0; RightA < BondDimensionRight; RightA++)
    {
      for (int RightC = 0; RightC <  BondDimensionRight;  RightC++)
	{
	  for (int  RightB = 0;  RightB < this->MPOBondDimension;  RightB++)
	    { 
	      L(RightA,RightB,RightC) = 0;
	      for (int LeftC = 0;  LeftC < BondDimensionLeft; LeftC++)
		{ 
		  for (int i = 0; i < this->PhysicalDimension; i++)
		    {
		      L(RightA,RightB,RightC) += M[i](LeftC,RightC) * A[i](RightA,RightB,LeftC);
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
    }
  
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
  for (int RightB = 0; RightB < this->MPOBondDimension ; RightB++)
    {
      for (int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
	{
	  for(int RightC = 0; RightC < BondDimensionRight; RightC++)
	    {
	      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		{
		   B[i](LeftA,RightB,RightC) +=  RightR(RightA,RightB,RightC) * M[i](LeftA,RightA);
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
  

  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      int MPOIndiceDown = this->GetIndiceDownFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceLeft = this->GetIndiceLeftFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceUp = this->GetIndiceUpFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceRight = this->GetIndiceRightFromTensorIndex(this->IndexValues[i]);
      
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
	      R(LeftA,LeftB,LeftC) = 0;
	      for (int RightC = 0;  RightC < BondDimensionRight; RightC++)
		{ 
		  for (int i = 0; i < this->PhysicalDimension; i++)
		    {
		      R(LeftA,LeftB,LeftC) += M[i](LeftC,RightC) * A[i](LeftA,LeftB,RightC);
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
      B[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
    }



  for (int i = 0; i < this->PhysicalDimension; i++)
    {
  for (int RightB = 0; RightB < this->MPOBondDimension ; RightB++)
    {
      for (int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
	{
	  for(int RightC = 0; RightC < BondDimensionRight; RightC++)
	    {
	      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		{
		  B[i](LeftA,RightB,RightC) +=  RightR(RightA,RightB,RightC) * vSource[(long int) this->PhysicalDimension*(BondDimensionRight * LeftA +RightA)+i];
		}
	    }
	}
    }
   }

  Tensor3<double> * A =  new  Tensor3<double>  [this->PhysicalDimension];
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      A[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
    }
 
 
 for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      int MPOIndiceDown = this->GetIndiceDownFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceLeft = this->GetIndiceLeftFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceUp = this->GetIndiceUpFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceRight = this->GetIndiceRightFromTensorIndex(this->IndexValues[i]);


      for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	{
	  for (int RightC = 0;  RightC < BondDimensionRight;  RightC++)
	    {
	      A[MPOIndiceUp](LeftA, MPOIndiceLeft,RightC) +=  B[MPOIndiceDown](LeftA,MPOIndiceRight,RightC) * this->ElementsValues[i];
	    }
	}
   } 

  delete [] B;
  int LastComponent = firstComponent + nbrComponent;

  gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout <<"First Part " << Dt << "s" << endl;
 gettimeofday (&TotalStartingTime, 0);
  for(int Index =  firstComponent; Index < LastComponent ;Index++)
  {
       for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	{
      for (int RightB = 0;  RightB < this->MPOBondDimension;  RightB++)
	{
            vDestination[Index] +=  LeftL(LeftA,RightB, Index/(BondDimensionRight* this->PhysicalDimension) ) * A[Index%this->PhysicalDimension](LeftA, RightB, (Index/this->PhysicalDimension) % BondDimensionRight);
        }
    }
 }
gettimeofday (&TotalEndingTime, 0);
Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
 cout <<"Second Part " << Dt << "s" << endl;
 delete [] A;
 return vDestination;
}
