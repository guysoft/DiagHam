#include "AbstractMPOperatorOBC.h"
#include "Tensor/Tensor3.h"
#include <iostream>

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
  return this->HilbertSpace->GetHilbertSpaceDimension();
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
  cout <<"I am in AbstractMPOperatorOBC"<<endl;
  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  
  Tensor3<double> & LeftL = this->Site->GetPreviousL();
  
  Tensor3<double> * B =  new  Tensor3<double> [this->PhysicalDimension];
  RealMatrix * M = this->Site->GetM();
   for (int i = 0; i < this->PhysicalDimension; i++)
    {
      cout <<"M[i] = "<< M[i]<<endl;
      B[i] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft);
    }
  
  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      int MPOIndiceDown = this->GetIndiceDownFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceLeft = this->GetIndiceLeftFromTensorIndex(this->IndexValues[i]);
      for (int RightA = 0; RightA < BondDimensionRight; RightA++)
	{
	  for(int LeftC = 0; LeftC < BondDimensionLeft; LeftC++)
	    {
	      for(int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
		{
		  B[MPOIndiceDown](RightA,MPOIndiceLeft,LeftC) +=  LeftL(LeftA,MPOIndiceLeft,LeftC)* M[i](LeftA,RightA);
		}
	    }
	}
    }
  
  Tensor3<double> * A = new Tensor3<double>[this->PhysicalDimension];
  
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      A[i] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft);
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
      cout <<"M[i] = "<< M[i]<<endl;
      B[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
    }
  
  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      int MPOIndiceDown = this->GetIndiceDownFromTensorIndex(this->IndexValues[i]);
      int MPOIndiceRight = this->GetIndiceRightFromTensorIndex(this->IndexValues[i]);
      for (int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
	{
	  for(int RightC = 0; RightC < BondDimensionRight; RightC++)
	    {
	      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		{
		  //cout <<"M[i](LeftA,RightA)"<<M[i](LeftA,RightA)<<endl;
                  cout <<" RightR(RightA,MPOIndiceRight,RightC) "<<RightR(RightA,MPOIndiceRight,RightC)<<endl;
		  B[MPOIndiceDown](LeftA,MPOIndiceRight,RightC) +=  RightR(RightA,MPOIndiceRight,RightC)* M[i](LeftA,RightA);
		  //cout <<"B[MPOIndiceDown](LeftA,MPOIndiceRight,RightC) "<<B[MPOIndiceDown](LeftA,MPOIndiceRight,RightC)<<endl;
		}
	    }
	}
    }

  for (int i = 0; i < this->PhysicalDimension; i++)
    { 
      B[i].PrintTensor();
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




/*
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      B[i] = new Tensor3<T>(site->BondDimensionRight,site->BondDimensionLeft,this->MPOBondDimension);
      for (int c = 0; c < this->BondDimensionLeft; c++)
	{
	  for(int previousa = 0; previousa < this->BondDimensionLeft; previousa++)
	    {
	      for(int previousb = 0; previousb < this->MPOperatorOBC->GetMPODimension(); previousb++)
		{
		  for(int previousc = 0; previousc < this->BondDimensionLeft; previousc++)
		    B[i](c,previousa,previousb) +=  LeftL(previousa,previousb,previousc)* (this->M[i] (previousc,c));
		}
	    }
	}
    }
  
  Tensor3<double> * A = new Tensor3<double>[this->PhysicalDimension];
  
  for (int i = 0; i<this->PhysicalDimension; i++)
    {
      A[i] = new Tensor3<double>(this->BondDimensionRight,this->BondDimensionLeft,this->MPOperatorOBC->GetMPODimension(this->SitePosition));
      
      for (int a = 0; a < this->BondDimensionLeft; a++)
	{
	  for (int previousa = 0; previousa < this->BondDimensionLeft; previousa++)
	    {
	      for (int b = 0; b < this->MPOperatorOBC->GetMPODimension(this->SitePosition); b++)
		{
		  for (int previousb = 0; previousb < this->MPOperatorOBC->GetMPODimension(this->SitePosition-1); previousb++)
		    {
		      for(int sigma2 = 0; sigma2 < this->PhysicalDimension; sigma2++)
			{
			  A[i](a,previousa,b) +=  B[sigma2](c , previousa,previousb) * this->MPOperatorOBC->GetTensorElements(this->SitePosition,i,sigma2, previousb, b);
			}
		    }
		}
		}
	    }
    }
    
    for (int i = 0; i<this->PhysicalDimension; i++)
    {
    delete B[i];
    }
  delete B;
  
  for (int a = 0; a < this->BondDimensionLeft; a++)
    {
      for (int c = 0; c < this->BondDimensionLeft; c++)
	{
	  for (int b = 0; b < this->MPOperatorOBC->GetMPODimension(this->SitePosition); b++)
	    { 
	      this->L(a,b,c) = 0;
	      for (int previousa = 0; previousa < this->BondDimensionLeft; previousa++)
		{ 
		  for (int i = 0; i<this->PhysicalDimension; i++)
		    {
		      this->L(a,b,c)  += Conj(this->M[i](previousa,a))* A[i](a,previousa,b);
		    }
		}
	    }
	}
    }

  for (int i = 0; i<this->PhysicalDimension; i++)
     {
       delete A[i];
     }
  delete [] A;
  
  return &NewL;
}

*/
