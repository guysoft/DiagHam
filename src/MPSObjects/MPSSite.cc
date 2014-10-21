#include "MPSSite.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "AbstractMPOperatorOBC.h"
#include "GeneralTools/GarbageFlag.h"

using std::endl;
using std::cout;

MPSSite::MPSSite()
{
  this->PhysicalDimension = 0;
  this->SiteOnLeft = 0;
  this->SiteOnRight = 0;
  this->SitePosition = 0;
  this->M = 0;
  this->L = 0;
  this->R = 0;
}

MPSSite::MPSSite(unsigned int sitePosition, unsigned int physicalDimension, MPSSite * siteOnLeft, MPSSite * siteOnRight , unsigned int bondDimension, AbstractMPOperatorOBC * mPOperator)
{
  this->PhysicalDimension = physicalDimension;
  this->SiteOnLeft = siteOnLeft;
  this->SiteOnRight = siteOnRight;
  this->SitePosition = sitePosition;
  this->OperatorToBeMinimized = mPOperator;
  this->Flag.Initialize();
  this->M = new RealMatrix[this->PhysicalDimension];
  this->L = 0;
  this->R = 0;
}

MPSSite::~MPSSite()
{
  //for (int i = 0 ; i < this->PhysicalDimension ; i++)
    //delete this->M[i];

  delete [] this->M;
  delete this->L;
  delete this->R;
}



// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

MPSSite & MPSSite::operator = (const MPSSite & site) 
{
  this->PhysicalDimension = site.PhysicalDimension;
  this->SiteOnLeft = site.SiteOnLeft;
  this->SiteOnRight = site.SiteOnRight;
  this->SitePosition = site.SitePosition;
  this->OperatorToBeMinimized = site.OperatorToBeMinimized;
  this->BondDimensionLeft = site.BondDimensionLeft;
  this->BondDimensionRight = site.BondDimensionRight;
  delete this->M;
  this->M = new RealMatrix[this->PhysicalDimension];
  for(int i= 0; i <this->PhysicalDimension; i++)
    this->M[i] = site.M[i];

  this->L = site.L;
  this->R = site.R;
  return *this;
}


void MPSSite::InitializeLeft(RealMatrix * newA)
{
  delete this->M;
  this->M = newA;
  this->BondDimensionLeft = this->M[0].GetNbrRow();
  this->BondDimensionRight = this->M[0].GetNbrColumn();
  delete this->L;
  this->L = new Tensor3<double>(this->BondDimensionRight,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionRight,true);
  cout <<this->BondDimensionLeft<<" "<< this->BondDimensionRight<<" "<<this->PhysicalDimension<<endl;
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeL(*this->L);
  cout <<"I have to print"<<endl;
  this->L->PrintTensor();
  cout <<"I have finish printing"<<endl;
}

void MPSSite::SetBondDimension(int bondDimensionLeft, int bondDimensionRight)
{
  this->BondDimensionLeft = bondDimensionLeft;
  this->BondDimensionRight = bondDimensionRight;
  cout << "Setting site Dimension = " <<this->SitePosition << " " <<this->BondDimensionLeft <<" " <<this->BondDimensionRight<<endl;;
}

void MPSSite::UpdateFromVector(RealVector & psi)
{
  delete this->M;
  this->M = new RealMatrix [this->PhysicalDimension];
  
  for(int i = 0; i < this->PhysicalDimension; i++)
    {
      this->M[i] = RealMatrix(this->BondDimensionLeft,this->BondDimensionRight, true);
    }
  
  /*
  for (int i = 0; i < psi.Dimension; i++)
    {
      this->M[i%d](i/(this->PhysicalDimension *this->this->BondDimensionRight),i/(this->PhysicalDimension)%this->this->BondDimensionRight) = psi[i];
      
      }
  */
  
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      for (int j = 0; j < this->BondDimensionLeft; j++)
	{
	  for (int k = 0; k < this->BondDimensionRight; j++)
	    {
	      this->M[i](j,k) = psi[(long int)this->PhysicalDimension*(this->BondDimensionRight*j+k) +i]; 
	    }
	}
    }
}


void MPSSite::InitializeRight(RealMatrix * newB)
{
  delete this->M;
  this->M = newB;
  this->BondDimensionLeft = this->M[0].GetNbrRow();
  this->BondDimensionRight = this->M[0].GetNbrColumn();
  delete this->R;
  this->R = new Tensor3<double>(this->BondDimensionLeft,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionLeft,true);
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeR(*this->R);
  this->R->PrintTensor();
}

bool MPSSite::CheckLeftNormalization()
{
  RealMatrix Result(this->BondDimensionLeft,this->BondDimensionLeft,true);
  for (int i = 0 ; i< this->PhysicalDimension ; i++)
    {
      RealMatrix Tmp = this->M[i].DuplicateAndTranspose();
      Result += Tmp * (this->M[i]);
    }
  
  for(int i = 0 ; i < this->BondDimensionLeft; i++)
    {
      if (Result(i,i) != 1 )
	  {
	    return false;
	  }
      for(int j = 1 ; j < this->BondDimensionLeft ; j++)
	{
	  if (Result(i,j) != 0 )
	    {
	      return false;
	    }
	}
    }
  return true;
}


bool MPSSite::CheckRightNormalization()
{
  RealMatrix Result(this->BondDimensionRight,this->BondDimensionRight,true);
  for (int i = 0 ; i < this->PhysicalDimension ; i++)
    {
      RealMatrix Tmp = this->M[i].DuplicateAndTranspose ();
      Result += (this->M[i])* Tmp;
    }
  
  for(int i = 0 ; i < this->BondDimensionRight; i++)
    {
      if (Result(i,i) != 1 )
	{
	  return false;
	}
      for(int j= 1;j< this->BondDimensionRight; j++)
	{
	  if (Result(i,j) != 0 )
	    {
	      return false;
	    }
	}
      
    }
  return true;
}

//can be used only if all matrices on the right sites are right-normalized
void MPSSite::BringMInRightCanonicalForm()
{
 cout << "BringMInRightCanonicalForm() for site" << this->SitePosition<<endl;
  RealMatrix TmpMatrix (this->BondDimensionLeft,this->BondDimensionRight * this->PhysicalDimension, true);
  
  for(int i = 0 ; i < this->BondDimensionLeft ; i++)
    {
      for(int p = 0 ; p < this->PhysicalDimension ; p++)
	{
	  for(int j = 0; j < this->BondDimensionRight; j++)
	    { 
	      TmpMatrix(i, this->PhysicalDimension * j + p) = this->M[p](i,j);
	    }
	}
    }
  cout << TmpMatrix<<endl;
  RealMatrix U,V;
  RealDiagonalMatrix SingularValues;
  TmpMatrix.SingularValueDecomposition(U,SingularValues,V,false);
  cout <<"SingularValues  = "<< SingularValues<<endl;
  cout <<U <<endl;
  cout <<V <<endl;
  U = U * SingularValues;
  cout <<  SingularValues<< endl;
  cout <<"After Multiplication"<<endl;
  cout <<U<<endl;;
 // U.Multiply(SingularValues);
  cout <<"After Multiplication"<<endl;
  for(int i =0 ; i <  this->PhysicalDimension; i++)
    {
      for(int j = 0 ; j <  this->BondDimensionLeft; j++)
	{
	  for(int k = 0 ; k <  this->BondDimensionRight; k++)
	    {
	      this->M[i](j,k) = V(this->PhysicalDimension * k + i,j);
	    }
	}
      cout <<"check multiplication order in MPSSite::BringMInRightCanonicalForm()"<<endl;
      this->SiteOnLeft->M[i] = this->SiteOnLeft->M[i]*U;
       }
  delete this->R;
  this->R = new Tensor3<double> (this->BondDimensionLeft,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionLeft) ;
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeR(*this->R);
}



// can be used only if all matrices on the right sites are right-normalized
// check the result

void MPSSite::BringMInRightCanonicalFormCareful()
{
  this->BringMInRightCanonicalForm();
  
if(this->CheckRightNormalization())
  ;
  else 
    cout <<"Right Normalization issue in BringMinRightCanonicalFormCareful()"<<endl;
}


// can be used only if all matrices on the left sites are left-normalized

void MPSSite::BringMInLeftCanonicalForm()
{
  RealMatrix TmpMatrix (this->PhysicalDimension*this->BondDimensionLeft,this->BondDimensionRight, true);
  
  for(int i = 0; i < this->BondDimensionLeft; i++)
    {
      for(int p = 0; p <   this->PhysicalDimension ; p++)
	{
	  for(int j = 0; j < this->BondDimensionRight; j++)
	    { 
	      TmpMatrix(this->PhysicalDimension*i + p,j) = this->M[p](i,j);
	    }
	}
    }
  
  RealMatrix U,V;
  RealDiagonalMatrix SingularValues;
  
  TmpMatrix.SingularValueDecomposition(U,SingularValues,V);
  
  V = V.Transpose()*SingularValues;
  for(int i = 0 ; i <  this->PhysicalDimension; i++)
    {
      for(int j = 0 ; j <  this->BondDimensionLeft; j++)
	{
	  for(int k = 0 ; k <  this->BondDimensionRight; k++)
	    {
	      this->M[i](j,k) = U(this->PhysicalDimension*j + i,k);
	    }
	}

      cout <<"check multiplication order in MPSSite::BringMInLeftCanonicalForm()"<<endl;
      
      this->SiteOnRight->M[i] =  V * this->SiteOnRight->M[i];
    }

  delete this->L;
  this->L = new Tensor3<double> (this->BondDimensionRight,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionRight) ;
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeL(*this->L);
}


// can be used only if all matrices on the left sites are left-normalized
// check the result

void MPSSite::BringMInLeftCanonicalFormCareful()
{
  this->BringMInLeftCanonicalForm();
  if(this->CheckLeftNormalization())
    ;
  else 
    cout <<"Left Normalization issue in BringMinLeftCanonicalFormCareful()"<<endl;
}


RealVector & MPSSite::GetMatrixInVectorForm()
{
  RealVector MatrixInVectorForm ((long int) this->PhysicalDimension*this->BondDimensionLeft*this->BondDimensionRight,true);
  for (int i = 0 ; i < this->PhysicalDimension ; i++)
    {
      for (int j = 0 ; j <  this->BondDimensionLeft ; j++)
	{
	  for (int k = 0 ; k <  this->BondDimensionRight ; k++)
	    {
	      MatrixInVectorForm[(long int) this->PhysicalDimension*( this->BondDimensionRight*j +k )+i ] = this->M[i](j,k); 
	    }
	}
    }
  return MatrixInVectorForm;
}


void MPSSite::InitializeWithRandomMatrices()
{
//   cout <<" start initialization i = "<<this->SitePosition <<endl;;
//   cout <<this->BondDimensionLeft<< " " << this->BondDimensionRight<<endl;
   for (int i = 0; i < this->PhysicalDimension; i++)
   {
    this->M[i] = RealMatrix(this->BondDimensionLeft,this->BondDimensionRight, true);
    for (int j = 0; j < this->BondDimensionLeft ; j++)
     for (int k = 0; k < this->BondDimensionRight ; k++)
      this->M[i](j,k) = ((double) rand() / (RAND_MAX) - 0.5);
   }
}
