#include <cmath>
#include "MPSSite.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "AbstractMPOperatorOBC.h"
#include "GeneralTools/GarbageFlag.h"

using std::endl;
using std::abs;
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
//  cout <<this->BondDimensionLeft<<" "<< this->BondDimensionRight<<" "<<this->PhysicalDimension<<endl;
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeL(*this->L);
  this->L->PrintTensor();
}

void MPSSite::SetBondDimension(int bondDimensionLeft, int bondDimensionRight)
{
  this->BondDimensionLeft = bondDimensionLeft;
  this->BondDimensionRight = bondDimensionRight;
  cout << "Setting site Dimension = " <<this->SitePosition << " " <<this->BondDimensionLeft <<" " <<this->BondDimensionRight<<endl;;
}

void MPSSite::UpdateFromVector(RealVector * psi)
{
//  cout <<"Enter UpdateFromVector(RealVector * psi)"<<endl;
//  cout <<*psi<<endl;
  delete [] this->M;
  this->M = new RealMatrix [this->PhysicalDimension];
  
  for(int i = 0; i < this->PhysicalDimension; i++)
    {
      this->M[i] = RealMatrix(this->BondDimensionLeft,this->BondDimensionRight, true);
    }
//  cout << this->PhysicalDimension<<" "<<this->BondDimensionLeft<<" "<<
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      for (int j = 0; j < this->BondDimensionLeft; j++)
	{
	  for (int k = 0; k < this->BondDimensionRight; k++)
	    {
	      this->M[i](j,k) = (*psi)[(long int)this->PhysicalDimension*(this->BondDimensionRight*j+k) +i]; 
	    }
	}
//         cout <<" update M"<<this->M[i]<<endl;
//	 cout <<"this->M[ " <<i<< "]"<< this->M[i]<<endl;
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
  RealMatrix Result(this->BondDimensionRight,this->BondDimensionRight,true);
  for (int i = 0 ; i< this->PhysicalDimension ; i++)
    {
      RealMatrix Tmp = this->M[i].DuplicateAndTranspose();
      Result += Tmp * (this->M[i]);
    }
  for(int i = 0 ; i < this->BondDimensionRight; i++)
    {
      if ( abs(Result(i,i)-1) > 1e-13  )
	  {
	    return false;
	  }
      for(int j = i+1 ; j < this->BondDimensionRight ; j++)
	{
	  if (abs(Result(i,j)) > 1e-13)
	    {
	      return false;
	    }
	}
    }
  return true;
}


bool MPSSite::CheckRightNormalization()
{
  RealMatrix Result(this->BondDimensionLeft,this->BondDimensionLeft,true);
  for (int i = 0 ; i < this->PhysicalDimension ; i++)
    {
      RealMatrix Tmp = this->M[i].DuplicateAndTranspose();
      Result += (this->M[i])* Tmp;
    }
//  cout <<Result<<endl;
  for(int i = 0 ; i < this->BondDimensionLeft; i++)
    {
      if ( abs(Result(i,i)-1) > 1e-13  )
	{
	  return false;
	}
      for(int j = i+1;j< this->BondDimensionLeft; j++)
	{
	  if (abs(Result(i,j)) > 1e-13)
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
// cout << "BringMInRightCanonicalForm() for site" << this->SitePosition<<endl;
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
//  cout << TmpMatrix<<endl;
  RealMatrix U,V;
  RealDiagonalMatrix SingularValues;

  TmpMatrix.SingularValueDecomposition(U,SingularValues,V,false);
//  cout <<"SingularValues  = "<< SingularValues<<endl;
//  cout <<U <<endl;
// cout <<V <<endl;
  U = U * SingularValues;
//  cout <<U<<endl;;

  for(int i =0 ; i <  this->PhysicalDimension; i++)
    {
      this->M[i].ClearMatrix();
      for(int j = 0 ; j <  this->BondDimensionLeft; j++)
	{
	  for(int k = 0 ; k <  this->BondDimensionRight; k++)
	    {
	      this->M[i](j,k) = V(j,this->PhysicalDimension * k + i);
	    }
	}
//      cout <<"this->M[ " <<i<< "]"<<endl<< this->M[i]<<endl;
  //    cout <<"check multiplication order in MPSSite::BringMInRightCanonicalForm()"<<endl;
       this->SiteOnLeft->M[i] = this->SiteOnLeft->M[i]*U;
       }
  delete this->R;
  this->R = new Tensor3<double> (this->BondDimensionLeft,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionLeft,true) ;
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
//  cout <<"check singular value decomposition in void MPSSite::BringMInLeftCanonicalForm()"<<endl;
//  cout <<  TmpMatrix<<endl;
  TmpMatrix.SingularValueDecomposition(U,SingularValues,V,false);
//  cout << SingularValues<<endl;
//  cout <<U<<endl;
//  cout <<V<<endl;
  V.Transpose();
  V = SingularValues*V;
//  cout <<V<<endl;
  for(int i = 0 ; i <  this->PhysicalDimension; i++)
    {
      this->M[i].ClearMatrix();
      for(int j = 0 ; j <  this->BondDimensionLeft; j++)
	{
	  for(int k = 0 ; k <  this->BondDimensionRight; k++)
	    {
	      this->M[i](j,k) = U(this->PhysicalDimension*j + i,k);
	    }
	}

//      cout <<"check multiplication order in MPSSite::BringMInLeftCanonicalForm()"<<endl;
//      cout<<" this->SiteOnRight->M[i]"<<endl;
//      cout<< this->SiteOnRight->M[i]<<endl; 
      this->SiteOnRight->M[i] =  V * this->SiteOnRight->M[i];
//      cout<<" this->SiteOnRight->M[i]"<<endl;
//      cout<< this->SiteOnRight->M[i]<<endl;
    }

  delete this->L;
  this->L = new Tensor3<double> (this->BondDimensionRight,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionRight,true);
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


void MPSSite::GetMatrixInVectorForm(RealVector *& resultInvector)
{
  delete resultInvector;
  resultInvector = new RealVector((long int) this->PhysicalDimension*this->BondDimensionLeft*this->BondDimensionRight,true);
//  cout <<"print Matrix in GetMatrixInVectorForm"<<endl;
//  cout << this->M[0]<<" " <<this->M[1];
  for (int i = 0 ; i < this->PhysicalDimension ; i++)
    {
      for (int j = 0 ; j <  this->BondDimensionLeft ; j++)
	{
	  for (int k = 0 ; k <  this->BondDimensionRight ; k++)
	    {
	      (*resultInvector)[(long int) this->PhysicalDimension*( this->BondDimensionRight*j +k )+i ] = this->M[i](j,k); 
	    }
	}
    }
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
