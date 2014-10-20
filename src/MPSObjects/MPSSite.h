
#ifndef _MPSSITE_H
#define _MPSSITE_H

#include "AbstractMPOperatorOBC.h"
#include "Tensor/Tensor3.h"
#include "Matrix/RealMatrix.h"

class AbstractMPOperatorOBC;

class MPSSite 
{
 protected:
   unsigned int PhysicalDimension;
   unsigned int BondDimensionLeft;
   unsigned int BondDimensionRight;
   unsigned int SitePosition;
   
   RealMatrix * M;  

   Tensor3<double> * L;
   Tensor3<double> * R;
   
   MPSSite * SiteOnLeft;
   MPSSite * SiteOnRight;


   MPSSite * TableAllSite;

   AbstractMPOperatorOBC * OperatorToBeMinimized;
   
 public:

   MPSSite();
   MPSSite(unsigned int sitePosition, unsigned int physicalDimension, MPSSite * siteOnLeft, MPSSite * siteOnRight , unsigned int bondDimension, AbstractMPOperatorOBC * mPOperator);

   ~MPSSite();

   // assignement (without duplicating datas)
   //
   // M = matrix to copy
   // return value = reference on modified matrix
   MPSSite & operator = (const MPSSite & site);

   
   bool CheckLeftNormalization();
   
   bool CheckRightNormalization();
   //   bool ComputeL();
   //bool ComputeR();
   
   void InitializeWithRandomMatrices();
   
   void InitializeLeft(RealMatrix * newA);
   void InitializeRight(RealMatrix * newB);
   void UpdateFromVector(RealVector & psi);
   RealVector & GetMatrixInVectorForm();          
   void BringMInLeftCanonicalForm();
   void BringMInRightCanonicalForm();
   void BringMInLeftCanonicalFormCareful();
   void BringMInRightCanonicalFormCareful();
   void SetBondDimension(int bondDimensionLeft, int bondDimensionRight);

   inline unsigned int GetBondDimensionRight() const
     {return this-> BondDimensionRight;}
   
   inline unsigned int GetBondDimensionLeft() const
     {return this-> BondDimensionLeft;}
   
   inline Tensor3<double> & GetPreviousL ()const
     { return (*this->SiteOnLeft->L);}

   inline Tensor3<double> & GetNextR ()const
     { return (*this->SiteOnRight->R);}
   inline RealMatrix * GetM() const
     { return this->M;}
   inline unsigned int GetSitePosition() const
      {return this->SitePosition;}
};

#endif
