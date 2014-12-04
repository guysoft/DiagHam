#ifndef _ABSTRACTMPSSITE_H
#define _ABSTRACTMPSSITE_H

#include "AbstractMPOperatorOBC.h"
#include "Tensor/Tensor3.h"
#include "GeneralTools/GarbageFlag.h"

class AbstractMPOperatorOBC;

class AbstractMPSSite 
{
 protected:

   // garbage flag to avoid data duplication
   GarbageFlag Flag;
   AbstractMPSSite * SiteOnLeft;
   AbstractMPSSite * SiteOnRight;

   unsigned int PhysicalDimension;
   unsigned int BondDimensionLeft;
   unsigned int BondDimensionRight;
   unsigned int MaxBondDimension;
   unsigned int SitePosition;
   
   AbstractMPOperatorOBC * OperatorToBeMinimized;
   
 public:

   AbstractMPSSite();
   AbstractMPSSite(unsigned int sitePosition, unsigned int physicalDimension, unsigned int bondDimension, AbstractMPOperatorOBC * mPOperator);

   ~AbstractMPSSite();

   // assignement (without duplicating datas)
   //
   // M = matrix to copy
   // return value = reference on modified matrix
   virtual AbstractMPSSite & operator = (const AbstractMPSSite & site);

   
   virtual bool CheckLeftNormalization();
   virtual bool CheckRightNormalization();
   
   virtual void InitializeWithRandomMatrices() = 0;
   
//   void InitializeLeft(RealMatrix * newA);
//   void InitializeRight(RealMatrix * newB);
//   void UpdateFromVector(RealVector * psi);
//   void GetMatrixInVectorForm(RealVector *& resultInvector );          

   virtual void BringMInLeftCanonicalForm();
   virtual void BringMInRightCanonicalForm();
   virtual void BringMInLeftCanonicalFormCareful();
   virtual void BringMInRightCanonicalFormCareful();
   void SetBondDimension(int bondDimensionLeft, int bondDimensionRight);
   inline void SetRightDimension(int bondDimensionRight);
   inline void SetLeftDimension(int bondDimensionLeft);

   //void ComputeDensityMatrixRight();
//   void ComputeDensityMatrixLeft();

   //void SymmetricUpdateOfTwoSites(AbstractMPSSite * leftSite , AbstractMPSSite * rightSite, RealVector * psi, RealDiagonalMatrix & singularValues );

   inline unsigned int GetBondDimensionRight() const
     {return this-> BondDimensionRight;}
   
   inline unsigned int GetBondDimensionLeft() const
     {return this-> BondDimensionLeft;}
   
   inline unsigned int GetSitePosition() const
      {return this->SitePosition;}
};




inline void AbstractMPSSite::SetRightDimension(int bondDimensionRight)
{
  this->BondDimensionRight = bondDimensionRight;
  this->SiteOnRight->BondDimensionLeft = bondDimensionRight;
}


inline void AbstractMPSSite::SetLeftDimension(int bondDimensionLeft)
{
  this->BondDimensionLeft = bondDimensionLeft;
  this->SiteOnLeft->BondDimensionRight = bondDimensionLeft;
}

inline void AbstractMPSSite::SetBondDimension(int bondDimensionLeft, int bondDimensionRight)
{
  this->BondDimensionLeft = bondDimensionLeft;
  this->BondDimensionRight = bondDimensionRight;
}


#endif
