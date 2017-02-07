////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                      Class author Cecile Repellin                          //
//                                                                            //
//                                                                            //
//            class of bond energy operator for spin-pseudospin model         //
//                                                                            //
//                        last modification : 16/01/2017                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "Operator/BondEnergySpinPseudospinOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/ComplexVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"


using std::cout;
using std::endl;


// constructor from default datas
//
// chain = pointer to the Hilbert space
// nbrSpin = number of spins

BondEnergySpinPseudospinOperator::BondEnergySpinPseudospinOperator(Spin1_2ChainWithPseudospinAnd2DTranslation* chain, int xMomentum, int nbrSpinX, int yMomentum, int nbrSpinY, int siteIndex1, int siteIndex2, int bondIndex)
{
  this->Chain = chain;
  this->NbrSpinX = nbrSpinX;
  this->NbrSpinY = nbrSpinY;
  this->NbrSpin = this->NbrSpinX * this->NbrSpinY;
  this->XMomentum = xMomentum;
  this->YMomentum = yMomentum;
  
  this->SiteIndex1 = siteIndex1;
  this->SiteIndex2 = siteIndex2;
  this->BondIndex = bondIndex;
  
  // projection of kagome onto the s = 1/2 states on each triangle
  this->PseudospinCouplingElements = new double[3];    
  this->PseudospinCouplingElements[0] = 0.0;
  this->PseudospinCouplingElements[1] = 1.0/sqrt(3.0);
  this->PseudospinCouplingElements[2] = -1.0/sqrt(3.0);
  
  this-> PseudospinDiagCouplingElements = new double*[3];
  for (int i = 0; i < 3; ++i)
    this->PseudospinDiagCouplingElements[i] = new double[2];
  this->PseudospinDiagCouplingElements[0][0] = 1.0;
  this->PseudospinDiagCouplingElements[0][1] = -1.0/3.0;
  
  this->PseudospinDiagCouplingElements[1][0] = 0.0;
  this->PseudospinDiagCouplingElements[1][1] = 2.0/3.0;
  
  this->PseudospinDiagCouplingElements[2][0] = 0.0;
  this->PseudospinDiagCouplingElements[2][1] = 2.0/3.0;
  
  this->EvaluateExponentialFactors();
  
}

// destructor
//

BondEnergySpinPseudospinOperator::~BondEnergySpinPseudospinOperator()
{
}

// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* BondEnergySpinPseudospinOperator::Clone ()
{
  return new BondEnergySpinPseudospinOperator(this->Chain, this->XMomentum, this->NbrSpinX, this->YMomentum, this->NbrSpinY, this->SiteIndex1, this->SiteIndex2, this->BondIndex);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void BondEnergySpinPseudospinOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Chain = (Spin1_2ChainWithPseudospinAnd2DTranslation*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* BondEnergySpinPseudospinOperator::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int BondEnergySpinPseudospinOperator::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex BondEnergySpinPseudospinOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, 
						     long firstComponent, long nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  Complex Element = 0.0;
  int pos;
  int pos2;
  double coef2;
  double coef;
  double TmpCoefficient;
  
  int nbrTranslationsX;
  int nbrTranslationsY;
  
  int nbrTranslationsX1;
  int nbrTranslationsY1;
  int Tmp1;
  int Tmp2;
  if (this->BondIndex == 0)      //AC
  {
    Tmp1 = 0;
    Tmp2 = 2;
  }
  if (this->BondIndex == 1)      //AB
  {
    Tmp1 = 0;
    Tmp2 = 1;
  }
  if (this->BondIndex == 2)      //BC
  {
    Tmp1 = 1;
    Tmp2 = 2;
  }
  
  
  for (int i = (int) firstComponent; i < dim; ++i)
    {	 
	Complex TmpValue = V2[i] * 0.5;
	
	TmpCoefficient = this->Chain->SziSzj(this->SiteIndex1, this->SiteIndex2, i);
	TmpCoefficient *= this->Chain->JDiagonali(this->SiteIndex1, i, this->PseudospinDiagCouplingElements[Tmp2]);
	TmpCoefficient *= this->Chain->JDiagonali(this->SiteIndex2, i, this->PseudospinDiagCouplingElements[Tmp1]);
	Element += Conj(V1[i]) * TmpCoefficient * V2[i];
	
	pos = this->Chain->SmiSpjJiJj(this->SiteIndex2, this->SiteIndex1, i, this->PseudospinDiagCouplingElements[Tmp1], this->PseudospinDiagCouplingElements[Tmp2], coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  Element += Conj(V1[pos]) * TmpValue * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	
	pos = this->Chain->SmiSpjJoffiJj(this->SiteIndex2, this->SiteIndex1, i, this->PseudospinDiagCouplingElements[Tmp2], coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  Element += Conj(V1[pos]) * TmpValue * coef * this->PseudospinCouplingElements[Tmp1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	
	pos = this->Chain->SmiSpjJiJoffj(this->SiteIndex2, this->SiteIndex1, i, this->PseudospinDiagCouplingElements[Tmp1], coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  Element += Conj(V1[pos]) *  TmpValue * coef * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	
	pos = this->Chain->SmiSpjJoffiJoffj(this->SiteIndex2, this->SiteIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  Element += Conj(V1[pos]) *  TmpValue * coef * this->PseudospinCouplingElements[Tmp2] * this->PseudospinCouplingElements[Tmp1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
	
	pos = this->Chain->SmiSpjJiJj(this->SiteIndex1, this->SiteIndex2, i, this->PseudospinDiagCouplingElements[Tmp2], this->PseudospinDiagCouplingElements[Tmp1], coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  Element += Conj(V1[pos]) *  TmpValue * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	
	pos = this->Chain->SmiSpjJoffiJj(this->SiteIndex1, this->SiteIndex2, i, this->PseudospinDiagCouplingElements[Tmp1], coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  Element += Conj(V1[pos]) *  TmpValue * coef * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	
	pos = this->Chain->SmiSpjJiJoffj(this->SiteIndex1, this->SiteIndex2, i, this->PseudospinDiagCouplingElements[Tmp2], coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  Element += Conj(V1[pos]) *   TmpValue * coef * this->PseudospinCouplingElements[Tmp1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	
	pos = this->Chain->SmiSpjJoffiJoffj(this->SiteIndex1, this->SiteIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  Element += Conj(V1[pos]) *   TmpValue * coef * this->PseudospinCouplingElements[Tmp1] * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		

	coef = this->Chain->SziSzj(this->SiteIndex1, this->SiteIndex2, i);
	TmpCoefficient = coef *  this->Chain->JDiagonali(this->SiteIndex2, i, this->PseudospinDiagCouplingElements[Tmp1]) * this->PseudospinCouplingElements[Tmp2];
	if (TmpCoefficient != 0.0)
	{
	  pos = this->Chain->JOffDiagonali(this->SiteIndex1, i, coef2, nbrTranslationsX, nbrTranslationsY);
	  if (pos!= dim)
	    Element += Conj(V1[pos]) *   2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	}
	
	TmpCoefficient = coef *  this->Chain->JDiagonali(this->SiteIndex1, i, this->PseudospinDiagCouplingElements[Tmp2]) * this->PseudospinCouplingElements[Tmp1];
	if (TmpCoefficient != 0.0)
	{
	  pos = this->Chain->JOffDiagonali(this->SiteIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
	  if (pos!= dim)
	    Element += Conj(V1[pos]) *   2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	}  
	  
	pos2 = this->Chain->JoffiJoffj(this->SiteIndex1, this->SiteIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
	if (pos2!= dim)
	  Element += Conj(V1[pos2]) *   2.0 * TmpValue * coef * coef2 * this->PseudospinCouplingElements[Tmp1] * this->PseudospinCouplingElements[Tmp2] *  this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
    }
  return Element;
}

// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& BondEnergySpinPseudospinOperator::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
						    int firstComponent, int nbrComponent)
{
  int dim = (int) (firstComponent + nbrComponent);
  Complex TmpCoefficient;
  int pos;
  int pos2;
  double coef2;
  double coef;
  
  int nbrTranslationsX;
  int nbrTranslationsY;
  
  int nbrTranslationsX1;
  int nbrTranslationsY1;
  int Tmp1;
  int Tmp2;
  if (this->BondIndex == 0)      //AC
  {
    Tmp1 = 0;
    Tmp2 = 2;
  }
  if (this->BondIndex == 1)      //AB
  {
    Tmp1 = 0;
    Tmp2 = 1;
  }
  if (this->BondIndex == 2)      //BC
  {
    Tmp1 = 1;
    Tmp2 = 2;
  }
  for (int i = (int) firstComponent; i < dim; ++i)
    {
	Complex TmpValue = vSource[i] * 0.5;
// 	TmpIndex1 = this->GetLinearizedIndex(j - this->Offset, k - 1);
// 	TmpIndex2 = this->GetLinearizedIndex(j, k);
      
	TmpCoefficient = this->Chain->SziSzj(this->SiteIndex1, this->SiteIndex2, i);
	TmpCoefficient *= this->Chain->JDiagonali(this->SiteIndex1, i, this->PseudospinDiagCouplingElements[Tmp2]);
	TmpCoefficient *= this->Chain->JDiagonali(this->SiteIndex2, i, this->PseudospinDiagCouplingElements[Tmp1]);
	vDestination[i] += TmpCoefficient * vSource[i];
	
	pos = this->Chain->SmiSpjJiJj(this->SiteIndex2, this->SiteIndex1, i, this->PseudospinDiagCouplingElements[Tmp1], this->PseudospinDiagCouplingElements[Tmp2], coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  vDestination[pos] += TmpValue * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	
	pos = this->Chain->SmiSpjJoffiJj(this->SiteIndex2, this->SiteIndex1, i, this->PseudospinDiagCouplingElements[Tmp2], coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[Tmp1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	
	pos = this->Chain->SmiSpjJiJoffj(this->SiteIndex2, this->SiteIndex1, i, this->PseudospinDiagCouplingElements[Tmp1], coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	
	pos = this->Chain->SmiSpjJoffiJoffj(this->SiteIndex2, this->SiteIndex1, i, coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[Tmp2] * this->PseudospinCouplingElements[Tmp1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		
	
	pos = this->Chain->SmiSpjJiJj(this->SiteIndex1, this->SiteIndex2, i, this->PseudospinDiagCouplingElements[Tmp2], this->PseudospinDiagCouplingElements[Tmp1], coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  vDestination[pos] += TmpValue * coef * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	
	pos = this->Chain->SmiSpjJoffiJj(this->SiteIndex1, this->SiteIndex2, i, this->PseudospinDiagCouplingElements[Tmp1], coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	
	pos = this->Chain->SmiSpjJiJoffj(this->SiteIndex1, this->SiteIndex2, i, this->PseudospinDiagCouplingElements[Tmp2], coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[Tmp1] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	
	pos = this->Chain->SmiSpjJoffiJoffj(this->SiteIndex1, this->SiteIndex2, i, coef, nbrTranslationsX, nbrTranslationsY);
	if (pos != dim)
	  vDestination[pos] += TmpValue * coef * this->PseudospinCouplingElements[Tmp1] * this->PseudospinCouplingElements[Tmp2] * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
		

	coef = this->Chain->SziSzj(this->SiteIndex1, this->SiteIndex2, i);
	TmpCoefficient = coef *  this->Chain->JDiagonali(this->SiteIndex2, i, this->PseudospinDiagCouplingElements[Tmp1]) * this->PseudospinCouplingElements[Tmp2];
	if (TmpCoefficient != 0.0)
	{
	  pos = this->Chain->JOffDiagonali(this->SiteIndex1, i, coef2, nbrTranslationsX, nbrTranslationsY);
	  if (pos!= dim)
	    vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	}
	
	TmpCoefficient = coef *  this->Chain->JDiagonali(this->SiteIndex1, i, this->PseudospinDiagCouplingElements[Tmp2]) * this->PseudospinCouplingElements[Tmp1];
	if (TmpCoefficient != 0.0)
	{
	  pos = this->Chain->JOffDiagonali(this->SiteIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
	  if (pos!= dim)
	    vDestination[pos] += 2.0 * TmpValue * TmpCoefficient * coef2 * this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
	}  
	  
	pos2 = this->Chain->JoffiJoffj(this->SiteIndex1, this->SiteIndex2, i, coef2, nbrTranslationsX, nbrTranslationsY);
	if (pos2!= dim)
	  vDestination[pos2] += 2.0 * TmpValue * coef * coef2 * this->PseudospinCouplingElements[Tmp1] * this->PseudospinCouplingElements[Tmp2] *  this->ExponentialFactors[nbrTranslationsX][nbrTranslationsY];
    }
  return vDestination;
}


// evaluate all exponential factors
//   

void BondEnergySpinPseudospinOperator::EvaluateExponentialFactors()
{
  this->ExponentialFactors = new Complex*[this->NbrSpinX];
  for (int i = 0; i < this->NbrSpinX; ++i)
    { 
      this->ExponentialFactors[i] = new Complex[this->NbrSpinY];
      for (int j = 0; j < this->NbrSpinY; ++j)
	{ 
	  this->ExponentialFactors[i][j] = Phase(2.0 * M_PI * ((this->XMomentum * ((double) i) / ((double) this->NbrSpinX))
							       + (this->YMomentum * ((double) j) / ((double) this->NbrSpinY))));
	}
    }
}

