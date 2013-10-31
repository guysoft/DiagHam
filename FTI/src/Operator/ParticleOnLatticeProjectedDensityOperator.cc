////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                   Class author: Cecile Repellin                            //
//                                                                            //
//                                                                            //
//         class of particle on lattice projected density operator            //
//                                                                            //
//                        last modification : 17/10/2013                      //
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


#include "config.h"
#include "Operator/ParticleOnLatticeProjectedDensityOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

  
using std::cout;
using std::endl;

//constructor
//
// particleSource = Hilbert space associated to the original system
// particleDestination = Hilbert space associated to the destination system
// tightBindingModel = tight binding model
// kx = total momentum in x-direction of the density operator
// ky = total momentum in y-direction of the density operator
ParticleOnLatticeProjectedDensityOperator::ParticleOnLatticeProjectedDensityOperator(ParticleOnSphere* particleSource, ParticleOnSphere* particleDestination, Abstract2DTightBindingModel* tightBindingModel, int kx, int ky)
{
  this->ParticleSource = particleSource;
  this->ParticleDestination = particleDestination;
  this->TightBindingModel = tightBindingModel;
  this->Kx=kx;
  this->Ky=ky;
  
  this->ParticleSource->SetTargetSpace(this->ParticleDestination);
}



//destructor
ParticleOnLatticeProjectedDensityOperator::~ParticleOnLatticeProjectedDensityOperator()
{
}


// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractOperator* ParticleOnLatticeProjectedDensityOperator::Clone()
{
  return new ParticleOnLatticeProjectedDensityOperator(*this);
}


// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeProjectedDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->ParticleSource = (ParticleOnSphere*) hilbertSpace;
}


// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeProjectedDensityOperator::GetHilbertSpace ()
{
  return this->ParticleSource;
}


// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeProjectedDensityOperator::GetHilbertSpaceDimension ()
{
  return this->ParticleSource->GetHilbertSpaceDimension();
}


// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnLatticeProjectedDensityOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								    int firstComponent, int nbrComponent)
{
  int Index;
  int Last = firstComponent + nbrComponent;
  double Coefficient;
  Complex oneBodyBasisCoefficient;
  int NbrSiteX = this->TightBindingModel->GetNbrSiteX();
  int NbrSiteY = this->TightBindingModel->GetNbrSiteY();
  
  for (int kx = 0; kx < NbrSiteX; ++kx)
  {
    for (int ky = 0; ky < NbrSiteY; ++ky)
    {
      oneBodyBasisCoefficient = 0;
      int LinearizedK = this->TightBindingModel->GetLinearizedMomentumIndex(kx, ky);
      int LinearizedKPlusQ = this->TightBindingModel->GetLinearizedMomentumIndexSafe(kx + this->Kx, ky + this->Ky);
      ComplexMatrix& LocalBasisK = this->TightBindingModel->GetOneBodyMatrix(LinearizedK);
      ComplexMatrix& LocalBasisKPlusQ = this->TightBindingModel->GetOneBodyMatrix(LinearizedKPlusQ);
//       cout << kx << " " << ky << " " << this->Kx << " " << this->Ky << " " << LinearizedK << " " << LinearizedKPlusQ << endl;
      for (int alpha = 0; alpha < this->TightBindingModel->GetNbrBands(); ++alpha)
      {
	oneBodyBasisCoefficient += Conj(LocalBasisK[0][alpha]) * LocalBasisKPlusQ[0][alpha];
// 	oneBodyBasisCoefficient += Conj(LocalBasisK[alpha][0]) * LocalBasisKPlusQ[alpha][0];
      }
      cout <<  kx << " " << ky << " " << ((kx + this->Kx) % NbrSiteX) << " " << ((ky + this->Ky) % NbrSiteY) << " : " << ((oneBodyBasisCoefficient.Re)*(oneBodyBasisCoefficient.Re) + (oneBodyBasisCoefficient.Im)*(oneBodyBasisCoefficient.Im)) << endl;
      for (int i = firstComponent; i < Last; ++i)
      {
	Index = this->ParticleSource->AdA(i, LinearizedK, LinearizedKPlusQ, Coefficient);
	if (Index < this->ParticleDestination->GetHilbertSpaceDimension())
	{
// 	  cout << kx << " " << ky << " " << ((kx + this->Kx) % NbrSiteX) << " " << ((ky + this->Ky) % NbrSiteY) << " : ";
// 	  this->ParticleSource->PrintState(cout, i) << " -> " ;
// 	  this->ParticleDestination->PrintState(cout, Index) << endl;
	  vDestination[Index] += oneBodyBasisCoefficient*Coefficient*vSource[i]/(NbrSiteX * NbrSiteY);
	}
      }
    }
  }
  double NormDestination = vDestination.Norm();
  cout << "Norm = " << NormDestination << endl;
  vDestination *= 1/NormDestination;
  return vDestination;
}
  


