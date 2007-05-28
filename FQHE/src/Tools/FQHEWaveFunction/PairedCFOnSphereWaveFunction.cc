////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                  //
//                                                                            //
//                                                                            //
//           class implementing a paired CF wave function on the sphere          //
//                                                                            //
//                        last modification : 18/05/2007                      //
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
#include "PairedCFOnSphereWaveFunction.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "Vector/RealVector.h"
#include "MathTools/FactorialCoefficient.h"

#include <iostream>

using std::cout;
using std::endl;


// default constructor
//

PairedCFOnSphereWaveFunction::PairedCFOnSphereWaveFunction()
{
  this->NbrParticles = 0;  
}

  // constructor
  //
  // nbrParticles = number of particles
  // nbrLandauLevel = number of Landau levels filled with composite fermions
  // nbrEffectiveFlux = number of flux quanta of the magnetic monopole field experienced by CF's
  // MooreReadCoefficient = prefactor of singular 1/z term in pair-wave function
  // CFCoefficients = prefactors of CF orbitals in shells 0, 1, 2, ... , nbrLandauLevel-1
  // correctPrefactors = flag that enables the correction of prefactors to adopt the conventions of previous code
  // jastrowPower = power to which the Jastrow factor has to be raised


PairedCFOnSphereWaveFunction::PairedCFOnSphereWaveFunction(int nbrParticles, int nbrLandauLevels, int nbrEffectiveFlux,
							   double MooreReadCoefficient, double * givenCFCoefficients,
							   bool correctPrefactors, int jastrowPower)
{
  this->NbrParticles = nbrParticles;
  this->NbrLandauLevels = nbrLandauLevels;
  this->AbsEffectiveFlux = abs(nbrEffectiveFlux);
  this->Orbitals = new JainCFOnSphereOrbitals(nbrParticles, nbrLandauLevels, nbrEffectiveFlux,jastrowPower);
  this->MooreReadCoefficient=MooreReadCoefficient;
  this->CFCoefficients= new double [NbrLandauLevels];
  for (int i=0; i<NbrLandauLevels; ++i) this->CFCoefficients[i] = givenCFCoefficients[i];
  this->ElementNorm=1.0;
  this->Slater = new ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->Flag.Initialize();
  this->Ji = new Complex[this->NbrParticles];

  this->gAlpha = new Complex*[this->NbrLandauLevels];
  for (int i=0; i< this->NbrLandauLevels; ++i)
    gAlpha[i]= new Complex[NbrParticles*NbrParticles];
  
  if (correctPrefactors)
    {
      int p=Orbitals->GetJastrowPower();
      FactorialCoefficient Coef;
      for (int n=0; n<NbrLandauLevels; ++n)
	{
	  Coef.SetToOne();
	  Coef.PartialFactorialDivide(nbrEffectiveFlux+p*(nbrParticles-1)+2,nbrEffectiveFlux+2*p*(nbrParticles-1)+1);
	  Coef.PartialFactorialMultiply(p*(nbrParticles-1)+n+2,2*p*(nbrParticles-1)+n+1);
	  this->CFCoefficients[n]*=Coef.GetNumericalValue()*Coef.GetNumericalValue();
	  //cout << "Correction["<<n<<"]="<<Coef.GetNumericalValue()<<"from: " << nbrEffectiveFlux+p*(nbrParticles-1)+2<<","<<nbrEffectiveFlux+2*p*(nbrParticles-1)+1<<","<< p*(nbrParticles-1)+n+2<<","<<2*p*(nbrParticles-1)+n+1<<endl;
	}
    }  
}

// copy constructor
//
// function = reference on the wave function to copy

PairedCFOnSphereWaveFunction::PairedCFOnSphereWaveFunction(const PairedCFOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->NbrLandauLevels = function.NbrLandauLevels;
  this->AbsEffectiveFlux = function.AbsEffectiveFlux;
  this->Flag = function.Flag;
  this->Orbitals = function.Orbitals;
  this->MooreReadCoefficient=function.MooreReadCoefficient;
  this->CFCoefficients=function.CFCoefficients;
  this->ElementNorm=function.ElementNorm;
  
  this->Slater = new ComplexSkewSymmetricMatrix(this->NbrParticles);
  this->Ji = new Complex[this->NbrParticles];
  this->gAlpha = new Complex*[this->NbrLandauLevels];
  for (int i=0; i< this->NbrLandauLevels; ++i)
    gAlpha[i]= new Complex[NbrParticles*NbrParticles];
  
}

// destructor
//

PairedCFOnSphereWaveFunction::~PairedCFOnSphereWaveFunction()
{
  if ( (this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete Orbitals;
      delete [] CFCoefficients;
    }
  for (int i=0; i< this->NbrLandauLevels; ++i) delete [] gAlpha[i];
  delete [] gAlpha;
  delete [] Ji;
  delete Slater;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PairedCFOnSphereWaveFunction::Clone ()
{
  return new PairedCFOnSphereWaveFunction(*this);
}


double fsgn(int x) // for calculating (-1)^x
{
  if (x&1) return -1.0;
  else return 1.0;
}


// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PairedCFOnSphereWaveFunction::operator ()(RealVector& x)
{
  int i, j, offset, alpha;
  this->OrbitalValues = (*Orbitals)(x);
  Complex tmp;
  // evaluate single particle Jastrow factors
  for (i=0;i<this->NbrParticles;i++)
    {
      Ji[i]=1.0;
      for(j=0;j<i;j++) Ji[i] *= Orbitals->JastrowFactorElement(i,j);
      for(j=i+1;j<NbrParticles;j++) Ji[i] *= Orbitals->JastrowFactorElement(i,j);
    }  
  // evaluate sums over orbitals m for each LL:
  for (i=0;i<this->NbrParticles;i++)
    for(j=0;j<i;j++)
      {
	alpha=0;
	for (int n=0;n<this->NbrLandauLevels;n++)
	  {
	    tmp=0.0;
	    offset=2*n*(n+this->AbsEffectiveFlux+1)+this->AbsEffectiveFlux;	    
	    for (int m2=-AbsEffectiveFlux-2*n; m2<=AbsEffectiveFlux+2*n;m2+=2)
	      {
		//offset-alpha gives Phi[] with -m 
		tmp+=fsgn((m2+AbsEffectiveFlux)/2)*OrbitalValues[alpha][i]*OrbitalValues[offset-alpha][j];
		//cout << "matching up " << alpha << " with " << offset-alpha<<" sign: "<<fsgn((m2+AbsEffectiveFlux)/2)<<endl;
		alpha++;
	      }
	    this->gAlpha[n][i*this->NbrParticles+j] = tmp;
	  }
      }
  // initialize Slater determinant (or Pfaffian matrix)
  for (int i=0;i<this->NbrParticles;++i)
    {
      for(int j=0;j<i;++j)
	{
	  tmp=0.0;
	  for (int n=0; n<this->NbrLandauLevels; ++n)
	    tmp+=this->CFCoefficients[n]*this->gAlpha[n][i*this->NbrParticles+j];	    
	  
	  Slater->SetMatrixElement(i,j, this->ElementNorm*this->Ji[i]*this->Ji[j]
				   *(MooreReadCoefficient/Orbitals->JastrowFactorElement(i,j) + tmp));
	}
    }  
  //cout << *Slater << endl;
  return Slater->Pfaffian();
}

// normalize the wave-function to one for the given particle positions
// x = point where the function has to be evaluated
void PairedCFOnSphereWaveFunction::AdaptNorm(RealVector& x)
{
  double det=Norm((*this)(x));
  while ((det<.1)||(det>50.0))
    {
      //cout <<"N'="<< this->ElementNorm << " det="<<det<<endl;
      if (det>1e300) 
	this->ElementNorm*= pow((double)1.0e-300,(double)2.0/this->NbrParticles);
      else if (det==0.0) 
	this->ElementNorm*= pow((double)1.0e300,(double)2.0/this->NbrParticles);
      else 
	this->ElementNorm*=pow(det,(double)-2.0/this->NbrParticles);
      det=Norm((*this)(x));
      //cout <<"N'="<< this->ElementNorm << endl;
    }
}
  
