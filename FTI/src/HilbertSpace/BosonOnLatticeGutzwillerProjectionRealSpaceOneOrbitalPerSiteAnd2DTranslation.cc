////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of bosons on lattice			              //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Antoine Sterdyniak                    //
//                                                                            //
//                        last modification : 12/09/2014                      //
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
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"

#include <math.h>
#include <cstdlib>
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation::BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// nbrSite = number of sites
// xMomentum = momentum sector in the x direction
// maxYMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction 
// memory = amount of memory granted for precalculations

BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation::BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation (int nbrBosons,  int lx, int ly, int xMomentum, int  maxXMomentum,
										      int yMomentum, int maxYMomentum, unsigned long memory):
 BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation (nbrBosons,lx*ly,xMomentum,maxXMomentum,yMomentum,maxYMomentum,memory), Lx(lx),Ly(ly)
{  
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation::BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation(const BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation & bosons) :  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation(bosons)
{
this->Lx = bosons.Lx;
this->Ly = bosons.Ly;
}

// destructor
//

BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation::~BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation ()
{
}


// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation & BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation::operator = (const BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation& bosons)
{
  BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::operator = (bosons);
  this->Lx = bosons.Lx;
  this->Ly = bosons.Ly;
  return (*this);
}


// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation::Clone()
{
  return new BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation(*this);
}






void BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation::GetCompositeFermionWavefunction(ComplexVector & trialState, ComplexMatrix & jastrowEigenVecs,ComplexMatrix & cFEigenVecs, double phaseTranslationX)
{
  int Nx = this->Lx / this->MaxXMomentum;
  int Ny = this->Ly / this->MaxYMomentum;
  int PositionX;
  int PositionY;
  Complex ** ExponentialFactors = new Complex*[this->MaxXMomentum];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    { 
      ExponentialFactors[i] = new Complex[this->MaxYMomentum];
      for (int j = 0; j < this->MaxYMomentum; ++j)
	{ 
	  ExponentialFactors[i][j] = Phase(2.0 * M_PI * ((this->XMomentum * ((double) i) / ((double) this->MaxXMomentum))
							       + (this->YMomentum * ((double) j) / ((double) this->MaxYMomentum))));
	}
    }

#ifdef __LAPACK__
  ComplexLapackDeterminant SlaterCF(NbrBosons);
  ComplexLapackDeterminant SlaterJastrow(NbrBosons);
#else
  ComplexMatrix SlaterCF(NbrBosons, NbrBosons);
  ComplexMatrix SlaterJastrow(NbrBosons, NbrBosons);
#endif

 int NbrTranslation;
 unsigned long * TemporaryState = new unsigned long [this->NbrBosons];
 
 for(int i = 0; i < this->HilbertSpaceDimension ; i++)
 {
    cout << this->StateDescription[i]<<" " <<this->NbrStateInOrbit[i]<<endl;
    NbrTranslation = 0;
    unsigned long  TmpStateDescription  = this->StateDescription[i];
    unsigned long  TmpStateDescription2  = this->StateDescription[i];
    int m=0;
    do{
      TmpStateDescription2 = TmpStateDescription;
        int n=0;
	do{  
     
          this->ConvertToMonomial(TmpStateDescription2,TemporaryState);
          for (int p = 0; p < NbrBosons; ++p)
          {
         	 cout <<TemporaryState[p]<<" ";
     	  }
	  cout <<endl;
          this->GetPositionSum(TemporaryState, PositionX, PositionY);

  for (int p = 0; p < NbrBosons; ++p)
  {

   for (int q = 0; q < NbrBosons; ++q)
	{
		  // need to consider proper ordering of matrix elements
		  // in Hilbert-space, largest quantum number q corresponds to position 0!

     		  SlaterCF.SetMatrixElement(p,q,cFEigenVecs[p][TemporaryState[q]]);
 		  SlaterJastrow.SetMatrixElement(p,q,jastrowEigenVecs[p][TemporaryState[q]]);
	}	      
       }
        
        trialState[i] +=  ExponentialFactors[n][m] * SlaterCF.Determinant() * SlaterJastrow.Determinant()* Phase(1.0*PositionY*n*phaseTranslationX);
        ++NbrTranslation;
        this->ApplySingleXTranslation(TmpStateDescription2);
       ++n;        
	} while((n < this->MaxXMomentum) && (TmpStateDescription2 != TmpStateDescription));
      this->ApplySingleYTranslation(TmpStateDescription);      
       ++m;  
    }while ((m < this->MaxYMomentum) && (TmpStateDescription !=  this->StateDescription[i]));

   
  trialState[i] /= sqrt(this->NbrStateInOrbit[i]);
  if(NbrTranslation != this->NbrStateInOrbit[i])
   {
	cout <<"Wrong Number of Translation for state !"  << i <<endl;
        cout <<NbrTranslation<< " " << this->NbrStateInOrbit[i]<<endl;
  }
   
}
   for (int i = 0; i < this->MaxXMomentum; ++i)
    { 
      delete [] ExponentialFactors[i];
    }
   delete []  ExponentialFactors;
}



