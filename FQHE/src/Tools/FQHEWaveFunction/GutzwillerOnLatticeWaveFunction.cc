////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                       Copyright (C) 2009 Gunnar Möller                     //
//                                                                            //
//                                                                            //
//           class for calculation of a Gutzwiller state on the lattice       //
//                                                                            //
//                        last modification : 27/10/2009                      //
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


#include "GutzwillerOnLatticeWaveFunction.h"

#include "Hamiltonian/AbstractQHEOnLatticeHamiltonian.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>
using std::cout;
using std::endl;

// constructor
// nbrParticles = particles in the condensate (should match space)
// hardCore = flag indicating whether double occupations may occur or whether hardcore particles are present
// space = target-space of many-body state
// variationalParameters = initial set of trial parameters
GutzwillerOnLatticeWaveFunction::GutzwillerOnLatticeWaveFunction(int nbrParticles, bool hardCore, ParticleOnLattice *space, ComplexVector *variationalParameters)
{
  this->NbrParticles = nbrParticles;
  this->Space = space;
  this->NbrSites = Space->GetNbrSites();
  this->HardCoreFlag = hardCore;
  if (this->HardCoreFlag)
    this->NbrVariationalParameters = 2*this->NbrSites;
  else
    this->NbrVariationalParameters = (this->NbrParticles+1)*this->NbrSites;
  if (variationalParameters!=NULL)
    {
      this->VariationalParameters = ComplexVector(*variationalParameters,true);
      if (variationalParameters->GetVectorDimension()!=NbrVariationalParameters)
	{
	  cout << "Attention, inconsistent number of variational parameters"<<endl;
	  this->VariationalParameters.Resize(NbrVariationalParameters);
	}
    }
  else
    this->VariationalParameters = ComplexVector(NbrVariationalParameters);
  
  this->Dim = Space->GetHilbertSpaceDimension();
  this->Hamiltonian=NULL;
  this->Architecture=NULL;
}

// destructor
GutzwillerOnLatticeWaveFunction::~GutzwillerOnLatticeWaveFunction()
{
}

// get Many-Body state
// return = resultingState
ComplexVector & GutzwillerOnLatticeWaveFunction::GetGutzwillerWaveFunction()
{
  this->TargetVector.Resize(Dim);
  this->TargetVector.ClearVector();

  // call main recursion
  this->Product(NbrSites-1, 0, 0x0ul, 1.0);
  
  this->TargetVector/=this->TargetVector.Norm();
//   cout <<"Test norm: "<<TargetVector.Norm()<<endl;
  return this->TargetVector;
}



// main recursion to calculate State \prod_i (\sum \chi_i + \psi_i a^\dagger_i) |state>
// nextQ = value quantum number in next operator to be applied
// nbrBosons = number of bosons already in state
// state = state to be acted upon
// prefactor = previous coefficients applied to state
// in last stage of recursion, writes to this->TargetVector
void GutzwillerOnLatticeWaveFunction::Product (int nextQ, int nbrBosons, unsigned long state, Complex prefactor)
{
  int Index;
  unsigned long ResultingState;
  double AdFactor;
  if (nextQ>0)
    {
      int NbrPlaced = 0;
      ResultingState = state;
      while (nbrBosons+NbrPlaced<this->NbrParticles)
	{
	  ResultingState = Space->Ad(ResultingState,nextQ,AdFactor);
	  ++NbrPlaced;
	  Product(nextQ-1, nbrBosons+NbrPlaced, ResultingState, prefactor*VariationalParameters[nextQ+NbrPlaced*NbrSites]);
	}
      Product(nextQ-1, nbrBosons, state, prefactor*VariationalParameters[nextQ]);
    }
  else
    {
      int NbrPlaced = 0;
      ResultingState = state;
      while (nbrBosons+NbrPlaced<this->NbrParticles)
	{
	  ResultingState = Space->Ad(ResultingState,nextQ,AdFactor);
	  ++NbrPlaced;
	}
      if ((Index=Space->CarefulFindStateIndex(ResultingState,-1))<Dim)
	{
	  TargetVector[Index]+= prefactor*VariationalParameters[nextQ+NbrPlaced*NbrSites];
	}
    }
}
