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
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "Tools/NewUnconstrainedOptimizsation.h"
#include <iostream>
#include <sys/time.h>

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
    {
      this->MaxOccupation = 1;
      this->NbrVariationalParameters = 3*this->NbrSites;
    }
  else
    {
      this->MaxOccupation = NbrParticles;
      this->NbrVariationalParameters = (2*this->NbrParticles+1)*this->NbrSites;
    }
  if (variationalParameters!=NULL)
    {
      this->VariationalParameters = RealVector(*variationalParameters,true);
      if (variationalParameters->GetVectorDimension()!=NbrVariationalParameters)
	{
	  cout << "Attention, inconsistent number of variational parameters"<<endl;
	  this->VariationalParameters.Resize(NbrVariationalParameters);
	}
    }
  else
    this->VariationalParameters = RealVector(NbrVariationalParameters);
  timeval RandomTime;
  gettimeofday (&(RandomTime), 0);
  this->RandomNumbers = new NumRecRandomGenerator(RandomTime.tv_sec);
  this->Dim = Space->GetHilbertSpaceDimension();
  this->Hamiltonian=NULL;
  this->Architecture=NULL;
}

// destructor
GutzwillerOnLatticeWaveFunction::~GutzwillerOnLatticeWaveFunction()
{
  delete this->RandomNumbers;
}

// get Many-Body state
// return = resultingState
ComplexVector & GutzwillerOnLatticeWaveFunction::GetGutzwillerWaveFunction()
{
  this->TargetVector.Resize(Dim);
  this->TargetVector.ClearVector();

  // call main recursion
  this->Product(NbrSites-1, 0, 0x0ul, 1.0);
  if (TargetVector.Norm()==0.0)
    cout << "Attention, obtained state with zero norm"<<endl;
  this->TargetVector/=this->TargetVector.Norm();
  //   cout <<"Test norm: "<<TargetVector.Norm()<<endl;
  return this->TargetVector;
}


void GutzwillerOnLatticeWaveFunction::SetVariationalParameters(RealVector &variationalParameters)
{
  if (variationalParameters.GetVectorDimension()!=NbrVariationalParameters)
    {
      cout << "Attention, inconsistent number of variational parameters"<<endl;
      variationalParameters.Resize(NbrVariationalParameters);
    }
  this->VariationalParameters.Copy(variationalParameters);
}

// set parameters to a random initial distribution (random phase)
void GutzwillerOnLatticeWaveFunction::SetToRandomPhase()
{
  Complex TmpC;
  
  for (int i=0; i<NbrSites; ++i)
    this->VariationalParameters[i]=1.0;
  for (int i=0; i<NbrSites; ++i)
    {
      TmpC=Polar((double)NbrParticles/NbrSites,2.0*M_PI*RandomNumbers->GetRealRandomNumber());
      this->VariationalParameters[NbrSites+2*i]=TmpC.Re;
      this->VariationalParameters[NbrSites+2*i+1]=TmpC.Im;
      for (int n=2; n<MaxOccupation; ++n)
	for (int i=0; i<NbrSites; ++i)
	  {
	    this->VariationalParameters[(2*n-1)*NbrSites+2*i]=TmpC.Re/n/n;
	    this->VariationalParameters[(2*n-1)*NbrSites+2*i+1]=TmpC.Im/n/n;
	  }
    }
}

// define a Hamiltonian to enable immediate evaluation of the energy
void GutzwillerOnLatticeWaveFunction::SetHamiltonian(AbstractQHEOnLatticeHamiltonian *hamiltonian)
{
  this->Hamiltonian=hamiltonian;
}
  
// define an architecture to enable multi-processor operations
void GutzwillerOnLatticeWaveFunction::SetArchitecture(AbstractArchitecture *architecture)
{
  this->Architecture=architecture;
}

// get expectation value of the energy
double GutzwillerOnLatticeWaveFunction::GetEnergy()
{
  if (this->Hamiltonian==NULL)
    {
      cout << "Please define a Hamiltonian first" << endl;
      exit(1);
    }
  if (this->Architecture==NULL)
    {
      cout << "Please define an Architecture first" << endl;
      exit(1);
    }
  this->GetGutzwillerWaveFunction();
  ComplexVector TmpState(this->Space->GetHilbertSpaceDimension());
  VectorHamiltonianMultiplyOperation Operation (this->Hamiltonian, &(this->TargetVector), &TmpState);
  Operation.ApplyOperation(this->Architecture);
  return Real(TargetVector*TmpState);
}




// main recursion to calculate State \prod_i (\sum \chi_i + \psi_i a^\dagger_i) |state>
// nextQ = value quantum number in next operator to be applied
// nbrBosons = number of bosons already in state
// state = state to be acted upon
// prefactor = previous coefficients applied to state
// in last stage of recursion, writes to this->TargetVector
void GutzwillerOnLatticeWaveFunction::Product (int nextQ, int nbrBosons, unsigned long state, Complex prefactor)
{
  // cout << "Calling: Product ("<<nextQ<<", "<< nbrBosons<<", "<<state<<", "<< prefactor<<")"<<endl;
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
	  Product(nextQ-1, nbrBosons+NbrPlaced, ResultingState, prefactor* Complex(VariationalParameters[nextQ+2*NbrPlaced*NbrSites],
										   VariationalParameters[nextQ+2*NbrPlaced*NbrSites]+1));
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
	  TargetVector[Index]+= prefactor*Complex(VariationalParameters[nextQ+2*NbrPlaced*NbrSites],
						  VariationalParameters[nextQ+2*NbrPlaced*NbrSites]+1);
	}
    }
}

// target function for optimizer routine:
double GutzwillerOnLatticeWaveFunction::EvaluateEnergy(int nbrParameters, double *x)
{
  if (nbrParameters!=this->NbrVariationalParameters-1)
    {
      cout << "Unexpected error"<<endl;
    }
  for (int i=1; i<this->NbrVariationalParameters; ++i)
    this->VariationalParameters[i]=x[i];
  return this->GetEnergy();
}

// optimize wavefunction starting from present settings of VariationalParameters
// tolerance = final tolerance on the variational parameters
// maxIter = maximal number of function evaluations
//
double GutzwillerOnLatticeWaveFunction::Optimize(double tolerance, int maxIter)
{
  double InitialStepSize=1.0;
  int EffectiveNbrVariationalParameters = NbrVariationalParameters-1;
  cout << "Initial parameters:" << endl << this->VariationalParameters;
  int NbrPoints = 2 * EffectiveNbrVariationalParameters + 1, rnf;
  double Result;
  double *Work = new double[(NbrPoints+13)*(NbrPoints+EffectiveNbrVariationalParameters)
			    + 3*EffectiveNbrVariationalParameters*(EffectiveNbrVariationalParameters+3)/2 + 12];
  // passing parameter vector to optimizer as vector indexed from 1, not 0:
  double *x = &this->VariationalParameters[0];
  double (GutzwillerOnLatticeWaveFunction::*TargetFunction)(int, double*)=&GutzwillerOnLatticeWaveFunction::EvaluateEnergy;
  GutzwillerOnLatticeWaveFunction *TargetObject=this;
  Result = NewUOA::newuoa(EffectiveNbrVariationalParameters, NbrPoints, x, InitialStepSize,
			  tolerance, &rnf, maxIter, Work, TargetObject, TargetFunction);
  cout << "Final parameters:" << endl << this->VariationalParameters;
  delete [] Work;
  return Result;
  
}
