////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                      class of SkyrmionOnSphere state wave function                    //
//                                                                            //
//                        last modification : 20/04/2005                      //
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
#include "Tools/FQHEWaveFunction/SkyrmionOnSphereWaveFunction.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "Options/Options.h"

#include <iostream>
using std::cout;
using std::endl;

// constructor
//
// create object to be initialized from system options
//
SkyrmionOnSphereWaveFunction::SkyrmionOnSphereWaveFunction(AbstractArchitecture* architecture, OptionManager &manager, int nbrParticles, int totalLzMax, int totalLz, int totalSz, QHEWaveFunctionManager *waveFunctionManager, int basisType)
{
  this->Architecture=architecture;
  this->UseExact=true;
  this->NbrParticles=nbrParticles;
  int PolarizedParticles = manager.GetInteger("nbr-particles");
  this->PolarizedLzMax = manager.GetInteger("polarized-lzmax"); 
  this->PolarizedLz = manager.GetInteger("polarized-lz");  
  
  this->AnalyticPolarizedWaveFunction=NULL;
  
  if (manager.GetString("polarized-state")==0)
    {
      if (waveFunctionManager!=NULL)
	{
	  UseExact=false;
	  PolarizedParticles=nbrParticles;
	  AnalyticPolarizedWaveFunction = waveFunctionManager->GetWaveFunction();
	}
      else
	{
	  cout << "An exact state vector is required for Skyrmion::\"polarized state\""<<endl;
	  exit(-1);
	}
    }

  bool Statistics = true;
  
  if (UseExact)
    {
      PolarizedParticles = manager.GetInteger("nbr-particles"); 
      PolarizedLzMax = manager.GetInteger("polarized-lzmax"); 
      PolarizedLz = manager.GetInteger("polarized-lz");
      
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(manager.GetString("polarized-state"),
						       PolarizedParticles, PolarizedLzMax, PolarizedLz, Statistics) == false)
	{
	  cout << "error while retrieving system parameters for polarized state " <<
	    manager.GetString("polarized-state") << endl;
	  exit(-1);
	}
      
      if (PolarizedParticles!=nbrParticles)
	{
	  cout << "Error: polarized and exact state have to have the same number of particles";
	  exit(-1);
	}
      
      if (PolarizedLzMax>totalLzMax)
	{
	  cout << "Error: polarized state has to be at a lower flux than the exact state";
	  exit(-1);
	}
      
      if (PolarizedState.ReadVector (manager.GetString("polarized-state")) == false)
	{
	  cout << "can't open vector file " << manager.GetString("polarized-state") << endl;
	  exit(-1);      
	}
    }


  int NbrBosons = manager.GetInteger("nbr-particles"); 
  this->BosonLzMax = 0;
  this->BosonLz = 0;
  this->BosonSz = 0;
  bool SzSymmetrizedBasis = false;
  bool SzMinusParity = false;
  bool LzSymmetrizedBasis = false;
  bool LzMinusParity = false;
  Statistics = false;  

  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(manager.GetString("bosonic-state"), NbrBosons,
							   BosonLzMax, BosonLz, BosonSz,
							   SzSymmetrizedBasis, SzMinusParity, 
							   LzSymmetrizedBasis, LzMinusParity, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << manager.GetString("bosonic-state") << endl;
      exit(-1);
    }

  if (NbrBosons!=nbrParticles)
    {
      cout << "Error: bosonic and exact state have to have the same number of particles";
      exit(-1);
    }
  
  if (BosonLzMax==0)
    BosonLzMax = totalLzMax-PolarizedLzMax;
  else if (BosonLzMax!=totalLzMax-PolarizedLzMax)
    {
      cout << "Error: total flux has to match: BosonLzMax == TotalLzMax-PolarizedLzMax";
      exit(-1);
    }

  if (BosonLz==0)
    BosonLz = totalLz;
  else if (BosonLz!=totalLz)
    {
      cout << "Error: total angular momentum has to match: BosonLz == TotalLz";
      exit(-1);
    }

  if (BosonSz==0)
    BosonSz = totalSz;
  else if (BosonSz!=totalSz)
    {
      cout << "Error: total spin has to match: BosonSz == TotalSz";
      exit(-1);
    }
  
  this->PolarizedSpace=NULL;
#ifdef __64_BITS__
  if (PolarizedLzMax <= 63)
#else
    if (PolarizedLzMax <= 31)
#endif
      {	
	PolarizedSpace = new FermionOnSphere(NbrParticles, PolarizedLz, PolarizedLzMax);
      }
    else
#ifdef __128_BIT_LONGLONG__
      if (PolarizedLzMax <= 126)
#else
	if (PolarizedLzMax <= 62)
#endif
	  {	    
	    PolarizedSpace = new FermionOnSphereLong(NbrParticles, PolarizedLz, PolarizedLzMax);
	  }
	else
	  {
	    cout << "States of this polarized Hilbert space cannot be represented in a single word." << endl;
	    exit(-1);
	  }

  this->BosonicSpace=new BosonOnSphereWithSpin(NbrParticles, BosonLz, BosonLzMax, BosonSz);

  this->OneBodyBasis = new ParticleOnSphereFunctionBasis(totalLzMax, basisType);
}

// copy constructor
//
// function = reference on the wave function to copy

SkyrmionOnSphereWaveFunction::SkyrmionOnSphereWaveFunction(const SkyrmionOnSphereWaveFunction& function)
{
  this->NbrParticles=function.NbrParticles;
  this->PolarizedLzMax=function.PolarizedLzMax;
  this->PolarizedLz=function.PolarizedLz;
  this->PolarizedState=function.PolarizedState;
  this->BosonicState=function.BosonicState;
  this->PolarizedSpace=function.PolarizedSpace;
  this->BosonicSpace=function.BosonicSpace;
  this->OneBodyBasis=function.OneBodyBasis;
  this->UseExact=function.UseExact;
  this->AnalyticPolarizedWaveFunction=function.AnalyticPolarizedWaveFunction;

}

// destructor
//
SkyrmionOnSphereWaveFunction::~SkyrmionOnSphereWaveFunction()
{
  if (!UseExact)
    delete this->PolarizedSpace;
  delete this->BosonicSpace;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* SkyrmionOnSphereWaveFunction::Clone ()
{
  return new SkyrmionOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex SkyrmionOnSphereWaveFunction::operator ()(RealVector& x)
{
  Complex ValueSpin,ValuePolarized;
  if (UseExact)
    {
      QHEParticleWaveFunctionOperation Operation(PolarizedSpace, &PolarizedState, &x, this->OneBodyBasis, /* TimeCoherence */ -1);
      Operation.ApplyOperation(this->Architecture);      
      ValuePolarized = Operation.GetScalar();
    }
  else
    {
      ValuePolarized = (*(this->AnalyticPolarizedWaveFunction))(x);
    }
  QHEParticleWaveFunctionOperation Operation(BosonicSpace, &BosonicState, &x, this->OneBodyBasis, /* TimeCoherence */ -1);
  Operation.ApplyOperation(this->Architecture);      
  ValueSpin = Operation.GetScalar();
  
  return ValuePolarized*ValueSpin;
}


// add an option group containing all options related to the skyrmion wave functions
//
// manager = pointer to the option manager
void SkyrmionOnSphereWaveFunction::AddSkyrmionOptionGroup(OptionManager &manager, QHEWaveFunctionManager *wfManager)
{  
  OptionGroup* SkyrmionGroup = new OptionGroup ("skyrmion options");
  manager+=SkyrmionGroup;
  
  (*SkyrmionGroup) += new SingleStringOption  ('\n', "polarized-state", "file name of polarized fermionic reference wave function (if omitted using analytic function)",0);
  (*SkyrmionGroup) += new SingleIntegerOption  ('l', "polarized-lzmax", "total number of flux quanta (0 if it has to be guessed from input file name)", 0);  
  (*SkyrmionGroup) += new SingleIntegerOption  ('z', "polarized-lz", "twice the total lz value of the system (0 if it has to be guessed from input file name)", 0);
  (*SkyrmionGroup) += new SingleStringOption  ('\n', "bosonic-state", "file name of spinful bosonic part of wave function",0);

  if (wfManager!=NULL)
    wfManager->AddOptionGroup(&manager);
}
