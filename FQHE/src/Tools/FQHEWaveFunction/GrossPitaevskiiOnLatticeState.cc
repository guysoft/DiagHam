////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//                 class of quatum Hall hamiltonian associated                //
//   to particles with contact interactions on a lattice in magnetic field    //
//                                                                            //
//                      last modification : 13/02/2008                        //
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
#include "GrossPitaevskiiOnLatticeState.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "Tools/NewUnconstrainedOptimizsation.h"
#include "GeneralTools/StringTools.h"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sys/time.h>
using std::cout;
using std::endl;

using std::ostream;

// switch for debugging output:
//#define DEBUG_OUTPUT



// constructor for contact interactions on a square lattice
//
// particles = Hilbert space associated to the system
// nbrStates = number of quantum states
// oneParticleTerms = file describing single particle terms
// twoParticleTerms = file describing two-particle terms
GrossPitaevskiiOnLatticeState::GrossPitaevskiiOnLatticeState(int nbrStates, const char* oneParticleTerms, const char* twoParticleTerms, LatticePhases *latticeGeometry, RealVector *variationalParameters)
{
  if (oneParticleTerms!=NULL)
    {
      this->OneParticleTerms = new char[strlen(oneParticleTerms)+2];
      strcpy(this->OneParticleTerms,oneParticleTerms);
    }
  else this->OneParticleTerms=NULL;
  if (twoParticleTerms!=NULL)
    {
      this->TwoParticleTerms = new char[strlen(twoParticleTerms)+2];
      strcpy(this->TwoParticleTerms,twoParticleTerms);
    }
  else this->TwoParticleTerms=NULL;
  this->LatticeGeometry=latticeGeometry;
  this->NbrSites = nbrStates;
  this->EvaluateInteractionFactors();
  if (variationalParameters!=NULL)
    {
      this->VariationalParameters = RealVector(*variationalParameters,true);
      if (variationalParameters->GetVectorDimension()!=2*NbrSites)
	{
	  cout << "Attention, inconsistent number of variational parameters"<<endl;
	  this->VariationalParameters.Resize(2*NbrSites);
	}
    }
  else
    this->VariationalParameters = RealVector(2*NbrSites);
  VariationalCoefficients.Resize(2*NbrSites);
  timeval RandomTime;
  gettimeofday (&(RandomTime), 0);
  this->RandomNumbers = new NumRecRandomGenerator(RandomTime.tv_sec);
  this->ChemicalPotential=1.0;
}

// destructor
//
GrossPitaevskiiOnLatticeState::~GrossPitaevskiiOnLatticeState()
{
  if (NbrHoppingTerms>0)
    {
      delete [] this->HoppingTerms;
      delete [] this->KineticQi;
      delete [] this->KineticQf;
    }
  if (NbrInteractionFactors>0)
    {
      delete [] this->InteractionFactors;
      delete [] this->Q1Value;
      delete [] this->Q2Value;
      delete [] this->Q3Value;
      delete [] this->Q4Value;
    }
  if (OneParticleTerms!=NULL)
    delete [] this->OneParticleTerms;
  if (this->TwoParticleTerms!=NULL)
    delete [] this->TwoParticleTerms;

}


void GrossPitaevskiiOnLatticeState::SetVariationalParameters(RealVector &variationalParameters)
{
  if (variationalParameters.GetVectorDimension()!=2*NbrSites)
    {
      cout << "Attention, inconsistent number of variational parameters"<<endl;
      variationalParameters.Resize(2*NbrSites);
    }
  this->VariationalParameters.Copy(variationalParameters);
}

// set parameters to a random initial distribution (random phase)
// amplitude = amplitude determining the density
void GrossPitaevskiiOnLatticeState::SetToRandomPhase(double amplitude)
{
  Complex TmpC;
  
  for (int i=0; i<NbrSites; ++i)
    {
      this->VariationalParameters[2*i]=sqrt(amplitude);
      this->VariationalParameters[2*i+1]=2.0*M_PI*RandomNumbers->GetRealRandomNumber();
    }
}

// get expectation value of the energy
double GrossPitaevskiiOnLatticeState::GetEnergy()
{
  Complex Result=0.0;
  double rho;
  double IntegratedDensity=0.0;
  for (int i=0; i<NbrSites; ++i)
    {
      rho=VariationalParameters[2*i];
      rho*=rho;
      VariationalCoefficients[i]=Polar(rho,VariationalParameters[2*i+1]);
      IntegratedDensity+=rho*rho;
    }
  for (int i=0; i<NbrInteractionFactors; ++i)
    {
      Result+=(InteractionFactors[i]*Conj(VariationalCoefficients[Q1Value[i]])*Conj(VariationalCoefficients[Q2Value[i]])*
		   VariationalCoefficients[Q3Value[i]]*VariationalCoefficients[Q4Value[i]]);
    }
  for (int i=0; i<NbrHoppingTerms; ++i)
    {
      Result+=(HoppingTerms[i]*Conj(VariationalCoefficients[KineticQf[i]])*VariationalCoefficients[KineticQi[i]]);
    }
  //cout << "Interaction Energy="<<Result;
  Result+=this->ChemicalPotential*IntegratedDensity;
  //cout << " density="<<IntegratedDensity <<" chem pot" <<this->ChemicalPotential<<" final Energy="<<Result<<endl;
  return Result.Re;
}

// get the total number of particles corresponding to the last configuration
double GrossPitaevskiiOnLatticeState::GetNbrParticles()
{
  double rho, IntegratedDensity=0.0;
  for (int i=0; i<NbrSites; ++i)
    {
      rho=VariationalParameters[2*i];
      rho*=rho;
      IntegratedDensity+=rho*rho;
    }
  return IntegratedDensity;
}



// target function for optimizer routine:
// version for 1st parameter fixed to constant
double GrossPitaevskiiOnLatticeState::EvaluateEnergy(int nbrParameters, double *x)
{
  for (int i=0; i<2*NbrSites; ++i)
    if (this->VariationalParameters[i]!=x[i])
      this->VariationalParameters[i]=x[i];
  ++this->NbrEvaluations;
  if ((this->NbrEvaluations<20)
      || ((this->NbrEvaluations<1000)&&(this->NbrEvaluations%10==0))
      || (this->NbrEvaluations%100==0))
    cout << ".";
  cout.flush();
  return this->GetEnergy();
}

// optimize wavefunction starting from present settings of VariationalParameters
// tolerance = final tolerance on the variational parameters
// maxIter = maximal number of function evaluations
//
double GrossPitaevskiiOnLatticeState::Optimize(double tolerance, int maxIter)
{
  double InitialStepSize=1.0;
  int EffectiveNbrVariationalParameters = 2*NbrSites;
  cout << "Starting Optimization ";
  this->NbrEvaluations=0;
  int NbrPoints = 2 * EffectiveNbrVariationalParameters + 1, rnf;
  double Result;
  double *Work = new double[(NbrPoints+13)*(NbrPoints+EffectiveNbrVariationalParameters)
			    + 3*EffectiveNbrVariationalParameters*(EffectiveNbrVariationalParameters+3)/2 + 12];
  // passing parameter vector to optimizer as vector indexed from 1, not 0:
  double *x = &this->VariationalParameters[0];
  double (GrossPitaevskiiOnLatticeState::*TargetFunction)(int, double*)=&GrossPitaevskiiOnLatticeState::EvaluateEnergy;
  GrossPitaevskiiOnLatticeState *TargetObject=this;
  Result = NewUOA::newuoa(EffectiveNbrVariationalParameters, NbrPoints, x, InitialStepSize,
			  tolerance, &rnf, maxIter, Work, TargetObject, TargetFunction);
  cout << endl << "total: "<<NbrEvaluations<< " evaluations"<<endl;
  delete [] Work;
  return Result;
  
}


// evaluate all interaction factors
//   
void GrossPitaevskiiOnLatticeState::EvaluateInteractionFactors()
{  
  if (this->OneParticleTerms!=0)
    {
      if (!IsFile(this->OneParticleTerms))
	{
	  cout << "Could not read file with single particle interactions "<<this->OneParticleTerms<<endl;
	  exit(1);
	}
      MultiColumnASCIIFile Parser;
      if ((Parser.Parse(this->OneParticleTerms))&&((Parser.GetNbrColumns()==4)||(Parser.GetNbrColumns()==3)))
	{
	  this->NbrHoppingTerms = Parser.GetNbrLines();
	  this->HoppingTerms = new Complex[NbrHoppingTerms];
	  this->KineticQf = Parser.GetAsIntegerArray (0);
	  this->KineticQi = Parser.GetAsIntegerArray (1);
	  double *TmpRe = Parser.GetAsDoubleArray (2);
	  if (Parser.GetNbrColumns()==4)
	    {
	      double *TmpIm = Parser.GetAsDoubleArray (3);
	      for (int i=0; i<NbrHoppingTerms; ++i)
		this->HoppingTerms[i]=Complex(TmpRe[i],TmpIm[i]);
	      delete [] TmpIm; 
	    }
	  else
	    {
	      for (int i=0; i<NbrHoppingTerms; ++i)
		this->HoppingTerms[i]=Complex(TmpRe[i],0.0);
	    }
	  delete [] TmpRe;
	}
      else
	{
	  cout << "Error parsing single particle interactions "<<this->OneParticleTerms<<
	    " (3 or 4 columns required: Qf, Qi, Re(M), [Im(M)])"<<endl;
	  exit(1);
	}
    }
  else
    {
      this->NbrHoppingTerms = 0;
    }


  if ((this->LatticeGeometry!=NULL)&&(this->LatticeGeometry->HaveOneParticlePotentials()))
    {
      int OldNumberTerms=this->NbrHoppingTerms;
      this->NbrHoppingTerms+=LatticeGeometry->GetNbrLocalPotentials();
      int *NewQi=new int[this->NbrHoppingTerms];
      int *NewQf=new int[this->NbrHoppingTerms];
      Complex *NewHoppingTerms=new Complex[this->NbrHoppingTerms];
      if (OldNumberTerms>0)
	{
	  for (int i=0; i<OldNumberTerms; ++i)
	    {
	      NewQi[i] = KineticQi[i];
	      NewQf[i] = KineticQi[i];
	      NewHoppingTerms[i] = HoppingTerms[i];
	    }
	  delete[]KineticQi;
	  delete[]KineticQf;
	  delete[]HoppingTerms;
	}
      this->KineticQi=NewQi;
      this->KineticQf=NewQf;
      this->HoppingTerms=NewHoppingTerms;
      cout << "Adding one particle potentials in Hamiltonian"<<endl;
      int NbrPotentials;
      int *PotentialPositions;
      double *Potentials = LatticeGeometry->GetOneParticlePotentials(NbrPotentials, PotentialPositions);
      for (int n=0; n<NbrPotentials; ++n)
	{
	  KineticQi[OldNumberTerms] = PotentialPositions[n];
	  KineticQf[OldNumberTerms] = KineticQi[OldNumberTerms];
	  HoppingTerms[OldNumberTerms] = Potentials[n];
#ifdef DEBUG_OUTPUT
	  cout << "H["<<KineticQi[OldNumberTerms]<<"->"<<KineticQf[OldNumberTerms]<<"]="<<HoppingTerms[OldNumberTerms]<<endl;
#endif
	  ++OldNumberTerms;
	}
    }

  // read two-particle interactions
  if (this->TwoParticleTerms!=0)
    {
      if (!IsFile(this->TwoParticleTerms))
	{
	  cout << "Could not read file with 2-particle interactions "<<this->TwoParticleTerms<<endl;
	  exit(1);
	}
      MultiColumnASCIIFile Parser;
      if ((Parser.Parse(this->TwoParticleTerms))&&((Parser.GetNbrColumns()==5)||(Parser.GetNbrColumns()==6)))
	{
	  this->NbrInteractionFactors = Parser.GetNbrLines();
	  this->InteractionFactors = new Complex[NbrInteractionFactors];
	  this->Q1Value = Parser.GetAsIntegerArray (0);
	  this->Q2Value = Parser.GetAsIntegerArray (1);
	  this->Q3Value = Parser.GetAsIntegerArray (2);
	  this->Q4Value = Parser.GetAsIntegerArray (3);
	  double *TmpRe = Parser.GetAsDoubleArray (4);
	  if (Parser.GetNbrColumns()==6)
	    {
	      double *TmpIm = Parser.GetAsDoubleArray (5);
	      for (int i=0; i<NbrInteractionFactors; ++i)
		this->InteractionFactors[i]=Complex(TmpRe[i],TmpIm[i]);
	      delete [] TmpIm; 
	    }
	  else
	    {
	      for (int i=0; i<NbrInteractionFactors; ++i)
		this->InteractionFactors[i]=Complex(TmpRe[i],0.0);
	    }
	  delete [] TmpRe;

	  // test symmetries
	  int oldQ1=Q1Value[0], oldQ2=Q1Value[0];
	  bool HaveLargerQ1=false, HaveSmallerQ1=false;
	  bool HaveLargerQ3=false, HaveSmallerQ3=false;
	  int Pos=0;
	  while (Pos<NbrInteractionFactors)
	    {
	      while ((Pos<NbrInteractionFactors)&&(Q1Value[Pos]==oldQ1)&&(Q2Value[Pos]==oldQ2))
		{
		  if (Q1Value[Pos]>Q2Value[Pos]) HaveLargerQ1=true;
		  else if (Q1Value[Pos]<Q2Value[Pos]) HaveSmallerQ1=true;
		  if (Q3Value[Pos]>Q4Value[Pos]) HaveLargerQ3=true;
		  else if (Q3Value[Pos]<Q4Value[Pos]) HaveSmallerQ3=true;
		  ++Pos;
		}
	      if (Pos<NbrInteractionFactors)
		{
		  if (Q1Value[Pos]>Q2Value[Pos]) HaveLargerQ1=true;
		  else if (Q1Value[Pos]<Q2Value[Pos]) HaveSmallerQ1=true;
		  if (Q3Value[Pos]>Q4Value[Pos]) HaveLargerQ3=true;
		  else if (Q3Value[Pos]<Q4Value[Pos]) HaveSmallerQ3=true;
		  oldQ1=Q1Value[Pos];
		  oldQ2=Q2Value[Pos];
		  ++Pos;
		}
	    }
	  
	  bool Q12Symmetry=false, Q34Symmetry=false;
	  if (HaveSmallerQ1^HaveLargerQ1)
	    Q12Symmetry=true;
	  if (HaveSmallerQ3^HaveLargerQ3)
	    Q34Symmetry=true;
	  if ((Q12Symmetry)||(Q34Symmetry))
	    cout << "Assuming symmetry in";
	  if (Q12Symmetry) cout << " Q12";
	  if (Q34Symmetry) cout << " Q34";
	  cout<<endl;

	  if (Q12Symmetry)
	    for (int i=0; i<NbrInteractionFactors; ++i)
	      if (Q1Value[Pos]!=Q2Value[Pos]) InteractionFactors[i]*=2.0;
	  if (Q34Symmetry)
	    for (int i=0; i<NbrInteractionFactors; ++i)
	      if (Q3Value[Pos]!=Q4Value[Pos]) InteractionFactors[i]*=2.0;

	}
      else
	{
	  cout << "Error parsing two-particle interactions "<<this->TwoParticleTerms<<
	    " (5 or 6 columns required: Q1, Q2, Q3, Q4, Re(M), [Im(M)])"<<endl;
	  exit(1);
	}
    }
  else
    {
      // we have no general four-particle interactions:     
      this->NbrInteractionFactors=0;
      cout << "No two-particle interactions"<<endl;
    }
}
