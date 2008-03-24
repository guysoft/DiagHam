////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of Moore Read state wave function on sphere             //
//                                                                            //
//                        last modification : 19/09/2004                      //
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
#include "Tools/FQHEWaveFunction/AdvancedMooreReadOnSphereWaveFunction.h"
#include "MathTools/BinomialCoefficients.h"
#include "Vector/RealVector.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;


// constructor
//
// nbrParticlesPerCluster = number of particles per cluster (=N/2)
// fermionicStatistics = flag indicating whether the pfaffian should be multiplied by a squared Jastrow Factor
AdvancedMooreReadOnSphereWaveFunction::AdvancedMooreReadOnSphereWaveFunction(int nbrParticlesPerCluster, bool fermionicStatistics)
{
  this->NbrParticles = 2*nbrParticlesPerCluster;
  this->ClusterSize = nbrParticlesPerCluster;
  this->NbrClusters = 2;
  this->FermionicStatistics = fermionicStatistics;
  this->SpinorUCoordinates = new Complex[this->NbrParticles];
  this->SpinorVCoordinates = new Complex[this->NbrParticles];
  this->JastrowFactorElements = new Complex*[this->NbrParticles];
  this->JastrowFactorSquares = new Complex*[this->NbrParticles];
  for (int i=0; i<NbrParticles; ++i)
    {
      this->JastrowFactorElements[i] = new Complex[this->NbrParticles];
      this->JastrowFactorSquares[i] = new Complex[this->NbrParticles];
    }
  this->EvaluatePermutations();
  this->Flag.Initialize();
}

// copy constructor
//
// function = reference on the wave function to copy

AdvancedMooreReadOnSphereWaveFunction::AdvancedMooreReadOnSphereWaveFunction(const AdvancedMooreReadOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->ClusterSize = function.ClusterSize;
  this->NbrClusters = function.NbrClusters;
  this->FermionicStatistics = function.FermionicStatistics;
  this->Permutations = function.Permutations;
  this->NbrPermutations = function.NbrPermutations;
  this->WeightOfPermutations = function.WeightOfPermutations;
  this->JastrowFactorElements = function.JastrowFactorElements;
  this->JastrowFactorSquares = function.JastrowFactorSquares;
  this->Flag = function.Flag;
}

// destructor
//

AdvancedMooreReadOnSphereWaveFunction::~AdvancedMooreReadOnSphereWaveFunction()
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      for (unsigned i = 0; i < this->NbrPermutations; ++i)
	delete[] this->Permutations[i];
      delete[] this->Permutations;
      delete[] this->SpinorUCoordinates;
      delete[] this->SpinorVCoordinates;
      for (int i=0; i<NbrParticles; ++i)
	{
	  delete [] this->JastrowFactorSquares[i];
	  delete [] this->JastrowFactorElements[i];
	}
      delete [] this->JastrowFactorSquares;
      delete [] this->JastrowFactorElements;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* AdvancedMooreReadOnSphereWaveFunction::Clone ()
{
  return new AdvancedMooreReadOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex AdvancedMooreReadOnSphereWaveFunction::operator ()(RealVector& x)
{  
  // for (int i = 0; i < this->NbrParticles; ++i)
//     {
//       SpinorUCoordinates[i].Re = cos(0.5 * x[i << 1]);
//       SpinorUCoordinates[i].Im = SpinorUCoordinates[i].Re;
//       SpinorUCoordinates[i].Re *= cos(0.5 * x[1 + (i << 1)]);
//       SpinorUCoordinates[i].Im *= sin(0.5 * x[1 + (i << 1)]);
//       SpinorVCoordinates[i].Re = sin(0.5 * x[i << 1]);
//       SpinorVCoordinates[i].Im = SpinorVCoordinates[i].Re;
//       SpinorVCoordinates[i].Re *= cos(0.5 * x[1 + (i << 1)]);
//       SpinorVCoordinates[i].Im *= -sin(0.5 * x[1 + (i << 1)]);
//     }

  // CalculateSpinors
  double s,c;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = cos(0.5 * x[i << 1]);
      this->SpinorUCoordinates[i].Im = this->SpinorUCoordinates[i].Re;
      this->SpinorUCoordinates[i].Re *= (c=cos(0.5 * x[1 + (i << 1)]));
      this->SpinorUCoordinates[i].Im *= -(s=sin(0.5 * x[1 + (i << 1)]));
      this->SpinorVCoordinates[i].Re = sin(0.5 * x[i << 1]);
      this->SpinorVCoordinates[i].Im = this->SpinorVCoordinates[i].Re;
      this->SpinorVCoordinates[i].Re *= c;
      this->SpinorVCoordinates[i].Im *= s;
      //cout << "U["<<i<<"]="<<SpinorUCoordinates[i]<<", "<< "V["<<i<<"]="<<SpinorVCoordinates[i]<<endl;
    }

  return this->ComplexEvaluations();
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex AdvancedMooreReadOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{  
  // Import from spinors
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = uv.Re(2*i);
      this->SpinorUCoordinates[i].Im = uv.Im(2*i);
      this->SpinorVCoordinates[i].Re = uv.Re(2*i+1);
      this->SpinorVCoordinates[i].Im = uv.Im(2*i+1);
    }

  return this->ComplexEvaluations();
}  


// evaluate permutations required for the Moore-Read state evaluation
// using: two symmetric blocs, only permutations changing particles
//        between blocks are required

void AdvancedMooreReadOnSphereWaveFunction::EvaluatePermutations()
{
  double totalWeight=0.0;
  BinomialCoefficients bico(NbrParticles);
  this->NbrPermutations = this->ClusterSize+1;
  this->WeightOfPermutations = new double[NbrPermutations];
  this->Permutations = new unsigned*[NbrPermutations];
  for (unsigned i=0; i<NbrPermutations; ++i)
    this->Permutations[i] = new unsigned[NbrParticles];
  for (int i=0; i<NbrParticles; ++i)
    this->Permutations[0][i] = i;
  this->WeightOfPermutations[0] = 1.0/(double)bico(NbrParticles,ClusterSize);
  totalWeight+=WeightOfPermutations[0];
  for (int i=0; i<ClusterSize; ++i)
    {
      for (int j=0; j<NbrParticles; ++j)
	this->Permutations[i+1][j] = this->Permutations[i][j];
      this->Permutations[i+1][i] = i+ClusterSize;
      this->Permutations[i+1][i+ClusterSize] = i;
      this->WeightOfPermutations[i+1] = (double)bico(ClusterSize,i+1)*(double)bico(ClusterSize,i+1)/(double)bico(NbrParticles,ClusterSize);
      totalWeight+=WeightOfPermutations[i+1];
      cout << "Permutation "<<i+1<<": [ "<<Permutations[i+1][0];
      for (int k=1; k<NbrParticles; ++k) cout<<", "<<Permutations[i+1][k];
      cout << "] has weight "<<WeightOfPermutations[i+1]<<endl;
    }
  cout << "TotalWeight="<<totalWeight<<endl;
  return;
}

// perform complex part of calculations
// uses internal spinor coordinates as input
//
Complex AdvancedMooreReadOnSphereWaveFunction::ComplexEvaluations()
{
  Complex JA, JB, Tmp;

  for (int i = 0; i < this->NbrParticles; ++i)
    for (int j = 0; j < i; ++j)
      {
	Tmp = ((this->SpinorUCoordinates[i] * this->SpinorVCoordinates[j]) - (this->SpinorUCoordinates[j] * this->SpinorVCoordinates[i]));
	JastrowFactorElements[i][j] = Tmp;
	JastrowFactorElements[j][i] = -Tmp;
      }

  unsigned *TmpP;
  Complex Value(0.0,0.0);
  cout << "Evaluating function"<<endl;
  for (unsigned i=0; i<NbrPermutations; ++i)
    {
      TmpP=this->Permutations[i];
      JA=1.0;
      JB=1.0;
      for (int k=1; k<ClusterSize; ++k)
	for (int j=0; j<k; ++j)
	{
	  JA*=JastrowFactorElements[TmpP[k]][TmpP[j]];
	  cout << "JA*=JastrowFactorElements("<<TmpP[k]<<","<<TmpP[j]<<")"<<endl;
	  JB*=JastrowFactorElements[TmpP[k+ClusterSize]][TmpP[j+ClusterSize]];
	  cout << "JB*=JastrowFactorElements("<<TmpP[k+ClusterSize]<<","<<TmpP[j+ClusterSize]<<")"<<endl;
	}
      cout << "Contribution: "<<JA*JA*JB*JB<<endl;
      cout << "Weight:" <<WeightOfPermutations[i]<<endl;
      Value += WeightOfPermutations[i]*JA*JA*JB*JB;
    }
  cout << "Symmetric Part: "<<Value<<endl;
  if (this->FermionicStatistics)
    {
      for (int i = 1; i < this->NbrParticles; ++i)
	for (int j = 0; j < i; ++j)
	  Value *=  JastrowFactorElements[i][j];
    }
  cout << "Value ="<<Value<<endl;
  
  return Value;
}
