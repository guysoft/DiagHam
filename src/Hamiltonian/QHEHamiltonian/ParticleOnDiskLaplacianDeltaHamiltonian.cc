////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated to particles on a disk with         //
//                          laplacian delta interaction                       //
//                                                                            //
//                        last modification : 05/02/2004                      //
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


#include "Hamiltonian/QHEHamiltonian/ParticleOnDiskLaplacianDeltaHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"
#include "MathTools/FactorialCoefficient.h"
#include "Architecture/AbstractArchitecture.h"

#include <iostream>
#include <math.h>
#include <stdlib.h>


using std::cout;
using std::endl;
using std::ostream;


// constructor from default datas
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them

ParticleOnDiskLaplacianDeltaHamiltonian::ParticleOnDiskLaplacianDeltaHamiltonian(ParticleOnDisk* particles, int nbrParticles,
										 AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{
  this->Particles = particles;
  this->MaxMomentum = this->Particles->GetMaxLz();
  this->NbrLzValue = this->MaxMomentum + 1;
  this->NbrParticles = nbrParticles;
  this->FastMultiplicationFlag = false;
  this->Architecture = architecture;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  if (TmpMemory < 1024)
	    cout  << "fast = " <<  TmpMemory << "b ";
	  else
	    if (TmpMemory < (1 << 20))
	      cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	    else
	  if (TmpMemory < (1 << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    cout  << "fast = " << (TmpMemory >> 30) << "Gb ";
	  if (memory > 0)
	    {
	      this->EnableFastMultiplication();
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);
  this->HamiltonianShift = 0.0;
}

// destructor
//

ParticleOnDiskLaplacianDeltaHamiltonian::~ParticleOnDiskLaplacianDeltaHamiltonian() 
{
  delete[] this->InteractionFactors;
  delete[] this->M1Value;
  delete[] this->M2Value;
  delete[] this->M3Value;
  delete[] this->M4Value;
  if (this->FastMultiplicationFlag == true)
    {
      int ReducedDim = this->Particles->GetHilbertSpaceDimension() / this->FastMultiplicationStep;
      if ((ReducedDim * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
	++ReducedDim;
      for (int i = 0; i < ReducedDim; ++i)
	{
	  delete[] this->InteractionPerComponentIndex[i];
	  delete[] this->InteractionPerComponentCoefficient[i];
	}
      delete[] this->InteractionPerComponentIndex;
      delete[] this->InteractionPerComponentCoefficient;
      delete[] this->NbrInteractionPerComponent;
    }
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnDiskLaplacianDeltaHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnDisk*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void ParticleOnDiskLaplacianDeltaHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}
  
// evaluate all interaction factors
//   

void ParticleOnDiskLaplacianDeltaHamiltonian::EvaluateInteractionFactors()
{
  int Pos = 0;
  int m4;
  double* TmpCoefficient = new double [this->NbrLzValue * this->NbrLzValue * this->NbrLzValue];
  double MaxCoefficient = 0.0;

  if (this->Particles->GetParticleStatistic() == ParticleOnDisk::FermionicStatistic)
    {
     for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if ((m4 >= 0) && (m3 > m4))
		{
		  TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					 + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					 - this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					 - this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		  if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		    MaxCoefficient = fabs(TmpCoefficient[Pos]);
		  ++Pos;
		}
	    }
      this->NbrInteractionFactors = 0;
      this->M1Value = new int [Pos];
      this->M2Value = new int [Pos];
      this->M3Value = new int [Pos];
      this->M4Value = new int [Pos];
      this->InteractionFactors = new double [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;
      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 < m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if ((m4 >= 0) && (m3 > m4))
		{
		  if  (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
		      this->M1Value[this->NbrInteractionFactors] = m1;
		      this->M2Value[this->NbrInteractionFactors] = m2;
		      this->M3Value[this->NbrInteractionFactors] = m3;
		      this->M4Value[this->NbrInteractionFactors] = m4;
		      ++this->NbrInteractionFactors;
		    }
		  ++Pos;
		}
	    }
    }
  else
    {
      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if ((m4 >= 0) && (m3 > m4))
		{
		  if (m1 != m2)
		    {
		      TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					     + this->EvaluateInteractionCoefficient(m2, m1, m4, m3)
					     + this->EvaluateInteractionCoefficient(m1, m2, m4, m3)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		    }
		  else
		    TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					   + this->EvaluateInteractionCoefficient(m1, m2, m4, m3));
		  if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		    MaxCoefficient = fabs(TmpCoefficient[Pos]);
		  ++Pos;
		}
	      else
		if (m3 == m4)
		  {
		    if (m1 != m2)
		      TmpCoefficient[Pos] = (this->EvaluateInteractionCoefficient(m1, m2, m3, m4)
					     + this->EvaluateInteractionCoefficient(m2, m1, m3, m4));
		    else
		      TmpCoefficient[Pos] = this->EvaluateInteractionCoefficient(m1, m2, m3, m4);
		    if (MaxCoefficient < fabs(TmpCoefficient[Pos]))
		      MaxCoefficient = fabs(TmpCoefficient[Pos]);
		    ++Pos;
		  }
	    }
      this->NbrInteractionFactors = 0;
      this->M1Value = new int [Pos];
      this->M2Value = new int [Pos];
      this->M3Value = new int [Pos];
      this->M4Value = new int [Pos];
      this->InteractionFactors = new double [Pos];
      cout << "nbr interaction = " << Pos << endl;
      Pos = 0;
      MaxCoefficient *= MACHINE_PRECISION;
      for (int m1 = 0; m1 <= this->MaxMomentum; ++m1)
	for (int m2 = 0; m2 <= m1; ++m2)
	  for (int m3 = 0; m3 <= this->MaxMomentum; ++m3)
	    {
	      m4 = m1 + m2 - m3;
	      if ((m4 >= 0) && (m3 > m4))
		{
		  if (fabs(TmpCoefficient[Pos]) > MaxCoefficient)
		    {
		      this->InteractionFactors[this->NbrInteractionFactors] = TmpCoefficient[Pos];
		      this->M1Value[this->NbrInteractionFactors] = m1;
		      this->M2Value[this->NbrInteractionFactors] = m2;
		      this->M3Value[this->NbrInteractionFactors] = m3;
		      this->M4Value[this->NbrInteractionFactors] = m4;
		      ++this->NbrInteractionFactors;
		    }
		  ++Pos;
		}
	    }
    }
  cout << "nbr interaction = " << this->NbrInteractionFactors << endl;
  cout << "====================================" << endl;
  delete[] TmpCoefficient;
}

// evaluate the numerical coefficient  in front of the a+_m1 a+_m2 a_m3 a_m4 coupling term
//
// m1 = first index
// m2 = second index
// m3 = third index
// m4 = fourth index
// return value = numerical coefficient

double ParticleOnDiskLaplacianDeltaHamiltonian::EvaluateInteractionCoefficient(int m1, int m2, int m3, int m4)
{
  if ((m1 == m2) || (m3 == m4))
    return 0.0;
  FactorialCoefficient Coef;
  Coef.SetToOne();
  if (m2 > 1)
    {
      Coef.PartialFactorialMultiply(m1 + 1, m1 + m2 - 1);
      Coef.FactorialDivide(m2);
    }
  else
    {
      if (m2 == 0)
	Coef /= m1;	
    }
  if (m4 > 1)
    {
      Coef.PartialFactorialMultiply(m3 + 1, m3 + m4 - 1);
      Coef.FactorialDivide(m4);
    }
  else
    {
      if (m4 == 0)
	Coef /= m3;	
    }
  Coef.Power2Divide(2 * (m1 + m2));
  return (sqrt(Coef.GetNumericalValue()) * ((double) ((m2 - m1) * (m3 - m4))) / M_PI);
}

// Output Stream overload
//
// Str = reference on output stream
// H = Hamiltonian to print
// return value = reference on output stream

ostream& operator << (ostream& Str, ParticleOnDiskLaplacianDeltaHamiltonian& H) 
{
  RealVector TmpV2 (H.Particles->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Particles->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Particles->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Particles->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  for (int i = 0; i < H.Particles->GetHilbertSpaceDimension(); i++)
    {
      for (int j = 0; j < H.Particles->GetHilbertSpaceDimension(); j++)
	{
	  Str << TmpV[j][i] << "    ";
	}
      Str << endl;
    }
  return Str;
}

// Mathematica Output Stream overload
//
// Str = reference on Mathematica output stream
// H = Hamiltonian to print
// return value = reference on output stream

MathematicaOutput& operator << (MathematicaOutput& Str, ParticleOnDiskLaplacianDeltaHamiltonian& H) 
{
  RealVector TmpV2 (H.Particles->GetHilbertSpaceDimension(), true);
  RealVector* TmpV = new RealVector [H.Particles->GetHilbertSpaceDimension()];
  for (int i = 0; i < H.Particles->GetHilbertSpaceDimension(); i++)
    {
      TmpV[i] = RealVector(H.Particles->GetHilbertSpaceDimension());
      if (i > 0)
	TmpV2[i - 1] = 0.0;
      TmpV2[i] = 1.0;
      H.LowLevelMultiply (TmpV2, TmpV[i]);
    }
  Str << "{";
  for (int i = 0; i < (H.Particles->GetHilbertSpaceDimension() - 1); i++)
    {
      Str << "{";
      for (int j = 0; j < (H.Particles->GetHilbertSpaceDimension() - 1); j++)
	{
	  Str << TmpV[j][i] << ",";
	}
      Str << TmpV[H.Particles->GetHilbertSpaceDimension() - 1][i];
      Str << "},";
    }
  Str << "{";
  for (int j = 0; j < (H.Particles->GetHilbertSpaceDimension() - 1); j++)
    {
      Str << TmpV[j][H.Particles->GetHilbertSpaceDimension() - 1] << ",";
    }
  Str << TmpV[H.Particles->GetHilbertSpaceDimension() - 1][H.Particles->GetHilbertSpaceDimension() - 1];
  Str << "}}";
  return Str;
}

