////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of U(1) wave function obtained form symmetrization of a   //
//                        SU(K) wave function on sphere                       //
//                                                                            //
//                        last modification : 17/03/2008                      //
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
#include "MathTools/BinomialCoefficients.h"
#include "Tools/FQHEWaveFunction/FQHESphereSymmetrizedSUKToU1WaveFunction.h"
#include "Tools/FQHEWaveFunction/HalperinOnSphereWaveFunction.h"
#include "Vector/RealVector.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;


// group maximum size in bits
#ifdef __64_BITS__
#define GROUP_MAXSIZE 12
#else
#define GROUP_MAXSIZE 6
#endif

// constructor
//
// nbrParticles = number of particles
// kValue = number of particle per cluster
// sUKWavefunction = pointer to the base SU(K) wave function
// fermionFlag = true if the final state should be a fermionic state
 
FQHESphereSymmetrizedSUKToU1WaveFunction::FQHESphereSymmetrizedSUKToU1WaveFunction(int nbrParticles, int kValue, Abstract1DComplexFunctionOnSphere* sUKWavefunction, bool fermionFlag)
{
  this->NbrParticles = nbrParticles;
  this->KValue = kValue;
  this->NbrParticlesPerColor = this->NbrParticles / this->KValue;
  this->SUKWaveFunction = (Abstract1DComplexFunctionOnSphere*) (sUKWavefunction->Clone());
  this->EvaluatePermutations();
  this->FermionFlag = fermionFlag;
  this->Flag.Initialize();
}

// copy constructor
//
// function = reference on the wave function to copy

FQHESphereSymmetrizedSUKToU1WaveFunction::FQHESphereSymmetrizedSUKToU1WaveFunction(const FQHESphereSymmetrizedSUKToU1WaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->NbrParticlesPerColor = function.NbrParticlesPerColor;
  this->KValue = function.KValue;
  this->Permutations = function.Permutations;
  this->NbrPermutations = function.NbrPermutations;
  this->SUKWaveFunction = (Abstract1DComplexFunctionOnSphere*) (function.SUKWaveFunction->Clone());
  this->FermionFlag = function.FermionFlag;
  this->Flag = function.Flag;
}

// destructor
//

FQHESphereSymmetrizedSUKToU1WaveFunction::~FQHESphereSymmetrizedSUKToU1WaveFunction()
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      for (unsigned long i = 0; i < this->NbrPermutations; ++i)
	delete[] this->Permutations[i];
      delete[] this->Permutations;
    }
  delete this->SUKWaveFunction;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* FQHESphereSymmetrizedSUKToU1WaveFunction::Clone ()
{
  return new FQHESphereSymmetrizedSUKToU1WaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex FQHESphereSymmetrizedSUKToU1WaveFunction::operator ()(RealVector& x)
{  
  ComplexVector TmpUV (this->NbrParticles * 2);
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      TmpUV[2 * i].Re = cos(0.5 * x[i << 1]);
      TmpUV[2 * i].Im = TmpUV[2 * i].Re;
      TmpUV[2 * i].Re *= cos(0.5 * x[1 + (i << 1)]);
      TmpUV[2 * i].Im *= sin(0.5 * x[1 + (i << 1)]);
      TmpUV[(2 * i) + 1].Re = sin(0.5 * x[i << 1]);
      TmpUV[(2 * i) + 1].Im = TmpUV[(2 * i) + 1].Re;
      TmpUV[(2 * i) + 1].Re *= cos(0.5 * x[1 + (i << 1)]);
      TmpUV[(2 * i) + 1].Im *= -sin(0.5 * x[1 + (i << 1)]);
    }
  return this->CalculateFromSpinorVariables(TmpUV);
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)

Complex FQHESphereSymmetrizedSUKToU1WaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{  
//  HalperinOnSphereWaveFunction* BaseFunction = new HalperinOnSphereWaveFunction (NbrParticles >> 1, NbrParticles >> 1, 2, 2, 0);
  Complex TotalValue = 0.0;
  ComplexVector TmpUV (this->NbrParticles * 2);
  for (unsigned long i = 0ul; i < this->NbrPermutations; ++i)
    {
      int TotalIndex = 0;
      for (int k = 0; k < this->NbrParticlesPerColor; ++k)
	{
	  unsigned long Tmp = this->Permutations[i][k];
	  for (int j = 0; j < this->KValue; ++j)
	    {
	      unsigned long TmpIndex = ((Tmp >> (j * GROUP_MAXSIZE)) & 0x3ful) << 1;
	      int TmpIndex2 = TotalIndex + (j * this->NbrParticlesPerColor * 2);
	      TmpUV.Re(TmpIndex2) = uv.Re(TmpIndex);
	      TmpUV.Im(TmpIndex2) = uv.Im(TmpIndex);
	      TmpIndex2 += 1;
	      TmpIndex += 1ul;
	      TmpUV.Re(TmpIndex2) = uv.Re(TmpIndex);
	      TmpUV.Im(TmpIndex2) = uv.Im(TmpIndex);
	    }
	  TotalIndex += 2;
	}
//       int TotalIndex = 0;
//       unsigned long Tmp = this->ColorPermutations[i];
//       for (int k = 0; k < this->NbrParticles; ++k)
// 	{
// 	  unsigned long TmpIndex = ((Tmp >> (k * 4)) & 0xful) << 1;
// //	  cout << (TotalIndex >> 1) << "->" << (TmpIndex >> 1) << ", ";
// 	  TmpUV.Re(TmpIndex) = uv.Re(TotalIndex);
// 	  TmpUV.Im(TmpIndex) = uv.Im(TotalIndex);
// 	  ++TotalIndex;
// 	  ++TmpIndex;
// 	  TmpUV.Re(TmpIndex) = uv.Re(TotalIndex);
// 	  TmpUV.Im(TmpIndex) = uv.Im(TotalIndex);
// 	  ++TotalIndex;
// 	}
// //      cout << " = " << this->SUKWaveFunction->CalculateFromSpinorVariables(TmpUV) << " " << BaseFunction->CalculateFromSpinorVariables(TmpUV) << endl;
      TotalValue += this->SUKWaveFunction->CalculateFromSpinorVariables(TmpUV);
    }
  if (this->FermionFlag == true)
    {
      Complex TmpU;
      Complex TmpV;
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrParticles; ++i)
	{
	  TmpU = TmpUV[2 * i];
	  TmpV = TmpUV[2 * i + 1];
	  for (int j = i + 1; j < this->NbrParticles; ++j)
	    WaveFunction *=  ((TmpU * TmpUV[2 * j + 1]) - (TmpV * TmpUV[2 * j]));
	}
      TotalValue *= WaveFunction;
    }
  return TotalValue;
}  


// evaluate all permutations requested to symmetrize the SU(K) state
//

void FQHESphereSymmetrizedSUKToU1WaveFunction::EvaluatePermutations()
{
  unsigned long Fact = 2;
  for (unsigned long i = 3; i <= ((unsigned long) this->NbrParticles); ++i)
    Fact *= i;
  unsigned long** Perm = new unsigned long* [Fact];

  Perm[0] = new unsigned long[this->NbrParticlesPerColor];
  unsigned long* TmpPerm = Perm[0];
  int Shift = 0;   
  for (int k = 0; k < this->NbrParticlesPerColor; ++k) 
    {
      TmpPerm[k] = (unsigned long) 0;
      for (int i = 0; i < this->KValue; ++i)    
        TmpPerm[k] |= ((unsigned long) (i + Shift)) << (i * GROUP_MAXSIZE);
      Shift += this->KValue; 
    }
  for (unsigned long i = 1; i < Fact; ++i)
    {
      Perm[i] = new unsigned long[this->NbrParticlesPerColor];
      for (int k = 0; k < this->NbrParticlesPerColor; ++k)
        Perm[i][k] = TmpPerm[k];
    }

  int GroupId = 0;
  int PosInGroup = 0;
  for (int Pos = 0; Pos < (this->NbrParticles - 2); ++Pos)
    {
      unsigned long Step = (unsigned long) 1;
      unsigned long Lim = (unsigned long)(this->NbrParticles - Pos - 1);
      for (unsigned long i = 2; i <= Lim; ++i) 
	Step *= i;
      unsigned long j = 0;
      unsigned long Lim2 = Fact / (Step * (Lim + 1l));
      for (unsigned long k = 0; k < Lim2; ++k)
	{
	  j += Step;
	  Lim = j + Step;
	  int GroupId2 = GroupId;
	  int PosInGroup2 = PosInGroup + 1;
	  if (PosInGroup2 == this->KValue)
	    {
	      PosInGroup2 = 0;
	      ++GroupId2;
	    }
	  for (int i = Pos + 1; i < this->NbrParticles; ++i)    
	    {
	      unsigned long Mask1 = ((unsigned long) 0x3f) << (PosInGroup * GROUP_MAXSIZE);
	      unsigned long Mask2 = ((unsigned long) 0x3f) << (PosInGroup2 * GROUP_MAXSIZE);
	      unsigned long NegMask1 = ~Mask1;
	      unsigned long NegMask2 = ~Mask2;
	      int Shift = (PosInGroup2 - PosInGroup) * GROUP_MAXSIZE;   
	      if (Shift > 0)
		{
		  for (; j < Lim; ++j)
		    {
		      unsigned long Tmp1 = (Perm[j][GroupId] & Mask1) << Shift;
		      unsigned long Tmp2 = (Perm[j][GroupId2] & Mask2) >> Shift;
		      Perm[j][GroupId] &= NegMask1;
		      Perm[j][GroupId] |= Tmp2;
		      Perm[j][GroupId2] &= NegMask2;
		      Perm[j][GroupId2] |= Tmp1;
		    }
		}
	      else
		{
		  Shift *= -1;
		  for (; j < Lim; ++j)
		    {
		      unsigned long Tmp1 = (Perm[j][GroupId] & Mask1) >> Shift;
		      unsigned long Tmp2 = (Perm[j][GroupId2] & Mask2) << Shift;
		      Perm[j][GroupId] &= NegMask1;
		      Perm[j][GroupId] |= Tmp2;
		      Perm[j][GroupId2] &= NegMask2;
		      Perm[j][GroupId2] |= Tmp1;
		    }
		}
	      ++PosInGroup2;
	      if (PosInGroup2 == this->KValue)
		{
		  PosInGroup2 = 0;
		  ++GroupId2;
		}
	      Lim += Step;
	    } 
	}
      ++PosInGroup;
      if (PosInGroup == this->KValue)
	{
	  PosInGroup = 0;
	  ++GroupId;
	}
    }
  unsigned long Mask1 = ((unsigned long) 0x3f) << (PosInGroup * GROUP_MAXSIZE);
  unsigned long Mask2 = ((unsigned long) 0x3f) << ((PosInGroup + 1) * GROUP_MAXSIZE);
  unsigned long NegMask = ~(Mask1 | Mask2);
  unsigned long Tmp;
  this->NbrPermutations = Fact;
  for (unsigned long i = 2; i < ((unsigned long) this->NbrParticlesPerColor); ++i)
    this->NbrPermutations /= i;
  this->Permutations = new unsigned long* [this->NbrPermutations];
  this->NbrPermutations = 0;
  int j;
  bool TmpFlag;
  int ReducedNbrParticlesPerColor = this-> NbrParticlesPerColor - 1;
  for (unsigned long i = 0; i < Fact; ++i)
    {
      j = 0;
      TmpFlag = false;
      while ((j < ReducedNbrParticlesPerColor) && (TmpFlag == false))
	{
	  if ((Perm[i][j] & ((unsigned long) 0x3f)) > (Perm[i][j + 1] & ((unsigned long) 0x3f))) 
	    TmpFlag = true;
	  ++j;
	}
      if (TmpFlag == false)
	{
	  this->Permutations[this->NbrPermutations] = Perm[i];
	  ++this->NbrPermutations;
	}
      else
	delete[] Perm[i];
      ++i;
      Tmp = ((Perm[i][GroupId] & Mask1) << GROUP_MAXSIZE) | ((Perm[i][GroupId] & Mask2) >> GROUP_MAXSIZE);
      Perm[i][GroupId] &= NegMask;
      Perm[i][GroupId] |= Tmp;
      j = 0;
      TmpFlag = false;
      while ((j < ReducedNbrParticlesPerColor) && (TmpFlag == false))
	{
	  if ((Perm[i][j] & ((unsigned long) 0x3f)) > (Perm[i][j + 1] & ((unsigned long) 0x3f))) 
	    TmpFlag = true;
	  ++j;
	}
      if (TmpFlag == false)
	{
	  this->Permutations[this->NbrPermutations] = Perm[i];
	  ++this->NbrPermutations;
	}
      else
	delete[]  Perm[i];
    }
  delete[] Perm;

  return;
}

