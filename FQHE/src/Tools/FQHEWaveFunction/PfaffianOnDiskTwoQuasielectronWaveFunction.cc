////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of Pfaffian wave function with two quasiholes on disk          //
//                                                                            //
//                        last modification : 23/10/2008                      //
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
#include "Tools/FQHEWaveFunction/PfaffianOnDiskTwoQuasielectronWaveFunction.h"
#include "Vector/RealVector.h"
#include "Matrix/ComplexSkewSymmetricMatrix.h"
#include "GeneralTools/Endian.h"

#include <iostream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// group maximum size in bits
#ifdef __64_BITS__
#define GROUP_MAXSIZE 12
#else
#define GROUP_MAXSIZE 6
#endif


// constructor
//
// nbrParticles = number of particles
// zElectron1 = position of the first quasielectron
// zElectron2 = position of the second quasielectron (spherical coordinates, theta angle)
// fermions = flag indicating whether to calculate bosonic or fermionic pfaffian

PfaffianOnDiskTwoQuasielectronWaveFunction::PfaffianOnDiskTwoQuasielectronWaveFunction(int nbrParticles,  Complex zElectron1, Complex zElectron2, bool fermions)
{
  this->NbrParticles = nbrParticles;
  this->ZElectron1 = zElectron1;
  this->ZElectron2 = zElectron2;
  this->GaussianWeight = exp (-0.125 * (SqrNorm(this->ZElectron1) + SqrNorm(this->ZElectron2)));

  this->TmpPfaffian = new Complex* [this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    this->TmpPfaffian[i] = new Complex [this->NbrParticles];
  this->TmpIndexArray = new int [this->NbrParticles];
  this->TmpGaussianWeights = new Complex [2 * this->NbrParticles];
  this->FermionFlag = fermions;

  this->EvaluatePermutations();
  this->Flag.Initialize();
}


// constructor from data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// zElectron1 = position of the first quasielectron
// zElectron2 = position of the second quasielectron (spherical coordinates, theta angle)
// fermions = flag indicating whether to calculate bosonic or fermionic pfaffian

PfaffianOnDiskTwoQuasielectronWaveFunction::PfaffianOnDiskTwoQuasielectronWaveFunction(char* filename, Complex zElectron1, Complex zElectron2, bool fermions)
{
  this->NbrParticles = 0;
  this->ReadPermutations(filename);
  this->ZElectron1 = zElectron1;
  this->ZElectron2 = zElectron2;
  this->GaussianWeight = exp (-0.125 * (SqrNorm(this->ZElectron1) + SqrNorm(this->ZElectron2)));

  this->TmpPfaffian = new Complex* [this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    this->TmpPfaffian[i] = new Complex [this->NbrParticles];
  this->TmpIndexArray = new int [this->NbrParticles];
  this->TmpGaussianWeights = new Complex [2 * this->NbrParticles];
  this->FermionFlag = fermions;

  this->Flag.Initialize();
}


// copy constructor
//
// function = reference on the wave function to copy

PfaffianOnDiskTwoQuasielectronWaveFunction::PfaffianOnDiskTwoQuasielectronWaveFunction(const PfaffianOnDiskTwoQuasielectronWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->ZElectron1 = function.ZElectron1;
  this->ZElectron2 = function.ZElectron2;
  this->GaussianWeight = exp (-0.125 * (SqrNorm(this->ZElectron1) + SqrNorm(this->ZElectron2)));
  this->TmpPfaffian = new Complex* [this->NbrParticles];
  for (int i = 0; i < this->NbrParticles; ++i)
    this->TmpPfaffian[i] = new Complex [this->NbrParticles];
  this->TmpIndexArray = new int [this->NbrParticles];
  this->TmpGaussianWeights = new Complex [2 * this->NbrParticles];

  this->FermionFlag=function.FermionFlag;

  this->Permutations = function.Permutations;
  this->NbrPermutations = function.NbrPermutations;
  this->Flag = function.Flag;
}

// destructor
//

PfaffianOnDiskTwoQuasielectronWaveFunction::~PfaffianOnDiskTwoQuasielectronWaveFunction()
{
  for (int i = 0; i < this->NbrParticles; ++i)
    delete[] this->TmpPfaffian[i];
  delete[] this->TmpPfaffian;
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete[] this->Permutations;
    }
  delete[] this->TmpIndexArray;
  delete[] this->TmpGaussianWeights;
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* PfaffianOnDiskTwoQuasielectronWaveFunction::Clone ()
{
  return new PfaffianOnDiskTwoQuasielectronWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex PfaffianOnDiskTwoQuasielectronWaveFunction::operator ()(RealVector& x)
{
  Complex Tmp, T1, T2;
  Complex WaveFunction(1.0);
  Complex TmpZ;
  Complex TmpZElectron1;
  Complex TmpZElectron2;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      TmpZ.Re = x[i << 1];
      TmpZ.Im = x[1 + (i << 1)];
      this->TmpGaussianWeights[i << 1] = exp (0.25 * (Conj(this->ZElectron1) * TmpZ));
      this->TmpGaussianWeights[(i << 1) + 1] = exp (0.25 * (Conj(this->ZElectron1) * TmpZ));
      for (int j = i + 1; j < this->NbrParticles; ++j)
	{
	  Tmp.Re = TmpZ.Re - x[j << 1];
	  Tmp.Im = TmpZ.Im - x[1 + (j << 1)];
	  this->TmpPfaffian[i][j] = Tmp;
	  this->TmpPfaffian[j][i] = -Tmp;
	  WaveFunction *= Tmp;
	}
    }
  if (this->FermionFlag == false)
    WaveFunction = 1.0;

  int NbrParticlesPerColor = this->NbrParticles >> 1;
  Complex WaveFunctionSum = 0.0;
  for (unsigned long i = 0ul; i < this->NbrPermutations; ++i)
    {
      unsigned long TmpPerm = this->Permutations[i];
      unsigned long TmpPerm2 = TmpPerm >> (NbrParticlesPerColor << 2);
      TmpPerm &= (0x1ul << (NbrParticlesPerColor << 2)) - 0x1ul;
      Complex WaveFunctionPart1 = 0.0;
      for (int j = 0; j < NbrParticlesPerColor; ++ j)	
	{
	  Complex WaveFunctionPart12 = 0.0;	  
	  for (int k = 0; k < NbrParticlesPerColor; ++k)
	    TmpIndexArray[k] = (TmpPerm >> (k << 2)) & 0xful;
	  TmpPerm = (TmpPerm  >> 4) | ((TmpPerm & 0xful) << ((NbrParticlesPerColor - 1) << 2));

	  Complex* TmpArray = this->TmpPfaffian[TmpIndexArray[0]];
	  for (int k = 1; k < NbrParticlesPerColor; ++k)
	    {
	      Tmp = 1.0; 
	      for (int l = 1; l < k; ++l)
		Tmp *= TmpArray[TmpIndexArray[l]];	      
	      for (int l = k + 1; l < NbrParticlesPerColor; ++l)
		Tmp *= TmpArray[TmpIndexArray[l]];	      
	      WaveFunctionPart12 += Tmp;
	    }
	  for (int k = 1; k < NbrParticlesPerColor; ++k)
	    for (int l = k + 1; l < NbrParticlesPerColor; ++l)
	      {
		TmpZ = this->TmpPfaffian[TmpIndexArray[k]][TmpIndexArray[l]];
		WaveFunctionPart12 *= TmpZ;
		WaveFunctionPart12 *= TmpZ;
	      }
	  WaveFunctionPart12 *= this->TmpGaussianWeights[(TmpIndexArray[0] << 1)];
	  WaveFunctionPart1 += WaveFunctionPart12;
	}
      Complex WaveFunctionPart2 = 0.0;
      for (int j = 0; j < NbrParticlesPerColor; ++ j)
	{
	  Complex WaveFunctionPart22 = 0.0;	  
	  for (int k = 0; k < NbrParticlesPerColor; ++k)
	    TmpIndexArray[k] = (TmpPerm2 >> (k << 2)) & 0xful;
	  TmpPerm2 = (TmpPerm2  >> 4) | ((TmpPerm2 & 0xful) << ((NbrParticlesPerColor - 1) << 2));
	  Complex* TmpArray = this->TmpPfaffian[TmpIndexArray[0]];
	  for (int k = 1; k < NbrParticlesPerColor; ++k)
	    {
	      Tmp = 1.0; 
	      for (int l = 1; l < k; ++l)
		Tmp *= TmpArray[TmpIndexArray[l]];	      
	      for (int l = k + 1; l < NbrParticlesPerColor; ++l)
		Tmp *= TmpArray[TmpIndexArray[l]];	      
	      WaveFunctionPart22 += Tmp;
	    }
	  for (int k = 1; k < NbrParticlesPerColor; ++k)
	    for (int l = k + 1; l < NbrParticlesPerColor; ++l)
	      {
		TmpZ = this->TmpPfaffian[TmpIndexArray[k]][TmpIndexArray[l]];
		WaveFunctionPart22 *= TmpZ;
		WaveFunctionPart22 *= TmpZ;
	      }
	  WaveFunctionPart22 *= this->TmpGaussianWeights[(TmpIndexArray[0] << 1)];
	  WaveFunctionPart2 += WaveFunctionPart22;
	}      

      WaveFunctionSum += WaveFunctionPart1 * WaveFunctionPart2;
    }
  WaveFunction *= this->GaussianWeight;
  WaveFunction *= WaveFunctionSum;
  return WaveFunction;
}

// evaluate all permutations requested to symmetrize the state
//

void PfaffianOnDiskTwoQuasielectronWaveFunction::EvaluatePermutations()
{
  unsigned long Fact = 2;
  for (unsigned long i = 3; i <= ((unsigned long) this->NbrParticles); ++i)
    Fact *= i;

  int NbrParticlesPerColor = this->NbrParticles >> 1;
  unsigned long TmpNbrPermutation = Fact;  
  unsigned long FactNbrParticlesPerColor = 1;
  for (unsigned long i = 2; i <= ((unsigned long) NbrParticlesPerColor); ++i)
    FactNbrParticlesPerColor *= i;
  TmpNbrPermutation /= FactNbrParticlesPerColor;
  TmpNbrPermutation /= FactNbrParticlesPerColor;
  TmpNbrPermutation /= 2;
  
  this->Permutations = new unsigned long[TmpNbrPermutation];
  
  unsigned long TmpPerm =  0x0ul;
  unsigned long TmpPerm2 = 0x0ul;
  unsigned long TmpPerm3 = 0x0ul;
  unsigned long* TmpArrayPerm = new unsigned long [this->NbrParticles];
  for (int k = 0; k < this->NbrParticles; ++k) 
    {
      TmpPerm |= ((unsigned long) k) << (k << 2);
      TmpArrayPerm[k] = (unsigned long) k;
    }
  
  this->Permutations[0] = TmpPerm;
  TmpNbrPermutation = 1ul;
  Fact /= ((unsigned long) this->NbrParticles);
  Fact *= ((unsigned long) (this->NbrParticles - NbrParticlesPerColor + 1));
  for (unsigned long j = 1; j < Fact; ++j)
    {
      int Pos1 = this->NbrParticles - 1;
      while (TmpArrayPerm[Pos1 - 1] >= TmpArrayPerm[Pos1])
	--Pos1;
      --Pos1;
      int Pos2 = this->NbrParticles - 1;      
      while (TmpArrayPerm[Pos2] <= TmpArrayPerm[Pos1])
	--Pos2;
      unsigned long TmpIndex = TmpArrayPerm[Pos1];
      TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
      TmpArrayPerm[Pos2] = TmpIndex;
      Pos2 = this->NbrParticles - 1;   
      Pos1++;
      while (Pos1 < Pos2)
	{
	  TmpIndex = TmpArrayPerm[Pos1];
	  TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	  TmpArrayPerm[Pos2] = TmpIndex;
	  ++Pos1;
	  --Pos2;
	}
      bool Flag = true;
      int TmpIndex2 = 0;
      int k = 1;
      while ((k < NbrParticlesPerColor) && (Flag == true))
	{
	  if (TmpArrayPerm[TmpIndex2] > TmpArrayPerm[TmpIndex2 + 1])
	    Flag = false;
	  ++TmpIndex2;
	  ++k;
	}
      ++TmpIndex2;
      k = 1;
      while ((k < NbrParticlesPerColor) && (Flag == true))
	{
	  if (TmpArrayPerm[TmpIndex2] > TmpArrayPerm[TmpIndex2 + 1])
	    Flag = false;
	  ++TmpIndex2;
	  ++k;
	}
      if (Flag == true)
	{
	  TmpPerm2 = 0ul;
	  TmpPerm3 = 0ul;
	  for (int i = 0; i < NbrParticlesPerColor; ++i)
	    {
	      TmpPerm2 |= TmpArrayPerm[i] << (i << 2);
	      TmpPerm3 |= TmpArrayPerm[i + NbrParticlesPerColor] << (i << 2);
	    }
	  if (TmpPerm2 < TmpPerm3)
	    {
	      this->Permutations[TmpNbrPermutation] = TmpPerm2 | (TmpPerm3 << (NbrParticlesPerColor << 2));	      
	      ++TmpNbrPermutation;
	    }
	}
    }
  
  this->NbrPermutations = TmpNbrPermutation;

  return;
}


// get all permutations requested to symmetrize the SU(K) state from data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// return value = true if no error occured

bool PfaffianOnDiskTwoQuasielectronWaveFunction::ReadPermutations(char* filename)
{
  ifstream File;
  File.open(filename, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << filename << endl;
      return false;
    }
  ReadLittleEndian(File, this->NbrParticles);
  ReadLittleEndian(File, this->NbrPermutations);
  this->Permutations = new unsigned long[this->NbrPermutations];
  for (unsigned long i = 0; i < this->NbrPermutations; ++i)
    ReadLittleEndian(File, this->Permutations[i]);
  File.close();
  return true;
}

// write all permutations requested to symmetrize the SU(K) state to data file 
//
// filename = pointer to the file name that described the symmetrization procedure
// return value = true if no error occured

bool PfaffianOnDiskTwoQuasielectronWaveFunction::WritePermutations(char* filename)
{
  ofstream File;
  File.open(filename, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << filename << endl;
      return false;
    }
  WriteLittleEndian(File, this->NbrParticles);
  WriteLittleEndian(File, this->NbrPermutations);
  for (unsigned long i = 0; i < this->NbrPermutations; ++i)
    WriteLittleEndian(File, this->Permutations[i]);
  File.close();
  return true;
}

