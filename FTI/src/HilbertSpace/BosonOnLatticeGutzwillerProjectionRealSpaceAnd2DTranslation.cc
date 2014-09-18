////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of bosons hardcore on lattice                      //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Antoine Sterdyniak                    //
//                                                                            //
//                        last modification : 10/09/2014                      //
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
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "GeneralTools/StringTools.h"

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

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation ()
{
  this->NbrBosons = 0;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->MaxMomentum = 0;
  this->NbrSite = 0;
  this->MomentumModulo = 0;
  this->XMomentum = 0; 
  this->YMomentum = 0;
  this->StateShift = 2 * (this->MaxMomentum / this->MomentumModulo);
  this->MomentumIncrement = (this->NbrBosons * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = ((unsigned long) 1);
  this->MaximumSignLookUp = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->StateMaxMomentum = 0;  
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrBosons = number of fermions
// nbrSite = number of sites
// xMomentum = momentum sector in the x direction
// maxXMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction 
// memory = amount of memory granted for precalculations

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation (int nbrBosons, int nbrSite, int xMomentum, int  maxXMomentum,  int yMomentum, int maxYMomentum, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MaxXMomentum = maxXMomentum;
  this->MomentumModulo = this->MaxXMomentum;

  this->StateXShift = 1;
  this->MomentumIncrement = (this->NbrBosons * this->StateShift / 2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->XMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = this->NbrSite / this->MaxXMomentum;
  this->ComplementaryStateXShift = this->MaxMomentum - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;
//   cout << "this->MaxXMomentum=" << this->MaxXMomentum << endl;
//   cout << "this->XMomentum=" << this->XMomentum << endl;
//   cout << "this->StateXShift=" << this->StateXShift << endl;
//   cout << "this->ComplementaryStateXShift=" << this->ComplementaryStateXShift << endl;
//   cout << "this->XMomentumMask=" << hex << this->XMomentumMask << dec << endl;

  this->MaxYMomentum = maxYMomentum;
  this->YMomentum = yMomentum % this->MaxYMomentum;
  this->NbrYMomentumBlocks = this->MaxXMomentum;
  this->StateYShift = (this->NbrSite / (this->MaxXMomentum * this->MaxYMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (0x1ul << this->StateYShift) - 0x1ul;
  this->YMomentumBlockMask = (0x1ul << this->YMomentumBlockSize) - 0x1ul;  
//   cout << "this->MaxYMomentum=" << this->MaxYMomentum << endl;
//   cout << "this->NbrYMomentumBlocks=" << this->NbrYMomentumBlocks << endl;
//   cout << "this->StateYShift=" << this->StateYShift << endl;
//   cout << "this->YMomentumBlockSize=" << this->YMomentumBlockSize << endl;
//   cout << "this->ComplementaryStateYShift=" << this->ComplementaryStateYShift << endl;
//   cout << "this->YMomentumMask=" << hex << this->YMomentumMask << dec << endl;
//   cout << "this->YMomentumBlockMask=" << hex << this->YMomentumBlockMask << dec << endl;
	
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons);
  cout << "intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->LargeHilbertSpaceDimension  = this->GenerateStates();
      this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension > 0l)
	{
	  this->StateMaxMomentum = new int [this->LargeHilbertSpaceDimension];  
	  int CurrentMaxMomentum = this->MaxMomentum;
	  while (((this->StateDescription[0] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
	    --CurrentMaxMomentum;
	  this->StateMaxMomentum[0] = CurrentMaxMomentum;
	  for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
	    {
	      while (((this->StateDescription[i] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
		--CurrentMaxMomentum;
	      this->StateMaxMomentum[i] = CurrentMaxMomentum;
	    }
	  this->GenerateLookUpTable(memory);
	  
#ifdef __DEBUG__
	  long UsedMemory = 0;
	  UsedMemory += (long) this->HilbertSpaceDimension * (sizeof(unsigned long) + sizeof(int));
	  cout << "memory requested for Hilbert space = ";
	  if (UsedMemory >= 1024)
	    if (UsedMemory >= 1048576)
	      cout << (UsedMemory >> 20) << "Mo" << endl;
	    else
	      cout << (UsedMemory >> 10) << "ko" <<  endl;
	  else
	    cout << UsedMemory << endl;
	  UsedMemory = this->NbrMomentum * sizeof(int);
	  UsedMemory += this->NbrMomentum * this->LookUpTableMemorySize * sizeof(int);
	  cout << "memory requested for lookup table = ";
	  if (UsedMemory >= 1024)
	    if (UsedMemory >= 1048576)
	      cout << (UsedMemory >> 20) << "Mo" << endl;
	    else
	      cout << (UsedMemory >> 10) << "ko" <<  endl;
	  else
	    cout << UsedMemory << endl;
#endif
	}
    }
 this->SignLookUpTable = 0;
 this->NbrParticleLookUpTable = 0;
}

// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation(const BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation& bosons)
{
  this->NbrBosons = bosons.NbrBosons;  
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->NbrSite = bosons.NbrSite;

  this->MaxXMomentum = bosons.MaxXMomentum;
  this->XMomentum = bosons.XMomentum;
  this->StateXShift = bosons.StateXShift;
  this->ComplementaryStateXShift = bosons.ComplementaryStateXShift;
  this->XMomentumMask = bosons.XMomentumMask;
  this->MaxYMomentum = bosons.MaxYMomentum;
  this->YMomentum = bosons.YMomentum;
  this->NbrYMomentumBlocks = bosons.NbrYMomentumBlocks;
  this->StateYShift = bosons.StateYShift;
  this->YMomentumBlockSize = bosons.YMomentumBlockSize;
  this->ComplementaryStateYShift = bosons.ComplementaryStateYShift;
  this->YMomentumMask = bosons.YMomentumMask;
  this->YMomentumBlockMask = bosons.YMomentumBlockMask;  

  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrMomentum = bosons.NbrMomentum;
  this->MomentumModulo = bosons.MomentumModulo;
  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->MomentumMask = bosons.MomentumMask;

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;

  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;

  this->NbrParticleLookUpTable = bosons.NbrParticleLookUpTable;

  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

  this->Flag = bosons.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::~BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation& BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::operator = (const BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateMaxMomentum;

      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      delete[] this->NbrParticleLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  this->NbrBosons = bosons.NbrBosons;  
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->NbrSite = bosons.NbrSite;

  this->MaxXMomentum = bosons.MaxXMomentum;
  this->XMomentum = bosons.XMomentum;
  this->StateXShift = bosons.StateXShift;
  this->ComplementaryStateXShift = bosons.ComplementaryStateXShift;
  this->XMomentumMask = bosons.XMomentumMask;
  this->MaxYMomentum = bosons.MaxYMomentum;
  this->YMomentum = bosons.YMomentum;
  this->NbrYMomentumBlocks = bosons.NbrYMomentumBlocks;
  this->StateYShift = bosons.StateYShift;
  this->YMomentumBlockSize = bosons.YMomentumBlockSize;
  this->ComplementaryStateYShift = bosons.ComplementaryStateYShift;
  this->YMomentumMask = bosons.YMomentumMask;
  this->YMomentumBlockMask = bosons.YMomentumBlockMask;  

  this->MaxMomentum = bosons.MaxMomentum;
  this->NbrMomentum = bosons.NbrMomentum;
  this->MomentumModulo = bosons.MomentumModulo;
  this->MomentumIncrement = bosons.MomentumIncrement;
  this->StateShift = bosons.StateShift;
  this->ComplementaryStateShift = bosons.ComplementaryStateShift;
  this->MomentumMask = bosons.MomentumMask;

  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->StateDescription = bosons.StateDescription;
  this->StateMaxMomentum = bosons.StateMaxMomentum;

  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;

  this->NbrParticleLookUpTable = bosons.NbrParticleLookUpTable;

  this->RescalingFactors = bosons.RescalingFactors;
  this->NbrStateInOrbit = bosons.NbrStateInOrbit;

  this->Flag = bosons.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::Clone()
{
  return new BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation(*this);
}


// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  this->RawGenerateStates(this->NbrBosons, this->NbrSite - 1, 0l);
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  int NbrTranslationY;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, NbrTranslationY) == this->StateDescription[i]))
	{
	//cout <<  this->StateDescription[i] <<endl;
	  if (this->TestMomentumConstraint(this->StateDescription[i]) == true)
	    {
	      ++TmpLargeHilbertSpaceDimension;
	    }
	  else
	    {
	      this->StateDescription[i] = 0x0ul;
	    }
	}
      else
	{
	  this->StateDescription[i] = 0x0ul;
	}
    }
  //cout << "new dim = " << TmpLargeHilbertSpaceDimension << endl;
  unsigned long* TmpStateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];  
  this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
  TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (this->StateDescription[i] != 0x0ul)
	{
	  TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
	  this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] = this->FindOrbitSize(this->StateDescription[i]);
	  ++TmpLargeHilbertSpaceDimension;
	}	
    }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;
  return TmpLargeHilbertSpaceDimension;
}

// apply a^+_m_sigma a_n_sigma operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator including the orbital and the spin index
// n = index of the annihilation operator including the orbital and the spin index
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::AdA (int index, int m, int n, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long State = this->StateDescription[index];
  if ((State & (0x1ul << n)) == 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  State &= ~(0x1ul << n);
  if ((State & (0x1ul << m))!= 0)
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  State |= (0x1ul << m);
  return this->SymmetrizeAdAdResult(State, coefficient, nbrTranslationX, nbrTranslationY);
}
  
// apply a^+_m1_sigma a^+_m2_sigma operator to the state produced using AuAu method (without destroying it)
//
// m1 = first index for creation operator
// m2 = second index for creation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslationX = reference on the number of translations to applied in the x direction to the resulting state to obtain the return orbit describing state
// nbrTranslationY = reference on the number of translations in the y direction to obtain the canonical form of the resulting state
// return value = index of the destination state 

int BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::AdAd (int m1, int m2, double& coefficient, int& nbrTranslationX, int& nbrTranslationY)
{
  unsigned long TmpState = this->ProdATemporaryState;
  if (((TmpState & (0x1ul << m1)) != 0) || ((TmpState & (0x1ul << m2)) != 0) || (m1 == m2))
    {
      coefficient = 0.0;
      return this->HilbertSpaceDimension;
    }
  coefficient = 1.0;
  TmpState |= (0x1ul << m1);
  return this->SymmetrizeAdAdResult(TmpState, coefficient, nbrTranslationX, nbrTranslationY);
}


void BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation::GetCompositeFermionWavefunction(ComplexVector & trialState, ComplexMatrix & jastrowEigenVecs,ComplexMatrix & cFEigenVecs)
{

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

 NbrTranslation = 0;
 unsigned long  TmpStateDescription  = this->StateDescription[i];
 unsigned long  TmpStateDescription2  = this->StateDescription[i];
 for (int m = 0; (m < this->MaxYMomentum) && (TmpStateDescription !=  this->StateDescription[i]) ; ++m)
    {
      TmpStateDescription2 = TmpStateDescription;
      for (int n = 0; (n < this->MaxXMomentum) && (TmpStateDescription2 != TmpStateDescription) ; ++n)
	{
          int  TmpMaxMomentum = this->NbrSite;
          while ((( TmpStateDescription2 >> TmpMaxMomentum) & 0x1ul) == 0x0ul)
	       --TmpMaxMomentum;
          this->ConvertToMonomial(TmpStateDescription2,TemporaryState);

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
        trialState[i] +=  ExponentialFactors[n][m] * SlaterCF.Determinant() * SlaterJastrow.Determinant();
        ++NbrTranslation;
        this->ApplySingleXTranslation(TmpStateDescription2);      
	}
      this->ApplySingleYTranslation(TmpStateDescription);      
    }

   
  trialState[i] /= sqrt(this->NbrStateInOrbit[i]);
  if(NbrTranslation != this->NbrStateInOrbit[i])
   {
	cout <<"Wrong Number of Translation for state !"  << i <<endl;
  }

}
   for (int i = 0; i < this->MaxXMomentum; ++i)
    { 
      delete [] ExponentialFactors[i];
    }
   delete []  ExponentialFactors;
}



