////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on lattice with spin                   //
//       in real space with translation invariance in two directions          //
//        with a constraint on the mininum number of on-site singlet          //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                        last modification : 14/07/2016                      //
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
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets.h"
#include "HilbertSpace/FermionOnLatticeWithSpinRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpace.h"
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
#include "Architecture/ArchitectureOperation/FQHETorusParticleEntanglementSpectrumOperation.h"
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

FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets::FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets ()
{
  this->NbrFermions = 0;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalSpin = 0;
  this->MaxMomentum = 0;
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSite = 0;
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = 1;
  this->XMomentum = 0; 
  this->YMomentum = 0;
  this->NbrSitePerUnitCell = 0;
  this->StateShift = 2 * (this->MaxMomentum / this->MomentumModulo);
  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = ((unsigned long) 1);
  this->MaximumSignLookUp = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->StateHighestBit = 0;  
  this->LargeHilbertSpaceDimension = 0;
  this->RescalingFactors = 0;
  this->TargetSpace = this;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// xMomentum = momentum sector in the x direction
// maxXMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets::FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets (int nbrFermions, int minNbrSinglets, 
																  int nbrSite, int xMomentum, int maxXMomentum,
																  int yMomentum, int  maxYMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = false;
  this->TotalSpin = 0;
  this->MinNbrSinglets = minNbrSinglets; 
  this->NbrFermionsUp = 0;
  this->NbrFermionsDown = 0;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MaxXMomentum = maxXMomentum;
  this->MaxYMomentum = maxYMomentum;
  this->NbrSitePerUnitCell = this->NbrSite /  (this->MaxYMomentum * this->MaxXMomentum);
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = this->MaxXMomentum;

  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->XMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = 2 * (this->NbrSite / this->MaxXMomentum);
  this->ComplementaryStateXShift = (2 * this->MaxMomentum) - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;
  this->TargetSpace = this;

  this->MaxYMomentum =  maxYMomentum;
  this->YMomentum = yMomentum % this->MaxYMomentum;
  this->NbrYMomentumBlocks = (2 * this->NbrSite) / this->StateXShift;
  this->StateYShift = 2 * (this->NbrSite / (this->MaxYMomentum * this->MaxXMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (0x1ul << this->StateYShift) - 0x1ul;
  this->YMomentumBlockMask = (0x1ul << this->YMomentumBlockSize) - 0x1ul;  
  this->YMomentumFullMask = 0x0ul;
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
    {
      this->YMomentumFullMask |= this->YMomentumMask << (i *  this->YMomentumBlockSize);
    }
  this->ComplementaryYMomentumFullMask = ~this->YMomentumFullMask; 

  this->NbrFermionsParity = (~((unsigned long) this->NbrFermions)) & 0x1ul;


  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSite - 1, 0);
  cout << "intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->GenerateSignLookUpTable();
      this->LargeHilbertSpaceDimension  = this->GenerateStates();
      this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension > 0l)
	{
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
	}
#endif
    }
}



// basic constructor when Sz is preserved
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// totalSpin = twice the total spin value
// xMomentum = momentum sector in the x direction
// maxXMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction
// memory = amount of memory granted for precalculations

FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets::FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets (int nbrFermions, int minNbrSinglets, int totalSpin, 
																  int nbrSite, int xMomentum, int maxXMomentum,
																  int yMomentum, int  maxYMomentum, unsigned long memory)
{  
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->SzFlag = true;
  this->TotalSpin = totalSpin;
  this->MinNbrSinglets = minNbrSinglets; 
  this->NbrFermionsUp = (totalSpin + this->NbrFermions) >> 1;
  this->NbrFermionsDown = this->NbrFermions - this->NbrFermionsUp;
  this->NbrSite = nbrSite;
  this->MaxMomentum =  this->NbrSite;
  this->NbrMomentum = this->MaxMomentum + 1;
  this->MaxXMomentum = maxXMomentum;
  this->MaxYMomentum = maxYMomentum;
  this->NbrSitePerUnitCell = this->NbrSite /  (this->MaxYMomentum * this->MaxXMomentum);
  this->NbrFermionStates = 2 * this->NbrMomentum;
  this->MomentumModulo = this->MaxXMomentum;

  this->MomentumIncrement = (this->NbrFermions * this->StateShift/2) % this->MomentumModulo;
  this->ComplementaryStateShift = 2 * this->MaxMomentum - this->StateShift;
  this->MomentumMask = (0x1ul << this->StateShift) - 0x1ul;

  this->XMomentum = xMomentum % this->MaxXMomentum;
  this->StateXShift = 2 * (this->NbrSite / this->MaxXMomentum);
  this->ComplementaryStateXShift = (2 * this->MaxMomentum) - this->StateXShift;
  this->XMomentumMask = (0x1ul << this->StateXShift) - 0x1ul;

  this->MaxYMomentum =  maxYMomentum;
  this->YMomentum = yMomentum % this->MaxYMomentum;
  this->NbrYMomentumBlocks = (2 * this->NbrSite) / this->StateXShift;
  this->StateYShift = 2 * (this->NbrSite / (this->MaxYMomentum * this->MaxXMomentum));
  this->YMomentumBlockSize = this->StateYShift * this->MaxYMomentum;
  this->ComplementaryStateYShift = this->YMomentumBlockSize - this->StateYShift;
  this->YMomentumMask = (0x1ul << this->StateYShift) - 0x1ul;
  this->YMomentumBlockMask = (0x1ul << this->YMomentumBlockSize) - 0x1ul;  
  this->YMomentumFullMask = 0x0ul;
  for (int i = 0; i < this->NbrYMomentumBlocks; ++i)
    {
      this->YMomentumFullMask |= this->YMomentumMask << (i *  this->YMomentumBlockSize);
    }
  this->ComplementaryYMomentumFullMask = ~this->YMomentumFullMask; 

  this->NbrFermionsParity = (~((unsigned long) this->NbrFermions)) & 0x1ul;
  this->TargetSpace = this;

  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->NbrSite - 1, this->NbrFermionsUp, 0);
  cout << "intermediate Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->GenerateSignLookUpTable();
      this->LargeHilbertSpaceDimension  = this->GenerateStates();
      this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
      cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
      if (this->LargeHilbertSpaceDimension > 0l)
	{
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
	}
#endif
    }
}



// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets::FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets(const FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets& fermions)
{
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->NbrSite = fermions.NbrSite;
  this->NbrFermionsParity = fermions.NbrFermionsParity;
  this->MinNbrSinglets = fermions.MinNbrSinglets; 

  this->NbrSitePerUnitCell = fermions.NbrSitePerUnitCell;
  this->MaxXMomentum = fermions.MaxXMomentum;
  this->XMomentum = fermions.XMomentum;
  this->StateXShift = fermions.StateXShift;
  this->ComplementaryStateXShift = fermions.ComplementaryStateXShift;
  this->XMomentumMask = fermions.XMomentumMask;
  this->MaxYMomentum = fermions.MaxYMomentum;
  this->YMomentum = fermions.YMomentum;
  this->NbrYMomentumBlocks = fermions.NbrYMomentumBlocks;
  this->StateYShift = fermions.StateYShift;
  this->YMomentumBlockSize = fermions.YMomentumBlockSize;
  this->ComplementaryStateYShift = fermions.ComplementaryStateYShift;
  this->YMomentumMask = fermions.YMomentumMask;
  this->YMomentumBlockMask = fermions.YMomentumBlockMask;  
  this->YMomentumFullMask = fermions.YMomentumFullMask;
  this->ComplementaryYMomentumFullMask = fermions.ComplementaryYMomentumFullMask; 
  if (fermions.TargetSpace != &fermions)
    {
      this->TargetSpace = fermions.TargetSpace;
    }
  else
    {
      this->TargetSpace = this;
    }
  this->NbrFermionStates = fermions.NbrFermionStates;
  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
}

// destructor
//

FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets::~FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets& FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets::operator = (const FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateHighestBit;

      delete[] this->LookUpTableShift;
      for (int i = 0; i < this->NbrMomentum; ++i)
	delete[] this->LookUpTable[i];
      delete[] this->LookUpTable;

      delete[] this->SignLookUpTable;
      delete[] this->NbrParticleLookUpTable;

      for (int i = 1; i <= this->MaxMomentum ; ++i)
	delete[] this->RescalingFactors[i];
      delete[] this->RescalingFactors;
      delete[] this->NbrStateInOrbit;
    }
  this->NbrFermionsUp = fermions.NbrFermionsUp;
  this->NbrFermionsDown = fermions.NbrFermionsDown;
  this->NbrFermions = fermions.NbrFermions;  
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalSpin = fermions.TotalSpin;
  this->SzFlag = fermions.SzFlag;
  this->NbrSite = fermions.NbrSite;
  this->NbrFermionsParity = fermions.NbrFermionsParity;
  this->MinNbrSinglets = fermions.MinNbrSinglets; 

  this->NbrSitePerUnitCell = fermions.NbrSitePerUnitCell;
  this->MaxXMomentum = fermions.MaxXMomentum;
  this->XMomentum = fermions.XMomentum;
  this->StateXShift = fermions.StateXShift;
  this->ComplementaryStateXShift = fermions.ComplementaryStateXShift;
  this->XMomentumMask = fermions.XMomentumMask;
  this->MaxYMomentum = fermions.MaxYMomentum;
  this->YMomentum = fermions.YMomentum;
  this->NbrYMomentumBlocks = fermions.NbrYMomentumBlocks;
  this->StateYShift = fermions.StateYShift;
  this->YMomentumBlockSize = fermions.YMomentumBlockSize;
  this->ComplementaryStateYShift = fermions.ComplementaryStateYShift;
  this->YMomentumMask = fermions.YMomentumMask;
  this->YMomentumBlockMask = fermions.YMomentumBlockMask;  
  this->YMomentumFullMask = fermions.YMomentumFullMask;
  this->ComplementaryYMomentumFullMask = fermions.ComplementaryYMomentumFullMask; 
  if (fermions.TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;

  this->NbrFermionStates = fermions.NbrFermionStates;
  this->MaxMomentum = fermions.MaxMomentum;
  this->NbrMomentum = fermions.NbrMomentum;
  this->MomentumModulo = fermions.MomentumModulo;
  this->MomentumIncrement = fermions.MomentumIncrement;
  this->StateShift = fermions.StateShift;
  this->ComplementaryStateShift = fermions.ComplementaryStateShift;
  this->MomentumMask = fermions.MomentumMask;

  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->StateHighestBit = fermions.StateHighestBit;

  this->MaximumLookUpShift = fermions.MaximumLookUpShift;
  this->LookUpTableMemorySize = fermions.LookUpTableMemorySize;
  this->LookUpTableShift = fermions.LookUpTableShift;
  this->LookUpTable = fermions.LookUpTable;

  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->NbrParticleLookUpTable = fermions.NbrParticleLookUpTable;

  this->RescalingFactors = fermions.RescalingFactors;
  this->NbrStateInOrbit = fermions.NbrStateInOrbit;

  this->ReorderingSign = fermions.ReorderingSign;

  this->Flag = fermions.Flag;

  this->LargeHilbertSpaceDimension = (long) this->HilbertSpaceDimension;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets::Clone()
{
  return new FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets(*this);
}

// generate all states corresponding to the constraints
//
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets::GenerateStates()
{
  this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
  if (this->SzFlag == false)
    this->RawGenerateStates(this->NbrFermions, this->NbrSite - 1, 0, 0l);
  else
    this->RawGenerateStates(this->NbrFermions, this->NbrSite - 1, this->NbrFermionsUp, 0, 0l);
  long TmpLargeHilbertSpaceDimension = 0l;
  int NbrTranslationX;
  int NbrTranslationY;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if ((this->FindCanonicalForm(this->StateDescription[i], NbrTranslationX, NbrTranslationY) == this->StateDescription[i]))
	{
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
  if (TmpLargeHilbertSpaceDimension == 0l)
    return 0l;
  unsigned long* TmpStateDescription = new unsigned long [TmpLargeHilbertSpaceDimension];  
  this->NbrStateInOrbit = new int [TmpLargeHilbertSpaceDimension];
  this->ReorderingSign = new unsigned long [TmpLargeHilbertSpaceDimension];
  TmpLargeHilbertSpaceDimension = 0l;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      if (this->StateDescription[i] != 0x0ul)
	{
	  TmpStateDescription[TmpLargeHilbertSpaceDimension] = this->StateDescription[i];
	  this->NbrStateInOrbit[TmpLargeHilbertSpaceDimension] = this->FindOrbitSize(this->StateDescription[i]);
	  unsigned long& TmpSign = this->ReorderingSign[TmpLargeHilbertSpaceDimension];
	  TmpSign = 0x0ul;	  
	  int Index = 1;
	  unsigned long TmpState =  this->StateDescription[i];
	  for (int m = 0; m < this->MaxYMomentum; ++m)
	    {
	      unsigned long TmpState2 = TmpState;
	      for (int n = 1; n < this->MaxXMomentum; ++n)
		{
		  TmpSign |= (this->GetSignAndApplySingleXTranslation(TmpState2) << Index) ^ ((TmpSign & (0x1ul << (Index - 1))) << 1);
		  ++Index;
		}
	      TmpSign |= ((this->GetSignAndApplySingleYTranslation(TmpState) << Index) 
			  ^ ((TmpSign & (0x1ul << (Index - this->MaxXMomentum))) << this->MaxXMomentum));
	      ++Index;
	    }
	  ++TmpLargeHilbertSpaceDimension;
	}
    }
  delete[] this->StateDescription;
  this->StateDescription = TmpStateDescription;

  this->StateHighestBit = new int [TmpLargeHilbertSpaceDimension];  
  int CurrentMaxMomentum = (2 * this->MaxMomentum) + 1;
  while (((this->StateDescription[0] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
    --CurrentMaxMomentum;
  this->StateHighestBit[0] = CurrentMaxMomentum;
  for (long i = 1l; i < TmpLargeHilbertSpaceDimension; ++i)
    {
      while (((this->StateDescription[i] >> CurrentMaxMomentum) & 0x1ul) == 0x0ul)
	--CurrentMaxMomentum;
      this->StateHighestBit[i] = CurrentMaxMomentum;
    }
  return TmpLargeHilbertSpaceDimension;
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// nbrSinglets = number of on-site singlets
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets::RawGenerateStates(int nbrFermions, int currentSite, int nbrSinglets, long pos)
{
  if (nbrFermions == 0)
    {
      if (nbrSinglets >= this->MinNbrSinglets)
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else
	{
	  return pos;
	}
    }
  if ((currentSite < 0) || (nbrFermions < 0))
    return pos;
  if (nbrFermions == 1)
    {
       if (nbrSinglets >= this->MinNbrSinglets)
	{
	  for (int j = currentSite; j >= 0; --j)
	    {
	      this->StateDescription[pos] = 0x2ul << (j << 1);
	      ++pos;
	      this->StateDescription[pos] = 0x1ul << (j << 1);
	      ++pos;
	    }
	}
      return pos;
    }
  long TmpPos = this->RawGenerateStates(nbrFermions - 2, currentSite - 1, nbrSinglets + 1, pos);
  unsigned long Mask = 0x3ul << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->RawGenerateStates(nbrFermions - 1, currentSite - 1, nbrSinglets, pos);
  Mask = 0x2ul << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
   TmpPos = this->RawGenerateStates(nbrFermions - 1, currentSite - 1, nbrSinglets, pos);
   Mask = 0x1ul << (currentSite << 1);
   for (; pos < TmpPos; ++pos)
     this->StateDescription[pos] |= Mask;
  return this->RawGenerateStates(nbrFermions, currentSite - 1, nbrSinglets, pos);
}

// generate all states corresponding to the constraints
// 
// nbrFermions = number of fermions
// currentSite = current site index
// nbrSpinUp = number of fermions with spin up
// nbrSinglets = number of on-site singlets
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets::RawGenerateStates(int nbrFermions, int currentSite, int nbrSpinUp, int nbrSinglets, long pos)
{
  if ((nbrFermions == 0) && (nbrSpinUp == 0))
    {
      if (nbrSinglets >= this->MinNbrSinglets)
	{
	  this->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else
	{
	  return pos;
	}
    }
  if ((currentSite < 0) || (nbrFermions < 0) || (nbrSpinUp > nbrFermions) || (nbrSpinUp < 0))
    return pos;
  if (nbrFermions == 1)
    {
      if (nbrSinglets >= this->MinNbrSinglets)
	{
	  if (nbrSpinUp == 1)
	    {
	      for (int j = currentSite; j >= 0; --j)
		{
		  this->StateDescription[pos] = 0x2ul << (j << 1);
		  ++pos;
		}
	    }
	  else
	    {
	      for (int j = currentSite; j >= 0; --j)
		{
		  this->StateDescription[pos] = 0x1ul << (j << 1);
		  ++pos;
		}
	    }
	}
      return pos;
    }
  long TmpPos = this->RawGenerateStates(nbrFermions - 2, currentSite - 1, nbrSpinUp - 1, nbrSinglets + 1, pos);
  unsigned long Mask = 0x3ul << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->RawGenerateStates(nbrFermions - 1, currentSite - 1, nbrSpinUp - 1, nbrSinglets, pos);
  Mask = 0x2ul << ((currentSite) << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  TmpPos = this->RawGenerateStates(nbrFermions - 1, currentSite - 1, nbrSpinUp, nbrSinglets, pos);
  Mask = 0x1ul << (currentSite << 1);
  for (; pos < TmpPos; ++pos)
    this->StateDescription[pos] |= Mask;
  return this->RawGenerateStates(nbrFermions, currentSite - 1, nbrSpinUp, nbrSinglets, pos);   
}


// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// currentSite = current site index
// nbrSinglets = number of on-site singlets
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets::EvaluateHilbertSpaceDimension(int nbrFermions, int currentSite, int nbrSinglets)
{
  if (nbrFermions == 0)
    {
      if (nbrSinglets >= this->MinNbrSinglets)
	{
	  return 1l;
	}
      else
	{
	  return 0l;
	}
    }
  if ((currentSite < 0) || (nbrFermions < 0))
    return 0l;
  if (nbrFermions == 1)
    {
       if (nbrSinglets >= this->MinNbrSinglets)
	{
	  return ((long) (2 * (currentSite + 1)));
	}
      return 0l;
    }
  long TmpPos = this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentSite - 1, nbrSinglets + 1);
  TmpPos += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentSite - 1, nbrSinglets);
  TmpPos += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentSite - 1, nbrSinglets);
  TmpPos += this->EvaluateHilbertSpaceDimension(nbrFermions, currentSite - 1, nbrSinglets);
  return TmpPos;
}

// evaluate Hilbert space dimension with a fixed number of fermions with spin up
//
// nbrFermions = number of fermions
// currentSite = current site index
// nbrSpinUp = number of fermions with spin up
// nbrSinglets = number of on-site singlets
// return value = Hilbert space dimension

long FermionOnLatticeWithSpinRealSpaceAnd2DTranslationMinNbrSinglets::EvaluateHilbertSpaceDimension(int nbrFermions, int currentSite, int nbrSpinUp, int nbrSinglets)
{
  if ((nbrFermions == 0) && (nbrSpinUp == 0))
    {
      if (nbrSinglets >= this->MinNbrSinglets)
	{
	  return 1l;
	}
      else
	{
	  return 0l;
	}
    }
  if ((currentSite < 0) || (nbrFermions < 0) || (nbrSpinUp > nbrFermions) || (nbrSpinUp < 0))
    return 0l;
  if (nbrFermions == 1)
    {
      if (nbrSinglets >= this->MinNbrSinglets)
	{
	  return ((long) (currentSite + 1));
	}
      return 0l;
    }
  long TmpPos = this->EvaluateHilbertSpaceDimension(nbrFermions - 2, currentSite - 1, nbrSpinUp - 1, nbrSinglets + 1);
  TmpPos += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentSite - 1, nbrSpinUp - 1, nbrSinglets);
  TmpPos += this->EvaluateHilbertSpaceDimension(nbrFermions - 1, currentSite - 1, nbrSpinUp, nbrSinglets);
  TmpPos += this->EvaluateHilbertSpaceDimension(nbrFermions, currentSite - 1, nbrSpinUp, nbrSinglets);   
  return TmpPos;
}


