////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                         class of bosons on square lattice                  //
//                                  in momentum space                         //
//                                                                            //
//                        last modification : 16/09/2011                      //
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
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
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
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"

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

BosonOnSquareLatticeMomentumSpace::BosonOnSquareLatticeMomentumSpace ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// kxMomentum = momentum along the x direction
// kyMomentum = momentum along the y direction
// memory = amount of memory granted for precalculations

BosonOnSquareLatticeMomentumSpace::BosonOnSquareLatticeMomentumSpace (int nbrBosons, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxMomentum = kxMomentum;
  this->KyMomentum = kyMomentum;
  this->LzMax = this->NbrSiteX * this->NbrSiteY;
  this->NbrLzValue = this->LzMax + 1;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0);
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->FermionBasis = new FermionOnSphere(0, 0, this->LzMax);
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->FermionBasis->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      this->FermionBasis->StateLzMax = new int [this->LargeHilbertSpaceDimension];  
      this->LargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->NbrSiteX - 1, this->NbrSiteY - 1, 0, 0, this->LzMax + this->NbrBosons, 0l);
      int TmpLzMax = this->LzMax;
      for (long i = 0l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  while ((this->FermionBasis->StateDescription[i] >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  this->FermionBasis->StateLzMax[i] = TmpLzMax;
	}
      this->FermionBasis->GenerateLookUpTable(memory);      
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
      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->FermionBasis-> LookUpTableMemorySize * sizeof(int);
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

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnSquareLatticeMomentumSpace::BosonOnSquareLatticeMomentumSpace(const BosonOnSquareLatticeMomentumSpace& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

BosonOnSquareLatticeMomentumSpace::~BosonOnSquareLatticeMomentumSpace ()
{
}

// assignement (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnSquareLatticeMomentumSpace& BosonOnSquareLatticeMomentumSpace::operator = (const BosonOnSquareLatticeMomentumSpace& bosons)
{
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrSiteX = bosons.NbrSiteX;
  this->NbrSiteY = bosons.NbrSiteY;
  this->KxMomentum = bosons.KxMomentum;
  this->KyMomentum = bosons.KyMomentum;
  this->NbrLzValue = bosons.NbrLzValue;
  this->FermionBasis = (FermionOnSphere*) bosons.FermionBasis->Clone();
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnSquareLatticeMomentumSpace::Clone()
{
  return new BosonOnSquareLatticeMomentumSpace(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured
bool BosonOnSquareLatticeMomentumSpace::WriteHilbertSpace (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->HilbertSpaceDimension);
  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);
  WriteLittleEndian(File, this->NbrBosons);
  WriteLittleEndian(File, this->NbrSiteX);
  WriteLittleEndian(File, this->NbrSiteY);
  WriteLittleEndian(File, this->KxMomentum);
  WriteLittleEndian(File, this->KyMomentum);
  if (this->HilbertSpaceDimension != 0)
    {
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->FermionBasis->StateDescription[i]);
      for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->FermionBasis->StateLzMax[i]);
    }
  else
    {
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->FermionBasis->StateDescription[i]);
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	WriteLittleEndian(File, this->FermionBasis->StateLzMax[i]);
    }
  File.close();
  return true;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& BosonOnSquareLatticeMomentumSpace::PrintState (ostream& Str, int state)
{
  this->FermionToBoson(this->FermionBasis->StateDescription[state], this->FermionBasis->StateLzMax[state], this->TemporaryState, this->TemporaryStateLzMax);
  Str << "[";
  for (int i = 0; i < this->NbrLzValue; ++i)
    {
      if (this->TemporaryState[i] > 0)
	{
	  int TmpKx = i / this->NbrSiteY;
	  int TmpKy = i % this->NbrSiteY;
	  for (int j = 0; j < this->TemporaryState[i]; ++j)
	    Str << "(" << TmpKx << "," << TmpKy << ")";
	}
    }
  Str << "]";
  return Str;
}

// generate all states corresponding to the constraints
// 
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// currentFermionicPosition = current fermionic position within the state description
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored
  
long BosonOnSquareLatticeMomentumSpace::GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int currentFermionicPosition, long pos)
{

  if (nbrBosons < 0)
    return pos;
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	{
	  this->FermionBasis->StateDescription[pos] = 0x0ul;	  
	  return (pos + 1l);
	}
      else	
	return pos;
    }
  if (currentKx < 0)
    return pos;

  for (int k = nbrBosons; k > 0; --k)
    {
      long TmpPos = this->GenerateStates(nbrBosons - k, currentKx, currentKy - 1, currentTotalKx + (k * currentKx), currentTotalKy + (k * currentKy), currentFermionicPosition - k - 1, pos);
      unsigned long Mask = ((0x1ul << k) - 0x1ul) << (currentFermionicPosition - k - 1);
      for (; pos < TmpPos; ++pos)
	this->FermionBasis->StateDescription[pos] |= Mask;
    }
  return this->GenerateStates(nbrBosons, currentKx, currentKy - 1, currentTotalKx, currentTotalKy, currentFermionicPosition - 1, pos);
};


// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// currentKx = current momentum along x for a single particle
// currentKy = current momentum along y for a single particle
// currentTotalKx = current total momentum along x
// currentTotalKy = current total momentum along y
// return value = Hilbert space dimension

long BosonOnSquareLatticeMomentumSpace::EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy)
{
  if (nbrBosons < 0)
    return 0l;
  if (currentKy < 0)
    {
      currentKy = this->NbrSiteY - 1;
      currentKx--;
    }
  if (nbrBosons == 0)
    {
      if (((currentTotalKx % this->NbrSiteX) == this->KxMomentum) && ((currentTotalKy % this->NbrSiteY) == this->KyMomentum))
	return 1l;
      else	
	return 0l;
    }
  if (currentKx < 0)
    return 0l;
  long Count = 0;
  for (int k = nbrBosons; k >= 0; --k)
    Count += this->EvaluateHilbertSpaceDimension(nbrBosons - k, currentKx, currentKy - 1, currentTotalKx + (k * currentKx), currentTotalKy + (k * currentKy));
  return Count;
}


// evaluate a density matrix of a subsystem of the whole system described by a given ground state. The density matrix is only evaluated in a given momentum sector and fixed number of particles. 
// 
// subsytemSizeX = number of sites along the x direction that belong to the subsytem
// subsytemSizeY = number of sites along the y direction that belong to the subsytem
// subsytemStartX = x momentum marking the beginning of the rectangluar subsystem
// subsytemStartY = y momentum marking the beginning of the rectangluar subsystem
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = Kx momentum sector in which the entanglement matrix has to be evaluated 
// kySector = Ky momentum sector in which the entanglement matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = density matrix of the subsytem

HermitianMatrix BosonOnSquareLatticeMomentumSpace::EvaluatePartialDensityMatrixMomentumSpace (int subsytemSizeX, int subsytemSizeY, int subsytemStartX, int subsytemStartY, int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState)
{
  subsytemSizeX = 1 + (-1 + subsytemSizeX + this->NbrSiteX) % this->NbrSiteX; // enforce subsystemSize in [1,NbrSiteX]
  subsytemSizeY = 1 + (-1 + subsytemSizeY + this->NbrSiteY) % this->NbrSiteY;
  if ((subsytemSizeX != this->NbrSiteX) && (subsytemSizeY != this->NbrSiteY))
    {
      cout << "not implemented yet!" << endl;
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }

  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementarySubsytemSizeX;
  int ComplementarySubsytemSizeY;
  int ComplementarySubsytemStartX;
  int ComplementarySubsytemStartY;
  if (subsytemSizeX == this->NbrSiteX)
    {
      ComplementarySubsytemSizeX = this->NbrSiteX;
      ComplementarySubsytemStartX = subsytemStartX;
      ComplementarySubsytemSizeY = this->NbrSiteY - subsytemSizeY;
      ComplementarySubsytemStartY = (subsytemStartY + subsytemSizeY) % this->NbrSiteY;
    }
  else
    {
      ComplementarySubsytemSizeX = this->NbrSiteX - subsytemSizeX;
      ComplementarySubsytemStartX = (subsytemStartX + subsytemSizeX) % this->NbrSiteX;
      ComplementarySubsytemSizeY = this->NbrSiteY;
      ComplementarySubsytemStartY = subsytemStartY;
    }

  int ComplementaryKxMomentum = (this->KxMomentum - kxSector + this->NbrSiteX) % this->NbrSiteX;
  int ComplementaryKyMomentum = (this->KyMomentum - kySector + this->NbrSiteY) % this->NbrSiteY;

//   BosonOnSquareLatticeNonPeriodicMomentumSpace SubsytemSpace (nbrParticleSector, this->NbrSiteX, this->NbrSiteY, subsytemSizeX, subsytemSizeY, subsytemStartX, subsytemStartY, kxSector, kySector);
//   BosonOnSquareLatticeNonPeriodicMomentumSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSiteX, this->NbrSiteY, ComplementarySubsytemSizeX, ComplementarySubsytemSizeY, ComplementarySubsytemStartX, ComplementarySubsytemStartY, ComplementaryKxMomentum, ComplementaryKyMomentum);

//   if (SubsytemSpace.HilbertSpaceDimension>0 && ComplementarySpace.HilbertSpaceDimension>0)
//   {
//       cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension
//            << "; complementary Hilbert space dimension = " << ComplementarySpace.HilbertSpaceDimension << endl;
//       HermitianMatrix TmpDensityMatrix (SubsytemSpace.HilbertSpaceDimension, true);
//       int nbrNonZeroElements = this->EvaluatePartialDensityMatrixMomentumSpaceCore (0, ComplementarySpace.HilbertSpaceDimension, &ComplementarySpace, &SubsytemSpace, groundState, &TmpDensityMatrix);
//       if (nbrNonZeroElements>0)
//           return TmpDensityMatrix;
//   }

//   HermitianMatrix TmpDensityMatrixZero (SubsytemSpace.HilbertSpaceDimension, true);

  HermitianMatrix TmpDensityMatrixZero (1, true);
  return TmpDensityMatrixZero;
}

// core part of the evaluation density matrix momentum space partition calculation
// 
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnSquareLatticeMomentumSpace::EvaluatePartialDensityMatrixMomentumSpaceCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace, ComplexVector& groundState,  HermitianMatrix* densityMatrix)
{
//   BosonOnSquareLatticeMomentumSpace* TmpHilbertSpace =  (BosonOnSquareLatticeNonPeriodicMomentumSpace*) complementaryHilbertSpace;
//   BosonOnSquareLatticeMomentumSpace* TmpDestinationHilbertSpace =  (BosonOnSquareLatticeNonPeriodicMomentumSpace*) destinationHilbertSpace;
//   int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
//   int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
//   Complex* TmpStateCoefficient = new Complex [TmpDestinationHilbertSpace->HilbertSpaceDimension];
//   int MaxIndex = minIndex + nbrIndex;
//   long TmpNbrNonZeroElements = 0l;
  
//   for (; minIndex < MaxIndex; ++minIndex)    
//     {
//       int Pos = 0;
//       unsigned long TmpState = TmpHilbertSpace->StateDescription[minIndex];
//       for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
// 	{
// 	  unsigned long TmpState2 = TmpDestinationHilbertSpace->StateDescription[j];
// 	  if ((TmpState & TmpState2) == 0x0ul)
// 	    {
//  	      int TmpLzMax = this->LzMax;
// 	      unsigned long TmpState3 = TmpState | TmpState2;
// 	      while ((TmpState3 >> TmpLzMax) == 0x0ul)
// 		--TmpLzMax;
// 	      int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
// 	      if (TmpPos != this->HilbertSpaceDimension)
// 		{
// 		  double Coefficient = 1.0;
// 		  unsigned long Sign = 0x0ul;
// 		  int Pos2 = TmpDestinationHilbertSpace->LzMax;
// 		  while ((Pos2 > 0) && (TmpState2 != 0x0ul))
// 		    {
// 		      while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
// 			--Pos2;
// 		      TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
// #ifdef  __64_BITS__
// 		      TmpState3 ^= TmpState3 >> 32;
// #endif	
// 		      TmpState3 ^= TmpState3 >> 16;
// 		      TmpState3 ^= TmpState3 >> 8;
// 		      TmpState3 ^= TmpState3 >> 4;
// 		      TmpState3 ^= TmpState3 >> 2;
// 		      TmpState3 ^= TmpState3 >> 1;
// 		      Sign ^= TmpState3;
// 		      TmpState2 &= ~(0x1ul << Pos2);
// 		      --Pos2;
// 		    }
//  		  if ((Sign & 0x1ul) == 0x0ul)		  
//  		    Coefficient *= 1.0;
//  		  else
//  		    Coefficient *= -1.0;
// 		  TmpStatePosition[Pos] = TmpPos;
// 		  TmpStatePosition2[Pos] = j;
// 		  TmpStateCoefficient[Pos] = Coefficient;
// 		  ++Pos;
// 		}
// 	    }
// 	}
//       if (Pos != 0)
// 	{
// 	  ++TmpNbrNonZeroElements;
// 	  for (int j = 0; j < Pos; ++j)
// 	    {
// 	      int Pos2 = TmpStatePosition2[j];
// 	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]]) * TmpStateCoefficient[j];
// 	      for (int k = 0; k < Pos; ++k)
// 		if (TmpStatePosition2[k] >= Pos2)
// 		  {
// 		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
// 		  }
// 	    }
// 	}
//     }
//   delete[] TmpStatePosition;
//   delete[] TmpStatePosition2;
//   delete[] TmpStateCoefficient;
//  return TmpNbrNonZeroElements;
  return 0l;
}

// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix is only evaluated in a given momentum sector and fixed number of particles
// 
// subsytemSizeX = number of sites along the x direction that belong to the subsytem
// subsytemSizeY = number of sites along the y direction that belong to the subsytem
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = Kx momentum sector in which the entanglement matrix has to be evaluated 
// kySector = Ky momentum sector in which the entanglement matrix has to be evaluated 
// groundState = reference on the total system ground state
// return value = entanglement matrix of the subsytem
  
ComplexMatrix BosonOnSquareLatticeMomentumSpace::EvaluatePartialEntanglementMatrixMomentumSpace (int subsytemSizeX, int subsytemSizeY, int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState)
{
  ComplexMatrix TmpEntanglementMatrix (1,2, true);
  return TmpEntanglementMatrix;
}
  
// evaluate an entanglement matrix of a subsystem of the whole system described by a given ground state, using particle partition. The entanglement matrix is only evaluated in a given Lz sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = Kx momentum sector in which the entanglement matrix has to be evaluated 
// kySector = Ky momentum sector in which the entanglement matrix has to be evaluated 
// groundState = reference on the total system ground state
// removeBinomialCoefficient = remove additional binomial coefficient in case the particle entanglement matrix has to be used for real space cut
// return value = entanglement matrix of the subsytem (return a wero dimension matrix if the entanglement matrix is equal to zero)

ComplexMatrix BosonOnSquareLatticeMomentumSpace::EvaluatePartialEntanglementMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, bool removeBinomialCoefficient)
{
  ComplexMatrix TmpEntanglementMatrix;
  return TmpEntanglementMatrix;
}
 
// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given Lz sector.
// 
// nbrBosonSector = number of particles that belong to the subsytem 
// lzSector = Lz sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnSquareLatticeMomentumSpace::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, int kxSector, int kySector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      if ((kxSector == 0) && (kySector == 0))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  if (nbrParticleSector == this->NbrBosons)
    {
      if ((kxSector == this->KxMomentum) && (kySector == this->KyMomentum))
	{
	  HermitianMatrix TmpDensityMatrix(1, true);
	  TmpDensityMatrix(0, 0) = 1.0;
	  return TmpDensityMatrix;
	}
      else
	{
	  HermitianMatrix TmpDensityMatrix;
	  return TmpDensityMatrix;
	}
    }
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  int ComplementaryKxMomentum = (this->KxMomentum - kxSector) % this->NbrSiteX;
  int ComplementaryKyMomentum = (this->KyMomentum - kySector) % this->NbrSiteY;
  if (ComplementaryKxMomentum < 0)
    ComplementaryKxMomentum += this->NbrSiteX;
  if (ComplementaryKyMomentum < 0)
    ComplementaryKyMomentum += this->NbrSiteY;
  cout << "kx = " << this->KxMomentum << " " << kxSector << " " << ComplementaryKxMomentum << endl;
  cout << "ky = " << this->KyMomentum << " " << kySector << " " << ComplementaryKyMomentum << endl;
  BosonOnSquareLatticeMomentumSpace SubsytemSpace (nbrParticleSector, this->NbrSiteX, this->NbrSiteY, kxSector, kySector);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnSquareLatticeMomentumSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSiteX, this->NbrSiteY, ComplementaryKxMomentum, ComplementaryKyMomentum);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;


  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}


// core part of the evaluation density matrix particle partition calculation
// 
// minIndex = first index to consider in complementary Hilbert space
// nbrIndex = number of indices to consider in complementary Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e part B)
// destinationHilbertSpace = pointer to the destination Hilbert space (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnSquareLatticeMomentumSpace::EvaluatePartialDensityMatrixParticlePartitionCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
											     ComplexVector& groundState,  HermitianMatrix* densityMatrix)
{
  BosonOnSquareLatticeMomentumSpace* TmpHilbertSpace =  (BosonOnSquareLatticeMomentumSpace*) complementaryHilbertSpace;
  BosonOnSquareLatticeMomentumSpace* TmpDestinationHilbertSpace =  (BosonOnSquareLatticeMomentumSpace*) destinationHilbertSpace;
  int* TmpStatePosition = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int* TmpStatePosition2 = new int [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  Complex* TmpStateCoefficient = new Complex [TmpDestinationHilbertSpace->HilbertSpaceDimension];
  int MaxIndex = minIndex + nbrIndex;
  long TmpNbrNonZeroElements = 0l;
  BinomialCoefficients TmpBinomial (this->NbrBosons);
  double TmpInvBinomial = 1.0 / sqrt(TmpBinomial(this->NbrBosons, TmpDestinationHilbertSpace->NbrBosons));
  
//   for (; minIndex < MaxIndex; ++minIndex)    
//     {
//       int Pos = 0;
//       unsigned long TmpState = TmpHilbertSpace->StateDescription[minIndex];
//       for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
// 	{
// 	  unsigned long TmpState2 = TmpDestinationHilbertSpace->StateDescription[j];
// 	  if ((TmpState & TmpState2) == 0x0ul)
// 	    {
//  	      int TmpLzMax = this->LzMax;
// 	      unsigned long TmpState3 = TmpState | TmpState2;
// 	      while ((TmpState3 >> TmpLzMax) == 0x0ul)
// 		--TmpLzMax;
// 	      int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
// 	      if (TmpPos != this->HilbertSpaceDimension)
// 		{
// 		  double Coefficient = TmpInvBinomial;
// 		  unsigned long Sign = 0x0ul;
// 		  int Pos2 = TmpDestinationHilbertSpace->LzMax;
// 		  while ((Pos2 > 0) && (TmpState2 != 0x0ul))
// 		    {
// 		      while (((TmpState2 >> Pos2) & 0x1ul) == 0x0ul)
// 			--Pos2;
// 		      TmpState3 = TmpState & ((0x1ul << (Pos2 + 1)) - 1ul);
// #ifdef  __64_BITS__
// 		      TmpState3 ^= TmpState3 >> 32;
// #endif	
// 		      TmpState3 ^= TmpState3 >> 16;
// 		      TmpState3 ^= TmpState3 >> 8;
// 		      TmpState3 ^= TmpState3 >> 4;
// 		      TmpState3 ^= TmpState3 >> 2;
// 		      TmpState3 ^= TmpState3 >> 1;
// 		      Sign ^= TmpState3;
// 		      TmpState2 &= ~(0x1ul << Pos2);
// 		      --Pos2;
// 		    }
//  		  if ((Sign & 0x1ul) == 0x0ul)		  
//  		    Coefficient *= 1.0;
//  		  else
//  		    Coefficient *= -1.0;
// 		  TmpStatePosition[Pos] = TmpPos;
// 		  TmpStatePosition2[Pos] = j;
// 		  TmpStateCoefficient[Pos] = Coefficient;
// 		  ++Pos;
// 		}
// 	    }
// 	}
//       if (Pos != 0)
// 	{
// 	  ++TmpNbrNonZeroElements;
// 	  for (int j = 0; j < Pos; ++j)
// 	    {
// 	      int Pos2 = TmpStatePosition2[j];
// 	      Complex TmpValue = Conj(groundState[TmpStatePosition[j]]) * TmpStateCoefficient[j];
// 	      for (int k = 0; k < Pos; ++k)
// 		if (TmpStatePosition2[k] >= Pos2)
// 		  {
// 		    densityMatrix->AddToMatrixElement(Pos2, TmpStatePosition2[k], TmpValue * groundState[TmpStatePosition[k]] * TmpStateCoefficient[k]);
// 		  }
// 	    }
// 	}
//     }
  delete[] TmpStatePosition;
  delete[] TmpStatePosition2;
  delete[] TmpStateCoefficient;
  return TmpNbrNonZeroElements;
}

// find state index from an array
//
// stateDescription = array describing the state (stored as kx1,ky1,kx2,ky2,...)
// return value = corresponding index, -1 if an error occured

int BosonOnSquareLatticeMomentumSpace::FindStateIndexFromArray(int* stateDescription)
{
  for (int i = 0; i <= this->LzMax; ++i)
    this->TemporaryState[i] = 0l;
  for (int i = 0; i < this->NbrBosons; ++i)
    ++this->TemporaryState[((stateDescription[i << 1] * this->NbrSiteY) + stateDescription[1 + (i << 1)])];
  unsigned long TmpState = this->BosonToFermion(this->TemporaryState, this->LzMax);
  int TmpLzMax = this->FermionBasis->LzMax;
  while ((TmpState >> TmpLzMax) == 0x0ul)
    --TmpLzMax;
  return this->FermionBasis->FindStateIndex(TmpState, TmpLzMax);
}
