////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                    class of fermions on lattice in real space              //
//                                                                            //
//                        last modification : 09/09/2014                      //
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
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"
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

BosonOnLatticeGutzwillerProjectionRealSpace::BosonOnLatticeGutzwillerProjectionRealSpace ()
{
  this->NbrBosons = 0;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->NbrSite = 0;
  this->LzMax = this->NbrSite - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->HilbertSpaceDimension = 0;
  this->StateDescription = 0;
  this->StateLzMax = 0;  
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// memory = amount of memory granted for precalculations

BosonOnLatticeGutzwillerProjectionRealSpace::BosonOnLatticeGutzwillerProjectionRealSpace (int nbrBosons, int nbrSite, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons);
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if (this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      long TmpLargeHilbertSpaceDimension = this->GenerateStates(this->NbrBosons, this->NbrSite - 1, 0l);
      this->StateLzMax = new int [this->LargeHilbertSpaceDimension];  
      int CurrentLzMax = this->NbrLzValue;
      while (((this->StateDescription[0] >> CurrentLzMax) & 0x1ul) == 0x0ul)
	--CurrentLzMax;
      this->StateLzMax[0] = CurrentLzMax;
      for (long i = 1l; i < this->LargeHilbertSpaceDimension; ++i)
	{
	  while (((this->StateDescription[i] >> CurrentLzMax) & 0x1ul) == 0x0ul)
	    --CurrentLzMax;
	  this->StateLzMax[i] = CurrentLzMax;
 	}
      if (TmpLargeHilbertSpaceDimension != this->LargeHilbertSpaceDimension)
	{
	  cout << "error while generating the Hilbert space, " << TmpLargeHilbertSpaceDimension << " generated states, should be " << this->LargeHilbertSpaceDimension << endl;
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
      UsedMemory = this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
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
// fermions = reference on the hilbert space to copy to copy

BosonOnLatticeGutzwillerProjectionRealSpace::BosonOnLatticeGutzwillerProjectionRealSpace(const BosonOnLatticeGutzwillerProjectionRealSpace & bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->NbrSite = bosons.NbrSite;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;  
  this->SignLookUpTable = bosons.SignLookUpTable;
  this->SignLookUpTableMask = bosons.SignLookUpTableMask;
  this->MaximumSignLookUp = bosons.MaximumSignLookUp;
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
}

// destructor
//

BosonOnLatticeGutzwillerProjectionRealSpace::~BosonOnLatticeGutzwillerProjectionRealSpace ()
{
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnLatticeGutzwillerProjectionRealSpace& BosonOnLatticeGutzwillerProjectionRealSpace::operator = (const BosonOnLatticeGutzwillerProjectionRealSpace& bosons)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->StateLzMax;
    }
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
  this->NbrSite = bosons.NbrSite;
  this->NbrLzValue = bosons.NbrLzValue;
  this->StateDescription = bosons.StateDescription;
  this->StateLzMax = bosons.StateLzMax;
  this->MaximumLookUpShift = bosons.MaximumLookUpShift;
  this->LookUpTableMemorySize = bosons.LookUpTableMemorySize;
  this->LookUpTableShift = bosons.LookUpTableShift;
  this->LookUpTable = bosons.LookUpTable;  
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnLatticeGutzwillerProjectionRealSpace::Clone()
{
  return new BosonOnLatticeGutzwillerProjectionRealSpace(*this);
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnLatticeGutzwillerProjectionRealSpace::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix(0, 0) = 1.0;
      return TmpDensityMatrix;
    }
  if (nbrParticleSector == this->NbrBosons)
    {
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix(0, 0) = 1.0;
      return TmpDensityMatrix;
    }
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  BosonOnLatticeGutzwillerProjectionRealSpace SubsytemSpace (nbrParticleSector, this->NbrSite);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnLatticeGutzwillerProjectionRealSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSite);
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

// evaluate the orbital cut entanglement matrix. The entanglement matrix is only evaluated for fixed number of particles
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// groundState = reference on the total system ground state
// keptOrbitals = array of orbitals that have to be kept, should be sorted from the smallest index to the largest index 
// nbrKeptOrbitals = array of orbitals that have to be kept
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = entanglement matrix of the subsytem

ComplexMatrix BosonOnLatticeGutzwillerProjectionRealSpace::EvaluatePartialEntanglementMatrix (int nbrParticleSector, int nbrKeptOrbitals, int* keptOrbitals, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  if ((nbrParticleSector >  nbrKeptOrbitals) || 
      (ComplementaryNbrParticles > (this->NbrSite - nbrKeptOrbitals)))
    {
      ComplexMatrix TmpEntanglementMatrix;
      return TmpEntanglementMatrix;
    }
  if (nbrKeptOrbitals == 0)
    {
      if (nbrParticleSector == 0)
	{
	  ComplexMatrix TmpEntanglementMatrix(1, this->HilbertSpaceDimension, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      TmpEntanglementMatrix[i][0] = groundState[i];
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  if (nbrKeptOrbitals == this->NbrSite)
    {
      if (nbrParticleSector == this->NbrBosons)
	{
	  ComplexMatrix TmpEntanglementMatrix(this->HilbertSpaceDimension, 1, true);
	  for (int i = 0; i < this->HilbertSpaceDimension; ++i)
	    {
	      TmpEntanglementMatrix[0][i] = groundState[i];
	    }
	  return TmpEntanglementMatrix;
	}
      else
	{
	  ComplexMatrix TmpEntanglementMatrix;
	  return TmpEntanglementMatrix;
	}
    }
  this->KeptOrbitals = new int [nbrKeptOrbitals];
  for (int i = 0 ; i < nbrKeptOrbitals; ++i) 
    this->KeptOrbitals[i] = keptOrbitals[i];
  BosonOnLatticeGutzwillerProjectionRealSpace SubsytemSpace (nbrParticleSector, nbrKeptOrbitals);
  BosonOnLatticeGutzwillerProjectionRealSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSite - nbrKeptOrbitals);
  ComplexMatrix TmpEntanglementMatrix (SubsytemSpace.GetHilbertSpaceDimension(), ComplementarySpace.HilbertSpaceDimension, true);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;

  long TmpEntanglementMatrixZero = this->EvaluatePartialEntanglementMatrixCore(0, ComplementarySpace.HilbertSpaceDimension, &ComplementarySpace, &SubsytemSpace, groundState, &TmpEntanglementMatrix);
//   FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpEntanglementMatrix);
//   Operation.ApplyOperation(architecture);
//   if (Operation.GetNbrNonZeroMatrixElements() > 0)	
  if (TmpEntanglementMatrixZero > 0)
     return TmpEntanglementMatrix;
   else
     {
       ComplexMatrix TmpEntanglementMatrixZero;
       return TmpEntanglementMatrixZero;
     }    
}
  
// core part of the evaluation orbital cut entanglement matrix calculation
// 
// minIndex = first index to consider in source Hilbert space
// nbrIndex = number of indices to consider in source Hilbert space
// complementaryHilbertSpace = pointer to the complementary Hilbert space (i.e. part B)
// destinationHilbertSpace = pointer to the destination Hilbert space  (i.e. part A)
// groundState = reference on the total system ground state
// densityMatrix = reference on the density matrix where result has to stored
// return value = number of components that have been added to the density matrix

long BosonOnLatticeGutzwillerProjectionRealSpace::EvaluatePartialEntanglementMatrixCore (int minIndex, int nbrIndex, ParticleOnSphere* complementaryHilbertSpace,  ParticleOnSphere* destinationHilbertSpace,
								       ComplexVector& groundState, ComplexMatrix* entanglementMatrix)
{
  BosonOnLatticeGutzwillerProjectionRealSpace* TmpHilbertSpace = (BosonOnLatticeGutzwillerProjectionRealSpace*) complementaryHilbertSpace;
  BosonOnLatticeGutzwillerProjectionRealSpace* TmpDestinationHilbertSpace = (BosonOnLatticeGutzwillerProjectionRealSpace*) destinationHilbertSpace;
  long TmpNbrNonZeroElements = 0;
  int* TraceOutOrbitals = new int [this->LzMax - TmpDestinationHilbertSpace->LzMax];
  int MaxIndex = minIndex + nbrIndex;
  int TmpIndex = 0;
  for (int i = 0; i < this->LzMax; ++i)
    {
      if (SearchInArray<int>(i, this->KeptOrbitals, TmpDestinationHilbertSpace->LzMax) < 0)
	TraceOutOrbitals[TmpIndex++] = i;
    }
  for (; minIndex < MaxIndex; ++minIndex)    
    {
      int Pos = 0;
      unsigned long TmpStateCompact = TmpHilbertSpace->StateDescription[minIndex];
      unsigned long TmpState = 0x0ul;
      for (int i = 0 ; i < TmpHilbertSpace->LzMax; ++i)
	TmpState |= (TmpStateCompact >> i) << TraceOutOrbitals[i];
      for (int j = 0; j < TmpDestinationHilbertSpace->HilbertSpaceDimension; ++j)
	{
	  unsigned long TmpStateCompact2 = TmpDestinationHilbertSpace->StateDescription[j];
	  unsigned long TmpState2 = 0x0ul;
	  for (int i = 0 ; i < TmpDestinationHilbertSpace->LzMax; ++i)
	    TmpState2 |= ((TmpStateCompact2 >> i) & 0x1ul) << this->KeptOrbitals[i];
	  unsigned long TmpState3 = TmpState | TmpState2;
	  int TmpLzMax = (this->LzMax << 1) + 1; 
	  while ((TmpState3 >> TmpLzMax) == 0x0ul)
	    --TmpLzMax;
	  int TmpPos = this->FindStateIndex(TmpState3, TmpLzMax);
	  if (TmpPos != this->HilbertSpaceDimension)
	    {
	      entanglementMatrix->AddToMatrixElement(j, minIndex,groundState[TmpPos]);
	      ++TmpNbrNonZeroElements;
	    }
	}
    }
  delete[] TraceOutOrbitals;
  return TmpNbrNonZeroElements;
}


