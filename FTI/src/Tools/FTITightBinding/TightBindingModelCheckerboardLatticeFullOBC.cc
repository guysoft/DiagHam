////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the checkerboard lattice       //
//                      with full open boundary conditions                    //
//                                                                            //
//                        last modification : 17/02/2018                      //
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
#include "Tools/FTITightBinding/TightBindingModelCheckerboardLatticeFullOBC.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include <iostream>

using std::cout;
using std::endl;
using std::ostream;



// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// mus = sublattice chemical potential on A sites
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
// blochFormFlag = use the Bloch form instead of the the traditional form

TightBindingModelCheckerboardLatticeFullOBC::TightBindingModelCheckerboardLatticeFullOBC(int nbrSiteX, int nbrSiteY, double t1, double t2, double mus, 
											 AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->MuS = mus;
  this->NbrBands = this->NbrSiteX * this->NbrSiteY;
  this->NbrStatePerBand = 1;
  this->Architecture = architecture;
  
  this->EnergyBandStructure = new double[this->NbrBands];

  this->FindConnectedOrbitals();
  this->ComputeBandStructure();
}



// destructor
//

TightBindingModelCheckerboardLatticeFullOBC::~TightBindingModelCheckerboardLatticeFullOBC()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelCheckerboardLatticeFullOBC::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  this->CoreComputeBandStructure(minStateIndex, nbrStates);
}

// find the orbitals connected to those located at the origin unit cell
// 
  
void TightBindingModelCheckerboardLatticeFullOBC::FindConnectedOrbitals()
{
  Complex TmpPhaseFactor = this->NNHopping * Phase(M_PI * 0.25);
  Complex TmpConjPhaseFactor = this->NNHopping * Phase(-M_PI * 0.25);

  this->NbrConnectedOrbitals = new int [this->NbrBands];
  this->ConnectedOrbitalIndices = new int* [this->NbrBands];
  this->ConnectedOrbitalSpatialIndices = new int* [this->NbrBands];
  this->ConnectedOrbitalHoppingAmplitudes = new Complex* [this->NbrBands];
  for (int x = 0; x < this->NbrSiteX; ++x)
    {
      for (int y = 0; y < this->NbrSiteY; ++x)
	{
	  int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndex(x, y);
	  this->NbrConnectedOrbitals[TmpLinearizedCoordinate] = 0;
	}
    }
  for (int x = 0; x < this->NbrSiteX; ++x)
    {
      for (int y = 0; y < this->NbrSiteY; ++x)
	{
	  int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndex(x, y);
	  this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
	  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x + 1, y + 1) >= 0)
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
	    }
	  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x - 1, y + 1) >= 0)
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
	    }
	  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x, y + 1) >= 0)
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
	    }
	  if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x + 1, y) >= 0)
	    {
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
	    }
	  this->ConnectedOrbitalIndices[TmpLinearizedCoordinate] = new int[this->NbrConnectedOrbitals[TmpLinearizedCoordinate]];	  
	  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate] = new int[this->NbrConnectedOrbitals[TmpLinearizedCoordinate]];
	  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate] = new Complex[this->NbrConnectedOrbitals[TmpLinearizedCoordinate]];	  
	}
    }
  for (int x = 0; x < this->NbrSiteX; ++x)
    {
      for (int y = 0; y < this->NbrSiteY; ++x)
	{
	  if (((x + y) & 1) == 0)
	    {
	      int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndex(x, y);
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinate] = 0;
	      this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpLinearizedCoordinate;
	      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
	      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->MuS;
	      ++this->NbrConnectedOrbitals[TmpLinearizedCoordinate];
	      if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x + 1, y + 1) >= 0)
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x + 1, y + 1);
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->NextNNHopping;
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		}
	      if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x - 1, y + 1) >= 0)
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x - 1, y + 1);
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = -this->NextNNHopping;
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		}
	      if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x, y + 1) >= 0)
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x, y + 1);
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpPhaseFactor;
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		}
	      if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x + 1, y) >= 0)
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x + 1, y);
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpConjPhaseFactor;
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		}
	    }
	  else
	    {
	      int TmpLinearizedCoordinate = this->GetRealSpaceTightBindingLinearizedIndex(x, y);
	      this->NbrConnectedOrbitals[TmpLinearizedCoordinate] = 0;
	      this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpLinearizedCoordinate;
	      this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
	      this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = -this->MuS;
	      ++this->NbrConnectedOrbitals[TmpLinearizedCoordinate];
	      if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x + 1, y + 1) >= 0)
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x + 1, y + 1);
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = -this->NextNNHopping;
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		}
	      if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x - 1, y + 1) >= 0)
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x - 1, y + 1);
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->NextNNHopping ;
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		}
	      if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x, y + 1) >= 0)
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x, y + 1);
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpConjPhaseFactor;
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		}
	      if (this->GetRealSpaceTightBindingLinearizedIndexSafe(x + 1, y) >= 0)
		{
		  this->ConnectedOrbitalIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = this->GetRealSpaceTightBindingLinearizedIndex(x + 1, y);
		  this->ConnectedOrbitalSpatialIndices[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = 0;
		  this->ConnectedOrbitalHoppingAmplitudes[TmpLinearizedCoordinate][this->NbrConnectedOrbitals[TmpLinearizedCoordinate]] = TmpPhaseFactor;
		  this->NbrConnectedOrbitals[TmpLinearizedCoordinate]++;
		}
	    }
	}
    }
}

