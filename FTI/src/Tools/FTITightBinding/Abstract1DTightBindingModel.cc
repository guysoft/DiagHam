////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 1D tight binding model                  //
//                                                                            //
//                        last modification : 01/05/2012                      //
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
#include "Tools/FTITightBinding/Abstract1DTightBindingModel.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"

#include <iostream>
#include <fstream>


using std::ofstream;
using std::endl;
using std::cout;


// default constructor
//

Abstract1DTightBindingModel::Abstract1DTightBindingModel()
{
  this->EnergyBandStructure = 0;
  this->OneBodyBasis = 0;
  this->EmbeddingX = RealVector();
}

// destructor
//

Abstract1DTightBindingModel::~Abstract1DTightBindingModel()
{
  if (this->OneBodyBasis != 0)
    {
      delete[] this->OneBodyBasis;
    }
  if (this->EnergyBandStructure != 0)
    {
      for (int i = 0; i < this->NbrBands; ++i)
	{
	  delete[] this->EnergyBandStructure[i];
	}
      delete[] this->EnergyBandStructure;
    }
}

// write an header that describes the tight binding model
// 
// output = reference on the output stream
// return value  = reference on the output stream

ofstream& Abstract1DTightBindingModel::WriteHeader(ofstream& output)
{
  int Dimension = 1;
  int HeaderSize = (((this->NbrBands + 2) * Dimension) * sizeof(double)) + ((Dimension + 1) * sizeof(int));
  WriteLittleEndian(output, HeaderSize);
  WriteLittleEndian(output, Dimension);
  WriteLittleEndian(output, this->NbrSiteX);
  WriteLittleEndian(output, this->KxFactor);
  WriteLittleEndian(output, this->GammaX);
  if (this->EmbeddingX.GetVectorDimension() != this->NbrBands)
  {
      double Tmp = 0.0;
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, Tmp);
  }
  else
  {
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, this->EmbeddingX[i]);
  }
  return output; 
}

// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool Abstract1DTightBindingModel::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# kx";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  File << endl;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      File << kx; 
      for (int i = 0; i < this->NbrBands; ++i)
	File << " " << this->EnergyBandStructure[i][kx];
      File << endl;
    }
  File.close();
  return true;
}

// flatten the energy spectrum
//
// flatteningFactor = flattening factor applie to each band (the band is rescale with respect to its average value)
// bandShift = shift each band by bandShift times the band index

void Abstract1DTightBindingModel::FlattenBands(double flatteningFactor, double bandShift)
{
  double* TmpAverageEnergies = new double[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    TmpAverageEnergies[i] = 0.0;
  for (int i = 0; i < this->NbrBands; ++i)
    {
      for (int j = 0; j < this->NbrStatePerBand; ++j)      
	TmpAverageEnergies[i] += this->EnergyBandStructure[i][j];
      TmpAverageEnergies[i] /= (double)  this->NbrStatePerBand;
    }
  for (int i = 0; i < this->NbrBands; ++i)
    {
      for (int j = 0; j < this->NbrStatePerBand; ++j)      
	{
	  this->EnergyBandStructure[i][j] = (((this->EnergyBandStructure[i][j] - TmpAverageEnergies[i]) * flatteningFactor)
					     + TmpAverageEnergies[i] + (((double) i) * bandShift));
	}
    }
  delete[] TmpAverageEnergies;
}

// get all the energies of a band, sorted from the smallest to the largest
//
// energies = reference to the array where the energies will be stored (the allocation is done by the method)
// momenta = reference to the array where the linearized momenta associated to each energy will be stored (the allocation is done by the method)
// bandIndex = index of the band to consider

void Abstract1DTightBindingModel::GetEnergies(double*& energies, int*& momenta, int bandIndex)
{
  energies = new double[this->NbrStatePerBand];
  momenta = new int[this->NbrStatePerBand];
  for (int i = 0; i < this->NbrStatePerBand; ++i)
    {
      energies[i] = this->EnergyBandStructure[bandIndex][i];
      momenta[i] = i;
    }
  SortArrayUpOrdering<int>(energies, momenta, this->NbrStatePerBand);
}

// get all the energies, sorted from the smallest to the largest
//
// energies = reference to the array where the energies will be stored (the allocation is done by the method)
// momenta = reference to the array where the linearized momentum associated to each energy will be stored (the allocation is done by the method)
// bandIndices = reference to the array where the band index associated to each energy will be stored (the allocation is done by the method)

void Abstract1DTightBindingModel::GetAllEnergies(double*& energies, int*& momenta, int*& bandIndices)
{
  int TotalNbrStates = this->NbrStatePerBand * this->NbrBands;
  energies = new double[TotalNbrStates];
  momenta = new int[TotalNbrStates];
  bandIndices = new int[TotalNbrStates];
  TotalNbrStates= 0;
  for (int i = 0; i < this->NbrStatePerBand; ++i)
    {
      for (int j = 0; j < this->NbrBands; ++j)
	{
	  energies[TotalNbrStates] = this->EnergyBandStructure[j][i];
	  momenta[TotalNbrStates] = i;
	  bandIndices[TotalNbrStates] = j;
	  ++TotalNbrStates;
	}
    }
  SortArrayUpOrdering<int>(energies, momenta, bandIndices, TotalNbrStates);
}

// evaluate the two point correlation function in a given region
//
// maxX = x coordinate of the region upper right corner 
// maxY = y coordinate of the region upper right corner 
// occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
// nbrOccupiedMomenta = number of occupied momenta
// bandIndex = index of the band to consider
// return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)

HermitianMatrix Abstract1DTightBindingModel::EvaluateFullTwoPointCorrelationFunction(int maxX, int maxY, int* occupiedMomenta, int nbrOccupiedMomenta, int bandIndex)
{
  int TotalNbrSites = maxX * this->NbrBands;
  int TmpMomentumX;
  HermitianMatrix EntanglementHamiltonian(TotalNbrSites, true);
  cout << "warning, Abstract1DTightBindingModel::EvaluateFullTwoPointCorrelationFunction is not implemented" << endl;
  return EntanglementHamiltonian;
}

// evaluate the mixed two point correlation function in a given region, assuming translation invariance along one direction
//
// maxX = length along the borken translation direction of the region 
// ky = momentum along the translation invariant direction
// occupiedMomenta = array that gives all the occupied momenta (as linearized indices)
// bandIndices = array that gives the band index of each occupied state
// nbrOccupiedMomenta = number of occupied momenta
// return value = matrix where the values of the two point correlation function will be stored (using the linearized position index as entry)

HermitianMatrix Abstract1DTightBindingModel::EvaluateFullMixedTwoPointCorrelationFunctionWithK(int maxX, int ky, int* occupiedMomenta, int* bandIndices, int nbrOccupiedMomenta)
{
  int TmpNbrOrbitalPerUnitCell = this->NbrBands;
  int TotalNbrSites = maxX * TmpNbrOrbitalPerUnitCell;
  int TmpMomentumX;
  int TmpMomentumY;
  HermitianMatrix EntanglementHamiltonian(TotalNbrSites, true);
  int TmpNbrStates = 0;
  for (int i = 0; i < nbrOccupiedMomenta; ++i)
    {
      if (occupiedMomenta[i] == ky)
	{
	  ++TmpNbrStates;
	}
    }
  if (TmpNbrStates == 0)
    {
      return EntanglementHamiltonian;
    }

  int* TmpOccupiedMomenta = new int [TmpNbrStates];
  TmpNbrStates = 0;
  for (int i = 0; i < nbrOccupiedMomenta; ++i)
    {
      if (occupiedMomenta[i] == ky)
	{
	  TmpOccupiedMomenta[TmpNbrStates] = bandIndices[i];
	  ++TmpNbrStates;
	}
    }
  for (int TmpY1 = 0; TmpY1 < maxX; ++TmpY1)
    {	  
      for (int TmpOrbital1 = 0; TmpOrbital1 < TmpNbrOrbitalPerUnitCell; ++TmpOrbital1)
	{
	  int TmpReducedLinearizedIndex1 = TmpOrbital1 + (TmpY1 * TmpNbrOrbitalPerUnitCell);
	  for (int TmpY2 = 0; TmpY2 < maxX; ++TmpY2)
	    {	  
	      for (int TmpOrbital2 = 0; TmpOrbital2 < TmpNbrOrbitalPerUnitCell; ++TmpOrbital2)
		{
		  int TmpReducedLinearizedIndex2 = TmpOrbital2 + (TmpY2 * TmpNbrOrbitalPerUnitCell);
		  Complex Tmp = 0.0;
		  for (int i = 0; i < TmpNbrStates; ++i)
		    {
		      Tmp +=  Conj(this->OneBodyBasis[ky][TmpOccupiedMomenta[i]][TmpReducedLinearizedIndex1]) * this->OneBodyBasis[ky][TmpOccupiedMomenta[i]][TmpReducedLinearizedIndex2];
		    }		  
		  EntanglementHamiltonian.SetMatrixElement(TmpReducedLinearizedIndex1, TmpReducedLinearizedIndex2, Tmp);
		}
	    }
	}
    }
  delete[] TmpOccupiedMomenta;

  return EntanglementHamiltonian;
}

