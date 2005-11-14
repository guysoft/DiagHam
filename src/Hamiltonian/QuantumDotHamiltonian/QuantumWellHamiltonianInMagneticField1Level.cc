////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum well in mangetic field      //
//             resricted to one subband and one Landau level                  //
//                                                                            //
//                      last modification : 11/13/2005                        //
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
#include "Hamiltonian/QuantumDotHamiltonian/QuantumWellHamiltonianInMagneticField1Level.h"
 
#include <math.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


#define HBARE_M0 0.115767635
#define HBAR 1.05457168e-34
#define HBAR_E 6.582119138e-16
#define ECHARGE 1.60217653e-19
#define M0 9.1093826e-31


// constructor from default data
//
// xSize = system dimension in the x direction (in Angstrom unit)
// ySize = system dimension in the y direction (in Angstrom unit)
// zSize = system dimension in the z direction (in Angstrom unit)
// mass = effective mass in the x direction (in electron mass unit)
// bField = B field value (in Tesla)
// zEnergy = z confinement
// landauIndex = Landau level index
// mailleParameter =
// bandOffset = conduction band offset between GaAs and InAs
// inDopage = In/Ga dopage ratio (=x with Ga_(1-x) In_x As)
// potentialDescription = name of the file that contains the potential description (null if the potential has to be evaluated)

QuantumWellHamiltonianInMagneticField1Level::QuantumWellHamiltonianInMagneticField1Level(double xSize, double ySize, double zSize, double mass, double bField, double zEnergy,
											 int landauIndex, double mailleParameter, double bandOffset, double inDopage, 
											 char* potentialDescription)
{
  this->XSize = xSize;
  this->YSize = ySize;
  this->ZSize = zSize;
  this->Mass = mass;
  this->BField = bField;
  this->ZEnergy = zEnergy;
  this->LandauIndex = landauIndex;
  this->MailleParameter = mailleParameter;
  this->BandOffset = bandOffset;
  this->InDopage = inDopage;

  this->MagneticLength = 1.0e10 * sqrt(HBAR_E / this->BField);
  cout << this->MagneticLength << endl;
  this->CyclotronEnergy = HBARE_M0 * this->BField / this->Mass;
  cout << this->CyclotronEnergy << endl;
  this->LandauDegeneracy = (int) ((this->XSize * this->YSize) / (2.0 * M_PI * this->MagneticLength * this->MagneticLength));
  cout << this->LandauDegeneracy << endl;
  this->NbrCells = (int) ((4.0 * this->XSize * this->YSize * this->ZSize) / (this->MailleParameter * this->MailleParameter * this->MailleParameter));
  cout << this->NbrCells << endl;
  this->GaXPosition = new double [this->NbrCells];
  this->GaYPosition = new double [this->NbrCells];
  this->GaZPosition = new double [this->NbrCells];
  this->InXPosition = new double [this->NbrCells];
  this->InYPosition = new double [this->NbrCells];
  this->InZPosition = new double [this->NbrCells];
  HermitianMatrix TmpHamiltonian(this->LandauDegeneracy, true);
  this->Hamiltonian = TmpHamiltonian;
  this->NbrXCells = (int)  (2.0 * this->XSize / this->MailleParameter);
  this->NbrYCells = (int)  (2.0 * this->YSize / this->MailleParameter);
  this->NbrZCells = (int)  (this->ZSize / this->MailleParameter);
  this->Potential = new BinaryThreeDConstantCellPotential(NbrXCells, NbrYCells, NbrZCells);
  if (potentialDescription == 0)
    {
      int Threshold = (int) (this->InDopage * ((double) RAND_MAX));
      double GaCoefficient = this->MailleParameter * this->MailleParameter * this->MailleParameter * 0.25 * (2.0 / (this->YSize * this->ZSize)) * this->BandOffset;
      double InCoefficient = GaCoefficient * (this->InDopage - 1.0);
      GaCoefficient *= this->InDopage;
      for (int k = 0; k < NbrZCells; ++k)
	for (int j = 0; j < NbrYCells; ++j)
	  for (int i = 0; i < NbrXCells; ++i)
	    {
	      if (rand() < Threshold)
		this->Potential->SetPotential(i, j, k, InCoefficient);
	      else
		this->Potential->SetPotential(i, j, k, GaCoefficient);
	    }
    }
  else
    {
      cout << "check " << endl;
      this->Potential->LoadBinaryPotential(potentialDescription);
    }
  this->EvaluateInteractionFactors();
}

// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

QuantumWellHamiltonianInMagneticField1Level::QuantumWellHamiltonianInMagneticField1Level(const QuantumWellHamiltonianInMagneticField1Level& hamiltonian)
{
}

// destructor
//

QuantumWellHamiltonianInMagneticField1Level::~QuantumWellHamiltonianInMagneticField1Level()
{
  delete this->Potential;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* QuantumWellHamiltonianInMagneticField1Level::Clone ()
{
  return new QuantumWellHamiltonianInMagneticField1Level(*this);
}


// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void QuantumWellHamiltonianInMagneticField1Level::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void QuantumWellHamiltonianInMagneticField1Level::ShiftHamiltonian (double shift)
{
}


// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex QuantumWellHamiltonianInMagneticField1Level::MatrixElement (RealVector& V1, RealVector& V2)
{
  Complex Tmp;
  return Tmp;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex QuantumWellHamiltonianInMagneticField1Level::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  Complex Tmp;
  return Tmp;
}

// store Hamiltonian into an hermitian matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix

HermitianMatrix& QuantumWellHamiltonianInMagneticField1Level::GetHamiltonian (HermitianMatrix& M)
{
  M = this->Hamiltonian;
  return M;
}
  
// evaluate all interaction factors
//   

void QuantumWellHamiltonianInMagneticField1Level::EvaluateInteractionFactors()
{
  double DiagonalTerm1 = this->ZEnergy + ((0.5 + ((double) this->LandauIndex)) * this->CyclotronEnergy);
  
  for (int i = 0; i < this->LandauDegeneracy; ++i)
    {      
      this->Hamiltonian.SetMatrixElement(i, i, DiagonalTerm1);
    }

  double XInc = this->MailleParameter * 0.5;
  double YInc = this->MailleParameter * 0.5;
  double ZInc = this->MailleParameter;
  double X = 0.5 * XInc;
  double Y = 0.5 * YInc;
  double Z = 0.5 * ZInc;
  double Coefficient;
  double KCoeffcient = 2.0 * M_PI / this->YSize;
  double XCoeffcient = this->MagneticLength * this->MagneticLength * KCoeffcient;
  double LandauPrefactor = pow(M_PI * this->MagneticLength * this->MagneticLength, -0.25);  
  for (int k = 0; k < this->NbrZCells; ++k)
    {
      Y = 0.5 * YInc;
      double ZPartValue = sin (M_PI * Z / this->ZSize);
      for (int j = 0; j < this->NbrYCells; ++j)
	{
	  X = 0.5 * XInc;
	  for (int i = 0; i < this->NbrXCells; ++i)
	    {
	      Coefficient = this->Potential->GetPotential(i, j, k);
	      for (int m = 0; m < this->LandauDegeneracy; ++m)
		{
		  double ShiftXM = (X / this->MagneticLength) - (this->MagneticLength * KCoeffcient * m);
		  double Landau11 = LandauPrefactor * exp (-0.25 * (ShiftXM * ShiftXM));
		  for (int n = m + 1; n < this->LandauDegeneracy; ++n)
		    {	
		      double ShiftXN = (X / this->MagneticLength) - (this->MagneticLength * KCoeffcient * n);
		      double Landau12 = LandauPrefactor * exp (-0.25 * (ShiftXN * ShiftXN));
		      Complex Tmp11 (cos(Y* KCoeffcient *((double) (n - m))), -sin (Y* KCoeffcient *((double) (n - m))));
		      Tmp11 *= Coefficient * ZPartValue * ZPartValue * Landau11 * Landau12;
		      this->Hamiltonian.AddToMatrixElement(m, n, Tmp11);
		    }
		  this->Hamiltonian.AddToMatrixElement(m, m, Coefficient * ZPartValue * ZPartValue * Landau11 * Landau11);
		}
	      X += XInc;
	    }
	  Y += YInc;	      
	}
      Z += ZInc;
    }
}

// save potential on disk
// 
// filename = name of the file (with path) where potential has to be saved
// return value = true if no error occured

bool QuantumWellHamiltonianInMagneticField1Level::SavePotential(char* filename)
{
  this->Potential->SaveBinaryPotential(filename);
  return true;
}
