////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 10/13/2005                        //
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
#include "Hamiltonian/QuantumDotHamiltonian/QuantumWellHamiltonianInMagneticField.h"
 
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
// zEnergy1 = z confinement in the first subband
// zEnergy2 = z confinement in the second subband
// landauIndex1 = Landau index of the first subband
// landauIndex2 = Landau index of the second subband
// mailleParameter =
// bandOffset = conduction band offset between GaAs and InAs
// inDopage = In/Ga dopage ratio (=x with Ga_(1-x) In_x As)

QuantumWellHamiltonianInMagneticField::QuantumWellHamiltonianInMagneticField(double xSize, double ySize, double zSize, double mass, double bField, double zEnergy1, double zEnergy2,
									     int landauIndex1, int landauIndex2, double mailleParameter, double bandOffset, double inDopage)
{
  this->XSize = xSize;
  this->YSize = ySize;
  this->ZSize = zSize;
  this->Mass = mass;
  this->BField = bField;
  this->ZEnergy1 = zEnergy1;
  this->ZEnergy2 = zEnergy2;
  this->LandauIndex1 = landauIndex1;
  this->LandauIndex2 = landauIndex2;
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
  HermitianMatrix TmpHamiltonian(this->LandauDegeneracy * 2, true);
  this->Hamiltonian = TmpHamiltonian;
  this->EvaluateInteractionFactors();
}

// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

QuantumWellHamiltonianInMagneticField::QuantumWellHamiltonianInMagneticField(const QuantumWellHamiltonianInMagneticField& hamiltonian)
{
}

// destructor
//

QuantumWellHamiltonianInMagneticField::~QuantumWellHamiltonianInMagneticField()
{
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* QuantumWellHamiltonianInMagneticField::Clone ()
{
  return new QuantumWellHamiltonianInMagneticField(*this);
}


// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void QuantumWellHamiltonianInMagneticField::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void QuantumWellHamiltonianInMagneticField::ShiftHamiltonian (double shift)
{
}


// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex QuantumWellHamiltonianInMagneticField::MatrixElement (RealVector& V1, RealVector& V2)
{
  Complex Tmp;
  return Tmp;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex QuantumWellHamiltonianInMagneticField::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  Complex Tmp;
  return Tmp;
}

// store Hamiltonian into an hermitian matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix

HermitianMatrix& QuantumWellHamiltonianInMagneticField::GetHamiltonian (HermitianMatrix& M)
{
  M = this->Hamiltonian;
  return M;
}
  
// evaluate all interaction factors
//   

void QuantumWellHamiltonianInMagneticField::EvaluateInteractionFactors()
{
  int Lim = 2 * this->LandauDegeneracy;
  double DiagonalTerm1 = this->ZEnergy1 + ((0.5 + ((double) this->LandauIndex1)) * this->CyclotronEnergy);
  double DiagonalTerm2 = this->ZEnergy2 + ((0.5 + ((double) this->LandauIndex2)) * this->CyclotronEnergy);
  
  for (int i = 0; i < Lim; i += 2)
    {      
      this->Hamiltonian.SetMatrixElement(i, i, DiagonalTerm1);
      this->Hamiltonian.SetMatrixElement(i + 1, i + 1, DiagonalTerm2);
    }

  double XInc = this->MailleParameter * 0.5;
  double YInc = this->MailleParameter * 0.5;
  double ZInc = this->MailleParameter;
  int NbrXCells = (int)  (2.0 * this->XSize / this->MailleParameter);
  int NbrYCells = (int)  (2.0 * this->YSize / this->MailleParameter);
  int NbrZCells = (int)  (this->ZSize / this->MailleParameter);
  double X = 0.5 * XInc;
  double Y = 0.5 * YInc;
  double Z = 0.5 * ZInc;
  int ProbabilityThreshold = (int) (((double) RAND_MAX) * this->InDopage);
  double GaCoefficient = this->MailleParameter * this->MailleParameter * this->MailleParameter * 0.25 * (2.0 / (this->YSize * this->ZSize)) * this->BandOffset;
  double InCoefficient = GaCoefficient * (this->InDopage - 1.0);
  GaCoefficient *= this->InDopage;
  double Coefficient;
  double KCoeffcient = 2.0 * M_PI / this->YSize;
  double XCoeffcient = this->MagneticLength * this->MagneticLength * KCoeffcient;
  double LandauPrefactor = pow(M_PI * this->MagneticLength * this->MagneticLength, -0.25);  
  int Threshold = (int) (this->InDopage * ((double) RAND_MAX));
  for (int k = 0; k < NbrZCells; ++k)
    {
      X = 0.5 * XInc;
      double ZPartValue1 = sin (M_PI * Z / this->ZSize);
      double ZPartValue2 = sin (2.0 * M_PI * Z / this->ZSize);
      for (int i = 0; i < NbrXCells; ++i)
	{
	  Y = 0.5 * YInc;
	  for (int j = 0; j < NbrYCells; ++j)
	    {
	      if (rand() < Threshold)
		Coefficient = InCoefficient;
	      else
		Coefficient = GaCoefficient;
	      for (int m = 0; m < this->LandauDegeneracy; ++m)
		{
		  double ShiftXM = (X / this->MagneticLength) - (this->MagneticLength * KCoeffcient * m);
		  double Landau11 = LandauPrefactor * exp (-0.25 * (ShiftXM * ShiftXM));
		  double Landau21 = Landau11 * M_SQRT1_2 * ((2.0 * ShiftXM * ShiftXM) - 1.0);
		  for (int n = m + 1; n < this->LandauDegeneracy; ++n)
		    {	
		      double ShiftXN = (X / this->MagneticLength) - (this->MagneticLength * KCoeffcient * n);
		      double Landau12 = LandauPrefactor * exp (-0.25 * (ShiftXN * ShiftXN));
		      double Landau22 = Landau12 * M_SQRT1_2 * ((2.0 * ShiftXN * ShiftXN) -1.0);
		      Complex Tmp11 (cos(Y* KCoeffcient *((double) (n - m))), sin (Y* KCoeffcient *((double) (n - m))));
		      Complex Tmp22 (Tmp11);
		      Complex Tmp12 (Tmp11);
		      Complex Tmp21 (Tmp11);
		      Tmp11 *= Coefficient * ZPartValue1 * ZPartValue1 * Landau11 * Landau12;
		      Tmp21 *= Coefficient * ZPartValue1 * ZPartValue2 * Landau12 * Landau22;
		      Tmp12 *= Coefficient * ZPartValue1 * ZPartValue2 * Landau11 * Landau22;
		      Tmp22 *= Coefficient * ZPartValue2 * ZPartValue2 * Landau22 * Landau21;
		      this->Hamiltonian.AddToMatrixElement(2 * m, 2 * n, Tmp11);
		      this->Hamiltonian.AddToMatrixElement(2 * m, 2 * n + 1, Tmp12);
		      this->Hamiltonian.AddToMatrixElement(2 * m + 1, 2 * n, Tmp21);
		      this->Hamiltonian.AddToMatrixElement(2 * m + 1, 2 * n + 1, Tmp22);
		    }
		  this->Hamiltonian.AddToMatrixElement(2 * m, 2 * m, Coefficient * ZPartValue1 * ZPartValue1 * Landau11 * Landau11);
		  this->Hamiltonian.AddToMatrixElement(2 * m, 2 * m + 1, Coefficient * ZPartValue1 * ZPartValue2 * Landau11 * Landau21);
		  this->Hamiltonian.AddToMatrixElement(2 * m + 1, 2 * m + 1, Coefficient * ZPartValue2 * ZPartValue2 * Landau21 * Landau21);		  
		}
	      Y += YInc;	      
	    }
	  X += XInc;
	}
      Z += ZInc;
    }
}

