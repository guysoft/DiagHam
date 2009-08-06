////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2009 Nicolas Regnault                  //
//                         class author: Gunnar Möller                        //
//                                                                            //
//                                                                            //
//                       class of Hamiltonian H=S*S on lattice                //
//                                                                            //
//                        last modification : 31/07/2009                      //
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
#include "SUNSpinOnLatticeQuadraticHamiltonian.h"
#include "HilbertSpace/GenericSUNSpinCollection.h"
#include "Hamiltonian/AbstractSUNSpinOnLatticeHamiltonian.h"
#include "Tools/LatticeConnections.h"
#include "GeneralTools/StringTools.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/SUNSpinPrecalculationOperation.h"



// constructor
// space = Hilbert space for problem
// lattice = class providing size and geometry / connections on lattice
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)
// precalculationFileName = option file name where precalculation can be read instead of reevaluting them
//
SUNSpinOnLatticeQuadraticHamiltonian::SUNSpinOnLatticeQuadraticHamiltonian(GenericSUNSpinCollection *space, LatticeConnections *lattice, AbstractArchitecture* architecture, long memory, char* precalculationFileName)
{ 
  this->Spins = space;
  this->Lattice = lattice;
  this->NbrSpins = Lattice->GetNbrSites();
  this->FastMultiplicationFlag = false;
  this->HamiltonianShift = 0.0;
  this->Architecture = architecture;
  this->EvaluateInteractionTerms();
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;
  this->DiskStorageFlag = false;
  
  if (precalculationFileName == 0)
    {
      if (memory > 0)
	{
	  long TmpMemory = this->FastMultiplicationMemory(memory);
	  cout  << "fast = ";
	  PrintMemorySize(cout, TmpMemory)<<endl;
	  if (memory > 0)
	    {
	      this->EnableFastMultiplication();
	    }
	}
    }
  else
    this->LoadPrecalculation(precalculationFileName);

}

// destructor
//
SUNSpinOnLatticeQuadraticHamiltonian::~SUNSpinOnLatticeQuadraticHamiltonian()
{
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian
//SUNSpinOnLatticeQuadraticHamiltonian::SUNSpinOnLatticeQuadraticHamiltonian* Clone ()
//{
//}


// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int SUNSpinOnLatticeQuadraticHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Spins->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SUNSpinOnLatticeQuadraticHamiltonian::ShiftHamiltonian (double shift)
{
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex SUNSpinOnLatticeQuadraticHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Spins->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex SUNSpinOnLatticeQuadraticHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> SUNSpinOnLatticeQuadraticHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> SUNSpinOnLatticeQuadraticHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}




// evaluate all interaction factors
//   
void SUNSpinOnLatticeQuadraticHamiltonian::EvaluateInteractionTerms()
{
  this->HaveComplexInteractions = false;
  int *Partners;
  int NbrPartners;
  int TmpNbrPermutationTerms=0;
  for (int s=0; s<NbrSpins; ++s)
    {      
      this->Lattice->GetPartners(s, Partners, NbrPartners);
      TmpNbrPermutationTerms+=NbrPartners;
    }
  this->NbrPermutationTerms=TmpNbrPermutationTerms;
  this->PermutationPrefactors = new double[this->NbrPermutationTerms];
  this->PermutationI = new int[this->NbrPermutationTerms];
  this->PermutationJ = new int[this->NbrPermutationTerms];
  int Pos=0;
  for (int s=0; s<NbrSpins; ++s)
    {      
      this->Lattice->GetPartners(s, Partners, NbrPartners);
      for (int p=0; p<NbrPartners; ++p)
	{
	  this->PermutationI[Pos]=s;
	  this->PermutationJ[Pos]=Partners[p];
	  this->PermutationPrefactors[Pos]=1.0;
	  ++Pos;
	}
    }
}
