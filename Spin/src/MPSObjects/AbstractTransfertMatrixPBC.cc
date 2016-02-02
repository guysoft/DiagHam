#include "AbstractTransfertMatrixPBC.h"
#include "Tensor/Tensor3.h"
#include "Matrix/RealSymmetricMatrix.h"
#include <iostream>
#include <sys/time.h>

using std::cout;
using std::endl;

AbstractTransfertMatrixPBC::AbstractTransfertMatrixPBC()
{
  this->NbrNonZeroElements = 0;
  this->PhysicalDimension = 0;
  this->MPOBondDimension = 0;
  this->IndiceBottomNonZeroTensorElementTopLeft = 0;
  this->IndiceRightNonZeroTensorElementTopLeft = 0;
  this->NbrNonZeroTensorElementTopLeft = 0;
  this->ValuesNonZeroTensorElementTopLeft = 0;
}

AbstractTransfertMatrixPBC::~AbstractTransfertMatrixPBC()
{
  for (int i = 0; i < this->MPOBondDimension; i++)
    {
      for(int j = 0; j < this->PhysicalDimension; j++)
	{
	  delete [] this->IndiceBottomNonZeroTensorElementTopLeft[i][j];
	  delete [] this->IndiceRightNonZeroTensorElementTopLeft[i][j];
	  delete [] this->ValuesNonZeroTensorElementTopLeft[i][j];
	}
      delete [] this->IndiceBottomNonZeroTensorElementTopLeft[i];
      delete [] this->IndiceRightNonZeroTensorElementTopLeft[i];
      delete [] this->ValuesNonZeroTensorElementTopLeft[i];
      delete [] this->NbrNonZeroTensorElementTopLeft[i];
    }
}


// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void AbstractTransfertMatrixPBC::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->HilbertSpace = (AbstractSpinChain * )hilbertSpace;
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* AbstractTransfertMatrixPBC::GetHilbertSpace ()
{
  return this->HilbertSpace;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractTransfertMatrixPBC::GetHilbertSpaceDimension ()
{
  return this->HilbertSpace->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void AbstractTransfertMatrixPBC::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}


// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored
RealVector& AbstractTransfertMatrixPBC::LowLevelMultiply(RealVector& vSource, RealVector& vDestination)
{ 
  return this->LowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}

// multiply a vector by the current hamiltonian and store result in another vector
// low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// return value = reference on vectorwhere result has been stored
ComplexVector& AbstractTransfertMatrixPBC::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination)
{ 
  return this->LowLevelMultiply(vSource, vDestination, 0, this->GetHilbertSpaceDimension());
}


RealVector& AbstractTransfertMatrixPBC::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent +  nbrComponent;
  for (int i = firstComponent ; i  <  LastComponent; i++)
    {
      int Dim = 0;
      int * IndiceLeft = this->HilbertSpace->GetBosonicOccupation(i);
      for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension;  IndiceTop++)
	{	  
	  Dim += this->EvaluateNbrResultingState(IndiceTop,IndiceLeft,this->HilbertSpace->GetSpinChainLength()-1,IndiceTop);
	}
      double * ResultingCoefficient = new double [Dim];
      unsigned long * ResultingIndex = new unsigned long[Dim];
      for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension;  IndiceTop++)
	{	  
	  this->GenerateResultingStateAndCoefficient(IndiceTop,IndiceLeft,this->HilbertSpace->GetSpinChainLength()-1,IndiceTop,ResultingCoefficient,ResultingIndex,0ul);
	}
      for (int p = 0; p < Dim; p++)
	{
	  vDestination[this->HilbertSpace->FindStateIndex(ResultingIndex[p])] = ResultingCoefficient[p]*vSource[i];
	}
    }
  return vDestination;
}


int AbstractTransfertMatrixPBC::EvaluateNbrResultingState(int indiceTop,int * indiceLeft,int chainSize, int lastIndice)
{
  int Tmp = 0;
  if(chainSize==0)
    {
      for (int i = 0; i < this->NbrNonZeroTensorElementTopLeft[indiceTop][indiceLeft[0]]; i++)
	{
	  if(this->IndiceBottomNonZeroTensorElementTopLeft[indiceTop][indiceLeft[0]][i] ==  lastIndice ) 
	    Tmp++;
	}
      return Tmp;
    }
  
  for (int i = 0; i < this->NbrNonZeroTensorElementTopLeft[indiceTop][indiceLeft[chainSize]]; i++)
    {
      Tmp+=this->EvaluateNbrResultingState(IndiceBottomNonZeroTensorElementTopLeft[indiceTop][indiceLeft[chainSize]][i],indiceLeft,chainSize-1,lastIndice);
    }
  return Tmp;
}


void AbstractTransfertMatrixPBC::InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile)
{
  this->NbrNonZeroElements = tensorElementsFile.GetNbrLines();
  int* IndexDown  = tensorElementsFile.GetAsIntegerArray (0);
  int* IndexUp = tensorElementsFile.GetAsIntegerArray (1);
  int* IndexLeft = tensorElementsFile.GetAsIntegerArray (2);
  int* IndexRight = tensorElementsFile.GetAsIntegerArray (3);
  double * ElementsValues = tensorElementsFile.GetAsDoubleArray (4);
  int TmpPhysicalDimension = 0;
  int TmpMPODimension = 0;
  unsigned long MemoryCost = 0;  
  cout <<"Nbr Non zero Elements = "<<  this->NbrNonZeroElements<<endl;
  
  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      if (IndexDown[i] > TmpPhysicalDimension)
	TmpPhysicalDimension = IndexDown[i];
      if (IndexLeft[i] > TmpMPODimension)
	TmpMPODimension = IndexLeft[i];
    }
  
  this->PhysicalDimension = TmpPhysicalDimension+1;
  this->MPOBondDimension =  TmpMPODimension+1;
  
  this->NbrNonZeroTensorElementTopLeft = new int * [this->MPOBondDimension];
  this->IndiceBottomNonZeroTensorElementTopLeft = new int ** [this->MPOBondDimension];
  this->IndiceRightNonZeroTensorElementTopLeft = new int ** [this->MPOBondDimension];
  MemoryCost+=sizeof(int*)*this->MPOBondDimension + 2*sizeof(int**)*this->MPOBondDimension;
  for (int i = 0; i < this->MPOBondDimension; i++)
    {
      this->NbrNonZeroTensorElementTopLeft[i] = new int [this->PhysicalDimension];
      this->IndiceBottomNonZeroTensorElementTopLeft[i] = new int * [this->PhysicalDimension];
      this->IndiceRightNonZeroTensorElementTopLeft[i] = new int * [this->PhysicalDimension];
      this->ValuesNonZeroTensorElementTopLeft[i] = new double * [this->PhysicalDimension];
      MemoryCost+=this->PhysicalDimension*(2*sizeof(int*)+sizeof(double*)+sizeof(int));
      for(int j = 0; j < this->PhysicalDimension; j++)
	{
	  this->NbrNonZeroTensorElementTopLeft[i][j]=0;
	}
    }
  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      this->NbrNonZeroTensorElementTopLeft[IndexUp[i]][IndexLeft[i]]++;
    } 
  
  for (int i = 0; i < this->MPOBondDimension; i++)
    {
      for(int j = 0; j < this->PhysicalDimension; j++)
	{
	  this->IndiceBottomNonZeroTensorElementTopLeft[i][j] = new int [this->NbrNonZeroTensorElementTopLeft[i][j]];
	  this->IndiceRightNonZeroTensorElementTopLeft[i][j] =  new int [this->NbrNonZeroTensorElementTopLeft[i][j]];
	  this->ValuesNonZeroTensorElementTopLeft[i][j] = new double [this->NbrNonZeroTensorElementTopLeft[i][j]];
	  this->NbrNonZeroTensorElementTopLeft[i][j]=0;
	  MemoryCost+= (2*sizeof(int)+ sizeof(double))*this->NbrNonZeroTensorElementTopLeft[i][j];
	}
    }
  
  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      this->IndiceBottomNonZeroTensorElementTopLeft[IndexUp[i]][IndexLeft[i]][this->NbrNonZeroTensorElementTopLeft[IndexUp[i]][IndexLeft[i]]] = IndexDown[i];
      this->IndiceRightNonZeroTensorElementTopLeft[IndexUp[i]][IndexLeft[i]][this->NbrNonZeroTensorElementTopLeft[IndexUp[i]][IndexLeft[i]]] = IndexRight[i];
      this->ValuesNonZeroTensorElementTopLeft[IndexUp[i]][IndexLeft[i]][this->NbrNonZeroTensorElementTopLeft[IndexUp[i]][IndexLeft[i]]] =  ElementsValues[i];
      this->NbrNonZeroTensorElementTopLeft[IndexUp[i]][IndexLeft[i]]++;
    }
  
  cout <<" Physical Dimension = " <<  this->PhysicalDimension<<endl;
  cout <<" MPO Dimension = " <<  this->MPOBondDimension <<endl;
  cout <<"Memory Cost = "<<MemoryCost <<endl;
  delete [] IndexDown;
  delete [] IndexUp;
  delete [] IndexLeft;
  delete [] IndexRight;
  delete [] ElementsValues;
}


int AbstractTransfertMatrixPBC::GenerateResultingStateAndCoefficient(int indiceTop, int * indiceLeft, int chainSize, int lastIndice, double * coefArray, unsigned long * stateArray, unsigned long pos)
{
  if(chainSize==0)
    {
      for (int i = 0; i < this->NbrNonZeroTensorElementTopLeft[indiceTop][indiceLeft[0]]; i++)
	{
	  if(this->IndiceBottomNonZeroTensorElementTopLeft[indiceTop][indiceLeft[0]][i] ==  lastIndice ) 
	    {
	      stateArray[pos] = this->HilbertSpace->EncodeSiteState(this->IndiceRightNonZeroTensorElementTopLeft[indiceTop][indiceLeft[0]][i],0);
	      coefArray[pos]  = this->ValuesNonZeroTensorElementTopLeft[indiceTop][indiceLeft[0]][i];
	      
	      pos++;
	    }
	}
      return pos;
    }
  int TmpPos = pos;
  
  for (int i = 0; i < this->NbrNonZeroTensorElementTopLeft[indiceTop][indiceLeft[chainSize]]; i++)
    {
      TmpPos +=this->GenerateResultingStateAndCoefficient(this->IndiceBottomNonZeroTensorElementTopLeft[indiceTop][indiceLeft[chainSize]][i],indiceLeft, chainSize-1, lastIndice, coefArray, stateArray, pos);
      for(;pos <TmpPos;pos++)
	{
	  stateArray[pos] |= this->HilbertSpace->EncodeSiteState(this->IndiceRightNonZeroTensorElementTopLeft[indiceTop][indiceLeft[chainSize]][i],chainSize);
	  coefArray[pos] *= this->ValuesNonZeroTensorElementTopLeft[indiceTop][indiceLeft[chainSize]][i];
	}
    }
  return pos;
}
