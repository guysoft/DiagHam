#include "ComplexPEPSTransfertMatrixPBC.h"

#include <iostream>
#include <sys/time.h>
#include "GeneralTools/ArrayTools.h"

#include "Matrix/RealDiagonalMatrix.h"


using std::cout;
using std::endl;

ComplexPEPSTransfertMatrixPBC::ComplexPEPSTransfertMatrixPBC ()
{
  this->BoundaryMatrix = 0; 
}

/*
ComplexPEPSTransfertMatrixPBC::ComplexPEPSTransfertMatrixPBC(MultiColumnASCIIFile & tensorElementsFile, AbstractArchitecture * architecture)
{
  this->InitializeTensorsElements(tensorElementsFile);
  this->Architecture = architecture;
  this->TemporaryArray = 0;
  this->Weigth1 = 0;
  this->Weigth2 = 0;
  this->NewHilbertSpace1 = 0;
  this->NewHilbertSpace2 = 0;
  this->OldHilbertSpace = 0;
}
*/


ComplexPEPSTransfertMatrixPBC::ComplexPEPSTransfertMatrixPBC(MultiColumnASCIIFile & tensorElementsFile, RealDiagonalMatrix * boundaryMatrix, AbstractArchitecture * architecture)
{
  this->InitializeTensorsElements(tensorElementsFile);
  this->Architecture = architecture;
  this->TemporaryArray = 0;
  this->TmpVector1 = 0;
  this->TmpVector2 = 0;
  this->StartVector = 0;
  this->EndVector = 0;
  this->ChainLength = 0;
  this->BoundaryMatrix = boundaryMatrix; 
}


ComplexPEPSTransfertMatrixPBC::ComplexPEPSTransfertMatrixPBC(MultiColumnASCIIFile & tensorElementsFile, AbstractArchitecture * architecture)
{
  this->InitializeTensorsElements(tensorElementsFile);
  this->Architecture = architecture;
  this->TemporaryArray = 0;
  this->TmpVector1 = 0;
  this->TmpVector2 = 0;
  this->StartVector = 0;
  this->EndVector = 0;
  this->ChainLength = 0;
  this->BoundaryMatrix = 0; 
}



ComplexPEPSTransfertMatrixPBC::~ComplexPEPSTransfertMatrixPBC()
{
  delete this->BoundaryMatrix;
  for (int i = 0; i < this->MPOBondDimension *  this->MPOBondDimension; i++)
    {
      for(int j = 0; j < this->PhysicalDimension * this->PhysicalDimension; j++)
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
  delete [] this->IndiceBottomNonZeroTensorElementTopLeft;
  delete [] this->IndiceRightNonZeroTensorElementTopLeft;
  delete [] this->ValuesNonZeroTensorElementTopLeft;
  delete [] this->NbrNonZeroTensorElementTopLeft;
  delete this->TmpVector1; delete  this->TmpVector2; delete  this->StartVector; delete this->EndVector;
  this->NbrNonZeroTensorElementTopLeft = 0;
}


void ComplexPEPSTransfertMatrixPBC::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{

  this->HilbertSpace = (AbstractSpinChain * )hilbertSpace;
  if ( this->ChainLength != this->HilbertSpace ->GetSpinChainLength() )
    {
      this->ChainLength =  this->HilbertSpace ->GetSpinChainLength();
  
      delete this->TmpVector1,   this->TmpVector2,   this->StartVector,   this->EndVector;
      delete [] this->PowerD;

      this->PowerD = new int[this->ChainLength+2];
      this->PowerD[0] = 1;
      for(int i =1; i <=this->ChainLength+1; i++)
	this->PowerD[i] = this->PowerD[i-1] * (this->MPOBondDimension * this->MPOBondDimension);

      this->TmpVector1 = new ComplexVector (this->PowerD[this->ChainLength+1],true);
      this->TmpVector2 = new ComplexVector (this->PowerD[this->ChainLength+1],true);
      this->EndVector =  new ComplexVector (this->PowerD[this->ChainLength],true);
      this->StartVector =  new ComplexVector (this->PowerD[this->ChainLength],true);
    }
  if(  this->BoundaryMatrix == 0)
    {
      this->BoundaryMatrix = new RealDiagonalMatrix( this->PowerD[1] ,true);
      this->BoundaryMatrix->SetToIdentity();
    }
}


/*
ComplexVector& ComplexPEPSTransfertMatrixPBC::LowLevelAddMultiply(ComplexVector & vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent +  nbrComponent;
  int * IndiceLeftBra = new int[this->ChainLength];
  int * IndiceLeftKet = new int[this->ChainLength];
  for (int i = firstComponent ; i  <  LastComponent; i++)
    {
      int Dim = 0;
      ((AbstractDoubledSpinChain * )this->HilbertSpace)->GetBosonicOccupation(i,IndiceLeftBra,IndiceLeftKet);
      for (int k =0; k < this->ChainLength;k++)
	{
	  this->TemporaryArray[k] = this->GetCommonIndexFromBraAndKetIndices(IndiceLeftBra[k],IndiceLeftKet[k]);
	}

      for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension;  IndiceTop++)
	{	  
	  Dim += this->EvaluateNbrResultingState(IndiceTop,this->ChainLength-1,IndiceTop);
	}

      Complex * ResultingCoefficient = new Complex [Dim];
      unsigned long * ResultingIndexBra = new unsigned long[Dim];
      unsigned long * ResultingIndexKet = new unsigned long[Dim];
      unsigned long Tmp=0;
      for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension ;  IndiceTop++)
	{	  
	  Tmp=this->GenerateResultingStateAndCoefficient(IndiceTop,this->ChainLength-1,IndiceTop,ResultingCoefficient,ResultingIndexBra,ResultingIndexKet,Tmp);
	}

      for (int p = 0; p < Dim; p++)
	{
	  vDestination[((AbstractDoubledSpinChain * )this->HilbertSpace)->FindStateIndex(ResultingIndexBra[p],ResultingIndexKet[p])] += ResultingCoefficient[p]*vSource[i];
	}
      delete [] ResultingCoefficient; delete [] ResultingIndexBra;  delete [] ResultingIndexKet; 
    }
  delete [] IndiceLeftBra;
  delete [] IndiceLeftKet;
  return vDestination;
}
*/

void ComplexPEPSTransfertMatrixPBC::InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile)
{
  this->NbrNonZeroElements = tensorElementsFile.GetNbrLines();
  int* IndexVertical  = tensorElementsFile.GetAsIntegerArray (0);
  int* IndexLeft = tensorElementsFile.GetAsIntegerArray (1);
  int* IndexUp = tensorElementsFile.GetAsIntegerArray (2);
  int* IndexRight = tensorElementsFile.GetAsIntegerArray (3);
  int* IndexDown  = tensorElementsFile.GetAsIntegerArray (4);

  Complex * ElementsValues = tensorElementsFile.GetAsComplexArray (5);
  int TmpPhysicalDimension = 0;
  int TmpMPODimension = 0;
  int VerticalDimension = 0;
  unsigned long MemoryCost = 0;  
  cout <<"Nbr Non zero Elements = "<<  this->NbrNonZeroElements<<endl;
  
  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      if (IndexVertical[i] >  VerticalDimension)
	{
	  VerticalDimension = IndexVertical[i];
	}
      if (IndexDown[i] > TmpPhysicalDimension)
	TmpPhysicalDimension = IndexDown[i];
      if (IndexLeft[i] > TmpMPODimension)
	TmpMPODimension = IndexLeft[i];
    }
  
  this->PhysicalDimension = TmpPhysicalDimension+1;
  this->MPOBondDimension =  TmpMPODimension+1;
  VerticalDimension++;
  
  this->NbrNonZeroTensorElementTopLeft = new int * [this->MPOBondDimension*this->MPOBondDimension];
  this->IndiceBottomNonZeroTensorElementTopLeft = new int ** [this->MPOBondDimension*this->MPOBondDimension];
  this->IndiceRightNonZeroTensorElementTopLeft = new int ** [this->MPOBondDimension*this->MPOBondDimension];
  this->ValuesNonZeroTensorElementTopLeft = new Complex ** [this->MPOBondDimension*this->MPOBondDimension];
  
  MemoryCost+= (sizeof(int*) + sizeof(Complex**) + 2*sizeof(int**))*this->MPOBondDimension*this->MPOBondDimension;

  for (int i = 0; i < this->MPOBondDimension * this->MPOBondDimension; i++)
    {
      this->NbrNonZeroTensorElementTopLeft[i] = new int [this->PhysicalDimension*this->PhysicalDimension];
      this->IndiceBottomNonZeroTensorElementTopLeft[i] = new int * [this->PhysicalDimension*this->PhysicalDimension];
      this->IndiceRightNonZeroTensorElementTopLeft[i] = new int * [this->PhysicalDimension*this->PhysicalDimension];
      this->ValuesNonZeroTensorElementTopLeft[i] = new Complex * [this->PhysicalDimension*this->PhysicalDimension];
      MemoryCost+=this->PhysicalDimension*this->PhysicalDimension*(2*sizeof(int*)+sizeof(double*)+sizeof(int));
      for(int j = 0; j < this->PhysicalDimension * this->PhysicalDimension; j++)
	{
	  this->NbrNonZeroTensorElementTopLeft[i][j] = 0;
	}
    }

  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      for(int j = 0 ; j < this->NbrNonZeroElements; j++)
	{
	  if (IndexVertical[i]==IndexVertical[j]) 
	    this->NbrNonZeroTensorElementTopLeft[GetCommonIndexFromBraAndKetIndices(IndexUp[i],IndexUp[j])][GetCommonIndexFromBraAndKetIndices(IndexLeft[i],IndexLeft[j])]++;
	} 
    }
  
  for (int i = 0; i < this->MPOBondDimension * this->MPOBondDimension; i++)
    {
      for(int j = 0; j < this->PhysicalDimension * this->PhysicalDimension; j++)
	{
	  this->IndiceBottomNonZeroTensorElementTopLeft[i][j] = new int [this->NbrNonZeroTensorElementTopLeft[i][j]];
	  this->IndiceRightNonZeroTensorElementTopLeft[i][j] =  new int [this->NbrNonZeroTensorElementTopLeft[i][j]];
	  this->ValuesNonZeroTensorElementTopLeft[i][j] = new Complex [this->NbrNonZeroTensorElementTopLeft[i][j]];
	  MemoryCost+= (2*sizeof(int)+ sizeof(Complex))*this->NbrNonZeroTensorElementTopLeft[i][j];
	  this->NbrNonZeroTensorElementTopLeft[i][j]=0;
	}
    }
  
  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      for(int j = 0 ; j < this->NbrNonZeroElements; j++)
	{
	  if (IndexVertical[i]==IndexVertical[j]) 
	    {
	      int IndiceUp = GetCommonIndexFromBraAndKetIndices(IndexUp[i],IndexUp[j]);
	      int IndiceLeft = GetCommonIndexFromBraAndKetIndices(IndexLeft[i],IndexLeft[j]);
	      this->IndiceBottomNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft][this->NbrNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft]] = GetCommonIndexFromBraAndKetIndices(IndexDown[i],IndexDown[j]);
	      this->IndiceRightNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft][this->NbrNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft]] = GetCommonIndexFromBraAndKetIndices(IndexRight[i],IndexRight[j]);
	      this->ValuesNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft][this->NbrNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft]] =  ElementsValues[i]*Conj(ElementsValues[j]);
	      this->NbrNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft]++;
	    }
	}
    }
  
  cout <<" Physical Dimension = " <<  this->PhysicalDimension<<endl;
  cout <<" MPO Dimension = " <<  this->MPOBondDimension <<endl;
  cout <<"Memory Cost = "<<MemoryCost <<endl;
  delete [] IndexVertical;
  delete [] IndexDown;
  delete [] IndexUp;
  delete [] IndexLeft;
  delete [] IndexRight;
  delete [] ElementsValues;
}

/*
int ComplexPEPSTransfertMatrixPBC::GenerateResultingStateAndCoefficient(int indiceTop, int chainSize, int lastIndice, Complex * coefArray, unsigned long * stateArrayBra, unsigned long * stateArrayKet, unsigned long pos)
{
  unsigned int IndexBra, IndexKet;
  if(chainSize==0)
    { 
      for (int i = 0; i < this->NbrNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[0]]; i++)
	{
	  if(this->IndiceBottomNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[0]][i] ==  lastIndice ) 
	    {
	      this->GetBraAndKetIndexFromCommonIndex(this->IndiceRightNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[0]][i], IndexBra, IndexKet);
	      stateArrayBra[pos] = ((AbstractDoubledSpinChain * )this->HilbertSpace)->EncodeSiteStateBra(IndexBra,0);
	      stateArrayKet[pos] = ((AbstractDoubledSpinChain * )this->HilbertSpace)->EncodeSiteStateKet(IndexKet,0);
	      coefArray[pos]  = this->ValuesNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[0]][i];
	      pos++;
	    }
	}
      return pos;
    }
  int TmpPos = pos;
  
  for (int i = 0; i < this->NbrNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[chainSize]]; i++)
    {
      TmpPos =this->GenerateResultingStateAndCoefficient(this->IndiceBottomNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[chainSize]][i], chainSize-1, lastIndice, coefArray, stateArrayBra, stateArrayKet, pos);
      for(;pos <TmpPos;++pos)
	{
	  this->GetBraAndKetIndexFromCommonIndex(this->IndiceRightNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[chainSize]][i], IndexBra, IndexKet);
	  stateArrayBra[pos] |= ((AbstractDoubledSpinChain * )this->HilbertSpace)->EncodeSiteStateBra(IndexBra,chainSize);
	  stateArrayKet[pos] |= ((AbstractDoubledSpinChain * )this->HilbertSpace)->EncodeSiteStateKet(IndexKet,chainSize);
	  coefArray[pos] *= this->ValuesNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[chainSize]][i];
	}
    }
  return pos;
}
*/
/*

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* ComplexPEPSTransfertMatrixPBC::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent +  nbrComponent;
  int * IndiceLeftBra = new int[this->ChainLength];
  int * IndiceLeftKet = new int[this->ChainLength];
  for (int i = firstComponent ; i  <  LastComponent; i++)
    {
      int Dim = 0;
      ((AbstractDoubledSpinChain * )this->HilbertSpace)->GetBosonicOccupation(i,IndiceLeftBra,IndiceLeftKet);
      for (int k =0; k < this->ChainLength;k++)
	{
	  this->TemporaryArray[k] = this->GetCommonIndexFromBraAndKetIndices(IndiceLeftBra[k],IndiceLeftKet[k]);
	}

      for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension;  IndiceTop++)
	{	  
	  Dim += this->EvaluateNbrResultingState(IndiceTop,this->ChainLength-1,IndiceTop);
	}

      Complex * ResultingCoefficient = new Complex [Dim];
      unsigned long * ResultingIndexBra = new unsigned long[Dim];
      unsigned long * ResultingIndexKet = new unsigned long[Dim];
      unsigned long Tmp=0;
      for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension ;  IndiceTop++)
	{	  
	  Tmp=this->GenerateResultingStateAndCoefficient(IndiceTop,this->ChainLength-1,IndiceTop,ResultingCoefficient,ResultingIndexBra,ResultingIndexKet,Tmp);
	}

      for (int p = 0; p < Dim; p++)
	{
	  int Index = ((AbstractDoubledSpinChain * )this->HilbertSpace)->FindStateIndex(ResultingIndexBra[p],ResultingIndexKet[p]);
	  for(int t=0; t < nbrVectors; t++)
	    vDestinations[t][Index] += ResultingCoefficient[p]*vSources[t][i];
	}
      delete [] ResultingCoefficient; delete [] ResultingIndexBra;  delete [] ResultingIndexKet; 
    }
  delete [] IndiceLeftBra;
  delete [] IndiceLeftKet;
  return vDestinations;  
}

*/

void  ComplexPEPSTransfertMatrixPBC::PrintTensorElements()
{
  cout <<"#Tensor IndiceLeft IndiceTop  IndicexRight IndiceBottom Values" <<endl;
  for(int IndiceLeft=0; IndiceLeft < this->PhysicalDimension * this->PhysicalDimension; IndiceLeft++)
    {
      for(int IndiceTop =0; IndiceTop < this->MPOBondDimension *this->MPOBondDimension ; IndiceTop++)
	{
	  for (int i = 0; i < this->NbrNonZeroTensorElementTopLeft[IndiceTop][IndiceLeft]; i++)
	    {
	      cout <<IndiceLeft<< " "<< IndiceTop  <<" "<< this->IndiceRightNonZeroTensorElementTopLeft[IndiceTop][IndiceLeft][i]<< " "<< this->IndiceBottomNonZeroTensorElementTopLeft[IndiceTop][IndiceLeft][i]<<" "<<this->ValuesNonZeroTensorElementTopLeft[IndiceTop][IndiceLeft][i]<<endl;
	    }
	}
    }
}


/*
int ComplexPEPSTransfertMatrixPBC::LowLevelAddMultiplyOnFirstSite(ComplexVector & vSource, unsigned long * oldHilbertSpace, unsigned long * newHilbertSpace, Complex * newWeigth, int oldHilbertSpaceDimension, int topIndice)
{
  int NbrElement =0;
  for(int i = 0; i< oldHilbertSpaceDimension;i++)
    {
      for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[topIndice][oldHilbertSpace[i]%this->PowerD[1]]; p++)
	{
	  NbrElement+=SearchInArrayAndSetWeight(this->GetNewIndexFromOldIndex(this->PowerD[1]*oldHilbertSpace[i], oldHilbertSpace[i]%this->PowerD[1], this->IndiceRightNonZeroTensorElementTopLeft[topIndice][oldHilbertSpace[i]%this->PowerD[1]][p],0,this->IndiceBottomNonZeroTensorElementTopLeft[topIndice][oldHilbertSpace[i]%this->PowerD[1]][p],1) ,newHilbertSpace,newWeigth, NbrElement,this->ValuesNonZeroTensorElementTopLeft[topIndice][oldHilbertSpace[i]%this->PowerD[1]][p]* vSource[i]);
	}
    }
  return NbrElement; 
}



int ComplexPEPSTransfertMatrixPBC::LowLevelAddMultiplyOnAnySite(unsigned long * oldHilbertSpace, Complex * oldWeigth, unsigned long * newHilbertSpace, Complex * newWeigth, int oldHilbertSpaceDimension, int position)
{
  int NbrElement =0;
  for(int i = 0; i< oldHilbertSpaceDimension;i++)
    {
      for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][(oldHilbertSpace[i]/this->PowerD[position+1])%this->PowerD[1]] ; p++)
	{
	  NbrElement+=SearchInArrayAndSetWeight(this->GetNewIndexFromOldIndex(oldHilbertSpace[i], (oldHilbertSpace[i]/this->PowerD[position+1])%this->PowerD[1], this->IndiceRightNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][(oldHilbertSpace[i]/this->PowerD[position+1])%this->PowerD[1]][p], oldHilbertSpace[i]%this->PowerD[1], this->IndiceBottomNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][(oldHilbertSpace[i]/this->PowerD[position+1])%this->PowerD[1]][p],position+1)

	    ,newHilbertSpace,newWeigth, NbrElement,this->ValuesNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][(oldHilbertSpace[i]/this->PowerD[position+1])%this->PowerD[1]][p] *  oldWeigth[i]);

	}
    }
  return NbrElement;
}



ComplexVector & ComplexPEPSTransfertMatrixPBC::LowLevelAddMultiplyOnLastSite(unsigned long * oldHilbertSpace, Complex * oldWeigth, ComplexVector & vDestination, int oldHilbertSpaceDimension, int topValue)
{
for(int i = 0; i< oldHilbertSpaceDimension;i++)
    {
    for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][oldHilbertSpace[i]/this->PowerD[this->ChainLength]] ; p++)
	{
	  if (this->IndiceBottomNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][oldHilbertSpace[i]/this->PowerD[this->ChainLength]][p] == topValue)
	    {
	      int Index = ((AbstractDoubledSpinChain * )this->HilbertSpace)->FindStateIndexFromLinearizedIndex(this->GetNewIndexFromOldIndex(oldHilbertSpace[i],oldHilbertSpace[i]/this->PowerD[this->ChainLength], this->IndiceRightNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][oldHilbertSpace[i]/this->PowerD[this->ChainLength]][p], oldHilbertSpace[i]%this->PowerD[1],0,this->ChainLength)/this->PowerD[1]);
	      
	      if ( Index < ((AbstractDoubledSpinChain * )this->HilbertSpace)->GetHilbertSpaceDimension() )
		{
		  vDestination[Index] +=  oldWeigth[i] * this->ValuesNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][oldHilbertSpace[i]/this->PowerD[this->ChainLength]][p];
		}
	    }
	}
    }
  return vDestination;
}


ComplexVector& ComplexPEPSTransfertMatrixPBC::LowLevelAddMultiply(ComplexVector & vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  Complex * TmpPointorWeigth;
  unsigned long * TmpPointorHilbertSpace;
  int WorkingDimension = 0;
  for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension;  IndiceTop++)
    {	  
      
      WorkingDimension = this->LowLevelAddMultiplyOnFirstSite(vSource,this->OldHilbertSpace, this->NewHilbertSpace1, this->Weigth1, ((AbstractDoubledSpinChain * )this->HilbertSpace)->GetHilbertSpaceDimension(),  IndiceTop);
      for(int Position = 1; Position <this->ChainLength - 1;Position++)
	{
	  WorkingDimension = this->LowLevelAddMultiplyOnAnySite(NewHilbertSpace1, this->Weigth1,   this->NewHilbertSpace2, this->Weigth2, WorkingDimension, Position);
	  TmpPointorHilbertSpace = this->NewHilbertSpace2;
	  this->NewHilbertSpace2 =  this->NewHilbertSpace1;
	  this->NewHilbertSpace1 = TmpPointorHilbertSpace;
	  TmpPointorWeigth =  this->Weigth1;
	  this->Weigth1 = this->Weigth2;
	  this->Weigth2 = TmpPointorWeigth;
	}
      this->LowLevelAddMultiplyOnLastSite(this->NewHilbertSpace1, this->Weigth1, vDestination, WorkingDimension, IndiceTop);
    }
  
  return vDestination;
}




int ComplexPEPSTransfertMatrixPBC::LowLevelAddMultiplyOnFirstSite(ComplexVector & vSource, unsigned long * oldHilbertSpace, unsigned long * newHilbertSpace, Complex * newWeigth, int oldHilbertSpaceDimension, int topIndice)
{
  int NbrElement =0;
  for(int i = 0; i< oldHilbertSpaceDimension;i++)
    {
      for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[topIndice][oldHilbertSpace[i]%this->PowerD[1]]; p++)
	{
	  NbrElement+=SearchInArrayAndSetWeight(this->GetNewIndexFromOldIndex(this->PowerD[1]*oldHilbertSpace[i], oldHilbertSpace[i]%this->PowerD[1], this->IndiceRightNonZeroTensorElementTopLeft[topIndice][oldHilbertSpace[i]%this->PowerD[1]][p],0,this->IndiceBottomNonZeroTensorElementTopLeft[topIndice][oldHilbertSpace[i]%this->PowerD[1]][p],1) ,newHilbertSpace,newWeigth, NbrElement,this->ValuesNonZeroTensorElementTopLeft[topIndice][oldHilbertSpace[i]%this->PowerD[1]][p]* vSource[i]);
	}
    }
  return NbrElement; 
}



int ComplexPEPSTransfertMatrixPBC::LowLevelAddMultiplyOnAnySite(unsigned long * oldHilbertSpace, Complex * oldWeigth, unsigned long * newHilbertSpace, Complex * newWeigth, int oldHilbertSpaceDimension, int position)
{
  int NbrElement =0;
  for(int i = 0; i< oldHilbertSpaceDimension;i++)
    {
      for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][(oldHilbertSpace[i]/this->PowerD[position+1])%this->PowerD[1]] ; p++)
	{
	  NbrElement+=SearchInArrayAndSetWeight(this->GetNewIndexFromOldIndex(oldHilbertSpace[i], (oldHilbertSpace[i]/this->PowerD[position+1])%this->PowerD[1], this->IndiceRightNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][(oldHilbertSpace[i]/this->PowerD[position+1])%this->PowerD[1]][p], oldHilbertSpace[i]%this->PowerD[1], this->IndiceBottomNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][(oldHilbertSpace[i]/this->PowerD[position+1])%this->PowerD[1]][p],position+1)

,newHilbertSpace,newWeigth, NbrElement,this->ValuesNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][(oldHilbertSpace[i]/this->PowerD[position+1])%this->PowerD[1]][p] *  oldWeigth[i]);

	}
    }
  return NbrElement;
}


ComplexVector & ComplexPEPSTransfertMatrixPBC::LowLevelAddMultiplyOnLastSite(unsigned long * oldHilbertSpace, Complex * oldWeigth, ComplexVector & vDestination, int oldHilbertSpaceDimension, int topValue)
{
  for(int i = 0; i< oldHilbertSpaceDimension;i++)
    {
      for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][oldHilbertSpace[i]/this->PowerD[this->ChainLength]] ; p++)
	{
	  if (this->IndiceBottomNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][oldHilbertSpace[i]/this->PowerD[this->ChainLength]][p] == topValue)
	    {
	      int Index = ((AbstractDoubledSpinChain * )this->HilbertSpace)->FindStateIndexFromLinearizedIndex(this->GetNewIndexFromOldIndex(oldHilbertSpace[i],oldHilbertSpace[i]/this->PowerD[this->ChainLength], this->IndiceRightNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][oldHilbertSpace[i]/this->PowerD[this->ChainLength]][p], oldHilbertSpace[i]%this->PowerD[1],0,this->ChainLength)/this->PowerD[1]);
	      
	      if ( Index < ((AbstractDoubledSpinChain * )this->HilbertSpace)->GetHilbertSpaceDimension() )
		{
		  vDestination[Index] +=  oldWeigth[i] * this->ValuesNonZeroTensorElementTopLeft[oldHilbertSpace[i]%this->PowerD[1]][oldHilbertSpace[i]/this->PowerD[this->ChainLength]][p];
		}
	    }
	}
    }
  return vDestination;
}
*/

ComplexVector& ComplexPEPSTransfertMatrixPBC::LowLevelAddMultiply(ComplexVector & vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  ComplexVector * TmpPointorVector;
  this->EndVector->ClearVector();
  ((AbstractDoubledSpinChain * )this->HilbertSpace)->ConvertToGeneralSpace(vSource,(*this->StartVector));
  for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension;  IndiceTop++)
    {	  
      this->LowLevelAddMultiplyOnFirstSite(IndiceTop);
      for(int Position = 1; Position <this->ChainLength - 1;Position++)
	{
	  this->LowLevelAddMultiplyOnAnySite(Position);
	  TmpPointorVector = this->TmpVector1;
	  this->TmpVector1 = this->TmpVector2;
	  this->TmpVector2 = TmpPointorVector;
	  this->TmpVector2->ClearVector();
	}
      this->LowLevelAddMultiplyOnLastSite(IndiceTop);
      this->TmpVector1->ClearVector();
    }
  ((AbstractDoubledSpinChain * )this->HilbertSpace)->AddConvertFromGeneralSpace((*this->EndVector),vDestination);
  return vDestination;
}

void ComplexPEPSTransfertMatrixPBC::LowLevelAddMultiplyOnFirstSite(int topIndice)
{
  for(int i = 0; i< this->StartVector->GetVectorDimension();i++) 
    {
      if ( Norm((*this->StartVector)[i]) > 1e-13 )
	{
	  for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[topIndice][i%this->PowerD[1]]; p++)
	    {
	      (*this->TmpVector1)[this->GetNewIndexFromOldIndex(this->PowerD[1]*i,i%this->PowerD[1], this->IndiceRightNonZeroTensorElementTopLeft[topIndice][i%this->PowerD[1]][p],0,this->IndiceBottomNonZeroTensorElementTopLeft[topIndice][i%this->PowerD[1]][p],1)]+=this->ValuesNonZeroTensorElementTopLeft[topIndice][i%this->PowerD[1]][p]*  (*this->StartVector)[i];
	    }
	}
    }
}


void ComplexPEPSTransfertMatrixPBC::LowLevelAddMultiplyOnAnySite(int position)
{
  for(int i = 0; i< this->TmpVector1->GetVectorDimension();i++) 
    {
      if ( Norm((*this->TmpVector1)[i]) > 1e-13 )
	{
	  for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[i%this->PowerD[1]][(i/this->PowerD[position+1])%this->PowerD[1]] ; p++)
	    {
	      (*this->TmpVector2)[this->GetNewIndexFromOldIndex(i, (i/this->PowerD[position+1])%this->PowerD[1], this->IndiceRightNonZeroTensorElementTopLeft[i%this->PowerD[1]][(i/this->PowerD[position+1])%this->PowerD[1]][p], i%this->PowerD[1], this->IndiceBottomNonZeroTensorElementTopLeft[i%this->PowerD[1]][(i/this->PowerD[position+1])%this->PowerD[1]][p],position+1)] += this->ValuesNonZeroTensorElementTopLeft[i%this->PowerD[1]][(i/this->PowerD[position+1])%this->PowerD[1]][p] *  (*this->TmpVector1)[i];
	    }
	}
    }
}
  

void ComplexPEPSTransfertMatrixPBC::LowLevelAddMultiplyOnLastSite(int topValue)
{
  double BoundaryValue;
  this->BoundaryMatrix->GetMatrixElement(topValue,topValue,BoundaryValue);
  for(int i = 0; i< this->TmpVector1->GetVectorDimension();i++) 
    {
      if ( Norm((*this->TmpVector1)[i]) > 1e-13 )
	{
	  for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[i%this->PowerD[1]][i/this->PowerD[this->ChainLength]] ; p++)
	    {
	      if (this->IndiceBottomNonZeroTensorElementTopLeft[i%this->PowerD[1]][i/this->PowerD[this->ChainLength]][p] == topValue)
		{
		  (*this->EndVector)[this->GetNewIndexFromOldIndex(i,i/this->PowerD[this->ChainLength], this->IndiceRightNonZeroTensorElementTopLeft[i%this->PowerD[1]][i/this->PowerD[this->ChainLength]][p],i%this->PowerD[1],0,this->ChainLength)/this->PowerD[1]] +=  (*this->TmpVector1)[i] * this->ValuesNonZeroTensorElementTopLeft[i%this->PowerD[1]][i/this->PowerD[this->ChainLength]][p] * BoundaryValue;
		  
		}
	    }
	}
    }
}


