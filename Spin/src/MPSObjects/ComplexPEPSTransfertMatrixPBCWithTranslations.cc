#include "ComplexPEPSTransfertMatrixPBCWithTranslations.h"

#include <iostream>
#include <sys/time.h>

using std::cout;
using std::endl;

ComplexPEPSTransfertMatrixPBCWithTranslations::ComplexPEPSTransfertMatrixPBCWithTranslations ()
{
}

ComplexPEPSTransfertMatrixPBCWithTranslations::ComplexPEPSTransfertMatrixPBCWithTranslations(MultiColumnASCIIFile & tensorElementsFile,AbstractArchitecture * architecture): ComplexPEPSTransfertMatrixPBC (tensorElementsFile,architecture)
{
}

ComplexPEPSTransfertMatrixPBCWithTranslations::~ComplexPEPSTransfertMatrixPBCWithTranslations()
{
}

/*
ComplexVector& ComplexPEPSTransfertMatrixPBCWithTranslations::LowLevelAddMultiply(ComplexVector & vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  cout <<"inside ComplexVector& ComplexPEPSTransfertMatrixPBCWithTranslations::LowLevelAddMultiply(ComplexVector & vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)"<<endl;
  int LastComponent = firstComponent +  nbrComponent;
  int * IndiceLeftBra = new int[this->ChainLength];
  int * IndiceLeftKet = new int[this->ChainLength];

  for (int i = firstComponent ; i  <  LastComponent; i++)
    {
      int Dim = 0;
      ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->GetBosonicOccupation(i,IndiceLeftBra,IndiceLeftKet);
      for (int k =0; k < this->HilbertSpace->GetSpinChainLength();k++)
	{
	  this->TemporaryArray[k]= this->GetCommonIndexFromBraAndKetIndices(IndiceLeftBra[k],IndiceLeftKet[k]);
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
      int NbrTranslation;
      unsigned long TmpStateDescritionKet, TmpStateDescritionBra;
      for (int p = 0; p < Dim; p++)
	{
	  ((AbstractDoubledSpinChainWithTranslations * ) this->HilbertSpace)->FindCanonicalForm(ResultingIndexBra[p],ResultingIndexKet[p],TmpStateDescritionBra ,TmpStateDescritionKet, NbrTranslation);
	  
	  int Index = ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->FindStateIndex(TmpStateDescritionBra,TmpStateDescritionKet);
	  if ( Index < ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->GetHilbertSpaceDimension() )
	    {
	      vDestination[Index] += ResultingCoefficient[p]*vSource[i]* this->ExponentialFactors[NbrTranslation] * ((AbstractDoubledSpinChainWithTranslations *) this->HilbertSpace)->GetRescalingFactor(i,Index); 
	    }
	}
      delete [] ResultingCoefficient; delete [] ResultingIndexBra;  delete [] ResultingIndexKet; 
    }
  delete [] IndiceLeftBra;
  delete [] IndiceLeftKet;
  return vDestination;
}
*/



void ComplexPEPSTransfertMatrixPBCWithTranslations::SetHilbertSpace(AbstractHilbertSpace * hilbertSpace)
{
  this->HilbertSpace = (AbstractSpinChain *) hilbertSpace;
  this->XMomentum= ((AbstractDoubledSpinChainWithTranslations *) this->HilbertSpace)->GetMomentum();
  
  if (  this->ChainLength != this->HilbertSpace->GetSpinChainLength() )
    {
      delete this->TmpVector1; 
      delete this->TmpVector2;
      delete this->EndVector;
      delete this->StartVector;
      this->ChainLength = this->HilbertSpace->GetSpinChainLength();
      this->MaxXMomentum= this->HilbertSpace->GetSpinChainLength();
      
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
  unsigned long MemoryCost =  (2*this->PowerD[this->ChainLength+1] + 2*this->PowerD[this->ChainLength])*sizeof(Complex);
  cout <<"Memory Cost " <<MemoryCost<<endl;
}

/*
ComplexVector* ComplexPEPSTransfertMatrixPBCWithTranslations::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent +  nbrComponent;
  int * IndiceLeftBra = new int[this->ChainLength];
  int * IndiceLeftKet = new int[this->ChainLength];

  for (int i = firstComponent ; i  <  LastComponent; i++)
    {
      int Dim = 0;
      ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->GetBosonicOccupation(i,IndiceLeftBra,IndiceLeftKet);
      for (int k =0; k < this->ChainLength ; k++)
	{
	  this->TemporaryArray[k]= this->GetCommonIndexFromBraAndKetIndices(IndiceLeftBra[k],IndiceLeftKet[k]);
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
      int NbrTranslation;
      unsigned long TmpStateDescritionKet, TmpStateDescritionBra;
      Complex TmpCoef;
      for (int p = 0; p < Dim; p++)
	{
	  ((AbstractDoubledSpinChainWithTranslations * ) this->HilbertSpace)->FindCanonicalForm(ResultingIndexBra[p],ResultingIndexKet[p],TmpStateDescritionBra ,TmpStateDescritionKet, NbrTranslation);
	  
	  int Index = ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->FindStateIndex(TmpStateDescritionBra,TmpStateDescritionKet);
	  if ( Index < ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->GetHilbertSpaceDimension() )
	    {
	      TmpCoef =  0;//ResultingCoefficient[p]*this->ExponentialFactors[NbrTranslation] * ((AbstractDoubledSpinChainWithTranslations *) this->HilbertSpace)->GetRescalingFactor(i,Index);
	      for(int t=0; t < nbrVectors; t++)
		vDestinations[t][Index] += TmpCoef *vSources[t][i];
	    }
	}
      delete [] ResultingCoefficient; delete [] ResultingIndexBra;  delete [] ResultingIndexKet; 
    }
  delete [] IndiceLeftBra;
  delete [] IndiceLeftKet;
  return vDestinations;
}
*/


ComplexVector& ComplexPEPSTransfertMatrixPBCWithTranslations::LowLevelAddMultiply(ComplexVector & vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  ComplexVector * TmpPointorVector;
  this->EndVector->ClearVector();
  ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->ConvertToGeneralSpaceWithMomentum(vSource,(*this->StartVector));
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
  ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->AddConvertFromGeneralSpaceWithMomentum((*this->EndVector),vDestination);
  return vDestination;
}
