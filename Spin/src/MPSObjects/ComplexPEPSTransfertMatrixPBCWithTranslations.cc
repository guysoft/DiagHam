#include "ComplexPEPSTransfertMatrixPBCWithTranslations.h"

#include <iostream>
#include <sys/time.h>

using std::cout;
using std::endl;

ComplexPEPSTransfertMatrixPBCWithTranslations::ComplexPEPSTransfertMatrixPBCWithTranslations ()
{
  this->ExponentialFactors=0;
}

ComplexPEPSTransfertMatrixPBCWithTranslations::ComplexPEPSTransfertMatrixPBCWithTranslations(MultiColumnASCIIFile & tensorElementsFile,AbstractArchitecture * architecture): ComplexPEPSTransfertMatrixPBC (tensorElementsFile,architecture)
{
  this->ExponentialFactors=0;
}

ComplexPEPSTransfertMatrixPBCWithTranslations::~ComplexPEPSTransfertMatrixPBCWithTranslations()
{
  delete [] this->ExponentialFactors;
}


ComplexVector& ComplexPEPSTransfertMatrixPBCWithTranslations::LowLevelAddMultiply(ComplexVector & vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
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



// evaluate all exponential factors
//   

void ComplexPEPSTransfertMatrixPBCWithTranslations::EvaluateExponentialFactors()
{
  this->ExponentialFactors = new Complex[this->MaxXMomentum];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    { 
      this->ExponentialFactors[i] = Phase(2.0 * M_PI * ((this->XMomentum * ((double) i) / ((double) this->MaxXMomentum))));
    }
}

void ComplexPEPSTransfertMatrixPBCWithTranslations::SetHilbertSpace(AbstractHilbertSpace * hilbertSpace)
{
  this->HilbertSpace = (AbstractSpinChain *) hilbertSpace;
  this->XMomentum= ((AbstractDoubledSpinChainWithTranslations *) this->HilbertSpace)->GetMomentum();
  this->ChainLength = this->HilbertSpace->GetSpinChainLength();
  this->MaxXMomentum= this->HilbertSpace->GetSpinChainLength();
  delete [] this->ExponentialFactors;
  delete [] this->TemporaryArray;
  this->TemporaryArray = new int[this->ChainLength];
  this->EvaluateExponentialFactors();
}


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
	      TmpCoef =  ResultingCoefficient[p]*this->ExponentialFactors[NbrTranslation] * ((AbstractDoubledSpinChainWithTranslations *) this->HilbertSpace)->GetRescalingFactor(i,Index);
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
