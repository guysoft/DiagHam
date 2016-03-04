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

/*
inline unsigned long ComplexPEPSTransfertMatrixPBCWithTranslations::GetNewIndexFromOldIndex(unsigned long oldIndex, int oldPhysicalSpin, int newPhysicalSpin, int oldVirtualSpin, int newVirtualSpin, int position)
{ 
  return  oldIndex + (oldPhysicalSpin -  newPhysicalSpin) *this->PowerD[position] + newVirtualSpin - oldVirtualSpin;
}


int ComplexPEPSTransfertMatrixPBCWithTranslations::LowLevelAddMultiplyOnFirstSite(ComplexVector & vSource, unsigned long * oldHilbertSpace, unsigned long * newHilbertSpace, Complex * newWeigth, int oldHilbertSpaceDimension, int topIndice)
{
  int NbrElement =0;
  for(int i = 0; i< oldHilbertSpaceDimension;i++)
    {
      for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[topIndice][oldHilbertSpaceElement[i]%this->PowedD[1]]; p++)
	{
	  NbrElement+=SearchInArrayAndSetWeight(this->PowedD[1]*this->GetNewIndexFromOldIndex(oldHilbertSpace[i], oldHilbertSpace[i]%this->PowedD[1], this->IndiceRightNonZeroTensorElementTopLeft[indiceTop][oldHilbertSpace[i]%this->PowedD[1]][p],0,this->IndiceBottomNonZeroTensorElementTopLeft[indiceTop][oldHilbertSpace[i]%this->PowedD[1]][p],1),newHilbertSpace,newWeigth, NbrElement,this->ValuesNonZeroTensorElementTopLeft[indiceTop][oldHilbertSpace[i]%this->PowedD[1]][p]* vSource[i]);
	}
    }
  return NbrElement; 
}



int ComplexPEPSTransfertMatrixPBCWithTranslations::LowLevelAddMultiplyOnAnySite(unsigned long * oldHilbertSpace, Complex * oldWeigth, unsigned long * newHilbertSpace, Complex * newWeigth, int oldHilbertSpaceDimension, int position)
{
  int NbrElement =0;
  for(int i = 0; i< oldHilbertSpaceDimension;i++)
    {
      for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[oldHilbertSpaceElement[i]%this->PowedD[1]][(oldHilbertSpaceElement[i]/this->PowedD[position+1])%this->PowedD[1]] ; p++)
	{
	  NbrElement+=SearchInArrayAndSetWeight(this->GetNewIndexFromOldIndex(oldHilbertSpace[i], (oldHilbertSpaceElement[i]/this->PowedD[position+1])%this->PowedD[1], this->IndiceRightNonZeroTensorElementTopLeft[oldHilbertSpaceElement[i]%this->PowedD[1]][(oldHilbertSpaceElement[i]/this->PowedD[position+1])%this->PowedD[1]][p], oldHilbertSpaceElement[i]%this->PowedD[1], this->IndiceBottomNonZeroTensorElementTopLeft[oldHilbertSpaceElement[i]%this->PowedD[1]][(oldHilbertSpaceElement[i]/this->PowedD[position+1])%this->PowedD[1]][p],position),newHilbertSpace,newWeigth, NbrElement,this->ValuesNonZeroTensorElementTopLeft[oldHilbertSpaceElement[i]%this->PowedD[1]][(oldHilbertSpaceElement[i]/this->PowedD[position+1])%this->PowedD[1]][p] *  oldWeigth[i]);
	}
    }
  return NbrElement;
}


ComplexVector & ComplexPEPSTransfertMatrixPBCWithTranslations::LowLevelAddMultiplyOnLastSite(unsigned long * oldHilbertSpace, Complex * oldWeigth, ComplexVector & vDestination, int oldHilbertSpaceDimension, int topValue)
{

  int NbrTranslation;
  for(int i = 0; i< oldHilbertSpaceDimension;i++)
    {
      for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[oldHilbertSpaceElement[i]%this->PowedD[1]][oldHilbertSpaceElement[i]/this->PowedD[this->ChainLength]] ; p++)
	{
	  if (this->IndiceBottomNonZeroTensorElementTopLeft[[oldHilbertSpaceElement[i]%this->PowedD[1]][oldHilbertSpaceElement[i]/this->PowedD[this->ChainLength]][p] == topValue)
	    {
	      NbrTranslation = 0;
	      int Index = ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->FindStateIndexFromLinearizedIndexAndNumberTranslation(this->GetNewIndexFromOldIndex(oldHilbertSpace[i],oldHilbertSpaceElement[i]/this->PowedD[this->ChainLength], this->IndiceRightNonZeroTensorElementTopLeft[oldHilbertSpaceElement[i]%this->PowedD[1]][oldHilbertSpaceElement[i]/this->PowedD[this->ChainLength]], oldHilbertSpaceElement[i]%this->PowedD[1],0,this->ChainLength)/this->PowedD[1]),NbrTranslation);
	      
	      if ( Index < ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->GetHilbertSpaceDimension() )
		{
		  vDestination[Index] +=  oldWeigth[i] * this->ValuesNonZeroTensorElementTopLeft[oldHilbertSpaceElement[i]%this->PowedD[1]][oldHilbertSpaceElement[i]/this->PowedD[this->ChainLength])][p]*this->ExponentialFactors[NbrTranslation] *  ((AbstractDoubledSpinChainWithTranslations *) this->HilbertSpace)->GetRescalingFactor(i,Index);
	      }
	}
    }
  return vDestination;
}



ComplexVector& ComplexPEPSTransfertMatrixPBCWithTranslations::LowLevelAddMultiply(ComplexVector & vSource, ComplexVector& vDestination)
{
  int TmpDimension = 1;

  for(int i=0; i <=this->ChainLength;i++)
    TmpDimension *= this->MPOBondDimension * this->MPOBondDimension;

  unsigned long * OldHilbertSpace = new unsigned long [  ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->GetHilbertSpaceDimension()];
  ((AbstractDoubledSpinChain * )this->HilbertSpace)->GetChainDescriptionInCondensedForm(OldHilbertSpace);
  Complex * Weigth1 = new Complex [TmpDimension];
  Complex * Weigth2 = new Complex [TmpDimension];
  unsigned long * NewHilbertSpace1 = new Complex [TmpDimension];
  unsigned long * NewHilbertSpace2 = new Complex [TmpDimension];
  Complex * TmpPointorWeigth;
  unsigned long * TmpPointorHilbertSpace;
  int WorkingDimension = 0;
  for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension;  IndiceTop++)
    {	  
      WorkingDimension = this->LowLevelAddMultiplyOnFirstSite(vSource,OldHilbertSpace, NewHilbertSpace1, Weigth1, ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->GetHilbertSpaceDimension(),  IndiceTop);
      for(int Position = 1; Position <this->ChainLength - 1;Position++)
	{
	  WorkingDimension = this->LowLevelAddMultiplyOnAnySite(NewHilbertSpace1, Weigth1,  NewHilbertSpace2, Weigth2, WorkingDimension, Position);
	  TmpPointorHilbertSpace = NewHilbertSpace2;
	  NewHilbertSpace2 = NewHilbertSpace1;
	  NewHilbertSpace1 = TmpPointorHilbertSpace;
	  TmpPointorWeigth = Weigth1;
	  Weigth1 =  Weigth2;
	  Weigth2 =  TmpPointorWeigth;
	}
      this->LowLevelAddMultiplyOnLastSite(NewHilbertSpace1, Weigth1, vDestination, WorkingDimension, IndiceTop);
    }
  
  delete [ ] Weigth1, Weigth2,NewHilbertSpace1,NewHilbertSpace2;
  return vDestination;
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
