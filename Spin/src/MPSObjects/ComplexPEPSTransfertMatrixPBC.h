#ifndef _COMPLEXPEPSTRANSFERTMATRIXPBC_H
#define _COMPLEXPEPSTRANSFERTMATRIXPBC_H

#include "MPSObjects/AbstractTransfertMatrixPBC.h"
#include "MPSObjects/AbstractMPSSite.h"
#include "Tensor/Tensor3.h"
#include "HilbertSpace/AbstractDoubledSpinChain.h"
#include "Architecture/AbstractArchitecture.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

class AbstractMPSSite;

class ComplexPEPSTransfertMatrixPBC : public  AbstractTransfertMatrixPBC
{
 protected:
  
  Complex *** ValuesNonZeroTensorElementTopLeft;
  
 public:

  ComplexPEPSTransfertMatrixPBC ();
  ComplexPEPSTransfertMatrixPBC(MultiColumnASCIIFile & tensorElementsFile,AbstractArchitecture * architecture);
  ~ComplexPEPSTransfertMatrixPBC ();  
  
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  //  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent);
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent);

  // multiply a set of vectors by the current hamiltonian for a given range of indices 
  // and add result to another set of vectors, low level function (no architecture optimization)
  //
  // vSources = array of vectors to be multiplied
  // vDestinations = array of vectors at which result has to be added
  // nbrVectors = number of vectors that have to be evaluated together
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = pointer to the array of vectors where result has been stored
  virtual ComplexVector* LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent);
   
  virtual int GenerateResultingStateAndCoefficient(int indiceTop, int chainSize, int lastIndice, Complex * coefArray, unsigned long * stateArrayBra,  unsigned long * stateArrayKet, unsigned long pos);
  
  virtual void PrintTensorElements();
  
  inline  int GetBondDimension() const {return this->MPOBondDimension; }; 

 protected:
  
  inline void GetBraAndKetIndexFromCommonIndex(unsigned int communIndex,unsigned int & braIndex ,unsigned int & ketIndex);
  inline unsigned int GetCommonIndexFromBraAndKetIndices(unsigned int braIndex, unsigned int ketIndex );


  virtual void InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile);
};


inline void ComplexPEPSTransfertMatrixPBC::GetBraAndKetIndexFromCommonIndex(unsigned int communIndex,unsigned int & braIndex ,unsigned int & ketIndex)
{
  braIndex = communIndex%this->MPOBondDimension;
  ketIndex = communIndex/this->MPOBondDimension;
}


inline unsigned int  ComplexPEPSTransfertMatrixPBC::GetCommonIndexFromBraAndKetIndices(unsigned int braIndex, unsigned int ketIndex )
{
  return ketIndex  * this->MPOBondDimension + braIndex;
}

#endif

