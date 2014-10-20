#ifndef _ABSTRACTMPOPERATOROBC_H
#define _ABSTRACTMPOPERATOROBC_H

#include "Hamiltonian/AbstractHamiltonian.h"
#include "MPSSite.h"
#include "Tensor/Tensor3.h"
#include "HilbertSpace/AbstractHilbertSpace.h"

class MPSSite;

class AbstractMPOperatorOBC : public AbstractHamiltonian
{
 protected:
  unsigned int NbrNonZeroElements;
  double * ElementsValues;
  unsigned int * IndexValues;
  unsigned int PhysicalDimension;
  unsigned int MPOBondDimension;
  unsigned int NbrSites;
  AbstractHilbertSpace* HilbertSpace;
  MPSSite * Site;
  
  AbstractMPOperatorOBC ();
  ~AbstractMPOperatorOBC ();
  
  virtual void InitializeTensorsElements() = 0;

  // global shift to apply to the diagonal matrix elements
  double HamiltonianShift;

 
 public:
  
  virtual void ComputeL(Tensor3<double> & L);
  
  virtual void ComputeR(Tensor3<double> & R);
  

  // set Hilbert space
  //
  // hilbertSpace = pointer to Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // set site to be acted on
  //
  // site = pointer to the siteto use 
  void SetSite (MPSSite* site);

  
  // get Hilbert space on which Hamiltonian acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();
  
  // return dimension of Hilbert space where Hamiltonian acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();
  
  // shift Hamiltonian from a given energy
  //
  // shift = shift value
  void ShiftHamiltonian (double shift);
  
  inline int GetMPODimension() const {return  MPOBondDimension;};
 protected:
  
  inline unsigned int GetIndiceDownFromTensorIndex(unsigned int tensorIndex);

  inline unsigned int GetIndiceUpFromTensorIndex(unsigned int tensorIndex);
  
  inline unsigned int GetIndiceLeftFromTensorIndex(unsigned int tensorIndex);
  
  inline unsigned int GetIndiceRightFromTensorIndex(unsigned int tensorIndex);
  
  inline unsigned int GetTensorIndexFromAllIndices(unsigned int indexDown, unsigned int indexUp, unsigned int indexLeft, unsigned int indexRight);

};

inline unsigned int AbstractMPOperatorOBC::GetIndiceDownFromTensorIndex(unsigned int tensorIndex)
{
  return  tensorIndex%this->PhysicalDimension;    
}


inline unsigned int AbstractMPOperatorOBC::GetIndiceUpFromTensorIndex(unsigned int tensorIndex)
{
  return  (tensorIndex/this->PhysicalDimension)%this->PhysicalDimension;    
}


inline unsigned int AbstractMPOperatorOBC::GetIndiceRightFromTensorIndex(unsigned int tensorIndex)
{
  return  (tensorIndex/(this->PhysicalDimension*this->PhysicalDimension)% this-> MPOBondDimension);    
}




inline unsigned int AbstractMPOperatorOBC::GetIndiceLeftFromTensorIndex(unsigned int tensorIndex)
{
  return  (tensorIndex/(this->PhysicalDimension*this->PhysicalDimension *this->MPOBondDimension));    
}


inline unsigned int AbstractMPOperatorOBC::GetTensorIndexFromAllIndices(unsigned int indexDown, unsigned int indexUp, unsigned int indexLeft, unsigned int indexRight)
{
  return  ((indexLeft * this-> MPOBondDimension  + indexRight)*this->PhysicalDimension +indexUp)*this->PhysicalDimension + indexDown;
}

#endif
