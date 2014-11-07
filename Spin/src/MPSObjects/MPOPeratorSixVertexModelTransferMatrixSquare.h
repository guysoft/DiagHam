#ifndef _MPOperatorSixVertexModelTransferMatrixSquare_H
#define _MPOperatorSixVertexModelTransferMatrixSquare_H

#include "MPSObjects/AbstractMPOperatorOBC.h"


class MPOPeratorSixVertexModelTransferMatrixSquare : public AbstractMPOperatorOBC
{
 protected:
  double *  LeftVector;
  double * RightVector;
  
 public:

  MPOPeratorSixVertexModelTransferMatrixSquare();
  MPOPeratorSixVertexModelTransferMatrixSquare(int nbrSites);
  ~MPOPeratorSixVertexModelTransferMatrixSquare();
  
  virtual void InitializeTensorsElements();
  void ComputeL(Tensor3<double> & L);
  void ComputeR(Tensor3<double> & R);
  void PrintTensorElements();
  
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent);

  virtual RealVector & LowLevelMultiplyTwoSites(RealVector& vSource, RealVector& vDestination,  int firstComponent, int nbrComponent);
  virtual RealVector & LowLevelMultiplyOneSite(RealVector& vSource, RealVector& vDestination,  int firstComponent, int nbrComponent);
};

#endif
