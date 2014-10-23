#ifndef _MPOperatorSixVertexModelTransferMatrixSquare_H
#define _MPOperatorSixVertexModelTransferMatrixSquare_H

#include "MPSObjects/AbstractMPOperatorOBC.h"


class MPOPeratorSixVertexModelTransferMatrixSquare : public AbstractMPOperatorOBC
{
 protected:
  double *  LeftVector;
  double * RightVector;
  
 public:
  
  MPOPeratorSixVertexModelTransferMatrixSquare(int nbrSites =0 );
  ~MPOPeratorSixVertexModelTransferMatrixSquare();
  
  virtual void InitializeTensorsElements();
  void ComputeL(Tensor3<double> & L);
  void ComputeR(Tensor3<double> & R);
  void PrintTensorElements();
  
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent);

};

#endif
