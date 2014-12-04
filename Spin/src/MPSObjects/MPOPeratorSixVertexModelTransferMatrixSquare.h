#ifndef _MPOperatorSixVertexModelTransferMatrixSquare_H
#define _MPOperatorSixVertexModelTransferMatrixSquare_H

#include "MPSObjects/RealMPOperatorOBC.h"


class MPOPeratorSixVertexModelTransferMatrixSquare : public RealMPOperatorOBC
{
 protected:

 public:

  MPOPeratorSixVertexModelTransferMatrixSquare();
  MPOPeratorSixVertexModelTransferMatrixSquare(int nbrSites);
  ~MPOPeratorSixVertexModelTransferMatrixSquare();
  
  virtual void InitializeTensorsElements();
  void PrintTensorElements();
};

#endif
