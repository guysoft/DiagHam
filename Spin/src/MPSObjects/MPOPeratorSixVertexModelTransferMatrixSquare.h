#ifndef _MPOperatorSixVertexModelTransferMatrixSquare_H
#define _MPOperatorSixVertexModelTransferMatrixSquare_H

#include "MPSObjects/AbstractMPOperatorOBC.h"


class MPOPeratorSixVertexModelTransferMatrixSquare : public AbstractMPOperatorOBC
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
