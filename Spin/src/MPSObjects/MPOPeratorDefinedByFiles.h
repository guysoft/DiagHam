#ifndef _MPOPeratorDefinedByFiles_H
#define _MPOPeratorDefinedByFiles_H

#include "MPSObjects/AbstractMPOperatorOBC.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

class MPOPeratorDefinedByFiles : public AbstractMPOperatorOBC
{
 protected:

  
 public:

  MPOPeratorDefinedByFiles();
  MPOPeratorDefinedByFiles(int nbrSites, MultiColumnASCIIFile tensorElementsFile, MultiColumnASCIIFile boundaryVectorsFile); 
  ~MPOPeratorDefinedByFiles();
  
  virtual void InitializeTensorsElements(MultiColumnASCIIFile tensorElementsFile);
  virtual void InitializeBoundaryVectors(MultiColumnASCIIFile boundaryVectorsFile);


  void PrintTensorElements();
 
};

#endif
