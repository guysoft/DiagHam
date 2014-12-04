#ifndef _REALMPOPeratorDefinedByFiles_H
#define _REALMPOPeratorDefinedByFiles_H

#include "MPSObjects/RealMPOperatorOBC.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

class RealMPOPeratorDefinedByFiles : public RealMPOperatorOBC
{
 protected:

  
 public:

  RealMPOPeratorDefinedByFiles();
  RealMPOPeratorDefinedByFiles(int nbrSites, MultiColumnASCIIFile & tensorElementsFile, MultiColumnASCIIFile & boundaryVectorsFile); 
  ~RealMPOPeratorDefinedByFiles();
  
  virtual void InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile);
  virtual void InitializeBoundaryVectors(MultiColumnASCIIFile & boundaryVectorsFile);

};

#endif
