#ifndef _DMRGFINITESIZEREALOBCMAINTASK_H 
#define _DMRGFINITESIZEREALOBCMAINTASK_H

#include "AbstractMPOperatorOBC.h"
#include "Architecture/AbstractArchitecture.h"

class DMRGFiniteSizeRealOBCMainTask 
{
 protected:
  const int NbrSites;
  
  MPSSite * LatticeSite;
  AbstractArchitecture * Architecture;
  
  AbstractMPOperatorOBC * MPOperator;
  
  int NbrSweep;
  int MaximumBondDimension;
  
 public:
  
  DMRGFiniteSizeRealOBCMainTask(MPSSite * latticeSite, AbstractMPOperatorOBC * mPOperator, int nbrSites, int NbrSweep,int MaximumBondDimension,  AbstractArchitecture * architecture);
  virtual ~DMRGFiniteSizeRealOBCMainTask();
  void RunAlgorithm();
  
 protected:
  void InitializeLattice();
  RealVector & OptimizeUsingLanczosLanczosAlgorithm (int siteIndex);
};

#endif
