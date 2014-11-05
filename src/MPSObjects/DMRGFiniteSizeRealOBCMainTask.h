#ifndef _DMRGFINITESIZEREALOBCMAINTASK_H 
#define _DMRGFINITESIZEREALOBCMAINTASK_H

#include "AbstractMPOperatorOBC.h"
#include "Architecture/AbstractArchitecture.h"
#include "LanczosAlgorithm/LanczosManager.h"

class DMRGFiniteSizeRealOBCMainTask 
{
 protected:
  const int NbrSites;
  
  MPSSite * LatticeSite;
  AbstractArchitecture * Architecture;
  LanczosManager *  AlgorithmManager;
  AbstractMPOperatorOBC * MPOperator;
  
  int NbrSweep;
  int MaximumBondDimension;
  double PreviousEnergy;
 public:
  
  DMRGFiniteSizeRealOBCMainTask(MPSSite * latticeSite, AbstractMPOperatorOBC * mPOperator, int nbrSites, int NbrSweep,int MaximumBondDimension,  AbstractArchitecture * architecture, LanczosManager* lanczos);
  virtual ~DMRGFiniteSizeRealOBCMainTask();
  void RunAlgorithm();
  
 protected:
  void InitializeLattice();
  void InitializeLatticeUsingIDMRG();
  void OptimizeUsingLanczosLanczosAlgorithm (int siteIndex);
  void TwoSiteOptimizationUsingLanczosLanczosAlgorithm (MPSSite * leftSite , MPSSite * rightSite);
};

#endif
