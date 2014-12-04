#ifndef _DMRGFINITESIZECOMPLEXOBCMAINTASK_H 
#define _DMRGFINITESIZECOMPLEXOBCMAINTASK_H

#include "ComplexMPOperatorOBC.h"
#include "Architecture/AbstractArchitecture.h"
#include "LanczosAlgorithm/LanczosManager.h"

class DMRGFiniteSizeComplexOBCMainTask 
{
 protected:
  const int NbrSites;
  
  ComplexMPSSite * LatticeSite;
  AbstractArchitecture * Architecture;
  LanczosManager *  AlgorithmManager;
  AbstractMPOperatorOBC * MPOperator;
  
  int NbrSweep;
  int MaximumBondDimension;
  double PreviousEnergy;
 public:
  
  DMRGFiniteSizeComplexOBCMainTask(ComplexMPSSite * latticeSite, AbstractMPOperatorOBC * mPOperator, int nbrSites, int NbrSweep,int MaximumBondDimension,  AbstractArchitecture * architecture, LanczosManager* lanczos);
  virtual ~DMRGFiniteSizeComplexOBCMainTask();
  void RunAlgorithm();
  
 protected:
  void InitializeLattice();
  void InitializeLatticeUsingIDMRG();
  void OptimizeUsingLanczosLanczosAlgorithm (int siteIndex);
  void TwoSiteOptimizationUsingLanczosLanczosAlgorithm (ComplexMPSSite * leftSite , ComplexMPSSite * rightSite, RealDiagonalMatrix & singularValues);
};

#endif
