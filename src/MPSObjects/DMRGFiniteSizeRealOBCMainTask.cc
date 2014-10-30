#include "DMRGFiniteSizeRealOBCMainTask.h"
#include "MPSSite.h"
#include <iostream>
#include <sys/time.h>
#include "Architecture/AbstractArchitecture.h"
#include "LanczosAlgorithm/BasicLanczosAlgorithmWithGroundState.h"
#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"

using std::cout;
using std::endl;

DMRGFiniteSizeRealOBCMainTask::DMRGFiniteSizeRealOBCMainTask(MPSSite * latticeSite, AbstractMPOperatorOBC * mPOperator, int nbrSites, int nbrSweep, int maximumBondDimension,  AbstractArchitecture * architecture) :  NbrSites(nbrSites)
{
  this->LatticeSite = latticeSite;
  this->MPOperator = mPOperator;
  this->NbrSweep = nbrSweep;
  this->MaximumBondDimension = maximumBondDimension;
  this->Architecture = architecture;
}


DMRGFiniteSizeRealOBCMainTask::~DMRGFiniteSizeRealOBCMainTask(){ }

void DMRGFiniteSizeRealOBCMainTask::RunAlgorithm()
{
  this->InitializeLattice();

  for(int CurrentSweep = 0; CurrentSweep <  this->NbrSweep; CurrentSweep++)
    {
      cout <<"Sweep : "<<CurrentSweep<<endl;
      
      for (int i = 0; i <  this->NbrSites - 1; i++)
	{
	  cout <<"From left to right "<<" site "<<i <<endl;
          this->MPOperator->SetSite(&this->LatticeSite[i]);
	  this->OptimizeUsingLanczosLanczosAlgorithm (i);
	  this->LatticeSite[i].BringMInLeftCanonicalFormCareful();  
	}
     for (int i =  this->NbrSites - 1; i >  0; i--)
	{
	  cout <<"From right to left,"<<" site "<<i <<endl;
	  this->MPOperator->SetSite(&this->LatticeSite[i]);
	  this->OptimizeUsingLanczosLanczosAlgorithm (i);
	  this->LatticeSite[i].BringMInRightCanonicalFormCareful();
	}
    }
}




void DMRGFiniteSizeRealOBCMainTask::InitializeLattice()
{
  for (int i = 0 ; i < NbrSites ; i++) 
    {
      LatticeSite[i].InitializeWithRandomMatrices();
    }

  for (int i = NbrSites - 1 ; i>0 ; i--) 
    {
      LatticeSite[i].BringMInRightCanonicalFormCareful(); 
    }
}



void DMRGFiniteSizeRealOBCMainTask::OptimizeUsingLanczosLanczosAlgorithm (int siteIndex)
{
  if (this->MPOperator->GetHilbertSpaceDimension() < 500 )
    {
    RealSymmetricMatrix HRep (this->MPOperator->GetHilbertSpaceDimension(), true);
    this->MPOperator->GetHamiltonian(HRep);
    
  //  cout <<"matrix form of the Ham" <<HRep<<endl;
    if (this->MPOperator->GetHilbertSpaceDimension() > 1)
     {
#ifdef __LAPACK__
		  RealDiagonalMatrix TmpDiag (this->MPOperator->GetHilbertSpaceDimension());
                  RealMatrix Q(this->MPOperator->GetHilbertSpaceDimension(), this->MPOperator->GetHilbertSpaceDimension());
		  HRep.LapackDiagonalize(TmpDiag, Q);
		  RealVector TmpEigenvector(this->MPOperator->GetHilbertSpaceDimension(),true);

                  this->MPOperator->LowLevelMultiply(Q[0], TmpEigenvector);

                  cout <<"TmpDiag[0] = " <<TmpDiag[0]<<endl;

                  cout << " Highest Energy = " <<(TmpEigenvector * Q[0]) << " " << endl;

         	  this->LatticeSite[siteIndex].UpdateFromVector(&Q[0]);


#endif

}
}
else
{
  RealVector * TmpVector = 0;
  BasicLanczosAlgorithmWithGroundState LanczosAlgorithm (this->Architecture);
  this->LatticeSite[siteIndex].GetMatrixInVectorForm(TmpVector);
  //cout <<(*TmpVector)<<endl;
  LanczosAlgorithm.InitializeLanczosAlgorithm(*TmpVector);
  double GroundStateEnergy;
  double Precision = 1.0;
  double PreviousLowest = 1e50;
  double Lowest = PreviousLowest;
  int CurrentNbrIterLanczos = 0;
  LanczosAlgorithm.SetHamiltonian(this->MPOperator);
  cout << "Run Lanczos Algorithm" << endl;
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  timeval TotalCurrentTime;
  double Dt;
  gettimeofday (&(TotalStartingTime), 0);
  int StartTimeSecond = TotalStartingTime.tv_sec;
  LanczosAlgorithm.RunLanczosAlgorithm(3);
  CurrentNbrIterLanczos = 4;
  RealTriDiagonalSymmetricMatrix TmpMatrix;  
  int CurrentTimeSecond = TotalCurrentTime.tv_sec;
  while ((LanczosAlgorithm.TestConvergence() == false)&&( CurrentNbrIterLanczos < 5 ))
  {
     ++CurrentNbrIterLanczos;
     LanczosAlgorithm.RunLanczosAlgorithm(1);
     TmpMatrix.Copy(LanczosAlgorithm.GetDiagonalizedMatrix());
     TmpMatrix.SortMatrixUpOrder();
    Lowest = TmpMatrix.DiagonalElement(0);
	  Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
	  PreviousLowest = Lowest; 
    cout << (TmpMatrix.DiagonalElement(0)) << " " << Lowest << " " << Precision << " ";
	  gettimeofday (&(TotalEndingTime), 0);
	  CurrentTimeSecond = TotalEndingTime.tv_sec;
    Dt = (double) (TotalEndingTime.tv_sec - TotalCurrentTime.tv_sec) + 
		((TotalEndingTime.tv_usec - TotalCurrentTime.tv_usec) / 1000000.0);		      
	      cout << "(" << Dt << " s for step " << CurrentNbrIterLanczos << ")";
	      TotalCurrentTime.tv_usec = TotalEndingTime.tv_usec;
	      TotalCurrentTime.tv_sec = TotalEndingTime.tv_sec;
	  cout << endl;

  }
      GroundStateEnergy = Lowest;
      cout << endl;
      cout << (TmpMatrix.DiagonalElement(0)) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	   << CurrentNbrIterLanczos << endl;
      RealVector GroundState = LanczosAlgorithm.GetGroundState();
      //cout <<"GroundState = " <<GroundState<<endl;
      this->LatticeSite[siteIndex].UpdateFromVector(&GroundState);
}

}
