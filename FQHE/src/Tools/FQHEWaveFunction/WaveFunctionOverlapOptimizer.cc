#include "WaveFunctionOverlapOptimizer.h"
#include "Matrix/RealMatrix.h"
#include <iostream>

using std::cout;
using std::endl;

// number of orders of magnitude of parameter values to scan
#define START_DIFFERENTIAL 0.0001
#define DIFFERENTIAL_SPREAD  7

WaveFunctionOverlapOptimizer::WaveFunctionOverlapOptimizer( Abstract1DComplexTrialFunction *trialState, char *historyFileName, int nbrParticles, bool excludeLastParameter, int maxPoints)
{
  this->NbrParticles = nbrParticles;  
  this->History = new MCHistoryRecord(historyFileName, 2*nbrParticles /* could add additional observables here */);
  this->Positions.Resize(2*nbrParticles);
  this->TrialState = trialState;
  this->NbrParameters = this->TrialState->GetNbrParameters();
  this->LastParameterExcluded = excludeLastParameter; // do not optimize last trial parameter
  if (excludeLastParameter)
    this->EffectiveNbrParameters = this->NbrParameters-1;
  else
    this->EffectiveNbrParameters = this->NbrParameters;
  this->Gradient.Resize(this->NbrParameters);
  this->StepDirection.Resize(this->NbrParameters);
  this->MaxPoints = maxPoints;
  this->MaxParameters = maxPoints*(1+2*this->EffectiveNbrParameters);
  if (this->MaxParameters<1+2*DIFFERENTIAL_SPREAD*this->EffectiveNbrParameters)
    this->MaxParameters = 1+2*DIFFERENTIAL_SPREAD*this->EffectiveNbrParameters;
  this->Differentials = new double[this->NbrParameters];
  this->NormObservation = new double[MaxParameters];
  this->OverlapObservation = new Complex[MaxParameters];
  this->NewParameters = new double*[MaxParameters];
  for (int i=0; i<MaxParameters; ++i) this->NewParameters[i]= new double[this->NbrParameters];
  this->InitialParameters = new double[this->NbrParameters];
  double *tmpDs = this->TrialState->GetTrialParameters();
  for (int i=0; i<this->NbrParameters; ++i) this->InitialParameters[i]=tmpDs[i];
  this->typicalSA=0.0;
  this->typicalTV=0.0;
  this->typicalWF=0.0;
  this->NormTrialObs = NULL;
  this->OverlapObs = NULL;
  
  // get an idea of the average amplitude of the different numbers in HistoryFile:
  int sampleCount;
  double SamplingAmplitude;
  Complex ExactValue;
  int toCheck = History->GetProjectedSamples();
  if (toCheck > 1000) toCheck=1000;
  int averageTypical=toCheck/2;
  int initialSkip=toCheck/2;
  bool HaveMoreHistory=true;
  for (int i=0; (HaveMoreHistory&&(i<initialSkip)); ++i)
    HaveMoreHistory=History->GetMonteCarloStep(sampleCount, SamplingAmplitude, &(Positions[0]), ExactValue);
  if (HaveMoreHistory)
    {
      HaveMoreHistory=true;
      for (int i=0;(HaveMoreHistory&&(i<averageTypical)); ++i)
	{
	  HaveMoreHistory=History->GetMonteCarloStep(sampleCount, SamplingAmplitude, &(Positions[0]), ExactValue);
	  typicalSA+=SamplingAmplitude;
	  typicalWF+=Norm(ExactValue);
	  typicalTV+=Norm((*TrialState)(Positions));
	}
      this->typicalSA/=averageTypical;
      this->typicalWF/=averageTypical;
      this->typicalTV/=averageTypical;
    }
  else
    {
      this->typicalSA=1.0;
      this->typicalWF=1.0;
      this->typicalTV=1.0;
    }
  
  History->RewindHistory();
  WeightedRealObservable NormExactObs(256);
  // calculate norm of exact wavefunction, once and for all:
  while (History->GetMonteCarloStep(sampleCount, SamplingAmplitude, &(Positions[0]), ExactValue))
	{
	  ExactValue /= typicalWF;
	  SamplingAmplitude /= typicalSA;
	  NormExactObs.Observe(SqrNorm(ExactValue)/SamplingAmplitude,(double)sampleCount);
	}
  this->NormExactWF = NormExactObs.Average();
  this->ErrorNormExactWF = NormExactObs.ErrorEstimate();
}

WaveFunctionOverlapOptimizer::~WaveFunctionOverlapOptimizer()
{
  delete History;
  for (int i=0; i<MaxParameters; ++i) delete [] this->NewParameters[i];  
  delete [] this->NewParameters;
  delete [] this->NormObservation;
  delete [] this->OverlapObservation;
  delete [] this->InitialParameters;
  delete [] this->Differentials;
}


// calculate Trial wavefunction for parameters that have been entered into NewParameters.
// the dimension of ManyValues indicates the number of parameters.
//
void WaveFunctionOverlapOptimizer::EvaluateTrialOverlaps()
{
  int sampleCount;
  double SamplingAmplitude;
  Complex Tmp, ExactValue;
  int nbrParametersInEvaluation = this->ManyValues.GetVectorDimension();
  if (this->NormTrialObs!=NULL) delete this->NormTrialObs;
  this->NormTrialObs = new WeightedRealVectorObservable(nbrParametersInEvaluation,256);
  if (this->OverlapObs != NULL) delete this->OverlapObs;
  this->OverlapObs = new WeightedComplexVectorObservable(nbrParametersInEvaluation,256);
  History->RewindHistory();
  int count =0;
  while ( History->GetMonteCarloStep(sampleCount, SamplingAmplitude, &(this->Positions[0]), ExactValue))
    {
      TrialState->GetForManyParameters(this->ManyValues, this->Positions, this->NewParameters);
      SamplingAmplitude /= typicalSA;
      ExactValue /= typicalWF;
      for (int i=0; i<nbrParametersInEvaluation; ++i)
	{
	  Tmp=ManyValues[i]/typicalTV;
	  NormObservation[i] = SqrNorm(Tmp)/SamplingAmplitude;
	  OverlapObservation[i] = Conj(Tmp)*ExactValue/SamplingAmplitude;
	}
      NormTrialObs->Observe(NormObservation,(double)sampleCount);
      OverlapObs->Observe(OverlapObservation,(double)sampleCount);
      count++;
    }
  for (int i=0; i<nbrParametersInEvaluation; ++i)
    {
      OverlapObservation[i] = OverlapObs->Average(i) / sqrt(NormTrialObs->Average(i)*this->NormExactWF);
      NormObservation[i] = SqrNorm(OverlapObservation[i]);
    }
}


void WaveFunctionOverlapOptimizer::DetermineGradientAndDifferentials(double *parameters)
{  
  this->ManyValues.Resize(1+2*DIFFERENTIAL_SPREAD*this->EffectiveNbrParameters);
  // set parameters for which the overlap shall be calculated:
  for (int s=0; s<1+2*DIFFERENTIAL_SPREAD*EffectiveNbrParameters; ++s)
    for (int i=0; i<NbrParameters; ++i)
      {
	this->NewParameters[s][i]=parameters[i];
	this->Differentials[i]=0.0;
	this->Gradient[i]=0.0;
      }
  int s=1;  
  for (int p=0; p<EffectiveNbrParameters; ++p)
    {
      double dA=START_DIFFERENTIAL;
      for (int i=0; i<DIFFERENTIAL_SPREAD; ++i)
	{
	  this->NewParameters[s++][p]+=dA;
	  this->NewParameters[s++][p]-=dA;
	  dA*=10.0;
	}
    }
  // define internal function with the instructions to calculate averages for the given parameters in ManyValues, NewParameters!  
  this->EvaluateTrialOverlaps();
  // OverlapObservation and NormObservation now contain data with the overlaps and their squares.
  this->MinDifferential=0.5;
  double maxOverlap;
  double *ProposedStepLength = new double[EffectiveNbrParameters];
  bool differentialToBeFound;
  for (int p=0; p<EffectiveNbrParameters; ++p)
    {
      double dA=START_DIFFERENTIAL;
      maxOverlap=NormObservation[0];      
      this->Differentials[p]=0.5;
      ProposedStepLength[p] = 0.5;
      differentialToBeFound=true;
      for (int i=0; i<DIFFERENTIAL_SPREAD; ++i)
	{	  
	  if (this->NormObservation[1+2*DIFFERENTIAL_SPREAD*p+2*i]>maxOverlap)
	    {
	      maxOverlap=this->NormObservation[1+2*DIFFERENTIAL_SPREAD*p+2*i];
	      //cout << "New maximum "<< maxOverlap<<" found for p="<<p<<" dA="<<dA<< endl;
	      ProposedStepLength[p] = dA;
	    }
	  if (this->NormObservation[2+2*DIFFERENTIAL_SPREAD*p+2*i]>maxOverlap)
	    {
	      maxOverlap=this->NormObservation[2+2*DIFFERENTIAL_SPREAD*p+2*i];
	      //cout << "New maximum "<< maxOverlap<<" found for p="<<p<<" dA="<<-dA<< endl;
	      ProposedStepLength[p] = dA;
	    }
	  if ((differentialToBeFound)&&(fabs(this->NormObservation[1+2*DIFFERENTIAL_SPREAD*p+2*i]
					     - this->NormObservation[1+2*DIFFERENTIAL_SPREAD*p+2*i+1]) /
					this->NormObservation[0] > 1e-4))
	    {
	      this->Gradient[p]=(this->NormObservation[1+2*DIFFERENTIAL_SPREAD*p+2*i] - this->NormObservation[1+2*DIFFERENTIAL_SPREAD*p+2*i+1]) / dA;
	      this->Differentials[p]=dA;
	      differentialToBeFound=false;
	      if (dA<this->MinDifferential) this->MinDifferential=dA;
	    }
	  dA*=10.0;
	}      
    }
  // do not allow StepLength to be bigger than 5 times the average differential, to start with.
  double maxGradient=0.0;
  int iMax=0;
  for (int i=0; i<EffectiveNbrParameters; ++i)
    if (fabs(Gradient[i])>maxGradient) { maxGradient = fabs(Gradient[i]); iMax=i; }
  this->StepLength = ProposedStepLength[iMax]*3.0/(double)this->MaxPoints;
  //if(this->StepLength > 100.0 * Differentials[iMax]) this->StepLength = 100.0 * Differentials[iMax];
  delete [] ProposedStepLength;
}

void WaveFunctionOverlapOptimizer::CalculateLinearSeries(RealVector &startParameters, RealVector &overlaps, RealMatrix &gradients)
{
  this->ManyValues.Resize(this->MaxParameters);
  // set parameters for which the overlap shall be calculated:
  int StatesPerPoint = 1+2*this->EffectiveNbrParameters;
  RealVector TmpParameters=startParameters;
  RealVector TmpStep= this->StepLength*this->StepDirection;
  for (int p=0; p<MaxPoints; ++p)
    {
      for (int i=0; i<NbrParameters; ++i)
	{
	  this->NewParameters[p*StatesPerPoint][i]=TmpParameters[i];
	  for (int j=0; j<EffectiveNbrParameters; ++j)
	    {
	      this->NewParameters[p*StatesPerPoint+2*j+1][i]=TmpParameters[i];
	      this->NewParameters[p*StatesPerPoint+2*j+2][i]=TmpParameters[i];
	    }
	}
      for (int j=0; j<EffectiveNbrParameters; ++j)
	{
	  this->NewParameters[p*StatesPerPoint+2*j+1][j]+=Differentials[j];
	  this->NewParameters[p*StatesPerPoint+2*j+2][j]-=Differentials[j];
	}
      TmpParameters+=TmpStep;
    }
  // define internal function with the instructions to calculate averages for the given parameters in ManyValues, NewParameters!  
  this->EvaluateTrialOverlaps();
  // OverlapObservation and NormObservation now contain data with the overlaps and their squares.
  for (int p=0; p<MaxPoints; ++p)
    {
      overlaps[p]=NormObservation[p*StatesPerPoint];
      for (int j=0; j<EffectiveNbrParameters; ++j)
	gradients[p][j]=(NormObservation[p*StatesPerPoint+2*j+1]-NormObservation[p*StatesPerPoint+2*j+2])/Differentials[j];
      for (int j=EffectiveNbrParameters; j<NbrParameters; ++j) gradients[p][j]=0.0;
    }
}
  
double WaveFunctionOverlapOptimizer::GetMaximumSqrOverlap(RealVector &optimalParameters, Complex &Overlap, double toleranceFinal, double toleranceIteration)
{
  // variables that are updated in the loop below:
  RealVector presentParameters(this->InitialParameters, this->NbrParameters);
  RealVector overlaps(this->MaxPoints);
  RealMatrix gradients(this->NbrParameters, this->MaxPoints);
  int StatesPerPoint = 1+2*this->EffectiveNbrParameters;
 
  // initial overlap and gradient of overlap, here:
  this->DetermineGradientAndDifferentials(InitialParameters);
  InitialSqrOverlap = NormObservation[0];
  cout << "Initial overlap = " << InitialSqrOverlap << endl;
  cout << "Initial parameters: " << presentParameters<<endl;
  cout << "Gradient is: " << Gradient.Norm() << endl;
  cout << "chose stepLength: "<<StepLength<<endl;
  this->StepDirection = this->Gradient;
  this->StepDirection.Normalize();
  cout << "StepDirection is: " << StepDirection;
  // ultimately, start a loop here:
  while ((this->Gradient.Norm() > toleranceFinal) && (StepLength > 1e-10))
    {
      // calculate a series of points (and gradients) along StepDirection, and spaced with StepLength
    
      this->CalculateLinearSeries(presentParameters, overlaps, gradients);
      double lastOverlap=overlaps[0];
      double scalarProduct=StepDirection*gradients[0]/gradients[0].Norm();
      int point=1;
      while ((point < MaxPoints) && (overlaps[point]>lastOverlap) && (scalarProduct>toleranceIteration))
	{	
	  scalarProduct = StepDirection*gradients[point]/gradients[point].Norm();
	  lastOverlap = overlaps[point];
	  point ++;
	}
      point--;
      cout << "Overlap "<<overlaps[point] <<" at point " <<point+1 <<"/"<<MaxPoints<<", Gradient: " << gradients[point].Norm() << endl;
      // update the new search direction and StepLength:
      if (point==0)
	{
	  // same direction, but smaller StepLength;
	  StepLength/=MaxPoints;
	  this->Gradient = gradients[0];
	  this->StepDirection = this->Gradient;
	  this->StepDirection.Normalize();
	  for (int i=0; i<this->EffectiveNbrParameters; ++i) Differentials[i]/=MaxPoints;
	  // cout << "New direction was not promising." << endl;
// 	  cout << "New direction:" << StepDirection << endl;
// 	  cout << "New StepLength:" << StepLength << endl;
	}
      else 
	{
	  // maybe separate case if continued until the last point?
	  this->StepDirection = (this->Gradient*((double)(gradients[point].SqrNorm()/this->Gradient.SqrNorm())));
	  this->StepDirection += gradients[point];
	  this->StepDirection.Normalize();
	  if (point<MaxPoints-1) // if we found a minimum before the end of our series -> reduce StepLength
	    this->StepLength*=0.5;
	  else // otherwise, increase StepLength slightly...
	    this->StepLength*=1.5;
	  this->Gradient = gradients[point];
	  for (int i=0; i<this->NbrParameters; ++i) presentParameters[i]=this->NewParameters[point*StatesPerPoint][i];
	  cout << "New parameters: " << presentParameters<<endl;
// 	  cout << "New direction:" << StepDirection << endl;
// 	  cout << "New StepLength:" << StepLength << endl;
	}      
    }
  return 0.0;
}



// double OverlapError(WeightedComplexObservable &ScalarProduct, WeightedRealObservable &NormObs1, WeightedRealObservable &NormObs2)
// //double OverlapError(ComplexObservable &ScalarProduct, RealObservable &NormObs1, RealObservable &NormObs2)
// {
//   double norm1 = NormObs1.Average();
//   double norm2 = NormObs2.Average();
//   double norm = sqrt(norm1*norm2);
//   double prod = Norm(ScalarProduct.Average());
//   return sqrt( DSQR(ScalarProduct.ErrorEstimate()/norm) + DSQR(prod*NormObs1.ErrorEstimate()*NormObs2.Average()/2.0/norm/norm/norm) + DSQR(prod*NormObs2.ErrorEstimate()*NormObs1.Average()/2.0/norm/norm/norm));
// }
