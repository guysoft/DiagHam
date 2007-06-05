////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2007 Gunnar Möller                  //
//                                                                            //
//                                                                            //
//           class of Jain composite fermion wave function on sphere          //
//                      with filled (pseudo) Landau levels                    //
//                                                                            //
//                        last modification : 18/05/2007                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#ifndef WAVEFUNCTIONOVERLAPOPTIMIZER_H
#define WAVEFUNCTIONOVERLAPOPTIMIZER_H


#include "config.h"
#include "MathTools/Complex.h"
#include "MathTools/NumericalAnalysis/Abstract1DComplexTrialFunction.h"
#include "MCObservables/WeightedRealObservable.h"
#include "MCObservables/WeightedRealVectorObservable.h"
#include "MCObservables/WeightedComplexVectorObservable.h"
#include "MCObservables/MCHistoryRecord.h"
#include <fstream>
using std::ofstream;

class WaveFunctionOverlapOptimizer
{
 protected:
  int NbrParticles;
  int NbrParameters;
  int EffectiveNbrParameters;
  Abstract1DComplexTrialFunction *TrialState;
  RealVector Positions;
  RealVector Gradient;
  ComplexVector ManyValues;
  double *InitialParameters;
  double InitialSqrOverlap;
  double **NewParameters;
  double *NormObservation;
  double *Differentials;
  double MinDifferential;
  double StepLength;
  Complex *OverlapObservation;
  int MaxPoints;
  int MaxParameters;  
  double NormExactWF;
  double ErrorNormExactWF;
  double OutlierLimit;
  WeightedRealVectorObservable *NormTrialObs;
  WeightedComplexVectorObservable *OverlapObs;
  MCHistoryRecord *History;
  bool LastParameterExcluded;
  double typicalSA;
  double typicalTV;
  double typicalWF;
  ofstream LogFile;
  
 public:

  WaveFunctionOverlapOptimizer( Abstract1DComplexTrialFunction *trialState, char *historyFileName, int nbrParticles, bool excludeLastParameter = true, int maxPoints = 50, char* logFileName = NULL);
  ~WaveFunctionOverlapOptimizer();
  
  double GetMaximumSqrOverlap(RealVector &optimalParameters, Complex &Overlap,
			      double toleranceFinal=1e-6, double toleranceIteration =0.01);

 private:

  void EvaluateTrialOverlaps();
  void DetermineGradientAndDifferentials(double *parameters);
  void CalculateLinearSeries(RealVector &startParameters, RealVector &stepDirection, RealVector &overlaps, RealMatrix &gradients);
  
};

#endif // WAVEFUNCTIONOVERLAPOPTIMIZER_H
