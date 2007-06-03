////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//         Copyright (C) 2007 Gunnar Moeller               //
//                                                                            //
//                                                                            //
//           class for storing the history of a Monte-Carlo overlap calculation  //
//                                                                            //
//                        last modification : 28/05/2007                      //
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


#ifndef MCHISTORYRECORD_H
#define MCHISTORYRECORD_H

#include <fstream>
#include "AbstractMCHistoryData.h"
#include "GeneralTools/List.h"
#include "Vector/RealVector.h"
#include "MathTools/Complex.h"

using std::fstream;

class MCHistoryRecord
{
 protected:
  int LastSampleCount;
  int TotalSampleCount;
  int TotalRecordCount;
  int ProjectedStepNum;

  int NbrPositions;
  
  int NumAdditionalData;
  AbstractMCHistoryData **AdditionalData;
  unsigned SkipAdditional;

  ofstream LogFile;
  ifstream HistoryFile;

  int RecordMode;
  std::streampos StartPos;
 public:

  enum
    {
      Recording = 0x01,
      Reading = 0x02
    };

  // constructors
  // for recording mode
  // projectedSamplex: expected number of MC steps
  // nbrPositions: Number of coordinates for each call of wavefunction (2*N for a two-D system)
  MCHistoryRecord(int projectedSamples, int nbrPositions, char* exactFile, char* samplingDescriptor, char* fileName, List<AbstractMCHistoryData> *additionalData=NULL);
    
  // for reading mode
  MCHistoryRecord(char *Input, int nbrPositions, List<AbstractMCHistoryData> *additionalData=NULL);

  // destructor -> automatically closes LogFile
  ~MCHistoryRecord();
  
  // record rejected step - to be called for rejected Microsteps
  void RecordRejectedStep();

  // record accepted step - to be called for each accepted step, or at every step to be written to file
  bool RecordAcceptedStep( double samplingAmplitude, RealVector &positions, Complex &valueExact);

  // read one MC sample back from file, gives back the parameters in call of RecordAcceptedStep
  // sampleCount additionally gives the multiplicity of each record
  bool GetMonteCarloStep( int &sampleCount, double & samplingAmplitude, double *positions, Complex &valueExact);

  // rewind in reading mode:
  void RewindHistory();

  // get projected Samples
  int GetProjectedSamples() {return ProjectedStepNum;}
};

#endif // MCHISTORYRECORD_H
