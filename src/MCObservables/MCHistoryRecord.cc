#include "MCHistoryRecord.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ListIterator.h"
#include <iostream>

using std::ios;
using std::cout;
using std::endl;

// constructors
// for recording mode
MCHistoryRecord::MCHistoryRecord(int projectedSamples, int nbrPositions, char* exactFile, char* samplingDescriptor, char* fileName, List<AbstractMCHistoryData> *additionalData)
{
  this->RecordMode = MCHistoryRecord::Recording;
  this->LastSampleCount=0;
  this->TotalSampleCount=0;
  this->TotalRecordCount=0;
  this->ProjectedStepNum=projectedSamples;
  this->NbrPositions = nbrPositions;
  LogFile.open(fileName, ios::binary | ios::out);  
  WriteLittleEndian(LogFile, this->ProjectedStepNum);
  WriteLittleEndian(LogFile, this->NbrPositions);  
  if (additionalData != NULL)
    {
      this->NumAdditionalData=additionalData->GetNbrElement();
      WriteLittleEndian(LogFile, this->NumAdditionalData);
      this->AdditionalData=new AbstractMCHistoryData*[NumAdditionalData];
      AbstractMCHistoryData *Data;
      int i=0;
      unsigned tmpU;
      for (ListIterator<AbstractMCHistoryData> LI(*additionalData); (Data=LI())!=NULL;++i)
	{
	  tmpU=Data->GetHistoryDataType();
	  WriteLittleEndian(LogFile,tmpU);
	  tmpU=Data->GetSizeOnDisk();
	  WriteLittleEndian(LogFile,tmpU);
	  this->AdditionalData[i]=Data;
	}
    }
  else
    {
      this->NumAdditionalData=0;
      this->AdditionalData=NULL;
      WriteLittleEndian(LogFile, this->NumAdditionalData);
    }
}


// for reading mode
MCHistoryRecord::MCHistoryRecord(char *Input, int nbrPositions, List<AbstractMCHistoryData> *additionalData)
{
  this->RecordMode = MCHistoryRecord::Reading;
  HistoryFile.open(Input, ios::binary | ios::in);
  this->TotalSampleCount=0;
  ReadLittleEndian(HistoryFile, this->ProjectedStepNum);
  ReadLittleEndian(HistoryFile, this->NbrPositions);
  if (this->NbrPositions != nbrPositions)
    {
      cout << "Mismatch of nbrPositions in reading MCHistoryRecord" << endl;
      exit (-1);
    }
  int numAdditional;
  ReadLittleEndian(HistoryFile,numAdditional);
  this->AdditionalData=NULL;
  if (numAdditional > 0)
    {
      cout << "Additional data present in MCHistoryRecord "<< Input << endl;
      if (additionalData != NULL)
	{
	  if (numAdditional != additionalData->GetNbrElement())
	    {
	      cout << "Mismatch of the number of additional MCHistoryRecord's you wish to read" << endl;
	      exit (1);
	    }	  
	  this->NumAdditionalData=numAdditional;
	  this->AdditionalData=new AbstractMCHistoryData*[NumAdditionalData];
	  this->SkipAdditional=0;
	  AbstractMCHistoryData *Data;
	  int i=0;
	  for (ListIterator<AbstractMCHistoryData> LI(*additionalData); (Data=LI())!=NULL;++i)
	    {
	      unsigned check, size;
	      ReadLittleEndian(HistoryFile,check);
	      ReadLittleEndian(HistoryFile,size);
	      if (check != Data->GetHistoryDataType())
		{
		  cout << "Mismatch of field number "<<i+1<<" for the additional MCHistoryRecord's you wish to read" << endl;
		  exit (1);
		}
	      this->AdditionalData[i]=Data;
	    }
	}
      else  // data present, but not requested.
	{
	  cout << "Additional data present in file is being discarded." << endl;
	  this->NumAdditionalData=0;
	  this->SkipAdditional=0;
	  unsigned check,skip;
	  AbstractMCHistoryData *Data;
	  for (ListIterator<AbstractMCHistoryData> LI(*additionalData); (Data=LI())!=NULL;)
	    {
	      ReadLittleEndian(HistoryFile,check);
	      ReadLittleEndian(HistoryFile,skip);
	      this->SkipAdditional+=skip;
	    }
	}
    }
  else
    { // no additional data:
      this->NumAdditionalData=0;
      this->SkipAdditional=0;
    }
  // first entry from writing is lastSampleCount that was zero...
  ReadLittleEndian(HistoryFile, this->LastSampleCount);
  if (this->LastSampleCount!=0)
    {
      cout << "Problem with header of History record " <<Input << endl;
      exit(2);
    }
}


MCHistoryRecord::~MCHistoryRecord()
{
  if (this->RecordMode & MCHistoryRecord::Recording)
    {
      // write last multiplicity
      WriteLittleEndian(LogFile,this->LastSampleCount);
      cout << "Total " << TotalRecordCount << " History records written." << endl;
      char signature = 'e'; // signature for the END
      WriteLittleEndian(LogFile,signature);
      LogFile.flush();
      LogFile.close();
      if (AdditionalData != 0) delete [] AdditionalData;
    }
  else if (this->RecordMode & MCHistoryRecord::Reading)
    {
      HistoryFile.close();
      if (AdditionalData != 0) delete [] AdditionalData;
    }
}

// record rejected step - to be called for rejected Microsteps
void MCHistoryRecord::RecordRejectedStep()
{
  this->LastSampleCount++;
  this->TotalSampleCount++;
}

// record accepted step - to be called for each accepted step, or at every step to be written to file
bool MCHistoryRecord::RecordAcceptedStep( double samplingAmplitude, RealVector &positions, Complex &valueExact)
{
  this->TotalSampleCount++;
  this->TotalRecordCount++;
  // first: write multiplicity of last configuration (is part of prior record)
  WriteLittleEndian(LogFile,this->LastSampleCount);
  // then: write info about new positions (starting new record)
  char signature = 'b'; // signature for a new Block
  WriteLittleEndian(LogFile,signature);
  WriteLittleEndian(LogFile,samplingAmplitude);
  for (int i = 0; i < this->NbrPositions; ++i)
    WriteLittleEndian(LogFile, positions[i]);
  WriteLittleEndian(LogFile,valueExact);
  if (NumAdditionalData>0)
    {
      for (int i=0; i<NumAdditionalData; ++i) 
	this->AdditionalData[i]->WriteMCHistoryData();
    }
  this->LastSampleCount=1; // reset counter for present position
  return true;
}

bool MCHistoryRecord::GetMonteCarloStep( int &sampleCount, double &samplingAmplitude, double *positions, Complex &valueExact)
{
  if ( ! HistoryFile.eof())
    {
      char signature; 
      ReadLittleEndian(HistoryFile,signature);
      if (signature == 'b') // recognized a new block
	{
	  ReadLittleEndian(HistoryFile,samplingAmplitude);
	  for (int i = 0; i < this->NbrPositions; ++i)
	    ReadLittleEndian(HistoryFile, positions[i]);
	  ReadLittleEndian(HistoryFile,valueExact);
	  if (NumAdditionalData>0)
	    {
	      for (int i=0; i<NumAdditionalData; ++i) 
		this->AdditionalData[i]->ReadMCHistoryData();
	    }
	  else if (SkipAdditional>0)
	    {	      
	      HistoryFile.seekg(SkipAdditional,ios::cur);
	    }
	}
      else if (signature == 'e')
	{
	  // cout << "End of file detected!" << endl;
	  return false;
	}
      ReadLittleEndian(HistoryFile,this->LastSampleCount);
      this->TotalSampleCount+=this->LastSampleCount;
      sampleCount=this->LastSampleCount;
      return true;
    }
  return false;
}
