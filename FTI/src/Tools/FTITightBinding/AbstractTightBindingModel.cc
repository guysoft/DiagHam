////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of abstract tight binding model                    //
//                                                                            //
//                        last modification : 01/05/2012                      //
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


#include "config.h"
#include "Tools/FTITightBinding/AbstractTightBindingModel.h"
#include "Architecture/ArchitectureOperation/FTIComputeBandStructureOperation.cc"
#include "GeneralTools/Endian.h"
#include "GeneralTools/OrderedList.h"
#include <sys/time.h>
#include <fstream>


using std::ofstream;
using std::ios;
using std::endl;

// constructor
//

AbstractTightBindingModel::AbstractTightBindingModel()
{
    this->Inversion = ComplexMatrix();
}

// default constructor
//

AbstractTightBindingModel::AbstractTightBindingModel()
{
  this->Architecture = 0;
}


// destructor
//

AbstractTightBindingModel::~AbstractTightBindingModel()
{
}

// write an ASCII header that describes the tight binding model
// 
// output = reference on the output stream
// commentChar = optional ASCII character to add in front of each header line
// return value  = reference on the output stream

ostream& AbstractTightBindingModel::WriteASCIIHeader(ostream& output, char commentChar)
{
  return output;
}

// write an header that describes the tight binding model
// 
// output = reference on the output stream
// return value  = reference on the output stream

ofstream& AbstractTightBindingModel::WriteHeader(ofstream& output)
{
  int HeaderSize = 0;
  WriteLittleEndian(output, HeaderSize);
  return output; 
}

// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool AbstractTightBindingModel::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# index";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  File << endl;
  for (int kx = 0; kx < this->NbrStatePerBand; ++kx)
    {
      File << kx; 
      for (int i = 0; i < this->NbrBands; ++i)
	File << " " << this->GetEnergy(i, kx);
      File << endl;
    }
  File.close();
  return true;
}

// compute the spread of a given band
// 
// band = band index
// return value = band spread 

double AbstractTightBindingModel::ComputeBandSpread(int band)
{
  if ((band < 0) || (band > this->NbrBands))
    {
      return -1.0;
    }
  double MinBand = this->GetEnergy(band, 0);
  double MaxBand = this->GetEnergy(band, 0);
  for (int i = 1; i < this->NbrStatePerBand; ++i)
    {
      if (this->GetEnergy(band, i) > MaxBand)
	{
	  MaxBand = this->GetEnergy(band, i);
	}
      else
	{
	  if (this->GetEnergy(band, i) < MinBand)
	    {
	      MinBand = this->GetEnergy(band, i);
	    }
	}
    }
  return (MaxBand - MinBand);
}

// compute the direct gap between two bands
// 
// band1 = index of the lower band 
// band2 = index of the upper band (-1 if it has to be set to band1 + 1)
// return value =  direct band gap

double AbstractTightBindingModel::ComputeDirectBandGap(int band1, int band2)
{
  if ((band1 < 0) || (band1 > this->NbrBands))
    {
      return -1.0;
    }
  if (band2 < 0)
    band2 = band1 + 1;
  if ((band2 < 0) || (band2 > this->NbrBands) || (band2 <= band1))
    {
      return -1.0;
    }
  double Gap = this->GetEnergy(band2, 0) - this->GetEnergy(band1, 0);
  for (int i = 1; i < this->NbrStatePerBand; ++i)
    {
      double TmpGap = this->GetEnergy(band2, i) - this->GetEnergy(band1, i);
      if (TmpGap < Gap)
	Gap = TmpGap;
    }
  return Gap;
}

// compute the direct gap between two bands
// 
// band1 = index of the lower band 
// band2 = index of the upper band (-1 if it has to be set to band1 + 1)
// return value =  direct band gap

double AbstractTightBindingModel::ComputeIndirectBandGap(int band1, int band2)
{
  if ((band1 < 0) || (band1 > this->NbrBands))
    {
      return -1.0;
    }
  if (band2 < 0)
    band2 = band1 + 1;
  if ((band2 < 0) || (band2 > this->NbrBands) || (band2 <= band1))
    {
      return -1.0;
    }
  double Max1=this->GetEnergy(band1, 0);
  double Min2=this->GetEnergy(band2, 0);
  for (int i = 1; i < this->NbrStatePerBand; ++i)
    {
      double Tmp1 = this->GetEnergy(band1, i);
      double Tmp2 = this->GetEnergy(band2, i);
      if (Tmp1 > Max1)
	Max1 = Tmp1;
      if (Tmp2 < Min2)
	Min2 = Tmp2;
    }
  return Min2-Max1;
}


class SingleParticleStateData
{
public:
  SingleParticleStateData();
  SingleParticleStateData(double e, double b, double i);
  ~SingleParticleStateData(){};
  double Energy;
  double BandIndex;
  double StateIndex;

  friend bool operator<(const SingleParticleStateData &d1, const SingleParticleStateData &d2);
  friend bool operator>(const SingleParticleStateData &d1, const SingleParticleStateData &d2);
  friend bool operator==(const SingleParticleStateData &d1, const SingleParticleStateData &d2);

  friend ostream& operator<<(ostream &str, const SingleParticleStateData &d);
  

};

SingleParticleStateData::SingleParticleStateData()
{
  this->Energy = 0.0;
  this->BandIndex = -1;
  this->StateIndex = -1;
}

SingleParticleStateData::SingleParticleStateData(double e, double b, double i)
{
  this->Energy = e;
  this->BandIndex = b;
  this->StateIndex = i;
}

bool operator<(const SingleParticleStateData &d1, const  SingleParticleStateData &d2)
{
  return (d1.Energy < d2.Energy);
}

bool operator>(const SingleParticleStateData &d1, const SingleParticleStateData &d2)
{
  return (d1.Energy > d2.Energy);
}

bool operator==(const SingleParticleStateData &d1, const SingleParticleStateData &d2)
{
  return (d1.Energy == d2.Energy);
}

ostream &operator<<(ostream &str, const SingleParticleStateData &d)
{
  str <<d.BandIndex<<"_"<<d.StateIndex<<" ("<<d.Energy<<")";
  return str;
}



// compute the ground state energy for a number of fermions filling the band 
// 
// nbrFermions = number of particles in the band structure
// bands = number of bands used in groundstate configuration
// return value =  total groundstate energy
double AbstractTightBindingModel::ComputeGroundstateEnergy(int nbrFermions, int &bands, bool verbose)
{
  bands = -1;
  if (nbrFermions<1) return 0.0;

  if (nbrFermions>this->NbrBands*this->NbrStatePerBand) return 0.0;
  OrderedList<SingleParticleStateData> AllStates(false); // ordered list, which does not eliminate duplicates

  for (int b=0; b<this->NbrBands; ++b)
    for (int i=0; i<this->NbrStatePerBand; ++i)
      AllStates += SingleParticleStateData(this->GetEnergy(b, i),b,i);

  bands = 0;
  double Energy=0.0;

  if ((verbose)&&(false))
    {
      cout << "Total number of states: "<<AllStates.GetNbrElement()<<endl;
      for (int n=0; n<AllStates.GetNbrElement(); ++n)
	{
	  SingleParticleStateData D=AllStates[n];
	  cout << "State "<<n<<": "<<D.BandIndex<<"_"<<D.StateIndex<<" (E="<<D.Energy<<")\n";
	}
    }

  for (int n=0; n<nbrFermions; ++n)
    {
      SingleParticleStateData D=AllStates[n];
      if (verbose) cout << n+1 <<"th occupied state: "<<D.BandIndex<<"_"<<D.StateIndex<<" (E="<<D.Energy<<")\n";
      Energy+=D.Energy;
      if (D.BandIndex>bands)
	bands=D.BandIndex;
    }
  
  ++bands;
  return Energy;
}


// return the energy of the lowest energy single-particle state
// 
double AbstractTightBindingModel::SingleParticleGroundstateEnergy()
{
  double Minimum = this->GetEnergy(0,0);
  for (int b=0; b<this->NbrBands; ++b)
    for (int i=0; i<this->NbrStatePerBand; ++i)
      if (this->GetEnergy(b, i) < Minimum)
	Minimum = this->GetEnergy(b, i);
  return Minimum;
}


  
  


// write the full band structure information in a binary file
//
// fileName = name of the output file 
// return value = true if no error occured  

bool AbstractTightBindingModel::WriteBandStructure(char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  WriteLittleEndian(File, this->NbrBands);
  WriteLittleEndian(File, this->NbrStatePerBand);
  if (this->Inversion.GetNbrRow() == 0)
  {
      for (int i = 0; i < this->NbrBands; ++i)
          for (int j = 0; j < this->NbrBands; ++j)
          {
              Complex TmpInversion = (i == j);
              WriteLittleEndian(File, TmpInversion);
          }
  }
  else
  {
      for (int i = 0; i < this->NbrBands; ++i)
          for (int j = 0; j < this->NbrBands; ++j)
          {
              Complex TmpInversion = this->Inversion[i][j];
              WriteLittleEndian(File, TmpInversion);
          }
  }
  this->WriteHeader(File);
  for (int kx = 0; kx < this->NbrStatePerBand; ++kx)
    for (int i = 0; i < this->NbrBands; ++i)
      {
	double Tmp = this->GetEnergy(i, kx);
	WriteLittleEndian(File, Tmp);
      }
  if (this->HaveOneBodyBasis() == true)
    {
      for (int kx = 0; kx < this->NbrStatePerBand; ++kx)
	this->GetOneBodyMatrix(kx).WriteMatrix(File);
    }
  File.close();
  return true;
}

// write the full band structure information in an ASCII file
//
// fileName = name of the output file 
// return value = true if no error occured  

bool AbstractTightBindingModel::WriteBandStructureASCII(char* fileName)
{
  ofstream File;
  this->WriteASCIIHeader(File, '#');
  File << "# index";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  for (int i = 0; i < this->NbrBands; ++i)
    for (int j = 0; j < this->NbrBands; ++j)
      File <<  "    U_{" << i << ", " << j << "}";
  File << endl;
  Complex Tmp;
  for (int LinearizedMomentumIndex = 0; LinearizedMomentumIndex < this->NbrStatePerBand; ++LinearizedMomentumIndex)
    {
      File << LinearizedMomentumIndex; 
      for (int i = 0; i < this->NbrBands; ++i)
	File << " " <<  this->GetEnergy(i, LinearizedMomentumIndex);
      for (int i = 0; i < this->NbrBands; ++i)
	for (int j = 0; j < this->NbrBands; ++j)
	  {
	    this->GetOneBodyMatrix(LinearizedMomentumIndex).GetMatrixElement(i, j, Tmp);
	    File <<  "    " << Tmp;
	  }
      File << endl;
    }

  File.close();
  return true;
}

// compute the band structure
//

void AbstractTightBindingModel::ComputeBandStructure()
{
  timeval TotalStartingTime;
  gettimeofday (&TotalStartingTime, 0);
  FTIComputeBandStructureOperation Operation (this);
  Operation.ApplyOperation(this->Architecture);
  timeval TotalEndingTime;
  gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout << "One-body diagonalization done in " << Dt << " s" << endl;
}

