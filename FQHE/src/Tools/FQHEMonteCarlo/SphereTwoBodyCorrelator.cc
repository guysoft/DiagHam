////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//        class for a binned correlation function on the sphere geometry      //
//                                                                            //
//                        last modification : 19/10/2009                      //
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
#include "SphereTwoBodyCorrelator.h"
#include <iostream>

using std::ios;
using std::ios_base;
using std::endl;


SphereTwoBodyCorrelator::SphereTwoBodyCorrelator()
{
  this->Bins=0;
}

// constructor
// resolution = total number of bins
// highres = number of points in high resolution interval at small r
// range =  ranger over which high resolution is implemented
SphereTwoBodyCorrelator::SphereTwoBodyCorrelator(int nbrFlux, int resolution, int highres, int range, bool printLength)
{
  this->Bins=resolution+highres-range+1;
  this->Resolution=resolution;
  this->Highres=highres;
  this->Range=range;
  this->Measures = 0.0;
  this->Correlations=new double[Bins];
  for (int j=0;j<Bins;++j)
    Correlations[j]=0.0;
  this->NbrFlux=nbrFlux;
  this->NbrParticles=0;
  this->Radius = sqrt(0.5*(double)NbrFlux); // the radius is also the inverse magnetic length
  this->PrintLength=printLength;
}
  
// destructor
SphereTwoBodyCorrelator::~SphereTwoBodyCorrelator()
{
  if (Bins!=0)
    delete [] Correlations;
}

// call to make an observation
// weight = relative weight of this sample
void SphereTwoBodyCorrelator::RecordValue(double weight)
{
  double Rij,phi,tmp;
  int index;
  this->Measures+=weight;
  for( int i=1;i<this->NbrParticles;++i)
    {
      for(int j=0; j<i;++j)
	{
	  // set Rij=sin(\theta_ij/2)
	  Rij=Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);  
	  phi=2.0*asin(Rij);
	  index= (int)((1.0-(tmp=cos(phi)))/2.0*Resolution);
	  if (index>=Range)
	    this->Correlations[index+Highres-Range+1]+=weight;  
	  else
	    {
	      index= (int)((1.0-tmp)/2.0*Highres*Resolution/Range);
	      this->Correlations[index+1]+=weight;
	    }
	}
    }
}

// print legend to the given stream
// all = flag indicating whether to print all, or shortened information
void SphereTwoBodyCorrelator::PrintLegend(std::ostream &output, bool all)
{
  output << "# r\tg(r)"<<endl;
}

// print status to the given stream
// all = flag indicating whether to print all, or shortened information
void SphereTwoBodyCorrelator::PrintStatus(std::ostream &output, bool all)
{
  // no action, for now
}

// print formatted data suitable for plotting
// ouput = the target stream
void SphereTwoBodyCorrelator::WriteDataFile(std::ostream &output)
{
  if (output.flags() & ios_base::binary)
    {
      // write as binary file
      this->WriteBinaryData(output);
    }
  else
    {
      // write as textfile

      double Units;
      if (this->PrintLength)
	{
	  Units=this->Radius;
	  output << "# r\tg(r)\n";
	}
      else
	{
	  Units=1.0;
	  output << "# phi\tg(phi)\n";
	}
      double Normalization=Measures/(Resolution*(double)Highres/Range)*NbrParticles*NbrParticles;
      for (int i=0;i<=Highres;i++)
	output << Units*acos(1-(double)i/(Resolution*(double)Highres/Range)*2.0)<<"\t"
	       << this->Correlations[i]/Normalization << endl;
      Normalization=Measures/Resolution*NbrParticles*NbrParticles;
      for (int i=1;i<=Resolution-Range; ++i)
	output << Units*acos(1-(double)(i+Range)/Resolution*2.0) <<"\t"
	       << this->Correlations[i]/Normalization << endl;
    }
}


// write binary data 
// ouput = the target stream
void SphereTwoBodyCorrelator::WriteBinaryData(std::ostream &output)
{
  std::cout << "Need to implement binary write"<<endl;
}

// set particle collection that the observable operates on
// system = particle collection
void SphereTwoBodyCorrelator::SetParticleCollection(AbstractParticleCollection *system)
{
  this->System = (ParticleOnSphereCollection*) system;
  this->NbrParticles = System->GetNbrParticles();
  this->System->GetSpinorCoordinates(SpinorUCoordinates, SpinorVCoordinates);
}

