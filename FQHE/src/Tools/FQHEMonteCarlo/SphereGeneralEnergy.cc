#include "SphereGeneralEnergy.h"
#include "GeneralTools/ConfigurationParser.h"

#include <cmath>
using std::sqrt;
using std::cout;
using std::endl;


double dsqrarg1;   // shift declaration to nrutil.cc to suppress warnings if unused in module that includes nrutil.h
#define DSQR(a) ((dsqrarg1=(a)) == 0.0 ? 0.0 : dsqrarg1*dsqrarg1)


// default constructor
SphereGeneralEnergy::SphereGeneralEnergy()
{
  this->NbrFlux=0;
}

// constructor
// nbrFlux = Number of Flux piercing sphere
// parameters = file describing parameters of the interaction
SphereGeneralEnergy::SphereGeneralEnergy(int nbrFlux, const char* parameters)
{
  this->Type=AbstractObservable::RealObservable;
  this->NbrFlux = nbrFlux;
  this->Values = new WeightedRealObservable();
  this->NbrObservations=0;

  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(parameters) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      exit(-1);
    }
  int TmpNbrParameters;
  if (InteractionDefinition.GetAsDoubleArray("Parameters", ' ', this->Coefficients, TmpNbrParameters) == false)
    {
      cout << "Parameters are not defined or has a wrong value in " << parameters << endl;
      exit(-1);
    }
  if (InteractionDefinition.GetAsSingleInteger("NumCoefficients", this->NbrParameters) == false)
    {
      cout << "NumCoefficients are not defined or has a wrong value in " << parameters << endl;
      exit(-1);
    }  
  if (NbrParameters!=TmpNbrParameters)
    {
      cout << "Values not consistent in " << parameters << endl;
      exit(-1);
    }
  int TmpNphi;
  if (InteractionDefinition.GetAsSingleInteger("Nphi", TmpNphi) == false)
    {
      cout << "Nphi is not defined or has a wrong value in " << parameters << endl;
      exit(-1);
    }
  cout << "Using TmpNphi="<<TmpNphi<<endl;
  this->Radius = sqrt(0.5*(double)TmpNphi); // the radius is also the inverse magnetic length

  if (this->NbrFlux==0)
    this->NbrFlux=TmpNphi;
}


// destructor
SphereGeneralEnergy::~SphereGeneralEnergy()
{
  if (NbrFlux>0)
    {
      delete Values;
    }
}

// call to make an observation
void SphereGeneralEnergy::RecordValue(double weight)
{
  int N = this->NbrParticles;
  ++NbrObservations;
  double rst, dij, sum=0.0;
  for (int i=1;i<N;i++)
    {
      for(int j=0;j<i;j++)
	{
	  dij = 2.0*Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	  rst = this->Coefficients[0]/ dij;
	  double p = this->Coefficients[this->NbrParameters-1];
	  for (int k=this->NbrParameters-2; k>0; --k)
	    {
	      p=p*dij + this->Coefficients[k];
	    }
	  rst+=p;
	  sum+=rst;
	}
    }
  this->Values->Observe(sum/Radius, weight);
}

// old version:
// // call to make an observation
// void SphereGeneralEnergy::RecordValue(double weight)
// {
//   int N = this->NbrParticles;
//   ++NbrObservations;
//   double rst, dij, sum=0.0;
//   for (int i=1;i<N;i++)
//     {
//       for(int j=0;j<i;j++)
// 	{
// 	  dij = 2.0*Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
// 	  rst = 1.0/ dij;
// 	  double p = this->Coefficients[this->NbrParameters-1];
// 	  for (int k=this->NbrParameters-2; k>-1; --k)
// 	    {
// 	      p=p*dij + this->Coefficients[k];
// 	    }
// 	  rst+=p;
// 	  sum+=rst;
// 	}
//     }
//   this->Values->Observe(sum/Radius, weight);
// }

// print legend to the given stream
void SphereGeneralEnergy::PrintLegend(std::ostream &output, bool all)
{
  if (all)
    {
      output << "E\t+/-";
    }
  else
    {
      output << "E\t+/-";
    }
}

// print status to the given stream
void SphereGeneralEnergy::PrintStatus(std::ostream &output, bool all)
{
  if (NbrObservations>0)
    {
      if (all)
	{
	  output << this->Values->Average()<<"\t"<<this->Values->ErrorEstimate();
	}
      else
	{
	  int tmp=output.precision();
	  output.precision(6);
	  output << this->Values->Average()<<"\t"<<this->Values->ErrorEstimate();
	  output.precision(tmp);
	}
    }
}

// print formatted data suitable for plotting
// ouput = the target stream
void SphereGeneralEnergy::WriteDataFile(std::ostream &output)
{
  output << "#  E  \t err  "<<endl;
  output << this->Values->Average()
	 <<"\t"<<this->Values->ErrorEstimate()<<endl;  
}

// set particle collection that the observable operates on
// system = particle collection
void SphereGeneralEnergy::SetParticleCollection(AbstractParticleCollection *system)
{
  if (system->GetCollectionType()!=AbstractParticleCollection::OnSphereCollection)
    {
      cout << "Need a particle collection on the sphere for SphereCoulombEnergy"<<endl;
      exit(1);
    }
  this->System = (ParticleOnSphereCollection*) system;
  this->NbrParticles = System->GetNbrParticles();
  this->System->GetSpinorCoordinates(SpinorUCoordinates, SpinorVCoordinates);
}
