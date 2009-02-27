#include "SphereWithSpinGeneralEnergy.h"
#include "GeneralTools/ConfigurationParser.h"

#include <cmath>
using std::sqrt;
using std::cout;
using std::endl;


double dsqrarg1;   // shift declaration to nrutil.cc to suppress warnings if unused in module that includes nrutil.h
#define DSQR(a) ((dsqrarg1=(a)) == 0.0 ? 0.0 : dsqrarg1*dsqrarg1)


// default constructor
SphereWithSpinGeneralEnergy::SphereWithSpinGeneralEnergy()
{
  this->NbrFlux=0;
}

// constructor
// nbrUp = Number of particles with up spin
// nbrFlux = Number of Flux piercing sphere
// parametersInter = file describing parameters of the interaction (inter-spin)
// parametersIntra = file describing parameters of the interaction (intra-spin)
// parametersIntra2 = file describing parameters of the interaction (intra-spin, other spin-channel)
SphereWithSpinGeneralEnergy::SphereWithSpinGeneralEnergy(int nbrUp, int nbrFlux, const char* parametersInter, const char* parametersIntra, const char* parametersIntra2)
{
  this->Type=AbstractObservable::RealObservable;
  this->NbrUp = nbrUp;
  this->NbrFlux = nbrFlux;
  this->Values = new WeightedRealObservable();
  this->NbrObservations=0;

  ConfigurationParser InteractionDefinition;
  if (InteractionDefinition.Parse(parametersInter) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      exit(-1);
    }
  int TmpNbrParameters;
  if (InteractionDefinition.GetAsDoubleArray("Parameters", ' ', this->CoefficientsInter, TmpNbrParameters) == false)
    {
      cout << "Parameters are not defined or has a wrong value in " << parametersInter << endl;
      exit(-1);
    }
  if (InteractionDefinition.GetAsSingleInteger("NumCoefficients", this->NbrParametersInter) == false)
    {
      cout << "NumCoefficients are not defined or has a wrong value in " << parametersInter << endl;
      exit(-1);
    }  
  if (NbrParametersInter!=TmpNbrParameters)
    {
      cout << "Values not consistent in " << parametersInter << endl;
      exit(-1);
    }
  int TmpNphi;
  if (InteractionDefinition.GetAsSingleInteger("Nphi", TmpNphi) == false)
    {
      cout << "Nphi is not defined or has a wrong value in " << parametersInter << endl;
      exit(-1);
    }
  if (this->NbrFlux==0)
    this->NbrFlux=TmpNphi;
  if (NbrFlux!=TmpNphi)
    {
      cout << "Interspin interactions have inconsistent flux value"<<endl;
      exit(-1);
    }
  this->Radius = sqrt(0.5*(double)TmpNphi); // the radius is also the inverse magnetic length


  // Intra-spin interactions:
  if (InteractionDefinition.Parse(parametersIntra) == false)
    {
      InteractionDefinition.DumpErrors(cout) << endl;
      exit(-1);
    }
  if (InteractionDefinition.GetAsDoubleArray("Parameters", ' ', this->CoefficientsIntra, TmpNbrParameters) == false)
    {
      cout << "Parameters are not defined or has a wrong value in " << parametersIntra << endl;
      exit(-1);
    }
  if (InteractionDefinition.GetAsSingleInteger("NumCoefficients", this->NbrParametersIntra) == false)
    {
      cout << "NumCoefficients are not defined or has a wrong value in " << parametersIntra << endl;
      exit(-1);
    }  
  if (NbrParametersIntra!=TmpNbrParameters)
    {
      cout << "Values not consistent in " << parametersIntra << endl;
      exit(-1);
    }
  if (InteractionDefinition.GetAsSingleInteger("Nphi", TmpNphi) == false)
    {
      cout << "Nphi is not defined or has a wrong value in " << parametersIntra << endl;
      exit(-1);
    }
  if (NbrFlux!=TmpNphi)
    {
      cout << "Intraspin interactions have inconsistent flux value"<<endl;
      exit(-1);
    }

  if (parametersIntra2!=NULL) // have asymmetric interactions?
    {
      this->HaveIntra2=true;
      // Intra-spin interactions:
      if (InteractionDefinition.Parse(parametersIntra2) == false)
	{
	  InteractionDefinition.DumpErrors(cout) << endl;
	  exit(-1);
	}
      int TmpNbrParameters;
      if (InteractionDefinition.GetAsDoubleArray("Parameters", ' ', this->CoefficientsIntra2, TmpNbrParameters) == false)
	{
	  cout << "Parameters are not defined or has a wrong value in " << parametersIntra2 << endl;
	  exit(-1);
	}
      if (InteractionDefinition.GetAsSingleInteger("NumCoefficients", this->NbrParametersIntra2) == false)
	{
	  cout << "NumCoefficients are not defined or has a wrong value in " << parametersIntra2 << endl;
	  exit(-1);
	}  
      if (NbrParametersIntra2!=TmpNbrParameters)
	{
	  cout << "Values not consistent in " << parametersIntra2 << endl;
	  exit(-1);
	}
      if (InteractionDefinition.GetAsSingleInteger("Nphi", TmpNphi) == false)
	{
	  cout << "Nphi is not defined or has a wrong value in " << parametersIntra2 << endl;
	  exit(-1);
	}
      if (NbrFlux!=TmpNphi)
	{
	  cout << "Second  intraspin interactions have inconsistent flux value"<<endl;
	  exit(-1);
	}
    }
  else
    {
      this->HaveIntra2=false;
      CoefficientsIntra2=CoefficientsIntra;
      NbrParametersIntra2=NbrParametersIntra;
    }
    
}


// destructor
SphereWithSpinGeneralEnergy::~SphereWithSpinGeneralEnergy()
{
  if (NbrFlux>0)
    {
      delete Values;
      delete [] CoefficientsInter;
      delete [] CoefficientsIntra;
      if (this->HaveIntra2)
	delete [] CoefficientsIntra2;
    }
}

// call to make an observation
void SphereWithSpinGeneralEnergy::RecordValue(double weight)
{
  int N = this->NbrParticles;
  ++NbrObservations;
  double rst, dij, sum=0.0;
  for (int i=1;i<this->NbrUp;i++)
    {
      for(int j=0;j<i;j++)
	{
	  dij = 2.0*Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	  rst = 1.0/ dij;
	  double p = this->CoefficientsIntra[this->NbrParametersIntra-1];
	  for (int k=this->NbrParametersIntra-2; k>-1; --k)
	    {
	      p=p*dij + this->CoefficientsIntra[k];
	    }
	  rst+=p;
	  sum+=rst;
	}
    }
  for (int i=this->NbrUp+1;i<N;i++)
    {
      for(int j=this->NbrUp;j<i;j++)
	{
	  dij = 2.0*Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	  rst = 1.0/ dij;
	  double p = this->CoefficientsIntra2[this->NbrParametersIntra2-1];
	  for (int k=this->NbrParametersIntra2-2; k>-1; --k)
	    {
	      p=p*dij + this->CoefficientsIntra2[k];
	    }
	  rst+=p;
	  sum+=rst;
	}
    }
  for (int i=0;i<this->NbrUp;i++)
    {
      for(int j=this->NbrUp;j<N;j++)
	{
	  dij = 2.0*Norm(SpinorUCoordinates[i]*SpinorVCoordinates[j]-SpinorUCoordinates[j]*SpinorVCoordinates[i]);
	  rst = 1.0/ dij;
	  double p = this->CoefficientsInter[this->NbrParametersInter-1];
	  for (int k=this->NbrParametersInter-2; k>-1; --k)
	    {
	      p=p*dij + this->CoefficientsInter[k];
	    }
	  rst+=p;
	  sum+=rst;
	}
    }
  this->Values->Observe(sum/Radius, weight);
}

// print legend to the given stream
void SphereWithSpinGeneralEnergy::PrintLegend(std::ostream &output, bool all)
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
void SphereWithSpinGeneralEnergy::PrintStatus(std::ostream &output, bool all)
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
void SphereWithSpinGeneralEnergy::WriteDataFile(std::ostream &output)
{
  output << "#  E  \t err  "<<endl;
  output << this->Values->Average()
	 <<"\t"<<this->Values->ErrorEstimate()<<endl;  
}

// set particle collection that the observable operates on
// system = particle collection
void SphereWithSpinGeneralEnergy::SetParticleCollection(AbstractParticleCollection *system)
{
  if (system->GetCollectionType()!=AbstractParticleCollection::OnSphereCollection)
    {
      cout << "Need a particle collection on the sphere for SphereCoulombEnergy"<<endl;
      exit(1);
    }
  this->System = (ParticleOnSphereCollection*) system;
  this->NbrParticles = System->GetNbrParticles();
  if (this->NbrParticles<this->NbrUp)
    {
      cout << "The number of particles needs to be superiour to the assigned value of up-spins in SphereWithSpinGeneralEnergy" << endl;
      exit(1);
    }
  this->System->GetSpinorCoordinates(SpinorUCoordinates, SpinorVCoordinates);
}
