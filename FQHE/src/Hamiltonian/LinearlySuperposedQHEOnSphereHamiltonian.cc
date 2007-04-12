#include "LinearlySuperposedQHEOnSphereHamiltonian.h"
#include <iostream>

using std::cout;
using std::endl;

// creates a Hamiltonian that consists of a superposition of Number Hamiltonians with an array of coefficients and
// Hamiltonians given as a list
LinearlySuperposedQHEOnSphereHamiltonian::LinearlySuperposedQHEOnSphereHamiltonian( int number, double *coefficients, AbstractQHEOnSphereHamiltonian **hamiltonians)
{
  this->NumHamiltonians=number;
  this->LinearCoefficients=coefficients;
  this->ListHamiltonians=hamiltonians;
  
  // Hilbert space associated to the system
  this->Particles=hamiltonians[0]->Particles;
  for (int n=1; n<NumHamiltonians;++n)
    if (hamiltonians[n]->Particles!=this->Particles)
      {
	cout << "Hamiltonians in LinearlySuperposedQHEOnSphereHamiltonian should be constructed on the same Hilbert space" << endl;
	exit(5);
      }

  this->VectorSum.Resize(this->Particles->GetHilbertSpaceDimension());
  this->VectorSummand.Resize(this->Particles->GetHilbertSpaceDimension());
  // number of particles
  this-> NbrParticles=hamiltonians[0]->NbrParticles;
  this->LzMax=hamiltonians[0]->LzMax;
  this->NbrLzValue=hamiltonians[0]->NbrLzValue;
}

LinearlySuperposedQHEOnSphereHamiltonian::~LinearlySuperposedQHEOnSphereHamiltonian()
{
  for (int n=0; n<NumHamiltonians;++n) delete ListHamiltonians[n];
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void LinearlySuperposedQHEOnSphereHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  for (int n=0; n<NumHamiltonians;++n) ListHamiltonians[n]->SetHilbertSpace(hilbertSpace);
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void LinearlySuperposedQHEOnSphereHamiltonian::ShiftHamiltonian (double shift)
{
  ListHamiltonians[0]->ShiftHamiltonian (shift);
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored
RealVector& LinearlySuperposedQHEOnSphereHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)
{
  for (int n=0; n<NumHamiltonians;++n)
    {
      this->VectorSummand.ClearVector();
      ListHamiltonians[n]->LowLevelAddMultiply(vSource, VectorSummand, firstComponent, nbrComponent);
      VectorSummand*=LinearCoefficients[n];
      vDestination+=VectorSummand;
    }
  return vDestination;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& LinearlySuperposedQHEOnSphereHamiltonian::LowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
						   int firstComponent, int nbrComponent)
{
  for (int n=0; n<NumHamiltonians;++n)
    {
      this->VectorSummand.ClearVector();
      ListHamiltonians[n]->LowLevelAddMultiplyPartialFastMultiply(vSource, VectorSummand, firstComponent, nbrComponent);
      VectorSummand*=LinearCoefficients[n];
      vDestination+=VectorSummand;
    }
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using disk storage option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& LinearlySuperposedQHEOnSphereHamiltonian::LowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
					   int firstComponent, int nbrComponent)
{
  for (int n=0; n<NumHamiltonians;++n)
    {
      this->VectorSummand.ClearVector();
      ListHamiltonians[n]->LowLevelAddMultiplyDiskStorage(vSource, VectorSummand, firstComponent, nbrComponent);
      VectorSummand*=LinearCoefficients[n];
      vDestination+=VectorSummand;
    }
  return vDestination;
}



// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* LinearlySuperposedQHEOnSphereHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
					int firstComponent, int nbrComponent)
{
  for (int n=0; n<NumHamiltonians;++n)
    {
      for (int v=0;v<nbrVectors;++v)
	{
	  this->VectorSummand.ClearVector();
	  ListHamiltonians[n]->LowLevelAddMultiply(vSources[v], VectorSummand, firstComponent, nbrComponent);
	  VectorSummand*=LinearCoefficients[n];
	  vDestinations[v]+=VectorSummand;
	}
    }
  return vDestinations;
}


  
// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* LinearlySuperposedQHEOnSphereHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
							   int firstComponent, int nbrComponent)
{
  for (int n=0; n<NumHamiltonians;++n)
      {
      for (int v=0;v<nbrVectors;++v)
	{
	  this->VectorSummand.ClearVector();
	  ListHamiltonians[n]->LowLevelAddMultiplyPartialFastMultiply(vSources[v], VectorSummand, firstComponent, nbrComponent);
	  VectorSummand*=LinearCoefficients[n];
	  vDestinations[v]+=VectorSummand;
	}
    }
  return vDestinations;
}
  
// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using disk storage option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored
  
RealVector* LinearlySuperposedQHEOnSphereHamiltonian::LowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
						   int firstComponent, int nbrComponent)
{
  for (int n=0; n<NumHamiltonians;++n)
    {
      for (int v=0;v<nbrVectors;++v)
	{
	  this->VectorSummand.ClearVector();
	  ListHamiltonians[n]->LowLevelAddMultiplyDiskStorage(vSources[v], VectorSummand, firstComponent, nbrComponent);
	  VectorSummand*=LinearCoefficients[n];
	  vDestinations[v]+=VectorSummand;
	}
    }
  return vDestinations;
}

// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long LinearlySuperposedQHEOnSphereHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  long rst=0;
  long tmp=0;
  long remainingMemoryAllowance=allowedMemory;
  for (int n=0; n<NumHamiltonians;++n)
    {
      tmp+=ListHamiltonians[n]->FastMultiplicationMemory(remainingMemoryAllowance);
      remainingMemoryAllowance-=tmp;
      rst+=tmp;
    }
  return rst;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element
  
long LinearlySuperposedQHEOnSphereHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  long rst=0;
  for (int n=0; n<NumHamiltonians;++n)
    rst+=ListHamiltonians[n]->PartialFastMultiplicationMemory(firstComponent, lastComponent);
  return rst;
}
  

// enable fast multiplication algorithm
//

void LinearlySuperposedQHEOnSphereHamiltonian::EnableFastMultiplication()
{
  for (int n=0; n<NumHamiltonians;++n) ListHamiltonians[n]->EnableFastMultiplication();
}

// enable fast multiplication algorithm using on disk cache 
//
// fileName = prefix of the name of the file where temporary matrix elements will be stored
  
void LinearlySuperposedQHEOnSphereHamiltonian::EnableFastMultiplicationWithDiskStorage(char* fileName)
{
  char *tmpString = new char[strlen(fileName)+10];
  for (int n=0; n<NumHamiltonians;++n)
    {
      sprintf(tmpString,"%s-H%d",fileName,n);
      ListHamiltonians[n]->EnableFastMultiplicationWithDiskStorage(tmpString);
    }
  delete [] tmpString;
}

// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs
  
bool LinearlySuperposedQHEOnSphereHamiltonian::SavePrecalculation (char* fileName)
{
  char *tmpString = new char[strlen(fileName)+10];
  bool rst=true;
  for (int n=0; n<NumHamiltonians;++n)
    {
      sprintf(tmpString,"%s-H%d",fileName,n);
      rst=rst && ListHamiltonians[n]->SavePrecalculation(tmpString);
    }
  delete [] tmpString;
  return rst;
}

// load precalculations from a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be read
// return value = true if no error occurs

bool LinearlySuperposedQHEOnSphereHamiltonian::LoadPrecalculation (char* fileName)
{
  char *tmpString = new char[strlen(fileName)+10];
  bool rst=true;
  for (int n=0; n<NumHamiltonians;++n)
    {
      sprintf(tmpString,"%s-H%d",fileName,n);
      rst = (rst && ListHamiltonians[n]->LoadPrecalculation(tmpString));
    }
  delete [] tmpString;
  return rst;
}


// evaluate all interaction factors
//   
void LinearlySuperposedQHEOnSphereHamiltonian::EvaluateInteractionFactors()
{
  for (int n=0; n<NumHamiltonians;++n) ListHamiltonians[n]->EvaluateInteractionFactors();
}
