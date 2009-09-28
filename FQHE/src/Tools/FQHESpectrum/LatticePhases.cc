////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2009 Nicolas Regnault                  //
//                         class author: Gunnar M�ller                        //
//                                                                            //
//                                                                            //
//  class that defines sites and nearest neighbor relationships on a lattice  //
//                                                                            //
//                        last modification : 05/26/2009                      //
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

#include "LatticePhases.h"
#include "Options/Options.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"

#include <iostream>
using std::cout;
using std::endl;

#include <cstdlib>
#include <cstring>
#include <cmath>
using std::fabs;

// generate the object using options from Option Manager
//
LatticePhases::LatticePhases()
{
  if (LatticePhases::Options==NULL)
    {
      cout << "Define the OptionManager, first, before creating any LatticePhases"<<endl;
      exit(1);
    }
  
  ConfigurationParser LatticeDefinition;
  if (LatticeDefinition.Parse(this->Options->GetString("lattice-definition")) == false)
    {
      LatticeDefinition.DumpErrors(cout) << endl;
      exit(-1);
    }
  if (LatticeDefinition["Descriptor"] == NULL)
    {
      cout << "Attention, 'Descriptor' is not defined, unnamed lattice geometry!" << endl;
      Descriptor = new char[10];
      sprintf(Descriptor,"unnamed");
    }
  else
    {
      this->Descriptor = new char[strlen(LatticeDefinition["Descriptor"])+1];
      strcpy(this->Descriptor, LatticeDefinition["Descriptor"]);
    }
  if ((LatticeDefinition.GetAsSingleInteger("NbrSites", NbrSitesPerCell) == false) || (NbrSitesPerCell <= 0))
    {
      cout << "NbrSites is not defined or has invalid value" << endl;
      exit(-1);
    }
  if ((LatticeDefinition.GetAsSingleInteger("Dimension", Dimension) == false) || (Dimension <= 0))
    {
      cout << "Dimension is not defined or has invalid value" << endl;
      exit(-1);
    }
  int TmpDimension;
  PeriodicRep = LatticePhases::Options->GetIntegers("cells-repeat",TmpDimension);
  if (TmpDimension==0)
    if ((LatticeDefinition.GetAsIntegerArray("PeriodicRepeat",',',PeriodicRep, TmpDimension) == false))
      {
	cout << "PeriodicRepeat is not defined or has invalid value, simulating a single unit cell" << endl;
	PeriodicRep = new int[Dimension];
	TmpDimension=Dimension;
	for (int i=0; i<Dimension; ++i) PeriodicRep[i]=1;
      }
  if (TmpDimension!=Dimension)
    {
      cout << "PeriodicRepeat does not have the right number of components (separator: ',')" << endl;
      exit(-1);
    }
  NbrSites=NbrSitesPerCell;
  for (int i=0; i<Dimension; ++i)
    NbrSites*=PeriodicRep[i];
  LatticeVectors.Resize(Dimension, Dimension);
  char *FieldName = new char[255];
  for (int i=0; i<Dimension; ++i)
    {
      sprintf(FieldName,"LatticeVector%d",i);
      int NbrComponents;
      double *Components;
      if (LatticeDefinition.GetAsDoubleArray(FieldName, ',', Components, NbrComponents) == false)
	{
	  cout << "error while parsing "<<FieldName<< " in " << this->Options->GetString("lattice-definition") << endl;
	  exit(-1);     
	}
      if (Dimension!=NbrComponents)
	{
	  cout << "Lattice Vectors need to have the dimension of the lattice!"<<endl;
	  exit(-1);
	}
      for (int j=0; j<Dimension; ++j)
	{
	  LatticeVectors[i][j]=Components[j];
	}
      delete [] Components;
    }
  double Rescale=1.0;
  if (Options->GetBoolean("normalize-lattice"))
    {
      double Area = LatticeVectors.Determinant();
      Rescale = pow((double)NbrSitesPerCell/Area,1.0/(double)Dimension);
      for (int i=0; i<Dimension; ++i)
	LatticeVectors[i]*=Rescale;
      Area = LatticeVectors.Determinant();
    }
  
  
  SubLatticeVectors.Resize(Dimension, NbrSitesPerCell);
  for (int i=0; i<NbrSitesPerCell; ++i)
    {
      sprintf(FieldName,"SubLatticeVector%d",i);
      int NbrComponents;
      double *Components;
      if (LatticeDefinition.GetAsDoubleArray(FieldName, ',', Components, NbrComponents) == false)
	{
	  cout << "error while parsing "<<FieldName<< " in " << this->Options->GetString("lattice-definition") << endl;
	  exit(-1);     
	}      
      if (Dimension!=NbrComponents)
	{
	  cout << "SubLattice Vectors need to have the dimension of the lattice!"<<endl;
	  exit(-1);
	}
      for (int j=0; j<Dimension; ++j)
	{
	  SubLatticeVectors[i][j]=Components[j];
	}
      delete [] Components;
      if (fabs(Rescale-1.0)>1e-10)
	SubLatticeVectors[i]*=Rescale;
    }
  
  // determine connectivity, and tunnelling phases
  // two alternate methods: indicate individual tunnelling phases, or give a gauge choice
  this->HaveGauge=false;
  if ((LatticeDefinition["UseGauge"]!=NULL)&&
      ( (strcmp(LatticeDefinition["UseGauge"],"yes")>0) || (strcmp(LatticeDefinition["UseGauge"],"YES")>0)
	|| (strcmp(LatticeDefinition["UseGauge"],"true")>0) || (strcmp(LatticeDefinition["UseGauge"],"TRUE")>0) ))
    {
      this->HaveGauge=true;
      if (LatticeDefinition.GetAsSingleDouble ("GaugeAxx", this->GaugeAxx)==false)
	this->GaugeAxx=0.0;
      if (LatticeDefinition.GetAsSingleDouble ("GaugeAxy", this->GaugeAxy)==false)
	this->GaugeAxy=0.0;
      if (LatticeDefinition.GetAsSingleDouble ("GaugeAyx", this->GaugeAyx)==false)
	this->GaugeAyx=0.0;
      if (LatticeDefinition.GetAsSingleDouble ("GaugeAyy", this->GaugeAyy)==false)
	this->GaugeAyy=0.0;
      this->AbsBField=GaugeAxy-GaugeAyx;
      if (fabs(AbsBField)<1e-14)
	{
	  cout << "If using gauge mode, the magnetic field B=GaugeAxy - GaugeAyx needs to be non-zero" << endl;
	  exit(-1);
	}       
      if (this->Options->GetBoolean("normalize-lattice"))
	{
	  this->AbsBField=1.0/this->AbsBField;
	  this->GaugeAxx*=this->AbsBField;
	  this->GaugeAxy*=this->AbsBField;
	  this->GaugeAyx*=this->AbsBField;
	  this->GaugeAyy*=this->AbsBField;
	  this->AbsBField=1.0;
	}
      cout << "Gauge used : A= ("<<GaugeAxx<<"*x +"<<GaugeAyx<<"*y) e_x"
	   << " + ("<<GaugeAxy << "*x +"<<GaugeAyy<<"*y) e_y, "
	   << "field strength B="<<this->AbsBField<<endl;
    }
	cout << "Attention, the code is currently not functional for gauges involving both Axy and Ayx!"<<endl;
	
  for (int i=0; i<Dimension; ++i)
    cout << "LatticeVector["<<i<<"]="<<endl<<LatticeVectors[i];
  for (int i=0; i<NbrSitesPerCell; ++i)
    cout << "SubLatticeVector["<<i<<"]="<<endl<<SubLatticeVectors[i];
  char ***NeighborString;
  int NbrPairs;
  int *NbrValues;
  RealSymmetricMatrix NeighborsInCellMatrix(NbrSitesPerCell, true);
  RealAntisymmetricMatrix TunnellingPhaseMatrix(NbrSitesPerCell, true);
  if (LatticeDefinition["NeighborsInCell"]!=NULL)
    {
      if (LatticeDefinition.GetAsStringMultipleArray ("NeighborsInCell", '|', ',', NeighborString, NbrPairs, NbrValues)==false)
	{
	  cout << "error while parsing NeighborsInCell in " << this->Options->GetString("lattice-definition") << endl;
	  exit(-1);
	}
      for (int p=0; p<NbrPairs; ++p)
	{
	  int s1, s2;
	  if (NbrValues[p]<2)
	    {
	      cout << "error while decoding NeighborsInCell in " << this->Options->GetString("lattice-definition") << endl;
	      cout << "Indicate paires of neighboring sites separated by commas and different pairs by bars: "
		   << "NeighborsInCell = s1,s2[,phaseA12=0.0] | s3, s4[,phaseA34=0.0] | ..."
		   << "Phases can either be indicated explitly, or will be deduced from gauge if defined"
		   << " or assumed to be one, otherwise"<<endl;
	      
	      exit(-1);
	    }
	  s1=strtod(NeighborString[p][0], NULL);
	  s2=strtod(NeighborString[p][1], NULL);
	  if ((s1<0)||(s1>=NbrSitesPerCell))
	    {
	      cout << "Attention: pair index "<<s1<<" out of range, ignoring pair ("<<s1<<", "<<s2<<")."<<endl;
	    }
	  else
	    {
	      if ((s2<0)||(s2>=NbrSitesPerCell))
		cout << "Attention: pair index "<<s2<<" out of range, ignoring pair ("<<s1<<", "<<s2<<")."<<endl;
	      else
		{
		  NeighborsInCellMatrix(s1,s2)=1.0;
		  // determine tunnelling phase
		  if (NbrValues[p]<3)
		    {
		      TunnellingPhaseMatrix(s1,s2)=0.0;
		    }
		  else		    
		    TunnellingPhaseMatrix(s1,s2)=strtod(NeighborString[p][2], NULL);
		}
	    }
	}
      for (int i=0; i<NbrPairs; ++i)
	{
	  for (int j=0; j<NbrValues[i]; ++j)
	    delete [] NeighborString[i][j];
	  delete [] NeighborString[i];
	}
      delete [] NbrValues;
      delete [] NeighborString;
    }
  cout << "NeighborsInCell="<<endl<<NeighborsInCellMatrix;
  if (LatticeDefinition.GetAsStringMultipleArray ("NeighborCells", '|', ',', NeighborString, NbrPairs, NbrValues)==false)
    {
      cout << "error while parsing NeighborCells in " << this->Options->GetString("lattice-definition") << endl;
      exit(-1);
    }
  this->NbrNeighborCells = NbrPairs;
  this->NeighborCells = new int*[NbrPairs];
  for (int p=0; p<NbrNeighborCells; ++p)
    {
      if (NbrValues[p]!=Dimension)
	{
	  cout << "error while decoding NeighborCells in " << this->Options->GetString("lattice-definition") << endl;
	  cout << "Indicate coordinates of neighboring cells as vectors of length Dimension separated by commas"<<endl
	       << "Separate multiple neighboring cells by bars: "
	       << "NeighborCells = v_11,...,v_1d | ... | v_k1, v_kd | ..."<<endl;
	  exit(-1);
	}
      this->NeighborCells[p] = new int[Dimension];
      for (int i=0; i<Dimension; ++i)
	this->NeighborCells[p][i] = strtod(NeighborString[p][i], NULL);
    }
  for (int i=0; i<NbrPairs; ++i)
    {
      for (int j=0; j<NbrValues[i]; ++j)
	delete [] NeighborString[i][j];
      delete [] NeighborString[i];
    }
  delete [] NeighborString;
  delete [] NbrValues;
  RealMatrix **NeighborsAcrossBoundary = new RealMatrix*[NbrNeighborCells];
  RealMatrix **PhasesAcrossBoundary = new RealMatrix*[NbrNeighborCells];
  for (int d=0; d<NbrNeighborCells; ++d)
    {
      sprintf(FieldName,"NeighborsAcrossBoundary%d",NeighborCells[d][0]);
      for (int i=1; i<Dimension; ++i)
	sprintf(FieldName,"%s_%d", FieldName, NeighborCells[d][i]);
      if (LatticeDefinition.GetAsStringMultipleArray (FieldName, '|', ',', NeighborString, NbrPairs, NbrValues)==false)
	{
	  cout << "could not parse "<<FieldName<<" in " << this->Options->GetString("lattice-definition")
	       << ", no connections added."<< endl;
	}
      NeighborsAcrossBoundary[d] = new RealMatrix(NbrSitesPerCell, NbrSitesPerCell, true);
      PhasesAcrossBoundary[d] = new RealMatrix(NbrSitesPerCell, NbrSitesPerCell, true);
      for (int p=0; p<NbrPairs; ++p)
	{
	  int s1, s2;
	  if (NbrValues[p]!=2)
	    {
	      cout << "error while decoding "<<FieldName<<" in " << this->Options->GetString("lattice-definition") << endl;
	      cout << "Indicate paires of neighboring sites separated by commas and different pairs by bars: "
		   << FieldName <<" = s1,s2[,phaseA12] | s3, s4 [,phaseA34] | ..."
		   << "Phases can either be indicated explitly, or will be deduced from gauge if defined"
		   << " or assumed to be one, otherwise"<<endl;
	      exit(-1);
	    }
	  s1=strtod(NeighborString[p][0], NULL);
	  s2=strtod(NeighborString[p][1], NULL);
	  if ((s1<0)||(s1>=NbrSitesPerCell))
	    {
	      cout << "Attention: pair index "<<s1<<" out of range, ignoring pair ("<<s1<<", "<<s2<<")."<<endl;
	    }
	  else
	    {
	      if ((s2<0)||(s2>=NbrSitesPerCell))
		cout << "Attention: pair index "<<s2<<" out of range, ignoring pair ("<<s1<<", "<<s2<<")."<<endl;
	      else
		{
		  NeighborsAcrossBoundary[d]->SetMatrixElement(s1, s2, 1.0);
		  // determine tunnelling phase
		  if (NbrValues[p]<3)
		    {
		      PhasesAcrossBoundary[d]->SetMatrixElement(s1, s2, 0.0);
		    }
		  else
		    PhasesAcrossBoundary[d]->SetMatrixElement(s1,s2,strtod(NeighborString[p][2], NULL));
		}
	    }
	}
      for (int i=0; i<NbrPairs; ++i)
	{
	  for (int j=0; j<NbrValues[i]; ++j)
	    delete [] NeighborString[i][j];
	  delete [] NeighborString[i];
	}
      delete [] NbrValues;
      delete [] NeighborString;
      cout << FieldName<<"="<<endl<<*NeighborsAcrossBoundary[d];
    }  
  
  this->Neighbors = new int*[NbrSites];
  this->TunnellingPhases = new double*[NbrSites];
  this->NbrNeighbors = new int[NbrSites];
  for (int i=0; i<NbrSites; ++i)
    this->NbrNeighbors[i] = 0;
  int *TmpNeighbors = new int[NbrNeighborCells*NbrSites];
  double *TmpPhases = new double[NbrNeighborCells*NbrSites];
  
  this->NbrCells = this->PeriodicRep[0];
  for (int d=1; d<Dimension; ++d)
    this->NbrCells *= this->PeriodicRep[d];

  int *CellCoordinates = new int[Dimension];
  int *CellCoordinates2 = new int[Dimension];
  int *Translation3 = new int[Dimension];
  for (int c=0; c<NbrCells; ++c)
    {
      this->GetCellCoordinates(c, CellCoordinates);
      cout << "Cell "<<c<<":"<< CellCoordinates[0]<<", "<<CellCoordinates[1]<<endl;
      int Site1, Site2, Site3;
      for (int i=0; i<NbrSitesPerCell; ++i)
	{
	  Site1 = this->GetSiteNumber(c, i);
	  cout << "Site 1="<<Site1<<endl;
	  for (int j=0; j<NbrSitesPerCell; ++j)
	    {
	      Site2 = this->GetSiteNumber(c, j);
	      
	      if (NeighborsInCellMatrix(i,j)>0.0)
		{		  
		  TmpNeighbors[NbrNeighbors[Site1]]=Site2;
		  if (this->HaveGauge)
		    TmpPhases[NbrNeighbors[Site1]] = GetTunnellingPhaseFromGauge(Site1, Site2);
		  else
		    TmpPhases[NbrNeighbors[Site1]] = TunnellingPhaseMatrix(Site1,Site2);
		  cout << "Neighbors "<<Site1<<"->"<<Site2<<" with phase "<<TmpPhases[NbrNeighbors[Site1]]<<endl;
		  ++NbrNeighbors[Site1];
		}
	      for (int d=0; d<NbrNeighborCells; ++d)
		{
		  if ((*(NeighborsAcrossBoundary[d]))(i,j)!=0.0)
		    {
		      for (int k=0; k<Dimension; ++k)
			{
			  CellCoordinates2[k]=CellCoordinates[k]+NeighborCells[d][k];
			  cout << "CellCoordinates["<<k<<"]="<<CellCoordinates[k]<< ", "
			       << "NeighborCells["<<d<<", "<<k<<"]="<<NeighborCells[d][k]<< ", "
			       << "CellCoordinates2["<<k<<"]="<<CellCoordinates2[k]<<endl;
			}
		      Site3 = this->GetSiteNumber(CellCoordinates2, j, Translation3);
		      cout << "Translation3= ["<<Translation3[0]<<", "<<Translation3[1]<<"]"<<endl;
		      TmpNeighbors[NbrNeighbors[Site1]]=Site3;
		      if (this->HaveGauge)
			{
			  cout << "Evaluating translation phase:"<<endl;
			  TmpPhases[NbrNeighbors[Site1]] = GetTunnellingPhaseFromGauge(Site1, Site3, Translation3);
			}
		      else
			TmpPhases[NbrNeighbors[Site1]]=(*PhasesAcrossBoundary[d])(Site1,Site3);
		      cout << "additional neighbors "<<Site1<<"->"<<Site3<<" from NeigborCell "<<d<<" at "<<
			CellCoordinates2[0]<<", "<<CellCoordinates2[1]<<", "<<j<<" : Site 3="<<Site3
			   <<" with translation "<<Translation3[0];
		      for (int i=1; i<Dimension; ++i) cout << " "<<Translation3[i];
		      cout << " and phase "<<TmpPhases[NbrNeighbors[Site1]]<<endl;
		      
		      ++NbrNeighbors[Site1];
		    }
		}	      
	    }
	  if (NbrNeighbors[Site1]>0)
	    {
	      Neighbors[Site1] = new int[NbrNeighbors[Site1]];
	      TunnellingPhases[Site1] = new double[NbrNeighbors[Site1]];
	    }
	  else Neighbors[Site1] = NULL;
	  for (int k=0; k<NbrNeighbors[Site1]; ++k)
	    {
	      Neighbors[Site1][k] = TmpNeighbors[k];
	      TunnellingPhases[Site1][k] = TmpPhases[k];
	    }
	}
    }

  cout << "LatticePhases created"<<endl;

  for (int d=0; d<NbrNeighborCells; ++d)
    {
      delete NeighborsAcrossBoundary[d];
      delete PhasesAcrossBoundary[d];
    }
  delete [] NeighborsAcrossBoundary;
  delete [] PhasesAcrossBoundary;
  delete [] CellCoordinates;
  delete [] CellCoordinates2;
  delete [] Translation3;
  delete [] FieldName;
  delete [] TmpNeighbors;
  delete [] TmpPhases;
}

// destructor
//
LatticePhases::~LatticePhases()
{
  if (NbrSites!=0)
    {
      delete [] this->PeriodicRep;
      for (int i=0; i<NbrSites; ++i)
	{
	  if (this->Neighbors[i]!=NULL)
	    delete [] this->Neighbors[i];
	  if (this->TunnellingPhases[i]!=NULL)
	    delete [] this->TunnellingPhases[i];
	}      
      delete [] this->Neighbors;
      delete [] this->TunnellingPhases;
      delete [] this->NbrNeighbors;
      for (int i=0; i<NbrNeighborCells; ++i)
	delete [] this->NeighborCells[i];
      delete [] this->NeighborCells;
      delete [] this->Descriptor;
    }
}


// get cell coordinates given the number of the unit cell
// nbrCell = cell to be looked up
// cellCoordinates = resulting coordinates, has to be reserved prior to call
void LatticePhases::GetCellCoordinates(int nbrCell, int *cellCoordinates)
{
  int Divisor=1;
  while (nbrCell<0) nbrCell+=NbrCells;
  for (int i=0; i<Dimension; ++i)
    {
      cellCoordinates[i] = (nbrCell/Divisor)%this->PeriodicRep[i];
      Divisor*=this->PeriodicRep[i];
    }
}

// get cell coordinates given the number of the unit cell
// nbrSite = cell to be looked up
// cellCoordinates = resulting coordinates, has to be reserved prior to call
// sublattice = resulting sublattice
void LatticePhases::GetSiteCoordinates(int nbrSite, int *cellCoordinates, int &sublattice)
{
  while (nbrSite<0) nbrSite+=NbrSites;
  sublattice = nbrSite%NbrSitesPerCell;
  int Divisor=NbrSitesPerCell;
  for (int i=0; i<Dimension-1; ++i)
    {
      cellCoordinates[i] = (nbrSite/Divisor)%this->PeriodicRep[i];
      Divisor*=this->PeriodicRep[i];
    }
  cellCoordinates[Dimension-1] = (nbrSite/Divisor)%this->PeriodicRep[Dimension-1];
}


// get number of a site in cell nbrCell
// cellCoordinates = coordinates of cell to be addressed
// sublattice = sublattice index
int LatticePhases::GetSiteNumber(int *cellCoordinates, int sublattice)
{
  int Result=this->Periodize(cellCoordinates[Dimension-1],Dimension-1);
  for (int i=Dimension-2; i>-1; --i)
    {
      Result*=this->PeriodicRep[i];
      Result+=this->Periodize(cellCoordinates[i],i);
    }
  Result*=NbrSitesPerCell;
  Result+=sublattice;
  return Result%NbrSites;
}

// get number of a site in cell nbrCell, and return translation vector back into the simulation cell
// cellCoordinates = coordinates of cell to be addressed
// sublattice = sublattice index
// translation = vector of tranlation back into simulation cell
int LatticePhases::GetSiteNumber(int *cellCoordinates, int sublattice, int *translation)
{
  cout << "Periodizing entry"<<Dimension-1<<endl;
  int Result=this->Periodize(cellCoordinates[Dimension-1], Dimension-1, translation[Dimension-1]);
  for (int i=Dimension-2; i>-1; --i)
    {
      Result*=this->PeriodicRep[i];
      Result+=this->Periodize(cellCoordinates[i], i, translation[i]);
    }
  Result*=NbrSitesPerCell;
  Result+=sublattice;
  return Result%NbrSites;
}


// request address of partners of site
// nbrSite = number of site whose partners to request
// nbrNeighbors = number of partners found
// Neighbors = array to partner sites
// phases = values of phase for tunnelling matrix element
void LatticePhases::GetNeighbors(int nbrSite, int &nbrNeighbors, int * &neighbors, double * &phases)
{
  if ((nbrSite>-1)&&(nbrSite<NbrSites))
    {
      neighbors = this->Neighbors[nbrSite];
      nbrNeighbors = this->NbrNeighbors[nbrSite];
      phases = this->TunnellingPhases[nbrSite];
    }
  else
    {
      nbrNeighbors = 0;
      neighbors = NULL;
      phases = NULL;
    }
}

// get total number of hopping terms
int LatticePhases::GetNbrHoppingTerms()
{
  int sum=0;
  for (int i=0; i<NbrSites; ++i)
    sum += this->NbrNeighbors[i];
  return sum;
}

// get a string describing the lattice geometry
// 
char *LatticePhases::GeometryString()
{
  char *rst = new char[strlen(this->Descriptor)+20];
  sprintf(rst,"%s_%d", this->Descriptor, this->PeriodicRep[0]);
  for (int i=1; i<Dimension; ++i)
    sprintf(rst,"%sx%d", rst, this->PeriodicRep[i]);
  return rst;
}


// add an option group containing all options related to the LatticeGeometry options
//
// manager = pointer to the option manager
void LatticePhases::AddOptionGroup(OptionManager* manager)
{
  LatticePhases::Options = manager;
  OptionGroup* LatticeGroup  = new OptionGroup ("lattice options");
  (*(LatticePhases::Options)) += LatticeGroup;

  (*LatticeGroup) += new SingleStringOption ('L', "lattice-definition", "File defining the geometry of the lattice");
  (*LatticeGroup) += new MultipleIntegerOption ('C', "cells-repeat", "number of times unit cell is repeated in the x-, y-,..., dim- directions of the lattice (overriding default given in definition)", ',');
  (*LatticeGroup) += new BooleanOption ('\n', "normalize-lattice", "normalize unit cell area to #of sites, and field strength to one");
}



OptionManager* LatticePhases::Options=NULL;

int MyRound(double a) {
return int(a + 0.5);
}

// simple sort algorithm
// array = integer array to be sorted
// length = length of array
void LatticePhases::ArraySort(int* array, int length)
{
  int inc = MyRound(length/2.0);
  int tmpI;
  while (inc > 0)
    {
      for (int i = inc; i< length; ++i)
	{
	  tmpI = array[i];
	  int j = i;
	  while ((j>=inc) && ( array[j-inc] < tmpI) )
	    {
	      array[j] = array[j - inc];
	      j = j - inc;
	    }
	  array[j] = tmpI;
	}
      inc = round(inc / 2.2);
    }
}

// calculate the tunnelling phase between two given sites from the gauge
// s1 = start site
// s2 = end site
// cellTranslation = indicating whether translation across a boundary ocurred
// return = relative phase
double LatticePhases::GetTunnellingPhaseFromGauge(int s1, int s2, int *cellTranslation)
{
  if (this->HaveGauge)
    {
      double Result=0.0;
      // calculate site coordinates
      int *S1Coordinates = new int[this->Dimension];
      int *S2Coordinates = new int[this->Dimension];
      int S1Sublattice, S2Sublattice;
      this->GetSiteCoordinates(s1, S1Coordinates, S1Sublattice);
      this->GetSiteCoordinates(s2, S2Coordinates, S2Sublattice);
      RealVector Position1(this->Dimension,true);
      RealVector Position2(this->Dimension,true);      
      RealVector Translation(this->Dimension,true);
      for (int i=0; i<Dimension; ++i)
	{
	  Position1.AddLinearCombination((double)S1Coordinates[i],LatticeVectors[i]);
	  Position2.AddLinearCombination((double)S2Coordinates[i],LatticeVectors[i]);
	  if (cellTranslation!=NULL)
	    Translation.AddLinearCombination((double)cellTranslation[i],LatticeVectors[i]);
	}      
      Position2.AddLinearCombination(-1.0,Translation);
      RealVector CellPosition2(this->Dimension);
      CellPosition2.Copy(Position2);
      Position1.AddLinearCombination(1.0,SubLatticeVectors[S1Sublattice]);
      Position2.AddLinearCombination(1.0,SubLatticeVectors[S2Sublattice]);
      CellPosition2.AddLinearCombination(1.0,SubLatticeVectors[S2Sublattice]);
      
      cout << "Position1="<<endl<<Position1;
      cout << "Position2="<<endl<<Position2;
      cout << "Translation="<<endl<<Translation;
      if (this->Dimension==2)
	{
	  Result += 0.5*GaugeAxx*(Position2[0]*Position2[0]-Position1[0]*Position1[0]);
	  Result += 0.5*GaugeAyx*(Position2[0]-Position1[0])*(Position2[1]+Position1[1]);
	  Result += 0.5*GaugeAxy*(Position2[1]-Position1[1])*(Position2[0]+Position1[0]);
	  Result += 0.5*GaugeAyy*(Position2[1]*Position2[1]-Position1[1]*Position1[1]);
	  // xxx: check signs and prefactors here!
	  cout << "Raw phase="<<Result;
	  if (Translation.SqrNorm()>1e-14)
	    {
	      if (Translation[0]>0.0)
		{
		  // translate in x-direction first:
		  double MagneticTranslation =(GaugeAxy*Translation[0]+GaugeAyy*Translation[1])*CellPosition2[1]; // (CellPosition2[1]+0.5*Translation[1]);
		  // then translate in y-direction
		  MagneticTranslation+=(GaugeAxx*Translation[0]+GaugeAyx*Translation[1])*(CellPosition2[0]+Translation[0]); // (CellPosition2[0]+Translation[0]+0.5*Translation[0]);
		  
		  Result += MagneticTranslation;
		  cout << ", magnetic translation: "<<MagneticTranslation;
		}
	      else
		{
		  // translate in y-direction first:
		  double MagneticTranslation =(GaugeAxx*Translation[0]+GaugeAyx*Translation[1])*CellPosition2[0]; // (CellPosition2[0]+0.5*Translation[0]);
		  // translate in x-direction first:
		  MagneticTranslation+= (GaugeAxy*Translation[0]+GaugeAyy*Translation[1])*(CellPosition2[1]+Translation[1]); // (CellPosition2[1]+Translation[1]+0.5*Translation[1]);
		  Result += MagneticTranslation;
		  cout << ", magnetic translation: "<<MagneticTranslation;
		}
	    }
	  cout << ", after corrections"<<Result<<endl;
	}
      else
	{
	  cout << "Need to define LatticePhases::GetTunnellingPhaseFromGauge for dimension d>2"<<endl;
	}
      delete [] S1Coordinates;
      delete [] S2Coordinates;
      return Result;
    }
  else
    return 0.0;
}
