////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2009 Nicolas Regnault                  //
//                         class author: Gunnar Möller                        //
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

#include "LatticeConnections.h"
#include "Options/Options.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>
using std::cout;
using std::endl;



// generate the object using options from Option Manager
//
LatticeConnections::LatticeConnections()
{
  if (LatticeConnections::Options==NULL)
    {
      cout << "Define the OptionManager, first, before creating any LatticeConnections"<<endl;
      exit(1);
    }

  ConfigurationParser LatticeDefinition;
  if (LatticeDefinition.Parse(this->Options->GetString("lattice-definition")) == false)
    {
      LatticeDefinition.DumpErrors(cout) << endl;
      exit(-1);
    }
  if ((Descriptor = LatticeDefinition["Descriptor"] ) == NULL)
    {
      cout << "Attention, 'Descriptor' is not defined, unnamed lattice geometry!" << endl;
      Descriptor = new char[10];
      sprintf(Descriptor,"unnamed");
    }
  if ((LatticeDefinition.GetAsSingleInteger("NbrSites", NbrSitesPerCell) == false) || (NbrSitesPerCell <= 0))
    {
      cout << "NbrSites is not defined or as a wrong value" << endl;
      exit(-1);
    }
  if ((LatticeDefinition.GetAsSingleInteger("Dimension", Dimension) == false) || (Dimension <= 0))
    {
      cout << "Dimension is not defined or as a wrong value" << endl;
      exit(-1);
    }
  int TmpDimension;
  PeriodicRep = LatticeConnections::Options->GetIntegers("cells-repeat",TmpDimension);
  if (TmpDimension==0)
    if ((LatticeDefinition.GetAsIntegerArray("PeriodicRepeat",',',PeriodicRep, TmpDimension) == false))
      {
	cout << "PeriodicRepeat is not defined or as a wrong value, simulating a single unit cell" << endl;
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
  SubLatticeVectors.Resize(NbrSitesPerCell, Dimension);
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
    }
  // determine connectivity:
  char ***NeighborString;
  int NbrPairs;
  int *NbrValues;
  RealSymmetricMatrix NeighborsInCellMatrix(NbrSitesPerCell, true);
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
	  if (NbrValues[p]!=2)
	    {
	      cout << "error while decoding NeighborsInCell in " << this->Options->GetString("lattice-definition") << endl;
	      cout << "Indicate paires of neighboring sites separated by commas and different pairs by bars: "
		   << "NeighborsInCell = s1,s2 | s3, s4 | ..."<<endl;
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
	      else NeighborsInCellMatrix(s1,s2)=1.0;
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
  cout << "NeighborsInCell="<<NeighborsInCellMatrix;
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
      for (int p=0; p<NbrPairs; ++p)
	{
	  int s1, s2;
	  if (NbrValues[p]!=2)
	    {
	      cout << "error while decoding "<<FieldName<<" in " << this->Options->GetString("lattice-definition") << endl;
	      cout << "Indicate paires of neighboring sites separated by commas and different pairs by bars: "
		   << FieldName <<" = s1,s2 | s3, s4 | ..."<<endl;
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
	      else NeighborsAcrossBoundary[d]->SetMatrixElement(s1, s2, 1.0);
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
      cout << FieldName<<"="<<*NeighborsAcrossBoundary[d]<<endl;
    }  
  
  this->Neighbors = new int*[NbrSites];
  this->NbrNeighbors = new int[NbrSites];
  for (int i=0; i<NbrSites; ++i)
    this->NbrNeighbors[i] = 0;
  int *TmpNeighbors = new int[NbrNeighborCells*NbrSites];
  
  this->NbrCells = this->PeriodicRep[0];
  for (int d=1; d<Dimension; ++d)
    this->NbrCells *= this->PeriodicRep[d];

  int *CellCoordinates = new int[Dimension];
  int *CellCoordinates2 = new int[Dimension];
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
	      cout << "Site 2="<<Site2;
	      if (NeighborsInCellMatrix(i,j)>0.0)
		{
		  cout << "... is neighbor"<<endl;
		  TmpNeighbors[NbrNeighbors[Site1]]=Site2;
		  ++NbrNeighbors[Site1];
		}
	      else cout << "... not neighbor"<<endl;
	      for (int d=0; d<NbrNeighborCells; ++d)
		{
		  if ((*(NeighborsAcrossBoundary[d]))(i,j)!=0.0)
		    {
		      for (int k=0; k<Dimension; ++k)
			CellCoordinates2[k]=CellCoordinates[k]+NeighborCells[d][k];
		      Site3 = this->GetSiteNumber(CellCoordinates2, j);
		      cout << "additional neighbor from NeigborCell "<<d<<" at "<<
			CellCoordinates2[0]<<", "<<CellCoordinates2[1]<<", "<<j<<" : Site 3="<<Site3<<endl;
		      TmpNeighbors[NbrNeighbors[Site1]]=Site3;
		      ++NbrNeighbors[Site1];
		    }
		}	      
	    }
	  if (NbrNeighbors[Site1]>0)
	    Neighbors[Site1] = new int[NbrNeighbors[Site1]];
	  else Neighbors[Site1] = NULL;
	  for (int k=0; k<NbrNeighbors[Site1]; ++k)
	    Neighbors[Site1][k] = TmpNeighbors[k];
	}
    }

  this->Partners = new int*[NbrSites];
  this->NbrPartners = new int[NbrSites];
  for (int i=0; i<NbrSites; ++i)
    this->NbrPartners[i]=0;
  
  for (int i=0; i<NbrSites; ++i)
    {
      cout <<  "Neighbors["<<i<<"] = "<<Neighbors[i][0];
      for (int k=1; k<NbrNeighbors[i]; ++k) cout <<" "<<Neighbors[i][k];
      cout << endl;
      this->ArraySort(Neighbors[i], NbrNeighbors[i]);
      cout <<  "sorted = "<<Neighbors[i][0];
      for (int k=1; k<NbrNeighbors[i]; ++k) cout <<" "<<Neighbors[i][k];
      cout << endl;
      int j=0; 
      while((j<NbrNeighbors[i])&&(Neighbors[i][j]>=i)) ++j;
      this->NbrPartners[i] = j;
      if (j>0)
	{
	  this->Partners[i] = new int[this->NbrPartners[i]];
	  for (int k=0; k<NbrPartners[i]; ++k)
	    this->Partners[i][k]=this->Neighbors[i][k];
	}
      else this->Partners[i] = NULL;
      if (NbrPartners[i]>0)
	{
	  cout <<  "Partners["<<i<<"] = "<<Partners[i][0];
	  for (int k=1; k<NbrPartners[i]; ++k) cout <<" "<<Partners[i][k];
	  cout << endl;
	}
      else cout << "no partners"<<endl;
    }

  cout << "LatticeConnections created"<<endl;

  for (int d=0; d<NbrNeighborCells; ++d)
    delete NeighborsAcrossBoundary[d];
  delete [] NeighborsAcrossBoundary;
  delete [] CellCoordinates;
  delete [] CellCoordinates2;
  delete [] FieldName;
  delete [] TmpNeighbors;
}

// destructor
//
LatticeConnections::~LatticeConnections()
{
  if (NbrSites!=0)
    {
      delete [] this->PeriodicRep;
      for (int i=0; i<NbrSites; ++i)
	{
	  if (this->Neighbors[i]!=NULL)
	    delete [] this->Neighbors[i];
	  if (this->Partners[i]!=NULL)
	    delete [] this->Partners[i];
	}      
      delete [] this->Neighbors;
      delete [] this->Partners;
      delete [] this->NbrNeighbors;
      delete [] this->NbrPartners;
      for (int i=0; i<NbrNeighborCells; ++i)
	delete [] this->NeighborCells[i];
      delete [] this->NeighborCells;
      delete [] this->Descriptor;
    }
}


// get cell coordinates given the number of the unit cell
// nbrCell = cell to be looked up
// cellCoordinates = resulting coordinates, has to be reserved prior to call
void LatticeConnections::GetCellCoordinates(int nbrCell, int *cellCoordinates)
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
void LatticeConnections::GetSiteCoordinates(int nbrSite, int *cellCoordinates, int &sublattice)
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
int LatticeConnections::GetSiteNumber(int *cellCoordinates, int sublattice)
{
  int Result=this->Periodize(cellCoordinates[Dimension-1],Dimension-1);
  for (int i=Dimension-1; i>-1; --i)
    {
      Result*=this->PeriodicRep[i];
      Result+=this->Periodize(cellCoordinates[i],i);
    }
  Result*=NbrSitesPerCell;
  Result+=sublattice;
  return Result%NbrSites;
}


// request address of partners of site
// nbrSite = number of site whose partners to request
// partners = array to partner sites
// nbrPartners = number of partners found
void LatticeConnections::GetPartners(int nbrSite, int * &partners, int &nbrPartners)
{
  if ((nbrSite>0)&&(nbrSite<NbrSites))
    {
      partners = this->Partners[nbrSite];
      nbrPartners = this->NbrPartners[nbrSite];
    }
  else
    {
      partners = NULL;
      nbrPartners = 0;
    }
}

// get a string describing the lattice geometry
// 
char *LatticeConnections::GeometryString()
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
void LatticeConnections::AddOptionGroup(OptionManager* manager)
{
  LatticeConnections::Options = manager;
  OptionGroup* LatticeGroup  = new OptionGroup ("lattice options");
  (*(LatticeConnections::Options)) += LatticeGroup;

  (*LatticeGroup) += new SingleStringOption  ('L', "lattice-definition", "File defining the geometry of the lattice");
  (*LatticeGroup) += new MultipleIntegerOption  ('c', "cells-repeat", "number of times unit cell is repeated in the x-, y-,..., dim- directions of the lattice (overriding default given in definition)", ',');
}



OptionManager* LatticeConnections::Options=NULL;

int round(double a) {
return int(a + 0.5);
}

// simple sort algorithm
// array = integer array to be sorted
// length = length of array
void LatticeConnections::ArraySort(int* array, int length)
{
  int inc = round(length/2.0);
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
