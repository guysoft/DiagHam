////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2006 Gunnar Moeller                       //
//                                                                            //
//                                                                            //
//              class for elementary factor in expansion of CF Orbitals       //
//                                                                            //
//                        last modification : 17/04/2006                      //
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

#include "DerivativeProduct.h"
#include "DerivativeProductFactor.h"
#include "SumDerivativeProduct.h"
#include "GeneralTools/ListIterator.h"
#include "GeneralTools/List.h"

#include <iostream>
using std::cout;
using std::endl;

DerivativeProduct::DerivativeProduct()
{
  this->CFOrbitals=NULL;
  this->PreFactor=1.0;
}

DerivativeProduct::DerivativeProduct( const DerivativeProductFactor &toWrap)
{
  this->CFOrbitals=toWrap.CFOrbitals;
  this->PreFactor=1.0;
  this->ProductFactors+=toWrap;
}

DerivativeProduct::DerivativeProduct(const DerivativeProduct &toCopy)
{
  this->CFOrbitals=toCopy.CFOrbitals;
  this->PreFactor=toCopy.PreFactor;
  this->ProductFactors=toCopy.ProductFactors;
}

DerivativeProduct::DerivativeProduct(DerivativeProduct &Reference, List<DerivativeProductFactor> PriorFactors,
		    List<DerivativeProductFactor> LaterFactors)
{
  this->CFOrbitals=Reference.CFOrbitals;
  this->PreFactor=Reference.PreFactor;
  this->ProductFactors=PriorFactors;
  this->ProductFactors.Link(LaterFactors);
  //cout << "Nbr Elements after Link: " << this->ProductFactors.GetNbrElement();
}

DerivativeProduct::~DerivativeProduct()
{
}

bool DerivativeProduct::isNonZero()
{
  DerivativeProductFactor *Factor;
  for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
    if(Factor->isZero()) return false;
  return true;
}

Complex DerivativeProduct::getValue()
{
  return Complex();
}

void DerivativeProduct::Simplify()
{
  DerivativeProductFactor *Factor;
  for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
    {
      this->PreFactor*=Factor->PreFactor;
      Factor->PreFactor=1.0;
    }
}


SumDerivativeProduct DerivativeProduct::Derivative(int DeriveU, int DeriveV)
{
  SumDerivativeProduct result(this->CFOrbitals);
  List<DerivativeProductFactor> PriorFactors;
  List<DerivativeProductFactor> LaterFactors;
  int elements=this->ProductFactors.GetNbrElement();
  //cout <<endl<< "Derivative of " << *this << endl;
  for (int pos=0; pos < elements; ++pos)
    {
      DerivativeProduct derivative = this->ProductFactors[pos].Derivative(DeriveU,DeriveV);
      //cout << "Element " << pos << "= " << this->ProductFactors[pos]<< endl;
      //cout << "its der. " << derivative << endl;      
      if (derivative.isNonZero())
	{
	  PriorFactors.DeleteList();
	  if (pos>0) PriorFactors=Extract(this->ProductFactors, 0, pos);
	  LaterFactors.DeleteList();
	  if (pos+1 < elements) LaterFactors=Extract(this->ProductFactors, pos+1);      
	  DerivativeProduct tmp(*this,PriorFactors,LaterFactors);
	  //cout << "other terms: " << tmp << endl;
	  tmp*=derivative;
	  //cout << "Full term: " << tmp << endl;
	  result +=tmp;
	}
    }
  return result;
}


/*DerivativeProduct& DerivativeProduct::operator*= (DerivativeProductFactor &toMultiply)
{
  DerivativeProductFactor *Factor;
  bool toBeInserted=true;
  int pos=0;
  if (toMultiply.isScalar()) // simple number to be multiplied...
    {
      this->PreFactor*=toMultiply.PreFactor;
      return *this;
    }
  for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL; ++pos)
    {
      if ( (*Factor) < toMultiply )
	{
	  if (Factor->isScalar())
	    Factor->Multiply(toMultiply);
	  else
	    this->ProductFactors.Insert(toMultiply,pos);
	  toBeInserted=false;
	  break;
	}
      else if ( Factor->Multiply(toMultiply) )
	{
	  toBeInserted=false;
	  break;
	}
    }
  if (toBeInserted) this->ProductFactors+=toMultiply;
  return *this;
}
*/

DerivativeProduct& DerivativeProduct::operator*= (DerivativeProductFactor &toMultiply)
{
  DerivativeProductFactor *Factor;
  bool toBeInserted=true;
  if (this->ProductFactors.GetNbrElement()>0)
    {
      for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
	{
	  if ( Factor->Multiply(toMultiply) )
	    {
	      toBeInserted=false;
	      break;
	    }
	}
    }
  if (toBeInserted) this->ProductFactors+=toMultiply;
  //cout << "before simplify: " << *this << endl;
  this->Simplify();
  //cout << "after simplify:  " << *this << endl;;
  return *this;
}

DerivativeProduct& DerivativeProduct::operator*= (const DerivativeProduct &toMultiply)
{
  DerivativeProductFactor *Factor;
  if (this->ProductFactors.GetNbrElement()>0)
    for (ListIterator<DerivativeProductFactor> LI(toMultiply.ProductFactors); (Factor=LI())!=NULL; )
      {
	(*this)*=(*Factor);
      }
  else *this = toMultiply;
  return *this;
}

/*
bool DerivativeProduct::operator ^ (DerivativeProduct &other)
{
  if (this->ProductFactors.GetNbrElement() != other.ProductFactors.GetNbrElement())
    return false;
  else
    {
      //cout <<"this  "<< this->ProductFactors;
      //cout << "other "<< other.ProductFactors;
      this->ProductFactors.UpOrder();
      other.ProductFactors.UpOrder();
      //cout <<"this  "<< this->ProductFactors;
      //cout << "other "<< other.ProductFactors;
      ListIterator<DerivativeProductFactor> OtherLI(other.ProductFactors);
      DerivativeProductFactor *Factor;
      DerivativeProductFactor *OtherFactor;
      for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
	{
	  OtherFactor=OtherLI();
	  if ( (*Factor) !=  (*OtherFactor) ) return false;
	}
      return true;
    }
}
*/

bool DerivativeProduct::operator ^ (DerivativeProduct &other)
{
  if (this->ProductFactors.GetNbrElement() != other.ProductFactors.GetNbrElement())
    return false;
  else
    {
      //cout <<"this  "<< this->ProductFactors;
      //cout << "other "<< other.ProductFactors;
      this->ProductFactors.UpOrder();
      other.ProductFactors.UpOrder();
      //cout <<"this  "<< this->ProductFactors;
      //cout << "other "<< other.ProductFactors;
      ListIterator<DerivativeProductFactor> OtherLI(other.ProductFactors);
      DerivativeProductFactor *Factor;
      DerivativeProductFactor *OtherFactor;
      for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
	{
	  OtherFactor=OtherLI();
	  if ( (*Factor) !=  (*OtherFactor) ) return false;
	}
      return true;
    }
}

bool DerivativeProduct::operator < (DerivativeProduct &other)
{
  if (this->ProductFactors.GetNbrElement() < other.ProductFactors.GetNbrElement())
    return true;
  else
    {     
      ListIterator<DerivativeProductFactor> OtherLI(other.ProductFactors);
      DerivativeProductFactor *Factor;
      DerivativeProductFactor *OtherFactor;
      for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
	{
	  OtherFactor=OtherLI();
	  //cout << "Comparing factors " << *Factor << " and " << *OtherFactor << endl;
	  if ( (*Factor) <  (*OtherFactor) ) return true;
	}
      return false;
    }
}

bool DerivativeProduct::operator > (DerivativeProduct &other)
{
  if (this->ProductFactors.GetNbrElement() > other.ProductFactors.GetNbrElement())
    return true;
  else
    {     
      ListIterator<DerivativeProductFactor> OtherLI(other.ProductFactors);
      DerivativeProductFactor *Factor;
      DerivativeProductFactor *OtherFactor;
      for (ListIterator<DerivativeProductFactor> LI(this->ProductFactors); (Factor=LI())!=NULL;)
	{
	  OtherFactor=OtherLI();
	  //cout << "Comparing factors " << *Factor << " and " << *OtherFactor << endl;
	  if ( (*Factor) >  (*OtherFactor) ) return true;
	}
      return false;
    }
}

ostream& operator << (ostream& str, DerivativeProduct& D)
{
  DerivativeProductFactor *Factor;
  if (D.ProductFactors.GetNbrElement()>=1)
    {
      if (D.PreFactor!=1.0) str << D.PreFactor << "*" ;
      ListIterator<DerivativeProductFactor> LI(D.ProductFactors);
      Factor=LI();
      str << (*Factor);
      for (; (Factor=LI())!=NULL;)
	str <<"*"<< (*Factor);    
    }
  else str << "1.0";
  str << " ";
  return str;
}
