////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                     Copyright (C) 2006 Gunnar Moeller                      //
//                                                                            //
//                                                                            //
//           class for elementary factor in expansion of CF Orbitals          //
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

#include "SumDerivativeProduct.h"
#include "GeneralTools/ListIterator.h"

#include <iostream>
using std::cout;
using std::endl;


SumDerivativeProduct::SumDerivativeProduct()
{
  this->CFOrbitals=NULL;
}

SumDerivativeProduct::SumDerivativeProduct(JainCFOnSphereOrbitals *CFOrbitals)
{
  this->CFOrbitals=CFOrbitals;
}

SumDerivativeProduct::SumDerivativeProduct(const DerivativeProductFactor &toWrap)
{
  this->CFOrbitals=toWrap.CFOrbitals;
  this->Summands+=DerivativeProduct(toWrap);
}

SumDerivativeProduct::SumDerivativeProduct(const DerivativeProduct &toWrap)
{
  this->CFOrbitals=toWrap.CFOrbitals;
  this->Summands+=toWrap;
}

SumDerivativeProduct::SumDerivativeProduct(const SumDerivativeProduct &toCopy)
{
  this->CFOrbitals=toCopy.CFOrbitals;
  this->Summands=toCopy.Summands;
}

SumDerivativeProduct::~SumDerivativeProduct()
{
}

Complex SumDerivativeProduct::getValue(int particle)
{
  Complex Sum(0.0);
  DerivativeProduct *Product;
  for (ListIterator<DerivativeProduct> LI(this->Summands); (Product=LI())!=NULL; )
    Sum+=Product->getValue(particle);
  return Sum;
}

Complex* SumDerivativeProduct::getValues()
{
  Complex *result, *CPtr, *TmpSummand;
  DerivativeProduct *Product;
  ListIterator<DerivativeProduct> LI(this->Summands);
  Product=LI();
  result = Product->getValues();
  while((Product=LI())!=NULL)
    {
      TmpSummand=Product->getValues();
      CPtr=result;
      for (int i=0; i<CFOrbitals->GetNbrParticles(); ++i)
	*(CPtr++) *= *(TmpSummand++);
    }
  return result;
}


void SumDerivativeProduct::TestHighestPowers()
{
  DerivativeProduct *Product;
  for (ListIterator<DerivativeProduct> LI(this->Summands); (Product=LI())!=NULL; )
    Product->TestHighestPowers();
}

SumDerivativeProduct SumDerivativeProduct::Derivative(int DeriveU, int DeriveV)
{
  SumDerivativeProduct result(this->CFOrbitals);
  DerivativeProduct *Product;
  //  cout << "Derivative of SumDerivativeProduct " << *this << endl;
  for (ListIterator<DerivativeProduct> LI(this->Summands); (Product=LI())!=NULL; )
    result += Product->Derivative(DeriveU, DeriveV);
  return result;
}

SumDerivativeProduct& SumDerivativeProduct::operator = (const SumDerivativeProduct& Assign)
{
  this->CFOrbitals=Assign.CFOrbitals;
  this->Summands=Assign.Summands;
  //this->Flag=Assign.Flag;
  //this->TmpSum=Assign.TmpSum;
  return *this;
}

SumDerivativeProduct& SumDerivativeProduct::operator*=(const DerivativeProduct &toMultiply)
{
  DerivativeProduct *Product;
  for (ListIterator<DerivativeProduct> LI(this->Summands); (Product=LI())!=NULL; )
    (*Product)*=toMultiply;
  return *this;
}

SumDerivativeProduct& SumDerivativeProduct::operator*= (const DerivativeProductFactor &toMultiply)
{
  DerivativeProduct *Product;
  for (ListIterator<DerivativeProduct> LI(this->Summands); (Product=LI())!=NULL; )
    (*Product)*=toMultiply;
  return *this;
}

SumDerivativeProduct& SumDerivativeProduct::operator+= ( DerivativeProduct &toAdd)
{
  DerivativeProduct *Summand;
  bool toBeInserted=true;
  int pos=0;
  //cout << "Adding " << toAdd << " to Sum " << *this << endl;
  for (ListIterator<DerivativeProduct> LI(this->Summands); (Summand=LI())!=NULL; ++pos)
    {
      if ( *Summand ^ toAdd ) // may be added together
	{
	  toBeInserted=false;
	  toAdd.Simplify();
	  Summand->Simplify();
	  Summand->PreFactor+=toAdd.PreFactor;
	}
    }
  if (toBeInserted) this->Summands+=toAdd;
  //cout << "Result: " << *this<< endl;
  return *this;
}

SumDerivativeProduct& SumDerivativeProduct::operator+= (const SumDerivativeProduct &toAdd)
{
  DerivativeProduct *Summand;
  for (ListIterator<DerivativeProduct> LI(toAdd.Summands); (Summand=LI())!=NULL;)
    {
      (*this)+=(*Summand);
    }
  this->Summands.UpOrder();  // clean up things a little...
  return *this;
}

SumDerivativeProduct& SumDerivativeProduct::operator+= (const DerivativeProductFactor &toAdd)
{
  (*this)+=DerivativeProduct(toAdd);
  return *this;
}

ostream& operator << (ostream& str, SumDerivativeProduct& S)
{
  DerivativeProduct *Summand;
  if (S.Summands.GetNbrElement()>=1)
    {
      ListIterator<DerivativeProduct> LI(S.Summands);
      Summand=LI();
      str << (*Summand);
      for (; (Summand=LI())!=NULL;)
	str <<"+" << (*Summand);    
    }
  else str << "0.0";
  str << " ";
  return str;
}
