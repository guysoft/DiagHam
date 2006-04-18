////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                       Copyright (C) 2006 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//            class for elementary factor in expansion of CF Orbitals         //
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

#ifndef DERIVATIVE_PRODUCT
#define DERIVATIVE_PRODUCT

#include "config.h"
#include "GeneralTools/List.h"
#include "MathTools/Complex.h"
#include "DerivativeProductFactor.h"
#include "JainCFOnSphereOrbitals.h"

class SumDerivativeProduct;

class DerivativeProduct
{
 private:  
  friend class JainCFOnSphereOrbitals;
  friend class SumDerivativeProduct;
  
 public:
  DerivativeProduct();
  DerivativeProduct( const DerivativeProductFactor &toWrap);
  DerivativeProduct( const DerivativeProduct &toCopy);
  DerivativeProduct(DerivativeProduct &Reference, List<DerivativeProductFactor> PriorFactors,
		    List<DerivativeProductFactor> LaterFactors);
  
  ~DerivativeProduct();

  void Simplify();

  bool isNonZero();

  Complex getValue();
  
  SumDerivativeProduct Derivative( int DeriveU, int DeriveV=0);
  
  DerivativeProduct& operator*= (const DerivativeProduct &toMultiply);
  DerivativeProduct& operator*= (DerivativeProductFactor &toMultiply);

  bool operator ^ ( DerivativeProduct &other);
  bool operator < ( DerivativeProduct &other);
  bool operator > ( DerivativeProduct &other);

  // Output Stream overload
  //
  // str = reference on output stream
  // D = DerivativeProductFactor
  // return value = referenceint GetNumSites(){return this->NSites;} on output stream
  friend ostream& operator << (ostream& str, DerivativeProduct& D);

  
 protected:
  JainCFOnSphereOrbitals *CFOrbitals;
  List<DerivativeProductFactor> ProductFactors;
  double PreFactor;
};


#endif //DERIVATIVE_PRODUCT
