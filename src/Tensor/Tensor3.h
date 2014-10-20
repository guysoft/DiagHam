
#ifndef _TENSOR3_H
#define _TENSOR3_H

#include<iostream>

using std::cout;
using std::endl;

template <class T> class Tensor3
{
 private :
  unsigned int FirstDimension;
  unsigned int SecondDimension;
  unsigned int ThirdDimension;
  T * TensorElements;
  
 public :
  
  Tensor3();
  Tensor3(int firstDimension,int secondDimension, int thirdDimension); 
  Tensor3(int firstDimension,int secondDimension, int thirdDimension,bool initiateFlag);
  // copy constructor (without duplicating datas)
  //
  Tensor3(const Tensor3 & tensor3); 
  ~Tensor3(){ delete [] this->TensorElements;}
  
  // assignement (without duplicating datas)
  //
  // M = matrix to copy
  // return value = reference on modified matrix
  Tensor3 & operator = (const Tensor3 & tensor3);
  
  void PrintTensor();
  
  T & operator()(unsigned int index1, unsigned int index2, unsigned int index3); 
};

template <class T>
Tensor3<T>::Tensor3()
{
  this->FirstDimension = 0;
  this->SecondDimension = 0;
  this->ThirdDimension = 0;
  this->TensorElements = 0;
}

template <class T>
Tensor3<T>::Tensor3(int firstDimension,int secondDimension, int thirdDimension)
{
  this->FirstDimension = firstDimension;
  this->SecondDimension = secondDimension;
  this->ThirdDimension = thirdDimension;
  this->TensorElements = new T [this->FirstDimension* this->SecondDimension* this->ThirdDimension];
}

template <class T>
Tensor3<T>::Tensor3(int firstDimension,int secondDimension, int thirdDimension,bool initiateFlag)
{
  this->FirstDimension = firstDimension;
  this->SecondDimension = secondDimension;
  this->ThirdDimension = thirdDimension;
  this->TensorElements = new T [this->FirstDimension* this->SecondDimension* this->ThirdDimension];
  if (initiateFlag)
    {
      for (int i = 0; i< this->FirstDimension* this->SecondDimension* this->ThirdDimension ; i++)
	this->TensorElements[i] = 0;
    }
}


template <class T>
Tensor3<T>::Tensor3(const Tensor3<T> & tensor3)
{
  this->FirstDimension = tensor3.FirstDimension;
  this->SecondDimension = tensor3.SecondDimension;
  this->ThirdDimension = tensor3.ThirdDimension;
  this->TensorElements = tensor3.TensorElements;
}


// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

template <class T>
Tensor3<T> & Tensor3<T>::operator = (const Tensor3<T> & tensor3) 
{
  this->FirstDimension = tensor3.FirstDimension;
  this->SecondDimension = tensor3.SecondDimension;
  this->ThirdDimension = tensor3.ThirdDimension;
  this->TensorElements = tensor3.TensorElements;
  return *this;
}
  

template <class T>
inline T & Tensor3<T>::operator() (unsigned int index1, unsigned int index2, unsigned int index3) 
{ 
  return TensorElements[index1 + index2*this->FirstDimension+this->FirstDimension*this->SecondDimension*index3];
} 

template <class T>
void Tensor3<T>::PrintTensor()
{
  for(int i=0; i< this->FirstDimension; i++)
    {
      for(int k=0; k< this->SecondDimension; k++)
	{
	  for(int p=0; p< this->ThirdDimension; p++)
	    {
	      cout <<i << " "<<k<< " "<<p<<" " << (*this)(i,k,p) << endl;
	    }
	}
    }
}

#endif
