#include "MPOPeratorDefinedByFiles.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

using std::cout;
using std::endl;
 


MPOPeratorDefinedByFiles::MPOPeratorDefinedByFiles()
{
}



MPOPeratorDefinedByFiles::MPOPeratorDefinedByFiles(int nbrSites, MultiColumnASCIIFile tensorElementsFile, MultiColumnASCIIFile boundaryVectorsFile)
{
  this->InitializeTensorsElements(tensorElementsFile);
  this->NbrSites = nbrSites;
  this->InitializeBoundaryVectors(boundaryVectorsFile);
}



MPOPeratorDefinedByFiles::~MPOPeratorDefinedByFiles()
{
  delete [] LeftVector;
  delete [] RightVector;
  delete [] ElementsValues;
  delete [] IndexValues;    
}


void MPOPeratorDefinedByFiles::InitializeTensorsElements(MultiColumnASCIIFile tensorElementsFile)
{

 this->NbrNonZeroElements = tensorElementsFile.GetNbrLines();

 int* IndexDown  = tensorElementsFile.GetLineAsIntegerArray (0);
 int* IndexUp = tensorElementsFile.GetLineAsIntegerArray (1);
 int* IndexLeft = tensorElementsFile.GetLineAsIntegerArray (2);
 int* IndexRight = tensorElementsFile.GetLineAsIntegerArray (3);
 this->ElementsValues = tensorElementsFile.GetAsDoubleArray (4);
 int TmpPhysicalDimension = 0;
 int TmpMPODimension = 0;


  cout <<"Nbr Non zero Elements = "<<  this->NbrNonZeroElements<<endl;
  this->ElementsValues = new double [this->NbrNonZeroElements];
  this->IndexValues = new unsigned int[this->NbrNonZeroElements];

  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
 {
  if (IndexDown[i] > TmpPhysicalDimension)
     TmpPhysicalDimension = IndexDown[i];
  if (IndexLeft[i] > TmpMPODimension)
     TmpMPODimension = IndexLeft[i];
   this->IndexValues[i] = this->GetTensorIndexFromAllIndices(IndexDown[i],IndexUp[i], IndexLeft[i] ,IndexRight[i]);
}

  this->PhysicalDimension = TmpPhysicalDimension+1;
  this->MPOBondDimension =  TmpMPODimension+1;
  delete [] IndexDown;
  delete [] IndexUp;
  delete [] IndexLeft;
  delete [] IndexRight;
}


void MPOPeratorDefinedByFiles::InitializeBoundaryVectors(MultiColumnASCIIFile boundaryVectorsFile)
{
  if ( boundaryVectorsFile.GetNbrLines() !=   this->MPOBondDimension)
  {
    cout <<"error size of the boundary vectors are not compatible with the one" <<endl;
  }
  else
  {
    if ( boundaryVectorsFile.GetNbrLines() == 2  )
{
     this->LeftVector = boundaryVectorsFile.GetLineAsDoubleArray (0); 
     this->RightVector = boundaryVectorsFile.GetLineAsDoubleArray (1);
}
else
{
   if ( boundaryVectorsFile.GetNbrLines() == 1  )
{
     this->LeftVector = boundaryVectorsFile.GetLineAsDoubleArray (0); 
     this->RightVector = boundaryVectorsFile.GetLineAsDoubleArray (0);
}
}
  }

}


void  MPOPeratorDefinedByFiles::PrintTensorElements()
{
  cout <<"#Tensor index indexDown indexUp indexLeft indexRight Check Index Values" <<endl;
  unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);      
      int Tmp = GetTensorIndexFromAllIndices( MPOIndiceDown,  MPOIndiceUp,  MPOIndiceLeft,  MPOIndiceRight);
      cout << this->IndexValues[i] <<" "<<MPOIndiceDown<< " "<< MPOIndiceUp<< " "<< MPOIndiceLeft<< " "<<MPOIndiceRight<<" " <<Tmp<<" "<< this->ElementsValues[i]<<endl;
    }
}
