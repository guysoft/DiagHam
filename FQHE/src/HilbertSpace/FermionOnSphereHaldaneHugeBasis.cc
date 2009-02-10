////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class of fermions on sphere using the Haldane basis            //
//                                                                            //
//                        last modification : 06/07/2006                      //
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


#include "config.h"
#include "HilbertSpace/FermionOnSphereHaldaneHugeBasis.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/ListIterator.h"
#include "GeneralTools/Endian.h"


#include <math.h>
#include <fstream>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// basic constructor
// 
// nbrFermions = number of fermions
// totalLz = twice the momentum total value
// lzMax = twice the maximum Lz value reached by a fermion
// maxFileSize = maximum file size (in MBytes)
// memory = amount of memory granted for precalculations
// referenceState = array that describes the reference state to start from
// symmetricFlag = indicate if a symmetric basis has to be used (only available if totalLz = 0)
// fullDimension = provide the full (i.e. without squeezing) Hilbert space dimension (0 if it has to be computed)

FermionOnSphereHaldaneHugeBasis::FermionOnSphereHaldaneHugeBasis (int nbrFermions, int totalLz, int lzMax, unsigned long maxFileSize, int* referenceState, unsigned long memory,
								  bool symmetricFlag, long fullDimension)
{
  this->TargetSpace = this;
  this->NbrFermions = nbrFermions;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->TotalLz = totalLz;
  this->LzMax = lzMax;
  this->NbrLzValue = this->LzMax + 1;
  this->SizeLimit = maxFileSize << 17;
  this->HilbertSpaceDimension = 0;
  this->ReferenceState = 0x0l;
  this->SymmetricFlag = symmetricFlag;
  this->MaximumSignLookUp = 16;
  int ReferenceStateLzMax = 0;
  for (int i = 0; i <= this->LzMax; ++i)
    {
      this->ReferenceState |= ((unsigned long) (referenceState[i] & 1)) << i;
      if (referenceState[i] == 1)
	ReferenceStateLzMax = i;
    }
  this->Flag.Initialize();

  if (fullDimension == 0l)
    this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrFermions, this->LzMax, this->TotalLz);
  else 
    this->LargeHilbertSpaceDimension = fullDimension;
  cout << "total dimension = " << this->LargeHilbertSpaceDimension << endl;
  ListIterator<unsigned long> IterFileSizes(this->FileSizes);
  unsigned long* TmpFileSize;
  int MinCommonLz = this->LzMax;
  while ((TmpFileSize = IterFileSizes()))
    {
      if (MinCommonLz > ((int) ((*TmpFileSize) & 0xffl)))
	MinCommonLz = ((*TmpFileSize) & 0xffl);
    }  
  int MaxPartialNbrFermions = this->LzMax - MinCommonLz;  
  this->StateHighestLzShift = MinCommonLz + 1;
  this->LargeHilbertSpaceDimension = 0l;
  int ShiftedTotalLz = (this->TotalLz + (this->NbrFermions * this->LzMax)) >> 1;
  this->PartialHilbertSpaceDimension = 0l;
  this->NbrFiles = 0;
  unsigned long LargestFile = 0;
  unsigned long DiskSpace = 0;
  for (int i = MaxPartialNbrFermions; i >= 0; --i)
    {
      int TmpMin = ((i * (i - 1)) >> 1);
      int TmpMax = ((MaxPartialNbrFermions - 1) * i) - ((i * (i - 1)) >> 1);
      int ShiftedTotalLz2 = ShiftedTotalLz - (i * this->StateHighestLzShift);
      int ShiftedTotalLz3  = ShiftedTotalLz2 - (((this->NbrFermions - i) * ((2 * MinCommonLz) - (this->NbrFermions - i) + 1)) >> 1);
      int ShiftedTotalLz4 = ShiftedTotalLz2 - (((this->NbrFermions - i) * (this->NbrFermions - i - 1)) >> 1);
      if (TmpMin < ShiftedTotalLz3)
	TmpMin = ShiftedTotalLz3;
      if (TmpMax > ShiftedTotalLz4)	
	TmpMax = ShiftedTotalLz4;
      if (TmpMax >= TmpMin)
	this->NbrFiles += (TmpMax - TmpMin) + 1;
    }
  this->StateDescriptionFileSizes = new unsigned long[this->NbrFiles];
  this->NbrFiles = 0;
  for (int i = MaxPartialNbrFermions; i >= 0; --i)
    {
      int TmpMin = ((i * (i - 1)) >> 1);
      int TmpMax = ((MaxPartialNbrFermions - 1) * i) - ((i * (i - 1)) >> 1);
      int ShiftedTotalLz2 = ShiftedTotalLz - (i * this->StateHighestLzShift);
      int ShiftedTotalLz3  = ShiftedTotalLz2 - (((this->NbrFermions - i) * ((2 * MinCommonLz) - (this->NbrFermions - i) + 1)) >> 1);
      int ShiftedTotalLz4 = ShiftedTotalLz2 - (((this->NbrFermions - i) * (this->NbrFermions - i - 1)) >> 1);
      if (TmpMin < ShiftedTotalLz3)
	TmpMin = ShiftedTotalLz3;
      if (TmpMax > ShiftedTotalLz4)	
	TmpMax = ShiftedTotalLz4;
      for (; TmpMax >= TmpMin; --TmpMax)
	{
          unsigned long TmpPartialDimension1 = this->ShiftedEvaluateHilbertSpaceDimension2(i, MaxPartialNbrFermions - 1, TmpMax);	  
	  unsigned long TmpPartialDimension2 = this->ShiftedEvaluateHilbertSpaceDimension2(this->NbrFermions - i, MinCommonLz, ShiftedTotalLz2 - TmpMax);
	  if (TmpPartialDimension1 == 0l)
	    TmpPartialDimension1 = 1l;
	  if (LargestFile < TmpPartialDimension2)
	    LargestFile = TmpPartialDimension2;
	  this->LargeHilbertSpaceDimension += TmpPartialDimension1 * TmpPartialDimension2;
	  this->PartialHilbertSpaceDimension += TmpPartialDimension1;
	  DiskSpace += TmpPartialDimension2;
	  this->StateDescriptionFileSizes[this->NbrFiles] = TmpPartialDimension2;
	  ++this->NbrFiles;
	}
    }
  cout << "total dimension = " << this->LargeHilbertSpaceDimension << "  (" << this->PartialHilbertSpaceDimension << ")" << endl;
  cout << "total requested disk space = " << (DiskSpace >> 17) <<  " Mbytes splitted in " << this->NbrFiles << " files " << endl;
  cout << "largest file = " << (LargestFile >> 7) << " KBytes (" << LargestFile<< " elements)" << endl;
  if (memory > (DiskSpace << 3))
    {
      cout << "disk storage is useless, keep everything in memory" << endl;
      this->NoDiskFlag = true;
    }
  else
    this->NoDiskFlag = false;
  this->StateDescription = new unsigned long [this->PartialHilbertSpaceDimension];
  this->StateDescriptionFileIndex = new int [this->PartialHilbertSpaceDimension];
  this->StateDescriptionIndexShift = new long [this->PartialHilbertSpaceDimension];
  this->PartialHilbertSpaceDimension = 0l;
  this->NbrFiles = 0;
  this->StateDescriptionIndexShift[0] = 0l;
  this->StateHighestLzToIndex = new int [1 << ((this->LzMax - MinCommonLz) + 2)];
  for (int i = MaxPartialNbrFermions; i >= 0; --i)
    {
      long TmpPartialHilbertSpaceDimension;
      int TmpMin = ((i * (i - 1)) >> 1);
      int TmpMax = ((MaxPartialNbrFermions - 1) * i) - ((i * (i - 1)) >> 1);
      int ShiftedTotalLz2 = ShiftedTotalLz - (i * this->StateHighestLzShift);
      int ShiftedTotalLz3  = ShiftedTotalLz2 - (((this->NbrFermions - i) * ((2 * MinCommonLz) - (this->NbrFermions - i) + 1)) >> 1);
      int ShiftedTotalLz4 = ShiftedTotalLz2 - (((this->NbrFermions - i) * (this->NbrFermions - i - 1)) >> 1);
      if (TmpMin < ShiftedTotalLz3)
	TmpMin = ShiftedTotalLz3;
      if (TmpMax > ShiftedTotalLz4)	
	TmpMax = ShiftedTotalLz4;
      for (; TmpMax >= TmpMin; --TmpMax)	  
	{
	  TmpPartialHilbertSpaceDimension = this->RawGenerateStates(i, MaxPartialNbrFermions - 1, MaxPartialNbrFermions - 1, TmpMax, this->PartialHilbertSpaceDimension);
	  for (; this->PartialHilbertSpaceDimension < TmpPartialHilbertSpaceDimension; ++this->PartialHilbertSpaceDimension)
	    this->StateDescriptionFileIndex[this->PartialHilbertSpaceDimension] = this->NbrFiles;
	  ++this->NbrFiles;	  
	}      
    }

  SortArrayDownOrdering(this->StateDescription, this->StateDescriptionFileIndex, this->PartialHilbertSpaceDimension);
  this->StateDescriptionIndexShift[0] = 0l;
  this->StateHighestLzToIndex[this->StateDescription[0]] = 0;
  for (long i = 1; i < this->PartialHilbertSpaceDimension; ++i)
    {
      this->StateDescriptionIndexShift[i] = this->StateDescriptionIndexShift[i - 1] + this->StateDescriptionFileSizes[this->StateDescriptionFileIndex[i - 1]];
      this->StateHighestLzToIndex[this->StateDescription[i]] = i;
    }
  this->HighestLzStateMask = ~((0x1l << this->StateHighestLzShift) - 1);
  this->NbrBuffers = memory / (LargestFile << 3);
  if (this->NbrBuffers < 2)
    this->NbrBuffers = 2;
  this->StateDescriptionBuffers = new unsigned long* [this->NbrBuffers];
  this->StateDescriptionLzSectorBuffers = new unsigned long* [this->NbrBuffers];
  for (int i = 0; i < this->NbrBuffers; ++i)
    {
      this->StateDescriptionBuffers[i] = new unsigned long[LargestFile];
      this->StateDescriptionLzSectorBuffers[i] = new unsigned long [this->StateHighestLzShift + 1];
    }
  this->StateDescriptionFileNames = new char*[this->NbrFiles];
  int CurrentFileIndex = 0;
  unsigned long* TmpStateDescription = this->StateDescription;
  for (int i = MaxPartialNbrFermions; i >= 0; --i)
    {
      int TmpMin = ((i * (i - 1)) >> 1);
      int TmpMax = ((MaxPartialNbrFermions - 1) * i) - ((i * (i - 1)) >> 1);
      int ShiftedTotalLz2 = ShiftedTotalLz - (i * this->StateHighestLzShift);
      int ShiftedTotalLz3  = ShiftedTotalLz2 - (((this->NbrFermions - i) * ((2 * MinCommonLz) - (this->NbrFermions - i) + 1)) >> 1);
      int ShiftedTotalLz4 = ShiftedTotalLz2 - (((this->NbrFermions - i) * (this->NbrFermions - i - 1)) >> 1);
      if (TmpMin < ShiftedTotalLz3)
	TmpMin = ShiftedTotalLz3;
      if (TmpMax > ShiftedTotalLz4)	
	TmpMax = ShiftedTotalLz4;
      for (; TmpMax >= TmpMin; --TmpMax)	  
	{
	  this->StateDescription = this->StateDescriptionBuffers[0];
	  this->RawGenerateStates(this->NbrFermions - i, MinCommonLz, MinCommonLz, ShiftedTotalLz2 - TmpMax, 0);
	  this->FindLzMaxSectors(this->StateDescription, this->StateDescriptionFileSizes[CurrentFileIndex], this->StateDescriptionLzSectorBuffers[0], this->StateHighestLzShift - 1);
	  this->StateDescriptionFileNames[CurrentFileIndex] = new char[256];
	  sprintf (this->StateDescriptionFileNames[CurrentFileIndex], "fermions_sphere_n_%d_2s_%d.%d.hb", this->NbrFermions, this->LzMax, CurrentFileIndex);	  
	  ofstream File; 
	  File.open(this->StateDescriptionFileNames[CurrentFileIndex], ios::binary | ios::out);
	  for (int i = 0; i <= this->StateHighestLzShift; ++i)
	    WriteLittleEndian(File, this->StateDescriptionLzSectorBuffers[0][i]);
	  for (unsigned long i = 0; i < this->StateDescriptionFileSizes[CurrentFileIndex]; ++i)
	    WriteLittleEndian(File, this->StateDescription[i]);
	  File.close();
	  ++CurrentFileIndex;
	}      
    }
  this->StateDescription = TmpStateDescription;
  this->BufferIndices = new int [this->NbrBuffers];
  this->BufferAges = new int [this->NbrBuffers];
  for (int i = 1; i < this->NbrBuffers; ++i)
    {
      this->BufferIndices[i] = -1;
      this->BufferAges[i] = this->NbrBuffers + 1;
    }
  this->BufferIndices[0] = this->NbrFiles - 1;
  this->BufferAges[0] = 1;
  this->FileToBuffer = new int [this->NbrFiles];
  for (int i = 0; i < this->NbrFiles; ++i)
    this->FileToBuffer[i] = -1;
  this->FileToBuffer[this->NbrFiles - 1] = 0;

#ifdef  __64_BITS__
   unsigned long ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 6) + 1l;
#else
   unsigned long ReducedHilbertSpaceDimension = (this->LargeHilbertSpaceDimension >> 5) + 1l;
#endif
   this->KeepStateFlag = new unsigned long [ReducedHilbertSpaceDimension];
   for (unsigned long i = 0; i < ReducedHilbertSpaceDimension; ++i)
     this->KeepStateFlag[i] = 0x0l;

   int MaxSweeps = (this->NbrFermions * (this->NbrFermions - 1)) >> 1;  
   this->TmpGeneratedStates =  new unsigned long [MaxSweeps * 1000];
   this->TmpGeneratedStatesLzMax = new int [MaxSweeps * 1000];
   long Memory = 0l;

   unsigned long TmpIndex = this->FindStateIndex(this->ReferenceState);
#ifdef  __64_BITS__
   this->KeepStateFlag[TmpIndex >> 6] = 0x1l << (TmpIndex & 0x3f);
#else
   this->KeepStateFlag[TmpIndex >> 5] = 0x1l << (TmpIndex & 0x1f);
#endif
   this->GenerateStates(ReferenceStateLzMax, this->ReferenceState, 1l, Memory);

   if ((this->SymmetricFlag == false) || (this->TotalLz != 0))
     {
       unsigned long NewHilbertSpaceDimension = 0;
       unsigned long TmpKeepStateFlag;
       int TmpNbrOne[] = {  
	 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 
	 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
	 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
	 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
	 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
	 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
	 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
	 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
	 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 
	 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
	 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
	 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
	 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 
	 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
	 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 
	 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8};
       for (unsigned long i = 0; i < ReducedHilbertSpaceDimension; ++i)
	 {
	   TmpKeepStateFlag = this->KeepStateFlag[i];
	   NewHilbertSpaceDimension += TmpNbrOne[TmpKeepStateFlag & 0xffl];
	   NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 8) & 0xffl];
	   NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 16) & 0xffl];
	   NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 24) & 0xffl];
#ifdef  __64_BITS__
	   NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 32) & 0xffl];
	   NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 40) & 0xffl];
	   NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 48) & 0xffl];
	   NewHilbertSpaceDimension += TmpNbrOne[(TmpKeepStateFlag >> 56) & 0xffl];      
#endif
	 }
       cout << "Haldane space dimension = " << NewHilbertSpaceDimension << endl;
       this->FullLargeHilbertSpaceDimension =  this->LargeHilbertSpaceDimension;
       this->LargeHilbertSpaceDimension = NewHilbertSpaceDimension;
       if (this->LargeHilbertSpaceDimension >= (1l << 31))
	 this->HilbertSpaceDimension = 0;
       else
	 this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
     }
   else
     {
       cout << "symmetrizing..." << endl;
#ifdef __64_BITS__
       this->InvertShift = 32 - ((this->LzMax + 1) >> 1);
#else
       this->InvertShift = 16 - ((this->LzMax + 1 ) >> 1);
#endif
       if ((this->LzMax & 1) == 0)
	 this->InvertUnshift = this->InvertShift - 1;
       else
	 this->InvertUnshift = this->InvertShift;
       unsigned long NewHilbertSpaceDimension = 0;
       unsigned long NewTotalHilbertSpaceDimension = 0;
       long TotalIndex = 0;
       long LocalIndex = 0;
       unsigned long TmpStateHighestLz = this->StateDescription[0] << this->StateHighestLzShift;
       int TmpFileIndex = this->StateDescriptionFileIndex[this->StateHighestLzToIndex[this->StateDescription[0]]];
       long TmpLimit = this->StateDescriptionFileSizes[TmpFileIndex];
       int TmpBufferIndex = this->FileToBuffer[TmpFileIndex];
       if (TmpBufferIndex < 0)
	 TmpBufferIndex = this->LoadLowestLzBuffer(TmpFileIndex);
       unsigned long* TmpStateDescriptionBuffers = this->StateDescriptionBuffers[TmpBufferIndex];
       for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	 {
#ifdef  __64_BITS__
	   unsigned long& TmpKeepStateFlag = this->KeepStateFlag[i >> 6];	   
#else
	   unsigned long& TmpKeepStateFlag = this->KeepStateFlag[i >> 5];	   
#endif	   
	   if (this->IsCanonicalState(TmpStateHighestLz | TmpStateDescriptionBuffers[LocalIndex]))
	     {
	       ++NewTotalHilbertSpaceDimension;
#ifdef  __64_BITS__
	       if (((TmpKeepStateFlag >> (i & 0x3ful)) & 0x1ul) != 0ul)
		 ++NewHilbertSpaceDimension;
	       else
		 TmpKeepStateFlag &= ~(0x1ul << (i & 0x3ful));
#else
	       if (((TmpKeepStateFlag >> (i & 0x1ful)) & 0x1ul) != 0ul)
		 ++NewHilbertSpaceDimension;
	       else
		 TmpKeepStateFlag &= ~(0x1ul << (i & 0x1ful));
#endif
	     }
	   ++LocalIndex;
	   if (LocalIndex == TmpLimit)
	     {
	       ++TotalIndex;
	       if (TotalIndex < this->PartialHilbertSpaceDimension)
		 {
		   TmpStateHighestLz = this->StateDescription[TotalIndex] << this->StateHighestLzShift;
		   TmpFileIndex = this->StateDescriptionFileIndex[this->StateHighestLzToIndex[this->StateDescription[TotalIndex]]];
		   TmpLimit = this->StateDescriptionFileSizes[TmpFileIndex];
		   TmpBufferIndex = this->FileToBuffer[TmpFileIndex];
		   if (TmpBufferIndex < 0)
		     TmpBufferIndex = this->LoadLowestLzBuffer(TmpFileIndex);
		   TmpStateDescriptionBuffers = this->StateDescriptionBuffers[TmpBufferIndex];
		   LocalIndex = 0ul;
		 }
	     }
	 }
       cout << "total symmetric space dimension = " << NewTotalHilbertSpaceDimension << endl;
       cout << "Haldane symmetric space dimension = " << NewHilbertSpaceDimension << endl;
       this->FullLargeHilbertSpaceDimension =  this->LargeHilbertSpaceDimension;
       this->LargeHilbertSpaceDimension = NewHilbertSpaceDimension;
       if (this->LargeHilbertSpaceDimension >= (1l << 31))
	 this->HilbertSpaceDimension = 0;
       else
	 this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
     }

  cout << "final Hilbert space dimension " << this->LargeHilbertSpaceDimension << endl;
}

// constructor from a binary file that describes the Hilbert space
//
// fileName = name of the binary file
// memoryHilbert = amount of memory granted to store the Hilbert space (in Mbytes)
// memory = amount of memory allowed for precalculations

FermionOnSphereHaldaneHugeBasis::FermionOnSphereHaldaneHugeBasis(char* fileName, long memoryHilbert, long memory)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      this->HilbertSpaceDimension = 0;
      return;
    }

  ReadLittleEndian(File, this->HilbertSpaceDimension);
  ReadLittleEndian(File, this->LargeHilbertSpaceDimension);
  cout << this->HilbertSpaceDimension << " " << this->LargeHilbertSpaceDimension << endl;
  ReadLittleEndian(File, this->NbrFermions);
  ReadLittleEndian(File, this->LzMax);
  ReadLittleEndian(File, this->TotalLz);
  ReadLittleEndian(File, this->ReferenceState);
  if ((this->LargeHilbertSpaceDimension << 3) < (memoryHilbert << 20))
    {
      this->StateDescription = new unsigned long [this->LargeHilbertSpaceDimension];
      for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
	ReadLittleEndian(File, this->StateDescription[i]);
    }
  else 
    this->StateDescription = 0;
  File.close();

#ifdef __64_BITS__
  this->InvertShift = 32 - ((this->LzMax + 1) >> 1);
#else
  this->InvertShift = 16 - ((this->LzMax + 1 ) >> 1);
#endif
  if ((this->LzMax & 1) == 0)
    this->InvertUnshift = this->InvertShift - 1;
  else
    this->InvertUnshift = this->InvertShift;

  this->TargetSpace = this;
  this->NbrLzValue = this->LzMax + 1;
  this->MaximumSignLookUp = 16;
  this->IncNbrFermions = this->NbrFermions + 1;
  this->Flag.Initialize();  

  if (this->StateDescription != 0)
    {
      this->GenerateLookUpTable(memory);
#ifdef __DEBUG__
      unsigned long UsedMemory = 0l;
      UsedMemory += ((unsigned long) this->LargeHilbertSpaceDimension) * (sizeof(unsigned long) + sizeof(int));
      UsedMemory += this->NbrLzValue * sizeof(int);
      UsedMemory += this->NbrLzValue * this->LookUpTableMemorySize * sizeof(int);
      UsedMemory +=  (1 << this->MaximumSignLookUp) * sizeof(double);
      cout << "memory requested for Hilbert space = ";
      if (UsedMemory >= 1024)
	if (UsedMemory >= 1048576)
	  cout << (UsedMemory >> 20) << "Mo" << endl;
	else
	  cout << (UsedMemory >> 10) << "ko" <<  endl;
      else
	cout << UsedMemory << endl;
#endif
    }
}



// copy constructor (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy

FermionOnSphereHaldaneHugeBasis::FermionOnSphereHaldaneHugeBasis(const FermionOnSphereHaldaneHugeBasis& fermions)
{
  this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->KeepStateFlag = fermions.KeepStateFlag;
}

// destructor
//

FermionOnSphereHaldaneHugeBasis::~FermionOnSphereHaldaneHugeBasis ()
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
      for (int i = 0; i < this->NbrBuffers; ++i)
	delete[] this->StateDescriptionBuffers[i];
      delete[] this->StateDescriptionBuffers;
      for (int i = 0; i < this->NbrFiles; ++i)
	delete[] this->StateDescriptionFileNames[i];
      delete[] this->StateDescriptionFileNames;
      delete[] this->StateDescriptionFileSizes;
      delete[] this->KeepStateFlag;
    }
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

FermionOnSphereHaldaneHugeBasis& FermionOnSphereHaldaneHugeBasis::operator = (const FermionOnSphereHaldaneHugeBasis& fermions)
{
  if ((this->HilbertSpaceDimension != 0) && (this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete[] this->StateDescription;
      delete[] this->SignLookUpTable;
      delete[] this->SignLookUpTableMask;
    }
  if (this->TargetSpace != &fermions)
    this->TargetSpace = fermions.TargetSpace;
  else
    this->TargetSpace = this;
  this->NbrFermions = fermions.NbrFermions;
  this->IncNbrFermions = fermions.IncNbrFermions;
  this->TotalLz = fermions.TotalLz;
  this->HilbertSpaceDimension = fermions.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = fermions.LargeHilbertSpaceDimension;
  this->StateDescription = fermions.StateDescription;
  this->LzMax = fermions.LzMax;
  this->NbrLzValue = fermions.NbrLzValue;
  this->Flag = fermions.Flag;
  this->SignLookUpTable = fermions.SignLookUpTable;
  this->SignLookUpTableMask = fermions.SignLookUpTableMask;
  this->MaximumSignLookUp = fermions.MaximumSignLookUp;
  this->KeepStateFlag = fermions.KeepStateFlag;
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* FermionOnSphereHaldaneHugeBasis::Clone()
{
  return new FermionOnSphereHaldaneHugeBasis(*this);
}

// save Hilbert space description to disk
//
// fileName = name of the file where the Hilbert space description has to be saved
// return value = true if no error occured

bool FermionOnSphereHaldaneHugeBasis::WriteHilbertSpace (char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out);
  if (!File.is_open())
    {
      cout << "can't open the file: " << fileName << endl;
      return false;
    }
  WriteLittleEndian(File, this->HilbertSpaceDimension);
  WriteLittleEndian(File, this->LargeHilbertSpaceDimension);
  WriteLittleEndian(File, this->NbrFermions);
  WriteLittleEndian(File, this->LzMax);
  WriteLittleEndian(File, this->TotalLz);
  WriteLittleEndian(File, this->ReferenceState);

  long TotalIndex = 0;
  long LocalIndex = 0;
  unsigned long TmpStateHighestLz = this->StateDescription[0] << this->StateHighestLzShift;
  int TmpFileIndex = this->StateDescriptionFileIndex[this->StateHighestLzToIndex[this->StateDescription[0]]];
  long TmpLimit = this->StateDescriptionFileSizes[TmpFileIndex];
  int TmpBufferIndex = this->FileToBuffer[TmpFileIndex];
  if (TmpBufferIndex < 0)
    TmpBufferIndex = this->LoadLowestLzBuffer(TmpFileIndex);
  unsigned long* TmpStateDescriptionBuffers = this->StateDescriptionBuffers[TmpBufferIndex];
  long Count = 0;;
  for (long i = 0; i < this->FullLargeHilbertSpaceDimension; ++i)
    {
#ifdef  __64_BITS__
      unsigned long& TmpKeepStateFlag = this->KeepStateFlag[i >> 6];	   
#else
      unsigned long& TmpKeepStateFlag = this->KeepStateFlag[i >> 5];	   
#endif	   
#ifdef  __64_BITS__
      if (((TmpKeepStateFlag >> (i & 0x3ful)) & 0x1ul) != 0ul)
	{
	  unsigned long TmpState = TmpStateHighestLz | TmpStateDescriptionBuffers[LocalIndex];
	  WriteLittleEndian(File, TmpState);
	  Count++;
	}
#else
      if (((TmpKeepStateFlag >> (i & 0x1ful)) & 0x1ul) != 0ul)
	{
	  unsigned long TmpState = TmpStateHighestLz | TmpStateDescriptionBuffers[LocalIndex];
	  WriteLittleEndian(File, TmpState);
	}
#endif
      ++LocalIndex;
      if (LocalIndex == TmpLimit)
	{
	  ++TotalIndex;
	  if (TotalIndex < this->PartialHilbertSpaceDimension)
	    {
	      TmpStateHighestLz = this->StateDescription[TotalIndex] << this->StateHighestLzShift;
	      TmpFileIndex = this->StateDescriptionFileIndex[this->StateHighestLzToIndex[this->StateDescription[TotalIndex]]];
	      TmpLimit = this->StateDescriptionFileSizes[TmpFileIndex];
	      TmpBufferIndex = this->FileToBuffer[TmpFileIndex];
	      if (TmpBufferIndex < 0)
		TmpBufferIndex = this->LoadLowestLzBuffer(TmpFileIndex);
	      TmpStateDescriptionBuffers = this->StateDescriptionBuffers[TmpBufferIndex];
	      LocalIndex = 0ul;
	    }
	}
    }  
  File.close();
  cout << "Nbr Saved states = " << Count << endl;
  return true;
}

// return a list of all possible quantum numbers 
//
// return value = pointer to corresponding quantum number

List<AbstractQuantumNumber*> FermionOnSphereHaldaneHugeBasis::GetQuantumNumbers ()
{
  List<AbstractQuantumNumber*> TmpList;
  TmpList += new SzQuantumNumber (this->TotalLz);
  return TmpList;
}

// return quantum number associated to a given state
//
// index = index of the state
// return value = pointer to corresponding quantum number

AbstractQuantumNumber* FermionOnSphereHaldaneHugeBasis::GetQuantumNumber (int index)
{
  return new SzQuantumNumber (this->TotalLz);
}

// extract subspace with a fixed quantum number
//
// q = quantum number value
// converter = reference on subspace-space converter to use
// return value = pointer to the new subspace

AbstractHilbertSpace* FermionOnSphereHaldaneHugeBasis::ExtractSubspace (AbstractQuantumNumber& q, 
							SubspaceSpaceConverter& converter)
{
  return 0;
}

// set a different target space (for all basic operations)
//
// targetSpace = pointer to the target space

void FermionOnSphereHaldaneHugeBasis::SetTargetSpace(ParticleOnSphere* targetSpace)
{
  this->TargetSpace = (FermionOnSphereHaldaneHugeBasis*) targetSpace;
}

// return Hilbert space dimension of the target space
//
// return value = Hilbert space dimension

int FermionOnSphereHaldaneHugeBasis::GetTargetHilbertSpaceDimension()
{
  return this->TargetSpace->HilbertSpaceDimension;
}

// apply a^+_m1 a^+_m2 a_n1 a_n2 operator to a given state (with m1+m2=n1+n2)
//
// index = index of the state on which the operator has to be applied
// m1 = first index for creation operator
// m2 = second index for creation operator
// n1 = first index for annihilation operator
// n2 = second index for annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereHaldaneHugeBasis::AdAdAA (int index, int m1, int m2, int n1, int n2, double& coefficient)
{
  return 0;
}

// apply Prod_i a^+_mi Prod_i a_ni operator to a given state (with Sum_i  mi= Sum_i ni)
//
// index = index of the state on which the operator has to be applied
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereHaldaneHugeBasis::ProdAdProdA (int index, int* m, int* n, int nbrIndices, double& coefficient)
{
  return 0;
}

// apply Prod_i a_ni operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next ProdA call
//
// index = index of the state on which the operator has to be applied
// n = array containg the indices of the annihilation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// return value =  multiplicative factor 

double FermionOnSphereHaldaneHugeBasis::ProdA (int index, int* n, int nbrIndices)
{
  return 0;
}

// apply Prod_i a^+_mi operator to the state produced using ProdA method (without destroying it)
//
// m = array containg the indices of the creation operators (first index corresponding to the leftmost operator)
// nbrIndices = number of creation (or annihilation) operators
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 

int FermionOnSphereHaldaneHugeBasis::ProdAd (int* m, int nbrIndices, double& coefficient)
{
  return 0;
}

// apply a^+_m a_m operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation and annihilation operator
// return value = coefficient obtained when applying a^+_m a_m

double FermionOnSphereHaldaneHugeBasis::AdA (int index, int m)
{
  return 0;
}

// apply a^+_m a_n operator to a given state 
//
// index = index of the state on which the operator has to be applied
// m = index of the creation operator
// n = index of the annihilation operator
// coefficient = reference on the double where the multiplicative factor has to be stored
// return value = index of the destination state 
  
int FermionOnSphereHaldaneHugeBasis::AdA (int index, int m, int n, double& coefficient)
{
  return 0;
}

// find state index
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

long FermionOnSphereHaldaneHugeBasis::FindStateIndex(unsigned long stateDescription)
{
  unsigned long TmpHighestLz = stateDescription >> this->StateHighestLzShift;
  int TmpBufferIndex =  this->StateHighestLzToIndex[TmpHighestLz];
  TmpHighestLz = this->StateDescriptionIndexShift[TmpBufferIndex];
  int TmpFileIndex = this->StateDescriptionFileIndex[TmpBufferIndex];
  TmpBufferIndex = this->FileToBuffer[TmpFileIndex];
  if (TmpBufferIndex < 0)
    TmpBufferIndex = this->LoadLowestLzBuffer(TmpFileIndex);

  stateDescription &= ~this->HighestLzStateMask;
  int TmpLzMax = this->StateHighestLzShift;
  while (((stateDescription >> TmpLzMax) & 0x1l) == 0l)
    --TmpLzMax;
  long PosMax = this->StateDescriptionLzSectorBuffers[TmpBufferIndex][TmpLzMax + 1];
  long PosMin = this->StateDescriptionLzSectorBuffers[TmpBufferIndex][TmpLzMax] - 1l;
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long* TmpStateDescriptionBuffers = this->StateDescriptionBuffers[TmpBufferIndex];
  unsigned long CurrentState = TmpStateDescriptionBuffers[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	PosMax = PosMid;
      else
	PosMin = PosMid;
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = TmpStateDescriptionBuffers[PosMid];
    }
  if (CurrentState == stateDescription)
    return (TmpHighestLz + PosMid);
  else
    return (TmpHighestLz + PosMin);
}

// load a lowest Lz part file into memory
//
// fileIndex = index of the file to read
// return file = index of the buffer that has been loaded

int FermionOnSphereHaldaneHugeBasis::LoadLowestLzBuffer(int fileIndex)
{
  int Index = 0;  
  while (this->BufferAges[Index] < this->NbrBuffers)
    ++Index;
  if (this->BufferIndices[Index] >= 0)
    this->FileToBuffer[this->BufferIndices[Index]] = -1;
  unsigned long* TmpBuffer = this->StateDescriptionBuffers[Index];
  unsigned long* TmpBuffer2 = this->StateDescriptionLzSectorBuffers[Index];
  unsigned long TmpSize = this->StateDescriptionFileSizes[fileIndex];
  ifstream File;
  File.open(this->StateDescriptionFileNames[fileIndex], ios::binary | ios::in);
  for (int i = 0; i <= this->StateHighestLzShift; ++i)
    ReadLittleEndian(File, TmpBuffer2[i]);
  for (unsigned long i = 0; i < TmpSize; ++i)
    ReadLittleEndian(File, TmpBuffer[i]);
  File.close();
  this->BufferAges[Index] = 0;
  this->BufferIndices[Index] = fileIndex;
  for (int i = 0; i < this->NbrBuffers; ++i)
    ++this->BufferAges[i];
  this->FileToBuffer[this->BufferIndices[Index]] = Index;
  return Index;
}

// find state index assuming the whole Hilbert space is stored in memory
//
// stateDescription = unsigned integer describing the state
// lzmax = maximum Lz value reached by a fermion in the state
// return value = corresponding index

long FermionOnSphereHaldaneHugeBasis::FindStateIndexMemory(unsigned long stateDescription, int lzmax)
{
  long PosMax = stateDescription >> this->LookUpTableShift[lzmax];
  long PosMin = this->LookUpTable[lzmax][PosMax];
  PosMax = this->LookUpTable[lzmax][PosMax + 1];
  long PosMid = (PosMin + PosMax) >> 1;
  unsigned long CurrentState = this->StateDescription[PosMid];
  while ((PosMax != PosMid) && (CurrentState != stateDescription))
    {
      if (CurrentState > stateDescription)
	{
	  PosMax = PosMid;
	}
      else
	{
	  PosMin = PosMid;
	} 
      PosMid = (PosMin + PosMax) >> 1;
      CurrentState = this->StateDescription[PosMid];
    }
  if (CurrentState == stateDescription)
    return PosMid;
  else
    return PosMin;
}

// print a given State
//
// Str = reference on current output stream 
// state = ID of the state to print
// return value = reference on current output stream 

ostream& FermionOnSphereHaldaneHugeBasis::PrintState (ostream& Str, int state)
{
  unsigned long TmpState = this->StateDescription[state];
  for (int i = 0; i < this->NbrLzValue; ++i)
    Str << ((TmpState >> i) & ((unsigned long) 0x1)) << " ";
  return Str;
}

// generate all states corresponding to the constraints
// 
// lzMax = momentum maximum value for a fermion in the state
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereHaldaneHugeBasis::GenerateStates(int lzMax, unsigned long referenceState, long pos, long& memory)
{
  int MaxSweeps = (this->NbrFermions * (this->NbrFermions - 1)) >> 1;  
  unsigned long* TmpGeneratedStates2 = this->TmpGeneratedStates + (MaxSweeps * memory);
  int* TmpLzMax = this->TmpGeneratedStatesLzMax  + (MaxSweeps  * memory);  
  memory += 1;
  int TmpCurrentLzMax = 2;
  int TmpCurrentLzMax2;
  int TmpMax = lzMax - 1;
  int NbrEntries = 0;
  unsigned long TmpReferenceState;
  
  while (TmpCurrentLzMax < TmpMax)
    {
      while ((TmpCurrentLzMax < TmpMax) && (((referenceState >> TmpCurrentLzMax) & 0x3l) != 0x2l))
	++TmpCurrentLzMax;
      if (TmpCurrentLzMax < TmpMax)
	{
	  TmpReferenceState = (referenceState & ~(0x3l << TmpCurrentLzMax)) | (0x1l << TmpCurrentLzMax);
	  TmpCurrentLzMax2 = TmpCurrentLzMax - 2;
	  while (TmpCurrentLzMax2 >= 0)
	    {
	      while ((TmpCurrentLzMax2 >= 0) && (((referenceState >> TmpCurrentLzMax2) & 0x3l) != 0x1l))
		--TmpCurrentLzMax2;
	      if (TmpCurrentLzMax2 >= 0)
		{
		  TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0x3l << TmpCurrentLzMax2)) | (0x2l << TmpCurrentLzMax2);
		  TmpLzMax[NbrEntries] = lzMax;
		  ++NbrEntries;
		  --TmpCurrentLzMax2;
		}	      
	    }
	  ++TmpCurrentLzMax;
	}
    }
  if (((referenceState >> TmpCurrentLzMax) & 0x3l) == 0x2l)
    {
      TmpReferenceState = (referenceState & ~(0x3l << TmpCurrentLzMax)) | (0x1l << TmpCurrentLzMax);
      TmpCurrentLzMax2 = TmpCurrentLzMax - 2;
      while (TmpCurrentLzMax2 >= 0)
	{
	  while ((TmpCurrentLzMax2 >= 0) && (((referenceState >> TmpCurrentLzMax2) & 0x3l) != 0x1l))
	    --TmpCurrentLzMax2;
	  if (TmpCurrentLzMax2 >= 0)
	    {
	      TmpGeneratedStates2[NbrEntries] = (TmpReferenceState & ~(0x3l << TmpCurrentLzMax2)) | (0x2l << TmpCurrentLzMax2);
	      TmpLzMax[NbrEntries] = lzMax - 1;
	      ++NbrEntries;
	      --TmpCurrentLzMax2;
	    }
	}      
    }

  unsigned long TmpIndex;
  int NbrNewEntries = 0;
  for (int i = 0; i < NbrEntries; ++i)
    {
      TmpIndex = this->FindStateIndex(TmpGeneratedStates2[i]);
#ifdef __64_BITS__
      if ((this->KeepStateFlag[TmpIndex >> 6] >> (TmpIndex & 0x3f)) & 0x1l)
	{
	  TmpGeneratedStates2[i] = 0x0l;
	}
      else
	{
	  this->KeepStateFlag[TmpIndex >> 6] |= 0x1l << (TmpIndex & 0x3f);
	  ++NbrNewEntries;
	}
#else
      if ((this->KeepStateFlag[TmpIndex >> 5] >> (TmpIndex & 0x1f)) & 0x1l)
	{
	  TmpGeneratedStates2[i] = 0x0l;
	}
      else
	{
	  this->KeepStateFlag[TmpIndex >> 5] |= 0x1l << (TmpIndex & 0x1f);
	  ++NbrNewEntries;
	}      
#endif
    }

  if (NbrNewEntries > 0)
    for (int i = 0; i < NbrEntries; ++i)
      if (TmpGeneratedStates2[i] != 0x0l)
	pos = this->GenerateStates(TmpLzMax[i], TmpGeneratedStates2[i], pos, memory);

  memory -= 1;
  return pos;
}


// generate all states (i.e. all possible skew symmetric polynomials with fixed Lz)
// 
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// currentLzMax = momentum maximum value for fermions that are still to be placed
// totalLz = momentum total value
// pos = position in StateDescription array where to store states
// return value = position from which new states have to be stored

long FermionOnSphereHaldaneHugeBasis::RawGenerateStates(int nbrFermions, int lzMax, int currentLzMax, int totalLz, long pos)
{
  if ((nbrFermions == 0) && (totalLz == 0))
    {
      this->StateDescription[pos] = 0l;
      return pos + 1l;
    }
  if ((nbrFermions == 0) || (totalLz < 0) || (currentLzMax < (nbrFermions - 1)))
    return pos;
  int LzTotalMax = ((2 * currentLzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (totalLz > LzTotalMax)
    return pos;
  if ((nbrFermions == 1) && (currentLzMax >= totalLz))
    {
      this->StateDescription[pos] = (0x1l << totalLz);
      return pos + 1l;
    }
  if (LzTotalMax == totalLz)
    {
      unsigned long Mask = 0l;
      for (int i = currentLzMax - nbrFermions + 1; i <= currentLzMax; ++i)
	Mask |= (((unsigned long) 1) << i);
      this->StateDescription[pos] = Mask;
      return pos + 1l;
    }

  int ReducedCurrentLzMax = currentLzMax - 1;
  long TmpPos = this->RawGenerateStates(nbrFermions - 1, lzMax, ReducedCurrentLzMax, totalLz - currentLzMax, pos);
  unsigned long Mask = ((unsigned long) 1) << currentLzMax;
  for (long i = pos; i < TmpPos; i++)
    this->StateDescription[i] |= Mask;
  if (lzMax == currentLzMax)
    return this->RawGenerateStates(nbrFermions, ReducedCurrentLzMax, ReducedCurrentLzMax, totalLz, TmpPos);
  else
    return this->RawGenerateStates(nbrFermions, lzMax, ReducedCurrentLzMax, totalLz, TmpPos);
}

// generate look-up table associated to current Hilbert space
// 
// memory = memory size that can be allocated for the look-up table

void FermionOnSphereHaldaneHugeBasis::GenerateLookUpTable(unsigned long memory)
{
  // evaluate look-up table size
  memory /= (sizeof(int*) * this->NbrLzValue);
  this->MaximumLookUpShift = 1;
  while (memory > 0)
    {
      memory >>= 1;
      ++this->MaximumLookUpShift;
    }
  if (this->MaximumLookUpShift > this->NbrLzValue)
    this->MaximumLookUpShift = this->NbrLzValue;
  this->LookUpTableMemorySize = 1 << this->MaximumLookUpShift;

  // construct  look-up tables for searching states
  this->LookUpTable = new long* [this->NbrLzValue];
  this->LookUpTableShift = new int [this->NbrLzValue];
  for (int i = 0; i < this->NbrLzValue; ++i)
    this->LookUpTable[i] = new long [this->LookUpTableMemorySize + 1];
  int CurrentLzMax = this->LzMax;
  unsigned long CurrentState = this->StateDescription[0];
  while (((CurrentState >> CurrentLzMax) & 0x1ul) == 0x0ul)
    --CurrentLzMax;
  int TmpLzMax = CurrentLzMax;
  long* TmpLookUpTable = this->LookUpTable[CurrentLzMax];
  if (CurrentLzMax < this->MaximumLookUpShift)
    this->LookUpTableShift[CurrentLzMax] = 0;
  else
    this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
  int CurrentShift = this->LookUpTableShift[CurrentLzMax];
  unsigned long CurrentLookUpTableValue = this->LookUpTableMemorySize;
  unsigned long TmpLookUpTableValue = CurrentState >> CurrentShift;
  while (CurrentLookUpTableValue > TmpLookUpTableValue)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = 0;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[CurrentLookUpTableValue] = 0;
  for (long i = 0; i < this->LargeHilbertSpaceDimension; ++i)
    {
      CurrentState = this->StateDescription[i];
      while (((CurrentState >> TmpLzMax) & 0x1ul) == 0x0ul)
	--TmpLzMax;
      if (CurrentLzMax != TmpLzMax)
	{
	  while (CurrentLookUpTableValue > 0)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[0] = i;
	  /*	  for (unsigned long j = 0; j <= this->LookUpTableMemorySize; ++j)
	    cout << TmpLookUpTable[j] << " ";
	    cout << endl << "-------------------------------------------" << endl;*/
 	  CurrentLzMax = TmpLzMax;
	  TmpLookUpTable = this->LookUpTable[CurrentLzMax];
	  if (CurrentLzMax < this->MaximumLookUpShift)
	    this->LookUpTableShift[CurrentLzMax] = 0;
	  else
	    this->LookUpTableShift[CurrentLzMax] = CurrentLzMax + 1 - this->MaximumLookUpShift;
	  CurrentShift = this->LookUpTableShift[CurrentLzMax];
	  TmpLookUpTableValue = CurrentState >> CurrentShift;
	  CurrentLookUpTableValue = this->LookUpTableMemorySize;
	  while (CurrentLookUpTableValue > TmpLookUpTableValue)
	    {
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	      --CurrentLookUpTableValue;
	    }
	  TmpLookUpTable[CurrentLookUpTableValue] = i;
	}
      else
	{
	  TmpLookUpTableValue = CurrentState >> CurrentShift;
	  if (TmpLookUpTableValue != CurrentLookUpTableValue)
	    {
	      while (CurrentLookUpTableValue > TmpLookUpTableValue)
		{
		  TmpLookUpTable[CurrentLookUpTableValue] = i;
		  --CurrentLookUpTableValue;
		}
	      TmpLookUpTable[CurrentLookUpTableValue] = i;
	    }
	}
    }
  while (CurrentLookUpTableValue > 0)
    {
      TmpLookUpTable[CurrentLookUpTableValue] = this->LargeHilbertSpaceDimension - 1l;
      --CurrentLookUpTableValue;
    }
  TmpLookUpTable[0] = this->LargeHilbertSpaceDimension - 1l;
  // look-up tables for evaluating sign when applying creation/annihilation operators
  int Size = 1 << this->MaximumSignLookUp;
  this->SignLookUpTable = new double [Size];
  int Count;
  int TmpNbr;
  for (int j = 0; j < Size; ++j)
    {
      Count = 0;
      TmpNbr = j;
      while (TmpNbr != 0)
	{
	  if (TmpNbr & 0x1)
	    ++Count;
	  TmpNbr >>= 1;
	}
      if (Count & 1)
	this->SignLookUpTable[j] = -1.0;
      else
	this->SignLookUpTable[j] = 1.0;
    }
#ifdef __64_BITS__
  this->SignLookUpTableMask = new unsigned long [128];
  for (int i = 0; i < 48; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 48; i < 64; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 48);
  for (int i = 64; i < 128; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#else
  this->SignLookUpTableMask = new unsigned long [64];
  for (int i = 0; i < 16; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0xffff;
  for (int i = 16; i < 32; ++i)
    this->SignLookUpTableMask[i] = ((unsigned long) 0xffff) >> (i - 16);
  for (int i = 32; i < 64; ++i)
    this->SignLookUpTableMask[i] = (unsigned long) 0;
#endif
}

// evaluate Hilbert space dimension
//
// nbrFermions = number of fermions
// lzMax = momentum maximum value for a fermion
// totalLz = momentum total value
// return value = Hilbert space dimension

long FermionOnSphereHaldaneHugeBasis::EvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  return this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax, (totalLz + nbrFermions * lzMax) >> 1);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

long FermionOnSphereHaldaneHugeBasis::ShiftedEvaluateHilbertSpaceDimension(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (lzMax < (nbrFermions - 1)))
    return 0l;
  int LzTotalMax = ((2 * lzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (LzTotalMax < totalLz)
    return 0l;
  if ((nbrFermions == 1) && (lzMax >= totalLz))
    return 1l;
  if (LzTotalMax == totalLz)
    return 1l;
  long Tmp1 = this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax);
  long Tmp2 = this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz);  
  if ((Tmp1 < this->SizeLimit) && (Tmp2 < this->SizeLimit) && ((Tmp1 + Tmp2) > this->SizeLimit))
    {
      this->FileSizes += ((((unsigned long) (nbrFermions - 1)) << 8) | ((unsigned long) (lzMax - 1)) | (((unsigned long) (totalLz - lzMax)) << 16));
      this->FileSizes += ((((unsigned long) nbrFermions) << 8) | ((unsigned long) (lzMax - 1)) | (((unsigned long) totalLz) << 16));
    }
  return (Tmp1 + Tmp2);
}

// evaluate Hilbert space dimension with shifted values for lzMax and totalLz
//
// nbrFermions = number of fermions
// lzMax = two times momentum maximum value for a fermion plus one 
// totalLz = momentum total value plus nbrFermions * (momentum maximum value for a fermion + 1)
// return value = Hilbert space dimension

long FermionOnSphereHaldaneHugeBasis::ShiftedEvaluateHilbertSpaceDimension2(int nbrFermions, int lzMax, int totalLz)
{
  if ((nbrFermions == 0) || (totalLz < 0)  || (lzMax < (nbrFermions - 1)))
    return 0l;
  int LzTotalMax = ((2 * lzMax - nbrFermions + 1) * nbrFermions) >> 1;
  if (LzTotalMax < totalLz)
    return 0l;
  if ((nbrFermions == 1) && (lzMax >= totalLz))
    return 1l;
  if (LzTotalMax == totalLz)
    return 1l;
  return (this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions - 1, lzMax - 1, totalLz - lzMax) +
	  this->ShiftedEvaluateHilbertSpaceDimension(nbrFermions, lzMax - 1, totalLz));  
}

// find position of each sector with a given maximum Lz
// 
// stateDescription = array that contains state description
// dimension = dimension of the stateDescription array
// lzSectors = array that where the position of each sector with a given maximum Lz will be stored
// lzMax = maximum momentum value for a fermion

void FermionOnSphereHaldaneHugeBasis::FindLzMaxSectors(unsigned long* stateDescription, long dimension, unsigned long* lzSectors, int lzMax)
{
  for (int i = 0; i <= lzMax; ++i)
    lzSectors[i] = 0xffffffffl;
  unsigned long TmpStateDescription;
  int CurrentLzMax = lzMax;
  int TmpLzMax;
  long Position = 0l;
  lzSectors[lzMax + 1] = 0l;
  while (Position < dimension)
    {
      TmpLzMax = CurrentLzMax;
      TmpStateDescription = stateDescription[Position];
      while (((TmpStateDescription >> TmpLzMax) & 0x1l) == 0x0l)
	--TmpLzMax;
      if (CurrentLzMax != TmpLzMax)
	{
	  lzSectors[CurrentLzMax] = Position;
	  lzSectors[TmpLzMax + 1] = Position;
	  CurrentLzMax = TmpLzMax;
	}
      ++Position;
    }
  lzSectors[CurrentLzMax] = dimension;
}

// evaluate wave function in real space using a given basis and only for agiven range of components
//
// state = vector corresponding to the state in the Fock basis
// position = vector whose components give coordinates of the point where the wave function has to be evaluated
// basis = one body real space basis to use
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = wave function evaluated at the given location

Complex FermionOnSphereHaldaneHugeBasis::EvaluateWaveFunction (RealVector& state, RealVector& position, AbstractFunctionBasis& basis,
							       int firstComponent, int nbrComponent)
{
  Complex Value;
  return Value;
}

// initialize evaluation of wave function in real space using a given basis and only for a given range of components and
//
// timeCoherence = true if time coherence has to be used

void FermionOnSphereHaldaneHugeBasis::InitializeWaveFunctionEvaluation (bool timeCoherence)
{
}
  
