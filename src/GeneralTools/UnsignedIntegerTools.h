////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2008 Gunnar Moller                   //
//                                                                            //
//                                                                            //
//          various unsigned long tools related to QHE Hilbert spaces         //
//                                                                            //
//                        last modification : 23/01/2006                      //
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


#ifndef UNSIGNED_INTEGER_TOOLS
#define UNSIGNED_INTEGER_TOOLS


// generate the smallest word with a given number of bits
//
// numBits = number of bits
// return value = smallest word
unsigned long smallestOne (int numBits);

// return next bigger word with same number of bits set.
//
// i = word from which next bigger word has to be built
// return value = next bigger word

unsigned long nextone (unsigned long i);

// generate a word with a given number of bits right of highest Bit,
// starting to number bits from the right with 1
//
// numBits = number of bits to write
// highestBit = poistion of the highest bit
// return value = corresponding word
unsigned long biggestOne (int numBits, int highestBit);

// return next smaller word with same number of bits set.
//
// i = word from which next smaller word has to be built
// return value = next smaller word
unsigned long lastone (unsigned long i);

// count bits in word
//
// i = word to test
// return value = number of bits
int bitcount (unsigned long i);


// count bits in word
//
// i = word to test
// return value = number of bits
inline double FastBitCount (unsigned long uLong)
{
#ifdef  __64_BITS__
  uLong ^= uLong >> 32;
#endif
  uLong ^= uLong >> 16;
  uLong ^= uLong >> 8;
  uLong ^= uLong >> 4;
  uLong ^= uLong >> 2;
  uLong ^= uLong << 1;
  uLong &= 0x2ul;
  return (1.0 - ((double) uLong));
}

// generates the smallest word with a given number of bits
//
// numBits = number of bits
// return value = smallest word

inline unsigned long smallestOne (int numBits)
{
  return ( (0x1ul << numBits) - 1); 
}

// generate a word with a given number of bits right of highest Bit,
// starting to number bits from the right with 1
//
//
// numBits = number of bits to write
// highestBit = poistion of the highest bit
// return value = corresponding word

inline unsigned long biggestOne (int numBits, int highestBit)
{
  if (numBits>highestBit) return 0ul;
  return ( ( ( (0x1ul << (highestBit-1)) - 1) | (0x1ul << (highestBit-1)) ) ^ ((0x1ul << (highestBit-numBits)) -1));
}

#endif
