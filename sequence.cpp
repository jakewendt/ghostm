/* The MIT License

   Copyright (c) 2011 Akiyama_Laboratory , Tokyo Institute of Technology.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be 
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
   SOFTWARE.
*/

#include <stdint.h>

#include "common.h"
#include "sequence.h"

// TODO: amino and DNA

// DNA number
const uint8_t Sequence::kDNAToNumber[] = {
    // ASCII code to Number
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  //                                            -
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5, 4, 4,
  //
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  //@  A  B  C   D  E  F  G   H  I  J  K   L  M  N  O
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  //P  Q  R  S   T  U  V  W   X  Y  Z
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
  //'  a  b  c   d  e  f  g   h  i  j  k   l  m  n  o
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4,
  //p  q  r  s   t  u  v  w   x  y  z
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,

    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};


// amino number
const uint8_t Sequence::kProteinToNumber[] = {
    // ASCII code to Number. 25 is sequence end
     23,23,23,23, 23,23,23,23, 23,23,23,23, 23,23,23,23,
     23,23,23,23, 23,23,23,23, 23,23,23,23, 23,23,23,23,
  //                                  *         -
     23,23,23,23, 23,23,23,23, 23,23,24,23, 23, 23,23,23,
  //
     23,23,23,23, 23,23,23,23, 23,23,23,23, 23,23,23,23,
  //  @  A  B  C   D  E  F  G   H  I  J  K   L  M  N  O
     23, 0,20, 4,  3, 6,13, 7,  8, 9,21,11, 10,12, 2,23,
  //  P  Q  R  S   T  U  V  W   X  Y  Z
     14, 5, 1,15, 16,23,19,17, 23,18,22,23, 23,23,23,23,
  //  '  a  b  c   d  e  f  g   h  i  j  k   l  m  n  o
     23, 0,20, 4,  3, 6,13, 7,  8, 9,21,11, 10,12, 2,23,
  //  p  q  r  s   t  u  v  w   x  y  z
     14, 5, 1,15, 16,23,19,17, 23,18,22,23, 23,23,23,23,
     23,23,23,23, 23,23,23,23, 23,23,23,23, 23,23,23,23,
     23,23,23,23, 23,23,23,23, 23,23,23,23, 23,23,23,23,
     23,23,23,23, 23,23,23,23, 23,23,23,23, 23,23,23,23,
     23,23,23,23, 23,23,23,23, 23,23,23,23, 23,23,23,23,
     23,23,23,23, 23,23,23,23, 23,23,23,23, 23,23,23,23,
     23,23,23,23, 23,23,23,23, 23,23,23,23, 23,23,23,23,
     23,23,23,23, 23,23,23,23, 23,23,23,23, 23,23,23,23,
     23,23,23,23, 23,23,23,23, 23,23,23,23, 23,23,23,23
};

// return length of sequence to number. you need memory free.
uint8_t *Sequence::ToNumber
(
    // no parameter
)
{
  uint8_t *array;
  int size = sequence_.size();
  array = new uint8_t[size];

  if (type_ == PROTEIN) {
    for (int i = 0; i < size; ++i) {
      array[i] = Sequence::kProteinToNumber[(int)sequence_.at(i)];
    }
  } else if (type_ == DNA) {
    for (int i = 0; i < size; ++i) {
          array[i] = Sequence::kDNAToNumber[(int)sequence_.at(i)];
        }
  }

  return array;
}
