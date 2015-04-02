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

#ifndef INDEX_H_
#define INDEX_H_

#include <iostream>
#include <stdint.h>

#include "common.h"

using namespace std;

class Index {
private:
  uint32_t seed_;
  uint32_t keys_count_length_;
  uint32_t positions_length_;
  uint32_t *keys_count_;
  uint32_t *positions_;

public:
  uint32_t GetSeed
  (
      // no parameter
  )
  {
    return seed_;
  }

  uint32_t *GetKeysCount
  (
      // no parameter
  )
  {
    return keys_count_;
  }

  uint32_t *GetAllPositions
  (
      // no parameter
  )
  {
    return positions_;
  }

  uint32_t GetKeysCountLength
  (
      // no parameter
  )
  {
    return keys_count_length_;
  }

  uint32_t GetPositionsLength
  (
      // no parameter
  )
  {
    return positions_length_;
  }

  // return key
  uint32_t GetKey
  (
      uint8_t sequence[]
  )
  {
    uint32_t i;
    uint32_t seed;
    uint32_t key;
    for (i = 0, seed = seed_, key = 0; seed != 0; ++i, seed >>= 1) {
      if (seed & 1) {
        key = key << CHARACTER_SIZE;
        key = key | sequence[i];
      }
    }
    return key;
  }

  // get positions and positions length
  // return positions
  uint32_t *GetPositions
  (
      uint8_t sequence[], // read sequence
      uint32_t *positions_length // store positions length
  )
  {
    uint32_t key = GetKey(sequence);
    *positions_length = keys_count_[key + 1] - keys_count_[key];
    return &positions_[keys_count_[key]];
  }

  uint32_t GetSeedLength
  (
      // no parameter
  )
  {
    uint32_t seed;
    uint32_t length;
    for (length = 0, seed = seed_; seed != 0; seed >>= 1, ++length)
      ;
    return length;
  }

  uint32_t GetSeedWeight
  (
      // no parameter
  )
  {
    uint32_t seed;
    uint32_t weight;
    for (weight = 0, seed = seed_; seed != 0; seed >>= 1) {
      if (seed & 1) {
        ++ weight;
      }
    }
    return weight;
  }

  // constructor
  Index
  (
      uint32_t seed,
      uint32_t keys_count_length,
      uint32_t positions_length,
      uint32_t keys_count[],
      uint32_t positions[]
  ){
    seed_ = seed;
    keys_count_length_ = keys_count_length;
    positions_length_ = positions_length;
    keys_count_ =keys_count;
    positions_ = positions;
  }

  // destructor
  virtual ~Index
  (
      //  no parameter
  ){
    delete [] keys_count_;
    delete [] positions_;
  }
};

#endif /* INDEX_H_ */
