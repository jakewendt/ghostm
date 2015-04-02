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

#ifndef DB_H_
#define DB_H_

#include <string>
#include <stdint.h>
#include <limits.h>
#include "index.h"

using namespace std;

class DB {
private:
  string db_prefix_;
  string *names_;
  uint32_t number_sequences_;
  uint32_t sequences_length_;
  uint32_t *positions_;
  uint8_t *sequences_;
  Index *index_;

  void SetNames
  (
      // no paramater
  );

  void SetPositions
  (
      // no parameter
  );

  void SetIndex
  (
      // no parameter
  );

  void SetSequences
  (
      // no parameter
  );

public:
  uint32_t GetNumberSequences
  (
      // no parameter
  )
  {
    return number_sequences_;
  }

  uint32_t GetSequencesLength
  (
      // no patameter
  )
  {
    return sequences_length_;
  }

  string *GetAllNames
  (
      // no parameter
  )
  {
    if (names_ == NULL) {
      SetNames();
    }
    return names_;
  }


  string GetName
  (
      uint32_t id
  )
  {
    if (names_ == NULL) {
      SetNames();
    }
    return names_[id];
  }

  uint32_t GetID
  (
      uint32_t position
  )
  {
    if (positions_ == NULL) {
      SetPositions();
    }

    if (positions_[number_sequences_ - 1]  <= position && position < sequences_length_) {
      return number_sequences_ - 1;
    }

    uint32_t left = 0;
    uint32_t right = number_sequences_ - 2;
    uint32_t mid;
    while (left <= right) {
      mid = (left + right)/2;
      if (positions_[mid]  <= position && position < positions_[mid + 1]) {
        return mid;
      } else if (positions_[mid] < position) {
        left = mid + 1;
      } else {
        right = mid - 1;
      }
    }

    // no hit
    return UINT_MAX;
  }

  uint32_t *GetAllPositions
  (
      // no parameter
  )
  {
    if (positions_ == NULL) {
      SetPositions();
    }
    return positions_;
  }

  uint32_t GetPosition
  (
      uint32_t id
  )
  {
    if (positions_ == NULL) {
      SetPositions();
    }
    return positions_[id];
  }

  Index *GetIndex
  (
      // no parameter
  )
  {
    if (index_ == NULL) {
      SetIndex();
    }
    return index_;
  }

  uint8_t *GetSequences
  (
      // no parameter
  )
  {
    if (sequences_ == NULL) {
      SetSequences();
    }
    return sequences_;
  }


  DB
  (
      string db_prefix,
      uint32_t number_sequences,
      uint32_t sequences_length
  )
  {
    db_prefix_ = db_prefix;
    number_sequences_ = number_sequences;
    sequences_length_ = sequences_length;
    names_ = NULL;
    positions_ = NULL;
    sequences_ = NULL;
    index_ = NULL;
  }
  virtual ~DB
  (
      // no parameter
  )
  {
    if (names_ != NULL) {
      delete [] names_;
    }

    if (positions_ != NULL) {
      delete [] positions_;
    }

    if (index_ != NULL) {
      delete index_;
    }

    if (sequences_ != NULL) {
      delete [] sequences_;
    }
  }
};

#endif /* DB_H_ */
