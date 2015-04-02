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

#ifndef QUERY_H_
#define QUERY_H_

#include <string>
#include <stdint.h>

using namespace std;

class Query {
private:
  uint32_t number_sequences_;
  uint32_t sequence_length_;
  string *names_;
  string query_prefix_;
  uint8_t *sequences_;

  void SetNames
  (
      // no parameter
  );

  void SetSequences
  (
      // no parameter
  );

public:
  Query
  (
      string query_prefix,
      uint32_t number_sequences,
      uint32_t sequence_length
  )
  {
    query_prefix_ = query_prefix;
    number_sequences_ = number_sequences;
    sequence_length_ = sequence_length;
    names_ = NULL;
    sequences_ = NULL;
  }
  virtual ~Query
  (
      // no parameter
  )
  {
    if (names_ != NULL) {
      delete [] names_;
    }

    if (sequences_ != NULL) {
      delete [] sequences_;
    }
  }

  uint32_t GetNumberSequences
  (
      // no parameter
  )
  {
    return number_sequences_;
  }

  uint32_t GetSequenceLength
  (
      // no parameter
  )
  {
    return sequence_length_;
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

};

#endif /* QUERY_H_ */
