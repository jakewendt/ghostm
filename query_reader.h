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

#ifndef QUERY_READER_H_
#define QUERY_READER_H_

#include <string>

using namespace std;

class Query;

class QueryReader {
public:
  // constructor
  QueryReader
  (
      string query_prefix
  );

  Query *Read
  (
      uint32_t id
  );

  Query *Read
  (
      // no parameter
  );

  uint32_t GetNextId() {
    return next_id_;
  }

  uint32_t GetMaxLengthConcatenatedSequence() {
    return max_length_sequence_;
  }

  uint32_t GetMaxNumberSequences() {
    return max_number_sequeces_;
  }

  void Close()
  {}

  // use default destructor
  //virtual ~QueryReader();

private:
  string query_prefix_;
  uint32_t division_;
  uint32_t max_length_sequence_;
  uint32_t max_number_sequeces_;
  uint32_t next_id_;
};

#endif /* QUERY_READER_H_ */
