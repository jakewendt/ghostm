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

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <stdint.h>
#include "query.h"
#include "query_reader.h"

using namespace std;
QueryReader::QueryReader(string query_prefix)
:query_prefix_(query_prefix), division_(0),
 max_length_sequence_(0),max_number_sequeces_(0),next_id_(0)
{
  ifstream in;
  string filename = query_prefix_ + ".inf";
  in.open(filename.c_str(), ios::binary);
  if (in) {
    in.read((char *) &division_, sizeof(division_));
    in.read((char *) &max_length_sequence_, sizeof(max_length_sequence_));
    in.read((char *) &max_number_sequeces_, sizeof(max_number_sequeces_));
    in.close();
  }
}

Query *QueryReader::Read
(
 uint32_t id
)
{
  if (id < division_) {
    uint32_t number_sequences;
    uint32_t sequence_length;
    stringstream query_prefix;
    query_prefix.str("");
    query_prefix <<  query_prefix_ << "_" << id;
    string filename = query_prefix.str() + ".inf";
    ifstream in;
    in.open(filename.c_str(), ios::binary);
    if (!in) {
      return NULL;
    }
    in.read((char *) &number_sequences, sizeof(uint32_t));
    in.read((char *) &sequence_length, sizeof(uint32_t));
    in.close();
    next_id_ = id + 1;
    return new Query(query_prefix.str(), number_sequences, sequence_length);
  } else {
    return NULL;
  }
}

Query *QueryReader::Read
(
    // no parameter
)
{
  if (next_id_ < division_) {
    uint32_t number_sequences;
    uint32_t sequence_length;
    stringstream query_prefix;
    query_prefix.str("");
    query_prefix <<  query_prefix_ << "_" << next_id_;
    string filename = query_prefix.str() + ".inf";
    ifstream in;
    in.open(filename.c_str(), ios::binary);
    if (!in) {
      return NULL;
    }
    in.read((char *) &number_sequences, sizeof(uint32_t));
    in.read((char *) &sequence_length, sizeof(uint32_t));
    in.close();
    ++next_id_;
    return new Query(query_prefix.str(), number_sequences, sequence_length);
  } else {
    return NULL;
  }
}
