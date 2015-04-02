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
#include "db.h"
#include "db_reader.h"

DBReader::DBReader
(
    string db_prefix
)
:db_prefix_(db_prefix), division_(0),seed_(0), max_length_concatenated_sequence_(0), sum_length_(0), next_id_(0)
{
  ifstream in;
  string filename = db_prefix_ + ".inf";
  in.open(filename.c_str(), ios::binary);
  if (in) {
    in.read((char *) &division_, sizeof(division_));
    in.read((char *) &seed_, sizeof(seed_));
    in.read((char *) &max_length_concatenated_sequence_, sizeof(max_length_concatenated_sequence_));
    in.read((char *) &sum_length_, sizeof(sum_length_));
    in.close();
  }
}

DB *DBReader::Read
(
    // no parameter
)
{
  if (next_id_ < division_) {
    uint32_t number_sequences;
    uint32_t sequences_length;
    stringstream prefix;
    prefix.str("");
    prefix <<  db_prefix_ << "_" << next_id_;
    string file_name = prefix.str() + ".inf";
    ifstream in;
    in.open(file_name.c_str(), ios::binary);
    if (!in) {
      return NULL;
    }
    in.read((char *) &number_sequences, sizeof(uint32_t));
    in.read((char *) &sequences_length, sizeof(uint32_t));
    in.close();
    ++next_id_;
    return new DB(prefix.str(), number_sequences, sequences_length);
  } else {
    return NULL;
  }
}
