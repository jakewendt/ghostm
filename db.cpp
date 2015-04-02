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

#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>
#include "index.h"
#include "db.h"

using namespace std;


void DB::SetNames
(
    // no parameter
)
{
  string file_name = db_prefix_ + ".nam";
  ifstream in;
  string line;

  names_ = new string[number_sequences_];

  in.open(file_name.c_str());

  uint32_t i;
  for (i = 0; i < number_sequences_ && !in.eof(); ++i) {
    std::getline(in, line);
    names_[i] = line;
  }

  if (i < number_sequences_) {
    cerr << "warning : couldn't read all sequence names" << endl;
  }

  in.close();
}


void DB::SetPositions
(
    void
)
{
  string file_name = db_prefix_ + ".pos";
  ifstream in;

  in.open(file_name.c_str(), ios::binary);

  positions_ = new uint32_t[number_sequences_];
  in.read((char *)positions_, sizeof(uint32_t)*number_sequences_);

  in.close();
}

void DB::SetSequences
(
    // no parameter
)
{
  string file_name = db_prefix_ + ".seq";
  ifstream in;

  in.open(file_name.c_str(), ios::binary);

  sequences_ = new uint8_t[sequences_length_];
  in.read((char *)sequences_, sizeof(uint8_t)*sequences_length_);

  in.close();
}

void DB::SetIndex
(
    // no parameter
)
{
  uint32_t seed;
  uint32_t keys_count_length;
  uint32_t *keys_count;
  uint32_t *positions;
  uint32_t positions_length;

  string file_name = db_prefix_ + ".ind";
  ifstream in;

  in.open(file_name.c_str(), ios::binary);
  if (!in) {
    return; // can't open file
  }
  in.read((char *) &seed, sizeof(uint32_t));
  in.read((char *) &keys_count_length, sizeof(uint32_t));
  in.read((char *) &positions_length, sizeof(uint32_t));
  keys_count = new uint32_t[keys_count_length];
  positions = new uint32_t[positions_length];
  in.read((char *) keys_count, sizeof(uint32_t)*keys_count_length);
  in.read((char *) positions, sizeof(uint32_t)*positions_length);

  in.close();
  index_ = new Index(seed, keys_count_length, positions_length, keys_count, positions);
}
