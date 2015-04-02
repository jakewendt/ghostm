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
#include "query.h"

using namespace std;

void Query::SetNames
(
    // no parameter
)
{
  string file_name = query_prefix_ + ".nam";
  ifstream in;
  string line;

  names_ = new string[number_sequences_];

  in.open(file_name.c_str());
  if (!in) {
    return;
  }
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

void Query::SetSequences
(
    // no parameter
)
{

  string file_name = query_prefix_ + ".seq";
  ifstream in;

  in.open(file_name.c_str(), ios::binary);
  if (!in) {
    // can't open file
    return;
  }
  sequences_ = new uint8_t[number_sequences_*sequence_length_];
  in.read((char *)sequences_, sizeof(uint8_t)*number_sequences_*sequence_length_);
  in.close();
}
