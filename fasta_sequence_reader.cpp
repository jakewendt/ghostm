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

#include <fstream>
#include <string>

#include "sequence.h"
#include "fasta_sequence_reader.h"

// create sequence
// return Sequence pointer, or if file is eof, return false.
Sequence* FASTASequenceReader::Read
(
    // no parameter
)
{
  std::string name;     // got name from file
  std::string sequence = "";  // got sequence form file
  std::string line;
  size_t position;


  if (!in_.eof()) {
    line = last_line_;

    // then jump to the next header line
    while (!in_.eof() && (line.length() == 0 || line.at(0) != '>')) {
      std::getline(in_, line);
    }
    if (in_.eof()) return NULL;

    // set name
    if (line.at(line.length() - 1) == '\r') {
      line = line.substr(0, line.length() - 1);
    }
    position = line.find_first_not_of("> ");
    if (position != std::string::npos) {
      name = line.substr(position);
    }

    // set sequence
    while (!in_.eof()) {
      std::getline(in_, line);
      if (line.length() != 0) {
        if (line.at(0) != '>') {
          if (line.at(line.length() - 1) == '\r') {
            line = line.substr(0, line.length() - 1);
          }
          if (line.at(line.length() - 1) == '+') {
            line = line.substr(0, line.length() - 1);
          }
          sequence += line;
        } else {
          break;
        }
      }
    }

    last_line_ = line;
    return new Sequence(name, sequence, sequence_type_);
  } else {
    return NULL;
  }
}
