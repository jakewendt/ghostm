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

#ifndef FASTA_SEQUENCE_READER_H_
#define FASTA_SEQUENCE_READER_H_

#include <fstream>
#include <assert.h>
#include <string>
#include <ios>
#include "sequence.h"
#include "sequence_reader.h"

// create sequence from fasta file
class FASTASequenceReader : public SequenceReader
{
private:
  uint8_t sequence_type_;
  std::ifstream in_; // input file stream
  std::string last_line_; // last input line

public:
  // create sequence
  Sequence* Read
  (
      // no parameter
  );

  void ChangeType
  (
      uint8_t type
  )
  {
    sequence_type_ = type;
  }

  // constructor
  FASTASequenceReader
  (
      const char* file_name, // open file name
      const int type
  )
  {
    in_.open(file_name);
    if (!in_) {
      throw std::ios_base::failure(std::string("can't open a file. :") + std::string(file_name));
    }
    last_line_ = "";
    sequence_type_ = type;
  }

  // use default copy constructor
  // FASTASequenceReader(const FASTASequenceReader old_fasta_sequence_reader)

  // use default "="
  // FASTASequenceReader operator = (const FASTASequenceReader& old_fasta_sequence_reader)

  // destructor
  ~FASTASequenceReader
  (
      // no parameter
  )
  {
    in_.close();
  }
};

#endif /* FASTA_SEQUENCE_READER_H_ */
