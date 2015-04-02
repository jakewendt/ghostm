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

#ifndef DB_CREATOR_H_
#define DB_CREATOR_H_

#include <vector>
#include <string>
#include <stdint.h>
#include <limits.h>
#include "command.h"

using namespace std;

class Sequence;
class SequenceReader;
class Index;

typedef struct {
  uint32_t seed;
  uint32_t max_length_concatenated_sequence;
} DBCreatorOption;

class DBCreator : public Command{
private:

  int SetParameter
  (
      int argc,
      char *argv[],
      string *base_output_prefix_pointer,
      string *input_file_pointer,
      DBCreatorOption *option_pointer
  );


  vector<Sequence *> *ReadSequences
  (
      string input_file,
      DBCreatorOption option
  );

  uint8_t *ConvertSequences
  (
      vector<Sequence *> *sequences,
      uint32_t *sequences_data_length_pointer,
      uint32_t *positions_pointer[]
  );

// construct
  Index *ConstructIndex
  (
      vector<Sequence *> *sequences,
      uint8_t sequences_data[],
      uint32_t sequences_data_length,
      uint32_t positions[],
      DBCreatorOption option
  );

  int WriteInformation
  (
      string output_prefix,
      int division,
      uint32_t seed,
      uint32_t max_length_concatenated_sequence,
      uint64_t sum_length
  );

  int WriteSubInformation
  (
      string output_prefix,
      uint32_t number_sequences,
      uint32_t sequences_data_length
  );

  int WriteNames
  (
      string output_prefix,
      vector<Sequence *> *sequences
  );


  int WriteSequences
  (
      string output_prefix,
      uint8_t sequences_data[],
      uint32_t sequences_data_length
  );

  int WritePositions
  (
      string output_prefix,
      uint32_t positions[],
      uint32_t number_sequences
  );

  int WriteIndex
  (
      string output_prefix,
      Index *index
  );

  void DeleteSequences
  (
      vector<Sequence *> *sequences
  );


public:

  int Execute
  (
      int argc,
      char *argv[]
  );


  // use default constructor
  // DBCreator()

  // use default destructor
  //virtual ~DBCreator();
};

#endif /* DB_CREATOR_H_ */
