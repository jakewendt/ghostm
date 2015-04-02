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

#ifndef QUERY_CREATER_H_
#define QUERY_CREATER_H_

#include <vector>
#include <string>
#include <stdint.h>
#include "command.h"

using namespace std;

class Sequence;
class SequenceReader;

typedef struct {
  uint8_t type; // string type
  uint32_t length_sequence;
  uint32_t max_length_concatenated_sequence;
} QueryCreatorOption;

class QueryCreator : public Command{
private:
  static const uint8_t kCodonTable[];

  int SetParameter
  (
      int argc,
      char *argv[],
      string *base_output_prefix_pointer,
      string *input_file_pointer,
      QueryCreatorOption *option_pointer
  );

  vector<Sequence *> *ReadSequences
  (
      string input_file,
      QueryCreatorOption option
  );

  vector<Sequence *> *ConvertDNAToProtein
  (
      vector<Sequence *> *dna_sequences
  );

  int WriteInformation
  (
      string output_prefix,
      int division,
      uint32_t max_length_concatenated_sequence,
      uint32_t max_numbe_sequences
  );

  int WriteSubInformation
  (
      string output_prefix,
      uint32_t number_sequences,
      uint32_t length_seuqnece
  );

  int WriteNames
  (
      string output_prefix,
      vector<Sequence *> *sequences
  );

  int WriteSequences
  (
      string output_prefix,
      uint32_t length_seuqnece,
      vector<Sequence *> *sequences
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
  //QueryCreator()

  // use default destructor
  // virtual ~QueryCreator();
};

#endif /* QUERY_CREATER_H_ */
