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
#include <vector>
#include <string>
#include <sstream>
#include <getopt.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <string.h>
#include <stdlib.h>
#include <stdexcept>
#include "sequence.h"
#include "sequence_reader.h"
#include "common.h"
#include "index.h"
#include "db_creator.h"

using namespace std;

int DBCreator::SetParameter
(
    int argc,
    char *argv[],
    string *base_output_prefix_pointer,
    string *input_file_pointer,
    DBCreatorOption *option_pointer
)
{
  // init
  option_pointer->seed = ((1 << 4) - 1); // 1111
  option_pointer->max_length_concatenated_sequence = (1 << 27);// 128M
  //option_pointer->max_length_concatenated_sequence = ( 1 << 6 );

  int c;
  while ((c = getopt(argc, argv, "i:o:k:l:")) >= 0) {
    switch (c) {
    case 'i':
      *input_file_pointer = optarg;
      break;

    case 'o':
      *base_output_prefix_pointer = optarg;
      break;

    case 'k': // set seed weight
      option_pointer->seed = ((1 << atoi(optarg)) - 1);
      break;

    case 'l': // set max length. max_length_concatenated_sequence is Mbyte.
      option_pointer->max_length_concatenated_sequence =  atoi(optarg)*(1 << 20);
      break;

    default:
      throw invalid_argument("");
    }
  }

  return SUCCESS;
}

vector<Sequence *> *DBCreator::ReadSequences
(
    string input_file,
    DBCreatorOption option
)
{
  static Sequence *next_sequence = NULL;
  static SequenceReader *sequence_reader = NULL;
  uint32_t sum_length = 0;
  Sequence *sequence;
  vector<Sequence *> *sequences;
  sequences = new vector<Sequence *>();
  if (sequence_reader == NULL) {
    sequence_reader = CreateSequenceReader(input_file.c_str(), FASTA_PROTEIN);
  }
  if (next_sequence != NULL) {
    sum_length = next_sequence->GetSequence().length() + 1; // + sequence end
    if (sum_length > option.max_length_concatenated_sequence) {
      cerr << "error : too small max length." << endl;
      return NULL;
    }
    sequences->push_back(next_sequence);
    next_sequence = NULL;
  }

  while ((sequence = sequence_reader->Read())) {
    sum_length += sequence->GetSequence().length() + 1; // + sequence end
    if (sum_length > option.max_length_concatenated_sequence) {
      next_sequence = sequence;
      sum_length -= sequence->GetSequence().length() + 1;
      break;
    }
    sequences->push_back(sequence);
  }

  if (sequences->size() == 0) {
    delete sequence_reader;
  }
  cerr << "sum length : " << sum_length << endl;
  return sequences;
}

uint8_t *DBCreator::ConvertSequences
(
    vector<Sequence *> *sequences,
    uint32_t *sequences_data_length_pointer,
    uint32_t *positions_pointer[]

)
{
  uint8_t end = SEQUENCE_END;
  uint32_t length = 0;
  uint32_t *positions = new uint32_t[sequences->size()];

  // calculate sequences_data_length and set positions
  uint32_t sequences_data_length = 0;
  for (uint32_t i = 0; i < sequences->size(); ++i) {
    positions[i] = sequences_data_length;
    sequences_data_length += (sequences->at(i)->GetSequence().length() + 1);
  }

  // set sequences_data_
  uint8_t *sequences_data = new uint8_t[sequences_data_length];
  memset(sequences_data, BASE_X, sizeof(uint8_t)*sequences_data_length);
  uint8_t *buf;
  for (uint32_t i = 0; i < sequences->size(); ++i) {
    buf = sequences->at(i)->ToNumber();
    length = sequences->at(i)->GetSequence().length();
    memcpy(sequences_data + positions[i], buf, sizeof(uint8_t)*length);
    memcpy(sequences_data+ positions[i] + length, &end, sizeof(uint8_t));
    delete [] buf;
  }

  *sequences_data_length_pointer = sequences_data_length;
  *positions_pointer = positions;

  return sequences_data;
}

Index *DBCreator::ConstructIndex
(
    vector<Sequence *> *sequences,
    uint8_t sequences_data[],
    uint32_t sequences_data_length,
    uint32_t positions[],
    DBCreatorOption option
)
{
  Index base_index = Index(option.seed, 0, 0, NULL, NULL);

  uint32_t seed_length  = base_index.GetSeedLength();
  uint32_t seed_weight = base_index.GetSeedWeight();

  uint32_t keys_length = sequences_data_length;
  uint32_t *keys = new uint32_t[keys_length];
  for (uint32_t i = 0; i < keys_length; ++i) {
    keys[i] = UINT_MAX;
  }

  uint32_t keys_count_length = (uint32_t)pow((double) ALPHABET_SIZE, (double) seed_weight) + 1;
  uint32_t *keys_count = new uint32_t[keys_count_length];

  for (uint32_t i = 0; i < keys_count_length; ++i) {
    keys_count[i] = 0;
  }

  // set keys
  uint32_t key;
  for (uint32_t i = 0; i < sequences->size(); ++i) {
    if (sequences->at(i)->GetSequence().length() > seed_length) {
      for (uint32_t j = positions[i]; sequences_data[j + seed_length - 1] != SEQUENCE_END; ++j) {

	// check sequence character
	bool contain_x = false;
	for (uint32_t k = 0; k < seed_length; ++k) {
	  if (sequences_data[j + k] == BASE_X) {
	    contain_x = true;
	  }
	}

	if (!contain_x) {
	  key = base_index.GetKey(sequences_data + j);
	  keys[j] = key;
	  ++keys_count[key + 1];
	}
      }
    }
  }

  for (uint32_t i = 1; i < keys_count_length; ++i) {
    keys_count[i] = keys_count[i - 1] + keys_count[i];
  }

  uint32_t *ids = new uint32_t[keys_count[keys_count_length -1]];
  uint32_t *counts = new uint32_t[keys_count_length];
  memset(counts, 0, sizeof(uint32_t)*keys_count_length);

  for (uint32_t i = 0; i < sequences->size(); ++i) {
    for (uint32_t j = positions[i]; sequences_data[j] != SEQUENCE_END; ++j) {
      key = keys[j];
      if (key != UINT_MAX) {
        ids[keys_count[key] + counts[key]] = j;
        ++counts[key];
      }
    }
  }

  Index *index = new Index(option.seed, keys_count_length, keys_count[keys_count_length -1], keys_count, ids);

  delete [] keys;
  delete [] counts;

  return index;
}
int DBCreator::WriteInformation
(
    string output_prefix,
    int division,
    uint32_t seed,
    uint32_t max_length_concatenated_sequence,
    uint64_t sum_length
)
{
  ofstream out;
  string file_name = output_prefix;
  file_name = file_name + ".inf";
  out.open(file_name.c_str(), ios::binary);
  out.write((char *) &division, sizeof(division));
  out.write((char *) &seed, sizeof(seed));
  out.write((char *) &max_length_concatenated_sequence, sizeof(max_length_concatenated_sequence));
  out.write((char *) &sum_length, sizeof(sum_length));

  // for tsubame ////////////////////////////////
  for (uint32_t i = 0; i < 32; ++i) {
    out.write((char *) &division, sizeof(int));
  }
  /////////////////////////////////////////////
  return SUCCESS;
}

int DBCreator::WriteSubInformation
(
    string output_prefix,
    uint32_t number_sequences,
    uint32_t sequences_data_length
)
{
  ofstream out;
  string file_name = output_prefix;
  file_name = file_name + ".inf";
  out.open(file_name.c_str(), ios::binary);
  out.write((char *) &number_sequences, sizeof(uint32_t));
  out.write((char *) &sequences_data_length, sizeof(uint32_t));
  out.close();
  return SUCCESS;
}

int DBCreator::WriteNames
(
    string output_prefix,
    vector<Sequence *> *sequences
)
{
  ofstream out;
  string file_name = output_prefix;
  file_name = file_name + ".nam";
  out.open(file_name.c_str());
  for (vector<Sequence *>::iterator it = sequences->begin(); it != sequences->end(); ++it) {
    out << (*it)->GetName() << endl;
  }
  out.close();

  return SUCCESS;
}

int DBCreator::WriteSequences
(
    string output_prefix,
    uint8_t sequences_data[],
    uint32_t sequences_data_length
)
{
  ofstream out;
  string file_name = output_prefix;
  file_name = file_name + ".seq";
  out.open(file_name.c_str(), ios::binary);
  out.write((char *) sequences_data, sizeof(uint8_t)*sequences_data_length);
  out.close();
  return SUCCESS;
}

int DBCreator::WritePositions
(
    string output_prefix,
    uint32_t positions[],
    uint32_t number_sequences
)
{
  ofstream out;
  string file_name = output_prefix;
  file_name = file_name + ".pos";
  out.open(file_name.c_str(), ios::binary);
  out.write((char *) positions, sizeof(uint32_t)*number_sequences);
  out.close();
  return SUCCESS;
}

int DBCreator::WriteIndex
(
    string output_prefix,
    Index *index
)
{
  ofstream out;
  string file_name = output_prefix;
  file_name = file_name + ".ind";
  int seed = index->GetSeed();
  int keys_count_length = index->GetKeysCountLength();
  int positions_length = index->GetPositionsLength();
  out.open(file_name.c_str(), ios::binary);
  out.write((char *) &seed, sizeof(uint32_t));
  out.write((char *) &keys_count_length, sizeof(uint32_t));
  out.write((char *) &positions_length, sizeof(uint32_t));
  out.write((char *) index->GetKeysCount(), sizeof(uint32_t)*keys_count_length);
  out.write((char *) index->GetAllPositions(), sizeof(uint32_t)*positions_length);
  out.close();
  return SUCCESS;
}


void DBCreator::DeleteSequences
(
    vector<Sequence *> *sequences
)
{
  for (vector<Sequence *>::iterator it = sequences->begin(); it != sequences->end(); ++it) {
    delete *it;
  }
  delete sequences;
}

int DBCreator::Execute
(
    int argc,
    char *argv[]
)
{
  clock_t start;
  stringstream output_prefix;
  string base_output_prefix = "";
  string input_file = "";
  DBCreatorOption option;

  vector<Sequence *> *sequences;
  uint8_t *sequences_data;
  uint32_t sequences_data_length;
  uint32_t *positions;
  uint64_t sum_length = 0;
  Index *index;

  if (SetParameter(argc, argv, &base_output_prefix, &input_file, &option) != SUCCESS) {
    return FAILURE;
  }

  // debug /////////////////////////
  cerr << "seed : "<<option.seed << endl;
  cerr << "max length :" << option.max_length_concatenated_sequence << endl;
  ////////////////////////////////

  for (int i = 0;; ++i) {
    // set output file prefix
    output_prefix.str("");
    output_prefix << base_output_prefix << "_" << i;

    // read sequences
    cerr << "[db] read sequences ..." << endl;
    start = clock();
    if ((sequences = ReadSequences(input_file, option)) == NULL) {
      return FAILURE;
    } else if (sequences->size() == 0) {
      cerr << "[db] red all sequences ..."<<  endl;
      cerr << "[db] division number : " << i <<  endl;
      if (WriteInformation(base_output_prefix, i, option.seed, option.max_length_concatenated_sequence, sum_length)  != SUCCESS) {
        return FAILURE;
      }
      break;
    }
    cerr << "number sequences : " << sequences->size() << endl;
    cerr << (float)(clock() - start) / (float)CLOCKS_PER_SEC <<"sec" << endl;

    // convert sequences
    cerr << "[db] convert sequences ..." << endl;
    start = clock();
    if ((sequences_data = ConvertSequences(sequences, &sequences_data_length, &positions)) == NULL) {
      return FAILURE;
    }
    cerr << (float)(clock() - start) / (float)CLOCKS_PER_SEC <<"sec" << endl;
    sum_length += sequences_data_length - sequences->size();

    // write info file
    cerr << "[db] create info file ..." << endl;
    start = clock();
    if (WriteSubInformation(output_prefix.str(), sequences->size(), sequences_data_length)  != SUCCESS) {
      return FAILURE;
    }
    cerr << (float)(clock() - start) / (float)CLOCKS_PER_SEC <<"sec" << endl;

    // write names file
    cerr << "[db] create names file ..." << endl;
    start = clock();
    if (WriteNames(output_prefix.str(), sequences)  != SUCCESS) {
      return FAILURE;
    }
    cerr << (float)(clock() - start) / (float)CLOCKS_PER_SEC <<"sec" << endl;

    // write sequences
    cerr << "[db] create sequences file ..." << endl;
    start = clock();
    if (WriteSequences(output_prefix.str(), sequences_data, sequences_data_length)  != SUCCESS) {
      return FAILURE;
    }
    cerr << (float)(clock() - start) / (float)CLOCKS_PER_SEC <<"sec" << endl;

    // write names file
    cerr << "[db] create positions file ..." << endl;
    start = clock();
    if (WritePositions(output_prefix.str(), positions, sequences->size())  != SUCCESS) {
      return FAILURE;
    }
    cerr << (float)(clock() - start) / (float)CLOCKS_PER_SEC <<"sec" << endl;

    // construct index
    cerr << "[db] create index file ..." << endl;
    start = clock();
    if ((index = ConstructIndex(sequences, sequences_data, sequences_data_length, positions, option))  == NULL) {
      return FAILURE;
    }
    // write index file
    if (WriteIndex(output_prefix.str(), index)  != SUCCESS) {
      return FAILURE;
    }
    cerr << (float)(clock() - start) / (float)CLOCKS_PER_SEC <<"sec" << endl;

    delete [] sequences_data;
    delete [] positions;
    DeleteSequences(sequences);
    delete index;

  }

  return SUCCESS;
}
