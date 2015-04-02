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
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include <stdexcept>
#include "common.h"
#include "sequence.h"
#include "sequence_reader.h"
#include "query_creator.h"

using namespace std;
const uint8_t QueryCreator::kCodonTable[] = {

    'k', // AAA
    'n', // AAC
    'k', // AAG
    'n', // AAT

    't', // ACA
    't', // ACC
    't', // ACG
    't', // ACT

    'r', // AGA
    's', // AGC
    'r', // AGG
    's', // AGT

    'i', // ATA
    'i', // ATC
    'm', // ATG
    'i', // ATT

    'q', // CAA
    'h', // CAC
    'q', // CAG
    'h', // CAT

    'p', // CCA
    'p', // CCC
    'p', // CCG
    'p', // CCT

    'r', // CGA
    'r', // CGC
    'r', // CGG
    'r', // CGT

    'l', // CTA
    'l', // CTC
    'l', // CTG
    'l', // CTT

    'e', // GAA
    'd', // GAC
    'e', // GAG
    'd', // GAT

    'a', // GCA
    'a', // GCC
    'a', // GCG
    'a', // GCT

    'g', // GGA
    'g', // GGC
    'g', // GGG
    'g', // GGT

    'v', // GTA
    'v', // GTC
    'v', // GTG
    'v', // GTT

    '*', // TAA
    'y', // TAC
    '*', // TAG
    'y', // TAT

    's', // TCA
    's', // TCC
    's', // TCG
    's', // TCT

    '*', // TGA
    'c', // TGC
    'w', // TGG
    'c', // TGT

    'l', // TTA
    'f', // TTC
    'l', // TTG
    'f', // TTT
    'x'  // other
};

int QueryCreator::SetParameter
(
    int argc,
    char *argv[],
    string *base_output_prefix_pointer,
    string *input_file_pointer,
    QueryCreatorOption *option_pointer
)
{
  // init
  option_pointer->max_length_concatenated_sequence = (1 << 27);
  option_pointer->length_sequence = 75;
  //option_pointer->max_length_concatenated_sequence = ( 1 << 6 );
  option_pointer->type = PROTEIN;

  int c;
  while ((c = getopt(argc, argv, "i:o:l:t:L:")) >= 0) {
    switch (c) {
    case 'i':
      *input_file_pointer = optarg;
      break;

    case 'o':
      *base_output_prefix_pointer = optarg;
      break;

    case 'l': // max length of one sequence
      option_pointer->length_sequence = atoi(optarg);
      break;

    case 'L': // set max_length_concatenated_sequence is Mbyte.
      option_pointer->max_length_concatenated_sequence =  atoi(optarg)*(1 << 20);
      break;

    case 't':
      if (strcmp(optarg, "d") == 0) {
        option_pointer->type = DNA;
      } else if (strcmp(optarg, "p")  == 0){
        option_pointer->type = PROTEIN;
      } else {
        throw invalid_argument("-t is not support " + string(optarg) + ".");
      }
      break;
    default:
      throw invalid_argument("");
    }
  }
  if (option_pointer->type == DNA) {
    option_pointer->max_length_concatenated_sequence =  option_pointer->max_length_concatenated_sequence/2;
    option_pointer->length_sequence /= 3;
    if (option_pointer->length_sequence >= MAX_QUERY_LENGTH) {
      cerr << "Warring: over upper limit of query length. Max is " << MAX_QUERY_LENGTH*3 << "."<<  endl;
      option_pointer->length_sequence = MAX_QUERY_LENGTH;
    }
  } else {
    if (option_pointer->length_sequence >= MAX_QUERY_LENGTH) {
      cerr << "Warring: over upper limit of query length. Max is "<< MAX_QUERY_LENGTH <<"."<<  endl;
      option_pointer->length_sequence = MAX_QUERY_LENGTH;
    }
  }

  return SUCCESS;
}

vector<Sequence *> *QueryCreator::ReadSequences
(
    string input_file,
    QueryCreatorOption option
)
{
  static Sequence *next_sequence = NULL;
  static SequenceReader *sequence_reader = NULL;
  uint32_t sum_length = 0;
  uint32_t sequnece_length = 0;
  Sequence *sequence;
  vector<Sequence *> *sequences;
  sequences = new vector<Sequence *>();
  if (sequence_reader == NULL) {
    if (option.type == DNA) {
      sequence_reader = CreateSequenceReader(input_file.c_str(), FASTA_DNA);
    } else {
      sequence_reader = CreateSequenceReader(input_file.c_str(), FASTA_PROTEIN);
    }
  }
  if (next_sequence != NULL) {
    sequnece_length = next_sequence->GetSequence().length();
    sum_length = sequnece_length;
    if (sum_length >option.max_length_concatenated_sequence) {
      cerr << "error : too small max length." << endl;
      return NULL;
    }

    sequences->push_back(next_sequence);
    next_sequence = NULL;
  }

  while ((sequence = sequence_reader->Read())) {
    sequnece_length = sequence->GetSequence().length();
    sum_length += sequnece_length;
    if (sum_length >option.max_length_concatenated_sequence) {
      next_sequence = sequence;
      break;
    }
    sequences->push_back(sequence);
  }

  if (sequences->size() == 0) {
    delete sequence_reader;
  }

  cerr << "sum length : " <<sum_length << endl;

  return sequences;
}

vector<Sequence *> *QueryCreator::ConvertDNAToProtein
(
    vector<Sequence *> *dna_sequences
)
{
  vector<Sequence *> *protein_sequences = new vector<Sequence *>();
  uint8_t mask = (1 << 2) - 1;
  uint32_t other_codon = 1 << 6;
  char c;
  stringstream ss;
  uint8_t *buf;
  uint8_t *dna_sequence[2];
  uint32_t dna_sequence_length = dna_sequences->at(0)->GetSequence().length();
  uint32_t max_protein_sequence_length = dna_sequence_length/3;

  for (uint32_t i  = 0; i < 2; ++i) {
    dna_sequence[i] = new uint8_t[dna_sequence_length];
  }

  for (uint32_t i = 0; i < dna_sequences->size(); ++i) {
    buf = dna_sequences->at(i)->ToNumber();
    for (uint32_t j = 0; j < dna_sequence_length; ++j) {
      c = buf[j];
      dna_sequence[0][j] = c;
      if (c > 3) {
        dna_sequence[1][dna_sequence_length - j - 1] = c;
      } else {
        c = (~c) & mask;
        dna_sequence[1][dna_sequence_length - j - 1] = c;
      }
    }
    delete [] buf;

    for (uint32_t j = 0; j < 2; ++j) {
      for (uint32_t offset = 0; offset < 3; ++offset) {
        ss.str("");
        uint32_t length = 0;
        bool stop = false;
        for (uint32_t k = offset + 2; k < dna_sequence_length; k += 3) {
          uint32_t codon = 0;
          for (int l = 2; l >= 0; --l) {
            if (dna_sequence[j][k - l] > 3) {
              codon = other_codon;
              break;
            }
            codon = codon << 2;
            codon = codon | dna_sequence[j][k - l];
          }
          switch (codon) {
          //start
          case 14: // AUG
            stop = false;
            break;

          // stop
          case 48: // TAA
          case 50: // TAG
          case 56: // TGA
            stop = true;
            break;
          }
          if (stop) {
            ss << '*';
          } else {
            ss << QueryCreator::kCodonTable[codon];
          }
          ++length;
        }
        for (uint32_t k = length; k < max_protein_sequence_length; ++k) {
          ss << '*';
        }

        protein_sequences->push_back(new Sequence(dna_sequences->at(i)->GetName(), ss.str(), PROTEIN));

      }
    }
  }
  for (uint32_t i = 0; i < 2; ++i) {
    delete [] dna_sequence[i];
  }

  return protein_sequences;
}

int QueryCreator::WriteInformation
(
    string output_prefix,
    int division,
    uint32_t max_length_sequence,
    uint32_t max_numbe_sequences
)
{
  ofstream out;
  string file_name = output_prefix;
  file_name = file_name + ".inf";
  out.open(file_name.c_str(), ios_base::binary);
  out.write((char *) &division, sizeof(division));
  out.write((char *) &max_length_sequence, sizeof(max_length_sequence));
  out.write((char *) &max_numbe_sequences, sizeof(max_numbe_sequences));
  // for tsubame ////////////////////////////////
  for (uint32_t i = 0; i < 32; ++i) {
    out.write((char *) &division, sizeof(division));
  }
  /////////////////////////////////////////////
  out.close();

  return SUCCESS;
}

int QueryCreator::WriteSubInformation
(
    string output_prefix,
    uint32_t number_sequences,
    uint32_t length_seuqnece
)
{
  ofstream out;
  string file_name = output_prefix;
  file_name = file_name + ".inf";
  out.open(file_name.c_str(), ios_base::binary);
  out.write((char *) &number_sequences, sizeof(uint32_t));
  out.write((char *) &length_seuqnece, sizeof(uint32_t));
  out.close();

  return SUCCESS;
}

int QueryCreator::WriteNames
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

int QueryCreator::WriteSequences
(
    string output_prefix,
    uint32_t length_seuqnece,
    vector<Sequence *> *sequences
)
{
  ofstream out;
  string file_name = output_prefix;
  uint8_t *sequence;
  vector<uint8_t> seq_data(length_seuqnece);

  file_name = file_name + ".seq";
  out.open(file_name.c_str(), ios::binary);
  for (vector<Sequence *>::iterator it = sequences->begin(); it != sequences->end(); ++it) {
    uint32_t seq_len = (*it)->GetSequence().length();
    sequence = (*it)->ToNumber();
    if (seq_len > length_seuqnece) {
      cerr << "warning : the length of sequence is over. " << (*it)->GetName() << endl;
      seq_len = length_seuqnece;
    }
    memcpy(&seq_data[0], sequence, sizeof(uint8_t)*seq_len);
    if (seq_len < length_seuqnece) {
      memset(&seq_data[seq_len], Sequence::kProteinToNumber[(int)'x'], length_seuqnece - seq_len);
    }
    
    //for (uint32_t i = 0; i < length_seuqnece; ++i) {
    //  cout << (int)seq_data[i] << " ";
    //}
    //cout << endl;
    out.write((char *) &seq_data[0], sizeof(uint8_t)*length_seuqnece);
    delete [] sequence;
  }
  out.close();
  return SUCCESS;
}

void QueryCreator::DeleteSequences
(
    vector<Sequence *> *sequences
)
{
  for (vector<Sequence *>::iterator it = sequences->begin(); it != sequences->end(); ++it) {
    delete *it;
  }
  delete sequences;
}



int QueryCreator::Execute
(
    int argc,
    char *argv[]
)
{
  clock_t start;
  stringstream output_prefix;

  string base_output_prefix = "";
  string input_file = "";
  QueryCreatorOption option;
  uint32_t max_number_sequences = 0;
  vector<Sequence *> *sequences;

  if (SetParameter(argc, argv, &base_output_prefix, &input_file, &option) != SUCCESS) {
    return FAILURE;
  }

  cerr << "length_sequence : " <<option.length_sequence << endl;
  cerr << "max_length_concatenated_sequence : " <<option.max_length_concatenated_sequence << endl;

  for (int i = 0;; ++i) {
    // set output file prefix
    output_prefix.str("");
    output_prefix << base_output_prefix << "_" << i;

    // read sequences
    cerr << "[query] read sequences ..." << endl;
    start = clock();
    if ((sequences = ReadSequences(input_file, option)) == NULL) {
      return FAILURE;
    } else if (sequences->size() == 0) {
      cerr << "[query] red all sequences ..."<<  endl;
      cerr << "[query] division number : " << i <<  endl;
      if (WriteInformation(base_output_prefix, i, option.length_sequence, max_number_sequences)  != SUCCESS) {
        return FAILURE;
      }
      break;
    }
    if (option.type == DNA) {
      vector<Sequence *> *dna_sequences = sequences;
      sequences = ConvertDNAToProtein(dna_sequences);
      DeleteSequences(dna_sequences);
    }

    cerr << "number sequences : " << sequences->size() << endl;
    cerr << (float)(clock() - start) / (float)CLOCKS_PER_SEC <<"sec" << endl;
    if (max_number_sequences < sequences->size()) {
      max_number_sequences = sequences->size();
    }
    // write info file
    cerr << "[query] create info file ..." << endl;
    start = clock();
    if (WriteSubInformation(output_prefix.str(), sequences->size(), option.length_sequence)  != SUCCESS) {
      return FAILURE;
    }
    cerr << (float)(clock() - start) / (float)CLOCKS_PER_SEC <<"sec" << endl;

    // write names file
    cerr << "[query] create names file ..." << endl;
    start = clock();
    if (WriteNames(output_prefix.str(), sequences)  != SUCCESS) {
      return FAILURE;
    }
    cerr << (float)(clock() - start) / (float)CLOCKS_PER_SEC <<"sec" << endl;

    // write sequences file
    cerr << "[query] create sequences file ..." << endl;
    start = clock();
    if (WriteSequences(output_prefix.str(), option.length_sequence, sequences)  != SUCCESS) {
      return FAILURE;
    }
    cerr << (float)(clock() - start) / (float)CLOCKS_PER_SEC <<"sec" << endl;

    DeleteSequences(sequences);
  }

  return SUCCESS;
}
