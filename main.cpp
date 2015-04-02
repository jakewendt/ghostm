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
#include <string.h>
#include <time.h>
#include "command.h"
#include "query_creator.h"
#include "db_creator.h"
#include "aligner.h"
#include "main.h"

using namespace std;
int usage() {
  cerr << "release\t" << kRelease << endl;
  cerr << "Command and Options" << endl;
  cerr << "db: convert a database fasta file to formatted GHOSTM database files " << endl;
  cerr << endl;
  cerr << "ghostm db [-i dbFastaFile] [-o dbName] [-k kSize] [-l chunkSize]" << endl;
  cerr << endl;
  cerr << "Options:" << endl;
  cerr << "(Required)" << endl;
  cerr << "-i STR    Protein sequences in FASTA format for a database" << endl;
  cerr << "-o STR    The name of database" << endl;
 cerr << "(Optional)" << endl;
  cerr << "-k INT    The size of K-mer\'s K.  [4]" << endl;
  cerr << "-l INT    Chunk size of the database (MB) [128]" << endl;
  cerr << endl;
  cerr << "                                  " << endl;
  cerr << "qry:  convert a query fasta file to formatted input query files" << endl;
  cerr << "" << endl;
  cerr << "ghostm qry [-i qryFastaFile] [-o qryName] [-l chunkSize] [-t queryType]" << endl;
  cerr << endl;
  cerr << "Options:" << endl;
  cerr << "(Required)" << endl;
  cerr << "-i STR    DNA or Protein sequences in FASTA format for query" << endl;
  cerr << "-o STR    The name of query" << endl;
  cerr << "(Optional)" << endl;
  cerr << "-l INT          Max query sequence length [75]" << endl;
  cerr << "-L INT          Chunk size of the query (MB) [128]" << endl;
  cerr << "-t STR    Sequence type of a query fasta file [p]" << endl;
  cerr << "  d   DNA" << endl;
  cerr << "  p   Amino acids" << endl;
  cerr << endl;
  cerr << "  " << endl;
  cerr << "aln:  Search homologs by using formatted queries and a database" << endl;
  cerr << "  " << endl;
  cerr << "ghostm aln [-b best][-D deviceId] [-l CandidatesSize] [-s skipSize] [-t threshold] [-r regionSize] [-e extendSize][-G openGap] [-E extendGap] [-M scoreMatrix] [?i queries] [-d databes] [-o output]" << endl;
  cerr << endl;
  cerr << "Options:" << endl;
  cerr << "(Required)" << endl;
  cerr << "-i STR    Input query file (must be formatted)" << endl;
  cerr << "-d STR    database file (must be formatted)" << endl;
  cerr << "-o STR    Output file" << endl;
 cerr << endl;
 cerr << "-D INT    GPU ID (run without a GPU if this option is not given)" << endl;
 cerr << "" << endl;
 cerr << "(Optional)" << endl;
 cerr << "-v        Verbose mode" << endl;
 cerr << "-M STR    Score matrix file[BLOSUM62]" << endl;
 cerr << "-G INT    Open gap penalty [11]" << endl;
 cerr << "-E INT    Extend gap penalty [1]" << endl;
 cerr << "-b INT    The number of the output for a query [10]" << endl;
 cerr << "-l INT    Maximun size of the candidates (MB) [128]" << endl;
 cerr << "-s INT    Skip number of query\'s K-mer [2]" << endl;
 cerr << "-r INT    The size of the regions [8]" << endl;
 cerr << "-t INT    Required minumum number of candidate seeds in a search region [2]" << endl;
 cerr << "-e INT    The width for extending an alignment region [2]" << endl;
 cerr << "-S INT		Start query chunk id [0]" << endl;	
 cerr << "-E INT		Last query chunk id [last chunk id]" << endl;
 return 1;
}

int main(int argc, char *argv[]) {

  Command *command;
  try {
    //clock_t start = clock();
    if (argc < 2) {
      return usage();
    } else if (strcmp(argv[1], "qry") == 0) {
      command = new QueryCreator();
    } else if (strcmp(argv[1], "db") == 0) {
      command = new DBCreator();
    } else if (strcmp(argv[1], "aln") == 0) {
      command = new Aligner();
    } else {
      cerr << "[main] unrecognized command " << argv[1] << endl;
      return 1;
    }
    if (command->Execute(argc - 1, argv + 1) == SUCCESS) {
      //cerr << endl << "Complete!" << endl;
    }
  } catch (exception &e) {
    cerr << e.what() << endl;
  }
  //cerr << argv[1] << "\t";
  //cerr << (float)(clock() - start) / (float)CLOCKS_PER_SEC << endl;
  return 0;
}
