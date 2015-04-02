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

#ifndef ALIGNER_H_
#define ALIGNER_H_

#include <fstream>
#include <string>
#include <list>
#include "statistics.h"
#include "alignment.h"
#include "command.h"

using namespace std;

class Query;
class DB;
class ScoreMatrix;

typedef struct {
  string output_file_name;
  string query_file_prefix;
  uint32_t start_query_file_id;
  uint32_t end_query_file_id;
  string db_file_prefix;
  uint32_t log_region_size;
  uint32_t shift_size;
  uint32_t threshold;
  uint32_t max_list_length;
  int open_gap;
  int extend_gap;
  ScoreMatrix *score_matrix;
  KarlinParameters statistics_parameters;
  uint32_t extend;
  uint32_t best;
  int device;
  int output_style;
  bool verbose;
} AlignerOption;

class Aligner : public Command{
public:
  //Aligner();
  //virtual ~Aligner();

  int Execute(int argc, char *argv[]);

private:
  uint32_t next_query_id_;
  list<Alignment> next_alignment_list_;

  int SetOption(int argc, char *argv[], AlignerOption &option);

  void SearchNext(Query *query,DB *db, AlignerOption &option, vector<Alignment> &alignment_list);

  void SearchNextCpu(Query *query,DB *db, AlignerOption &option, vector<Alignment> &alignment_list);

  void CalculateScore(Query *query,DB *db, vector<Alignment> &alignment_list, AlignerOption &option);

  void CalculateScoreCpu(Query *query,DB *db, vector<Alignment> &alignment_list, AlignerOption &option);

  void Merge(vector< vector<Alignment> > &result_list, vector<Alignment> &alignment_list, Query *query, DB *db, AlignerOption &option);

  void TraceBack(Query *query,DB *db, Alignment &alignment, AlignerOption &option);

  void WriteOutput(ostream &out, vector< vector<Alignment> > &result_list, Query *query, uint64_t db_length, const KarlinParameters &statistical_parameters);

  void WriteOutputV2(ostream &out, vector< vector<Alignment> > &result_list, Query *query);
  void WriteOutputV1(ostream &out, vector< vector<Alignment> > &result_list, Query *query);
};

#endif /* ALIGNER_H_ */
