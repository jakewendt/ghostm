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

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <getopt.h>
#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <string.h>
#include <stdint.h>
#include <time.h>
#include <math.h>
#include <stdexcept>
#include <limits.h>
#include "statistics.h"
#include "score_matrix.h"
#include "score_matrix_reader.h"
#include "query.h"
#include "query_reader.h"
#include "db.h"
#include "db_reader.h"
#include "common.h"
#include "alignment.h"
#include "aligner.h"
#include "aligner_gpu.h"

class AlignmentComp
{
public:
  bool operator()
  (
     const Alignment &a1,
     const Alignment &a2
  ) const
  {
    return (a1.GetScore() > a2.GetScore());
  }
};

int Aligner::Execute(int argc, char *argv[]) {
  AlignerOption option;
  SetOption(argc, argv,option);
  ofstream out(option.output_file_name.c_str());

  if(option.verbose){
    cout << "#     G H O S T M "<< endl;
    cout << "# * GPU-base HOmology Search Tool for Metagenomics *"<< endl << endl;
  }

  if (option.device != CPU) {
#ifdef NUMBER_GPUS
    InitGpu();
    {
      DBReader db_reader(option.db_file_prefix);
      QueryReader query_reader(option.query_file_prefix);
      if (CheckGpuMemory(db_reader.GetSeed(),option.shift_size, option.max_list_length, query_reader.GetMaxLengthConcatenatedSequence(), 
			  query_reader.GetMaxNumberSequences(), db_reader.GetMaxLengthConcatenatedSequence())) {
        throw invalid_argument(string("error:Out of GPU memory"));
      }
    }
    if(option.verbose){ 
      cout << "  Maximun size of the candidates (-l option): " << option.max_list_length 
	   << " (" << (option.max_list_length>>20) << "MB)"<< endl;
      cout << "  Init GPU device [" << option.device << "] ... ";
      cout.flush();
    }
    SetOptionGpu(option.max_list_length,  option.score_matrix->GetMatrix(), option.device);
    if(option.verbose){ cout << "ok." << endl; printGpuInfo(option.device); cout << endl;}
#endif
  }

  clock_t start;
  QueryReader query_reader(option.query_file_prefix);
  Query *query = NULL;
  if (option.start_query_file_id == UINT_MAX) {
    query = query_reader.Read();
  } else {
    query = query_reader.Read(option.start_query_file_id);
  }

  if (query != NULL) {
    while (query != NULL) {
      if (option.device != CPU) {
#ifdef NUMBER_GPUS
	SetQueryGpu(query->GetSequences(), query->GetNumberSequences(), query->GetSequenceLength());
#endif
      }

      vector< vector<Alignment> > result_list(query->GetNumberSequences());
      DBReader db_reader(option.db_file_prefix);
      DB *db = db_reader.Read();
      if (db != NULL) {
        while (db != NULL) {

	  if (option.device != CPU) {
#ifdef NUMBER_GPUS
	    Index *db_index = db->GetIndex();
	    SetDbGpu(db->GetSequences(), db->GetSequencesLength(), db_index->GetKeysCount(),
		     db_index->GetKeysCountLength(), db_index->GetAllPositions(), db_index->GetPositionsLength());
#endif
	  }

          next_query_id_ = 0;
          next_alignment_list_.clear();
          vector<Alignment> alignment_list;
	  while (1) {
	    if(option.verbose) cout << "|Search alignment candidates ... ";
	    start = clock();
	    SearchNext(query, db, option, alignment_list);
	    if (alignment_list.empty()) {
	      if(option.verbose) cout << "fin." << endl;
	      break;
	    }
	    if(option.verbose) cout << " " << (float)(clock() - start)/ (float)CLOCKS_PER_SEC << " sec." << endl;
#if 0
	    for (vector<Alignment>::iterator it = alignment_list.begin(); it != alignment_list.end(); ++it) {
	      cout <<"id     " << it->GetQueryId() << " start " << it->GetDbStart() << endl;
	    }
#endif
	    if(option.verbose) cout << "|Calculate scores ... ";
	    start = clock();
	    CalculateScore(query, db, alignment_list, option);
	    if(option.verbose) cout << (float)(clock() - start)/ (float)CLOCKS_PER_SEC << " sec." << endl;

#if 0
	    for (vector<Alignment>::iterator it = alignment_list.begin(); it != alignment_list.end(); ++it) {
	      cout <<"id     " << it->GetQueryId() << " start " << it->GetDbStart() << " end " << it->GetDbEnd()  << " score " << it->GetScore()<< endl;
	    }
#endif
	    
	    if(option.verbose) cout << "|Merge results ... ";
	    start = clock();
	    Merge(result_list, alignment_list, query, db, option);
	    if(option.verbose) cout << (float)(clock() - start)/ (float)CLOCKS_PER_SEC << "sec" << endl;
#if 0
	    cerr << "result_list size " << result_list.size() << endl;
	    for (uint32_t i = 0; i < query->GetNumberSequences(); ++i) {
	      cerr << "result_list[" << i << "] size " << result_list[i].size() << endl;
	      for (vector<Alignment>::iterator it = result_list[i].begin(); it != result_list[i].end(); ++it) {
		cerr <<"id     " << it->GetQueryId() << endl;
		cerr << "start " << it->GetDbStart() << endl;
		cerr << "end " << it->GetDbEnd() << endl;
	      }
	    }
#endif
	  }
	  delete db;
    db = db_reader.Read();
  }

	if(option.verbose) cout << "|Write results ... ";
	start = clock();
	DBReader db_reader(option.db_file_prefix);
	switch (option.output_style) {
	case 1 :
	  WriteOutputV1(out, result_list, query);
	  break;
	case 2 :
    WriteOutputV2(out, result_list, query);
    break;
	default :
	  WriteOutput(out, result_list, query, db_reader.GetSumDbLength(), option.statistics_parameters);
	  break;
	}
	if(option.verbose) cout <<  (float)(clock() - start)/ (float)CLOCKS_PER_SEC << " sec." << endl;

      } else {
        cerr << "[Aligner] error: don't find db file." << endl;
        db_reader.Close();
        delete query;
        break;
      }
      db_reader.Close();
      delete query;
      query = NULL;
      if (query_reader.GetNextId() <= option.end_query_file_id) {
	query = query_reader.Read();
      }
    }
  } else {
    cerr << "[Aligner] error: don't find query file." << endl;
  }

  if (option.device != CPU) {
#ifdef NUMBER_GPUS
    FreeGpu();
#endif
  }

  query_reader.Close();
  if ( option.score_matrix != NULL) {
    delete option.score_matrix;
  }
  out.close();
  if(option.verbose) cout << "Complete." << endl;

  return SUCCESS;
}

int Aligner::SetOption(int argc, char *argv[], AlignerOption &option) {
  // init
  option.output_file_name = "";
  option.query_file_prefix = "";
  option.db_file_prefix = "";
  option.log_region_size = 4;
  option.shift_size = 2;
  option.threshold = 2;
  option.max_list_length = 1 << 27; // 27 for 128MB
  option.open_gap = -11;
  option.extend_gap = -1;
  option.score_matrix = NULL;
  option.extend = 2;
  option.best = 10;
  option.device = CPU;
  option.start_query_file_id = UINT_MAX;
  option.end_query_file_id = UINT_MAX;
  option.output_style = 0;
  option.verbose = false;

  string score_matrix_file_name = "BLOSUM62";

  int c;
  while ((c = getopt(argc, argv, "b:d:D:e:E:G:i:l:M:o:r:s:t:S:L:y:v")) >= 0) {
    switch (c) {
    case 'b':
      option.best = atoi(optarg);
      break;

    case 'd':
      option.db_file_prefix = optarg;
      break;

    case 'D':
#ifdef NUMBER_GPUS
      option.device =  atoi(optarg);
#else
      throw invalid_argument(string("not support gpu."));
#endif
      break;

    case 'e':
      option.extend = atoi(optarg);
      break;

    case 'E':
      option.extend_gap = (-1)*atoi(optarg);
      break;

    case 'G':
      option.open_gap = (-1)*atoi(optarg);
      break;

    case 'i':
      option.query_file_prefix = optarg;
      break;

    case 'S':
      option.start_query_file_id = atoi(optarg);
      break;

    case 'L':
      option.end_query_file_id = atoi(optarg);
      break;

    case 'l':
      option.max_list_length = atoi(optarg)*(1 << 20);
      break;

    case 'M':
      score_matrix_file_name = optarg;
      break;

    case 'o':
      option.output_file_name = optarg;
      break;

    case 'r':
      option.log_region_size = (int)log2(atoi(optarg));
      if (option.log_region_size < 1) {
        option.log_region_size = 1;
      }
      break;

    case 's':
      option.shift_size = atoi(optarg);
      break;

    case 't':
      option.threshold = atoi(optarg);
      break;

    case 'y':
      option.output_style = atoi(optarg);
      break;

    case 'v':
      option.verbose = true;
      break;

    default:
      throw invalid_argument("");
    }
  }
  ScoreMatrixReader reader;
  option.score_matrix = reader.Read(score_matrix_file_name);

  if (option.score_matrix == NULL) {
    throw invalid_argument(string("can't find score matrix. :"));
  }

  if (option.output_style < 0 && 2 <= option.output_style) {
    throw invalid_argument(string("can't support output style"));
  }

  if (option.output_style == 0) {
    Statistics statistics;
    statistics.CalculateGappedKarlinParameters(*option.score_matrix, option.open_gap, option.extend_gap, &option.statistics_parameters);
  }
  return SUCCESS;
}

void Aligner::SearchNext(Query *query,DB *db, AlignerOption &option, vector<Alignment> &alignment_list) {
  if (!alignment_list.empty()) {
    alignment_list.clear();
  }
  if (option.device != CPU) {
#ifdef NUMBER_GPUS
    Index *db_index = db->GetIndex();
    vector<uint32_t> alignment_count_list(query->GetNumberSequences()+1);
    vector<uint32_t> starts(option.max_list_length);
    uint32_t query_count = SearchNextGpu(query->GetSequenceLength(), 
					 query->GetNumberSequences(), 
					 db_index->GetSeed(),
					 option.threshold, 
					 option.shift_size, option.log_region_size, option.max_list_length,
					 next_query_id_, &alignment_count_list[0], &starts[0]);
    if (query_count > 0) {
      //cerr << "  number query: " << query_count << endl;
      uint32_t number_alignments = alignment_count_list[query_count];
      alignment_list.resize(number_alignments);
      for (uint32_t i = 0, j = 1; i < number_alignments; ++i) {
	for (;alignment_count_list[j] <= i; ++j)
	  ;
	alignment_list[i].SetQueryId(j-1+next_query_id_);
	alignment_list[i].SetDbStart(starts[i]);
#if 0
	cerr << "id: " << alignment_list[i].GetQueryId() <<  " start: " << alignment_list[i].GetDbStart() << endl;
#endif
      }
      next_query_id_ += query_count;
    }
#endif
  } else {
    SearchNextCpu(query, db, option, alignment_list);
  }
}

void Aligner::SearchNextCpu(Query *query,DB *db, AlignerOption &option, vector<Alignment> &alignment_list) {
  uint32_t alignment_count = 0;
  if (next_query_id_ == query->GetNumberSequences()) {
      return;
  } else {
    copy(next_alignment_list_.begin(), next_alignment_list_.end(), back_inserter(alignment_list));
    next_alignment_list_.clear();
  }
  alignment_count = alignment_list.size();

  uint8_t *query_sequences = query->GetSequences();
  uint32_t query_sequence_length = query->GetSequenceLength();
  uint32_t number_query_sequences = query->GetNumberSequences();
  Index *db_index = db->GetIndex();
  uint32_t seed_length = db_index->GetSeedLength();

  uint32_t list_length = (query_sequence_length - seed_length)/option.shift_size + 1;
  uint32_t distance_list[list_length];
  uint32_t *positions_list[list_length];
  uint32_t positions_length_list[list_length];
  uint32_t id_list[list_length];

  // debug /////////////////////
  //uint32_t p;
  ////////////////////////////

  uint32_t d;
  uint32_t distance;
  uint32_t count;
  uint32_t next_distance;
  uint32_t next_count;

  uint32_t query_offset;
  uint32_t threshold = option.threshold - 1;

  for (uint32_t i = next_query_id_; i < number_query_sequences; ++i) {

    // init
    query_offset = i*query_sequence_length;
    for (uint32_t j = 0; j < list_length; ++j) {
      positions_list[j] = db_index->GetPositions(&query_sequences[query_offset + j*option.shift_size], &(positions_length_list[j]));

      // debug /////////////////////
      //d = positions_length_list[j];
      ////////////////////////////

      uint32_t k = 0;
      for (; k < positions_length_list[j] && positions_list[j][k] < j*option.shift_size; ++k)
        ;

      id_list[j] = k;
      distance_list[j] = UINT_MAX;
      if (k < positions_length_list[j]) {
        distance_list[j] = (positions_list[j][k] - j*option.shift_size) >> option.log_region_size;
        // debug /////////////////////////
        //p = positions_list[j][k];
        //d = sort_list[j].distance;
        ////////////////////////////
        ++id_list[j];
      }
    }

    // debug //////////////////////////
    //for (uint32_t j = 0; j < list_length; ++j) {
    //  cerr << sort_list[j].id << " " << sort_list[j].distance << endl;
    //}
    ////////////////////////////////////

    distance = 0;
    count = 0;
    while (1) {
      // check min distance
      next_count = 0;
      next_distance = distance_list[0];
      for (uint32_t j = 1; j < list_length; ++j) {
        if (distance_list[j] < next_distance) {
          next_distance = distance_list[j];
        }
      }
      if (next_distance == UINT_MAX) {
        break;
      }

      // count
      for (uint32_t j = 0; j < list_length; ++j) {
        if (next_distance == distance_list[j]) {
          ++next_count;
          distance_list[j] = UINT_MAX;
          for (uint32_t k = id_list[j]; k < positions_length_list[j]; ++k) {
            d = (positions_list[j][k] - j*option.shift_size) >> option.log_region_size;

            if (d != next_distance) {
              distance_list[j] = d;
              id_list[j] = k + 1;
              break;
            }
          }
        }
      }

      if ((next_distance - distance) == 1) {
        count += next_count;
      }

      if (count > threshold) {
        distance = distance << option.log_region_size;
        Alignment a;
        a.SetQueryId(i);
        a.SetDbStart(distance);
        next_alignment_list_.push_back(a);
        ++alignment_count;
      }
      count = next_count;
      distance = next_distance;

    }

    // last check
    if (count > threshold) {
      distance = distance << option.log_region_size;
      Alignment a;
      a.SetQueryId(i);
      a.SetDbStart(distance);
      next_alignment_list_.push_back(a);
      ++alignment_count;
    }


    if (alignment_count > option.max_list_length) {
      next_query_id_ = i + 1;
      return;
    }
    copy(next_alignment_list_.begin(), next_alignment_list_.end(), back_inserter(alignment_list));
    next_alignment_list_.clear();
  }

  next_query_id_ = number_query_sequences;
  return;
}

void Aligner::CalculateScore(Query *query,DB *db, vector<Alignment> &alignment_list, AlignerOption &option) {
  if (option.device != CPU) {
#ifdef NUMBER_GPUS
    vector<uint32_t> scores(option.max_list_length);
    vector<uint32_t> ends(option.max_list_length);;
    uint32_t base_search_length = query->GetSequenceLength() + 2*option.extend + 2*(1 << option.log_region_size);
    CalculateScoreGpu(db->GetSequencesLength(), query->GetSequenceLength(), 
		      alignment_list.size(), &scores[0], &ends[0], 
		      base_search_length, option.extend, option.open_gap, option.extend_gap);

    uint32_t number_alignments = alignment_list.size();
    for (uint32_t i = 0; i < number_alignments; ++i) {
      alignment_list[i].SetScore(scores[i]);
      alignment_list[i].SetDbEnd(ends[i]);

    }
#endif
  } else {
    CalculateScoreCpu(query, db, alignment_list, option);
  }
}

void Aligner::CalculateScoreCpu(Query *query,DB *db, vector<Alignment> &alignment_list, AlignerOption &option) {
  uint8_t *db_sequence = db->GetSequences();
  uint32_t db_length = db->GetSequencesLength();
  uint8_t *query_sequences = query->GetSequences();
  uint32_t query_sequence_length = query->GetSequenceLength();
  uint32_t base_search_length = query->GetSequenceLength() + 2*option.extend + 2*(1 << option.log_region_size);
  uint32_t offset = option.extend;
  int *score_matrix = option.score_matrix->GetMatrix();
  int open_gap = option.open_gap;
  int extend_gap = option.extend_gap;
  uint8_t db_character;
  uint32_t j, k;
  int score = 0;
  int local_score = 0; // score in cell
  int max_score = 0; // max score at the ticket
  uint32_t max_end = 0; // alignment end of max score at the ticket
  int db_offset;
  int query_offset;
  int score_matrix_offset;
  uint32_t db_search_length;
  uint32_t query_search_length;
  int dp_column[MAX_COLUMN_LENGTH];
  int insertion_column[MAX_COLUMN_LENGTH];
  int deletion_score;
  int temp_score; // stored score of previous cell

  query_search_length = query_sequence_length + 1;

  for (vector<Alignment>::iterator it = alignment_list.begin(); it != alignment_list.end(); ++it) {

    // init
    db_offset = it->GetDbStart() - offset;
    if (db_offset < 0) {
      db_offset = 0;
    }
    db_search_length = base_search_length;
    if (db_offset + db_search_length > db_length) {
      db_search_length = db_length - db_offset;
    }

    query_offset = it->GetQueryId()*query_sequence_length - 1;

    max_score = 0;
    for (j = 0; j < query_search_length; ++j) {
      dp_column[j] = 0;
      insertion_column[j] = 0;
    }
#if 0
    // debug ///////////////////////////////////
    cerr << "query_id:" << it->GetQueryId() << " db offset:" << db_offset << endl;
    fprintf(stderr, "      ");
    for (uint32_t x = 1; x < query_search_length; ++x) {
      fprintf(stderr, "%3d", query_sequences[query_offset + x]);
    }
    cerr << endl;

    fprintf(stderr, "   ");
    for (uint32_t x = 0; x < query_search_length; ++x) {
      fprintf(stderr, "%3d", dp_column[x]);
    }
    cerr << endl;
    //////////////////////////////////////////////////////
#endif
    // calculate score
    for (j = 0; j < db_search_length; ++j) {
      if ((db_character = db_sequence[db_offset + j]) != SEQUENCE_END) {
        score_matrix_offset = db_character*ALPHABET_SIZE;
        temp_score = 0;
        deletion_score = 0;
        for (k = 1; k < query_search_length; ++k) {
          local_score = 0;
          // match or mismatch
          score = temp_score + score_matrix[score_matrix_offset + query_sequences[query_offset + k]];
          if (score > 0) {
            local_score = score;
          }

          // insertion
          if (insertion_column[k] + extend_gap < dp_column[k] + open_gap) {
            insertion_column[k] = dp_column[k] + open_gap;
          } else {
            insertion_column[k] += extend_gap;
          }

          if (insertion_column[k] > local_score) {
            local_score = insertion_column[k];
          }

          // deletion
          score = dp_column[k - 1] + open_gap;
          if (deletion_score + extend_gap < dp_column[k - 1] + open_gap) {
            deletion_score = dp_column[k - 1] + open_gap;
          } else {
            deletion_score += extend_gap;
          }

          if (deletion_score > local_score) {
            local_score = deletion_score;
          }

          // update score column
          temp_score = dp_column[k];
          dp_column[k] = local_score;

          // update max score
          if (local_score >= max_score) {
            max_score = local_score;
            max_end = j;
          }
        }
#if 0
      // debug /////////////////////////////////////////////
          fprintf(stderr, "%3d", db_character);
          for (uint32_t x = 0; x < query_search_length; ++x) {
            fprintf(stderr, "%3d", dp_column[x]);
          }
          cerr << endl;
      //////////////////////////////////////////////////////
#endif
      } else {
        // reset score_column
        for (k = 0; k < query_search_length; ++k) {
          dp_column[k] = 0;
          insertion_column[k] = 0;
        }
      }
    }

    // set alignment result
    it->SetScore(max_score);
    it->SetDbEnd(db_offset + max_end);
    // debug /////////////////////
    //if (i == 37047) {
    //cerr << "max score:" << scores[i] << endl;
    //cerr << "start:" << starts[i]<< endl;
    //cerr << "end:" << ends[i] << endl;
    //cerr << endl;
    //}
    ////////////////////////////
  }
}

void Aligner::Merge(vector< vector<Alignment> > &result_list, vector<Alignment> &alignment_list, Query *query, DB *db, AlignerOption &option) {
  vector<Alignment>::iterator alignment_list_it = alignment_list.begin();
  vector<uint32_t> overlap(db->GetNumberSequences());
  vector<Alignment> l;
  for (vector<uint32_t>::iterator it = overlap.begin(); it != overlap.end(); ++it) {
    *it = UINT_MAX;
  }
  // debug ///////////////////////////////////////////////////////
  //cerr << "alignment list size: " << alignment_list.size() << endl;
  ////////////////////////////////////////////////////////////////
  string prev_query_name = query->GetName(0);
  for (uint32_t i = 0; i < query->GetNumberSequences(); ++i) {
    string query_name = query->GetName(i);
    if (prev_query_name != query_name) {
      uint32_t id = i - 1; 
      sort(l.begin(), l.end(), AlignmentComp());
      for (vector<Alignment>::iterator it = l.begin(); it != l.end(); ++it) {
	uint32_t db_id = it->GetDbId();
	if (db_id == UINT_MAX) {
	  db_id = db->GetID(it->GetDbEnd());
	  if (overlap[db_id] != id) {
	    overlap[db_id] = id;
	    TraceBack(query, db, *it, option);
	    uint32_t db_start = it->GetDbStart();
	    uint32_t db_end = it->GetDbEnd();
	    uint32_t db_position = db->GetPosition(db_id);
	    it->SetDbId(db_id);
	    it->SetDbName(db->GetName(db_id));
	    it->SetDbStart(db_start - db_position);
	    it->SetDbEnd(db_end - db_position);
	    result_list[id].push_back(*it);
	  }
	} else {
	  result_list[id].push_back(*it);
	}
	if (result_list[id].size() >= option.best) {
	  break;
	}
      }
      l.clear();
    }
    prev_query_name = query_name;
    vector<Alignment>::iterator start = alignment_list.begin();
    //cerr << "size: " << result_list[i].size() << endl;
    //cerr << "capacity: " <<result_list[i].capacity() << endl;
    for (; alignment_list_it != alignment_list.end(); ++alignment_list_it) {
      if (i != alignment_list_it->GetQueryId()) {
        break;
      }
      l.push_back(*alignment_list_it);
    }
    for (vector<Alignment>::iterator it = result_list[i].begin(); it != result_list[i].end(); ++it) {
      l.push_back(*it);
    }
    result_list[i].clear();
  }

  uint32_t last = query->GetNumberSequences() - 1;
  sort(l.begin(), l.end(), AlignmentComp());
  for (vector<Alignment>::iterator it = l.begin(); it != l.end(); ++it) {
    uint32_t db_id = it->GetDbId();
    if (db_id == UINT_MAX) {
      db_id = db->GetID(it->GetDbEnd());
      if (overlap[db_id] != last) {
	overlap[db_id] = last;
	TraceBack(query, db, *it, option);
	uint32_t db_start = it->GetDbStart();
	uint32_t db_end = it->GetDbEnd();
	uint32_t db_position = db->GetPosition(db_id);
	it->SetDbId(db_id);
	it->SetDbName(db->GetName(db_id));
	it->SetDbStart(db_start - db_position);
	it->SetDbEnd(db_end - db_position);
	result_list[last].push_back(*it);
      }
    } else {
      result_list[last].push_back(*it);
    }
    if (result_list[last].size() >= option.best) {
      break;
    }
  }
}

void Aligner::TraceBack(Query *query,DB *db, Alignment &alignment, AlignerOption &option) {
  uint8_t *db_sequence = db->GetSequences();
  uint8_t *query_sequences = query->GetSequences();
  uint32_t query_sequence_length = query->GetSequenceLength();
  uint32_t base_search_length = query->GetSequenceLength() + 2*option.extend *2*(1 << option.log_region_size);
  int *score_matrix = option.score_matrix->GetMatrix();
  int open_gap = option.open_gap;
  int extend_gap = option.extend_gap;
  uint8_t db_character;
  uint32_t j;
  int k;
  int score = 0;
  int local_score = 0; // score in cell
  int max_score = 0; // max score at the ticket
  uint32_t max_start = 0; // alignment start of max score at the ticket
  uint32_t db_offset;
  uint32_t query_offset;
  int score_matrix_offset;
  uint32_t db_search_length;
  int query_search_length;
  int dp_column[MAX_COLUMN_LENGTH];
  int insertion_column[MAX_COLUMN_LENGTH];
  int deletion_score;
  int temp_score; // stored score of previous cell
  uint32_t temp_match, new_match, ins_match, del_match;
  uint32_t temp_aln_len, new_aln_len, ins_aln_len, del_aln_len;
  uint32_t qry_aln_len[MAX_COLUMN_LENGTH], aln_match[MAX_COLUMN_LENGTH];
  uint32_t max_aln_match, max_qry_start, max_qry_aln_len;

  query_search_length = query->GetSequenceLength() + 1;
  // init
  db_offset = alignment.GetDbEnd();

  db_search_length = base_search_length;
  if (db_offset < db_search_length) {
    db_search_length = db_offset + 1;
  }
  query_offset = alignment.GetQueryId() * query_sequence_length;

  max_score = 0;
  for (k = 0; k < query_search_length; ++k) {
    dp_column[k] = 0;
    insertion_column[k]= 0;
    aln_match[k] = 0;
    qry_aln_len[k] = 0;
  }

  //////////////////////////////////////////////////////
  // calculate score
  /*
  cout << "======  ";
  for (k = query_search_length - 2; 0 <= k; --k) {
    fprintf(stdout, "%3d", query_sequences[query_offset + k]);
  }
  cout << endl;
  cout << "--------------------------------------------------------" << endl;
  */
  max_aln_match = max_qry_start = max_qry_aln_len = 0;
  for (j = 0; j < db_search_length; ++j) {
    if ((db_character = db_sequence[db_offset - j]) != SEQUENCE_END) {
      score_matrix_offset = db_character * ALPHABET_SIZE;
      temp_score = 0;
      deletion_score = 0;
      temp_match =  new_match = ins_match = del_match = 0;
      temp_aln_len =  new_aln_len = ins_aln_len = del_aln_len= 0;
      for (k = query_search_length - 2; 0 <= k; --k) {
        local_score = 0;
	new_match   = 0;
	new_aln_len     = 0;
        // match or mismatch
        score = temp_score + score_matrix[score_matrix_offset
            + query_sequences[query_offset + k]];
        if (score > 0) {
          local_score = score;
	  if(db_character == query_sequences[query_offset + k]){
	    new_match = temp_match + 1;
	  }else{
	    new_match = temp_match;
	  }

	  new_aln_len = temp_aln_len + 1;
	  
        }

        // insertion
        if (insertion_column[k] + extend_gap < dp_column[k] + open_gap) {
          insertion_column[k] = dp_column[k] + open_gap;
        } else {
          insertion_column[k] += extend_gap;
        }
	ins_match = aln_match[k];
	ins_aln_len = qry_aln_len[k] + 1;

        if (insertion_column[k] > local_score) {
          local_score = insertion_column[k];
	  new_match   = ins_match;
	  new_aln_len     = ins_aln_len;
        }

        // deletion
        score = dp_column[k + 1] + open_gap;
        if (deletion_score + extend_gap < dp_column[k + 1] + open_gap) {
          deletion_score = dp_column[k + 1] + open_gap;
        } else {
          deletion_score += extend_gap;
        }
	del_match      = aln_match[k+1];
	del_aln_len        = qry_aln_len[k+1] + 1;

        if (deletion_score > local_score) {
          local_score = deletion_score;
	  new_match   = del_match;
	  new_aln_len     = del_aln_len;
        }

        // update score column
        temp_score = dp_column[k];
        dp_column[k] = local_score;
	
	temp_match = aln_match[k];
	aln_match[k] = new_match;

	temp_aln_len = qry_aln_len[k];
	qry_aln_len[k] = new_aln_len;

        // update max score
        //if (local_score >= max_score) {
	if (local_score > max_score) { //extend if the score increases
          max_score = local_score;
          max_start = j;
	  max_aln_match = new_match;
	  max_qry_start = k;
	  max_qry_aln_len   = new_aln_len;
        }
      }
      // debug /////////////////////////////////////////////
      /*
      fprintf(stdout, "|%3d|", db_character);
      for (int x = query_search_length - 1; 0 <= x; --x) {
        fprintf(stdout, "%3d", dp_column[x]);
      }
      cout << endl;
      
      fprintf(stdout, "|***|");
      for (int x = query_search_length - 1; 0 <= x; --x) {
        fprintf(stdout, "%3d", aln_match[x]);
      }
      cout << endl;
      
      //fprintf(stdout, "|---|");
      //for (int x = query_search_length - 1; 0 <= x; --x) {
      //  fprintf(stdout, "%3d", qry_aln_len[x]);
      //}
      //cout << endl;
      */
      //////////////////////////////////////////////////////
    } else {
      break;
    }
  }
  //get real query seq length
  {
    int aln_len, db_aln_len;
    db_aln_len  = max_start+1;
    //qry_aln_len = max_qry_aln_len;
    //if(db_aln_len > qry_aln_len){ aln_len = qry_aln_len; }
    //else{ aln_len = db_aln_len; }
    aln_len = max_qry_aln_len;

    //printf("#SCORE:%d: MATCH:%d / LEN:%d => SEQID: %g\n", max_score, max_aln_match, aln_len,
    //   (float)max_aln_match/(float)aln_len );
    //printf("#[DB]    START:%d  END:%d: LEN:%d \n", db_offset-max_start, db_offset, db_aln_len);
  
    alignment.SetDbStart(db_offset - max_start);
    alignment.SetSeqId( (float)max_aln_match/(float)aln_len  );
    alignment.SetAlnLen( aln_len );
    alignment.SetAlnMatch( max_aln_match );
  }
}

void Aligner::WriteOutput(ostream &out, vector< vector<Alignment> > &result_list, Query *query, uint64_t db_length, const KarlinParameters &statistical_parameters) {
  uint8_t *query_sequences = query->GetSequences();
  for (uint32_t i = 0; i < query->GetNumberSequences(); ++i) {
    string query_name = query->GetName(i);
    uint32_t query_length = 0;
    uint32_t start = i*query->GetSequenceLength();
    uint32_t end = (i + 1)*query->GetSequenceLength() - 1;
    uint32_t offset = 0;
    for (offset = end; offset > start && query_sequences[offset] == BASE_X; --offset)
      ;
    query_length = offset - start + 1;
    uint64_t search_space = Statistics::CalculateSearchSpace(query_length, db_length);
    //out << "# Fields: Query id, Subject id, % identity, alignment length, matches, s. start, s. end, e-value, bit-score \n";
    for (vector<Alignment>::iterator it = result_list[i].begin(); it != result_list[i].end(); ++it) {
      float bit_score = Statistics::Nominal2Normalized(it->GetScore(), statistical_parameters);
      float e_value = Statistics::Nominal2EValue(it->GetScore(), search_space, statistical_parameters);
      out << query_name << "\t";
      out << it->GetDbName() << "\t";
      out << it->GetSeqId()*100 << "\t";
      out << it->GetAlnLen() << "\t";
      out << it->GetAlnMatch() << "\t";
      out << it->GetDbStart() + 1<< "\t";
      out << it->GetDbEnd() + 1 << "\t";
      out << e_value << "\t";
      out << bit_score << "\t";
      out << endl;
    }
  }
}


void Aligner::WriteOutputV2(ostream &out, vector< vector<Alignment> > &result_list, Query *query) {
  for (uint32_t i = 0; i < query->GetNumberSequences(); ++i) {
    string query_name = query->GetName(i);
    //out << "# Fields: Query id, Subject id, raw score, s. start, s. end, % identity, alignment length, matches\n";
    for (vector<Alignment>::iterator it = result_list[i].begin(); it != result_list[i].end(); ++it) {

      out << query_name << "\t";
      out << it->GetDbName() << "\t";
      out << it->GetScore() << "\t";
      out << it->GetDbStart() + 1<< "\t";
      out << it->GetDbEnd() + 1 << "\t";
      out << it->GetSeqId() << "\t";
      out << it->GetAlnLen() << "\t";
      out << it->GetAlnMatch() ;
      out << endl;
    }
  }
}

void Aligner::WriteOutputV1(ostream &out, vector< vector<Alignment> > &result_list, Query *query) {
  for (uint32_t i = 0; i < query->GetNumberSequences(); ++i) {
    string query_name = query->GetName(i);
    for (vector<Alignment>::iterator it = result_list[i].begin(); it != result_list[i].end(); ++it) {
      out << query_name << "\t";
      out << it->GetDbName() << "\t";
      out << it->GetScore() << "\t";
      out << it->GetDbStart() + 1<< "\t";
      out << it->GetDbEnd() + 1<< endl;
    }
  }
}


