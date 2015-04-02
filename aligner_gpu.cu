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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <cutil.h>
#include <stdint.h>
#include "common.h"
#include "aligner_gpu.h"

#define NUMBER_CONVERT_BLOCKS 128
#define NUMBER_CONVERT_THREADS 256

#define NUMBER_COUNT_BLOCKS 128
#define NUMBER_COUNT_THREADS 256

#define NUMBER_SET_BLOCKS 128
#define NUMBER_SET_THREADS 256

#define NUMBER_SCORE_BLOCKS 128
#define NUMBER_SCORE_THREADS 256

texture<int,  1, cudaReadModeElementType> score_matrix_texture;
texture<uint8_t, 1, cudaReadModeElementType> query_sequences_texture;

int g_device;
int *g_d_score_matrix;
uint32_t g_score_matrix_size;
uint8_t *g_d_query_sequences;
uint32_t g_query_sequences_size;
uint8_t *g_d_db_sequences;
uint32_t g_db_sequences_size;
uint32_t *g_d_keys_count;
uint32_t g_keys_count_size;
uint32_t *g_d_positions;
uint32_t g_positions_size;
uint32_t *g_d_keys;
uint32_t g_keys_size;
uint32_t *g_d_alignment_count_list;
uint32_t *g_alignment_count_list;
uint32_t g_alignment_count_list_size;
uint32_t g_start_query_id;
uint32_t g_query_count;
uint32_t *g_d_scores;
uint32_t g_scores_size;
uint32_t *g_d_starts;
uint32_t g_starts_size;
uint32_t *g_d_ends;
uint32_t g_ends_size;

__global__ void ConvertToKey
(
    uint8_t query_sequences[],
    uint32_t query_sequence_length,
    uint32_t number_query_sequences,
    uint32_t list_length,
    uint32_t seed,
    uint32_t seed_length,
    uint32_t shift_size,
    uint32_t keys[],
    uint32_t keys_length
)
{
  uint32_t i,j;
  uint32_t s;
  uint32_t key;
  uint32_t query_id;
  uint32_t query_offset;
  uint32_t list_id;
  uint32_t stride = gridDim.x*blockDim.x;;


  for (i = blockIdx.x*blockDim.x + threadIdx.x; i < keys_length; i += stride) {
    query_id = i/list_length;
    list_id = i%list_length;
    query_offset = (query_id*query_sequence_length) + (list_id*shift_size);
    for (j = 0, s = seed, key = 0; s != 0; ++j, s >>= 1) {
      if (s & 1) {
        key = key << CHARACTER_SIZE;
        key = key | query_sequences[query_offset + j];
      }
    }
    keys[i] = key;
  }

  return;
}

__global__ void CountQueryAlignment
(
    uint32_t number_query_sequences,
    uint32_t list_length,
    uint32_t keys[],
    uint32_t keys_count[],
    uint32_t positions[],
    uint32_t threshold,
    uint32_t shift_size,
    uint32_t log_region_size,
    uint32_t alignment_count_list[]
)
{
  uint32_t i,j,k,l;
  uint32_t d;
  uint32_t distance;
  uint32_t next_distance;
  uint32_t count;
  uint32_t next_count;
  uint32_t distance_list[MAX_LIST_SIZE];
  uint32_t positions_end_list[MAX_LIST_SIZE];
  uint32_t positions_id_list[MAX_LIST_SIZE];
  uint32_t count_list[MAX_LIST_SIZE];
  uint32_t sequence_offset;
  uint32_t number_alignment_list;
  uint32_t key;
  uint32_t keys_offset;
  uint32_t stride = gridDim.x*blockDim.x;


  // init
  for (i = 0; i < list_length; ++i) {
    distance_list[i] = UINT_MAX;
    positions_id_list[i] = 0;
    positions_end_list[i] = 0;
    count_list[i] = 0;
  }
  --threshold;

  for (i = blockIdx.x*blockDim.x + threadIdx.x; i < number_query_sequences; i += stride) {
    
    // init
    number_alignment_list = 0;
    keys_offset = i*list_length;
    for (j = 0; j < list_length; ++j) {
      sequence_offset = j*shift_size;
      key = keys[keys_offset + j];
      positions_id_list[j] = keys_count[key];
      positions_end_list[j] = keys_count[key + 1];

      for (k = positions_id_list[j]; k < positions_end_list[j] && positions[k] < sequence_offset; ++k)
        ;

      distance_list[j] = UINT_MAX;
      if (k < positions_end_list[j]) {
        distance_list[j] = (positions[k] - sequence_offset) >> log_region_size;
        ++k;
      }
      positions_id_list[j] = k;
      
    }

    distance = 0;
    count = 0;

    while (1) {
      // check min distance
      next_count = 1;
      next_distance = distance_list[0];
      count_list[0] = 0;
      for (j = 1; j < list_length; ++j) {
        if (distance_list[j] < next_distance) {
          next_distance = distance_list[j];
          next_count = 1;
          count_list[0] = j;
        } else if (distance_list[j] == next_distance) {
          count_list[next_count] = j;
          ++next_count;
        }
      }
      
      if (next_distance == UINT_MAX) {
        break;
      }
      
      // set distance list
      for (j = 0; j < next_count; ++j) {
        k = count_list[j];
        sequence_offset = k*shift_size;
        d = UINT_MAX;
        for (l = positions_id_list[k]; l < positions_end_list[k]; ++l) {
          d = (positions[l] - sequence_offset) >> log_region_size;
          if (d != next_distance) {
            ++l;
            break;
          }
	  d = UINT_MAX;
        }
	distance_list[k] = d;
	positions_id_list[k] = l;
      }

      if ((next_distance - distance) == 1) {
        count += next_count;
      }

      if (count > threshold) {
        ++number_alignment_list;
      }

      count = next_count;
      distance = next_distance;
      
    }
    
    // last check
    if (count > threshold) {
      ++number_alignment_list;
    }
    
    alignment_count_list[i] = number_alignment_list;
    
  }
  return;
}

__global__ void SetAlignmentList
(
    uint32_t start_query_id,
    uint32_t number_set_querys,
    uint32_t alignment_list[],
    uint32_t alignment_count_list[],
    uint32_t list_length,
    uint32_t keys[],
    uint32_t keys_count[],
    uint32_t positions[],
    uint32_t threshold,
    uint32_t shift_size,
    uint32_t log_region_size
)
{
  uint32_t i,j,k,l;
  uint32_t d;
  uint32_t distance;
  uint32_t next_distance;
  uint32_t count;
  uint32_t next_count;
  uint32_t distance_list[MAX_LIST_SIZE];
  uint32_t positions_end_list[MAX_LIST_SIZE];
  uint32_t positions_id_list[MAX_LIST_SIZE];
  uint32_t count_list[MAX_LIST_SIZE];
  uint32_t sequence_offset;
  uint32_t alignment_list_id;
  uint32_t alignment_list_end_id;
  uint32_t key;
  uint32_t keys_offset;
  uint32_t stride = gridDim.x*blockDim.x;


  // init
  for (i = 0; i < list_length; ++i) {
    distance_list[i] = UINT_MAX;
    positions_id_list[i] = 0;
    positions_end_list[i] = 0;
    count_list[i] = 0;
  }
  --threshold;

  for (i = blockIdx.x*blockDim.x + threadIdx.x; i < number_set_querys; i += stride) {
    
    // init
    alignment_list_id = alignment_count_list[i];
    alignment_list_end_id = alignment_count_list[i + 1];
    keys_offset = (start_query_id + i)*list_length;

    for (j = 0; j < list_length; ++j) {
      sequence_offset = j*shift_size;
      key = keys[keys_offset + j];
      positions_id_list[j] = keys_count[key];
      positions_end_list[j] = keys_count[key + 1];
      
      for (k = positions_id_list[j]; k < positions_end_list[j] && positions[k] < sequence_offset; ++k)
        ;
      
      distance_list[j] = UINT_MAX;
      if (k < positions_end_list[j]) {
        distance_list[j] = (positions[k] - sequence_offset) >> log_region_size;
        ++k;
      }
      positions_id_list[j] = k;
      
    }

    distance = 0;
    count = 0;
    while (alignment_list_id < alignment_list_end_id) {
      // check min distance
      next_count = 1;
      next_distance = distance_list[0];
      count_list[0] = 0;
      for (j = 1; j < list_length; ++j) {
        if (distance_list[j] < next_distance) {
          next_distance = distance_list[j];
          next_count = 1;
          count_list[0] = j;
        } else if (distance_list[j] == next_distance) {
          count_list[next_count] = j;
          ++next_count;
        }
      }
     
      // set distance list
      for (j = 0; j < next_count; ++j) {
        k = count_list[j];
        sequence_offset = k*shift_size;
        d = UINT_MAX;
        for (l = positions_id_list[k]; l < positions_end_list[k]; ++l) {
          d = (positions[l] - sequence_offset) >> log_region_size;
          if (d != next_distance) {
            ++l;
            break;
          }
	  d = UINT_MAX;
        }
	distance_list[k] = d;
	positions_id_list[k] = l;
      }
    

      if ((next_distance - distance) == 1) {
        count += next_count;
      }

      if (count > threshold) {
        alignment_list[alignment_list_id] = distance << log_region_size;
	++alignment_list_id;
      }

      count = next_count;
      distance = next_distance;
    }
  }
  return;
}

__global__ void CalculateScore
(
 uint8_t db_sequence[],
 uint32_t db_length, // db sequence length
 uint32_t query_sequence_length, // one query length
 uint32_t alignment_count_list[],
 uint32_t alignment_count_list_length,
 uint32_t number_alignment_list, // the number of alignment_list
 uint32_t scores[], // score array in alignment_list
 uint32_t starts[], // start array in alignment_list
 uint32_t ends[], // end array in alignment_list
 uint32_t start_query_id,
 uint32_t base_search_length,
 uint32_t offset,
 int open_gap,
 int extend_gap
)
{
  uint8_t db_character;
  uint32_t i, j, k;
  int score = 0;
  int local_score = 0; // score in cell
  int max_score = 0; // max score at the alignment
  uint32_t max_end = 0; // alignment end of max score at the alignment
  uint32_t alignment_count_id;
  int db_offset;
  int query_offset;
  int score_matrix_offset;
  uint32_t db_search_length;
  uint32_t query_search_length;
  int dp_column[MAX_COLUMN_LENGTH];
  int insertion_column[MAX_COLUMN_LENGTH];
  int deletion_score;
  int temp_score; // stored score of previous cell
  uint32_t stride = gridDim.x*blockDim.x;
  
  alignment_count_id = 0;
  alignment_count_list_length = alignment_count_list_length - 1;
  query_search_length = query_sequence_length + 1;

  
  for (i = blockIdx.x*blockDim.x + threadIdx.x; i< number_alignment_list; i += stride) {

    // init
    db_offset = starts[i] - offset;
    if (db_offset < 0) {
      db_offset = 0;
    }
    db_search_length = base_search_length;
    if (db_offset + db_search_length > db_length) {
      db_search_length = db_length - db_offset;
    }

    for (; alignment_count_id < alignment_count_list_length; ++alignment_count_id) {
     if (alignment_count_list[alignment_count_id]  <= i && i < alignment_count_list[alignment_count_id + 1]) {
       break;
     }
    }
    query_offset = (start_query_id + alignment_count_id)*query_sequence_length - 1;

    max_score = 0;
    for (j = 0; j < query_search_length; ++j) {
      dp_column[j] = 0;
      insertion_column[j] = 0;
    }


    // calculate score
    for (j = 0; j < db_search_length; ++j) {
      if ((db_character = db_sequence[db_offset + j]) != SEQUENCE_END) {
        score_matrix_offset = db_character*ALPHABET_SIZE;
        temp_score = 0;
        deletion_score = 0;
        for (k = 1; k < query_search_length; ++k) {
          local_score = 0;
          // match or mismatch
	  score = temp_score + tex1Dfetch(score_matrix_texture,
					  score_matrix_offset +
					  tex1Dfetch(query_sequences_texture, query_offset + k));
          if (score > local_score) {
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

      } else {
        // reset score_column
        for (k = 0; k < query_search_length; ++k) {
          dp_column[k] = 0;
          insertion_column[k] = 0;
        }
      }

    }

    // set alignment result
    scores[i] = (uint32_t)max_score;
    ends[i] = db_offset + max_end;
    //ends[i] = query_offset;


  }
  return;
}

size_t GetNeededGPUMemorySize
(
  uint32_t seed,
  uint32_t shift_size,
  uint32_t max_list_length,
  uint32_t max_query_length,
  uint32_t max_number_queries,
  uint32_t max_db_length
)
{
  uint32_t seed_weight = 0;
  for (uint32_t s = seed; s > 0; s = s >> 1) {
  	if (s&1) {
  		++seed_weight;
  	}
  }
  uint32_t seed_length = 0;
  for (uint32_t s = seed; s != 0; s >>= 1, ++seed_length)
    ;

  // candidates
  size_t scores_size = sizeof(uint32_t)*max_list_length;
  size_t starts_size = sizeof(uint32_t)*max_list_length;
  size_t ends_size = sizeof(uint32_t)*max_list_length;
  
  // query
  size_t query_sequences_size = sizeof(uint8_t)*max_query_length;
  size_t alignment_count_list_size = sizeof(uint32_t)*(max_number_queries + 1);
  uint32_t keys_length = max_number_queries*((max_query_length - seed_length)/shift_size + 1);
  size_t key_size = sizeof(uint32_t)*keys_length;
  
  // db
  size_t db_sequence_size = sizeof(uint8_t)*max_db_length;
  uint32_t keys_count_length = (uint32_t)pow((double) ALPHABET_SIZE, (double) seed_weight) + 1;
  size_t keys_count_size = sizeof(uint32_t)*keys_count_length;
  size_t positions_size = sizeof(uint32_t)*max_db_length;
  
  // other
  size_t score_matrix_size = sizeof(int)*ALPHABET_SIZE*ALPHABET_SIZE;
  
  return scores_size +
  	starts_size +
  	ends_size +
  	query_sequences_size +
  	alignment_count_list_size +
  	key_size +
  	db_sequence_size +
  	keys_count_size +
  	positions_size +
	score_matrix_size;	
  	
}

// if global memory is over, return 1. otherwise return 0.
int CheckGpuMemory
(
   uint32_t seed,
   uint32_t shift_size,
   uint32_t max_list_length,
   uint32_t max_query_length,
   uint32_t max_number_queries,
   uint32_t max_db_length
)
{
  cudaDeviceProp deviceProp;
  CUDA_SAFE_CALL(cudaGetDeviceProperties(&deviceProp, g_device));
  size_t global_memory_size =  GetNeededGPUMemorySize(seed, shift_size, max_list_length, max_query_length, max_number_queries, max_db_length);

#if 0
    fprintf(stderr, "memory size.       %llu bytes\n", global_memory_size);
    fprintf(stderr, "global memory size.%llu bytes\n", deviceProp.totalGlobalMem);
    fflush(stderr);
#endif


  if (global_memory_size > deviceProp.totalGlobalMem) {
    return 1;
  } else {
    return 0;
  }
}

int InitGpu ()
{
  g_device = 0;
  g_d_score_matrix = NULL;
  g_score_matrix_size = 0;
  g_d_query_sequences = NULL;
  g_query_sequences_size = 0;
  g_d_db_sequences = NULL;
  g_db_sequences_size = 0;
  g_d_keys_count = NULL;
  g_keys_count_size = 0;
  g_d_positions = NULL;
  g_positions_size = 0;
  g_d_keys = NULL;
  g_keys_size = 0;
  g_d_alignment_count_list = NULL;
  g_alignment_count_list = NULL;
  g_alignment_count_list_size = 0;
  g_query_count = 0;
  g_d_scores = NULL;
  g_scores_size = 0;
  g_d_starts = NULL;
  g_starts_size = 0;
  g_d_ends = NULL;
  g_ends_size = 0;

  return 0;
}

int SetOptionGpu
(
 uint32_t max_list_length,
 int score_matrix[],
 int device
)
{
  g_device = device;
  CUDA_SAFE_CALL(cudaSetDevice(g_device));

  // for > cuda 3.0 
  //cudaFuncSetCacheConfig(ConvertToKey, cudaFuncCachePreferL1);
  //cudaFuncSetCacheConfig(CountQueryAlignment, cudaFuncCachePreferL1);
  //cudaFuncSetCacheConfig(SetAlignmentList, cudaFuncCachePreferL1);
  //cudaFuncSetCacheConfig(CalculateScore, cudaFuncCachePreferL1);

  g_scores_size = sizeof(uint32_t)*max_list_length;
  g_starts_size = sizeof(uint32_t)*max_list_length;
  g_ends_size = sizeof(uint32_t)*max_list_length;

  CUDA_SAFE_CALL(cudaMalloc((void **)&g_d_scores, g_scores_size));
  CUDA_SAFE_CALL(cudaMalloc((void **)&g_d_starts, g_starts_size));
  CUDA_SAFE_CALL(cudaMalloc((void **)&g_d_ends, g_ends_size));

  g_score_matrix_size = sizeof(int)*ALPHABET_SIZE*ALPHABET_SIZE;
  
  // set score matrix to device
  CUDA_SAFE_CALL(cudaMalloc((void **)&g_d_score_matrix, g_score_matrix_size));
  CUDA_SAFE_CALL(cudaMemcpy(g_d_score_matrix, score_matrix,
			    g_score_matrix_size, cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaBindTexture(NULL, score_matrix_texture, g_d_score_matrix, g_score_matrix_size));
  return 0;
}

void printGpuInfo(int device){
  cudaDeviceProp deviceProp;
  cudaGetDeviceProperties(&deviceProp, device);
  fprintf(stdout, "  [GPU] device: \"%s\"\n", deviceProp.name);
  fprintf(stdout, "  [GPU] global memory size: %u bytes (%gMB)\n", deviceProp.totalGlobalMem, (deviceProp.totalGlobalMem/1048576.0) );
  //fprintf(stdout, "  [GPU] number of cores: %d\n", nGpuArchCoresPerSM[deviceProp.major] * deviceProp.multiProcessorCount);
}

int SetQueryGpu
(
 uint8_t sequences[], 
 uint32_t number_sequences,
 uint32_t sequence_length
 )
{
  uint32_t new_size;
  new_size = sizeof(uint8_t)*number_sequences*sequence_length;
  if (new_size > g_query_sequences_size) {
    if (g_d_query_sequences != NULL) {
      CUDA_SAFE_CALL(cudaFree(g_d_query_sequences));
    }
    CUDA_SAFE_CALL(cudaMalloc((void **)&g_d_query_sequences, new_size));
     g_query_sequences_size = new_size;
  }
  CUDA_SAFE_CALL(cudaMemcpy(g_d_query_sequences, sequences,
			    new_size, cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaBindTexture(NULL, query_sequences_texture, g_d_query_sequences, 
				 new_size));

  new_size = sizeof(uint32_t)*(number_sequences + 1);
  if (new_size > g_alignment_count_list_size) {
    if (g_d_alignment_count_list != NULL) {
      CUDA_SAFE_CALL(cudaFree(g_d_alignment_count_list));
    }
    if (g_alignment_count_list != NULL) {
      free(g_alignment_count_list);
    }
    CUDA_SAFE_CALL(cudaMalloc((void **)&g_d_alignment_count_list, new_size));
    g_alignment_count_list= (uint32_t *)malloc(sizeof(uint32_t)*(number_sequences + 1));
    if (g_alignment_count_list == NULL) {
      return 1;
    }
    g_alignment_count_list_size = new_size;
  }
  
  return 0;
}

int SetDbGpu
(
 uint8_t sequences[],
 uint32_t sequences_legnth,
 uint32_t keys_count[],
 uint32_t keys_count_length,
 uint32_t positions[],
 uint32_t positions_length
 )
{
  uint32_t new_size = sizeof(uint8_t)*sequences_legnth;
  if (new_size > g_db_sequences_size) {
    if (g_d_db_sequences != NULL) {
      CUDA_SAFE_CALL(cudaFree(g_d_db_sequences));
    }
    CUDA_SAFE_CALL(cudaMalloc((void **)&g_d_db_sequences, new_size));
    g_db_sequences_size = new_size;

  }
  CUDA_SAFE_CALL(cudaMemcpy(g_d_db_sequences, sequences, new_size, cudaMemcpyHostToDevice));
  
  new_size = sizeof(uint32_t)*keys_count_length;
  if (new_size > g_keys_count_size) {
    if (g_d_keys_count != NULL) {
      CUDA_SAFE_CALL(cudaFree(g_d_keys_count));
    }
    CUDA_SAFE_CALL(cudaMalloc((void **)&g_d_keys_count, new_size));
    g_keys_count_size = new_size;
  }
  CUDA_SAFE_CALL(cudaMemcpy(g_d_keys_count, keys_count, new_size, cudaMemcpyHostToDevice));
  
  new_size = sizeof(uint32_t)*positions_length;
  if (new_size > g_positions_size) {
    if (g_d_positions != NULL) {
      CUDA_SAFE_CALL(cudaFree(g_d_positions));
    }
    CUDA_SAFE_CALL(cudaMalloc((void **)&g_d_positions, new_size));
    g_positions_size = new_size;
  }
  CUDA_SAFE_CALL(cudaMemcpy(g_d_positions, positions, new_size, cudaMemcpyHostToDevice));

#if 0
    // debug ///////////////////////////////////////////////////
    fprintf(stderr, "keys_count\n");
    for (uint32_t i = 0; i < 10; ++i) {
      fprintf(stderr, "%u ", keys_count[i]);
    }
    fprintf(stderr, "\n");

    fprintf(stderr, "positions\n");
    for (uint32_t i = 0; i < 10; ++i) {
      fprintf(stderr, "%u ", positions[i]);
    }
    fprintf(stderr, "\n");
    /////////////////////////////////////////////////////////////
#endif
  
  return 0;
}

uint32_t SearchNextGpu
(
    uint32_t query_sequence_length,
    uint32_t number_query_sequences,
    uint32_t seed,
    uint32_t threshold,
    uint32_t shift_size,
    uint32_t log_region_size,
    uint32_t max_number_alignments,
    uint32_t start_query_id,
    uint32_t *alignment_count_list,
    uint32_t *starts
)
{
  uint32_t i,j;
  uint32_t s;
  uint32_t seed_length;
  uint32_t list_length;
  uint32_t keys_length;
  uint32_t query_count;
  uint32_t number_alignments;
  uint32_t new_size;

  // debug ////////////////////////////
  unsigned int timer;
  //float gpu_time;
  ///////////////////////////////////////

  if (start_query_id == number_query_sequences) {
    return 0;
  }
  g_start_query_id = start_query_id;

  //debug ////////////////////////////
  cutCreateTimer(&timer);
  ///////////////////////////////////////
  

  // set seed length
  for (seed_length = 0, s = seed; s != 0; s >>= 1, ++seed_length)
    ;
  
  list_length = (query_sequence_length - seed_length)/shift_size + 1;
  keys_length = number_query_sequences*list_length;


  new_size = sizeof(uint32_t)*keys_length;
  if (new_size > g_keys_size) {
    if (g_d_query_sequences != NULL) {
      CUDA_SAFE_CALL(cudaFree(g_d_keys));
    }
    CUDA_SAFE_CALL(cudaMalloc((void **)&g_d_keys, new_size));
    g_keys_size = new_size;
  }

  if (start_query_id == 0) { 
  //debug /////////////////////////
    cutResetTimer(timer);
    cutStartTimer(timer);
  ////////////////////////////////////


    ConvertToKey<<<NUMBER_CONVERT_BLOCKS, NUMBER_CONVERT_THREADS>>>(g_d_query_sequences, query_sequence_length, number_query_sequences, list_length, seed, seed_length, shift_size, g_d_keys, keys_length);
    
    // debug ////////////////////////////////////
    //cudaThreadSynchronize();
    //cutStopTimer(timer);
    //gpu_time = cutGetTimerValue(timer)*1.0e-03;
    //fprintf(stderr, "\n   convert into keys:  %9.3e [GPU]\n", gpu_time);
    /////////////////////////////////////////////
    
    
#if 0
    // debug ///////////////////////////////////////////////////
    fprintf(stderr, "keys\n");
    uint32_t *keys = (uint32_t *)malloc(sizeof(uint32_t)*keys_length);
    CUDA_SAFE_CALL(cudaMemcpy(keys, g_d_keys,
			      sizeof(uint32_t)*keys_length, cudaMemcpyDeviceToHost));
    for (i = 0; i < list_length; ++i) {
      fprintf(stderr, "%u ", keys[i]);
    }
    fprintf(stderr, "\n");
    free(keys);
    /////////////////////////////////////////////////////////////
#endif

#if 0
    // debug ///////////////////////////////////////////////////
    fprintf(stderr, "keys_count\n");
    uint32_t *keys_count = (uint32_t *)malloc(g_keys_count_size);
    CUDA_SAFE_CALL(cudaMemcpy(keys_count, g_d_keys_count,
			      g_keys_count_size, cudaMemcpyDeviceToHost));
    for (i = 0; i < 10; ++i) {
      fprintf(stderr, "%u ", keys_count[i]);
    }
    fprintf(stderr, "\n");
    free(keys_count);

    fprintf(stderr, "positions\n");
    uint32_t *positions = (uint32_t *)malloc(g_positions_size);
    CUDA_SAFE_CALL(cudaMemcpy(positions, g_d_positions,
			      g_positions_size, cudaMemcpyDeviceToHost));
    for (i = 0; i < 10; ++i) {
      fprintf(stderr, "%u ", positions[i]);
    }
    fprintf(stderr, "\n");
    free(positions);
    /////////////////////////////////////////////////////////////
#endif
    
    // count each query alignments
    alignment_count_list[0] = 0;
    //debug /////////////////////////
    cutResetTimer(timer);
    cutStartTimer(timer);
    /////////////////////////////////////
#if 0
    //debug /////////////////////////
    fprintf(stderr, "number_query_sequences %u \n", number_query_sequences);
    fprintf(stderr, "list_length %u \n", list_length);
    fprintf(stderr, "threshold %u \n", threshold);
    fprintf(stderr, "shift_size %u \n", shift_size);
    fprintf(stderr, "log_region_size %u \n", log_region_size);
    /////////////////////////////////////
#endif
    CountQueryAlignment<<<NUMBER_COUNT_BLOCKS, NUMBER_COUNT_THREADS>>>(number_query_sequences, list_length, 
								       g_d_keys, g_d_keys_count, g_d_positions,
								       threshold, shift_size, log_region_size,
								       g_d_alignment_count_list);

    // debug //////////////////
    //cudaThreadSynchronize();
    //cutStopTimer(timer);
    //gpu_time = (cutGetTimerValue(timer)*1.0e-03);
    //fprintf(stderr, "   count tickets: %9.3e [GPU]\n", gpu_time);
    ////////////////////////////
    g_alignment_count_list[0] = 0;
    CUDA_SAFE_CALL(cudaMemcpy(&g_alignment_count_list[1], g_d_alignment_count_list,
			      sizeof(uint32_t)*number_query_sequences, cudaMemcpyDeviceToHost));
    
#if 0
    // debug ///////////////////////////////////////////////////
    fprintf(stderr, "alignment count\n");
    for (i = 1; i <= /*number_query_sequences*/ 10; ++i) {
      fprintf(stderr, "%u ", alignment_count_list[i]);
      }
    fprintf(stderr, "\n");
    //return 0;
    /////////////////////////////////////////////////////////////
#endif
  }

  number_alignments = 0;
  alignment_count_list[0] = 0;
  for (i = start_query_id + 1, j = 1; 
       (number_alignments + g_alignment_count_list[i] < max_number_alignments) && 
	 (i <= number_query_sequences);
       ++i, ++j) {
    number_alignments += g_alignment_count_list[i];
    alignment_count_list[j] = number_alignments;
  }
  query_count = j - 1;

#if 0
  // debug ///////////////////////////////////////////////////
  fprintf(stderr, "start query id %u \n", start_query_id);
  fprintf(stderr, "query_count %u \n", query_count);
  fprintf(stderr, "alignment_count_list\n");
  for (i = 0; i < query_count + 1; ++i) {
    fprintf(stderr, "%u ", alignment_count_list[i]);
  }
  fprintf(stderr, "\n");
  /////////////////////////////////////////////////////////////
#endif

  CUDA_SAFE_CALL(cudaMemcpy(g_d_alignment_count_list, alignment_count_list, 
			    sizeof(uint32_t)*(query_count + 1), cudaMemcpyHostToDevice));
    

    //debug /////////////////////////
  cutResetTimer(timer);
  cutStartTimer(timer);
  ////////////////////////////////////
  
  SetAlignmentList<<<NUMBER_SET_BLOCKS, NUMBER_SET_THREADS>>>(start_query_id, query_count, g_d_starts, 
							      g_d_alignment_count_list, list_length, g_d_keys,
							      g_d_keys_count, g_d_positions, threshold, 
							      shift_size, log_region_size);

  // debug ////////////////////////////////////
  //cudaThreadSynchronize();
  //cutStopTimer(timer);
  //gpu_time = cutGetTimerValue(timer)*1.0e-03;
  //printf("   set tickets: %9.3e [GPU]\n", gpu_time);
  /////////////////////////////////////////////
  
  CUDA_SAFE_CALL(cudaMemcpy(starts, g_d_starts, sizeof(uint32_t)*number_alignments, cudaMemcpyDeviceToHost));
  
#if 0
  // debug ///////////////////////////////////////////////////
  fprintf(stderr, "starts\n");
  for (i = number_alignments - 10; i < number_alignments; ++i) {
  fprintf(stderr, "%u ", starts[i]);
  }
  fprintf(stderr, "\n");
  /////////////////////////////////////////////////////////////
#endif

  g_query_count = query_count;
  return query_count;
}

void  CalculateScoreGpu
(
 uint32_t db_length, // db sequence length
 uint32_t query_sequence_length, // one query length
 uint32_t number_alignment_list, // the number of alignment_list
 uint32_t scores[], // score array in alignment_list
 uint32_t ends[], // end array in alignment_list
 uint32_t base_search_length,
 uint32_t offset,
 int open_gap,
 int extend_gap
)
{
  dim3 dim_grid(NUMBER_SCORE_BLOCKS);
  dim3 dim_block(NUMBER_SCORE_THREADS);
#if 0
    // debug ///////////////////////////////////////////////////
    fprintf(stderr, "alignment_count_list_length\n");
    uint32_t alignment_count_list_length = g_query_count + 1;
    uint32_t *alignment_count_list = (uint32_t *)malloc(sizeof(uint32_t)*alignment_count_list_length);
    CUDA_SAFE_CALL(cudaMemcpy(alignment_count_list, g_d_alignment_count_list,
			      sizeof(uint32_t)*alignment_count_list_length, cudaMemcpyDeviceToHost));
    for (uint32_t i = 0; i < alignment_count_list_length; ++i) {
      fprintf(stderr, "%u ", alignment_count_list[i]);
    }
    fprintf(stderr, "\n");
    free(alignment_count_list);
    /////////////////////////////////////////////////////////////
#endif

#if 0
    //debug /////////////////////////
    fprintf(stderr, "query_sequence_length %u \n", query_sequence_length);
    fprintf(stderr, "alignment_count_list_length %u \n", g_query_count + 1);
    fprintf(stderr, "number_alignment_list %u \n", number_alignment_list);
    fprintf(stderr, "start_query_id %u \n", g_start_query_id);
    fprintf(stderr, "base_search_length %u \n", base_search_length);
    fprintf(stderr, "offset %u \n", offset);
    fprintf(stderr, "open_gap %d \n", open_gap);
    fprintf(stderr, "extend_gap %d \n", extend_gap);
    /////////////////////////////////////
#endif

  CalculateScore<<<dim_grid, dim_block>>>(g_d_db_sequences, db_length, query_sequence_length, 
					  g_d_alignment_count_list,
					  g_query_count + 1, number_alignment_list, 
					  g_d_scores, g_d_starts, g_d_ends, g_start_query_id, base_search_length, 
					  offset, open_gap, extend_gap);
  CUT_CHECK_ERROR("calculatealignment_listcore() execution failed.\n");

  // return score
  CUDA_SAFE_CALL(cudaMemcpy(scores, g_d_scores, g_scores_size, cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(ends, g_d_ends, g_ends_size, cudaMemcpyDeviceToHost));


  return;
}

int FreeGpu
(
 // no parameter
)
{
  // unbind texture
  CUDA_SAFE_CALL(cudaUnbindTexture(score_matrix_texture));
  CUDA_SAFE_CALL(cudaUnbindTexture(query_sequences_texture));

  // free gpu memory
  CUDA_SAFE_CALL(cudaFree(g_d_db_sequences));
  CUDA_SAFE_CALL(cudaFree(g_d_query_sequences));
  CUDA_SAFE_CALL(cudaFree(g_d_score_matrix));
  CUDA_SAFE_CALL(cudaFree(g_d_alignment_count_list));
  CUDA_SAFE_CALL(cudaFree(g_d_scores));
  CUDA_SAFE_CALL(cudaFree(g_d_starts));
  CUDA_SAFE_CALL(cudaFree(g_d_ends));
 
  // free memory
  free(g_alignment_count_list);
  return 0;
}


