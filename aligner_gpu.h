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

#ifndef ALIGNER_GPU_H_
#define ALIGNER_GPU_H_

#include <stdint.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C"{
#endif

  int InitGpu ();

  size_t GetNeededGPUMemorySize
  (
   uint32_t seed,
   uint32_t shift_size,
   uint32_t max_list_length,
   uint32_t max_query_length,
   uint32_t max_number_queries,
   uint32_t max_db_length
  );

  int CheckGpuMemory
  (
   uint32_t seed,
   uint32_t shift_size,
   uint32_t max_list_length,
   uint32_t max_query_length,
   uint32_t max_number_queries,
   uint32_t max_db_length
  );

  int SetOptionGpu
  (
   uint32_t max_list_length,
   int score_matrix[],
   int device
   );

  void printGpuInfo(int device);
  
  int SetQueryGpu
  (
   uint8_t sequences[], 
   uint32_t number_sequences,
   uint32_t sequence_length
   );

  int SetDbGpu
  (
   uint8_t sequences[],
   uint32_t sequences_legnth,
   uint32_t keys_count[],
   uint32_t keys_count_length,
   uint32_t positions[],
   uint32_t positions_length
   );

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
   );

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
   );

  int FreeGpu ();
    
  
    
#ifdef __cplusplus
}
#endif

#endif /* ALIGNER_GPU_H_ */













/* JAKE - copied from cutil.h */

#  define CUDA_SAFE_CALL_NO_SYNC( call) {                                    \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }

#  define CUDA_SAFE_CALL( call)     CUDA_SAFE_CALL_NO_SYNC(call);                                            \

    //! Check for CUDA error
#ifdef _DEBUG
#  define CUT_CHECK_ERROR(errorMessage) {                                    \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    err = CUT_DEVICE_SYNCHRONIZE();                                           \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    }
#else
#  define CUT_CHECK_ERROR(errorMessage) {                                    \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    }
#endif


