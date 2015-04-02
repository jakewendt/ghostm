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

#ifndef STATISTICS_H_
#define STATISTICS_H_

#include "score_matrix.h"
#include <vector>
#include <stdint.h>

using namespace std;

typedef struct KarlinParameters {
  float lambda;
  float K;
  float H;
} KarlinBlk;

class Statistics {
public:
  static float Nominal2Normalized(int normal_score, const KarlinParameters &paramters);
  static int Normalized2Nominal(float normalized_score,
      const KarlinParameters &paramters);
  static double Nominal2EValue(int normal_score, uint64_t search_space, const KarlinParameters &paramters);

  static uint64_t CalculateSearchSpace(uint32_t query_length, uint64_t db_length);

  Statistics();
  virtual ~Statistics();

  void CalculateUngappedIdealKarlinParameters(
      const ScoreMatrix &score_matrix, KarlinParameters *paramters);
  void CalculateUngappedKarlinParameters(uint8_t *query, uint32_t length,
      const ScoreMatrix &score_matrix, KarlinParameters *paramters);
  void CalculateGappedKarlinParameters(const ScoreMatrix &score_matrix,
      int open_gap, int extend_gap, KarlinParameters *parameters);
  void CalculateAlphaBeta(const ScoreMatrix &score_matrix, int open_gap, int extend_gap, double *alpha, double *beta);

private:
  void CalculateScoreProbabilities(const ScoreMatrix &score_matrix,
      double *score_frequency1, double *score_frequency2,
      double *score_probabilities);
  void Normalize(double *frequency, double norm);
  void ComposeSequenceFrequency(uint8_t *sequence, uint32_t length,
      double *frequency);

  vector<double> letter_prob_;
  uint8_t min_regular_letter_code_;
  uint8_t max_regular_letter_code_;
  uint8_t min_code_;
  uint8_t max_code_;

};

#endif /* STATISTICS_H_ */
