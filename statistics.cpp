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

#include "assert.h"
#include "statistics.h"
#include "score_matrix.h"
#include "karlin.h"
#include <iostream>
#include <vector>
#include <stdint.h>
#include <math.h>
#include <stdexcept>

using namespace std;

using namespace std;

float Statistics::Nominal2Normalized(int normal_score, const KarlinParameters &paramters) {
  float normalized_score;
  normalized_score = (((static_cast<float>(normal_score)*paramters.lambda) - log(paramters.K))/static_cast<float>(log(2.0)));
  return normalized_score;
}
int Statistics::Normalized2Nominal(float normalized_score, const KarlinParameters &paramters) {
  float normal_score;
  normal_score = (normalized_score*static_cast<float>(log(2.0)) + paramters.K)/paramters.lambda;
  return static_cast<int>(floor(normal_score));
}

double Statistics::Nominal2EValue(int normal_score, uint64_t search_space, const KarlinParameters &paramters) {
  double e_value;
  e_value = search_space*paramters.K*exp(static_cast<double>(-1.0*normal_score*paramters.lambda));
  return e_value;
}

uint64_t Statistics::CalculateSearchSpace(uint32_t query_length, uint64_t db_length) {
  return query_length*db_length;
}

Statistics::Statistics() :
    min_regular_letter_code_(0), max_regular_letter_code_(19), min_code_(0), max_code_(24) {

  min_regular_letter_code_ = 0;
  max_regular_letter_code_ = 19;
  min_code_ = 0;
  max_code_ = 24;
  letter_prob_.resize(max_code_ + 1);
  for (uint32_t i = 0; i <= max_code_; ++i) {
    letter_prob_[i] = 0.0;
  }

  // robinson prob
  letter_prob_[10] = 90.19;
  letter_prob_[0] = 78.05;
  letter_prob_[7] = 73.77;
  letter_prob_[15] = 71.20;
  letter_prob_[19] = 64.41;
  letter_prob_[6] = 62.95;
  letter_prob_[16] = 58.41;
  letter_prob_[11] = 57.44;
  letter_prob_[3] = 53.64;
  letter_prob_[14] = 52.03;
  letter_prob_[9] = 51.42;
  letter_prob_[1] = 51.29;
  letter_prob_[2] = 44.87;
  letter_prob_[5] = 42.64;
  letter_prob_[13] = 38.56;
  letter_prob_[18] = 32.16;
  letter_prob_[12] = 22.43;
  letter_prob_[8] = 21.99;
  letter_prob_[4] = 19.25;
  letter_prob_[17] = 13.30;
  Normalize(&letter_prob_[0], 1.0);
}

Statistics::~Statistics() {
}

void Statistics::CalculateUngappedIdealKarlinParameters(const ScoreMatrix &score_matrix, KarlinParameters *paramters) {

  int highest = score_matrix.GetHighestValue();
  int lowest = score_matrix.GetLowestValue();

  double score_probabilities0[highest - lowest + 1];
  double *score_probabilities = &score_probabilities0[-lowest];

  CalculateScoreProbabilities(score_matrix, &letter_prob_[0], &letter_prob_[0], score_probabilities);

  BlastKarlinBlkCalc(score_probabilities, lowest, highest, &(paramters->lambda), &(paramters->K), &(paramters->H));

}

void Statistics::CalculateUngappedKarlinParameters(uint8_t *sequence,
    uint32_t length, const ScoreMatrix &score_matrix, KarlinParameters *paramters) {

  double query_frequency[max_code_ + 1];

  ComposeSequenceFrequency(sequence, length, query_frequency);
  Normalize(query_frequency, 1.0);

  int highest = score_matrix.GetHighestValue();
  int lowest = score_matrix.GetLowestValue();

  double score_probabilities0[highest - lowest + 1];
  double *score_probabilities = &score_probabilities0[-lowest];
  CalculateScoreProbabilities(score_matrix, query_frequency, &letter_prob_[0], score_probabilities);

  BlastKarlinBlkCalc(score_probabilities, lowest, highest, &(paramters->lambda), &(paramters->K), &(paramters->H));

}


void Statistics::CalculateGappedKarlinParameters(const ScoreMatrix &score_matrix, int open_gap, int extend_gap, KarlinParameters *parameters) {
  if (score_matrix.GetName() == "BLOSUM62" && open_gap == -11 && extend_gap == -1) {
    parameters->lambda = 0.267;
    parameters->K = 0.041;
    parameters->H = 0.14;
  } else if (score_matrix.GetName() == "PAM30" && open_gap == -9 && extend_gap == -1){
    parameters->lambda = 0.294;
    parameters->K = 0.11;
    parameters->H = 0.61;
  } else {
    throw invalid_argument("error: not support score option");
  }
}

void Statistics::CalculateScoreProbabilities(const ScoreMatrix &score_matrix, double *score_frequency1, double *score_frequency2, double *score_probabilities) {
  int highest = score_matrix.GetHighestValue();
  int lowest = score_matrix.GetLowestValue();

  uint32_t number_letters = score_matrix.GetNumberLetters();

  for (int i = lowest; i <= highest; ++i) {
    score_probabilities[i] = 0.0;
  }

  int *matrix = score_matrix.GetMatrix();
  for (uint8_t c1 = min_code_; c1 <= max_code_; ++c1) {
    uint32_t offset = c1*number_letters;
    for (uint8_t c2 = min_code_; c2 <= max_code_; ++c2) {
      int score = matrix[offset + c2];
      score_probabilities[score] += score_frequency1[c1]*score_frequency2[c2];
    }
  }

  double sum_score = 0.0;
  for (int i = lowest; i <= highest; ++i) {
    if (score_probabilities[i] > 0.0) {
      sum_score += score_probabilities[i];
    }
  }
  if (sum_score <= 0.0) {
    throw std::invalid_argument("sum of score frequencies is 0");
  }
  for (int i = lowest; i <= highest; ++i) {
    if (score_probabilities[i] > 0.0) {
      score_probabilities[i] /= sum_score;
    }
  }
}

void Statistics::Normalize(double *frequency, double norm) {
  if (frequency == NULL || norm == 0.0) {
    assert(!"invalid arguments");
  }

  double sum = 0.0;
  for (uint8_t c = min_code_; c <= max_code_; ++c) {
    double p = frequency[c];
    if (p < 0.0) {
      assert(!"invalid frequencies");
    }
    sum += p;
  }

  if (sum <= 0.0) {
    throw std::invalid_argument("sum of frequencies is 0");
  }

  for (uint8_t c = min_code_; c <= max_code_; ++c) {
    frequency[c] /= sum;
    frequency[c] *= norm;
  }
}

void Statistics::ComposeSequenceFrequency(uint8_t *sequence, uint32_t length, double *frequency) {
  if (sequence == NULL || frequency == NULL) {
    assert(!"invalid arguments");
  }

  for (uint8_t c = min_code_; c <= max_code_; ++c) {
    frequency[c] = 0.0;
  }

  for (uint32_t i = 0; i < length; ++i) {
    uint8_t c = sequence[i];
    if (min_regular_letter_code_ <= c && c <= max_regular_letter_code_) {
      ++frequency[c];
    }
  }
  return;
}



