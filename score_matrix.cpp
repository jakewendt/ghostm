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

#include "score_matrix.h"
#include <stdint.h>

ScoreMatrix::ScoreMatrix()
:name_(""), matrix_(0), highest_value_(0), lowest_value_(0), number_letters_(0)
{}

ScoreMatrix::ScoreMatrix(string name, int *matrix, uint32_t number_letters)
:name_(name), matrix_(0), highest_value_(0), lowest_value_(0), number_letters_(number_letters)
{
  matrix_ = matrix;
  for (uint32_t i = 0; i < number_letters; ++i) {
    for (uint32_t j = 0; j < number_letters; ++j) {
      int s = matrix_[number_letters_*i + j];
      if (highest_value_ < s) {
        highest_value_ = s;
      }
      if (lowest_value_ > s) {
        lowest_value_ = s;
      }
    }
  }
}
