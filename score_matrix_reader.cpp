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
#include <string>
#include <list>
#include <string.h>
#include <stdlib.h>
#include <sstream>
#include "common.h"
#include "score_matrix.h"
#include "score_matrix_reader.h"
#include "sequence.h"

using namespace std;

const string ScoreMatrixReader::default_matrix_name_="BLOSUM62";
const string ScoreMatrixReader::default_matrix_ = 
"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *\nA  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0 -2 -1  0 -4\n R -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3 -1  0 -1 -4\n N -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3  3  0 -1 -4\n D -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3  4  1 -1 -4\n C  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1 -3 -3 -2 -4\n Q -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2  0  3 -1 -4\n E -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n G  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3 -1 -2 -1 -4\n H -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3  0  0 -1 -4\n I -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3 -3 -3 -1 -4\n L -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1 -4 -3 -1 -4\n K -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2  0  1 -1 -4\n M -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1 -3 -1 -1 -4\n F -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1 -3 -3 -1 -4\n P -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2 -2 -1 -2 -4\n S  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2  0  0  0 -4\n T  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0 -1 -1  0 -4\n W -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3 -4 -3 -2 -4\n Y -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1 -3 -2 -1 -4\n V  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4 -3 -2 -1 -4\n B -2 -1  3  4 -3  0  1 -1  0 -3 -4  0 -3 -3 -2  0 -1 -4 -3 -3  4  1 -1 -4\n Z -1  0  0  1 -3  3  4 -2  0 -3 -3  1 -1 -3 -1  0 -1 -3 -2 -2  1  4 -1 -4\n X  0 -1 -1 -1 -2 -1 -1 -1 -1 -1 -1 -1 -1 -1 -2  0  0 -2 -1 -1 -1 -1 -1 -4\n * -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4 -4  1\n";

ScoreMatrix* ScoreMatrixReader::Read(const string &filename) {
  ifstream in;
  string name = filename;
  string::size_type index = name.find_last_of('/');
  if (index != string::npos) {
    name = name.substr(index + 1, name.length() - index);
  }

  in.open(filename.c_str());
  if (in) {
    return Read(in, name);
    in.close();
  } else {
    istringstream default_in(default_matrix_);
    return Read(default_in, default_matrix_name_);
  }
}


list<string> ScoreMatrixReader::Split(string str, string delim) {
  list<string> result;
  string::size_type index;
  index = str.find_first_of(delim);
  while (index != string::npos) {
    if (index > 0) {
      result.push_back(str.substr(0, index));
    }
    str = str.substr(index + 1);
    index = str.find_first_of(delim);
  }
  if (str.length() > 0) {
    result.push_back(str);
  }
  return result;
}

 ScoreMatrix* ScoreMatrixReader::Read(istream &in, const string &name) {
  int number_letters = ALPHABET_SIZE;
  int matrix_size = number_letters * number_letters;
  int *matrix = new int[matrix_size];
  string line;
  int line_number = 0;
  char column_headings[number_letters];
  char row_headings[number_letters];
  for (int i = 0; i < matrix_size; ++i) {
    matrix[i] = 0;
  }

  while(!in.eof()) {
    std::getline(in, line);
    if (line[0] != '\0' && line[0] != '#' && line_number < SEQUENCE_END) {
      list<string> str_list = Split(line, " ");
      int i = 0;
      for (list<string>::iterator it = str_list.begin(); it != str_list.end() && i < SEQUENCE_END; ++it, ++i) {
        if (line_number == 0) {
          column_headings[i] = it->c_str()[0];
        } else if (i == 0) {
          row_headings[line_number - 1] = it->c_str()[0];
        } else {
          int value = atoi(it->c_str());
          matrix[Sequence::kProteinToNumber[(int)row_headings[line_number - 1]]*ALPHABET_SIZE +
                 Sequence::kProteinToNumber[(int)column_headings[i - 1]]] = value;
        }
      }
      ++line_number;
    }
  }

  return new ScoreMatrix(name, matrix, number_letters);
}
