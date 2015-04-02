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

#ifndef COMMON_H_
#define COMMON_H_

#define MAX_FILE_NAME_LENGTH 1024

#define ALPHABET_SIZE 32 // DNA is 4. Amino is 32(20 + X + SequenceEnd).
#define CHARACTER_SIZE 5 // one alphabet bit size

#define SEQUENCE_END 25
#define BASE_X 23

#define NUMBER_GPUS 4
#define MAX_COLUMN_LENGTH 128 // score column length at calculating score
#define MAX_LIST_SIZE MAX_COLUMN_LENGTH // key list size at creating tickets
#define MAX_QUERY_LENGTH (MAX_COLUMN_LENGTH - 1)

enum device {
  CPU = -1
};

#endif /* COMMON_H_ */
