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

#ifndef SEQUENCE_H_
#define SEQUENCE_H_

#include <string>
#include <assert.h>
#include <stdint.h>

enum sequence_type{
  DNA,
  PROTEIN
};

// sequence data class
class Sequence
{
public:
  static const uint8_t kDNAToNumber[];
  static const uint8_t kProteinToNumber[];
  // return sequence name
  std::string GetName
  (
      // no parameter
  ) const
  {
    return name_;
  }

  // return sequence
  std::string GetSequence
  (
      // no parameter
  ) const
  {
    return sequence_;
  }

  uint8_t GetType
  (
      // no parameter
  )
  {
    return type_;
  }

  // return number array of converted sequence. you need memory free.
  uint8_t *ToNumber
  (
      // no parameter
  );

  // constructor
  Sequence
  (
      const std::string name,    //sequence name
      const std::string sequence  //sequence
  )
  {
    name_ = name;
    sequence_ = sequence;
    type_ = PROTEIN;
  }

  Sequence
  (
      const std::string name,    //sequence name
      const std::string sequence,  //sequence
      const uint8_t type
  )
  {
    name_ = name;
    sequence_ = sequence;
    type_ = type;
  }

  // use default copy constructor
  // Sequence(const Sequence old_sequence)

  // use default "="
  // Sequence operator = (const Sequence& old_sequence)

  // use default destructor
  //virtual ~Sequence();

private:
  uint8_t type_;
  std::string name_;     // sequence name
  std::string sequence_;   // sequence

};

#endif /* SEQUENCE_H_ */
