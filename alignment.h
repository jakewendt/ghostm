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

#ifndef ALIGNMENT_H_
#define ALIGNMENT_H_

#include <string>
#include <stdint.h>
#include <limits.h>

using namespace std;

class Alignment {
public:
  Alignment()
  :query_id_(UINT_MAX), db_id_(UINT_MAX), db_name_(""), score_(0),
    db_start_(UINT_MAX), db_end_(UINT_MAX),
    aln_len_(UINT_MAX), aln_match_(UINT_MAX), seq_id_(0.0)
  {};
  virtual ~Alignment()
  {};

  uint32_t GetQueryId() const
  {
    return query_id_;
  }

  void SetQueryId(uint32_t query_id)
  {
    query_id_ = query_id;
  }

  uint32_t GetDbId() const
  {
    return db_id_;
  }

  void SetDbId(uint32_t db_id)
  {
    db_id_ = db_id;
  }

  string GetDbName() const
  {
    return db_name_;
  }

  void SetDbName(string db_name)
  {
    db_name_ = db_name;
  }

  uint32_t GetScore() const
  {
    return score_;
  }

  void SetScore( uint32_t score)
  {
    score_ = score;
  }

  uint32_t GetDbStart() const
  {
    return db_start_;
  }

  void SetDbStart( uint32_t db_start)
  {
    db_start_ = db_start;
  }

  uint32_t GetDbEnd() const
  {
    return db_end_;
  }

  void SetDbEnd( uint32_t db_end)
  {
    db_end_ = db_end;
  }

  float GetSeqId() const
  {
    return seq_id_;
  }

  void SetSeqId( float seq_id)
  {
    seq_id_ = seq_id;
  }

  uint32_t GetAlnLen() const
  {
    return aln_len_;
  }

  void SetAlnLen( uint32_t aln_len)
  {
    aln_len_ = aln_len;;
  }
  
  uint32_t GetAlnMatch() const
  {
    return aln_match_;
  }

  void SetAlnMatch( uint32_t aln_match)
  {
    aln_match_ = aln_match;
  }


private:
  uint32_t query_id_;
  uint32_t db_id_;
  string db_name_;
  uint32_t score_;
  uint32_t db_start_;
  uint32_t db_end_;
  uint32_t aln_len_;
  uint32_t aln_match_;
  float seq_id_;
};

#endif /* ALIGNMENT_H_ */
