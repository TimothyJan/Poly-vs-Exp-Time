
///////////////////////////////////////////////////////////////////////////////
// poly_exp.hpp
//
// Definitions for two algorithms that solve the Maximum Subarray Problem,
// and one algorithm that solves the Subset Sum Problem.
//
///////////////////////////////////////////////////////////////////////////////

#pragma once

#include <algorithm>
#include <functional>
#include <optional>
#include <vector>

//Tim Added
#include <numeric>  //for std::accumulate
#include <utility>  //for std::pair
#include <cmath>    //for std::pow exponent
#include <iostream>

namespace subarray {

// A summed_span represents a non-empty range of indices inside of a vector of
// ints, stored in a begin iterator and end iterator. The class also stores
// the sum of the ints in that range.
//
// Just like in the rest of the C++ standard library, the range includes all
// elements in [begin, end), or in other words the range includes begin, and all
// elements up to BUT NOT INCLUDING end itself.
class summed_span {
public:
  using iterator = std::vector<int>::const_iterator;

private:
  iterator begin_, end_;
  int sum_;

public:

  // Constructor, given the begin iterator, end iterator, and sum of elements
  // in the range. begin must come before end. sum must be the total of all
  // elements in the range. O(1) time.
  summed_span(iterator begin, iterator end, int sum)
  : begin_(begin), end_(end), sum_(sum) {
    assert(begin < end);
  }

  // Constructor, given only the begin and end iterators. The sum is computed
  // in O(n) time.
  summed_span(iterator begin, iterator end)
  : summed_span(begin, end, std::accumulate(begin, end, 0)) {}

  // Equality tests, two spans are equal when each of their iterators are equal.
  bool operator== (const summed_span& rhs) const {
    return (begin_ == rhs.begin_) && (end_ == rhs.end_);
  }

  // Accessors.
  const iterator& begin() const { return begin_; }
  const iterator& end  () const { return end_  ; }
  int sum() const { return sum_; }

  // Compute the number of elements in the span.
  size_t size() const { return end_ - begin_; }

  // Stream insertion operator, so this class is printable.
  friend std::ostream& operator<<(std::ostream& stream, const summed_span& rhs) {
    stream << "summed_span, size=" << rhs.size() << ", sum=" << rhs.sum();
    return stream;
  }
};

// Compute the maximum subarray of input; i.e. the non-empty contiguous span of
// elements with the maximum sum. input must be nonempty. This function uses an
// exhaustive search algorithm that takes O(n^3) time.
summed_span max_subarray_exh(const std::vector<int>& input) {
  assert(!input.empty());
  // TODO: 
  //Base Case
  if(input.size() == 1)
    return summed_span(input.begin(), input.end());

  //General
  int b = 0, e = 1;
  for(int i=0;i<=input.size()-1;i++){
    for(int j=i+1;j<=input.size();j++){
      if(accumulate(input.begin()+i,input.begin()+j,0) > accumulate(input.begin()+b,input.begin()+e,0))
        {
          b=i;
          e=j;
        }
    }
  }
  //  return summed_span(input.begin(), input.begin() + 1);
  return summed_span(input.begin()+b, input.begin()+e);
}

// Compute the maximum subarray using a decrease-by-half algorithm that takes
// O(n log n) time.
std::pair<int,int> maximum_subarray_crossing(const std::vector<int>& input, int low, int middle, int high){
  int left_sum = -2147483647, right_sum = -2147483647, sum = 0, b=0, e=1;
  for(int i = middle; i >= low; i--){
    sum += input[i];
    if(sum > left_sum){
      left_sum = sum;
      b=i;
    }
  }
  sum=0;
  for(int i = middle+1; i <= high; i++){
    sum += input[i];
    if(sum > right_sum){
      right_sum = sum;
      e=i;
    }
  }
  return std::make_pair(b,e+1);
}

std::pair<int,int> maximum_subarray_recursive(const std::vector<int>& input, int low, int high){
  if(low == high)
    return std::make_pair(low, low+1);
  int middle = (low+high)/2;
  std::pair<int,int> entirely_left = maximum_subarray_recursive(input, low, middle);
  std::pair<int,int> entirely_right = maximum_subarray_recursive(input, middle+1, high);
  std::pair<int,int> crossing = maximum_subarray_crossing(input, low, middle, high);
  //comparing which has the greatest sum
  if(accumulate(input.begin()+entirely_left.first,input.begin()+entirely_left.second,0) > accumulate(input.begin()+entirely_right.first,input.begin()+entirely_right.second,0) && accumulate(input.begin()+entirely_left.first,input.begin()+entirely_left.second,0)>accumulate(input.begin()+crossing.first,input.begin()+crossing.second,0))
    return entirely_left;
  else if(accumulate(input.begin()+entirely_right.first,input.begin()+entirely_right.second,0) > accumulate(input.begin()+entirely_left.first,input.begin()+entirely_left.second,0) && accumulate(input.begin()+entirely_right.first,input.begin()+entirely_right.second,0)>accumulate(input.begin()+crossing.first,input.begin()+crossing.second,0))
    return entirely_right;
  else
    return crossing;  
}

summed_span max_subarray_dbh(const std::vector<int>& input) {
  assert(!input.empty());
  // TODO: 
  //Base Case
  if(input.size() == 1)
    return summed_span(input.begin(), input.end());

  //General
  std::pair<int,int> b_e = maximum_subarray_recursive(input, 0, input.size()-1);

  // return summed_span(input.begin(), input.begin() + 1);
  return summed_span(input.begin()+b_e.first, input.begin()+b_e.second);
}

// Solve the subset sum problem: return a non-empty subset of input that adds
// up to exactly target. If no such subset exists, return an empty optional.
// input must not be empty, and must contain fewer than 64 elements.
// Note that the returned subset must never be empty, even if target == 0.
// This uses an exhaustive search algorithm that takes exponential O(n * 2^n)
// time.
std::optional<std::vector<int>>
subset_sum_exh(const std::vector<int>& input, int target) {

  assert(!input.empty());
  assert(input.size() < 64);

  // TODO: 
  std::vector<int> elements;
  std::vector<int> candidate;
  for(int i=0; i<input.size(); i++)
    elements.push_back(input[i]);
  for(uint64_t i=0; i<std::pow(2,elements.size()); i++){
    candidate.clear();
    for(uint64_t j=0; j<input.size(); j++){
      if(((i>>j)&1) == 1)  //issue
        candidate.push_back(elements[j]);
      if(candidate.size()>0 && std::accumulate(candidate.begin(),candidate.end(),0) == target)
        return candidate;
    }
  } 
  //return std::make_optional<std::vector<int>>();
  return std::nullopt;
 }

}
