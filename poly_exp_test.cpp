
///////////////////////////////////////////////////////////////////////////////
// poly_exp_test.hpp
//
// Unit tests for the functionality declared in poly_exp.hpp .
///////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <random>
#include <vector>

#include "gtest/gtest.h"

#include "poly_exp.hpp"

TEST(max_subarray_exh_SmallCases, max_subarray_exh_SmallCases) {
  { // small vectors, all positive, so include everything
    std::vector<int> one{1}, two{1, 2}, three{1, 2, 3};
    EXPECT_EQ(subarray::summed_span(one.begin(), one.end()),
              subarray::max_subarray_exh(one));
    EXPECT_EQ(subarray::summed_span(two.begin(), two.end()),
              subarray::max_subarray_exh(two));
    EXPECT_EQ(subarray::summed_span(three.begin(), three.end()),
              subarray::max_subarray_exh(three));
  }

  { // small vectors, one negative element, avoid that
    std::vector<int> two{-1, 2}, three{1, 2, -1}, five{1, 2, -9, 2, 2};
    EXPECT_EQ(subarray::summed_span(two.begin() + 1, two.begin() + 2),
              subarray::max_subarray_exh(two));
    EXPECT_EQ(subarray::summed_span(three.begin(), three.begin() + 2),
              subarray::max_subarray_exh(three));
    EXPECT_EQ(subarray::summed_span(five.begin() + 3, five.begin() + 5),
              subarray::max_subarray_exh(five));
  }

  { // all negative, chooses the one largest element
    std::vector<int> one{-1}, two{-1, -2}, three{-2, -1, -3}, four{-4, -2, -1, -3};
    EXPECT_EQ(subarray::summed_span(one.begin(), one.begin() + 1),
              subarray::max_subarray_exh(one));
    EXPECT_EQ(subarray::summed_span(two.begin(), two.begin() + 1),
              subarray::max_subarray_exh(two));
    EXPECT_EQ(subarray::summed_span(three.begin() + 1, three.begin() + 2),
              subarray::max_subarray_exh(three));
    EXPECT_EQ(subarray::summed_span(four.begin() + 2, four.begin() + 3),
              subarray::max_subarray_exh(four));
  }

  { // all tied and positive, picks everything
    std::vector<int> five{2, 2, 2, 2, 2};
    EXPECT_EQ(subarray::summed_span(five.begin(), five.begin() + 5),
              subarray::max_subarray_exh(five));
  }

  { // all tied and negative, picks one element
    std::vector<int> five{-2, -2, -2, -2, -2};
    auto result = subarray::max_subarray_exh(five);
    EXPECT_EQ(1, result.size());
    EXPECT_EQ(-2, result.sum());
  }
}

TEST(max_subarray_exh_LargeCases, max_subarray_exh_LargeCases) {
  { // nontrivial example from CLRS page 70
    std::vector<int> clrs{
      13, -3, -25, 20, -3, -16, -23, 18, 20, -7, 12, -5, -22, 15, -4, 7
    };
    EXPECT_EQ(subarray::summed_span(clrs.begin() + 7, clrs.begin() + 11),
              subarray::max_subarray_exh(clrs));
  }

  { // random instance
    std::vector<int> medium;
    std::mt19937 rng(0);
    std::uniform_int_distribution<> randint(-10, +10);
    for (unsigned i = 0; i < 1000; ++i) {
      medium.push_back(randint(rng));
    }
    ASSERT_EQ(1000, medium.size());
    ASSERT_EQ(164, subarray::max_subarray_exh(medium).sum());
  }
}

TEST(max_subarray_dbh_SmallCases, max_subarray_dbh_SmallCases) {
  { // small vectors, all positive, so include everything
    std::vector<int> one{1}, two{1, 2}, three{1, 2, 3};
    EXPECT_EQ(subarray::summed_span(one.begin(), one.end()),
              subarray::max_subarray_dbh(one));
    EXPECT_EQ(subarray::summed_span(two.begin(), two.end()),
              subarray::max_subarray_dbh(two));
    EXPECT_EQ(subarray::summed_span(three.begin(), three.end()),
              subarray::max_subarray_dbh(three));
  }

  { // small vectors, one negative element, avoid that
    std::vector<int> two{-1, 2}, three{1, 2, -1}, five{1, 2, -9, 2, 2};
    EXPECT_EQ(subarray::summed_span(two.begin() + 1, two.begin() + 2),
              subarray::max_subarray_dbh(two));
    EXPECT_EQ(subarray::summed_span(three.begin(), three.begin() + 2),
              subarray::max_subarray_dbh(three));
    EXPECT_EQ(subarray::summed_span(five.begin() + 3, five.begin() + 5),
              subarray::max_subarray_dbh(five));
  }

  { // all negative, chooses the one largest element
    std::vector<int> one{-1}, two{-1, -2}, three{-2, -1, -3}, four{-4, -2, -1, -3};
    EXPECT_EQ(subarray::summed_span(one.begin(), one.begin() + 1),
              subarray::max_subarray_dbh(one));
    EXPECT_EQ(subarray::summed_span(two.begin(), two.begin() + 1),
              subarray::max_subarray_dbh(two));
    EXPECT_EQ(subarray::summed_span(three.begin() + 1, three.begin() + 2),
              subarray::max_subarray_dbh(three));
    EXPECT_EQ(subarray::summed_span(four.begin() + 2, four.begin() + 3),
              subarray::max_subarray_dbh(four));
  }

  { // all tied and positive, picks everything
    std::vector<int> five{2, 2, 2, 2, 2};
    EXPECT_EQ(subarray::summed_span(five.begin(), five.begin() + 5),
              subarray::max_subarray_dbh(five));
  }

  { // all tied and negative, picks one element
    std::vector<int> five{-2, -2, -2, -2, -2};
    auto result = subarray::max_subarray_dbh(five);
    EXPECT_EQ(1, result.size());
    EXPECT_EQ(-2, result.sum());
  }
}

TEST(max_subarray_dbh_LargeCases, max_subarray_dbh_LargeCases) {
  { // nontrivial example from CLRS page 70
    std::vector<int> clrs{
      13, -3, -25, 20, -3, -16, -23, 18, 20, -7, 12, -5, -22, 15, -4, 7
    };
    EXPECT_EQ(subarray::summed_span(clrs.begin() + 7, clrs.begin() + 11),
              subarray::max_subarray_dbh(clrs));
  }

  { // random instance
    std::vector<int> medium;
    std::mt19937 rng(0);
    std::uniform_int_distribution<> randint(-10, +10);
    for (unsigned i = 0; i < 1000; ++i) {
      medium.push_back(randint(rng));
    }
    ASSERT_EQ(1000, medium.size());
    ASSERT_EQ(164, subarray::max_subarray_dbh(medium).sum());
  }
}

TEST(subset_sum_exh, subset_sum_exh) {
  // one element that is not target
  EXPECT_FALSE(subarray::subset_sum_exh({5}, 1));

  // one element that IS target
  {
    auto result = subarray::subset_sum_exh({5}, 5);
    ASSERT_TRUE(result.has_value());
    EXPECT_EQ(1, result->size());
    EXPECT_EQ(5, (*result)[0]);
  }

  // two positive elements that need to be used
  {
    auto result = subarray::subset_sum_exh({1, 3}, 4);
    ASSERT_TRUE(result.has_value());
    auto sorted = *result;
    std::sort(sorted.begin(), sorted.end());
    EXPECT_EQ(2, sorted.size());
    EXPECT_EQ(1, sorted[0]);
    EXPECT_EQ(3, sorted[1]);
  }

  // one positive and one negative that cancel out
  {
    auto result = subarray::subset_sum_exh({5, -2}, 3);
    ASSERT_TRUE(result.has_value());
    auto sorted = *result;
    std::sort(sorted.begin(), sorted.end());
    EXPECT_EQ(2, sorted.size());
    EXPECT_EQ(-2, sorted[0]);
    EXPECT_EQ(5, sorted[1]);
  }

  // doesn't pick empty set
  {
    EXPECT_FALSE(subarray::subset_sum_exh({1, 2, 3}, 0));
  }

  // small instances with no solution
  {
    EXPECT_FALSE(subarray::subset_sum_exh({2, -1}, 0));
    EXPECT_FALSE(subarray::subset_sum_exh({2, 4, 6}, 5));
    EXPECT_FALSE(subarray::subset_sum_exh({8, 2, -5, 3}, 1));
  }

  // Wikipedia https://en.wikipedia.org/wiki/Subset_sum_problem
  {
    auto result = subarray::subset_sum_exh({-7, -3, -2, 5, 8}, 0);
    ASSERT_TRUE(result.has_value());
    EXPECT_FALSE(result->empty());
    EXPECT_EQ(0, std::accumulate(result->begin(), result->end(), 0));
  }

  // CLRS page 1097
  {
    std::vector<int> input{
      1, 2, 7, 14, 49, 98, 343, 686, 2409, 2793, 16808, 17206, 117705, 117993
    };
    auto result = subarray::subset_sum_exh(input, 138457);
    ASSERT_TRUE(result.has_value());
    EXPECT_FALSE(result->empty());
    EXPECT_EQ(138457, std::accumulate(result->begin(), result->end(), 0));
  }
}
