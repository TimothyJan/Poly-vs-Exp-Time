///////////////////////////////////////////////////////////////////////////////
// poly_exp_timing.cpp
//
// Example code showing how to run each algorithm while measuring
// elapsed times precisely. You should modify this program to gather
// all of your experimental data.
//
///////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <cassert>
#include <iostream>
#include <random>
#include <vector>

#include "timer.hpp"

#include "poly_exp.hpp"

void print_bar() {
  std::cout << std::string(79, '-') << std::endl;
}

void print_int_vector(const std::vector<int>& v) {
  std::cout << "{";
  bool first = true;
  for (auto x : v) {
    if (first) {
      first = false;
    } else {
      std::cout << ", ";
    }
    std::cout << x;
  }
  std::cout << "}";
}

int main() {

  // Feel free to change these constants to suit your needs.
  const size_t n = 20,
               print_input_limit = 100,
               max_subarray_exh_limit = 10000,
               subset_sum_exh_limit = 28;

  assert(n > 0);

  std::vector<int> subarray_input, subset_sum_input;
  {
    std::mt19937 rng(0); // Use a hardcoded seed for reproducibility between runs.
    const unsigned billion = 1E9;
    std::uniform_int_distribution<> subarray_dist(-100, +100),
                                    subset_sum_dist(-billion, +billion);
    for (size_t i = 0; i < n; ++i) {
      subarray_input.push_back(subarray_dist(rng));
      subset_sum_input.push_back(subset_sum_dist(rng));
    }
  }
  assert(n == subarray_input.size());
  assert(n == subset_sum_input.size());

  Timer timer;
  double elapsed;

  print_bar();
  std::cout << "n = " << n << std::endl;
  if (n > print_input_limit) {
    std::cout << "(input too large to print)" << std::endl;
  } else {
    std::cout << "subarray_input: " << std::endl;
    print_int_vector(subarray_input);
    std::cout << std::endl;
    std::cout << "subset_sum_input: " << std::endl;
    print_int_vector(subset_sum_input);
    std::cout << std::endl;
  }

  print_bar();
  std::cout << "max_subarray_dc" << std::endl;
  {
    timer.reset();
    auto solution = subarray::max_subarray_dbh(subarray_input);
    elapsed = timer.elapsed();
    std::cout << "solution: " << solution << std::endl
              << "elapsed time=" << elapsed << " seconds" << std::endl;
  }

  print_bar();
  std::cout << "max_subarray_exh" << std::endl;
  if (n > max_subarray_exh_limit) {
    std::cout << "(skipped because n > " << max_subarray_exh_limit << ")" << std::endl;
  } else {
    timer.reset();
    auto solution = subarray::max_subarray_exh(subarray_input);
    elapsed = timer.elapsed();
    std::cout << "solution: " << solution << std::endl
              << "elapsed time=" << elapsed << " seconds" << std::endl;
  }

  print_bar();
  std::cout << "subset_sum_exh" << std::endl;
  if (n > subset_sum_exh_limit) {
    std::cout << "(skipped because n > " << subset_sum_exh_limit << ")" << std::endl;
  } else {
    timer.reset();
    auto solution = subarray::subset_sum_exh(subset_sum_input, 1);
    elapsed = timer.elapsed();
    std::cout << "solution:" << std::endl;
    if (!solution.has_value()) {
      std::cout << "(no solution)" << std::endl;
    } else {
      print_int_vector(*solution);
      std::cout << std::endl << "sum="
                << std::accumulate(solution->begin(), solution->end(), 0)
                << std::endl;
    }
    std::cout << "elapsed time=" << elapsed << " seconds" << std::endl;
  }

  print_bar();

  return 0;
}
