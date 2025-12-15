#pragma once
#include "flag.hpp"
#include <iostream>

using namespace std;

struct MultiplicationPartialResult {
public:
  int idx_multiplied;
  int idx_left;
  int idx_right;
  int good_maps;
  int all_maps;
};

std::istream &operator>>(std::istream &in,
                         MultiplicationPartialResult &partial_result);

struct ProjectionPartialResult {
public:
  int idx_projected;
  int idx_original;
  int good_maps;
  int all_maps;
};

std::istream &operator>>(std::istream &in,
                         ProjectionPartialResult &partial_result);

struct FlagProduct {
  int type_idx;
  int f1_idx;
  int f2_idx;
};

std::ostream &operator<<(std::ostream &stream, const FlagProduct &flag_product);

std::istream &operator>>(std::istream &in, FlagProduct &flag_product);

bool operator==(const FlagProduct &a, const FlagProduct &b);

struct FlagProductPartialResult {
  FlagProduct product;
  int count;
};

std::istream &operator>>(std::istream &in,
                         FlagProductPartialResult &partial_result);

const vector<ProjectionPartialResult> &
get_projection_cache(flag type, int vertices, int new_vertices,
                     int remaining_labeled);

// TODO: Unify multiplication and flag product cache. They do the same thing.
const vector<MultiplicationPartialResult> &
get_multiplication_cache(flag type, int size_left, int size_right);

const vector<FlagProductPartialResult> &
get_flag_product_cache(int Kn, int idx_in_list, int typesize);
