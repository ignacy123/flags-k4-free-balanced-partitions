#pragma once
#include "flag-calculator.hpp"
#include "flag.hpp"

inline flag get_type(const vector<flag_coeff> &left_side,
                     const vector<flag_coeff> &right_side,
                     const vector<flag_coeff> &random) {
  flag type;
  if (left_side.size() > 0) {
    left_side[0].g.get_type_subflag(type);
  } else if (right_side.size() > 0) {
    right_side[0].g.get_type_subflag(type);
  } else if (random.size() > 0) {
    random[0].g.get_type_subflag(type);
  }

  if (left_side.size() > 0) {
    flag left_type;
    left_side[0].g.get_type_subflag(left_type);
    assert(type.is_identical_to(left_type));
  } else if (right_side.size() > 0) {
    flag right_type;
    right_side[0].g.get_type_subflag(right_type);
    assert(type.is_identical_to(right_type));
  } else if (random.size() > 0) {
    flag random_type;
    random[0].g.get_type_subflag(random_type);
    assert(type.is_identical_to(random_type));
  }
  return type;
}

template <int root_size> struct CutComponents {
public:
  FlagVector<root_size, root_size + 2> fixed_left;
  FlagVector<root_size, root_size + 2> fixed_right;
  FlagVector<root_size, root_size + 1> denominator;
  FlagVector<root_size, root_size + 1> space_left;
  FlagVector<root_size, root_size + 1> space_right;
  FlagVector<root_size, root_size + 2> one_left;
  FlagVector<root_size, root_size + 2> one_right;
  FlagVector<root_size, root_size + 2> neither;
  FlagVector<root_size, root_size + 2> lower_bound;
};

template <int root_size>
CutComponents<root_size> prepare_cut_components(vector<flag_coeff> left_side,
                                                vector<flag_coeff> right_side,
                                                vector<flag_coeff> random,
                                                double bound) {
  flag type = get_type(left_side, right_side, random);
  FlagVector<root_size, root_size + 1> sum(type);
  for (auto from_left : left_side) {
    assert(from_left.g.have_same_type(type));
    assert(from_left.g.m_Theta == root_size);
    assert(from_left.g.m_vertices == root_size + 1);
    sum += from_left;
  }
  for (auto from_right : right_side) {
    assert(from_right.g.have_same_type(type));
    assert(from_right.g.m_Theta == root_size);
    assert(from_right.g.m_vertices == root_size + 1);
    sum += from_right;
  }
  for (auto from_random : random) {
    assert(from_random.g.have_same_type(type));
    assert(from_random.g.m_Theta == root_size);
    assert(from_random.g.m_vertices == root_size + 1);
    sum += from_random;
  }
  FlagVector<root_size, root_size + 1> ones(type, 1.);
  // Check if all possible combinations have been passed.
  assert(ones == sum);

  FlagVector<root_size, root_size + 2> both_left =
      get_edges_inside<root_size>(type, left_side);
  FlagVector<root_size, root_size + 2> both_right =
      get_edges_inside<root_size>(type, right_side);
  FlagVector<root_size, root_size + 2> neither =
      get_edges_inside<root_size>(type, random);
  FlagVector<root_size, root_size + 2> one_left =
      get_edges_between<root_size>(type, left_side, random);
  FlagVector<root_size, root_size + 2> one_right =
      get_edges_between<root_size>(type, right_side, random);

  FlagVector<root_size, root_size + 1> denominator =
      FlagVector<root_size, root_size + 1>::from_vector(type, random);
  FlagVector<root_size, root_size + 1> single_left =
      FlagVector<root_size, root_size + 1>::from_vector(type, left_side);
  FlagVector<root_size, root_size + 1> single_right =
      FlagVector<root_size, root_size + 1>::from_vector(type, right_side);
  FlagVector<root_size, root_size + 1> space_left = 1. / 2 - single_left;
  FlagVector<root_size, root_size + 1> space_right = 1. / 2 - single_right;

  return CutComponents<root_size>{
      both_left,
      both_right,
      denominator,
      space_left,
      space_right,
      one_left,
      one_right,
      neither,
      FlagVector<root_size, root_size + 2>(type, 2 * bound)};
}

template <int root_size, int final_size> struct CutInfo {
  FlagVector<root_size, final_size> left_side;
  FlagVector<root_size, final_size> right_side;
  FlagVector<root_size, final_size> lower_bound;
};
