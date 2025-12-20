#include "balanced-bipartition-utils.hpp"
#include "flag.hpp"

template <int root_size>
FlagVector<root_size, root_size + 3>
simplified_rooted_cut(vector<flag_coeff> left_side,
                      vector<flag_coeff> right_side,
                      vector<flag_coeff> random) {
  auto cut_components =
      prepare_cut_components<root_size>(left_side, right_side, random);

  auto result = (cut_components.fixed_left + cut_components.fixed_right -
                 cut_components.lower_bound) *
                cut_components.denominator;
  result += cut_components.space_left * cut_components.one_left;
  result += cut_components.space_right * cut_components.one_right;
  result += (cut_components.space_left + cut_components.space_right) *
            cut_components.neither;
  return result;
}

template <int root_size>
FlagVector<root_size, root_size + 4>
rooted_cut(vector<flag_coeff> left_side, vector<flag_coeff> right_side,
           vector<flag_coeff> random, double bound) {
  auto cut_components =
      prepare_cut_components<root_size>(left_side, right_side, random, bound);

  auto result = cut_components.fixed_left * cut_components.denominator *
                cut_components.denominator;
  result += cut_components.space_left * cut_components.one_left *
            cut_components.denominator;
  result += cut_components.fixed_right * cut_components.denominator *
            cut_components.denominator;
  result += cut_components.space_right * cut_components.one_right *
            cut_components.denominator;
  result += (cut_components.space_left * cut_components.space_left +
             cut_components.space_right * cut_components.space_right) *
            cut_components.neither;

  result -= cut_components.lower_bound * cut_components.denominator *
            cut_components.denominator;

  return result;
}
