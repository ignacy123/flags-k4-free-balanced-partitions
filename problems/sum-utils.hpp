#include "balanced-bipartition-utils.hpp"
#include "flag.hpp"

template <int root_size>
FlagVector<root_size, root_size + 3>
simplified_rooted_cut(vector<flag_coeff> left_side,
                      vector<flag_coeff> right_side, vector<flag_coeff> random,
                      double bound) {
  auto cut_halves = prepare_cut_halves_simplified<root_size>(
      left_side, right_side, random, bound);
  return cut_halves.left_side + cut_halves.right_side - cut_halves.lower_bound;
}

template <int root_size>
FlagVector<root_size, root_size + 4>
rooted_cut(vector<flag_coeff> left_side, vector<flag_coeff> right_side,
           vector<flag_coeff> random, double bound) {
  auto cut_halves =
      prepare_cut_halves<root_size>(left_side, right_side, random, bound);
  return cut_halves.left_side + cut_halves.right_side - cut_halves.lower_bound;
}
