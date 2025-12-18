#include "config.hpp"
#include "flag-calculator.hpp"
#include "flag.hpp"
#include "option-parser.hpp"
#include "problem.hpp"
#include "solver.hpp"
#include <cassert>
#include <vector>

double BOUND = 2. / 9;

flag get_type(const vector<flag> &left_side, const vector<flag> &right_side,
              const vector<flag> &random) {
  flag type;
  if (left_side.size() > 0) {
    left_side[0].get_type_subflag(type);
  } else if (right_side.size() > 0) {
    right_side[0].get_type_subflag(type);
  } else if (random.size() > 0) {
    random[0].get_type_subflag(type);
  }

  if (left_side.size() > 0) {
    flag left_type;
    left_side[0].get_type_subflag(left_type);
    assert(type.is_identical_to(left_type));
  } else if (right_side.size() > 0) {
    flag right_type;
    right_side[0].get_type_subflag(right_type);
    assert(type.is_identical_to(right_type));
  } else if (random.size() > 0) {
    flag random_type;
    random[0].get_type_subflag(random_type);
    assert(type.is_identical_to(random_type));
  }
  return type;
}

template <int root_size> struct CutComponents {
public:
  FlagVector<root_size, root_size + 2> fixed;
  FlagVector<root_size, root_size + 1> denominator;
  FlagVector<root_size, root_size + 1> space_left;
  FlagVector<root_size, root_size + 1> space_right;
  FlagVector<root_size, root_size + 2> one_left;
  FlagVector<root_size, root_size + 2> one_right;
  FlagVector<root_size, root_size + 2> neither;
};

template <int root_size>
CutComponents<root_size> prepare_cut_components(vector<flag> left_side,
                                                vector<flag> right_side,
                                                vector<flag> random) {
  flag type = get_type(left_side, right_side, random);
  FlagVector<root_size, root_size + 1> sum(type);
  for (flag from_left : left_side) {
    assert(from_left.have_same_type(type));
    assert(from_left.m_Theta == root_size);
    assert(from_left.m_vertices == root_size + 1);
    sum += from_left;
  }
  for (flag from_right : right_side) {
    assert(from_right.have_same_type(type));
    assert(from_right.m_Theta == root_size);
    assert(from_right.m_vertices == root_size + 1);
    sum += from_right;
  }
  for (flag from_random : random) {
    assert(from_random.have_same_type(type));
    assert(from_random.m_Theta == root_size);
    assert(from_random.m_vertices == root_size + 1);
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

  FlagVector<root_size, root_size + 2> fixed = both_left + both_right - BOUND;
  FlagVector<root_size, root_size + 1> denominator =
      FlagVector<root_size, root_size + 1>::from_vector(type, random);
  FlagVector<root_size, root_size + 1> single_left =
      FlagVector<root_size, root_size + 1>::from_vector(type, left_side);
  FlagVector<root_size, root_size + 1> single_right =
      FlagVector<root_size, root_size + 1>::from_vector(type, right_side);
  FlagVector<root_size, root_size + 1> space_left = 1. / 2 - single_left;
  FlagVector<root_size, root_size + 1> space_right = 1. / 2 - single_right;

  return CutComponents<root_size>{fixed,       denominator, space_left,
                                  space_right, one_left,    one_right,
                                  neither};
}

template <int root_size>
FlagVector<root_size, root_size + 3>
simplified_rooted_cut(vector<flag> left_side, vector<flag> right_side,
                      vector<flag> random) {
  auto cut_components =
      prepare_cut_components<root_size>(left_side, right_side, random);

  auto result = (cut_components.fixed * cut_components.denominator);
  result += cut_components.space_left * cut_components.one_left;
  result += cut_components.space_right * cut_components.one_right;
  result += (cut_components.space_left + cut_components.space_right) *
            cut_components.neither;
  return result;
}

template <int root_size>
FlagVector<root_size, root_size + 4> rooted_cut(vector<flag> left_side,
                                                vector<flag> right_side,
                                                vector<flag> random) {
  auto cut_components =
      prepare_cut_components<root_size>(left_side, right_side, random);

  auto result = (cut_components.fixed * cut_components.denominator *
                 cut_components.denominator);
  result += cut_components.space_left * cut_components.one_left *
            cut_components.denominator;
  result += cut_components.space_right * cut_components.one_right *
            cut_components.denominator;
  result += (cut_components.space_left * cut_components.space_left +
             cut_components.space_right * cut_components.space_right) *
            cut_components.neither;
  return result;
}

int main(int argc, char *argv[]) {
  ProblemConfig::instance().data_directory = "k4-free";
  ProblemConfig::instance().csdp_binary = "csdp-no-accelerate";
  ProblemConfig::instance().sdpa_binary = "sdpa";
  parse_options(argc, argv, ProblemConfig::instance());
  Problem problem(2 * 2 * 2 * 3 * 3 * 3 * 3 * 5, true);

  auto edge_vector = FlagVector<0, 2>(flag("2 0  2"));
  auto objective = (edge_vector * (2. / 3 - edge_vector)).project<0, 5>();
  objective += edge_vector * (2. / 9 - FlagVector<0, 3>(flag("3 0  2 2  2")));
  auto triangle_common_edge0 = FlagVector<3, 4>(flag("4 3  2 2 2  2 2  1"));
  auto triangle_common_edge1 = FlagVector<3, 4>(flag("4 3  2 2 1  2 2  2"));
  auto triangle_common_edge2 = FlagVector<3, 4>(flag("4 3  2 2 2  2 1  2"));
  objective += (1. - triangle_common_edge0 - triangle_common_edge1 -
                triangle_common_edge2)
                   .project<0, 5>();
  objective = 1. - objective;

  problem.add_objective(objective);

  flag edge_left("3 2  2 2  1");
  flag edge_right("3 2  2 1  2");
  flag edge_both("3 2  2 2  2");
  flag edge_neither("3 2  2 1  1");

  auto symmetric_cut_on_edge =
      rooted_cut<2>({edge_left}, {edge_right}, {edge_neither, edge_both});
  problem.add_constraint(symmetric_cut_on_edge);

  auto asymmetric_cut_on_edge =
      rooted_cut<2>({edge_left, edge_neither}, {edge_right}, {edge_both});
  problem.add_constraint(asymmetric_cut_on_edge);

  auto extra_cut_on_edge =
      rooted_cut<2>({edge_left}, {edge_both}, {edge_right, edge_neither});
  problem.add_constraint(extra_cut_on_edge);

  auto second_extra_cut_on_edge =
      rooted_cut<2>({edge_left, edge_neither}, {edge_both}, {edge_right});
  problem.add_constraint(second_extra_cut_on_edge);

  solve_sdp_for_problem(problem.get_constraints(), problem.get_objective());
}
