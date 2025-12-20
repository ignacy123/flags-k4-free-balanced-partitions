#include "config.hpp"
#include "flag-calculator.hpp"
#include "flag.hpp"
#include "option-parser.hpp"
#include "problem.hpp"
#include "solver.hpp"
#include "sum-utils.hpp"
#include <cassert>

double BOUND = 1. / 9;

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

  flag_coeff edge_left("3 2  2 2  1");
  flag_coeff edge_right("3 2  2 1  2");
  flag_coeff edge_both("3 2  2 2  2");
  flag_coeff edge_neither("3 2  2 1  1");

  auto symmetric_cut_on_edge = rooted_cut<2>({edge_left}, {edge_right},
                                             {edge_neither, edge_both}, BOUND);
  problem.add_constraint(symmetric_cut_on_edge);

  auto asymmetric_cut_on_edge = rooted_cut<2>({edge_left, edge_neither},
                                              {edge_right}, {edge_both}, BOUND);
  problem.add_constraint(asymmetric_cut_on_edge);

  auto extra_cut_on_edge = rooted_cut<2>({edge_left}, {edge_both},
                                         {edge_right, edge_neither}, BOUND);
  problem.add_constraint(extra_cut_on_edge);

  auto second_extra_cut_on_edge = rooted_cut<2>(
      {edge_left, edge_neither}, {edge_both}, {edge_right}, BOUND);
  problem.add_constraint(second_extra_cut_on_edge);

  solve_sdp_for_problem(problem.get_constraints(), problem.get_objective());
}
