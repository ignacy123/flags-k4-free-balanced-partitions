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
  ProblemConfig::instance().csdp_binary = "csdp";
  ProblemConfig::instance().sdpa_binary = "sdpa";
  ProblemConfig::instance().upper_bound = false;
  parse_options(argc, argv, ProblemConfig::instance());
  int scaling_factor = 2 * 2 * 2 * 2 * 3 * 3 * 3 * 3 * 5;
  Problem problem(scaling_factor, true);

  auto edge_vector = FlagVector<0, 2>(flag("2 0  2 2  2"));
  auto objective = (edge_vector * (2. / 3 - edge_vector)).project<0, 5>();
  objective +=
      edge_vector * (2. / 9 - FlagVector<0, 3>(flag("3 0  2 2 2  2 2  2")));
  auto triangle_common_edge0 =
      FlagVector<3, 4>(flag("4 3  2 2 2 0  2 2 2  2 2  1"));
  auto triangle_common_edge1 =
      FlagVector<3, 4>(flag("4 3  2 2 2 0  2 2 1  2 2  2"));
  auto triangle_common_edge2 =
      FlagVector<3, 4>(flag("4 3  2 2 2 0  2 2 2  2 1  2"));
  objective += (1. - triangle_common_edge0 - triangle_common_edge1 -
                triangle_common_edge2)
                   .project<0, 5>();
  objective = 1. - objective;

  problem.add_objective(objective);

  // Color 0 - any color.
  // Color 1 - vertices of low degree (<= 1/2), also called blue.
  // Color 2 - vertices of high degree (>= 1/2), also called red.

  auto degree_blue = 1. / 2 - (FlagVector<1, 2>(flag("2 1  1 0  2")));
  auto degree_red = FlagVector<1, 2>(flag("2 1  2 0  2")) - 1. / 2;
  problem.add_constraint(degree_blue);
  problem.add_constraint(degree_red);

  flag_coeff blue_vertex_connected("2 1  1 0  2");
  flag_coeff blue_vertex_disconnected("2 1  1 0  1");

  auto cut_on_blue_vertex = rooted_cut<1>({blue_vertex_connected}, {},
                                          {blue_vertex_disconnected}, BOUND);
  problem.add_constraint(cut_on_blue_vertex);

  flag_coeff red_edge_left("3 2  2 2 0  2 2  1");
  flag_coeff red_edge_right("3 2  2 2 0  2 1  2");
  flag_coeff red_edge_both("3 2  2 2 0  2 2  2");
  flag_coeff red_edge_neither("3 2  2 2 0  2 1  1");

  auto symmetric_cut_on_red_edge =
      rooted_cut<2>({red_edge_left}, {red_edge_right},
                    {red_edge_neither, red_edge_both}, BOUND);
  problem.add_constraint(symmetric_cut_on_red_edge);

  auto asymmetric_cut_on_red_edge =
      rooted_cut<2>({red_edge_left, red_edge_neither}, {red_edge_right},
                    {red_edge_both}, BOUND);
  problem.add_constraint(asymmetric_cut_on_red_edge);

  auto cut_det_both_red_edge =
      rooted_cut<2>({red_edge_left}, {red_edge_both},
                    {red_edge_right, red_edge_neither}, BOUND);
  problem.add_constraint(cut_det_both_red_edge);

  auto cut_extra_both_red_edge =
      rooted_cut<2>({red_edge_left, red_edge_neither}, {red_edge_both},
                    {red_edge_right}, BOUND);
  problem.add_constraint(cut_extra_both_red_edge);

  // flag cherry_none("4 3  2 2 2 0  2 1 1  2 1  1");
  // flag cherry_left_only("4 3  2 2 2 0  2 1 2  2 1  1");
  // flag cherry_middle_only("4 3  2 2 2 0  2 1 1  2 2  1");
  // flag cherry_right_only("4 3  2 2 2 0  2 1 1  2 1  2");
  // flag cherry_left_middle("4 3  2 2 2 0  2 1 2  2 2  1");
  // flag cherry_middle_right("4 3  2 2 2 0  2 1 1  2 2  2");
  // flag cherry_left_right("4 3  2 2 2 0  2 1 2  2 1  2");
  // flag cherry_all("4 3  2 2 2 0  2 1 2  2 2  2");

  // auto cut_on_cherry = simplified_rooted_cut<3>(
  //     {cherry_left_only, cherry_left_right, cherry_right_only},
  //     {cherry_none, cherry_middle_only, cherry_middle_right},
  //     {cherry_left_middle, cherry_all});
  // problem.add_constraint(cut_on_cherry);

  solve_sdp_for_problem(problem.get_constraints(), problem.get_objective());
}
