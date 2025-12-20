#include "balanced-bipartition-utils.hpp"
#include "flag-calculator.hpp"
#include "flag.hpp"
#include "option-parser.hpp"
#include "problem.hpp"
#include "solver.hpp"
#include <cassert>
#include <vector>

double BOUND = 425. / 6498;

template <int root_size>
CutInfo<root_size - 1, root_size + 4> casted_cut(vector<flag_coeff> left_side,
                                                 vector<flag_coeff> right_side,
                                                 vector<flag_coeff> random) {

  auto regular_cut_halves =
      prepare_cut_halves<root_size>(left_side, right_side, random, BOUND);

  return CutInfo<root_size - 1, root_size + 4>{
      regular_cut_halves.left_side
          .template project<root_size - 1, root_size + 4>(),
      regular_cut_halves.right_side
          .template project<root_size - 1, root_size + 4>(),
      regular_cut_halves.lower_bound
          .template project<root_size - 1, root_size + 4>(),
  };
}

CutInfo<0, 4> degree_cut(const vector<flag_coeff> &left,
                         const vector<flag_coeff> &random) {
  FlagVector<0, 1> sum;
  FlagVector<0, 1> ones(flag(), 1.);
  for (auto flag : left) {
    assert(flag.g.m_vertices == 1);
    sum += flag;
  }
  for (auto flag : random) {
    assert(flag.g.m_vertices == 1);
    sum += flag;
  }
  assert(sum == ones);
  auto both_left = get_edges_inside<0>(flag(), left);
  auto neither = get_edges_inside<0>(flag(), random);
  auto one_left = get_edges_between<0>(flag(), left, random);
  auto single_left = FlagVector<0, 1>::from_vector(flag(), left);
  auto space_left = 1. / 2 - single_left;
  FlagVector<0, 1> space_right(flag(), 1. / 2);
  auto random_vector = FlagVector<0, 1>::from_vector(flag(), random);
  assert(random_vector == space_left + space_right);

  auto left_side = both_left * random_vector * random_vector;
  left_side += one_left * space_left * random_vector;
  left_side += neither * space_left * space_left;

  auto right_side = neither * space_right * space_right;
  FlagVector<0, 2> lower_bound(flag(), 2 * BOUND);
  return CutInfo<0, 4>{left_side, right_side,
                       lower_bound * random_vector * random_vector};
}

int main(int argc, char *argv[]) {
  ProblemConfig::instance().data_directory = "wlike";
  ProblemConfig::instance().csdp_binary = "csdp";
  ProblemConfig::instance().sdpa_binary = "sdpa";
  ProblemConfig::instance().upper_bound = false;
  parse_options(argc, argv, ProblemConfig::instance());
  Problem problem;

  auto edge_vector = FlagVector<0, 2>(flag("2 0  0 0  2"));
  auto objective = 1. - edge_vector;

  problem.add_objective(objective);

  // Color 0 - any color.
  // Color 1 - low degree (<= 1 / 2).
  // Color 2 - high degree (>= 1 / 2).
  auto degree_blue = 1. / 2 - (FlagVector<1, 2>(flag("2 1  1 0  2")));
  auto degree_red = FlagVector<1, 2>(flag("2 1  2 0  2")) - 1. / 2;
  problem.add_constraint(degree_blue);
  problem.add_constraint(degree_red);

  // auto random_cut = FlagVector<0, 2>(flag("2 0  0 0  2")) - 8 * BOUND;
  // problem.add_constraint(random_cut);

  vector<flag_coeff> vertices_less, vertices_more;
  if (ProblemConfig::instance().case_number & 1) {
    cerr << "Assuming there's more vertices of low degree." << endl;
    vertices_more.push_back(flag_coeff("1 0  1"));
    vertices_less.push_back(flag_coeff("1 0  2"));
  } else {
    cerr << "Assuming there's more vertices of high degree." << endl;
    vertices_more.push_back(flag_coeff("1 0  2"));
    vertices_less.push_back(flag_coeff("1 0  1"));
  }
  problem.add_constraint(FlagVector<0, 1>::from_vector(flag(), vertices_more) -
                         1. / 2);
  auto cut_on_degrees = degree_cut(vertices_less, vertices_more);
  // if (ProblemConfig::instance().case_number >> 1 & 1) {
  //   cerr << "Assuming that in the cut degree there is more edges in the half
  //   "
  //           "with vertices of the more prevalent degree."
  //        << endl;
  //   problem.add_constraint(cut_on_degrees.right_side -
  //                          cut_on_degrees.left_side);
  //   problem.add_constraint(cut_on_degrees.right_side -
  //                          cut_on_degrees.lower_bound);
  // } else {
  //   cerr << "Assuming that in the cut degree there is more edges in the half
  //   "
  //           "with vertices of both kind."
  //        << endl;
  //   problem.add_constraint(cut_on_degrees.left_side -
  //                          cut_on_degrees.right_side);
  //   problem.add_constraint(cut_on_degrees.left_side -
  //                          cut_on_degrees.lower_bound);
  // }

  flag_coeff blue_vertex_connected_blue("2 1  1 1  2");
  flag_coeff blue_vertex_connected_red("2 1  1 2  2");
  flag_coeff blue_vertex_disconnected_blue("2 1  1 1  1");
  flag_coeff blue_vertex_disconnected_red("2 1  1 2  1");

  auto cut_info_blue_vertex = prepare_cut_halves<1>(
      {blue_vertex_connected_red}, {},
      {blue_vertex_disconnected_blue, blue_vertex_disconnected_red}, BOUND);
  problem.add_constraint(cut_info_blue_vertex.left_side -
                         cut_info_blue_vertex.right_side);
  problem.add_constraint(cut_info_blue_vertex.left_side -
                         cut_info_blue_vertex.lower_bound);

  flag_coeff red_vertex_disconnected_blue("2 1  2 1  1");
  flag_coeff red_vertex_disconnected_red("2 1  2 2  1");
  flag_coeff red_vertex_connected_blue("2 1  2 1  2");
  flag_coeff red_vertex_connected_red("2 1  2 2  2");

  auto first_cut_on_red_vertex = prepare_cut_halves<1>(
      {red_vertex_disconnected_blue, red_vertex_disconnected_red},
      {red_vertex_connected_red}, {red_vertex_connected_blue}, BOUND);
  // We don't need to care about the other case, because right side is a subset
  // of a triangle-free graph, which allows us to get a very strong bound even
  // with a simple application of Mantel's theorem.
  problem.add_constraint(first_cut_on_red_vertex.left_side -
                         first_cut_on_red_vertex.right_side);
  problem.add_constraint(first_cut_on_red_vertex.left_side -
                         first_cut_on_red_vertex.lower_bound);

  flag_coeff red_edge_left_only_blue("3 2  2 2 1  2 2  1");
  flag_coeff red_edge_left_only_red("3 2  2 2 2  2 2  1");
  flag_coeff red_edge_right_only_blue("3 2  2 2 1  2 1  2");
  flag_coeff red_edge_right_only_red("3 2  2 2 2  2 1  2");
  flag_coeff red_edge_neither_blue("3 2  2 2 1  2 1  1", 0.5);
  flag_coeff red_edge_neither_red("3 2  2 2 2  2 1  1", 0.5);
  flag_coeff red_edge_both_blue("3 2  2 2 1  2 2  2");
  flag_coeff red_edge_both_red("3 2  2 2 2  2 2  2");
  auto casted_cut_on_red_edge =
      casted_cut<2>({red_edge_left_only_blue, red_edge_left_only_red,
                     red_edge_neither_blue, red_edge_neither_red},
                    {red_edge_right_only_blue, red_edge_right_only_red,
                     red_edge_neither_blue, red_edge_neither_red},
                    {red_edge_both_blue, red_edge_both_red});
  problem.add_constraint(casted_cut_on_red_edge.left_side -
                         casted_cut_on_red_edge.right_side);
  problem.add_constraint(casted_cut_on_red_edge.left_side -
                         casted_cut_on_red_edge.lower_bound);

  solve_sdp_for_problem(problem.get_constraints(), problem.get_objective());
}
