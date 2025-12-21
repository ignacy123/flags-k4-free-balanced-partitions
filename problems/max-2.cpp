#include "balanced-bipartition-utils.hpp"
#include "flag-calculator.hpp"
#include "flag.hpp"
#include "option-parser.hpp"
#include "problem.hpp"
#include "solver.hpp"
#include <cassert>
#include <string>
#include <vector>

double BOUND = 5. / 72;

template <int root_size>
CutInfo<root_size - 1, root_size + 3>
simple_projected_cut(vector<flag_coeff> left_side,
                     vector<flag_coeff> right_side, vector<flag_coeff> random) {

  auto regular_cut_halves = prepare_cut_halves_simplified<root_size>(
      left_side, right_side, random, BOUND);

  return CutInfo<root_size - 1, root_size + 3>{
      regular_cut_halves.left_side
          .template project<root_size - 1, root_size + 3>(),
      regular_cut_halves.right_side
          .template project<root_size - 1, root_size + 3>(),
      regular_cut_halves.lower_bound
          .template project<root_size - 1, root_size + 3>(),
  };
}

CutInfo<1, 5> projected_edge_cut(int color, vector<int> other_colors) {

  FlagVector<1, 5> left_side(flag("1 1  " + to_string(color)));
  FlagVector<1, 5> right_side(flag("1 1  " + to_string(color)));
  FlagVector<1, 5> lower_bound(flag("1 1  " + to_string(color)));
  // Heavily relying on the linearity of averaging operator.
  for (int other_color : other_colors) {
    flag_coeff left_only("3 2  " + to_string(color) + " " +
                         to_string(other_color) + " 0  2 2  1");
    flag_coeff right_only("3 2  " + to_string(color) + " " +
                          to_string(other_color) + " 0  2 1  2");
    flag_coeff neither("3 2  " + to_string(color) + " " +
                           to_string(other_color) + " 0  2 1  1",
                       0.5);
    flag_coeff both("3 2  " + to_string(color) + " " + to_string(other_color) +
                    " 0  2 2  2");
    auto cut_halves = prepare_cut_halves_simplified<2>(
        {left_only, neither}, {right_only, neither}, {both}, BOUND);
    left_side += cut_halves.left_side.project<1, 5>();
    right_side += cut_halves.right_side.project<1, 5>();
    lower_bound += cut_halves.lower_bound.project<1, 5>();
  }

  return CutInfo<1, 5>{left_side, right_side, lower_bound};
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
  ProblemConfig::instance().data_directory = "k4-free";
  ProblemConfig::instance().csdp_binary = "csdp";
  ProblemConfig::instance().sdpa_binary = "sdpa";
  ProblemConfig::instance().upper_bound = false;
  parse_options(argc, argv, ProblemConfig::instance());
  Problem problem;

  auto edge_vector = FlagVector<0, 2>(flag("2 0  0 0  2"));
  auto objective = 1. - edge_vector;

  problem.add_objective(objective);

  // Color 0 - any color.
  // Color 1 - vertices of low degree (<= 1/2), also called blue. These vertices
  // have more edges on one side of the first cut that we consider, which only
  // works for vertices of low degree.
  // Color 2 - vertices of low degree (<= 1/2), also called cyan. These vertices
  // have more edges on one side of the first cut that we consider, which only
  // works for vertices of low degree.
  // Color 3 - vertices of high degree (>= 1/2), also called red. These vertices
  // have more edges on one side of the second cut that we consider, which only
  // works for vertices of high degree.
  // Color 4 - vertices of high degree (>= 1/2), also called magenta. These
  // vertices have more edges on one side of the second cut that we consider,
  // which only works for vertices of high degree.

  auto degree_blue = 1. / 2 - (FlagVector<1, 2>(flag("2 1  1 0  2")));
  auto degree_cyan = 1. / 2 - (FlagVector<1, 2>(flag("2 1  2 0  2")));
  auto degree_red = FlagVector<1, 2>(flag("2 1  3 0  2")) - 1. / 2;
  auto degree_magenta = FlagVector<1, 2>(flag("2 1  4 0  2")) - 1. / 2;
  problem.add_constraint(degree_blue);
  problem.add_constraint(degree_cyan);
  problem.add_constraint(degree_red);
  problem.add_constraint(degree_magenta);

  auto random_cut = FlagVector<0, 2>(flag("2 0  0 0  2")) - 8 * BOUND;
  problem.add_constraint(random_cut);

  vector<flag_coeff> vertices_less, vertices_more;
  if (ProblemConfig::instance().case_number & 1) {
    cerr << "Assuming there's more vertices of low degree." << endl;
    vertices_more.push_back(flag_coeff("1 0  1"));
    vertices_more.push_back(flag_coeff("1 0  2"));
    vertices_less.push_back(flag_coeff("1 0  3"));
    vertices_less.push_back(flag_coeff("1 0  4"));
  } else {
    cerr << "Assuming there's more vertices of high degree." << endl;
    vertices_more.push_back(flag_coeff("1 0  3"));
    vertices_more.push_back(flag_coeff("1 0  4"));
    vertices_less.push_back(flag_coeff("1 0  1"));
    vertices_less.push_back(flag_coeff("1 0  2"));
  }
  problem.add_constraint(FlagVector<0, 1>::from_vector(flag(), vertices_more) -
                         FlagVector<0, 1>::from_vector(flag(), vertices_less));
  auto cut_on_degrees = degree_cut(vertices_less, vertices_more);
  if (ProblemConfig::instance().case_number >> 1 & 1) {
    cerr << "Assuming that in the cut degree there is more edges in the half "
            "with vertices of the more prevalent degree."
         << endl;
    problem.add_constraint(cut_on_degrees.right_side -
                           cut_on_degrees.left_side);
    problem.add_constraint(cut_on_degrees.right_side -
                           cut_on_degrees.lower_bound);
  } else {
    cerr << "Assuming that in the cut degree there is more edges in the half "
            "with vertices of both kind."
         << endl;
    problem.add_constraint(cut_on_degrees.left_side -
                           cut_on_degrees.right_side);
    problem.add_constraint(cut_on_degrees.left_side -
                           cut_on_degrees.lower_bound);
  }

  flag_coeff blue_vertex_connected("2 1  1 0  2");
  flag_coeff blue_vertex_disconnected("2 1  1 0  1");
  flag_coeff cyan_vertex_connected("2 1  2 0  2");
  flag_coeff cyan_vertex_disconnected("2 1  2 0  1");

  auto cut_info_blue_vertex = prepare_cut_halves<1>(
      {blue_vertex_connected}, {}, {blue_vertex_disconnected}, BOUND);
  problem.add_constraint(cut_info_blue_vertex.left_side -
                         cut_info_blue_vertex.right_side);
  problem.add_constraint(cut_info_blue_vertex.left_side -
                         cut_info_blue_vertex.lower_bound);

  auto cut_info_cyan_vertex = prepare_cut_halves<1>(
      {cyan_vertex_connected}, {}, {cyan_vertex_disconnected}, BOUND);
  problem.add_constraint(cut_info_cyan_vertex.right_side -
                         cut_info_cyan_vertex.left_side);
  problem.add_constraint(cut_info_cyan_vertex.right_side -
                         cut_info_cyan_vertex.lower_bound);

  flag_coeff red_vertex_disconnected("2 1  3 0  1");
  flag_coeff red_vertex_connected_blue("2 1  3 1  2");
  flag_coeff red_vertex_connected_cyan("2 1  3 2  2");
  flag_coeff red_vertex_connected_red("2 1  3 3  2");
  flag_coeff red_vertex_connected_magenta("2 1  3 4  2");

  auto cut_on_red_vertex = prepare_cut_halves<1>(
      {red_vertex_disconnected},
      {red_vertex_connected_red, red_vertex_connected_magenta},
      {red_vertex_connected_blue, red_vertex_connected_cyan}, BOUND);
  // We don't need to care about the other case, because right side is a subset
  // of a triangle-free graph, which allows us to get a very strong bound even
  // with a simple application of Mantel's theorem.
  problem.add_constraint(cut_on_red_vertex.left_side -
                         cut_on_red_vertex.right_side);
  problem.add_constraint(cut_on_red_vertex.left_side -
                         cut_on_red_vertex.lower_bound);

  flag_coeff magenta_vertex_disconnected("2 1  4 0  1");
  flag_coeff magenta_vertex_connected_blue("2 1  4 1  2");
  flag_coeff magenta_vertex_connected_cyan("2 1  4 2  2");
  flag_coeff magenta_vertex_connected_red("2 1  4 3  2");
  flag_coeff magenta_vertex_connected_magenta("2 1  4 4  2");

  auto cut_on_magenta_vertex = prepare_cut_halves<1>(
      {magenta_vertex_disconnected},
      {magenta_vertex_connected_red, magenta_vertex_connected_magenta},
      {magenta_vertex_connected_blue, magenta_vertex_connected_cyan}, BOUND);
  // We don't need to care about the other case, because right side is a subset
  // of a triangle-free graph, which allows us to get a very strong bound even
  // with a simple application of Mantel's theorem.
  problem.add_constraint(cut_on_magenta_vertex.left_side -
                         cut_on_magenta_vertex.right_side);
  problem.add_constraint(cut_on_magenta_vertex.left_side -
                         cut_on_magenta_vertex.lower_bound);

  auto projected_edge_cut_on_red_vertex = projected_edge_cut(3, {3, 4});
  problem.add_constraint(projected_edge_cut_on_red_vertex.left_side -
                         projected_edge_cut_on_red_vertex.right_side);
  problem.add_constraint(projected_edge_cut_on_red_vertex.left_side -
                         projected_edge_cut_on_red_vertex.lower_bound);

  auto projected_edge_cut_on_magenta_vertex = projected_edge_cut(4, {3, 4});
  problem.add_constraint(projected_edge_cut_on_magenta_vertex.right_side -
                         projected_edge_cut_on_magenta_vertex.left_side);
  problem.add_constraint(projected_edge_cut_on_magenta_vertex.right_side -
                         projected_edge_cut_on_magenta_vertex.lower_bound);

  // flag_coeff rr_edge_left_only("3 2  3 3 0  2 2  1");
  // flag_coeff rr_edge_right_only("3 2  3 3 0  2 1  2");
  // flag_coeff rr_edge_neither("3 2  3 3 0  2 1  1", 0.5);
  // flag_coeff rr_edge_both("3 2  3 3 0  2 2  2");
  // auto projected_cut_on_rr_edge =
  //     projected_cut<2>({rr_edge_left_only, rr_edge_neither},
  //                   {rr_edge_right_only, rr_edge_neither}, {rr_edge_both});
  // problem.add_constraint(projected_cut_on_rr_edge.left_side -
  //                        projected_cut_on_rr_edge.right_side);
  // problem.add_constraint(projected_cut_on_rr_edge.left_side -
  //                        projected_cut_on_rr_edge.lower_bound);

  // flag_coeff mm_edge_left_only("3 2  4 4 0  2 2  1");
  // flag_coeff mm_edge_right_only("3 2  4 4 0  2 1  2");
  // flag_coeff mm_edge_neither("3 2  4 4 0  2 1  1", 0.5);
  // flag_coeff mm_edge_both("3 2  4 4 0  2 2  2");
  // auto projected_cut_on_mm_edge =
  //     projected_cut<2>({mm_edge_left_only, mm_edge_neither},
  //                   {mm_edge_right_only, mm_edge_neither}, {mm_edge_both});
  // problem.add_constraint(projected_cut_on_mm_edge.right_side -
  //                        projected_cut_on_mm_edge.left_side);
  // problem.add_constraint(projected_cut_on_mm_edge.right_side -
  //                        projected_cut_on_mm_edge.lower_bound);

  solve_sdp_for_problem(problem.get_constraints(), problem.get_objective());
}
