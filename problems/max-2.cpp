#include "balanced-bipartition-utils.hpp"
#include "config.hpp"
#include "flag-calculator.hpp"
#include "flag.hpp"
#include "option-parser.hpp"
#include "problem.hpp"
#include "solver.hpp"
#include <cassert>
#include <string>
#include <vector>

double BOUND = 0.067;

CutInfo<1, 5> cut_vertex_high_degree(int color, vector<int> colors_low_degree,
                                     vector<int> colors_high_degree) {
  flag_coeff disconnected("2 1  " + to_string(color) + " 0  1");

  vector<flag_coeff> connected_low_degree;
  vector<flag_coeff> connected_high_degree;
  for (int low_degree_color : colors_low_degree) {
    connected_low_degree.push_back(flag_coeff("2 1  " + to_string(color) + " " +
                                              to_string(low_degree_color) +
                                              "  2"));
  }
  for (int high_degree_color : colors_high_degree) {
    connected_high_degree.push_back(
        flag_coeff("2 1  " + to_string(color) + " " +
                   to_string(high_degree_color) + "  2"));
  }
  return prepare_cut_halves<1>({disconnected}, connected_high_degree,
                               connected_low_degree, BOUND);
}

CutInfo<1, 5> cut_vertex_low_degree(int color) {
  flag_coeff connected("2 1  " + to_string(color) + " 0  2");
  flag_coeff disconnected("2 1  " + to_string(color) + " 0  1");

  return prepare_cut_halves<1>({connected}, {}, {disconnected}, BOUND);
}

CutInfo<1, 5> projected_cut_high_degree_edge(int color,
                                             vector<int> other_colors) {

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

CutInfo<1, 6> projected_cut_low_degree_edge(int color,
                                            vector<int> other_colors) {

  FlagVector<1, 6> left_side(flag("1 1  " + to_string(color)));
  FlagVector<1, 6> right_side(flag("1 1  " + to_string(color)));
  FlagVector<1, 6> lower_bound(flag("1 1  " + to_string(color)));
  // Heavily relying on the linearity of averaging operator.
  for (int other_color : other_colors) {
    flag_coeff left_only("3 2  " + to_string(color) + " " +
                         to_string(other_color) + " 0  2 2  1");
    flag_coeff right_only("3 2  " + to_string(color) + " " +
                          to_string(other_color) + " 0  2 1  2");
    flag_coeff neither("3 2  " + to_string(color) + " " +
                       to_string(other_color) + " 0  2 1  1");
    flag_coeff both("3 2  " + to_string(color) + " " + to_string(other_color) +
                    " 0  2 2  2");
    auto cut_halves = prepare_cut_halves<2>({left_only, both}, {right_only},
                                            {neither}, BOUND);
    left_side += cut_halves.left_side.project<1, 6>();
    right_side += cut_halves.right_side.project<1, 6>();
    lower_bound += cut_halves.lower_bound.project<1, 6>();
  }

  return CutInfo<1, 6>{left_side, right_side, lower_bound};
}

CutInfo<0, 4> degree_cut(const vector<int> &vertices_less,
                         const vector<int> &vertices_more) {
  vector<flag_coeff> left;
  vector<flag_coeff> random;
  for (int color : vertices_less) {
    left.push_back(flag_coeff("1 0  " + to_string(color)));
  }
  for (int color : vertices_more) {
    random.push_back(flag_coeff("1 0  " + to_string(color)));
  }
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

FlagVector<0, 1> vertex_color_constraint(vector<int> vertices_more) {
  FlagVector<0, 1> result;
  for (int color : vertices_more) {
    result += flag_coeff("1 0  " + to_string(color));
  }
  return result - 1. / 2;
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

  vector<int> low_degree_vertices{1, 2};
  vector<int> high_degree_vertices{3};
  if (ProblemConfig::instance().case_number == 1) {
    high_degree_vertices.push_back(4);
  }
  for (int low_degree_vertex : low_degree_vertices) {
    problem.add_constraint(
        1. / 2 - FlagVector<1, 2>(
                     flag("2 1  " + to_string(low_degree_vertex) + " 0  2")));
  }
  for (int high_degree_vertex : high_degree_vertices) {
    problem.add_constraint(
        FlagVector<1, 2>(
            flag("2 1  " + to_string(high_degree_vertex) + " 0  2")) -
        1. / 2);
  }

  auto random_cut = FlagVector<0, 2>(flag("2 0  0 0  2")) - 8 * BOUND;
  problem.add_constraint(random_cut);

  FlagVector<0, 2> high_degree_edges;
  for (int i = 0; i < high_degree_vertices.size(); i++) {
    for (int j = i; j < high_degree_vertices.size(); j++) {
      high_degree_edges +=
          FlagVector<0, 2>("2 0  " + to_string(high_degree_vertices[i]) + " " +
                           to_string(high_degree_vertices[j]) + " 2");
    }
  }
  // Consider the total number of red vertices.
  // In other words, define the cutoff by with respect to the density of edges
  // between vertices of high degree.
  double high_degree_edges_cutoff = 0.1;
  if (ProblemConfig::instance().case_number == 1) {
    problem.add_constraint(high_degree_edges - high_degree_edges_cutoff);
  }
  if (ProblemConfig::instance().case_number == 2) {
    problem.add_constraint(high_degree_edges_cutoff - high_degree_edges);
  }

  vector<int> vertices_more, vertices_less;
  if (ProblemConfig::instance().case_number == 1 or
      ProblemConfig::instance().case_number == 2 or
      ProblemConfig::instance().case_number == 4) {
    cerr << "Assuming there's more vertices of low degree." << endl;
    vertices_more = low_degree_vertices;
    vertices_less = high_degree_vertices;
  }
  if (ProblemConfig::instance().case_number == 3 or
      ProblemConfig::instance().case_number == 5) {
    cerr << "Assuming there's more vertices of high degree." << endl;
    vertices_more = high_degree_vertices;
    vertices_less = low_degree_vertices;
  }
  problem.add_constraint(vertex_color_constraint(vertices_more));

  auto cut_on_degrees = degree_cut(vertices_less, vertices_more);
  if (ProblemConfig::instance().case_number == 4 or
      ProblemConfig::instance().case_number == 5) {
    cerr << "Assuming that in the cut degree there is more edges in the half "
            "with vertices of the more prevalent degree."
         << endl;
    problem.add_constraint(cut_on_degrees.right_side -
                           cut_on_degrees.left_side);
    problem.add_constraint(cut_on_degrees.right_side -
                           cut_on_degrees.lower_bound);
  }
  if (ProblemConfig::instance().case_number == 1 or
      ProblemConfig::instance().case_number == 2 or
      ProblemConfig::instance().case_number == 3) {
    cerr << "Assuming that in the cut degree there is more edges in the half "
            "with vertices of both kind."
         << endl;
    problem.add_constraint(cut_on_degrees.left_side -
                           cut_on_degrees.right_side);
    problem.add_constraint(cut_on_degrees.left_side -
                           cut_on_degrees.lower_bound);
  }

  if (ProblemConfig::instance().case_number == 2) {
    auto projected_edge_cut_on_blue_vertex =
        projected_cut_low_degree_edge(1, low_degree_vertices);
    problem.add_constraint(projected_edge_cut_on_blue_vertex.left_side -
                           projected_edge_cut_on_blue_vertex.right_side);
    problem.add_constraint(projected_edge_cut_on_blue_vertex.left_side -
                           projected_edge_cut_on_blue_vertex.lower_bound);

    auto projected_edge_cut_on_cyan_vertex =
        projected_cut_low_degree_edge(2, low_degree_vertices);
    problem.add_constraint(projected_edge_cut_on_cyan_vertex.right_side -
                           projected_edge_cut_on_cyan_vertex.left_side);
    problem.add_constraint(projected_edge_cut_on_cyan_vertex.right_side -
                           projected_edge_cut_on_cyan_vertex.lower_bound);
  } else {
    auto blue_vertex_cut = cut_vertex_low_degree(1);
    problem.add_constraint(blue_vertex_cut.left_side -
                           blue_vertex_cut.right_side);
    problem.add_constraint(blue_vertex_cut.left_side -
                           blue_vertex_cut.lower_bound);

    auto cyan_vertex_cut = cut_vertex_low_degree(2);
    problem.add_constraint(cyan_vertex_cut.right_side -
                           cyan_vertex_cut.left_side);
    problem.add_constraint(cyan_vertex_cut.right_side -
                           cyan_vertex_cut.lower_bound);
  }

  auto red_vertex_cut =
      cut_vertex_high_degree(3, low_degree_vertices, high_degree_vertices);
  // We don't need to care about the other case, because right side is a
  // subset of a triangle-free graph, which allows us to get a very strong
  // bound even with a simple application of Mantel's theorem.
  problem.add_constraint(red_vertex_cut.left_side - red_vertex_cut.right_side);
  problem.add_constraint(red_vertex_cut.left_side - red_vertex_cut.lower_bound);

  if (ProblemConfig::instance().case_number == 1) {
    auto magenta_vertex_cut =
        cut_vertex_high_degree(4, low_degree_vertices, high_degree_vertices);
    // We don't need to care about the other case, because right side is a
    // subset of a triangle-free graph, which allows us to get a very strong
    // bound even with a simple application of Mantel's theorem.
    problem.add_constraint(magenta_vertex_cut.left_side -
                           magenta_vertex_cut.right_side);
    problem.add_constraint(magenta_vertex_cut.left_side -
                           magenta_vertex_cut.lower_bound);

    auto projected_edge_cut_on_red_vertex =
        projected_cut_high_degree_edge(3, high_degree_vertices);
    problem.add_constraint(projected_edge_cut_on_red_vertex.left_side -
                           projected_edge_cut_on_red_vertex.right_side);
    problem.add_constraint(projected_edge_cut_on_red_vertex.left_side -
                           projected_edge_cut_on_red_vertex.lower_bound);

    auto projected_edge_cut_on_magenta_vertex =
        projected_cut_high_degree_edge(4, high_degree_vertices);
    problem.add_constraint(projected_edge_cut_on_magenta_vertex.right_side -
                           projected_edge_cut_on_magenta_vertex.left_side);
    problem.add_constraint(projected_edge_cut_on_magenta_vertex.right_side -
                           projected_edge_cut_on_magenta_vertex.lower_bound);
  }

  solve_sdp_for_problem(problem.get_constraints(), problem.get_objective());
}
