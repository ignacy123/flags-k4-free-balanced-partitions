#include "flag-calculator.hpp"
#include "flag.hpp"
#include "option-parser.hpp"
#include "problem.hpp"
#include "solver.hpp"
#include <cassert>
#include <vector>

double BOUND = 5. / 72;

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
      FlagVector<root_size, root_size + 2>(type, 2 * BOUND)};
}

template <int root_size, int final_size> struct CutInfo {
  FlagVector<root_size, final_size> left_side;
  FlagVector<root_size, final_size> right_side;
  FlagVector<root_size, final_size> lower_bound;
};

template <int root_size>
CutInfo<root_size, root_size + 4> regular_cut(vector<flag> left_side,
                                              vector<flag> right_side,
                                              vector<flag> random) {
  auto cut_components =
      prepare_cut_components<root_size>(left_side, right_side, random);
  auto left = cut_components.fixed_left * cut_components.denominator *
              cut_components.denominator;
  left += cut_components.space_left * cut_components.one_left *
          cut_components.denominator;
  left += cut_components.neither * cut_components.space_left *
          cut_components.space_left;
  auto right = cut_components.fixed_right * cut_components.denominator *
               cut_components.denominator;
  right += cut_components.space_right * cut_components.one_right *
           cut_components.denominator;
  right += cut_components.neither * cut_components.space_right *
           cut_components.space_right;
  return CutInfo<root_size, root_size + 4>{left, right,
                                           cut_components.lower_bound *
                                               cut_components.denominator *
                                               cut_components.denominator};
}

CutInfo<0, 4> degree_cut(vector<flag> left, vector<flag> random) {
  FlagVector<0, 1> sum;
  FlagVector<0, 1> ones(flag(), 1.);
  for (auto flag : left) {
    assert(flag.m_vertices == 1);
    sum += flag;
  }
  for (auto flag : random) {
    assert(flag.m_vertices == 1);
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
  int scaling_factor = 2 * 2 * 2 * 2 * 3 * 3 * 3 * 3 * 5;
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
  // auto degree_magenta = FlagVector<1, 2>(flag("2 1  4 0  2")) - 1. / 2;
  problem.add_constraint(degree_blue);
  problem.add_constraint(degree_cyan);
  problem.add_constraint(degree_red);
  // problem.add_constraint(degree_magenta);

  vector<flag> vertices_less, vertices_more;
  if (ProblemConfig::instance().case_number & 1) {
    vertices_more.push_back(flag("1 0  1"));
    vertices_more.push_back(flag("1 0  2"));
    vertices_less.push_back(flag("1 0  3"));
    // vertices_less.push_back(flag("1 0  4"));
  } else {
    vertices_more.push_back(flag("1 0  3"));
    // vertices_more.push_back(flag("1 0  4"));
    vertices_less.push_back(flag("1 0  1"));
    vertices_less.push_back(flag("1 0  2"));
  }
  problem.add_constraint(FlagVector<0, 1>::from_vector(flag(), vertices_more) -
                         FlagVector<0, 1>::from_vector(flag(), vertices_less));
  auto cut_on_degrees = degree_cut(vertices_less, vertices_more);
  if (ProblemConfig::instance().case_number >> 1 & 1) {
    problem.add_constraint(cut_on_degrees.right_side -
                           cut_on_degrees.left_side);
    problem.add_constraint(cut_on_degrees.right_side -
                           cut_on_degrees.lower_bound);
  } else {
    problem.add_constraint(cut_on_degrees.left_side -
                           cut_on_degrees.right_side);
    problem.add_constraint(cut_on_degrees.left_side -
                           cut_on_degrees.lower_bound);
  }

  flag blue_vertex_connected("2 1  1 0  2");
  flag blue_vertex_disconnected("2 1  1 0  1");
  flag cyan_vertex_connected("2 1  2 0  2");
  flag cyan_vertex_disconnected("2 1  2 0  1");

  auto cut_info_blue_vertex =
      regular_cut<1>({blue_vertex_connected}, {}, {blue_vertex_disconnected});
  problem.add_constraint(cut_info_blue_vertex.left_side -
                         cut_info_blue_vertex.right_side);
  problem.add_constraint(cut_info_blue_vertex.left_side -
                         cut_info_blue_vertex.lower_bound);

  auto cut_info_cyan_vertex =
      regular_cut<1>({cyan_vertex_connected}, {}, {cyan_vertex_disconnected});
  problem.add_constraint(cut_info_cyan_vertex.right_side -
                         cut_info_cyan_vertex.left_side);
  problem.add_constraint(cut_info_cyan_vertex.right_side -
                         cut_info_cyan_vertex.lower_bound);

  flag red_vertex_disconnected("2 1  3 0  1");
  flag red_vertex_connected_blue("2 1  3 1  2");
  flag red_vertex_connected_cyan("2 1  3 2  2");
  flag red_vertex_connected_red("2 1  3 3  2");
  // flag red_vertex_connected_magenta("2 1  3 4  2");

  auto cut_info_red_vertex =
      regular_cut<1>({red_vertex_disconnected}, {red_vertex_connected_red},
                     {red_vertex_connected_blue, red_vertex_connected_cyan});
  // We don't need to care about the other case, because right side is a subset
  // of a triangle-free graph, which allows us to get a very strong bound even
  // with a simple application of Mantel's theorem.
  problem.add_constraint(cut_info_red_vertex.left_side -
                         cut_info_red_vertex.right_side);
  problem.add_constraint(cut_info_red_vertex.left_side -
                         cut_info_red_vertex.lower_bound);

  solve_sdp_for_problem(problem.get_constraints(), problem.get_objective());
}
