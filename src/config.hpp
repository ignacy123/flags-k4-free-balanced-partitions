#pragma once
#include "utils.hpp"
#include <sstream>
#include <string>

using namespace std;

#ifndef V
#define V                                                                      \
  4 // Maximum number of vertices of flags. If you know you need smaller, use
    // smaller - will be faster and less memory consuming
#endif

#ifndef G_COLORED_EDGES
#define G_COLORED_EDGES 2
#endif

#define G_NOT_ALL_FLAGS_USED // Some extra warnings if there are unexpected
                             // flags. This should be enabled

#define G_USING_ZERO_AS_ANY_COLOR // If something has color 0, it is considered
                                  // to have anycolor or being uncolored.
                                  // Usefull when trying to find extensions of
                                  // some graph

#define G_BE_BRAVE_AND_IGNORE_SAFETY_ASSERTS // Will be slower, but if not
                                             // defined. But if you don't trust
                                             // the program, comment it :-)

#define USE_REDUCED_TYPES // 4 4 1:123  and 4 4 1:234 are considered redundant
                          // and only one is used. This you want to have on. Not
                          // having it is stupid unless you have ordered
                          // vertices

#ifdef G_COLORED_VERTICES
#define COLORS_VERTICES (G_COLORED_VERTICES + 1)
#endif
#ifdef G_COLORED_EDGES
#define COLORS_EDGES (G_COLORED_EDGES + 1)
#endif

#ifndef G_BLOW_UP_COLOR_EDGES
#define G_BLOW_UP_COLOR_EDGES 1
#endif

#ifndef G_PRECISION
#define G_PRECISION 16
#endif

struct ProblemConfig {
  string data_directory = "k4-free";
  string sdp_directory = "sum-bicol";
  string csdp_binary = "csdp-no-accelerate";
  string sdpa_binary = "sdpa";

  int largest_flag_to_show = V;

  bool use_product_linear_constraints = false;
  bool use_product_square_constraints = false;
  bool upper_bound = false;
  bool verbose_output = false;
  bool use_csdp = true;
  bool use_sdpa = false;
  bool skip_empty_flags_when_processing = true;

  static ProblemConfig &instance() {
    static ProblemConfig inst;
    return inst;
  }

  string get_sdp_problem_name() {
    stringstream filename;
    create_dir(ProblemConfig::instance().sdp_directory);
    create_dir("SDP");
    create_dir("SDP/" + ProblemConfig::instance().sdp_directory);
    filename << "SDP/" << ProblemConfig::instance().sdp_directory
             << "/problem.dat-s";
    return filename.str();
  }

  string get_sdp_solution_name() {
    stringstream filename;
    filename << get_sdp_problem_name() << ".result";
    return filename.str();
  }

  string get_sdp_extremal_construction_name() {
    stringstream filename;
    filename << get_sdp_problem_name() << ".extremal.txt";
    return filename.str();
  }
};
