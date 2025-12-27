#include <assert.h>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <errno.h>
#include <fstream>
#include <iostream>
#include <istream>
#include <limits.h>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#include <vector>

#include "cache.hpp"
#include "config.hpp"
#include "flag-calculator.hpp"
#include "flag-generator.hpp"
#include "flag.hpp"
#include "option-parser.hpp"
#include "solution-processor.hpp"
#include "utils.hpp"

#ifdef _USING_OMP_
#include <omp.h>
#endif

using namespace std;

int print_CSDP_simple_linear_constraints(
    ostream &ostr, const vector<FlagVector<0, V>> &constraints, int i,
    int matrixID, int blockID, bool print_blocks, bool print_products,
    int verbose_output) {

  int constraint_id = 0;
  for (int j = 0; j < (int)constraints.size(); j++) {
    constraint_id++;

    if (print_products) {
      double d = constraints[j].get_coefficients()[i];
      if (d != 0) {
        ostr.precision(G_PRECISION);
        ostr << matrixID << " " << blockID << " " << constraint_id << " "
             << constraint_id << " " << smart_round(d) << endl;
      }
    }
  }

  if (print_blocks) {
    ostr << -constraint_id << " ";
    if (verbose_output) {
      cerr << "Simple linear are used for " << constraint_id << "/"
           << constraints.size() << " constraints" << endl;
    }
  }

  return constraint_id;
}

int print_CSDP_constraints_header(const vector<FlagVector<0, V>> &constraints,
                                  ostream &ostr = cout,
                                  int verbose_output = 1) {
  int blocks = 0;

  // next we have linear constraints blocks

  if (constraints.size() != 0) {
    // Simple linear constraints
    int simple_constraints = print_CSDP_simple_linear_constraints(
        ostr, constraints, 0, 0, 0, true, false, verbose_output);
    if (simple_constraints != 0) {
      // block_ID_simple_linear_constraints = blocks+1;
      blocks++;
    }
  }

  return blocks;
}

int print_CSDP_constraints_blocks(ostream &ostr,
                                  const vector<FlagVector<0, V>> &constraints,
                                  int i, int matrixID, int current_csdp_block,
                                  int verbose_output = 1) {
  if (print_CSDP_simple_linear_constraints(ostr, constraints, i, matrixID,
                                           current_csdp_block, false, true,
                                           verbose_output) != 0) {
    current_csdp_block++;
  }

  return current_csdp_block;
}

int print_CSDP_additional_blocks_header(
    const vector<FlagVector<0, V>> &constraints, ostream &ostr,
    int verbose_output = 1) {
  int blocks = 0;

  if (verbose_output) {
    cerr << "Added " << blocks << " additional blocks." << endl;
  }
  return blocks;
}

// According to http://plato.asu.edu/ftp/sdpa_format.txt
//
// Solved program
//   max C.X
//   subject to A_1.X = b_1
//              A_2.X = b_2
//              A_mdim.X = b_mdim
//   where X is positive semidefinite matrix
//
// In our UPPER BOUND application
//
//      min    t
//      s.t.  D_i + [A.X]_i <= t
//
// We add a slack variable s_i to make the <= to = and move D_i (density -
// constant to the other side) and change to maximization:
//
//      max   -t
//      s.t.  [A.X]_i + s_i - t = -D_i
//
//   C.X is just entry 0 2 1 1 -1.0 which means that in the second block,
//   the first variable is with -1. So we are maximizing -t which corresponds to
//   minimizing t
//
// Now LOWER BOUND application
//
//      max    t
//      s.t.  D_i - [A.X]_i  >= t
//
// We add a slack variable s_i to make the >= to = and move D_i (density -
// constant to the other side) and change to maximization:
//
//      max   t
//      s.t.  [A.X]_i + s_i + t = D_i
//
// Finally, Additional constraints applications:
// Say we want to add  a*G >= c. Then we could add
//  k(a*G(sampled in H_i) -c) >= 0
//
// For lower bound we get (max t)
//     [A.X]_i + k*(a*G_i-c) + s_i + t = D_i
//
// For upper bound we get (min t)
//     [A.X]_i + k*(a*G_i-c) + s_i -t = D_i
//
//
// Notice that by definition, t >= 0
// We can write t = t1 - t2 instead and make t a free variable.

void print_CSDP_specific_part(const vector<FlagVector<0, V>> &constraints,
                              const FlagVector<0, V> objective,
                              bool upper_bound = true, ostream &ostr = cout,
                              int verbose_output = 1) {
  int mdim = 0;
  mdim = (int)get_unlabeled_flags(V).size();

  // hack
  int extra_onstraints = 0;
  // extra_onstraints +=
  // print_CSDP_product_linear_constraints_extra_lines_in_csdp(ostr, Kn, 0,
  // false, verbose_output)

  ostr << mdim + extra_onstraints
       << endl; //<<" =mdim" << endl; // number of constraint matrices

  stringstream ssblocks; // String stream for all blocks in the matrix
  int blocks = 0;        // number of blocks in ssblocks

  // First blocks are flags products:
  for (int f = 0; f < (int)get_flags().size(); f++)
    ssblocks << get_flags()[f].size() << " ";
  blocks = (int)get_flags().size();

  // next we have one diagonal block of slack variables
  ssblocks << -1 - mdim
           << " "; // get_number_of_slacks() but we don't go for it here...
  int objective_index = mdim + 1; // Objective variables come after slacks

  blocks++;

  int first_constraint_block = blocks + 1;

  blocks +=
      print_CSDP_constraints_header(constraints, ssblocks, verbose_output);

  int additional_blocks_offset = blocks + 1;
  blocks += print_CSDP_additional_blocks_header(constraints, ssblocks,
                                                verbose_output);

  ostr << blocks << endl;
  ostr << ssblocks.str() << endl;

  for (unsigned int i = 0; i < get_unlabeled_flags(V).size(); i++) {
    ostr.precision(G_PRECISION);
    if (upper_bound)
      ostr << -1 * smart_round(objective.get_coefficients()[i])
           << " "; // outputing b_i
    else
      ostr << smart_round(objective.get_coefficients()[i])
           << " "; // outputing b_i
  }
  ostr << endl;

  if (upper_bound)
    ostr << "0 " << (int)get_flags().size() + 1 << " " << objective_index << " "
         << objective_index << " -1" << endl; // Objective function
  else
    ostr << "0 " << (int)get_flags().size() + 1 << " " << objective_index << " "
         << objective_index << " 1" << endl; // Objective function

  // printing of subflags
  for (int i = 0; i < (int)get_unlabeled_flags(V).size(); i++) {
    // cout << "flag: "; c4free_subrgaphs[e][i].print();

    stringstream ss;

    int matrixID = i + 1;

    // variable for slack
    ss << matrixID << " " << (int)get_flags().size() + 1 << " " << matrixID
       << " " << matrixID << " 1" << endl;

    // variable for objective function
    if (upper_bound)
      ss << matrixID << " " << (int)get_flags().size() + 1 << " "
         << objective_index << " " << objective_index << " -1" << endl;
    else
      ss << matrixID << " " << (int)get_flags().size() + 1 << " "
         << objective_index << " " << objective_index << " 1" << endl;

    // int next_block =
    print_CSDP_constraints_blocks(ss, constraints, i, matrixID,
                                  first_constraint_block, verbose_output);
    ostr << ss.str();
  }

  ostr.flush();
}

void print_CSDP_flag_products_part(const vector<FlagVector<0, V>> &constraints,
                                   bool upper_bound = true,
                                   ostream &ostr = cout, int verbose_output = 1,
                                   int lowest_not_computed = 0) {

  for (int i = lowest_not_computed; i < (int)get_unlabeled_flags(V).size();
       i++) {

    stringstream ss;

    int matrixID = i + 1;
    // We go through all possible type sizes
    const flag &f = get_unlabeled_flags(V)[i];
    for (int typesize = f.m_vertices - 2; typesize >= 0; typesize -= 2) {
      const vector<FlagProductPartialResult> &cache =
          get_flag_product_cache(V, i, typesize);

      if (verbose_output) {
        cerr << "Flag product cache size: " << cache.size() << endl;
      }

      for (auto partial_result : cache) {
        ss << matrixID << " " << partial_result.product << " "
           << partial_result.count << endl;
      }
    }
    ostr << ss.str();
  }
}

void print_CSDP(const vector<FlagVector<0, V>> &constraints,
                const FlagVector<0, V> &objective, bool upper_bound = true,
                ostream &ostr = cout, int verbose_output = 1) {
  cerr << "Computing specific part of the SDP..." << endl;
  print_CSDP_specific_part(constraints, objective, upper_bound, ostr,
                           verbose_output);

  cerr << "Computing flag products part of the SDP..." << endl;
  print_CSDP_flag_products_part(constraints, upper_bound, ostr, verbose_output);
}

int solve_sdp_for_problem(const vector<FlagVector<0, V>> &constraints,
                          const FlagVector<0, V> &objective) {
  generate_unlabeled_flags_of_size(V, false, false, false,
                                   ProblemConfig::instance().verbose_output);

  // labeled flags
  if (!load_labeled_flags_from_file(V, true)) {

    cerr << "Generating labeled flags..." << endl;

    // Generating all labeled flags:
    int last_types = (int)get_flags().size();
    //// TODO: Make this parallel
    for (int i = 1; i <= V / 2; i++) {
      cerr << "Getting labeled flags of size " << V - i << ":" << V - 2 * i
           << endl;
      generate_labeled_flags(V - i, V - 2 * i,
                             ProblemConfig::instance().verbose_output, false);
      int gain = (int)get_flags().size() - last_types;
      cerr << "Got " << gain << " types ";
      for (int j = last_types; j < (int)get_flags().size(); j++)
        cerr << get_flags()[j].size() << " ";
      cerr << endl;
      last_types = (int)get_flags().size();
    }
    dump_labeled_flags(V);
  }
  cerr << "Labeled flag have " << (int)get_flags().size() << " types. Counts: ";
  for (int i = 0; i < (int)get_flags().size(); i++)
    cerr << " " << get_flags()[i].size();
  cerr << endl;
  if (get_flags().size() == 0) {
    cerr << "WARNING: No labeled flags. Computations might not work." << endl;
  }

  string output_file = ProblemConfig::instance().get_sdp_problem_name();
  string result_file = ProblemConfig::instance().get_sdp_solution_name();

  ofstream outfile;
  outfile.open(output_file.c_str(), ofstream::out);
  if (!outfile.good()) {
    cerr << "Failed opening file " << output_file << endl;
    return -1;
  }
  cerr << "Generating SDP program to " << output_file << endl;
  print_CSDP(constraints, objective, ProblemConfig::instance().upper_bound,
             outfile, ProblemConfig::instance().verbose_output);
  outfile.close();

  cerr << "Trying hack by setting MKL_DEBUG_CPU_TYPE=5" << endl;
  setenv("MKL_DEBUG_CPU_TYPE", "5", 1);

  int csdp_return_value = 1;
  if (system(NULL)) {
    // cout << "Command processor exists";
    // assert(use_sdpa == false);
    string command_test;

    cerr << "Solvers found: ";

    bool sdpa_available = false;
    command_test =
        "which  " + ProblemConfig::instance().sdpa_binary + " > /dev/null 2>&1";
    if (system(command_test.c_str())) {
      // Command doesn't exist...
    } else {
      // Command does exist, do something with it...
      sdpa_available = true;
      cerr << " sdpa ";
    }

    bool csdp_available = false;
    command_test =
        "which " + ProblemConfig::instance().csdp_binary + "  > /dev/null 2>&1";
    if (system(command_test.c_str())) {
      // Command doesn't exist...
    } else {
      // Command does exist, do something with it...
      csdp_available = true;
      cerr << " csdp ";
    }
    cerr << endl;

    stringstream system_command;
    string log_file = output_file;
  }

  pid_t csdp_PID = 0;
  csdp_PID = fork();

  if (csdp_PID == 0) {
    // clean memory! TODO

    if (!ProblemConfig::instance().use_csdp) {
      execlp("sdpa", "sdpa", "-ds", output_file.c_str(), "-o",
             result_file.c_str(), (char *)NULL);
    } else {
      execlp(ProblemConfig::instance().csdp_binary.c_str(), "csdp",
             output_file.c_str(), result_file.c_str(), (char *)NULL);
    }

    cerr << "Strugglig with executing: "
         << ProblemConfig::instance().csdp_binary.c_str() << endl;
    cerr << "Something went wrong: " << strerror(errno) << endl;
  } else {
    int csdo_status;
    pid_t tpid = wait(&csdo_status);
    if (tpid != csdp_PID) {
      cerr << "Something went wrong!" << endl;
    }
  }

  cerr << "Done solving SDP." << endl;
  cerr << "Trying to process solution..." << endl;

  if (!ProblemConfig::instance().use_csdp) {
    cerr << "Can only process CSDP solutions..." << endl;
    exit(0);
  }

  process_sdp_solution(ProblemConfig::instance());

  return 0;
}
