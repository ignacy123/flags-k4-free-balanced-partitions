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
                                   bool use_sdp_temp = false,
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
                ostream &ostr = cout, int verbose_output = 1,
                bool use_sdp_temp = false, bool sdp_temp_up_to_date = true) {
  cerr << "Computing specific part of the SDP..." << endl;
  print_CSDP_specific_part(constraints, objective, upper_bound, ostr,
                           verbose_output);

  if (!use_sdp_temp || sdp_temp_up_to_date == false) {
    cerr << "Computing flag products part of the SDP..." << endl;
    print_CSDP_flag_products_part(constraints, upper_bound, ostr,
                                  verbose_output, use_sdp_temp);
    return;
  }

  cerr << "Getting flag products for SDP..." << endl;

  stringstream filename;

  filename << filename_prefix() << "__n" << V << "_sdp_products.txt";

  // The sdp_product_file may not contain all products if the computation
  // was accidentally interrupted. We try to go from we stopped...
  // First find were we stopped

  int lowest_not_computed = 0;

  ifstream sdp_product_file;
  sdp_product_file.open(filename.str().c_str(), ifstream::in);
  while (sdp_product_file.good()) {
    int id, type, f1, f2, d;
    sdp_product_file >> id >> type >> f1 >> f2 >> d;
    if (!sdp_product_file.good())
      break; // if the read failed, do not use it...
    if (id > lowest_not_computed) {
      lowest_not_computed = id; // note that in file the ID is +1
    }
  }
  sdp_product_file.close();

  if (lowest_not_computed < (int)get_unlabeled_flags(V).size()) {
    if (lowest_not_computed != 0) {
      cerr << "Warning: file " << filename.str()
           << " contains only products up to " << lowest_not_computed
           << " out of " << (int)get_unlabeled_flags(V).size() << endl;
      cerr << "          We assume it happened by accident so the rows "
              "starting with "
           << lowest_not_computed << " will be deleted and recomputed" << endl;
      stringstream mv_command;
      mv_command << "rm -f " << filename.str() << "~ ;  mv " << filename.str()
                 << " " << filename.str() << "~";
      cerr << "          Executing " << mv_command.str() << endl;
      if (system(mv_command.str().c_str()) != 0) {
        cerr << " Execution failed" << endl;
        exit(1);
      }
      // We remove the last line since it may be incoplete and not catched by
      // the grep
      stringstream grep_command;
      grep_command << "sed \\$d " << filename.str() << "~ |   grep -v '^"
                   << lowest_not_computed << " ' " << " >" << filename.str();
      cerr << "          Executing " << grep_command.str() << endl;
      if (system(grep_command.str().c_str()) != 0) {
        cerr << " Execution failed" << endl;
        exit(1);
      }
      lowest_not_computed--;
    }

    ofstream sdp_product_file_out;
    sdp_product_file_out.open(filename.str().c_str(),
                              ofstream::out | ofstream::app);
    if (!sdp_product_file_out.good()) {
      cerr << "Failed creating file with sdp products " << filename.str()
           << endl;
      exit(1);
    }

    cerr << "Computing flag products of the SDP to file " << filename.str()
         << endl;
    print_CSDP_flag_products_part(constraints, upper_bound,
                                  sdp_product_file_out, verbose_output,
                                  use_sdp_temp, lowest_not_computed);
    sdp_product_file_out.close();
  }

  sdp_product_file.open(filename.str().c_str(), ifstream::in);
  if (!sdp_product_file.good()) {
    cerr << "Failed obtaining file with sdp products " << filename.str()
         << endl;
    exit(1);
  }

  cerr << "Copying flag products of the SDP from " << filename.str() << endl;
  ostr << sdp_product_file.rdbuf();
  sdp_product_file.close();
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

    if (ProblemConfig::instance().use_sdpa && sdpa_available) {
      log_file += ".sdpa.log";
      system_command << "'" << ProblemConfig::instance().sdpa_binary
                     << "' -ds '" << output_file << "' -o '" << result_file
                     << "'  2>&1  | tee  '" << log_file << "'";
    } else if (ProblemConfig::instance().use_csdp && csdp_available) {
      log_file += ".csdp.log";
      system_command << "'" << ProblemConfig::instance().csdp_binary << "'  '"
                     << output_file << "'  '" << result_file
                     << "'  2>&1  | tee  '" << log_file << "'";
    } else {
      cerr << "ERROR: No semidefinie programming solver found." << endl;
      return 1;
    }
    cerr << "Executing: " << system_command.str() << endl;
    csdp_return_value = system(system_command.str().c_str());

    // execlp(csdp_binary.c_str(),"csdp",output_file.c_str(),result_file.c_str(),(char
    // *)NULL);
  }

  if (csdp_return_value != 0) {
    cout << "Command processor doesn't exists or the system call failed.";

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
  }
  // execlp("sdpa","sdpa","-ds",output_file.c_str(),"-o",result_file.c_str(),(char
  // *)NULL); system("csdp",output_file.c_str());

  cerr << "Done solving SDP." << endl;
  cerr << "Trying to process solution..." << endl;

  if (!ProblemConfig::instance().use_csdp) {
    cerr << "Can only process CSDP solutions..." << endl;
    exit(0);
  }

  process_sdp_solution(ProblemConfig::instance());

  return 0;
}
