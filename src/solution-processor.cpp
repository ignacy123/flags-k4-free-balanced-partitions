#include "solution-processor.hpp"
#include "config.hpp"
#include "flag-calculator.hpp"
#include "flag-generator.hpp"
#include "flag.hpp"
#include "option-parser.hpp"
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>

void process_csdp_solution(istream &ist, ostream &ost, int largest_flag_to_show,
                           bool skip_if_empty) {

  double sum = 0;
  double sum_print = 0;
  double count_printed = 0;

  // TODO: Make sure this gets called automatically when we try to access the
  // unlabeled flags!
  generate_unlabeled_flags_of_size(V);
  FlagVector<0, V> densities;
  for (int i = 0; i < (int)get_unlabeled_flags(V).size(); i++) {
    double density;
    ist >> density;

    sum += density;

    count_printed++;
    sum_print += density;
    flag_and_coefficient fc;
    fc.g = get_unlabeled_flags(V)[i];
    fc.coefficient = density;
    densities += fc;
  }

  for (int j = largest_flag_to_show; j > 0; j--) {
    ost << "Vector of extremal densities on " << j << " vertices:" << endl;
    for (int i = 0; i < (int)get_unlabeled_flags(j).size(); i++) {
      flag_and_coefficient fc;
      fc.g = get_unlabeled_flags(j)[i];
      fc.coefficient = densities.density(j, i);
      fc.print(ost, skip_if_empty, to_string(i));
    }
    ost << endl;
  }
  ost << "Printed: " << count_printed << "/" << get_unlabeled_flags(V).size()
      << "  Sum of coefficients: " << sum_print << "/" << sum << endl;
}

void process_sdp_solution(ProblemConfig &problem_config) {
  ifstream results;
  string result_name = problem_config.get_sdp_solution_name();
  results.open(result_name.c_str(), ifstream::in);
  if (!results.good()) {
    cerr << "Failed opening file with CSDP result " << result_name << endl;
    exit(1);
  }
  cerr << "Processing CSDP solution file: " << result_name << endl;
  cout << endl;
  process_csdp_solution(results, cout, problem_config.largest_flag_to_show,
                        problem_config.skip_empty_flags_when_processing);
  results.close();
  results.open(result_name.c_str(), ifstream::in);
  if (!results.good()) {
    cerr << "Failed opening file with CSDP result " << result_name << endl;
    exit(1);
  }
  ofstream outfile;
  outfile.open(problem_config.get_sdp_extremal_construction_name().c_str(),
               ofstream::out);
  if (!outfile.good()) {
    cerr << "Failed opening file "
         << problem_config.get_sdp_extremal_construction_name() << endl;
    exit(1);
  }
  cerr << "Saving full densities of extremal construction to file: "
       << result_name << endl;
  process_csdp_solution(results, outfile, V, false);
  outfile.close();
}
