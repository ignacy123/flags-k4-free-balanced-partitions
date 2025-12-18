#include "option-parser.hpp"
#include "config.hpp"
#include <iostream>

using namespace std;
void extra_arguments(int extra, int processed, int argc) {
  if (processed + extra >= argc) {
    cerr << "Expected more arguments..." << endl;
    exit(1);
  }
}

void help() {
  cout << "This program is used to show densities of small graphs in the "
          "graphon attaining the proven bound."
       << endl;
  cout << "They are obtained from the solution to the dual problem, which by "
          "itself are densities of the largest flags used."
       << endl;
  cout << "However, these densities are difficult to read by themselves so it "
          "makes more sense to aggregate them to smaller flags, which is what "
          "this program does."
       << endl;
  cout << "Following flags are supported:" << endl;
  cout << "-l -- maximum size of the flag that should be actually printed"
       << endl;
  cout << "-s -- directory with the SDP problem and solution." << endl;
  cout << "-d -- directory with the flags." << endl;
  cout << "-c -- number of the case." << endl;
  cout
      << "-a -- using this switch blocks skipping flags equal to 0 (useful if "
         "you're planning to process the solution programmatically afterwards.."
      << endl;
}

void parse_options(int argc, char *argv[], ProblemConfig &problem_config) {
  int i = 1;
  while (i < argc) {
    if (string(argv[i]) == "-help") {
      help();
    }
    if (string(argv[i]) == "-l") {
      extra_arguments(1, i, argc);
      problem_config.largest_flag_to_show = atoi(argv[++i]);
    }
    if (string(argv[i]) == "-c") {
      extra_arguments(1, i, argc);
      problem_config.case_number = atoi(argv[++i]);
    }
    if (string(argv[i]) == "-s") {
      extra_arguments(1, i, argc);
      problem_config.sdp_directory = argv[++i];
    }
    if (string(argv[i]) == "-d") {
      extra_arguments(1, i, argc);
      problem_config.data_directory = argv[++i];
    }
    if (string(argv[i]) == "-a") {
      problem_config.skip_empty_flags_when_processing = false;
    }
    i++;
  }
}
