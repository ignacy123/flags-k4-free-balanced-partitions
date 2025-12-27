#include "config.hpp"
#include "flag-calculator.hpp"
#include "option-parser.hpp"
#include <cassert>
#include <string>

double BOUND = 0.0662;

int main(int argc, char *argv[]) {
  ProblemConfig::instance().data_directory = "graph-checks";
  ProblemConfig::instance().data_directory = "graph-checks-forbidden";
  ProblemConfig::instance().upper_bound = false;
  parse_options(argc, argv, ProblemConfig::instance());
  FlagVector<0, V> flags;
  flags *= 2;
}
