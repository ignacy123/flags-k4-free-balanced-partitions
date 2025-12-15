#include "utils.hpp"
#include <cassert>
#include <iostream>
#include <math.h>
#include <sstream>
#include <sys/stat.h>

string time_to_str(time_t time_taken) {

  int days_t = time_taken / 60 / 60 / 24;
  int hours_t = (time_taken % (60 * 60 * 24)) / 60 / 60;
  int min_t = (time_taken % (60 * 60)) / 60;
  int sec_t = time_taken % (60);

  stringstream ss;
  ss << days_t << "d+" << hours_t << "h+" << min_t << "m+" << sec_t << "s";

  return ss.str();
}

// Try if the number should be actuallly rounded
double smart_round(double d, double precision) {
  double dr = round(d);
  if (abs(dr - d) < precision)
    return dr;
  return d;
}

double g_smart_round_precision = 0.00000001;
double smart_round(double d) { return smart_round(d, g_smart_round_precision); }

int binomial(int n, int k) {
  int b = 1;
  for (int i = 0; i < k; ++i) {
    b *= (n - i);
    b /= (i + 1);
  }
  return b;
}

bool create_dir(string path) {
  if (mkdir(path.c_str(), 0755) != 0) {
    if (errno == EEXIST)
      return true;
    return false;
  }
  return true;
}

bool starts_with_flag_expression(string filename) {
  stringstream ss(filename);

  char ch;
  ss >> ch;

  return isdigit(ch) || ch == '-' || ch == '(' || ch == '.' || ch == '{' ||
         ch == '[';

  // Below is an old version, we need ( for flag calculator
  double number;
  ss >> number;
  return !ss.fail();
}
