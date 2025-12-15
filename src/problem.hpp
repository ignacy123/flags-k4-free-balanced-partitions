#pragma once

#include "flag-calculator.hpp"

class Problem {
public:
  Problem(int scaling_factor = 1, bool ensure_integers = false)
      : scaling_factor(scaling_factor), ensure_integers(ensure_integers) {};
  template <int theta, int vertices>
  void add_constraint(FlagVector<theta, vertices> constraint,
                      bool multiply = true) {
    if (multiply) {
      FlagVector<theta, theta + V - vertices> filling_flags(
          constraint.get_type());
      for (flag filling_flag : filling_flags.get_flags()) {
        FlagVector<theta, theta + V - vertices> filling_vector(filling_flag);
        FlagVector<theta, V> multiplied = constraint * filling_vector;
        add_constraint_raw(multiplied);
      }
    } else {
      add_constraint_raw(constraint);
    }
  }

  template <int theta, int vertices>
  void add_objective(FlagVector<theta, vertices> raw_objective) {
    FlagVector<0, vertices> unlabeled =
        raw_objective.template project<0, vertices>();
    FlagVector<0, V> projected_up = unlabeled.template project<0, V>();
    objective = projected_up;
  }

  const vector<FlagVector<0, V>> &get_constraints() {
    cerr << "There are " << constraints.size() << " constraints." << endl;
    return constraints;
  }
  const FlagVector<0, V> &get_objective() { return objective; }

private:
  template <int theta, int vertices>
  void add_constraint_raw(FlagVector<theta, vertices> constraint) {
    auto projected = constraint.template project<0, V>();
    projected *= scaling_factor;
    if (ensure_integers) {
      projected.ensure_coefficients_are_integers();
    }
    constraints.push_back(projected);
  }

  vector<FlagVector<0, V>> constraints;
  FlagVector<0, V> objective;
  int scaling_factor = 1;
  bool ensure_integers = false;
};

template <int theta, int vertices>
void show_vector(FlagVector<theta, vertices> flag_vector, string name) {
  auto unlabeled = flag_vector.template project<0, V>();
  cout << name << ": " << endl;
  unlabeled.debug_print(cout);
  cout << endl;
}
