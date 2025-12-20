#pragma once
#include "cache.hpp"
#include "flag-generator.hpp"
#include "flag.hpp"
#include "utils.hpp"
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>
#include <vector>

template <int theta, int vertices> class FlagVector {
public:
  FlagVector<theta, vertices>(const flag &blueprint = flag(),
                              double coefficient = 0.)
      : type([&]() {
          flag type;
          blueprint.get_type_subflag(type);
          return type;
        }()) {
    assert(type.m_Theta == theta);
    get_labeled_flags_of_one_type(vertices, type, flags);
    coefficients.assign(flags.size(), coefficient);
    if (blueprint.m_Theta == blueprint.m_vertices &&
        blueprint.m_vertices != vertices) {
      return;
    }
    if (blueprint.m_vertices != vertices) {
      cerr << "Instantiated vector with a flag that isn't a plain type or a "
              "compatible flag"
           << endl;
      exit(1);
    }
    (*this) += blueprint;
  };

  static FlagVector<theta, vertices> from_vector(flag type,
                                                 vector<flag> flags) {
    FlagVector<theta, vertices> result(type);
    for (int i = 0; i < flags.size(); i++) {
      assert(flags[i].m_vertices == vertices);
      assert(flags[i].have_same_type(type));
      result += flags[i];
    }
    return result;
  }

  static FlagVector<theta, vertices> load_from_file(string filename,
                                                    flag type) {
    FlagVector<theta, vertices> result(type);
    ifstream infile;
    infile.open(filename.c_str(), ifstream::in);
    if (!infile.good()) {
      cerr << "Failed opening file " << filename;
      exit(1);
    }

    FilteringIstream filteredinfile(infile);
    double coefficient;
    flag flag;
    while (true) {
      filteredinfile >> coefficient;
      if (filteredinfile.fail())
        break;
      flag.load_from_stream(filteredinfile, vertices, theta);
      flag_coeff fc;
      fc.g = flag;
      fc.coefficient = coefficient;
      cerr << fc << endl;
      FlagVector<theta, vertices> added(flag);
      added = added * coefficient;
      result += added;
    }

    return result;
  }

  void ensure_coefficients_are_integers() {
    double error = 0.000000000001;
    for (double coefficient : coefficients) {
      int rounded = round(coefficient);
      if (abs(rounded - coefficient) > error) {
        cerr << "Non-integer coefficient found: " << coefficient
             << ", exiting..." << endl;
        exit(1);
      }
    }
  }

  const vector<flag> &get_flags() { return flags; }
  const vector<double> &get_coefficients() const { return coefficients; }
  const flag &get_type() { return type; }

  bool operator==(const FlagVector<theta, vertices> &other) {
    if (!type.is_identical_to(other.type)) {
      return false;
    }
    if (coefficients != other.coefficients) {
      return false;
    }
    return true;
  }

  void debug_print(std::ostream &stream) {
    stream << flags.size() << endl;
    print(stream, false);
  }

  friend ostream &operator<<(std::ostream &stream,
                             const FlagVector<theta, vertices> &flag_vector) {
    flag_vector.print(stream);
    return stream;
  };

  void operator+=(double constant) {
    for (int i = 0; i < flags.size(); i++) {
      coefficients[i] += constant;
    }
  }

  void operator-=(double constant) {
    for (int i = 0; i < flags.size(); i++) {
      coefficients[i] -= constant;
    }
  }

  void operator*=(double constant) {
    for (int i = 0; i < flags.size(); i++) {
      coefficients[i] *= constant;
    }
  }

  void operator+=(const flag &flag) { add_flag_with_coefficient(flag); };

  void operator+=(const flag_coeff &fc) {
    add_flag_with_coefficient(fc.g, fc.coefficient);
  };

  friend FlagVector<theta, vertices> &
  operator+=(FlagVector<theta, vertices> &flag_vector,
             const FlagVector<theta, vertices> &other) {
    assert(flag_vector.type.have_same_type(other.type));
    for (int i = 0; i < flag_vector.flags.size(); i++) {
      flag_vector.coefficients[i] += other.coefficients[i];
    }
    return flag_vector;
  };

  friend FlagVector<theta, vertices> &
  operator-=(FlagVector<theta, vertices> &flag_vector,
             const FlagVector<theta, vertices> &other) {
    flag_vector += (other * (-1));
    return flag_vector;
  };

  template <int new_theta, int vertices_left, int vertices_right>
  friend FlagVector<new_theta, vertices_left + vertices_right - new_theta>
  operator*(const FlagVector<new_theta, vertices_left> left,
            const FlagVector<new_theta, vertices_right> right);

  template <int new_theta, int new_vertices>
  friend FlagVector<new_theta, new_vertices>
  operator*(const FlagVector<new_theta, new_vertices> flag_vector,
            double coefficient);

  friend FlagVector<theta, vertices>
  operator+(double constant, const FlagVector<theta, vertices> &flag_vector) {
    FlagVector<theta, vertices> result(flag_vector.type);
    for (int i = 0; i < flag_vector.coefficients.size(); i++) {
      result.coefficients[i] = constant + flag_vector.coefficients[i];
    }
    return result;
  }

  friend FlagVector<theta, vertices>
  operator+(const FlagVector<theta, vertices> &flag_vector, double constant) {
    return constant + flag_vector;
  }

  friend FlagVector<theta, vertices>
  operator-(double constant, const FlagVector<theta, vertices> &flag_vector) {
    return constant + (-flag_vector);
  };

  friend FlagVector<theta, vertices>
  operator-(const FlagVector<theta, vertices> &flag_vector, double constant) {
    return flag_vector + (-constant);
  };

  FlagVector<theta, vertices>
  operator+(const FlagVector<theta, vertices> &other) {
    FlagVector<theta, vertices> result(type);
    for (int i = 0; i < flags.size(); i++) {
      result.coefficients[i] = coefficients[i] + other.coefficients[i];
    }
    return result;
  };

  FlagVector<theta, vertices> operator-() const {
    FlagVector<theta, vertices> result(type);
    for (int i = 0; i < flags.size(); i++) {
      result.coefficients[i] = -coefficients[i];
    }
    return result;
  }

  FlagVector<theta, vertices>
  operator-(const FlagVector<theta, vertices> &other) {
    return (*this) + (other * double(-1.));
  };

  template <int remaining_labeled, int new_vertices>
  FlagVector<remaining_labeled, new_vertices> project() {
    assert(remaining_labeled <= theta);
    assert(new_vertices >= vertices);
    flag type_copy = type;
    type_copy.set_Theta(remaining_labeled);
    flag reduced_type;
    type_copy.get_type_subflag(reduced_type);
    FlagVector<remaining_labeled, new_vertices> result(reduced_type);
    const vector<ProjectionPartialResult> &projection_cache =
        get_projection_cache(type, vertices, new_vertices, remaining_labeled);
    for (auto partial_result : projection_cache) {
      double coeff =
          (double)partial_result.good_maps / (double)partial_result.all_maps;
      result.coefficients[partial_result.idx_projected] +=
          coeff * coefficients[partial_result.idx_original];
    }
    return result;
  }

  double density(int m_vertices, int idx, int remaining_labeled = 0) {
    double density = 0;
    vector<ProjectionPartialResult> projection_cache =
        get_projection_cache(type, m_vertices, vertices, remaining_labeled);
    for (auto partial_result : projection_cache) {
      if (partial_result.idx_original != idx)
        continue;
      density +=
          ((double)partial_result.good_maps / (double)partial_result.all_maps) *
          coefficients[partial_result.idx_projected];
    }
    return density;
  }

private:
  void add_flag_with_coefficient(const flag &flag, double coefficient = 1) {
    assert(type.have_same_type(flag));
    bool found_at_least_one = false;
    for (int i = 0; i < flags.size(); i++) {
      if (flags[i].is_isomorphic_to(flag)) {
        coefficients[i] += coefficient;
        found_at_least_one = true;
      }
    }
    if (!found_at_least_one) {
      cerr << "Error: No isomorphic flag found" << endl;
      exit(1);
    }
  };

  void print(std::ostream &stream, bool skip_if_empty = true) const {
    for (int i = 0; i < flags.size(); i++) {
      flag_coeff fc;
      fc.coefficient = coefficients[i];
      fc.g = flags[i];
      fc.print(stream, skip_if_empty);
    }
  }

  vector<flag> flags;
  vector<double> coefficients;
  flag type;
  template <int, int> friend class FlagVector;
};

template <int new_theta, int vertices_left, int vertices_right>
FlagVector<new_theta, vertices_left + vertices_right - new_theta>
operator*(const FlagVector<new_theta, vertices_left> left,
          const FlagVector<new_theta, vertices_right> right) {
  if (!left.type.is_identical_to(right.type)) {
    cerr << "Error: tried multiplying vectors of different type" << endl;
    exit(1);
  }
  string type = left.type.print("");
  int new_vertices = vertices_left + vertices_right - new_theta;
  FlagVector<new_theta, vertices_left + vertices_right - new_theta> result(
      left.type);
  const vector<MultiplicationPartialResult> &multiplication_cache =
      get_multiplication_cache(left.type, vertices_left, vertices_right);
  for (auto partial_result : multiplication_cache) {
    double coeff =
        (double)partial_result.good_maps / (double)partial_result.all_maps;
    result.coefficients[partial_result.idx_multiplied] +=
        coeff * left.coefficients[partial_result.idx_left] *
        right.coefficients[partial_result.idx_right];
  }

  return result;
};

template <int new_theta, int new_vertices>
FlagVector<new_theta, new_vertices>
operator*(const FlagVector<new_theta, new_vertices> flag_vector,
          double coefficient) {
  FlagVector<new_theta, new_vertices> result(flag_vector.type);
  for (int i = 0; i < flag_vector.flags.size(); i++) {
    result.coefficients[i] = flag_vector.coefficients[i] * coefficient;
  }
  return result;
};

pair<bool, flag> get_edge_between(flag flag_left, flag flag_right);

template <int root_size>
FlagVector<root_size, root_size + 2> get_edges_inside(flag type,
                                                      vector<flag> flags) {
  FlagVector<root_size, root_size + 2> result(type);
  for (int i = 0; i < flags.size(); i++) {
    for (int j = i; j < flags.size(); j++) {
      auto edge_between = get_edge_between(flags[i], flags[j]);
      if (edge_between.first) {
        result += edge_between.second;
      }
    }
  }
  return result;
}

template <int root_size>
FlagVector<root_size, root_size + 2>
get_edges_between(flag type, vector<flag> left, vector<flag> right) {
  FlagVector<root_size, root_size + 2> result(type);
  for (const flag flag_left : left) {
    for (const flag flag_right : right) {
      auto edge_between = get_edge_between(flag_left, flag_right);
      if (edge_between.first) {
        result += edge_between.second;
      }
    }
  }
  return result;
}
