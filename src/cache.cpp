#include "cache.hpp"
#include "flag-generator.hpp"
#include "utils.hpp"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unordered_map>

// Cache of operations on flags.
unordered_map<string, vector<MultiplicationPartialResult>> multiplication_cache;
unordered_map<string, vector<ProjectionPartialResult>> projection_cache;
unordered_map<string, vector<FlagProductPartialResult>> flag_product_cache;

std::istream &operator>>(std::istream &in,
                         MultiplicationPartialResult &partial_result) {

  in >> partial_result.idx_multiplied >> partial_result.idx_left >>
      partial_result.idx_right >> partial_result.good_maps >>
      partial_result.all_maps;

  return in;
}

std::istream &operator>>(std::istream &in,
                         ProjectionPartialResult &partial_result) {

  in >> partial_result.idx_projected >> partial_result.idx_original >>
      partial_result.good_maps >> partial_result.all_maps;
  return in;
}

std::ostream &operator<<(std::ostream &stream,
                         const FlagProduct &flag_product) {
  stream << flag_product.type_idx << " " << flag_product.f1_idx << " "
         << flag_product.f2_idx;
  return stream;
}

std::istream &operator>>(std::istream &in, FlagProduct &flag_product) {
  in >> flag_product.type_idx >> flag_product.f1_idx >> flag_product.f2_idx;
  return in;
}

std::istream &operator>>(std::istream &in,
                         FlagProductPartialResult &partial_result) {
  in >> partial_result.product >> partial_result.count;
  return in;
}

bool operator==(const FlagProduct &a, const FlagProduct &b) {
  return a.type_idx == b.type_idx && a.f1_idx == b.f1_idx &&
         a.f2_idx == b.f2_idx;
}

template <typename CacheType, typename CacheFileWriter>
void write_cache_to_file(string filename, CacheFileWriter writing_function) {
  ofstream outfile;
  outfile.open(filename.c_str(), ofstream::out);
  if (!outfile.good()) {
    cerr << "Failed opening file " << filename << endl;
    exit(1);
  }
  writing_function(outfile);
}

template <typename CacheType, typename CacheFileWriter>
const vector<CacheType> load_cache_from_file(string filename,
                                             string name_prefix,
                                             CacheFileWriter writing_function) {
  create_dir(name_prefix);
  ifstream infile;
  infile.open(filename.c_str(), ifstream::in);
  if (!infile.good()) {
    cerr << "Couldn't open cache file " << filename << ", generating..."
         << endl;
    write_cache_to_file<CacheType>(filename, writing_function);
  }

  infile.close();
  infile.open(filename.c_str(), ifstream::in);
  if (!infile.good()) {
    cerr << "Couldn't open filename " << filename
         << " even after generating cache, exiting..." << endl;
    exit(1);
  }

  string line;
  vector<CacheType> result;
  while (getline(infile, line)) {
    istringstream iss(line);
    CacheType partial_result;
    iss >> partial_result;
    result.push_back(partial_result);
  }
  return result;
}

template <typename CacheType, typename CacheFileWriter>
const vector<CacheType> &
get_cache(string filename, string name_prefix, CacheFileWriter writing_function,
          unordered_map<string, vector<CacheType>> &memory_cache) {
  if (memory_cache.find(filename) != memory_cache.end()) {
    return memory_cache[filename];
  }
  memory_cache[filename] =
      load_cache_from_file<CacheType>(filename, name_prefix, writing_function);
  return memory_cache[filename];
}

const vector<MultiplicationPartialResult> &
get_multiplication_cache(flag type, int size_left, int size_right) {
  int n = size_left + size_right - type.m_vertices;

  stringstream filename_stream;
  string type_name = type.print("");
  string name_prefix = filename_prefix("MULT");
  create_dir(name_prefix);
  filename_stream << name_prefix << "/n" << n << "__left" << size_left
                  << "__right" << size_right << "__" << type_name << ".txt";
  string filename = filename_stream.str();

  auto multiplication_cache_writing_function = [&](ofstream &outfile) {
    hydrate_multiplication_cache(type, size_left, size_right, outfile);
  };
  return get_cache(filename, name_prefix, multiplication_cache_writing_function,
                   multiplication_cache);
}

const vector<ProjectionPartialResult> &
get_projection_cache(flag type, int vertices, int new_vertices,
                     int remaining_labeled) {
  stringstream filename_stream;
  string type_name = type.print("");
  string name_prefix = filename_prefix("PROJ");
  create_dir(name_prefix);
  filename_stream << name_prefix << "/n" << vertices << "__new_n"
                  << new_vertices << "__new_theta" << remaining_labeled << "__"
                  << type_name << ".txt";
  string filename = filename_stream.str();

  auto projection_cache_writing_function = [&](ofstream &outfile) {
    hydrate_projection_cache(type, vertices, new_vertices, remaining_labeled,
                             outfile);
  };
  return get_cache(filename, name_prefix, projection_cache_writing_function,
                   projection_cache);
}

const vector<FlagProductPartialResult> &
get_flag_product_cache(int Kn, int idx_in_list, int typesize) {
  stringstream filename_stream;
  string name_prefix = filename_prefix("PROD");
  create_dir(name_prefix);
  filename_stream << name_prefix << "/n" << Kn << "__i" << idx_in_list
                  << "__typesize" << typesize << ".txt";
  string filename = filename_stream.str();
  auto flag_product_cache_writing_function = [&](ofstream &outfile) {
    hydrate_flag_product_cache(Kn, idx_in_list, typesize, outfile);
  };
  return get_cache(filename, name_prefix, flag_product_cache_writing_function,
                   flag_product_cache);
}
