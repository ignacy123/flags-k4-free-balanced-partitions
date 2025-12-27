#include "flag-generator.hpp"
#include "cache.hpp"
#include "config.hpp"
#include "flag.hpp"
#include "utils.hpp"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <sstream>
#include <string>
#include <unordered_map>

vector<flag> g_forbidden_subflags;

// This should make a cash of loaded flags and have the things a little faster
unordered_map<string, vector<flag>> g_labeled_flags_of_one_type_map;
vector<flag> g_unlabeled_flags[V + 1];
vector<flag>
    g_types[V + 1]; // types - it is just same as g_unlabeled_flags, but all
                    // vertices are labeled :-) Used for reduced types
vector<vector<flag>> g_flags; // flags to process using multiplication - every
                              // type&size is a separate list

vector<vector<flag>> &get_flags() { return g_flags; }
vector<flag> &get_unlabeled_flags(int Kn) { return g_unlabeled_flags[Kn]; }
vector<flag> &get_forbidden_subflags() {
  // This will always run if we don't use any forbidden subflags,
  // but I don't really care.
  if (g_forbidden_subflags.size() == 0)
    load_forbidden();
  return g_forbidden_subflags;
}

#ifdef G_COLORED_VERTICES
vector<int> g_vertex_color_pattern;
#endif

string filename_prefix(string type_of_file, string directory) {
  stringstream filename;

  if (directory == "") {
    directory = ProblemConfig::instance().data_directory;
  }

  create_dir(directory);
  filename << directory << "/" << type_of_file
#ifdef G_COLORED_VERTICES
           << "_vertices" << COLORS_VERTICES - 1
#else
#endif
#ifdef G_COLORED_EDGES
           << "_edges" << COLORS_EDGES - 1
#endif
#ifdef G_ORIENTED_EDGES
           << "oriented"
#endif
      ;
  return filename.str();
}

#ifdef G_COLORED_VERTICES
void add_color(int color) { g_vertex_color_pattern.push_back(color); }
#endif

bool load_flags_from_file(string filename, vector<flag> &flag_list,
                          int verbose_output = 0) {
  OPEN_FILE_SMARTLY_RETURN_FALSE_ON_FAIL(istr, filename);

  if (verbose_output)
    cerr << "Loading labeled flags from file " << filename << endl;

  flag f;
  while (f.load_from_stream((*istr), -1, -1)) {
    flag_list.push_back(f);
  }

  /// infile.close();

  return true;
}

bool dump_flags_to_file(string filename, vector<flag> &flag_list) {
  ofstream outfile;
  outfile.open(filename.c_str(), ofstream::out);
  if (!outfile.good()) {
    cerr << "Failed opening file " << filename << endl;
    return false;
  }

  cerr << "Writing flags to file " << filename << endl;

  for (unsigned int x = 0; x < flag_list.size(); x++) {
    outfile << flag_list[x].print() << endl;
  }

  outfile.close();

  return true;
}

bool load_flags_and_coefficients_from_file(string filename,
                                           vector<flag_coeff> &flag_list,
                                           int verbose_output) {

  OPEN_FILE_SMARTLY_RETURN_FALSE_ON_FAIL(istr, filename);

  if (verbose_output)
    cerr << "Loading labeled flags from file " << filename << endl;

  flag_coeff fc;
  (*istr) >> fc.coefficient;
  cerr << "Loaded coefficient " << fc.coefficient << endl;
  while ((*istr) && fc.g.load_from_stream((*istr), -1, -1)) {
    flag_list.push_back(fc);
    (*istr) >> fc.coefficient;
  }

  // infile.close();

  return true;
}

bool load_labeled_flags_from_file(int sizeKn, int verbose_output) {
  stringstream filename;
  filename << filename_prefix() << "__n" << sizeKn << "_labeled.txt";

  ifstream infileF;
  infileF.open(filename.str().c_str(), ifstream::in);
  if (!infileF.good()) {
    cerr << "Failed opening file with labeled flags " << filename.str() << endl;
    return false;
  }

  FilteringIstream infile(infileF);

  /*
      ifstream infile;
      infile.open (filename.str().c_str(), ifstream::in);
      if (!infile.good())
      {
          cerr << "Failed opening file with labeled flags " << filename.str() <<
     endl; return false;
      }
  */

  if (verbose_output) {
    cerr << "Loading labeled flags from file " << filename.str() << endl;
  }

  flag f;
  while (f.load_from_stream(infile, -1, -1)) {
    include_flag_in_list(f, get_flags());
  }

  infileF.close();

  return true;
}

void dump_flag_and_coefficient(const flag_coeff &fc,
                               bool use_smart_round = false) {
  if (fc.coefficient == 0)
    return;

  cout.precision(G_PRECISION);
  if (use_smart_round) {
    if (smart_round(fc.coefficient) == 0)
      return;
    cout << smart_round(fc.coefficient);
  } else
    cout << fc.coefficient;
  cout << "  " << fc.g.print() << endl;
}

void dump_flags_and_coefficients(vector<flag_coeff> &flag_list,
                                 bool use_smart_round = false) {
  for (unsigned int x = 0; x < flag_list.size(); x++) {
    dump_flag_and_coefficient(flag_list[x], use_smart_round);
  }
}

// final_flags means that the flags will not be used any further
// for genertin other flags. Important when combined with
// G_USE_FIRST_EDGE_COLOR_FOR_BLOWUP_ONLY
bool load_unlabeled_flags_from_file(
    int sizeKn, bool final_flags = true,
    bool remove_duplicates_while_loading = false,
    bool remove_forbidden_wile_loading = false) {
  stringstream filename;

  filename << filename_prefix() << "__n" << sizeKn << "_unlabeled.txt";

  ifstream infile;
  infile.open(filename.str().c_str(), ifstream::in);
  if (!infile.good()) {
    cerr << "Failed opening file with unlabeled flags " << filename.str()
         << endl;
    return false;
  }

  /*
  for (int i = 0; i < sizeKn; i++)
  {
      if (g_unlabeled_flags[i].size() > 0)
      {
          cerr << "Flags of given size already exists, cannot load again." <<
  endl; assert(0);
      }
  }
  */

  cerr << "Loading unlabeled flags from file " << filename.str() << endl;

  // For large instances, it really helps to reserve the space in containers
  //  before loading... so we read the file twice
  int sizes_cnt[V + 1], duplicates_cnt[V + 1], forbidden_cnt[V + 1];
  for (int i = 0; i < V + 1; i++) {
    sizes_cnt[i] = 0;
    duplicates_cnt[i] = 0;
    forbidden_cnt[i] = 0;
  }
  flag g;
  while (g.load_from_stream(infile, -1, 0)) {
    sizes_cnt[g.m_vertices]++;
  }

  for (int i = 0; i < V + 1; i++) {
    // cerr << "Expecting for size " << i << " flag number: " << sizes_cnt[i] <<
    // endl;
    g_unlabeled_flags[i].clear();
    g_unlabeled_flags[i].reserve(sizes_cnt[i]);
  }

  //    infile.seek(0);
  infile.close();
  infile.open(filename.str().c_str(), ifstream::in);

  int tested = 0;

  while (g.load_from_stream(infile, -1, 0)) {

    tested++;
    assert(sizeKn >= g.m_vertices);

    if (remove_forbidden_wile_loading && is_flag_forbidden(g)) {
      forbidden_cnt[g.m_vertices]++;
      continue;
    }

    if (remove_duplicates_while_loading) {
      volatile bool new_flag = true;

#pragma omp parallel for shared(new_flag)
      for (int i = 0; i < (int)g_unlabeled_flags[g.m_vertices].size(); i++) {

        if (new_flag == false)
          continue;

        if (g.is_isomorphic_to(g_unlabeled_flags[g.m_vertices][i])) {
#pragma omp atomic write
          new_flag = false;
          // This allows to break when a duplicate found - may be usefull to
          // check if the duplicates are indeed duplicates or if it is a mistake
          // in the program....
          // if (duplicates_cnt[g.m_vertices] == 2)
          //{
          //    cerr << "Found duplicates" << endl;
          //    cout << g.print() << endl;
          //    cout << g_unlabeled_flags[g.m_vertices][i].print() << endl;
          //    exit(0);
          //}
#ifndef _USING_OMP_
          break;
#endif
        }
      }
      if (!new_flag) {
        duplicates_cnt[g.m_vertices]++;
        continue;
      }
    }
    g_unlabeled_flags[g.m_vertices].push_back(g);
    // g.reverse_rotation_system();
    // g_unlabeled_flags[g.m_vertices].push_back(g);
    // cerr << g.print() << " " << g_unlabeled_flags[g.m_vertices].size() << "
    // tested " << tested << endl;
  }

  infile.close();

  for (int i = 0; i <= sizeKn; i++) {
    cerr << "Loaded # of unlabeled flags of size " << i << " is "
         << g_unlabeled_flags[i].size();
    if (remove_forbidden_wile_loading) {
      cerr << " and " << forbidden_cnt[i] << " forbiddens";
    }
    if (remove_duplicates_while_loading) {
      cerr << " and " << duplicates_cnt[i] << " duplicates";
    }
    cerr << endl;
  }

  return true;
}

bool dump_unlabeled_flags(int sizeKn) {
  stringstream filename;
  filename << filename_prefix() << "__n" << sizeKn << "_unlabeled.txt";

  ofstream outfile;
  outfile.open(filename.str().c_str(), ofstream::out);
  if (!outfile.good()) {
    cerr << "Failed opening file " << filename.str() << endl;
    return false;
  }

  cerr << "Writing unlabeled flags to file " << filename.str() << endl;

  for (int f = 0; f <= V + 1; f++) {
    for (unsigned int x = 0; x < g_unlabeled_flags[f].size(); x++) {
      outfile << g_unlabeled_flags[f][x].print() << endl;
    }
  }

  outfile.close();

  return true;
}

bool dump_labeled_flags(int sizeKn) {
  stringstream filename;
  filename << filename_prefix() << "__n" << sizeKn << "_labeled.txt";

  ofstream outfile;
  outfile.open(filename.str().c_str(), ofstream::out);
  if (!outfile.good()) {
    cerr << "Failed opening file " << filename.str() << endl;
    return false;
  }

  cerr << "Writing labeled flags to file " << filename.str() << endl;

  for (int f = 0; f < (int)get_flags().size(); f++) {
    for (unsigned int x = 0; x < get_flags()[f].size(); x++) {
      outfile << get_flags()[f][x].print() << endl;
    }
    if (f < (int)get_flags().size() - 1)
      outfile << endl;
  }

  outfile.close();

  return true;
}

void generate_all_unlabeled_subflags_from_size(int sizeKn, int verbose_output) {
  cerr << "Generating smaller unlabeled flags from size " << sizeKn << endl;

  int mapping[V + 1];

  for (int n = sizeKn; n > 0; n--) {
    // Delete the ones we are just about to start generating
    g_unlabeled_flags[n - 1].clear();

    // graph to be processed
    for (int id = 0; id < (int)g_unlabeled_flags[n].size(); id++) {
      // vertex to be skipped
      for (int skip = 0; skip < n; skip++) {
        int mapID = 0;
        for (int u = 0; u < n; u++) {
          if (u == skip)
            continue;
          mapping[mapID++] = u;
        }

        flag F;
        F.as_subflag(g_unlabeled_flags[n][id], mapping, n - 1, 0);

        volatile bool new_flag = true;

#pragma omp parallel for shared(new_flag)
        for (int i = 0; i < (int)g_unlabeled_flags[F.m_vertices].size(); i++) {
          if (new_flag == false)
            continue;

          if (F.is_isomorphic_to(g_unlabeled_flags[F.m_vertices][i])) {
#pragma omp atomic write
            new_flag = false;
#ifndef _USING_OMP_
            break;
#endif
          }
        }
        if (!new_flag)
          continue;

        g_unlabeled_flags[F.m_vertices].push_back(F);
      }
    }
    if (verbose_output) {
      cerr << "# of unlabeled flags of size " << n - 1 << " is "
           << g_unlabeled_flags[n - 1].size() << endl;
    }
  }
}

bool g_already_in_known_flags(flag &g, vector<flag> &flag_list) {
  for (unsigned int i = 0; i < flag_list.size(); i++) {
    if (flag_list[i].is_isomorphic_to(g)) {
      return true;
    }
  }
  return false;
}

void add_g_to_known_flags(flag &g, vector<flag> &flag_list) {
  flag_list.push_back(g);
  // if (g.m_3edges_cnt >= 20)
  {
    //    cerr << g.print() << endl;
  }
}

#ifdef G_COLORED_EDGES
void try_color_edge(flag &g, int u, int v, vector<flag> &flag_list) {
  if (v >= g.m_vertices || v == u) {
    u++;
    v = u + 1;
  }
  if (v >= g.m_vertices) {
    if (is_flag_forbidden(g))
      return;

    g.create_minlex_signature();

    if (!g_already_in_known_flags(g, flag_list))
      add_g_to_known_flags(g, flag_list);
    //		add_g_to_known_flags(flagh_list);
    return;
  }

  // Color only uncolored edges.....
  if (g.m_color_edge[u][v] == 0) {
    // Either color 1 us used for blow-up or for tournaments it is not used at
    // all since color 1 does not change orientation
    for (int color = 1; color < COLORS_EDGES; color++) {
      g.color_edge(u, v, color);
      try_color_edge(g, u, v + 1, flag_list);
#ifdef G_ORIENTED_EDGES
      if (color > G_ORIENTED_EDGES_UNORIENTED_COLORS) {
        g.color_edge(v, u, color);
        try_color_edge(g, u, v + 1, flag_list);
      }
#endif
    }
    g.color_edge(u, v, 0);
  } else {
    try_color_edge(g, u, v + 1, flag_list);
  }
}
#endif

void try_extensions_of_g_to_last_vertex_edges(flag &g,
                                              vector<flag> &flag_list) {

  int v = g.m_vertices;

#ifdef G_COLORED_EDGES
  try_color_edge(g, v - 1, 0, flag_list);
  return;
#endif

  v++; // this is to avoid warning of unused v
  assert(0);
}

void try_extensions_of_g_to_last_vertex(flag &g, vector<flag> &flag_list) {
#ifdef G_COLORED_VERTICES
  int v = g.m_vertices;

  if ((int)g_vertex_color_pattern.size() > v - 1) {
    // cerr << "USing it" << g_vertex_color_pattern[v-1] << endl;
    g.color_vertex(v - 1, g_vertex_color_pattern[v - 1]);
    try_extensions_of_g_to_last_vertex_edges(g, flag_list);
    return;
  }

  for (int c = 1; c < COLORS_VERTICES; c++) {
    // if (g.m_colored_vertices[c]+1 > 5)
    //{
    //     continue;
    // }
    g.color_vertex(v - 1, c);
    try_extensions_of_g_to_last_vertex_edges(g, flag_list);
  }
#else
  try_extensions_of_g_to_last_vertex_edges(g, flag_list);
#endif
}

void dump_vector_of_flags(vector<flag> *flag_lists, int previous_size,
                          int sizeKn) {
  stringstream filename;
  filename << filename_prefix() << "__n" << sizeKn << "_unlabeled_dump.txt";

  ofstream outfile;
  outfile.open(filename.str().c_str(), ofstream::out);
  if (!outfile.good()) {
    cerr << "Failed opening file " << filename.str() << endl;
    return;
  }

  // dumping of flag lists...
  for (int i = 0; i < previous_size; i++) {
    for (int f = 0; f < (int)flag_lists[i].size(); f++) {
      outfile << flag_lists[i][f].print() << endl;
    }
    outfile << endl;
  }

  outfile.close();
}
inline string get_filename_for_labeled_flags(int flag_size, const flag &type) {
  // Saving of the flags to a separate file and/or loading them
  return filename_prefix() + "__size" + to_string((long long)flag_size) +
         "_type" + type.print("") + ".txt";
}

bool labeled_flags_of_one_already_exist(int flag_size, const flag &type) {
  string filename = get_filename_for_labeled_flags(flag_size, type);

  // if it is already in the map, it exists
  if (g_labeled_flags_of_one_type_map.find(filename) !=
      g_labeled_flags_of_one_type_map.end())
    return true;

  vector<flag> flag_list;
  if (load_flags_from_file(filename, flag_list)) {
    g_labeled_flags_of_one_type_map.insert(make_pair(filename, flag_list));
    return true;
  }

  return false;
}

void generate_next_map(int type_size, int flag_size, int id,
                       vector<int> &current_map, vector<bool> &used_in_map,
                       vector<vector<int>> &mappings_to_try) {
  // mapping found
  if (id >= flag_size) {
    mappings_to_try.push_back(current_map);
    return;
  }

  // finish mapping of unlabeled vertices in any order
  if (id >= type_size) {
    for (int i = 0; i < flag_size; i++) {
      if (used_in_map[i] == false) {
        used_in_map[i] = true;
        current_map[id] = i;
        generate_next_map(type_size, flag_size, id + 1, current_map,
                          used_in_map, mappings_to_try);
        used_in_map[i] = false;
        return;
      }
    }
    assert(0);
  }

  // Try all possibilities for labeled vertices
  for (int i = 0; i < flag_size; i++) {
    if (used_in_map[i] == false) {
      used_in_map[i] = true;
      current_map[id] = i;
      generate_next_map(type_size, flag_size, id + 1, current_map, used_in_map,
                        mappings_to_try);
      used_in_map[i] = false;
    }
  }
}

void generate_maps(int type_size, int flag_size,
                   vector<vector<int>> &mappings_to_try) {
  vector<int> current_map;
  vector<bool> used_in_map;
  for (int i = 0; i < flag_size; i++) {
    current_map.push_back(i);
    used_in_map.push_back(false);
  }
  generate_next_map(type_size, flag_size, 0, current_map, used_in_map,
                    mappings_to_try);
}

// The following function megres dest_flags with src_flags and clears src_flags
void merge_flags_lists_parallel(vector<flag> &dest_flags,
                                vector<flag> &src_flags) {

  //    cout << flag_list.size() << endl;

  vector<flag> new_flags;
  new_flags.reserve(src_flags.size());

#pragma omp parallel for
  for (int k = 0; k < (int)src_flags.size(); k++) {
    if (!g_already_in_known_flags(src_flags[k], dest_flags)) {
#pragma omp critical(merging_lists)
      {
        new_flags.push_back(src_flags[k]);
      }
    }
  }
  // #ifdef DONT_USE_C11
  dest_flags.insert(dest_flags.end(), new_flags.begin(), new_flags.end());
  src_flags.clear();
  vector<flag>().swap(src_flags); // free the memory
                                  // #else
  //     dest_flags.insert(dest_flags.end(),
  //     std::make_move_iterator(new_flags.begin()),
  //     std::make_move_iterator(new_flags.end()) ); src_flags.clear();
  //     src_flags.shrink_to_fit(); // free the memory - does not work well on
  //     icc
  // #endif
}

void generate_unlabeled_flags_of_size(int i, int verbose_output,
                                      bool dump_unlabeled_while_generating) {

  int previous_size = (int)g_unlabeled_flags[i - 1].size();

  // vector <flag> flag_lists[previous_size]; // does not work on mac (llvm :-(
  // )
  vector<flag> *flag_lists;
  flag_lists = new vector<flag>[previous_size];
  if (flag_lists == NULL) {
    cerr << "Memory allocation failed" << endl;
    exit(1);
  }

  if ((int)g_unlabeled_flags[i - 1].size() < 1) {
    cerr << "No flags of size " << i - 1 << endl;
    cerr << "Generating flags of size " << i << " failed." << endl;
    exit(0);
  }

  //            #pragma omp parallel for
#pragma omp parallel for ordered schedule(dynamic)
  for (int j = 0; j < (int)g_unlabeled_flags[i - 1].size(); j++) {
    flag_lists[j].reserve(previous_size);

    flag g;
    g.set_vertices(i);
    g.copy_from(g_unlabeled_flags[i - 1][j]);
    // vector <flag> valid_extensions_of_g;
    // cout << "going for " << j << endl;
    try_extensions_of_g_to_last_vertex(g, flag_lists[j]);
    if (verbose_output) {
      cerr << "generated size " << i << ": " << j << "/"
           << (int)g_unlabeled_flags[i - 1].size() << " found "
           << flag_lists[j].size() << " flags " << endl;
    }
  }
  if (verbose_output) {
    cerr << "Flags generated. Merging part starts." << endl;
  }

  if (dump_unlabeled_while_generating)
    dump_vector_of_flags(flag_lists, previous_size, i);

  for (int power = 1; power < previous_size; power *= 2) {
    if (verbose_output) {
      cerr << "Trying " << power << " in " << previous_size << endl;
    }
    for (int j = 0; j < previous_size; j += 2 * power) {
      if (j + power < previous_size) {
        if (verbose_output) {
          cerr << "Merging  " << j << "<-" << j + power << "  of sizes "
               << flag_lists[j].size() << " " << flag_lists[j + power].size()
               << endl;
        }
        merge_flags_lists_parallel(flag_lists[j], flag_lists[j + power]);
        // flag_lists[j+power].clear();
      }
    }
    if (dump_unlabeled_while_generating)
      dump_vector_of_flags(flag_lists, previous_size, i);
  }

  g_unlabeled_flags[i].swap(flag_lists[0]);

  delete[] flag_lists;
}

void generate_unlabeled_flags_of_size(int Kn, bool force_generate_flags,
                                      bool remove_duplicates_while_loading,
                                      bool remove_forbidden_wile_loading,
                                      int verbose_output,
                                      bool dump_unlabeled_while_generating) {
  // Already loaded
  if (g_unlabeled_flags[Kn].size() > 0) {
    return;
  }

  // Getting unlabeled flags...
  if (force_generate_flags ||
      !load_unlabeled_flags_from_file(Kn, true, remove_duplicates_while_loading,
                                      remove_forbidden_wile_loading)) {
    // Try to reuse previous graphs...
    int loaded_size = 0;

    // blow-upping and generating would not work correctly...
    if (!force_generate_flags) {
      loaded_size = Kn - 1;
      cerr << "Trying to load file with smaller flags..." << endl;
      while (loaded_size > 0 &&
             !load_unlabeled_flags_from_file(loaded_size, false))
        loaded_size--;
    }

    if (loaded_size == 0) {
      cerr << "Unable to load any smaller graphs. Starting from empty graph."
           << endl;
      flag g_zero;
      g_zero.set_vertices_and_Theta(0, 0);
      assert(g_unlabeled_flags[0].size() == 0);
      g_unlabeled_flags[0].push_back(g_zero);
    }

    cerr << "Generating unlabeled flags of sizes " << loaded_size + 1 << " to "
         << Kn << endl;

    for (int i = loaded_size + 1; i <= Kn; i++) {
      generate_unlabeled_flags_of_size(i, verbose_output,
                                       dump_unlabeled_while_generating);
      cerr << "Unlabbeled flags of size " << i << ": "
           << g_unlabeled_flags[i].size() << endl;
    }

    dump_unlabeled_flags(Kn);
  }

  if ((int)g_unlabeled_flags[Kn].size() < 1) {
    cerr << "No flags of size " << Kn << endl;
    cerr << "Generating flags of size " << Kn << " failed." << endl;
    exit(0);
  }

  /*
  // And assert test
  for (int sz = 0; sz < Kn; sz++)
  {
      for (int i = 0; i < (int)g_unlabeled_flags[sz].size(); i++)
          for (int j = i+1; j < (int)g_unlabeled_flags[sz].size(); j++)
          {
              if
  (g_unlabeled_flags[sz][i].is_isomorphic_to(g_unlabeled_flags[sz][j]))
              {
                  cerr << endl;
                  cerr << "ERROR while generating flags of size " << Kn << " got
  for size " << sz << " identical flags" << endl; cerr << "i=" << i << ": " <<
  g_unlabeled_flags[sz][i].print() << " j=" << j << ": " <<
  g_unlabeled_flags[sz][j].print() << endl; assert(0);
              }

          }
  }
  */
}

// Not paralel version - very slow and not completely correct - in particular
// for ordered version.
void generate_labeled_flags_of_one_type(int flag_size, const flag &type,
                                        vector<flag> &flag_list) {

  cerr << "Generating labeled flags of size " << flag_size << " of type "
       << type.print() << endl;

  string filename = get_filename_for_labeled_flags(flag_size, type);

  // Make sure we have the unlabeled flags of right size loaded
  generate_unlabeled_flags_of_size(flag_size, false, false, false,
                                   ProblemConfig::instance().verbose_output);

  int to_label_size = (int)g_unlabeled_flags[flag_size].size();
  // cerr << "Type size " << type_size << " for " <<
  // g_unlabeled_flags[flag_size][i].print() << endl;
  int mapping[flag_size];
  flag F;

  int type_size = type.labeled_vertices_cnt();

  if (type_size != 0) {

    vector<vector<int>> mappings_to_try;
    generate_maps(type_size, flag_size, mappings_to_try);

    vector<flag> found_already;
    for (int i = 0; i < to_label_size; i++) {
      found_already.clear();
      // if (verbose_output)
      //      cerr << "Flag size "<<flag_size << " type size " << type_size << "
      //      labeling " << i << "/" << to_label_size << endl;

      // for (int j = 0; j < flag_size; j++) mapping[j]=j;

      // do {
      for (int mapID = 0; mapID < (int)mappings_to_try.size(); mapID++) {
        for (int j = 0; j < flag_size; j++)
          mapping[j] = mappings_to_try[mapID][j];

        F.as_subflag(g_unlabeled_flags[flag_size][i], mapping, flag_size,
                     type_size); // taking a subflag
        // cerr << "Found " << F.print() << " X " << F.have_same_type(type) << "
        // " << !find_flag_in_list(F, flag_list) << endl; if
        // (F.have_same_type(type) && find_flag_in_list(F, flag_list) == -1) if
        // (F.have_same_type(type) && find_flag_in_list_nonparalel(F, flag_list)
        // == -1)
        if (F.have_same_type(type) &&
            find_flag_in_list_nonparalel(F, found_already) == -1) {
          found_already.push_back(F);
          flag_list.push_back(F);
        }
        // cerr << "x";
      }
      //} while ( std::next_permutation(mapping,mapping+flag_size) );
      // cerr << endl;
      // cerr << "found_already.size()=" << found_already.size() << "
      // flag_list.size()=" << flag_list.size() << endl;
    }
  } else {
    // cerr << "Bubu " << flag_size << " "  <<
    // g_unlabeled_flags[flag_size].size() << endl; flag_list.erase();
    flag_list = g_unlabeled_flags[flag_size];
  }

  if (dump_flags_to_file(filename, flag_list)) {
    cerr << "Generated and written flags to " << filename << endl;
  }
}

// Not paralel version - very slow and not completely correct - in particular
// for ordered version.
void get_labeled_flags_of_one_type(int flag_size, const flag &type,
                                   vector<flag> &flag_list) {

  // cerr << "XXXXX" << endl;

  // Saving of the flags to a separate file and/or loading them
  string filename = get_filename_for_labeled_flags(flag_size, type);

  // umap.insert(make_pair("e", 2.718));

  if (g_labeled_flags_of_one_type_map.find(filename) ==
      g_labeled_flags_of_one_type_map.end()) {
#pragma omp critical(generating_special_labeled_flags)
    {
      if (g_labeled_flags_of_one_type_map.find(filename) ==
          g_labeled_flags_of_one_type_map.end()) {
        // cerr << "About to start loading " << filename << endl;
        if (load_flags_from_file(filename, flag_list)) {
          cerr << "Loaded " << flag_list.size() << " flags from " << filename
               << endl;
        } else {
          cerr << "Generating flags for file " << filename << endl;
          generate_labeled_flags_of_one_type(flag_size, type, flag_list);
        }

        g_labeled_flags_of_one_type_map.insert(make_pair(filename, flag_list));
      }
    }
  } else {
    // cerr << "Type already exists" << endl;
  }

  assert(g_labeled_flags_of_one_type_map.find(filename) !=
         g_labeled_flags_of_one_type_map.end());

  flag_list = g_labeled_flags_of_one_type_map[filename];
}

// vector<flag> g_types[V];
void create_types_of_size(int type_size) {
  if (g_types[type_size].size() != 0)
    return;

  // Look if there are already some flags with given types...
  for (int i = 0; i < (int)g_flags.size(); i++) {
    if (g_flags[i].size() == 0)
      continue;

    if (g_flags[i][0].labeled_vertices_cnt() == type_size) {
      flag type;
      g_flags[i][0].get_type_subflag(type);
      g_types[type_size].push_back(type);
    }
  }
  if (g_types[type_size].size() != 0)
    return;

  // if not, just label graphs of apropriate size
  g_types[type_size] = g_unlabeled_flags[type_size];
  generate_unlabeled_flags_of_size(type_size);

#pragma omp parallel for
  for (int i = 0; i < (int)g_types[type_size].size(); i++) {
    g_types[type_size][i].make_all_vertices_labeled();
    //
    cerr << "New type: " << g_types[type_size][i].print() << endl;
  }

  cerr << "Generated " << g_types[type_size].size() << " types of size "
       << type_size << endl;
}

int get_reduced_type_id(const flag &f) {
  int type_size = f.labeled_vertices_cnt();

  // cerr << "Testing size " << type_size << endl;

  for (int i = 0; i < (int)g_types[type_size].size(); i++) {
    if (f.have_same_type(g_types[type_size][i])) {
      return i;
    }
  }

  /*
  int found_type = -1;
  for (int i = 0; i < (int)g_types[type_size].size(); i++)
  {
      if (f.have_same_type(g_types[type_size][i]))
      {
          if (found_type == -1)
          {
              found_type=i;
          }
          else
          {
              cerr << "Type confusion! " << endl;
              cerr << f.print() << endl;
              cerr << g_types[type_size][i].print() << endl;
              cerr <<  g_types[type_size][found_type].print() << endl;
              exit(1);
          }
      }
  }
  return found_type;
  */
  return -1;
}

void get_types_of_size(int type_size, bool types_from_file) {
  //    create_types_of_size(type_size);
  //    return;

  if (!types_from_file) {
    create_types_of_size(type_size);
    cerr << g_types[type_size].size() << " types of size " << type_size
         << " generated" << endl;
    return;
  }

  string filename =
      filename_prefix() + "__type_" + to_string(type_size) + ".txt";

  if (!load_flags_from_file(filename, g_types[type_size])) {
    create_types_of_size(type_size);
    dump_flags_to_file(filename, g_types[type_size]);
  } else {
    cerr << g_types[type_size].size() << " types of size " << type_size
         << " loaded from file " << filename << endl;
  }
}

void generate_labeled_flags(int flag_size, int type_size, int verbose_output,
                            bool types_from_file) {
  // First make types
  get_types_of_size(type_size, types_from_file);

  // cerr << "B:EEEE  " << g_types[type_size].size() << endl;
  // cerr <<  "flag_size=" << flag_size << " type_size=" << type_size << endl;
  // cerr << "check_sensible_flags=" << check_sensible_flags << endl;

  int to_label_size = (int)g_unlabeled_flags[flag_size].size();

  //
  vector<vector<flag>> *flags;
  flags = new vector<vector<flag>>[to_label_size];

  int types = (int)g_types[type_size].size();

  // Avoid copying flags
  for (int i = 0; i < to_label_size; i++) {
    flags[i].resize(types);
  }

  cerr << "Generating labeled flags for Cauchy-Swartz in parallel " << endl;
#pragma omp parallel for ordered
  for (int i = 0; i < to_label_size; i++) {
    if (verbose_output)
      cerr << "Flag size " << flag_size << " type size " << type_size
           << " labeling " << i << "/" << to_label_size << endl;

    // cerr << "Type size " << type_size << " for " <<
    // g_unlabeled_flags[flag_size][i].print() << endl;
    //  Permutations for labeling... Could be improved to distinguish only the
    //  first type_size vertices.
    int mapping[flag_size];
    flag F;
    for (int j = 0; j < flag_size; j++)
      mapping[j] = j;

    // int found_good_type = 0;
    do {
      F.as_subflag(g_unlabeled_flags[flag_size][i], mapping, flag_size,
                   type_size); // taking a subflag
      // cerr << "Got one" << endl;
      int type_of_F = get_reduced_type_id(F);
      if (type_of_F == -1) {
        // We have each type only once -  3 3   1 1  2  and  3 3   1 2  1 are
        // different types but enough to use just one
        continue;
      }
      // found_good_type++;
      if (!g_already_in_known_flags(F, flags[i][type_of_F])) {
        flags[i][type_of_F].push_back(F);
      }
      // include_flag_in_list_if_new(F,flags[i],false);
    } while (std::next_permutation(mapping, mapping + flag_size));
    // cout << "Found good type " << found_good_type << endl;
  }

  // Todo: this might be also maybe paralelized - there is some paralelism
  // inside...
  for (int i = 0; i < to_label_size; i++) {
    for (int j = 0; j < (int)flags[i].size(); j++)
      for (int k = 0; k < (int)flags[i][j].size(); k++)
        include_flag_in_list_if_new(flags[i][j][k], g_flags);

    // merge_typed_flags_lists_parallel(g_flags, flags[i]);
  }

  for (int t = 0; t < types; t++) {
    int found_t = 0;
    for (int i = 0; i < to_label_size; i++) {
      found_t += (int)flags[i][t].size();
    }
    if (found_t == 0) {
      cerr << "Type " << t << " has no labeled flags" << endl;
    }
  }

  /*
      // Merging together to flags[0]
      cerr << "Merging flags  " << endl;
      for (int power = 1; power < to_label_size; power *= 2)
      {
  //        if (verbose_output)
          {
          //    cerr << "Trying " << power << " in " << to_label_size << endl;
          }
          for (int j = 0; j < to_label_size; j += 2*power)
          {
              if (j+power < to_label_size)
              {
                  //if (verbose_output)
                  {
                  //    cerr << "Merging  " << j << "<-" << j+power << "  of
  sizes " <<  flags[j].size() << " " << flags[j+power].size() << endl;
                  }
                  for (int t = 0; t < types; t++)
                  {
                      merge_flags_lists_parallel(flags[j][t],flags[j+power][t]);
                  }
              }
          }
      }

      // Merging to g_flags
      cerr << "Final merge " << g_flags.size() << " " << flags[0].size() <<
  endl; for (int t = 0; t < types; t++)
      {
          g_flags.push_back(flags[0][t]);
      }
  */
  // merge_typed_flags_lists_parallel(g_flags, flags[0]);

  // cerr << "delete[] " << endl;

  delete[] flags;

  // cerr << "done " << endl;
}

inline void add_g_to_flags_list_if_new(flag &g, vector<flag> &flag_list) {
  if (g_already_in_known_flags(g, flag_list))
    return;
  add_g_to_known_flags(g, flag_list);
}

void extensions_of_g_edges(flag &g, vector<flag> &flag_list) {
  // cerr << "Generating extensions of " << g.print() << endl;

#ifdef G_COLORED_EDGES
  try_color_edge(g, 0, 1, flag_list);
  return;
#endif
}

#ifdef G_COLORED_VERTICES
void try_color_vertex(flag &g, int u, vector<flag> &flag_list) {
  if (u >= g.m_vertices) {
#if defined(G_COLORED_EDGES) || defined(G_COLORED_3EDGES) ||                   \
    defined(G_3EDGES) || defined(G_4EDGES) || defined(G_ROOTED_3EDGES) ||      \
    defined(G_ROOTED_4EDGES)
    extensions_of_g_edges(g, flag_list);
    //        try_color_edge(g, 0,1,flag_list);
#else
    if (is_flag_forbidden(g))
      return;

    g.create_minlex_signature();

    if (!g_already_in_known_flags(g, flag_list))
      add_g_to_known_flags(g, flag_list);
#endif
    return;
  }

  // Color only uncolored vertices.....
  if (g.m_color_vertex[u] == 0) {
    for (int color = 1; color < COLORS_VERTICES; color++) {
      g.color_vertex(u, color);
      try_color_vertex(g, u + 1, flag_list);
    }
    g.color_vertex(u, 0);
  } else {
    try_color_vertex(g, u + 1, flag_list);
  }
}
#endif

void extensions_of_g(flag &g, vector<flag> &flag_list) {
#ifdef G_COLORED_VERTICES
  try_color_vertex(g, 0, flag_list);
#else
  extensions_of_g_edges(g, flag_list);
#endif
}

bool compare_flag_sizes(const flag &f1, const flag &f2) {
  return f1.m_vertices < f2.m_vertices;
}

void sort_flags_by_size(vector<flag> &flag_list) {
  sort(flag_list.begin(), flag_list.end(), compare_flag_sizes);
}

void load_forbidden() {
  string filename =
      filename_prefix(
          "forbidden",
          ProblemConfig::instance().get_forbidden_subgraphs_directory()) +
      ".txt";

  ifstream infile;
  infile.open(filename.c_str(), ifstream::in);
  if (!infile.good()) {
    cerr << "Failed opening file " << filename
         << " no forbidden structures are used" << endl;
    return;
  }

  FilteringIstream filteredinfile(infile);

  int flags_in_file = 0;

  flag h;
  while (h.load_from_stream(filteredinfile, -1, 0)) {
    flags_in_file++;
    // check if h is not already there
    if (!(g_already_in_known_flags(h, g_forbidden_subflags))) {
      g_forbidden_subflags.push_back(h);
      // cerr << "Loaded forbidden graph " << h.print() << endl;
    }
  }

  infile.close();

  // Idea is that it is faster to test smaller subflags than bigger
  // so lets first test smaller ones
  sort_flags_by_size(g_forbidden_subflags);

  cerr << "Loaded " << g_forbidden_subflags.size() << " forbidden graphs from "
       << filename << endl;
  if (flags_in_file != (int)g_forbidden_subflags.size()) {
    cerr << "File " << filename << " contains "
         << flags_in_file - g_forbidden_subflags.size() << " duplicate  graphs."
         << endl;
  }

  // for (int i = 0; i < g_forbidden_subflags.size(); i++)
  //{
  //     cerr << g_forbidden_subflags[i].print() << endl;
  // }
}

// The following function megres dest_flags with src_flags and clears src_flags
template <class data_type>
void merge_vectors_parallel(vector<data_type> &dest_vector,
                            vector<data_type> &src_vector) {
  vector<data_type> new_vector;
  new_vector.reserve(src_vector.size());

#pragma omp parallel for
  for (int k = 0; k < (int)src_vector.size(); k++) {
    bool already_exists = false;
    for (const auto &test : dest_vector) {
      if (test == src_vector[k]) {
        already_exists = true;
      }
    }
    if (already_exists == false) {
#pragma omp critical(merging_lists)
      {
        new_vector.push_back(src_vector[k]);
      }
    }
  }
  // #ifdef DONT_USE_C11
  dest_vector.insert(dest_vector.end(), new_vector.begin(), new_vector.end());
  src_vector.clear();
  vector<data_type>().swap(src_vector); // free the memory
                                        // #else
  //     dest_flags.insert(dest_flags.end(),
  //     std::make_move_iterator(new_flags.begin()),
  //     std::make_move_iterator(new_flags.end()) ); src_flags.clear();
  //     src_flags.shrink_to_fit(); // free the memory - does not work well on
  //     icc
  // #endif
}

/*
    void try_crossings_extension_of_g_by_adding_vertex(flag &g, vector<flag>
&flag_list)
{
    int n = g.m_vertices+1;
    int v = n-1;
    flag f;
    f.set_vertices_and_Theta(n,0);
    f.copy_from(g);
#ifdef CROSSING_EXTENSIONS_COLOR
    f.color_vertex(v, CROSSING_EXTENSIONS_COLOR);
#else
    #ifdef G_COLORED_VERTICES
    #error The color fo extension needs to be present
    #endif
#endif

    // Now we try to add the remaining crossings...
}
*/

void generate_all_types_containing_one_subtype(int big_type_size,
                                               const flag &type,
                                               vector<flag> &flag_list) {
  flag g;
  g.set_vertices_and_Theta(big_type_size, big_type_size);
  g.copy_from(type);

  extensions_of_g(g, flag_list);
}

void hydrate_flag_product_cache(int Kn, int idx_in_list, int typesize,
                                ofstream &outfile) {
  const flag &f = get_unlabeled_flags(Kn)[idx_in_list];
  int subset_split[V];
  int mapping_F1[V];
  int mapping_F2[V];
  flag F1, F2;

  // Try to find decomposition of vertices of f into 3 sets.
  // first is of size typesize and the other two are of sizes privatesize
  // Then typesize+privatesize creates one flag.

  vector<pair<FlagProduct, int>> found_pairs;

  int privatesize = (f.m_vertices - typesize) / 2;
  int totalsize = typesize + privatesize; // size of F1 and F2
  int label_f1 = typesize;
  int label_f2 = typesize + 1;

  for (int i = 0; i < typesize; i++)
    subset_split[i] = i;
  assert(typesize + 2 * privatesize <= V);
  for (int i = 0; i < privatesize; i++) {
    subset_split[typesize + i] = label_f1;
    subset_split[typesize + i + privatesize] =
        label_f2; // I don't know why the warning is here...
  }

  do {
    // std::cout << subset_split[0] << ' ' << subset_split[1] << ' ' <<
    // subset_split[2] << ' ' << subset_split[3] << '\n';

    // Construction of inverz mapping
    int inF1 = typesize;
    int inF2 = typesize;
    for (int i = 0; i < f.m_vertices; i++) {
      if (subset_split[i] < typesize) {
        mapping_F1[subset_split[i]] = i;
        mapping_F2[subset_split[i]] = i;
      }
      if (subset_split[i] == label_f1)
        mapping_F1[inF1++] = i;
      if (subset_split[i] == label_f2)
        mapping_F2[inF2++] = i;
    }

    F1.as_subflag(f, mapping_F1, totalsize, typesize);

    int typeID = get_flag_type_in_list(F1, get_flags());
    if (typeID == -2)
      continue;
#ifdef G_NOT_ALL_FLAGS_USED
    if (typeID == -1)
      continue;
#endif
    assert(typeID != -1);

    F2.as_subflag(f, mapping_F2, totalsize, typesize);

    int idF1 = find_flag_in_list(F1, get_flags()[typeID]);
    int idF2 = find_flag_in_list(F2, get_flags()[typeID]);
    if (idF1 > idF2)
      continue; // count only once when idF1 <= idF2
    if (idF1 == -1) {
#ifdef G_NOT_ALL_FLAGS_USED
      continue;
#endif
      cerr << "Missing " << F1.print() << endl;
    }
    assert(idF1 >= 0 && idF2 >= 0);

    FlagProduct flag_product{typeID + 1, idF1 + 1, idF2 + 1};

    bool known = false;
    for (int i = 0; i < (int)found_pairs.size(); i++) {
      if (found_pairs[i].first == flag_product) {
        found_pairs[i].second++;
        known = true;
        break;
      }
    }
    if (known == false) {
      found_pairs.push_back(make_pair(flag_product, 1));
    }

  } while (std::next_permutation(subset_split, subset_split + f.m_vertices));

  for (int i = 0; i < (int)found_pairs.size(); i++) {
    FlagProductPartialResult partial_result;
    partial_result.product = found_pairs[i].first;
    partial_result.count = found_pairs[i].second;

    outfile << partial_result.product << " " << partial_result.count << endl;
  }
}

void hydrate_projection_cache(flag type, int vertices, int new_vertices,
                              int remaining_labeled, ofstream &outfile) {
  vector<flag> flags;
  vector<flag> flags_smaller_type;
  flag type_copy = type;
  type_copy.set_Theta(remaining_labeled);
  flag reduced_type;
  type_copy.get_type_subflag(reduced_type);
  get_labeled_flags_of_one_type(vertices, type, flags);
  get_labeled_flags_of_one_type(new_vertices, reduced_type, flags_smaller_type);

  for (int i = 0; i < flags_smaller_type.size(); i++) {
    if (i % 1000 == 0) {
      cerr << "Calculating for " << i << " out of " << flags_smaller_type.size()
           << endl;
    }
    for (int j = 0; j < flags.size(); j++) {
      pair<int, int> proj_result = P_F1_IN_H(flags[j], flags_smaller_type[i]);
      if (proj_result.first == 0 or proj_result.second == 0)
        continue;
      outfile << i << " " << j << " " << proj_result.first << " "
              << proj_result.second << endl;
    }
  }
}

void hydrate_multiplication_cache(flag type, int size_left, int size_right,
                                  ofstream &outfile) {
  int n = size_left + size_right - type.m_vertices;

  vector<flag> flags_big;
  vector<flag> flags_left;
  vector<flag> flags_right;
  get_labeled_flags_of_one_type(n, type, flags_big);
  get_labeled_flags_of_one_type(size_left, type, flags_left);
  get_labeled_flags_of_one_type(size_right, type, flags_right);

  for (int i = 0; i < flags_big.size(); i++) {
    if (i % 1000 == 0) {
      cerr << "Calculating for " << i << " out of " << flags_big.size() << endl;
    }
    for (int j = 0; j < flags_left.size(); j++) {
      for (int k = 0; k < flags_right.size(); k++) {
        pair<int, int> mult_result =
            P_F1_F2_IN_labeled_H(flags_left[j], flags_right[k], flags_big[i]);
        if (mult_result.first == 0)
          continue;
        outfile << i << " " << j << " " << k << " " << mult_result.first << " "
                << mult_result.second << endl;
      }
    }
  }
}
