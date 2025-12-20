#pragma once
#include "config.hpp"
#include <string>
#include <utility>
#include <vector>

using namespace std;

// This is the main class that holds a flag.
class flag {
public:
  flag();
  flag(const std::string &str);
  flag(const flag &f);

  flag &operator=(const flag &f);
  // Nothing to do with isomorphisms.
  bool operator<(const flag &f) const;

  int m_vertices; // number of vertices
  int m_Theta;    // number of labeled vertices in case of not ordered. Binary
                  // string which are labeled in ordered
  int m_Theta_class; // Optional thing - same theta can have several classes -
                     // this should help with doing fancier rounding.
                     // m_Theta_class == 0 means any class. If not specified,
                     // any class is the default

  //	int m_id; // id of the original flag

#ifdef G_COLORED_VERTICES
  int m_color_vertex[V]; // color of a vertex (for bipartite complete flags?)
  int m_colored_vertices[COLORS_VERTICES];
#endif

#ifdef G_COLORED_EDGES
  int m_color_edge[V][V]; // adjacency matrix
  //	int m_colored_deg[V][COLORS_EDGES]; // number of edges from each vertex
  // of a particular color
  int m_colored_edges[COLORS_EDGES]; // number of edges of each color
//	vector<int> m_deg_sorted;
//	bool m_deg_sorted_valid;
#endif

  void make_all_vertices_labeled();

  bool is_labeled(const int v) const;

  // returns the number of labeled vertices
  int labeled_vertices_cnt() const;

  // First m_Theta vertices are considered labeled in that order
  void set_Theta(int Theta, int Theta_class = 0);

  int get_Theta();

  void set_vertices_and_Theta(int vertices, int Theta, int Theta_class = 0);

  // set the nymber of vertices and clears the graph
  void set_vertices(int vertices);

#ifdef G_COLORED_VERTICES
  void color_vertex(int u, int color);
#endif

#ifdef G_COLORED_EDGES
  void color_edge(int u, int v, int color);
#endif

#ifdef G_COLORED_VERTICES
  void permute_vertex_colors(const vector<int> &color_permutation);
#endif

#ifdef G_COLORED_EDGES
  // Makes a permutation of edges of this flag
  void permute_edge_colors(const vector<int> &color_permutation);
#endif

public:
  // This is the function to call in general program
  bool have_same_type(const flag &f) const;

  template <bool verbose_output = false>
  bool is_isomorphic_to_colorblind_colored_3edges(const flag &h) const;

  template <bool verbose_output = false>
  bool is_isomorphic_to_colorblind_colored_edges(const flag &h) const;

  template <bool verbose_output = false>
  bool is_isomorphic_to_colorblind_oriented_edges(const flag &h) const;

  template <bool verbose_output = false>
  bool is_isomorphic_to_colorblind_vertices(const flag &h) const;

  template <bool verbose_output = false>
  bool is_isomorphic_to_reverseblind_leftright_system(const flag &h) const;

  bool is_isomorphic_to_reverseblind_rotation_system(const flag &h) const;

  template <bool verbose_output = false>
  bool is_isomorphic_to(const flag &h) const;

  // perm is mapping vertices of h to vertices of this. pv is a candidate for
  // mapping of v, all vertices before v are already mapped.
  bool is_map_up_to_v_correct(int v, int pv, const int *perm,
                              const flag &h) const;

  // perm is saying how are vertices of h mapped to *this
  // it means color of x in h is the same as color of perm[x] in *this
  //
  //  v in h is mapped to perm[v] in *this
  //
  bool is_mapping_an_equality(int *perm, const flag &h) const;

  // Real isomorphism testing
  bool make_perm_isomorphic(int v, int *perm, bool *used, const flag &h) const;

  template <bool verbose_output = false>
  bool is_isomorphic_to_not_colorblind(const flag &h) const;

#if defined(G_3EDGES) || defined(G_4EDGES) || defined(G_ROOTED_3EDGES) ||      \
    defined(G_ROOTED_4EDGES) ||                                                \
    defined(G_MAYBE_ROOTED_KEDGES) // || (defined(G_COLORED_EDGES) &&
                                   // (G_COLORED_EDGES == 2))
  // perm is mapping vertices of h to vertices of this. pv is a candidate for
  // mapping of v, all vertices before v are already mapped.
  bool is_map_up_to_v_correct_notinduced(int v, int pv, const int *perm,
                                         const flag &h) const;

  // Real isomorphism testing
  bool make_perm_for_noninuced_subflags(int v, int *perm, bool *used,
                                        const flag &h) const;

  bool has_as_notinduced_subflag(const flag &h) const;
#else
  bool has_as_notinduced_subflag(const flag &h) const;
#endif

  bool is_identical_to(const flag &h, bool write_why_not = true) const;

  void create_minlex_signature_blind_3edges(const flag &h,
                                            string &lex_min) const;

  void create_minlex_signature_blind_edges(const flag &h,
                                           string &lex_min) const;

  void create_minlex_signature_oriented_edges(const flag &h,
                                              string &lex_min) const;

  void create_minlex_signature_blind_vertices(const flag &h,
                                              string &lex_min) const;

  void create_minlex_signature();

  void get_type_subflag(flag &f) const;

  // removes two vertices
  void remove_vertices(int u, int v);

  void remove_labeled_vertices();

  void as_subflag(const flag &h, const int *mapping, int vertices, int theta,
                  bool create_signature = true);

  void as_subflag_in_blowup(const flag &h, const int *mapping, int vertices,
                            int theta, bool create_signature = true);

  int e();

  int e3();

  bool contains_as_subflag_map_rest(const flag &g, int *mapping,
                                    int next_to_map, int first_available) const;

  // Count density of g as a subflags
  bool contains_as_subflag(const flag &g) const;

  void density_subflag_map_rest(const flag &g, int &total, int &good,
                                int *mapping, int next_to_map,
                                int first_available) const;

  // Count density of g as a subflags
  double density_subflag(const flag &g) const;

  int count_subflag(const flag &g) const;

  void count_labeled_copies_of_map_rest(const flag &g_labeled, int &good_maps,
                                        int *mapping, int *used,
                                        int next_to_map) const;

  int count_labeled_copies_of(const flag &g) const;

  void
  generate_subflags_of_size_n_rest(int n, vector<flag> &subflags, int *mapping,
                                   int next_to_map, int first_available,
                                   const vector<int> &available_to_map) const;

  // Could be written as a template with size of in_subflag - should be faster
  // an/or with arrays
  // if in_subrgaph is specified, then only subflags containing in_subrgaph
  // are generated
  int generate_subflags_of_size_n(
      int n, vector<flag> &subflags,
      const vector<int> &in_subflag = vector<int>()) const;

  // Copy data from g if the number of vertices of g is <= current number
  void copy_from(const flag &g);

  /*****************************************************************************************************************
   * Printing ****************/
  void print_for_human() const;

  template <typename T>
  string print_latex(bool use_label, const T &graph_label,
                     bool color_1_nonedge = false);

  string print(const string &delimeter = " ") const;

#ifdef G_COLORED_VERTICES
  void print_vertex_coloring(ostream &ostr, const string &delimiter) const;
#endif

#ifdef G_COLORED_EDGES
  void print_edge_coloring(ostream &ostr, const string &delimiter) const;
#endif

  /*****************************************************************************************************************
   * Loading ****************/
  void load_from_string(const char *str);

  bool load_from_stream(istream &stream, int assumed_vertices,
                        int assumed_theta);

  // Sequence of functions for checking if two flags have the same type //
private:
  // TODO: This sequence could be rewritten by using early check and no need
  // for creating copies of F - could be more effective.
  bool have_same_type_not_colorblind(const flag &f) const;
};

bool is_flag_forbidden(const flag &g, int verbose_output = 0);

class flag_coeff {
public:
  bool operator==(const flag_coeff &fc) const;

  flag g;
  double coefficient;
  void print(std::ostream &stream, bool skip_if_empty = true,
             string prefix = "") const;
};

std::ostream &operator<<(std::ostream &stream, const flag_coeff &fc);

int get_flag_type_in_list(const flag &f, vector<vector<flag>> &flags);

int find_flag_in_list(const flag &f, const vector<flag> &fv);

int find_flag_in_list_nonparalel(const flag &f, const vector<flag> &fv);

void include_flag_in_list(const flag &f, vector<vector<flag>> &flags,
                          int id = -3, bool ignore_reduced_types = true);

void include_flag_in_list_if_new(const flag &f, vector<vector<flag>> &flags,
                                 bool ignore_reduced_types = true);

pair<int, int> P_F1_IN_H(const flag &F1, const flag &H, bool density = true);

pair<int, int> P_F1_F2_IN_labeled_H(const flag &F1, const flag &F2,
                                    const flag &H);
