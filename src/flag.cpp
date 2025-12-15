#include "flag.hpp"
#include "config.hpp"
#include "flag-generator.hpp"
#include "utils.hpp"
#include <algorithm>
#include <array>
#include <assert.h>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;

// adds more precision in specifying the blow-up
static array<int, V> g_blow_up_color_edges = [] {
  static array<int, V> arr;
  for (int i = 0; i < V; i++)
    arr[i] = G_BLOW_UP_COLOR_EDGES;
  return arr;
}();

flag::flag() {
  m_Theta = 0;
  m_Theta_class = 0;
  m_vertices = 0;
}

flag::flag(const string &str) : flag() {
  stringstream ss(str);
  load_from_stream(ss, -1, -1);
}

void flag::make_all_vertices_labeled() { m_Theta = m_vertices; }

bool flag::is_labeled(const int v) const {
  assert(v < m_vertices);
  return v < m_Theta;
}

// returns the number of labeled vertices
int flag::labeled_vertices_cnt() const { return m_Theta; }

// First m_Theta vertices are considered labeled in that order
void flag::set_Theta(int Theta, int Theta_class) {
  // if !G_ORDERED_VERTICES, then the labeled vertices are first. Not true for
  // ordered.
  // cerr << m_vertices << " " << Theta << endl;
  if (Theta > m_vertices) {
    cerr << "FAIL: Attempting to load grap on " << m_vertices
         << " vertices with Theta=" << Theta << endl;
  }
  assert(Theta <= m_vertices);
  m_Theta = Theta;
  m_Theta_class = Theta_class;
}

int flag::get_Theta() { return m_Theta; }

void flag::set_vertices_and_Theta(int vertices, int Theta, int Theta_class) {
  set_vertices(vertices);
  set_Theta(Theta, Theta_class);
}

// set the nymber of vertices and clears the graph
void flag::set_vertices(int vertices) {
  assert(vertices <= V);
  m_vertices = vertices;

#ifdef G_COLORED_VERTICES
  for (int u = 0; u < V; u++) {
    m_color_vertex[u] = 0;
  }
  m_colored_vertices[0] = m_vertices;
  for (int c = 1; c < COLORS_VERTICES; c++)
    m_colored_vertices[c] = 0;
#endif

#ifdef G_COLORED_EDGES
  for (int u = 0; u < V; u++) {
    for (int v = 0; v < V; v++) {
      m_color_edge[u][v] = 0;
    }
  }

  m_colored_edges[0] = e();
  for (int c = 1; c < COLORS_EDGES; c++)
    m_colored_edges[c] = 0;
#endif
}

#ifdef G_COLORED_VERTICES
void flag::color_vertex(int u, int color) {
  if (color >= COLORS_VERTICES) {
    cerr << "Coloring vertex " << u << " with color " << color
         << " that is more than " << COLORS_VERTICES - 1 << endl;
    cerr << "we have " << m_vertices << " vertices and " << m_Theta
         << " m_Theta" << endl;
  }

  assert(u < m_vertices);
  assert(color >= 0);
  assert(color < COLORS_VERTICES);

  m_colored_vertices[m_color_vertex[u]]--;

  m_color_vertex[u] = color;

  m_colored_vertices[color]++;
}
#endif

#ifdef G_COLORED_EDGES
void flag::color_edge(int u, int v, int color) {

#if defined(G_ORIENTED_EDGES) || defined(G_REMOVE_ORIENTATION_WHEN_LOADING)
  if (color < 0) {
    color_edge(v, u, -color);
    return;
  }
#endif
  assert(u != v);
  assert(u < m_vertices);
  assert(v < m_vertices);
  assert(0 <= u);
  assert(0 <= v);
  assert(color >= 0);
  assert(color < COLORS_EDGES);

  int old_color = m_color_edge[u][v];
#ifdef G_ORIENTED_EDGES
  if (old_color < 0)
    old_color = -old_color;
#endif

  m_colored_edges[old_color]--;

#ifdef G_ORIENTED_EDGES
  if (color >
      G_ORIENTED_EDGES_UNORIENTED_COLORS) // defaul is 1, it is for not
                                          // oriented edges, non-edges....
  {
    m_color_edge[u][v] = color;
    m_color_edge[v][u] = -color;
  } else {
    m_color_edge[u][v] = m_color_edge[v][u] = color;
  }
#else
  m_color_edge[u][v] = m_color_edge[v][u] = color;
#endif

  m_colored_edges[color]++;
}
#endif

#ifdef G_COLORED_VERTICES
//    void permute_vertex_colors(int color_permutation[])
void flag::permute_vertex_colors(const vector<int> &color_permutation) {
  for (int v = 0; v < m_vertices; v++) {
    color_vertex(v, color_permutation[m_color_vertex[v]]);
  }
}
#endif

#ifdef G_COLORED_EDGES
// Makes a permutation of edges of this flag
// Maybe be rewritten for higher effectivity
void flag::permute_edge_colors(const vector<int> &color_permutation) {
  assert(color_permutation.size() >= COLORS_EDGES);

  for (int u = 0; u < m_vertices - 1; u++)
    for (int v = u + 1; v < m_vertices; v++) {
#ifdef G_ORIENTED_EDGES
      if (m_color_edge[u][v] < 0) {
        color_edge(u, v, -color_permutation[-m_color_edge[u][v]]);
      } else {
        color_edge(u, v, color_permutation[m_color_edge[u][v]]);
      }
#else
      color_edge(u, v, color_permutation[m_color_edge[u][v]]);
#endif
    }
}
#endif

bool flag::have_same_type_colorblind_colored_3edges(const flag &h) const {
  return have_same_type_colorblind_colored_edges(h);
}

bool flag::have_same_type_colorblind_colored_edges(const flag &h) const {
  return have_same_type_colorblind_oriented_edges(h);
}

bool flag::have_same_type_colorblind_oriented_edges(const flag &h) const {
  return have_same_type_colorblind_vertices(h);
}

bool flag::have_same_type_colorblind_vertices(const flag &h) const {
  return have_same_type_leftright_blind(h);
}

bool flag::have_same_type_leftright_blind(const flag &h) const {
  return have_same_type_rotation_reverse_blind(h);
}

bool flag::have_same_type_rotation_reverse_blind(const flag &h) const {
  return have_same_type_not_colorblind(h);
}

bool flag::have_same_type_not_colorblind(const flag &f) const {
  // if (m_Theta != f.m_Theta) return false;
  // if (m_Theta == 0) return true;

#ifdef G_COLORED_VERTICES
  // same vertex color
  for (int u = 0; u < m_Theta; u++) {
#ifdef G_USING_ZERO_AS_ANY_COLOR
    if (m_color_vertex[u] == 0 || f.m_color_vertex[u] == 0)
      continue;
#endif
    if (m_color_vertex[u] != f.m_color_vertex[u])
      return false;
  }
#endif

#ifdef G_COLORED_EDGES
  assert(m_Theta < V);
  for (int u = 0; u < m_Theta; u++)
    for (int v = u + 1; v < m_Theta; v++) {
#ifdef G_USING_ZERO_AS_ANY_COLOR
      if (m_color_edge[u][v] == 0 || f.m_color_edge[u][v] == 0)
        continue;
#endif
      if (m_color_edge[u][v] != f.m_color_edge[u][v])
        return false;
    }
#endif

  return true;
}

bool flag::have_same_type(const flag &f) const {
  if (m_Theta_class != 0 && f.m_Theta_class != 0 &&
      m_Theta_class != f.m_Theta_class)
    return false;
  if (m_Theta == 0 && f.m_Theta == 0)
    return true;
  if (labeled_vertices_cnt() != f.labeled_vertices_cnt())
    return false;

  return have_same_type_colorblind_colored_3edges(f);
}

template <bool verbose_output>
bool flag::is_isomorphic_to_colorblind_colored_3edges(const flag &h) const {
  return is_isomorphic_to_colorblind_colored_edges<verbose_output>(h);
}

template <bool verbose_output>
bool flag::is_isomorphic_to_colorblind_colored_edges(const flag &h) const {
  return is_isomorphic_to_colorblind_oriented_edges<verbose_output>(h);
}
template bool
flag::is_isomorphic_to_colorblind_colored_edges(const flag &h) const;

template <bool verbose_output>
bool flag::is_isomorphic_to_colorblind_oriented_edges(const flag &h) const {
  return is_isomorphic_to_colorblind_vertices<verbose_output>(h);
}

template <bool verbose_output>
bool flag::is_isomorphic_to_colorblind_vertices(const flag &h) const {
  return is_isomorphic_to_reverseblind_leftright_system<verbose_output>(h);
}

template <bool verbose_output>
bool flag::is_isomorphic_to_reverseblind_leftright_system(const flag &h) const {
  return is_isomorphic_to_reverseblind_rotation_system(h);
}

bool flag::is_isomorphic_to_reverseblind_rotation_system(const flag &h) const {
  return is_isomorphic_to_not_colorblind(h);
}

template <bool verbose_output>
bool flag::is_isomorphic_to(const flag &h) const {
  // Just a quick kill
  if (m_vertices != h.m_vertices) {
    if (verbose_output)
      cerr << "The number of vertices is different" << endl;
    return false;
  }
  if (labeled_vertices_cnt() != h.labeled_vertices_cnt()) {
    if (verbose_output)
      cerr << "The number of labeled vertices is different" << endl;
    return false;
  }

  if (m_Theta_class != 0 && h.m_Theta_class != 0 &&
      m_Theta_class != h.m_Theta_class) {
    if (verbose_output)
      cerr << "The m_Theta_class is different" << endl;
    return false;
  }

  // cerr << "Testing " << h.print() << " and " << print() << endl;

  return is_isomorphic_to_colorblind_colored_3edges<verbose_output>(h);
}
template bool flag::is_isomorphic_to<false>(const flag &h) const;

// perm is mapping vertices of h to vertices of this. pv is a candidate for
// mapping of v, all vertices before v are already mapped.
bool flag::is_map_up_to_v_correct(int v, int pv, const int *perm,
                                  const flag &h) const {

#ifdef G_COLORED_VERTICES
  // 0 as a joker color
  if ((m_color_vertex[pv] != h.m_color_vertex[v]) &&
      (m_color_vertex[pv] != 0) && (0 != h.m_color_vertex[v]))
    return false;
#endif

#ifdef G_COLORED_EDGES
  // check edges, 0 as joker
  for (int u = 0; u < v; u++)
    if (m_color_edge[perm[u]][pv] != h.m_color_edge[u][v] &&
        m_color_edge[perm[u]][pv] != 0 && h.m_color_edge[u][v] != 0)
      return false;
#endif

  return true;
}

// perm is saying how are vertices of h mapped to *this
// it means color of x in h is the same as color of perm[x] in *this
//
//  v in h is mapped to perm[v] in *this
//
bool flag::is_mapping_an_equality(int *perm, const flag &h) const {
#ifdef G_COLORED_EDGES
  //	check coloring
  for (int u = 0; u < m_vertices; u++) {
    for (int x = u + 1; x < m_vertices; x++) {
      if (m_color_edge[perm[u]][perm[x]] == 0 || h.m_color_edge[u][x] == 0)
        continue;
      if (m_color_edge[perm[u]][perm[x]] != h.m_color_edge[u][x])
        return false;
    }
  }
#endif

#ifdef G_COLORED_VERTICES
  for (int u = 0; u < m_vertices; u++) {
    if (m_color_vertex[perm[u]] == 0 || h.m_color_vertex[u] == 0)
      continue;
    if (m_color_vertex[perm[u]] != h.m_color_vertex[u])
      return false;
  }
#endif

  return true;
}

// Real isomorphism testing
bool flag::make_perm_isomorphic(int v, int *perm, bool *used,
                                const flag &h) const {
  // Test if the permutation is isomorphism
  if (v >= m_vertices) {
#ifndef G_BE_BRAVE_AND_IGNORE_SAFETY_ASSERTS
    assert(is_mapping_an_equality(perm, h) == true);
#endif

    /*
    #ifdef G_LEFTRIGHT
                            // Check crossings
                            // now crossing of pairs of edges
                            for (int u1 = 0; u1 < m_vertices; u1++)
                    for (int u2 = u1+1; u2 < m_vertices; u2++)
                    {
                        bool test_normal = true;
                        bool test_flip = true;
                        for (int v1 = 0; v1 < m_vertices; v1++)
                        {
                            if (u1 == v1 || u2 == v1) continue;

                            if (m_leftright[perm[u1]][perm[u2]][perm[v1]] !=
    h.m_leftright[u1][u2][v1]) test_normal=false; if
    (m_leftright[perm[u1]][perm[u2]][perm[v1]] !=
    -1*h.m_leftright[u1][u2][v1]) test_flip=false;
                        }
                        if (test_normal == false && test_flip == false)
                        {
                            //cerr << "Help needed" << endl;
                            return false;
                        }
                    }
    #endif
     */

    return true;
  }

  for (int pv = m_Theta; pv < m_vertices; pv++) {
    // check if used
    if (used[pv])
      continue;

    perm[v] = pv;
    if (!is_map_up_to_v_correct(v, pv, perm, h))
      continue;

    used[pv] = true;

    if (make_perm_isomorphic(v + 1, perm, used, h))
      return true;

    used[pv] = false;
  }

  return false;
}

template <bool verbose_output>
bool flag::is_isomorphic_to_not_colorblind(const flag &h) const {
  // if (labeled_vertices_cnt() > 0)
  //     cerr << "Starting testing not colorblind " << print() <<" and "<<
  //     h.print() << endl;

  // cerr << endl << "Testing " << endl << h.print() << endl << print() <<
  // endl;

  // quick check
  if (m_vertices != h.m_vertices) {
    // if (labeled_vertices_cnt() > 0)
    // cerr << "Wrong number of vertices " << endl;
    return false;
  }
  if (labeled_vertices_cnt() != h.labeled_vertices_cnt()) {
    // if (labeled_vertices_cnt() > 0)
    // cerr << "Wrong number of labeled vertices " << labeled_vertices_cnt()
    // << " " << h.labeled_vertices_cnt() << endl;
    return false;
  }

#ifdef G_COLORED_VERTICES
  for (int c = 1; c < COLORS_VERTICES; c++) {
    // cerr << "c=" << c << "   "  << m_colored_vertices[c] << " "  <<
    // m_colored_vertices[0] << " " << h.m_colored_vertices[c] << " "  <<
    // h.m_colored_vertices[0] << endl;

    if (m_colored_vertices[c] >
        h.m_colored_vertices[c] + h.m_colored_vertices[0])
      return false;
    if (m_colored_vertices[c] + m_colored_vertices[0] < h.m_colored_vertices[c])
      return false;
  }
#endif

#ifdef G_COLORED_EDGES
  for (int c = 0; c < COLORS_EDGES; c++) {
    // if (m_vertices==0)
    //     cerr << c << "  " << m_colored_edges[c] << " " <<
    //     h.m_colored_edges[c] << endl;
    if (m_colored_edges[c] > h.m_colored_edges[c] + h.m_colored_edges[0])
      return false;
    if (m_colored_edges[c] + h.m_colored_edges[0] < h.m_colored_edges[c])
      return false;
  }
#endif

  //        if (m_vertices == 0)
  //            cerr << "Testing not colorblind " << print() <<" and "<<
  //            h.print() << endl;

  // cerr << "Testing not colorblind " << print() <<" and "<< h.print() <<
  // endl;

  // Check for same type
  if (!have_same_type_not_colorblind(h))
    return false;

  // try check all permutations --- WASTING HERE!!!!!!!!
  int perm[V];
  bool used[V];

  assert(0 <= m_Theta);
  assert(m_Theta <= V);

  // Theta is preserved by the mapping
  for (int i = 0; i < m_Theta; i++) {
    perm[i] = i;
    used[i] = true;
  }

  for (int i = m_Theta; i < V; i++) {
    perm[i] = 0;
    used[i] = false;
  }

  /*
  #ifdef G_ROTATION_SYSTEM
          bool new_test = use_rotation_system_for_isomorphis(h);
          bool old_test = make_perm_isomorphic(m_Theta, perm, used, h);

          if (new_test != old_test)
          {
              cerr << "h=" << h.print() << endl;
              cerr << "g=" << print() << endl;
              cerr << "new=" << new_test << endl;
              cerr << "old=" << old_test << endl;
          }

          //assert(new_test == old_test);

          return new_test;
  #endif
  */
  bool rv = make_perm_isomorphic(m_Theta, perm, used, h);
  if (verbose_output == true && rv == true) {
    cerr << "Isomorphic with permutation ";
    for (int i = 0; i < m_vertices; i++) {
      cerr << i << "->" << perm[i] << " ";
    }
    cerr << endl;
    cerr << "This maps the second graph to the first one" << endl;
  }
  return rv;
}

#if defined(G_3EDGES) || defined(G_4EDGES) || defined(G_ROOTED_3EDGES) ||      \
    defined(G_ROOTED_4EDGES) ||                                                \
    defined(G_MAYBE_ROOTED_KEDGES) // || (defined(G_COLORED_EDGES) &&
                                   // (G_COLORED_EDGES == 2))
// perm is mapping vertices of h to vertices of this. pv is a candidate for
// mapping of v, all vertices before v are already mapped.
bool flag::is_map_up_to_v_correct_notinduced(int v, int pv, const int *perm,
                                             const flag &h) const {
  // cerr << "Testing " << v << " " << pv << endl;

#ifdef G_COLORED_VERTICES
  if ((m_color_vertex[pv] != h.m_color_vertex[v]) &&
      (m_color_vertex[pv] != 0) && (0 != h.m_color_vertex[v]))
    return false;
#endif

  return true;
}

// Real isomorphism testing
bool flag::make_perm_for_noninuced_subflags(int v, int *perm, bool *used,
                                            const flag &h) const {

  // Test if the permutation is isomorphism
  if (v >= h.m_vertices) {
    return true;
  }

  for (int pv = m_Theta; pv < m_vertices; pv++) {

    // check if used
    if (used[pv])
      continue;

    // cerr << "Trying v="<< v <<  " to " << pv << endl;
    perm[v] = pv;
    if (!is_map_up_to_v_correct_notinduced(v, pv, perm, h))
      continue;

    used[pv] = true;

    if (make_perm_for_noninuced_subflags(v + 1, perm, used, h))
      return true;

    used[pv] = false;
  }

  return false;
}

bool flag::has_as_notinduced_subflag(const flag &h) const {
  // cerr << "Frist check " << h.print() << " in " << print() << endl;

  // YUFEI

  // YUFEI

  // Check for same type
  if (!have_same_type(h))
    return false;

  // cerr << "Testing " << h.print() << " in " << print() << endl;

  // try check all permutations --- WASTING HERE!!!!!!!!
  int perm[V];
  bool used[V];

  // Theta is preserved by the mapping
  for (int i = 0; i < m_Theta; i++) {
    perm[i] = i;
    used[i] = true;
  }

  for (int i = m_Theta; i < V; i++) {
    perm[i] = 0;
    used[i] = false;
  }

  return make_perm_for_noninuced_subflags(m_Theta, perm, used, h);
}
#else
bool flag::has_as_notinduced_subflag(const flag &h) const {
  // Not implemented for other
  assert(0);
  return false;
}
#endif

bool flag::is_identical_to(const flag &h, bool write_why_not) const {
  if (m_vertices != h.m_vertices) {
    if (write_why_not)
      cerr << "Number of vertices is different: " << m_vertices << " vs "
           << h.m_vertices << endl;
    return false;
  }

  if (m_Theta != h.m_Theta) {
    if (write_why_not)
      cerr << "Thetas are different: " << m_Theta << " vs " << h.m_Theta
           << endl;
    return false;
  }

  int perm[] = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
                11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
  assert(m_vertices < 20);
  return is_mapping_an_equality(perm, h);
}

void flag::create_minlex_signature_blind_3edges(const flag &h,
                                                string &lex_min) const {
  string Hprint = h.print("");
  if (Hprint.compare(lex_min) < 0) {
    lex_min = Hprint;
  }
}

void flag::create_minlex_signature_blind_edges(const flag &h,
                                               string &lex_min) const {
  create_minlex_signature_blind_3edges(h, lex_min);
}

void flag::create_minlex_signature_oriented_edges(const flag &h,
                                                  string &lex_min) const {
  create_minlex_signature_blind_edges(h, lex_min);
}

void flag::create_minlex_signature_blind_vertices(const flag &h,
                                                  string &lex_min) const {
  create_minlex_signature_oriented_edges(h, lex_min);
}

void flag::create_minlex_signature() {}

void flag::get_type_subflag(flag &f) const {
  const int mapping[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
  f.as_subflag(*this, mapping, m_Theta, m_Theta);
  f.m_Theta_class = m_Theta_class;
}

// removes two vertices
void flag::remove_vertices(int u, int v) {
  int mapping[V];

  int id = 0;
  for (int i = 0; i < m_vertices; i++) {
    if (i == u || i == v)
      continue;
    mapping[id++] = i;
  }

  flag h;
  h.as_subflag(*this, mapping, m_vertices - 2, 0);

  *this = h;
}

void flag::remove_labeled_vertices() {
  if (m_Theta == 0)
    return;

  int mapping[V];

  for (int i = 0; i < m_vertices - m_Theta; i++) {
    mapping[i] = i + m_Theta;
  }

  flag h;
  h.as_subflag(*this, mapping, m_vertices - m_Theta, 0);

  *this = h;
}

void flag::as_subflag(const flag &h, const int *mapping, int vertices,
                      int theta, bool create_signature) {
  set_vertices_and_Theta(vertices, theta);

#ifdef G_COLORED_EDGES
  for (int x = 0; x < m_vertices - 1; x++)
    for (int y = x + 1; y < m_vertices; y++) {
      color_edge(x, y, h.m_color_edge[mapping[x]][mapping[y]]);
    }
#endif

#ifdef G_COLORED_VERTICES
  // copy vertices
  for (int x = 0; x < m_vertices; x++)
    color_vertex(x, h.m_color_vertex[mapping[x]]);
#endif
}

void flag::as_subflag_in_blowup(const flag &h, const int *mapping, int vertices,
                                int theta, bool create_signature) {
  set_vertices_and_Theta(vertices, theta);

#ifdef G_COLORED_EDGES
  for (int x = 0; x < m_vertices - 1; x++)
    for (int y = x + 1; y < m_vertices; y++) {
      if (mapping[x] != mapping[y]) {
        color_edge(x, y, h.m_color_edge[mapping[x]][mapping[y]]);
      } else {
        // color_edge(x,y,G_BLOW_UP_COLOR_EDGES);
        color_edge(x, y, g_blow_up_color_edges[mapping[x]]);
      }
    }
#endif

#ifdef G_COLORED_VERTICES
  // copy vertices
  for (int x = 0; x < m_vertices; x++)
    color_vertex(x, h.m_color_vertex[mapping[x]]);
#endif
}

int flag::e() { return (m_vertices * (m_vertices - 1)) / 2; }

int flag::e3() {
  return (m_vertices * (m_vertices - 1) * (m_vertices - 2)) / 6;
}

bool flag::contains_as_subflag_map_rest(const flag &g, int *mapping,
                                        int next_to_map,
                                        int first_available) const {
  if (next_to_map >= g.m_vertices) {
    flag h;
    h.as_subflag(*this, mapping, g.m_vertices, g.m_Theta);

    // This may happen!!
    //  assert(g.m_Theta == 0);

    // cerr << "Testing " << h.print() << " and " << g.print() << endl;

    if (h.is_isomorphic_to(g)) {
      return true;
    }

    return false;
  }

  if (first_available < g.m_Theta) {
    mapping[next_to_map] = first_available;
    if (contains_as_subflag_map_rest(g, mapping, next_to_map + 1,
                                     first_available + 1))
      return true;
  } else {
    for (int i = first_available; i < m_vertices; i++) {
      mapping[next_to_map] = i;
      if (contains_as_subflag_map_rest(g, mapping, next_to_map + 1, i + 1))
        return true;
    }
  }

  return false;
}

// Count density of g as a subflags
bool flag::contains_as_subflag(const flag &g) const {
  if (g.m_vertices > m_vertices)
    return false;

  int mapping[g.m_vertices];

  return contains_as_subflag_map_rest(g, mapping, 0, 0);
}

void flag::density_subflag_map_rest(const flag &g, int &total, int &good,
                                    int *mapping, int next_to_map,
                                    int first_available) const {
  if (next_to_map >= g.m_vertices) {
    flag h;
    h.as_subflag(*this, mapping, g.m_vertices, g.m_Theta);

    assert(g.m_Theta == 0);

    if (h.is_isomorphic_to(g))
      good++;

    total++;
    return;
  }

  for (int i = first_available; i < m_vertices; i++) {
    mapping[next_to_map] = i;
    density_subflag_map_rest(g, total, good, mapping, next_to_map + 1, i + 1);
  }
}

// Count density of g as a subflags
double flag::density_subflag(const flag &g) const {
  if (g.m_vertices > m_vertices)
    return 0;

  int mapping[g.m_vertices];

  int total = 0;
  int good = 0;

  density_subflag_map_rest(g, total, good, mapping, 0, 0);

  return (double)good / (double)total;
}

int flag::count_subflag(const flag &g) const {
  if (g.m_vertices > m_vertices)
    return 0;

  int mapping[g.m_vertices];

  int total = 0;
  int good = 0;

  density_subflag_map_rest(g, total, good, mapping, 0, 0);

  return good;
}

void flag::count_labeled_copies_of_map_rest(const flag &g_labeled,
                                            int &good_maps, int *mapping,
                                            int *used, int next_to_map) const {
  if (next_to_map >= g_labeled.m_vertices) {
    flag h;
    h.as_subflag(*this, mapping, g_labeled.m_vertices, g_labeled.m_Theta);

    assert(g_labeled.m_Theta == g_labeled.m_vertices);

    if (h.have_same_type(g_labeled))
      good_maps++;

    return;
  }

  for (int i = 0; i < m_vertices; i++) {
    if (used[i])
      continue;
    mapping[next_to_map] = i;
    used[i] = 1;
    count_labeled_copies_of_map_rest(g_labeled, good_maps, mapping, used,
                                     next_to_map + 1);
    used[i] = 0;
  }
}

int flag::count_labeled_copies_of(const flag &g) const {
  if (g.m_vertices > m_vertices)
    return 0;

  // This does not work correctly for unlabeled graphs
  assert(g.m_Theta == 0 && m_Theta == 0);

  // We use have_same_type to check correct mappings...
  flag g_labeled = g;
  g_labeled.m_Theta = g.m_vertices;

  int mapping[g.m_vertices];
  int used[m_vertices];
  for (int i = 0; i < m_vertices; i++) {
    used[i] = 0;
  }

  int good_maps = 0;

  count_labeled_copies_of_map_rest(g_labeled, good_maps, mapping, used, 0);

  return good_maps;
}

void flag::generate_subflags_of_size_n_rest(
    int n, vector<flag> &subflags, int *mapping, int next_to_map,
    int first_available, const vector<int> &available_to_map) const {
  // cerr << "Runnig n=" << n << " mapping=" << mapping[0] << "," <<
  // mapping[1] << "," << mapping[2] << " next_to_map=" << next_to_map <<
  // endl;

  if (next_to_map >= n) {
    flag h;
    h.as_subflag(*this, mapping, n, 0);
    bool h_is_new = true;
    for (int i = 0; i < (int)subflags.size(); ++i) {
      if (h.is_isomorphic_to(subflags[i])) {
        h_is_new = false;
        break;
      }
    }
    if (h_is_new)
      subflags.push_back(h);
    return;
  }

  // cerr << "first_available=" << first_available  << " (n-next_to_map)=" <<
  // (n-next_to_map) << endl; cerr <<
  // "(int)available_to_map.size()-(n-next_to_map)=" <<
  // (int)available_to_map.size()-(n-next_to_map) << endl;
  for (int i = first_available;
       i < (int)available_to_map.size() - (n - next_to_map - 1); i++) {
    // cerr << "i=" << i << endl;
    mapping[next_to_map] = available_to_map[i];
    generate_subflags_of_size_n_rest(n, subflags, mapping, next_to_map + 1,
                                     i + 1, available_to_map);
  }
}

// Could be written as a template with size of in_subflag - should be faster
// an/or with arrays
// if in_subrgaph is specified, then only subflags containing in_subrgaph are
// generated
int flag::generate_subflags_of_size_n(int n, vector<flag> &subflags,
                                      const vector<int> &in_subflag) const {
  if (n > m_vertices)
    return 0;

  // This does not work correctly for unlabeled graphs
  assert(m_Theta == 0);

  int mapping[n];

  vector<int> available_to_map;
  for (int i = 0; i < m_vertices; i++)
    if (find(in_subflag.begin(), in_subflag.end(), i) == in_subflag.end())
      available_to_map.push_back(i);

  for (int i = 0; i < (int)in_subflag.size(); i++)
    mapping[i] = in_subflag[i];

  // cerr << "Working on it... " << available_to_map.size() << endl;
  // cerr << "Starting with " << (int)in_subflag.size() << endl;

  generate_subflags_of_size_n_rest(n, subflags, mapping, (int)in_subflag.size(),
                                   0, available_to_map);

  return (int)subflags.size();
}

// Copy data from g if the number of vertices of g is <= current number
void flag::copy_from(const flag &g) {
  assert(m_vertices >= g.m_vertices);

#ifdef G_COLORED_VERTICES
  // copy vertices
  for (int x = 0; x < g.m_vertices; x++)
    color_vertex(x, g.m_color_vertex[x]);
#endif

#ifdef G_COLORED_EDGES
  for (int x = 0; x < g.m_vertices - 1; x++)
    for (int y = x + 1; y < g.m_vertices; y++) {
      color_edge(x, y, g.m_color_edge[x][y]);
    }
#endif
}

/*****************************************************************************************************************
 * Printing ****************/
void flag::print_for_human() const {
  // Prints the number of vertices, zero as not type and then part of the
  // diagonal matrix
  cout << m_vertices << " " << m_Theta << " ";

#ifdef G_COLORED_VERTICES
  print_vertex_coloring(cout, " ");
  cout << " ";
#endif

#ifdef G_COLORED_EDGES
  print_edge_coloring(cout, " ");
  cout << " ";
#endif

  cout << endl;
}

template <typename T>
string flag::print_latex(bool use_label, const T &graph_label,
                         bool color_1_nonedge) {
  if (m_vertices == 0) {
    return "";
  }

  stringstream ss;

  // Prints the number of vertices, zero as not type and then part of the
  // diagonal matrix
  ss << " %\n\\begin{tikzpicture}[flag_pic]";
  ss << "\\outercycle{" << m_vertices << "}{" << m_Theta << "}\n";

  const string edgestr = "--";

#ifdef G_COLORED_EDGES
  for (int u = 0; u < m_vertices; u++) {
    for (int v = u + 1; v < m_vertices; v++) {
#ifdef G_ORIENTED_EDGES
      if (abs(m_color_edge[u][v]) <= G_ORIENTED_EDGES_UNORIENTED_COLORS)
        ss << "\\draw[edge_color" << m_color_edge[u][v] << "] (x" << u << ")"
           << edgestr << "(x" << v << ");";
      else if (m_color_edge[u][v] > 0)
        ss << "\\draw[edge_color" << m_color_edge[u][v] << ", -latex] (x" << u
           << ")" << edgestr << "(x" << v << ");";
      else
        ss << "\\draw[edge_color" << -m_color_edge[u][v] << ", latex-] (x" << u
           << ")" << edgestr << "(x" << v << ");";
#else
      ss << "\\draw[edge_color" << m_color_edge[u][v] << "] (x" << u << ")"
         << edgestr << "(x" << v << ");";
#endif
    }
    ss << "  ";
  }
  ss << endl;
#else // Edges have no colors
#endif

#ifdef G_COLORED_VERTICES
  for (int i = 0; i < m_vertices; i++) {
    if (i < m_Theta)
      ss << "\\draw (x" << i << ") node[labeled_vertex_color"
         << m_color_vertex[i] << "]{};";
    else
      ss << "\\draw (x" << i << ") node[vertex_color" << m_color_vertex[i]
         << "]{};";
  }
  ss << endl;
#else
  for (int i = 0; i < m_vertices; i++) {
    if (i < m_Theta)
      ss << "\\draw (x" << i << ") node[labeled_vertex]{};";
    else
      ss << "\\draw (x" << i << ") node[unlabeled_vertex]{};";
  }
  ss << endl;
#endif
  for (int i = 0; i < m_vertices; i++) {
    ss << "\\labelvertex{" << i << "}";
  }

  if (use_label) {
    ss << "\\draw (labelpoint) node{" << graph_label << "};" << endl;
  }
  ss << "\\end{tikzpicture} " << endl;

  return ss.str();
}
template string flag::print_latex(bool use_label, const int &graph_label,
                                  bool color_1_nonedge);

string flag::print(const string &delimeter) const {
  stringstream ss;

  // Prints the number of vertices, zero as not type and then part of the
  // diagonal matrix
  ss << m_vertices << delimeter;
  if (m_Theta_class != 0)
    ss << m_Theta << "." << m_Theta_class << delimeter << delimeter;
  else
    ss << m_Theta << delimeter << delimeter;

#ifdef G_COLORED_VERTICES
  print_vertex_coloring(ss, delimeter);
  ss << delimeter;
#endif

#ifdef G_COLORED_EDGES
  for (int u = 0; u < m_vertices; u++) {
    for (int v = u + 1; v < m_vertices; v++)
      ss << delimeter << m_color_edge[u][v];
    ss << delimeter;
  }
#endif

  return ss.str();
}

#ifdef G_COLORED_VERTICES
void flag::print_vertex_coloring(ostream &ostr, const string &delimiter) const {
  for (int i = 0; i < m_vertices; i++) {
    ostr << delimiter << m_color_vertex[i];
  }
}
#endif

#ifdef G_COLORED_EDGES
void flag::print_edge_coloring(ostream &ostr, const string &delimiter) const {
  for (int u = 0; u < m_vertices; u++) {
    for (int v = u + 1; v < m_vertices; v++)
      ostr << delimiter << m_color_edge[u][v];
    ostr << delimiter;
  }
}
#endif

/*****************************************************************************************************************
 * Loading ****************/
void flag::load_from_string(const char *str) {
  if (strcmp(str, "cin") == 0) {
    load_from_stream(cin, -1, -1);
  } else {
    stringstream s(str);
    load_from_stream(s, -1, -1);
  }
}

bool flag::load_from_stream(istream &stream, int assumed_vertices,
                            int assumed_theta) {

  // int poistion = stream.tellg();
  // cerr << "load_from_stream  Rading: " << stream.tellg() <<  "; " <<
  // stream.rdbuf() << endl; stream.seekg(poistion);

  int vertices, theta, theta_class = 0; //, theta;
  string theta_dot_class; //  (theta loaded as string since it may be
                          //  THETA.CLASS)
  stream >> vertices;

  if (stream.fail())
    return false;

  // stream >> theta_dot_class;
  stream >> theta;

  if (stream.fail()) {
    return false;
  }

  int nextcharacter = stream.peek();
  if (nextcharacter == '.') {
#ifndef G_FLAG_PRODUCTS_SLOW_LOW_MEMORY
    cerr << "The program was not compiled with support for repeated flag "
            "types.\n";
    cerr << "Use -DG_FLAG_PRODUCTS_SLOW_LOW_MEMORY when compiling the program."
         << endl;
    assert(0);
#endif
    char ch;
    stream >> ch;
    stream >> theta_class;
  }
  if (nextcharacter == EOF) {
    stream.clear();
  }

  if (stream.fail()) {
    cerr << "Fail after . read" << endl;
    return false;
  }

  /*
//       ss << m_Theta << "." <<  m_Theta_class << delimeter << delimeter;
  string only_theta="", only_theta_class="";

  std::size_t found_dot = theta_dot_class.find(".");
  if (found_dot != std::string::npos)
  {
      only_theta = theta_dot_class.substr(0,found_dot);
      only_theta_class = theta_dot_class.substr(found_dot+1);
#ifndef G_FLAG_PRODUCTS_SLOW_LOW_MEMORY
      cerr << "The program was not compiled with support for repeated flag
types.\n"; cerr << "Use -DG_FLAG_PRODUCTS_SLOW_LOW_MEMORY when compiling the
program." << endl; assert(0); #endif
  }
  else
  {
      only_theta = theta_dot_class;
      only_theta_class = "";
  }
  theta = stol(only_theta);
  if (only_theta_class.length() > 0)
  {
      theta_class = stol(only_theta_class);
  }
  */

  if (assumed_theta != -1 && theta != assumed_theta) {
    cerr << "Err: assumed theta " << assumed_theta << " but loaded " << theta
         << " for " << vertices << " vertices " << endl;
    cerr << "Rest of the stream: " << stream.rdbuf() << endl;
  }
  assert(assumed_theta == -1 || theta == assumed_theta);
  assert(assumed_vertices == -1 || assumed_vertices == vertices);
  set_vertices_and_Theta(vertices, theta, theta_class);

  // cerr << "Set " <<vertices << " " << theta << " " << theta_class << endl;

#ifdef G_COLORED_VERTICES
  for (int i = 0; i < m_vertices; i++) {
    int color_v;
    stream >> color_v;
    color_vertex(i, color_v);
  }
#endif

#ifdef G_COLORED_EDGES
  int color;
  for (int u = 0; u < vertices; u++) {
    for (int v = u + 1; v < vertices; v++) {
      stream >> color;
      color_edge(u, v, color);
    }
  }
#endif

  if (!stream) {
    cerr << "Stream failed while loading flag " << print() << endl;
    assert(0);
  }

  // std::cout << " good()=" << stream.good() << endl;
  // std::cout << " eof()=" << stream.eof() << endl;
  // std::cout << " fail()=" << stream.fail() << endl;
  // std::cout << " bad()=" << stream.bad() << endl;
  // cerr << print() << endl;

  create_minlex_signature();

  return true;
}

flag::flag(const flag &f) { *this = f; }

flag &flag::operator=(const flag &f) {
  m_vertices = f.m_vertices;
  m_Theta = f.m_Theta;
  m_Theta_class = f.m_Theta_class;

#ifdef G_COLORED_VERTICES
  for (int i = 0; i < m_vertices; i++) {
    m_color_vertex[i] = f.m_color_vertex[i];
  }

  // memcpy(m_color_vertex, f.m_color_vertex, sizeof(int)*m_vertices);

  for (int i = 0; i < COLORS_VERTICES; i++) {
    m_colored_vertices[i] = f.m_colored_vertices[i];
  }
#endif

#ifdef G_COLORED_EDGES
  for (int i = 0; i < m_vertices; i++)
    for (int j = 0; j < m_vertices; j++) {
      m_color_edge[i][j] = f.m_color_edge[i][j];
    }

  for (int i = 0; i < COLORS_EDGES; i++) {
    m_colored_edges[i] = f.m_colored_edges[i];
  }

#endif

  return *this;
}
bool flag::operator<(const flag &f) const {
  string me = print();
  string other = f.print();
  return me < other;
}

bool is_flag_forbidden(const flag &g, int verbose_output) {
  // cerr << "Testing forbidden " << g.print() << endl;
  // verbose_output = 1;

#ifdef HACK_ALL_COLOR_FORBIDDEN
  if (g.m_vertices < 4) {
    return false;
  }

  // For each color we count how many times it is used on G
  // Also check that the graph is a complete graph
  int color_counts[HACK_ALL_COLOR_FORBIDDEN_COLORS];
  for (int c = 0; c < HACK_ALL_COLOR_FORBIDDEN_COLORS; c++) {
    color_counts[c] = 0;
  }

  for (int u = 0; u < g.m_vertices; u++)
    for (int v = u + 1; v < g.m_vertices; v++) {
      int label = (g.m_color_edge[u][v]) - 1;
      if (label == 0) {
        return false;
      }
      for (int c = 0; c < HACK_ALL_COLOR_FORBIDDEN_COLORS; c++) {
        color_counts[c] += (label >> c) & 1;
      }
    }

  int colors_used = 0;
  for (int c = 0; c < HACK_ALL_COLOR_FORBIDDEN_COLORS; c++) {
    if (color_counts[c] > 0)
      colors_used++;
  }

  if (colors_used != HACK_ALL_COLOR_FORBIDDEN_COLORS) {
    return false;
  }

  // Now we need to check if we can find an actual copy that has all colors. Not
  // just a fake for example by a random edge. We pick
  // HACK_ALL_COLOR_FORBIDDEN_COLORS edges and see if they are in the correct
  // color
  if (g.m_vertices >= 4) {
    /*
    int edge_vertices[][2] = { [0,1], [0,2], [0,3], [1,2], [1,3], [2,3] };
    for (e1 = 0; e1 < 6; e1++)
    {
        if (  ((g.m_color_edge[c1_u][c1_v]-1)>>0)& 1 != 1) continue;

    }
    */

    int colors_picked[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11};
    int colors_valid = HACK_ALL_COLOR_FORBIDDEN_COLORS;
    /// int colors_picked[]={0,0,1,1,2,2};
    // int colors_picked[]={0,0,0,1,1,2};
    // int colors_picked[]={0,0,0,0,1,2};
    // int colors_valid = 6;
    assert(HACK_ALL_COLOR_FORBIDDEN_COLORS < 11);

    int edges_used[] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int edges_u[V * V];
    int edges_v[V * V];
    int e = 0;
    for (int u = 0; u < g.m_vertices; u++)
      for (int v = u + 1; v < g.m_vertices; v++) {
        edges_u[e] = u;
        edges_v[e] = v;
        e++;
      }

    if (forbidden_rainbow_pattern_found(g, colors_picked, colors_valid,
                                        edges_used, edges_u, edges_v, 0)) {
      return true;
    }

    return false;
    /*

    for (int c1_u = 0; c1_u < g.m_vertices; c1_u++)
        for (int c1_v = c1_u+1; c1_v < g.m_vertices; c1_v++)
        {
            if (  (((g.m_color_edge[c1_u][c1_v]-1)>>0)& 1) != 1) continue;

            for (int c2_u = 0; c2_u < g.m_vertices; c2_u++)
                for (int c2_v = c2_u+1; c2_v < g.m_vertices; c2_v++)
                {
                    if (c1_u == c2_u && c1_v == c2_v) continue;
                    if (  (((g.m_color_edge[c2_u][c2_v]-1)>>1)& 1) != 1)
    continue;

                    for (int c3_u = 0; c3_u < g.m_vertices; c3_u++)
                        for (int c3_v = c3_u+1; c3_v < g.m_vertices; c3_v++)
                        {
                            if (c1_u == c3_u && c1_v == c3_v) continue;
                            if (c2_u == c3_u && c2_v == c3_v) continue;
                            if (  (((g.m_color_edge[c3_u][c3_v]-1)>>2)& 1) != 1)
    continue; for (int c4_u = 0; c4_u < g.m_vertices; c4_u++) for (int c4_v =
    c4_u+1; c4_v < g.m_vertices; c4_v++)
                                {
                                    if (c1_u == c4_u && c1_v == c4_v) continue;
                                    if (c2_u == c4_u && c2_v == c4_v) continue;
                                    if (c3_u == c4_u && c3_v == c4_v) continue;
                                    if (  (((g.m_color_edge[c4_u][c4_v]-1)>>3)&
    1) != 1) continue;

                                    if (HACK_ALL_COLOR_FORBIDDEN_COLORS == 4)
                                    {
                                        return true;
                                    }
                                    else
                                    {
                                        for (int c5_u = 0; c5_u < g.m_vertices;
    c5_u++) for (int c5_v = c5_u+1; c5_v < g.m_vertices; c5_v++)
                                            {
                                                if (c1_u == c5_u && c1_v ==
    c5_v) continue; if (c2_u == c5_u && c2_v == c5_v) continue; if (c3_u == c5_u
    && c3_v == c5_v) continue; if (c4_u == c5_u && c4_v == c5_v) continue; if (
    (((g.m_color_edge[c5_u][c5_v]-1)>>4)& 1) != 1) continue;

                                                if
    (HACK_ALL_COLOR_FORBIDDEN_COLORS == 5)
                                                {
                                                    return true;
                                                }
                                                else
                                                {
                                                    for (int c6_u = 0; c6_u <
    g.m_vertices; c6_u++) for (int c6_v = c6_u+1; c6_v < g.m_vertices; c6_v++)
                                                        {
                                                            if (c1_u == c6_u &&
    c1_v == c6_v) continue; if (c2_u == c6_u && c2_v == c6_v) continue; if (c3_u
    == c6_u && c3_v == c6_v) continue; if (c4_u == c6_u && c4_v == c6_v)
    continue; if (c5_u == c6_u && c5_v == c6_v) continue; if (
    (((g.m_color_edge[c6_u][c6_v]-1)>>5)& 1) != 1) continue;

                                                            return true;
                                                        }
                                                }
                                            }
                                    }
                                }
                        }
                }
        }
    */
  }

  return false;
#endif

  for (int i = 0; i < (int)get_forbidden_subflags().size(); i++) {
    // cerr << "Testing " << g.print() << " and " <<
    // g_forbidden_subflags[i].print() << endl;
    if (g.contains_as_subflag(get_forbidden_subflags()[i])) {
      if (verbose_output > 0) {
        cerr << "Forbidding " << g.print() << " because of "
             << get_forbidden_subflags()[i].print() << endl;
      }
      return true;
    }
  }
  return false;
}

bool flag_and_coefficient::operator==(const flag_and_coefficient &fc) const {
  if (fc.coefficient != coefficient)
    return false;
  return g.is_isomorphic_to(fc.g);
}

void flag_and_coefficient::print(std::ostream &stream, bool skip_if_empty,
                                 string prefix) const {
  double rounded = smart_round(coefficient);
  if (skip_if_empty && rounded == 0)
    return;
  if (prefix != "") {
    stream << left << setw(7) << prefix;
  }
  stream << fixed << setprecision(G_PRECISION) << smart_round(coefficient)
         << "  " << g.print() << endl;
};

std::ostream &operator<<(std::ostream &stream, const flag_and_coefficient &fc) {
  fc.print(stream);

  return stream;
}

int find_flag_in_list_nonparalel(const flag &f, const vector<flag> &fv) {
  for (int i = 0; i < (int)fv.size(); ++i) {
    if (f.is_isomorphic_to(fv[i])) {
      return i;
    }
  }
  return -1;
}

int find_flag_in_list(const flag &f, const vector<flag> &fv) {
  // Not worth paralelizing the search.
  if (fv.size() < 500) {
    return find_flag_in_list_nonparalel(f, fv);
  }

  // attempt to paralelize the search if it was not parallel yet
#ifdef _USING_OMP_
  if (omp_in_parallel() == false) {
    // cerr << "Trying in parallel" << endl;
    volatile int location = -1;

#pragma omp parallel for shared(location)
    for (int i = 0; i < (int)fv.size(); ++i) {
      if (location != -1)
        continue;
      if (f.is_isomorphic_to(fv[i])) {
        location = i;
      }
    }
    return location;
  }
#endif

  return find_flag_in_list_nonparalel(f, fv);
}

int get_flag_type_in_list(const flag &f, vector<vector<flag>> &flags) {
#ifdef USE_REDUCED_TYPES
  const int identity[] = {0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
                          11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,
                          22, 23, 24, 25, 26, 27, 28, 29, 30};
  flag f_type;
  f_type.as_subflag(f, identity, f.m_Theta, 0);
#endif

  // Find correct type
  for (int i = 0; i < (int)flags.size(); i++) {
    if (flags[i].size() == 0) {
      cerr << "Some bug" << endl;
    }
    if (f.have_same_type(flags[i][0]))
      return i;
#ifdef USE_REDUCED_TYPES
    if (f.m_Theta == flags[i][0].m_Theta &&
        f.m_Theta_class == flags[i][0].m_Theta_class) {
      flag g_flags_type;
      g_flags_type.as_subflag(flags[i][0], identity, flags[i][0].m_Theta, 0);
      if (g_flags_type.is_isomorphic_to(f_type))
        return -2;
    }
#endif
  }

  return -1;
}

void include_flag_in_list(const flag &f, vector<vector<flag>> &flags, int id,
                          bool ignore_reduced_types) {
  if (id == -3) {
    id = get_flag_type_in_list(f, flags);
  }
#ifdef USE_REDUCED_TYPES
  if (id == -2 && ignore_reduced_types)
    return;
#endif
  if (id == -1 || id == -2) {
    vector<flag> new_type_list;
    new_type_list.push_back(f);
    flags.push_back(new_type_list);
    return;
  }

  flags[id].push_back(f);
}

void include_flag_in_list_if_new(const flag &f, vector<vector<flag>> &flags,
                                 bool ignore_reduced_types) {
  int id = get_flag_type_in_list(f, flags);
#ifdef USE_REDUCED_TYPES
  if (id == -2 && ignore_reduced_types) {
    return;
  }
#endif
  if (id == -1 || id == -2) {
    //        cerr << "new type" << endl;
    include_flag_in_list(f, flags, id, ignore_reduced_types);
    return;
  }

  if (find_flag_in_list(f, flags[id]) != -1) {
    //        cerr << "duplicate" << endl;
    return;
  }
  //    cerr << "new for known type" << f.print() << endl;
  flags[id].push_back(f);
}
pair<int, int> P_F1_IN_H(const flag &F1, const flag &H, bool density) {

  assert(F1.m_vertices <= H.m_vertices);

  flag type;
  F1.get_type_subflag(type);

  if (H.labeled_vertices_cnt() > 0) {
    flag typeH;
    H.get_type_subflag(typeH);
    if (!type.contains_as_subflag(typeH))
      return make_pair(0, 1);
  }

  // Some easy checks...
#if defined(G_COLORED_EDGES) && !defined(G_COLORED_EDGES_BLIND)

  for (int c = 1; c < COLORS_EDGES; c++) {
    if (F1.m_colored_edges[c] > H.m_colored_edges[c])
      return make_pair(0, 1);
  }
#endif

#if defined(G_COLORED_VERTICES) && !defined(G_COLORED_VERTICES_BLIND)
  for (int c = 1; c < COLORS_VERTICES; c++) {
    if (F1.m_colored_vertices[c] > H.m_colored_vertices[c])
      return make_pair(0, 1);
  }
#endif

  int good_maps = 0;
  int all_maps = 0;

  // Map labeled vertices as 0,1,2,theta,-1,-1,-1,-4,-4,-4
  // where -1 indicates number of vertices in F1
  int vertex_split[2 * V]; // the 2*V instead of V is here for overflow warning

  int t = 0;
  for (int i = 0; i < type.m_vertices; i++)
    vertex_split[t++] = i;
  for (int i = type.m_vertices; i < F1.m_vertices; i++)
    vertex_split[t++] = -1;
  for (; t < H.m_vertices; t++)
    vertex_split[t] = -4;

  std::sort(vertex_split + H.labeled_vertices_cnt(),
            vertex_split + H.m_vertices);

  // cerr << "FFF" << endl;

  flag tmp;
  flag tmp_type;
  int mapping[V];
  int next_to_map; // temporary
  do {

    // cerr << "DDDD" << endl;
    // for (int i  = 0; i < F1.m_vertices; i++)
    //{
    //     cerr << i << "->" << mapping[i] << " ";
    // }
    // cerr << endl;

    // check if Theta correct
    for (int j = 0; j < H.m_vertices; j++) {
      if (vertex_split[j] < 0)
        continue;
      if (vertex_split[j] < type.m_vertices)
        mapping[vertex_split[j]] = j;
    }

    // Finish creating the mapping
    next_to_map = type.m_vertices;
    for (int j = 0; j < H.m_vertices; j++)
      if (vertex_split[j] == -1)
        mapping[next_to_map++] = j;

    tmp_type.as_subflag(H, mapping, type.m_vertices, type.m_Theta);
    tmp.as_subflag(H, mapping, F1.m_vertices, type.m_Theta);

#ifdef G_COLORED_VERTICES_SAMPLED_SEPARATELY_BY_COLORS
    // check if the labeling is color preserving
    bool color_preserved = true;
    for (int v = 0; v < F1.m_Theta; v++) {
      if (F1.m_color_vertex[v] != tmp_type.m_color_vertex[v]) {
        color_preserved = false;
        // cerr << "Color not preserved for labeled  v=" << v << endl;
        break;
      }
    }
    if (color_preserved == false) {
      continue;
    }

    // check if the map is color preserving
    color_preserved = true;
    for (int c = 1; c < COLORS_VERTICES; c++) {
      if (F1.m_colored_vertices[c] != tmp.m_colored_vertices[c]) {
        color_preserved = false;
        break;
      }
    }
    if (color_preserved == false)
      continue;
#endif

    all_maps++;

    // check if type is correct
    if (!tmp_type.is_isomorphic_to(type))
      continue;

    // check if F1 is correct
    tmp.as_subflag(H, mapping, F1.m_vertices, type.m_Theta);
    if (!tmp.is_isomorphic_to(F1))
      continue;

    good_maps++;

  } while (std::next_permutation(vertex_split + H.labeled_vertices_cnt(),
                                 vertex_split + H.m_vertices));

  // cerr << "Good/All: " << good_maps << "/" << all_maps << endl;
  if (all_maps == 0) {
    cerr << "all_maps is zero" << endl;
    cerr << "Calculating P(F1,H), where" << endl;
    cerr << "F1 = " << F1.print() << endl;
    cerr << " H = " << H.print() << endl;
    assert(0);
  }
  return make_pair(good_maps, all_maps);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// P(F1,F2,H)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// P(F1,F2,H)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// P(F1,F2,H)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// P(F1,F2,H)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// P(F1,F2,H)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// P(F1,F2,H)

// Rest of F2 is also piskced as a subset
void pick_F2_mapping(const flag &F1, const flag &F2, const flag &H,
                     int *mapping_F2, int next_to_map, bool *used_from_H,
                     int &good_maps) {
  if (next_to_map >= F2.m_vertices) {
    flag H_F2;
    H_F2.as_subflag(H, mapping_F2, F2.m_vertices, F2.m_Theta);
    if (!H_F2.is_isomorphic_to(F2))
      return;

    good_maps++;

    return;
  }

  int min_vertex = 0;
  if (next_to_map > F2.m_Theta)
    min_vertex = mapping_F2[next_to_map - 1] + 1;

  for (int v = min_vertex; v < H.m_vertices; v++) {
    if (used_from_H[v])
      continue;

    used_from_H[v] = true;
    mapping_F2[next_to_map] = v;
    pick_F2_mapping(F1, F2, H, mapping_F2, next_to_map + 1, used_from_H,
                    good_maps);
    used_from_H[v] = false;
  }
}

// Rest of F1 is picked only as a subset
void pick_F1_mapping(const flag &F1, const flag &F2, const flag &H,
                     int *mapping_F1, int *mapping_F2, int next_to_map,
                     bool *used_from_H, int &good_maps) {
  if (next_to_map >= F1.m_vertices) {
    // cout << "." << endl;

    // Check if the mapping is corect!
    flag H_F1;
    H_F1.as_subflag(H, mapping_F1, F1.m_vertices, F1.m_Theta);

    if (!H_F1.is_isomorphic_to(F1))
      return;

    pick_F2_mapping(F1, F2, H, mapping_F2, F2.m_Theta, used_from_H, good_maps);
    return;
  }

  int min_vertex = 0;
  if (next_to_map > F1.m_Theta)
    min_vertex = mapping_F1[next_to_map - 1] + 1;

  for (int v = min_vertex; v < H.m_vertices; v++) {
    if (used_from_H[v])
      continue;

    used_from_H[v] = true;
    mapping_F1[next_to_map] = v;
    pick_F1_mapping(F1, F2, H, mapping_F1, mapping_F2, next_to_map + 1,
                    used_from_H, good_maps);
    used_from_H[v] = false;
  }
}

// Theta mapping is picked as all permutations...
void pick_theta_mapping(const flag &F1, const flag &F2, const flag &H,
                        int *mapping_theta, int next_to_map, bool *used_from_H,
                        int &good_maps) {
  if (next_to_map >= F1.m_Theta) {
    int mapping_F1[F1.m_vertices];
    int mapping_F2[F2.m_vertices];
    for (int i = 0; i < F1.m_Theta; i++) {
      mapping_F2[i] = mapping_F1[i] = mapping_theta[i];
    }

    pick_F1_mapping(F1, F2, H, mapping_F1, mapping_F2, F1.m_Theta, used_from_H,
                    good_maps);
    return;
  }

  for (int v = 0; v < H.m_vertices; v++) {
    if (used_from_H[v])
      continue;

    if (!H.is_map_up_to_v_correct(next_to_map, v, mapping_theta, F1))
      continue;

    used_from_H[v] = true;
    mapping_theta[next_to_map] = v;
    pick_theta_mapping(F1, F2, H, mapping_theta, next_to_map + 1, used_from_H,
                       good_maps);
    used_from_H[v] = false;
  }
}

double P_F1_F2_IN_H(const flag &F1, const flag &F2, const flag &H,
                    bool return_only_good_maps = false) {
  assert(F1.m_Theta == F2.m_Theta);

#if defined(G_COLORED_EDGES) && !defined(G_COLORED_EDGES_BLIND)
  for (int c = 1; c < COLORS_EDGES; c++) {
    //		cout << "e" << endl;
    if (F1.m_colored_edges[c] > H.m_colored_edges[c])
      return 0;
    if (F2.m_colored_edges[c] > H.m_colored_edges[c])
      return 0;
  }
#endif

#ifdef G_COLORED_VERTICES
  for (int c = 1; c < COLORS_VERTICES; c++) {
    // cout << "v " << F1.m_colored_vertices[c] << " " <<
    // H.m_colored_vertices[c] << endl;
    if (F1.m_colored_vertices[c] > H.m_colored_vertices[c])
      return 0;
    if (F2.m_colored_vertices[c] > H.m_colored_vertices[c])
      return 0;
  }
#endif

  int mapping_theta[F1.m_Theta];
  bool used_from_H[H.m_vertices];

  for (int i = 0; i < H.m_vertices; i++) {
    used_from_H[i] = false;
  }

  int good_maps = 0;
  pick_theta_mapping(F1, F2, H, mapping_theta, 0, used_from_H, good_maps);

  if (return_only_good_maps)
    return good_maps;

  int all_maps = 1;
#ifdef G_COLORED_VERTICES_SAMPLED_SEPARATELY_BY_COLORS
  int colors_available[COLORS_VERTICES];
  int colors_needed_F1[COLORS_VERTICES];
  int colors_needed_F2[COLORS_VERTICES];
  colors_available[0] = 0;
  for (int c = 1; c < COLORS_VERTICES; c++) {
    colors_available[c] = H.m_colored_vertices[c];
    colors_needed_F1[c] = F1.m_colored_vertices[c];
    colors_needed_F2[c] = F2.m_colored_vertices[c];
  }
  for (int i = 0; i < F1.m_Theta; i++) {
    all_maps *= colors_available[F1.m_color_vertex[i]]--;
    colors_needed_F1[F1.m_color_vertex[i]]--;
    colors_needed_F2[F1.m_color_vertex[i]]--;
  }
  for (int c = 1; c < COLORS_VERTICES; c++) {
    all_maps *= binomial(colors_available[c], colors_needed_F1[c]);
    all_maps *= binomial(colors_available[c] - colors_needed_F1[c],
                         colors_needed_F2[c]);
  }
#else
  for (int i = 0; i < F1.m_Theta; i++)
    all_maps *= (H.m_vertices - i);
  all_maps *= binomial(H.m_vertices - F1.m_Theta, F1.m_vertices - F1.m_Theta);
  all_maps *=
      binomial(H.m_vertices - F1.m_vertices, F2.m_vertices - F1.m_Theta);
#endif
  //    if (good_maps > 0)
  //        cout << "good_maps/all_maps: " << good_maps << "/" << all_maps <<
  //        endl;

  return (double)good_maps / (double)all_maps;
}

pair<int, int> P_F1_F2_IN_labeled_H(const flag &F1, const flag &F2,
                                    const flag &H) {

  // cerr << F1.m_Theta << " " << F2.m_Theta << endl;
  assert(F1.m_Theta == F2.m_Theta);
  if (F1.m_vertices + F2.m_vertices - F1.labeled_vertices_cnt() >
      H.m_vertices) {
    cerr << F1.m_vertices << "+" << F2.m_vertices << "-"
         << F1.labeled_vertices_cnt() << " > " << H.m_vertices << endl;
  }
  assert(F1.m_vertices + F2.m_vertices - F1.labeled_vertices_cnt() <=
         H.m_vertices);
  assert(F1.labeled_vertices_cnt() >= H.labeled_vertices_cnt());

  flag type;
  F1.get_type_subflag(type);

  flag type_F2;
  F2.get_type_subflag(type_F2);

  if (type.is_isomorphic_to(type_F2) == false) {
    return make_pair(0, 1);
  }

  if (H.labeled_vertices_cnt() > 0) {
    flag typeH;
    H.get_type_subflag(typeH);
    if (!type.contains_as_subflag(typeH))
      return make_pair(0, 1);
  }

  // Some easy checks...
#if defined(G_COLORED_EDGES) && !defined(G_COLORED_EDGES_BLIND)
  for (int c = 1; c < COLORS_EDGES; c++) {
    if (F1.m_colored_edges[c] + F2.m_colored_edges[c] -
            type.m_colored_edges[c] >
        H.m_colored_edges[c])
      return make_pair(0, 1);
  }
#endif
#if defined(G_COLORED_VERTICES) && !defined(G_COLORED_VERTICES_BLIND)
  for (int c = 1; c < COLORS_VERTICES; c++) {
    if (F1.m_colored_vertices[c] + F2.m_colored_vertices[c] -
            type.m_colored_vertices[c] >
        H.m_colored_vertices[c])
      return make_pair(0, 1);
  }
#endif

  int good_maps = 0;
  int all_maps = 0;

  int vertex_split[V];

  int t = 0;
  for (int i = 0; i < type.m_vertices; i++)
    vertex_split[t++] = i;
  for (int i = type.m_vertices; i < F1.m_vertices; i++)
    vertex_split[t++] = -1;
  for (int i = type.m_vertices; i < F2.m_vertices; i++)
    vertex_split[t++] = -2;
  // for (int i = type.m_vertices; i < F3.m_vertices; i++) vertex_split[t++] =
  // -3;
  for (; t < H.m_vertices; t++)
    vertex_split[t] = -4;

  std::sort(vertex_split + H.labeled_vertices_cnt(),
            vertex_split + H.m_vertices);

  flag tmp_type;
  flag tmp_F1;
  flag tmp_F2;
  int mapping[V];
  int next_to_map; // temporary
  do {

    // check if Theta correct
    for (int j = 0; j < H.m_vertices; j++) {
      if (vertex_split[j] < 0)
        continue;
      if (vertex_split[j] < type.m_vertices)
        mapping[vertex_split[j]] = j;
    }
    tmp_type.as_subflag(H, mapping, type.m_vertices, type.m_Theta);

    // cerr << tmp.print() << " vs " << type.print();

    // check if F1 correct
    next_to_map = type.m_vertices;
    for (int j = 0; j < H.m_vertices; j++)
      if (vertex_split[j] == -1)
        mapping[next_to_map++] = j;
    tmp_F1.as_subflag(H, mapping, F1.m_vertices, type.m_Theta);

    // check if F2 correct
    next_to_map = type.m_vertices;
    for (int j = 0; j < H.m_vertices; j++)
      if (vertex_split[j] == -2)
        mapping[next_to_map++] = j;
    tmp_F2.as_subflag(H, mapping, F2.m_vertices, type.m_Theta);

#ifdef G_COLORED_VERTICES_SAMPLED_SEPARATELY_BY_COLORS
    // check if the labeling is color preserving
    bool color_preserved = true;
    for (int v = 0; v < F1.m_Theta; v++) {
      if (F1.m_color_vertex[v] != tmp_type.m_color_vertex[v]) {
        color_preserved = false;
        break;
      }
    }
    if (color_preserved == false)
      continue;

    // check if the map is color preserving
    color_preserved = true;
    for (int c = 1; c < COLORS_VERTICES; c++) {
      if (F1.m_colored_vertices[c] != tmp_F1.m_colored_vertices[c] ||
          F2.m_colored_vertices[c] != tmp_F2.m_colored_vertices[c]) {
        color_preserved = false;
        break;
      }
    }
    if (color_preserved == false)
      continue;
#endif

    all_maps++;

    if (!tmp_type.is_isomorphic_to(type))
      continue;
    if (!tmp_F1.is_isomorphic_to(F1))
      continue;
    if (!tmp_F2.is_isomorphic_to(F2))
      continue;

    good_maps++;

  } while (std::next_permutation(vertex_split + H.labeled_vertices_cnt(),
                                 vertex_split + H.m_vertices));

  // cerr << "Good/All: " << good_maps << "/" << all_maps << endl;

  return make_pair(good_maps, all_maps);
}
