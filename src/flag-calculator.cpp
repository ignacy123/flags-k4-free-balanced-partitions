#include "flag-calculator.hpp"
#include "flag.hpp"
#include <cassert>
#include <utility>

pair<bool, flag_coeff> get_edge_between(const flag_coeff &flag_left,
                                        const flag_coeff &flag_right) {
  assert(flag_left.g.have_same_type(flag_right.g));
  assert(flag_left.g.m_vertices == flag_left.g.m_Theta + 1);
  assert(flag_left.g.m_vertices == flag_right.g.m_vertices);

  flag raw_flag;
  raw_flag.set_vertices_and_Theta(flag_left.g.m_Theta + 2, flag_left.g.m_Theta);
  // Copy type.
  for (int i = 0; i < flag_left.g.m_Theta; i++) {
    for (int j = i + 1; j < flag_left.g.m_Theta; j++) {
      raw_flag.color_edge(i, j, flag_left.g.m_color_edge[i][j]);
    }
  }

#ifdef G_COLORED_VERTICES
  for (int i = 0; i < flag_left.g.m_Theta; i++) {
    raw_flag.color_vertex(i, flag_left.g.m_color_vertex[i]);
  }
  raw_flag.color_vertex(flag_left.g.m_Theta,
                        flag_left.g.m_color_vertex[flag_left.g.m_Theta]);
  raw_flag.color_vertex(flag_right.g.m_Theta + 1,
                        flag_right.g.m_color_vertex[flag_right.g.m_Theta]);
#endif

  for (int i = 0; i < flag_left.g.m_Theta; i++) {
    raw_flag.color_edge(i, flag_left.g.m_Theta,
                        flag_left.g.m_color_edge[i][flag_left.g.m_Theta]);
  }
  for (int i = 0; i < flag_right.g.m_Theta; i++) {
    raw_flag.color_edge(i, flag_right.g.m_Theta + 1,
                        flag_right.g.m_color_edge[i][flag_right.g.m_Theta]);
  }
  raw_flag.color_edge(flag_left.g.m_Theta, flag_left.g.m_Theta + 1, 2);
  flag_coeff result;
  result.g = raw_flag;
  result.coefficient = flag_left.coefficient * flag_right.coefficient;
  return pair<bool, flag_coeff>(!is_flag_forbidden(raw_flag), result);
}
