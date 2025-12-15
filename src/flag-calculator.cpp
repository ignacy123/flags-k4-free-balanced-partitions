#include "flag-calculator.hpp"
#include "flag.hpp"
#include <cassert>
#include <utility>

pair<bool, flag> get_edge_between(flag flag_left, flag flag_right) {
  assert(flag_left.have_same_type(flag_right));
  assert(flag_left.m_vertices == flag_left.m_Theta + 1);
  assert(flag_left.m_vertices == flag_right.m_vertices);

  flag result;
  result.set_vertices_and_Theta(flag_left.m_Theta + 2, flag_left.m_Theta);
  // Copy type.
  for (int i = 0; i < flag_left.m_Theta; i++) {
    for (int j = i + 1; j < flag_left.m_Theta; j++) {
      result.color_edge(i, j, flag_left.m_color_edge[i][j]);
    }
  }
  for (int i = 0; i < flag_left.m_Theta; i++) {
    result.color_edge(i, flag_left.m_Theta,
                      flag_left.m_color_edge[i][flag_left.m_Theta]);
  }
  for (int i = 0; i < flag_right.m_Theta; i++) {
    result.color_edge(i, flag_right.m_Theta + 1,
                      flag_right.m_color_edge[i][flag_right.m_Theta]);
  }
  result.color_edge(flag_left.m_Theta, flag_left.m_Theta + 1, 2);
  return pair<bool, flag>(!is_flag_forbidden(result), result);
}
