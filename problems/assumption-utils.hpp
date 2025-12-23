#include "flag-calculator.hpp"

template <int theta, int vertices>
FlagVector<0, 2 * vertices> regularity_assumption(flag blueprint) {
  cerr << "Regularity assumption " << blueprint.print() << endl;
  flag type;
  blueprint.get_type_subflag(type);
  type.remove_labeled_vertices();
  FlagVector<theta, vertices> flag_vector(blueprint);
  FlagVector<0, theta> root_vector(type);
  auto unlabeled = flag_vector.template project<0, vertices>();
  auto squared = flag_vector * flag_vector;
  return unlabeled * unlabeled -
         squared.template project<0, 2 * vertices - theta>() * root_vector;
}
