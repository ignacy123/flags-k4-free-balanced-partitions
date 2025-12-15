#pragma once

#include "flag-calculator.hpp"

int solve_sdp_for_problem(const vector<FlagVector<0, V>> &constraints,
                          const FlagVector<0, V> &objective);
