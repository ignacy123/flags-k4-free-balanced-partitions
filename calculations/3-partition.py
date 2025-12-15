from scipy.optimize import minimize, Bounds, differential_evolution
import numpy as np
from itertools import product


bounds = Bounds([1 / 2], [1])


def edges_inside(p: float):
    c_1 = 1 / 3 * p * p
    c_2_c_3 = 1 / 3 * p * (1 - p)
    r_1_r_2 = 1 / 3 * p
    total_left = c_2_c_3 + r_1_r_2
    space_first = 1 / 3 - c_1
    c_2 = space_first / total_left * c_2_c_3
    c_3 = c_2_c_3 - c_2
    r_1 = space_first / total_left * r_1_r_2
    r_2 = r_1_r_2 - r_1
    fixed_second = r_2 + c_3
    space_second = 1 / 3 - fixed_second
    total_remaining = 2 / 3 - fixed_second
    r_3_r_4 = 1 / 3 - r_1_r_2
    c_4_c_5 = 1 / 3 - c_1 - c_2_c_3
    l_1_l_2 = 1 / 3
    r_3 = space_second / total_remaining * r_3_r_4
    r_4 = r_3_r_4 - r_3
    c_4 = space_second / total_remaining * c_4_c_5
    c_5 = c_4_c_5 - c_4
    l_1 = space_second / total_remaining * l_1_l_2
    l_2 = l_1_l_2 - l_1

    raw_edges_first = c_1 * r_1 + c_2 * r_1
    raw_edges_second = (
        r_2 * c_3
        + r_2 * l_1
        + r_2 * c_4
        + c_3 * l_1
        + c_3 * r_3
        + l_1 * r_3
        + l_1 * c_4
        + r_3 * c_4
    )
    raw_edges_third = l_2 * r_4 + l_2 * c_5 + r_4 * c_5
    return -p * (raw_edges_first + raw_edges_second + raw_edges_third)


def grid_search(bounds: Bounds):
    precision = 40
    linspaces = []
    largest, largest_witness = 0, None
    for lower, upper in zip(bounds.lb, bounds.ub):
        linspaces.append(np.linspace(lower, upper, precision))

    for idx, starting_point in enumerate(product(*linspaces)):
        if idx % 10000 == 0:
            print(f"Iteration {idx}")
        result = minimize(
            edges_inside,
            starting_point,
            bounds=bounds,
            method="slsqp",
            options={"ftol": 1e-12, "maxiter": 1000, "eps": 1e-12},
        )
        if -result.fun > largest:
            print(f"New maximum: {-result.fun} at {result.x}")
            largest = -result.fun
            largest_witness = result

    return largest_witness


result = grid_search(bounds)
print(f"Minimize: worst partition {-result.fun} at {result.x}")

result = differential_evolution(
    edges_inside,
    bounds=bounds,
)
print(f"Differential evolution: worst partition {-result.fun} at {result.x}")
