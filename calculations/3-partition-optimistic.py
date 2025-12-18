# This is a test of how well a cut is going to perform on a sparsified balanced tripartite graph.
# The cut itself is based on an edge uv. We assume that there is at most 1 / 3 vertices connected to only to u, and at most 1 / 3 vertices connected only to v.
# That is a very optimistic assumption. We initialize 2 blobs by assigning vertices that are only connected to u and vertices that are only connected to v.
# The third blob is common neighborhood. The rest is distributed randomly.
from scipy.optimize import minimize, Bounds, differential_evolution
import numpy as np
from itertools import product


bounds = Bounds([1 / 2], [1])


def edges_inside(args: list[float]):
    p = args[0]
    c_1 = 1 / 3 * p * p
    l_1 = 1 / 3 * p
    c_2 = 1 / 3 * p * (1 - p)

    l_2_l_2_l_3 = 1 / 3 * (1 - p)
    c_3_c_3_c_4 = 1 / 3 - c_1 - 2 * c_2

    total_random = 2 * l_2_l_2_l_3 + c_3_c_3_c_4
    space_left = 1 / 3 - l_1 - c_2
    l_2 = space_left / total_random * l_2_l_2_l_3
    l_3 = l_2_l_2_l_3 - 2 * l_2
    c_3 = space_left / total_random * c_3_c_3_c_4
    c_4 = c_3_c_3_c_4 - 2 * c_3

    total_left = l_1 + 2 * l_2 + l_3
    total_center = c_1 + 2 * c_2 + 2 * c_3 + c_4
    left_side = l_1 + 2 * l_2 + c_2 + c_3
    center_side = c_1 + 2 * l_3 + c_4

    print(f"{total_left=}, {total_center=}, {left_side=}, {center_side=}")

    raw_edges_first = (l_1 + l_2) * (c_2 + c_3 + l_2) + l_2 * (c_2 + c_3)
    raw_edges_third = (c_1 + c_4) * (l_3 + l_3) + l_3 * l_3
    return -p * (raw_edges_first + raw_edges_first + raw_edges_third)


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
