# This is a test of how well a cut is going to perform on a sparsified balanced tripartite graph.
# The cut itself is based on a triangle uvw. We initialize the blobs with the common neighborhoods of the edges of the triangle.
# Then we fill the blob defined by N(vw) by the vertices connected only to v, only to w and to none of the three.
# This is more than 1/3, as this is the non-neighborhood of v (whose degree is at most 2/3).
# We then randomly assign the remaining vertices.
from scipy.optimize import minimize, Bounds, differential_evolution
import numpy as np
from itertools import product


bounds = Bounds([1 / 2], [1])


def edges_inside(args: list[float]):
    p = args[0]
    c_1 = 1 / 3 * p * p
    l_1 = 1 / 3 * p * p
    r_1 = 1 / 3 * p * p
    c_2_c_3_c_4 = 1 / 3 * (1 - p) * (1 - p) + 2 * 1 / 3 * p * (1 - p)
    l_2_l_3_l_4 = 1 / 3 * (1 - p) * p + 1 / 3 * (1 - p) * (1 - p)
    r_2_r_3_r_4 = 1 / 3 * (1 - p) * p + 1 / 3 * (1 - p) * (1 - p)
    l_5_l_6 = 1 / 3 * p * (1 - p)
    r_5_r_6 = 1 / 3 * p * (1 - p)

    space_top = 1 / 3 - 1 / 3 * p * p
    random_possibly_top = c_2_c_3_c_4 + l_2_l_3_l_4 + r_2_r_3_r_4
    c_2 = space_top / random_possibly_top * c_2_c_3_c_4
    l_2 = space_top / random_possibly_top * l_2_l_3_l_4
    r_2 = space_top / random_possibly_top * r_2_r_3_r_4

    c_3_c_4 = c_2_c_3_c_4 - c_2
    l_3_l_4 = l_2_l_3_l_4 - l_2
    r_3_r_4 = r_2_r_3_r_4 - r_2

    c_3 = c_3_c_4 / 2
    l_3 = l_3_l_4 / 2
    r_3 = r_3_r_4 / 2

    c_4 = c_3_c_4 / 2
    l_4 = l_3_l_4 / 2
    r_4 = r_3_r_4 / 2

    l_5 = l_5_l_6 / 2
    r_5 = r_5_r_6 / 2

    l_6 = l_5_l_6 / 2
    r_6 = r_5_r_6 / 2

    total_left = l_1 + l_2 + l_3 + l_4 + l_5 + l_6
    total_right = r_1 + r_2 + r_3 + r_4 + r_5 + r_6
    total_center = c_1 + c_2 + c_3 + c_4

    print(f"{total_left=}, {total_right=}, {total_center=}")

    raw_edges_first = (c_1 + c_2) * (l_2 + r_2) + l_2 * r_2
    raw_edges_second = (l_1 + l_3 + l_5) * (c_3 + r_3 + r_5) + (r_3 + r_5) * c_3
    raw_edges_third = (r_1 + r_5 + r_6) * (c_4 + r_4 + r_6) + (l_4 + l_6) * c_4
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
