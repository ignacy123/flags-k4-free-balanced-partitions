# This is a test of how well a cut is going to perform on a sparsified balanced tripartite graph.
# The cut itself is based on a triangle uvw. We make a lot of optimistic assumptions.
# We initialize the blobs by common neighborhoods of edges.
# For every vertex, we distribute vertices only seen by it to the neighboring blobs. The rest is distributed randomly.
from scipy.optimize import minimize, Bounds, differential_evolution
import numpy as np
from itertools import product


bounds = Bounds([1 / 2], [1])


def edges_inside(args: list[float]):
    p = args[0]

    correct = 1 / 3 * p * p
    correct += 1 / 3 * p * (1 - p)
    incorrect_one = 1 / 2 * 1 / 3 * p * (1 - p)

    space = 1 / 3 - correct - 2 * incorrect_one

    correct += 1 / 3 * space
    incorrect_one += 1 / 3 * space

    total = correct + 2 * incorrect_one
    raw_edges = correct * 2 * incorrect_one + incorrect_one * incorrect_one

    print(f"{total=}, {raw_edges=}")

    return -p * 3 * raw_edges


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
