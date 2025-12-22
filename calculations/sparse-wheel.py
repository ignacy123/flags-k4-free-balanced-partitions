from scipy.optimize import minimize, Bounds, differential_evolution
import numpy as np
from itertools import product


bounds = Bounds([0, 1 / 2], [1 / 2, 1])


def blue_edge_partition_asymmetric(
    center_size: float, p: float, print_details: bool = False
):
    blob_size = (1 - center_size) / 5
    random = 2 * blob_size * p
    space_left = 1 / 2 - center_size - 2 * p * blob_size
    space_right = 1 / 2 - 4 * blob_size * (1 - p) - blob_size
    prob_left = space_left / random
    prob_right = space_right / random

    one_left = 3 * blob_size * blob_size * p * p * p
    one_left += 2 * blob_size * p * center_size

    one_right = 3 * blob_size * blob_size * p * (1 - p) * p
    one_right += blob_size * blob_size * p * p

    edges_left = 2 * blob_size * p * center_size
    edges_left += one_left * prob_left

    edges_right = 3 * blob_size * blob_size * (1 - p) * (1 - p) * p
    edges_right += 2 * blob_size * blob_size * (1 - p) * p
    edges_right += one_right * prob_right

    if print_details:
        print(f"Blue edge partition, left: {edges_left}, right {edges_right}")

    return max(edges_left, edges_right)


def blue_edge_partition(center_size: float, p: float, print_details: bool = False):
    center_proportion_left = 1
    center_left = center_size * center_proportion_left
    center_right = center_size - center_left

    blob_size = (1 - center_size) / 5
    random = 4 * blob_size * (1 - p) + blob_size
    space_left = 1 / 2 - center_left - 2 * p * blob_size
    space_right = 1 / 2 - center_right - 2 * p * blob_size
    prob_left = space_left / random
    prob_right = space_right / random

    one_left = 3 * blob_size * blob_size * p * (1 - p) * p
    one_left += blob_size * blob_size * p * p
    one_left += 4 * blob_size * (1 - p) * center_left
    one_left += blob_size * center_left

    one_right = 3 * blob_size * blob_size * p * (1 - p) * p
    one_right += blob_size * blob_size * p * p
    one_right += 4 * blob_size * (1 - p) * center_right
    one_right += blob_size * center_right

    neither = 3 * blob_size * blob_size * (1 - p) * (1 - p) * p
    neither += 2 * blob_size * blob_size * (1 - p) * p

    edges_left = 2 * blob_size * p * center_left
    edges_left += one_left * prob_left
    edges_left += neither * prob_left * 1 / 2

    if print_details:
        print(f"{neither=}, {prob_left=}")

    edges_right = 2 * blob_size * p * center_right
    edges_right += one_right * prob_right
    edges_right += neither * prob_right

    if print_details:
        print(f"Blue edge partition, left: {edges_left}, right {edges_right}")

    return max(edges_left, edges_right)


def blue_vertex_partition(center_size: float, p: float, print_details: bool = False):
    center_proportion_left = 1
    center_left = center_size * center_proportion_left
    center_right = center_size - center_left

    blob_size = (1 - center_size) / 5
    random = 3 * blob_size + 2 * (1 - p) * blob_size
    space_right = 1 / 2 - center_right
    extra_size = random - space_right
    prob_left = extra_size / random
    prob_right = space_right / random

    edges_left = center_left * (1 / 2 - center_left)
    edges_left += 4 * blob_size * blob_size * p * p * prob_left
    edges_neither = 4 * blob_size * blob_size * p * (1 - p) + blob_size * blob_size * p
    edges_left += edges_neither * prob_left * prob_left

    edges_right = center_right * (1 / 2 - center_right)
    edges_right += edges_neither * prob_right * prob_right

    if print_details:
        print(f"Blue vertex partition, left: {edges_left}, right {edges_right}")

    return max(edges_left, edges_right)


def red_vertex_partition(center_size: float, p: float, print_details: bool = False):
    blob_size = (1 - center_size) / 5
    edges_wheel = blob_size * blob_size * 5 * p
    space_left = 1 / 2 - center_size
    prob_left = space_left / (1 - center_size)
    prob_right = (1 / 2) / (1 - center_size)

    edges_left = center_size * (1 / 2 - center_size)
    edges_left += edges_wheel * prob_left * prob_left

    edges_right = edges_wheel * prob_right * prob_right

    if print_details:
        print(f"Red vertex partition, left: {edges_left}, right {edges_right}")

    return max(edges_left, edges_right)


def random_partition(center_size: float, p: float, print_details: bool = False):
    blob_size = (1 - center_size) / 5

    total_edges = center_size * (1 - center_size)
    total_edges += p * blob_size * blob_size * 5

    if print_details:
        print(f"Random partition: {total_edges / 4}")

    return total_edges / 4


def degree(center_size: float, p: float):
    blob_size = (1 - center_size) / 5
    return center_size + 2 * blob_size * p


def partition(args: list[float]):
    center_size, p = args
    if degree(center_size, p) > 1 / 2:
        return 0
    on_blue_edge = blue_edge_partition(center_size, p)
    on_blue_vertex = blue_vertex_partition(center_size, p)
    on_red_vertex = red_vertex_partition(center_size, p)
    random = random_partition(center_size, p)

    return -min(on_blue_edge, on_red_vertex, random)


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
            partition,
            starting_point,
            bounds=bounds,
            method="slsqp",
            options={"ftol": 1e-12, "maxiter": 1000, "eps": 1e-12},
        )
        if -result.fun > largest:
            print(f"New maximum: {-result.fun} at {result.x}")
            print(f"Degree: {degree(result.x[0], result.x[1])}")
            blue_vertex_partition(result.x[0], result.x[1], True)
            blue_edge_partition(result.x[0], result.x[1], True)
            red_vertex_partition(result.x[0], result.x[1], True)
            random_partition(result.x[0], result.x[1], True)
            largest = -result.fun
            largest_witness = result

    return largest_witness


def main():
    result = grid_search(bounds)
    assert result
    print(f"Minimize: worst partition {-result.fun} at {result.x}")

    result = differential_evolution(
        partition,
        bounds=bounds,
    )
    print(f"Differential evolution: worst partition {-result.fun} at {result.x}")


if __name__ == "__main__":
    main()
