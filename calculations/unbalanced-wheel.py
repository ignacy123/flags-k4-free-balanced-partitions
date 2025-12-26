from scipy.optimize import minimize, Bounds, LinearConstraint
import numpy as np
from itertools import product


def half(x, a, b, c):
    return x * a * b + (1 / 2 - x * a - b - c) * (x * a + b + c)


def degree_partition(a, b, c, d, e):
    edges_in_c5 = a * b + b * c + c * d + d * e + e * a
    total_c5_size = a + b + c + d + e
    edges_between_c5_and_rest = (1 / 2 - total_c5_size) * total_c5_size

    return edges_in_c5 + edges_between_c5_and_rest


def red_edge_partition(a, b, c, d, e, print_info=False):
    # x = (
    #     -a
    #     + 2 * a**2
    #     + b
    #     - 2 * b**2
    #     - c
    #     + 4 * a * c
    #     + 2 * c**2
    #     + d
    #     - 4 * b * d
    #     - 2 * d**2
    #     - e
    #     + 2 * a * e
    #     + 4 * c * e
    #     + 2 * e**2
    # ) / (2 * a * (-1 + 2 * a + b + 2 * c + 2 * d + e))
    x = 1 / 2
    left = half(x, a, b, d)
    right = half(1 - x, a, e, c)
    middle_point = max(left, right)
    if print_info:
        print(f"{left=}, {right=}")
    # left_point = max(half(0, a, b, d), half(1, a, e, c))
    # right_point = max(half(1, a, b, d), half(0, a, e, c))

    return middle_point


def optimistic_cut(a, b, c, d, e):
    p_a = red_edge_partition(a, b, c, d, e)
    p_b = red_edge_partition(b, c, d, e, a)
    p_c = red_edge_partition(c, d, e, a, b)
    p_d = red_edge_partition(d, e, a, b, c)
    p_e = red_edge_partition(e, a, b, c, d)

    return min(p_a, p_b, p_c, p_d, p_e)


def projected_cut_by_vertex(a, b, c, d, e, print_info=False):
    left_side_with_e = half(1 / 2, c, d, a)
    right_side_with_e = half(1 / 2, c, b, e)
    left_side_with_b = half(1 / 2, d, c, a)
    right_side_with_b = half(1 / 2, d, e, b)

    averaged_left_side = (e * left_side_with_e + b * left_side_with_b) / (e + b)
    averaged_right_side = (e * right_side_with_e + b * right_side_with_b) / (e + b)

    if print_info:
        print(f"{averaged_left_side=}, {averaged_right_side=}")
        print(f"{left_side_with_b=}, {right_side_with_b=}")
        print(f"{left_side_with_e=}, {right_side_with_e=}")

    return averaged_left_side, averaged_right_side


def random_red_edge_partition(a, b, c, d, e, print_info=False):
    left_a, right_a = projected_cut_by_vertex(a, b, c, d, e, print_info)
    left_b, right_b = projected_cut_by_vertex(b, c, d, e, a, print_info)
    left_c, right_c = projected_cut_by_vertex(c, d, e, a, b, print_info)
    left_d, right_d = projected_cut_by_vertex(d, e, a, b, c, print_info)
    left_e, right_e = projected_cut_by_vertex(e, a, b, c, d, print_info)
    red = []
    magenta = []

    for size, left, right in [
        (a, left_a, right_a),
        (b, left_b, right_b),
        (c, left_c, right_c),
        (d, left_d, right_d),
        (e, left_e, right_e),
    ]:
        if left >= right:
            red.append((size, left))
        else:
            magenta.append((size, right))

    red_averaged = 0
    red_total_size = 0
    for size, value in red:
        red_total_size += size
        red_averaged += size * value

    if red_total_size > 0:
        red_averaged /= red_total_size
    else:
        red_averaged = 1

    magenta_averaged = 0
    magenta_total_size = 0
    for size, value in magenta:
        magenta_total_size += size
        magenta_averaged += size * value

    if magenta_total_size > 0:
        magenta_averaged /= magenta_total_size
    else:
        magenta_averaged = 1

    if print_info:
        print(f"{left_a=}, {left_b=}, {left_c=}, {left_d=}, {left_e=}")
        print(f"{right_a=}, {right_b=}, {right_c=}, {right_d=}, {right_e=}")
        print(f"{red_total_size=}, {magenta_total_size=}")
        print(f"{red_averaged=}, {magenta_averaged=}")

    return min(red_averaged, magenta_averaged)


def red_vertex_partition(owning_blob, first_opposite, second_opposite):
    deterministic_part_size = owning_blob + first_opposite + second_opposite
    edges_from_center = (1 / 2 - deterministic_part_size) * deterministic_part_size

    return first_opposite * second_opposite + edges_from_center


def random_red_vertex_partition(a, b, c, d, e, print_info=False):
    p_a = red_vertex_partition(a, c, d)
    p_b = red_vertex_partition(b, d, e)
    p_c = red_vertex_partition(c, e, a)
    p_d = red_vertex_partition(d, a, b)
    p_e = red_vertex_partition(e, b, c)

    if print_info:
        print(f"{p_a=}, {p_b=}, {p_c=}, {p_d=}, {p_e=}")

    return (a * p_a + b * p_b + c * p_c + d * p_d + e * p_e) / (a + b + c + d + e)


def random_cut(a, b, c, d, e):
    center_size = 1 - a - b - c - d - e
    red_edges = a * b + b * c + c * d + d * e + e * a
    total_edges = center_size * (1 - center_size) + red_edges

    return total_edges / 4


def tested_partition(args):
    a, b, c, d, e = args
    p_0 = degree_partition(a, b, c, d, e)
    edge_partition = random_red_edge_partition(a, b, c, d, e)
    vertex_partition = random_red_vertex_partition(a, b, c, d, e)
    random_partition = random_cut(a, b, c, d, e)
    optimistic_partition = optimistic_cut(a, b, c, d, e)

    return -min(p_0, optimistic_partition, vertex_partition, random_partition)


def print_info(args):
    a, b, c, d, e = args
    red_edge_density = a * b + b * c + c * d + d * e + e * a
    red_edge_density *= 2
    print(f"{red_edge_density=}")
    red_vertex_density = a + b + c + d + e
    blue_vertex_density = 1 - red_vertex_density
    print(f"{blue_vertex_density=}")
    print(f"{red_vertex_density=}")
    degrees = []
    for i in range(5):
        left_neighbor, right_neighbor = args[i - 1], args[(i + 1) % 5]
        degree = left_neighbor + right_neighbor + blue_vertex_density
        degrees.append(degree)

    print(f"{degrees=}")

    deg_partition = degree_partition(a, b, c, d, e)
    print(f"{deg_partition=}")

    edge_partition = random_red_edge_partition(a, b, c, d, e, True)
    print(f"{edge_partition=}")

    vertex_partition = random_red_vertex_partition(a, b, c, d, e, True)
    print(f"{vertex_partition=}")

    random_partition = random_cut(a, b, c, d, e)
    print(f"{random_partition=}")

    optimistic = optimistic_cut(a, b, c, d, e)
    print(f"{optimistic=}")

    print()


def grid_search(bounds: Bounds):
    precision = 4
    linspaces = []
    sanity = LinearConstraint([[1, 1, 1, 1, 1]], 1 / 10, 1 / 2)
    largest, largest_witness = 0, None
    for lower, upper in zip(bounds.lb, bounds.ub):
        linspaces.append(np.linspace(lower, upper, precision))

    for idx, starting_point in enumerate(product(*linspaces)):
        if idx % 10000 == 0:
            print(f"Iteration {idx}")
        result = minimize(
            tested_partition,
            starting_point,
            bounds=bounds,
            constraints=[sanity],
            options={"ftol": 1e-12, "maxiter": 1000, "eps": 1e-12},
        )
        if not result.success:
            continue
        if -result.fun > largest:
            print(f"New maximum: {-result.fun} at {result.x}")
            print_info(result.x)
            largest = -result.fun
            largest_witness = result

    return largest_witness


def main():
    lower_bound = 1 / 60
    bounds = Bounds([lower_bound] * 5, [1 / 2, 1 / 2, 1 / 2, 1 / 2, 1 / 2])
    grid_search(bounds)


if __name__ == "__main__":
    main()
