import itertools
from dataclasses import dataclass


@dataclass(frozen=True)
class Flag:
    size: int
    prob: float
    edges: dict[tuple[int, int], int]
    root_size: int


def read_flags(f):
    for line in f:
        data = line.split()
        prob = float(data[0])
        size = int(data[1])
        root_size = int(data[2])
        edges = {}
        for left, right in itertools.product(range(size), repeat=2):
            edges[(left, right)] = 0

        idx = 3
        for i in range(0, size - 1):
            for j in range(i + 1, size):
                edges[(i, j)] = int(data[idx]) - 1
                edges[(j, i)] = int(data[idx]) - 1
                idx += 1

        yield Flag(size, prob, edges, root_size)
