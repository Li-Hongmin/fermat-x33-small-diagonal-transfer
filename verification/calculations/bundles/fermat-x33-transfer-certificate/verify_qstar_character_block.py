#!/usr/bin/env python3
"""Check the orbit-average versus descended-exterior normalization.

On an order-five orbit, the source projector averages the five determinant
blocks simultaneously.  Descent of the exterior of the averaged H^1 vectors
has coefficient 1/5^4 in each block, whereas the source average has 1/5.
Thus the canonical descended exterior is 5^3 times the source vector, and
the raw tensor correspondence is 5^-3 times the descended exterior map.
This checks the finite scalar only; it is not a geometric proof.
"""

from fractions import Fraction
from itertools import permutations


DIMENSION = 4
ORBIT_LENGTH = 5


def permutation_sign(permutation):
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def tensor_index(indices):
    value = 0
    for index in indices:
        value = value * DIMENSION + index
    return value


def alt_basis():
    result = {}
    for permutation in permutations(range(DIMENSION)):
        result[tensor_index(permutation)] = Fraction(permutation_sign(permutation))
    return result


def add(left, right):
    result = dict(left)
    for index, coefficient in right.items():
        result[index] = result.get(index, Fraction(0)) + coefficient
    return {index: coefficient for index, coefficient in result.items() if coefficient}


def apply_tensor(matrix, vector):
    result = {}
    for source_index, source_coefficient in vector.items():
        digits = []
        remainder = source_index
        for _ in range(DIMENSION):
            digits.append(remainder % DIMENSION)
            remainder //= DIMENSION
        digits.reverse()
        for targets in __import__("itertools").product(range(DIMENSION), repeat=DIMENSION):
            coefficient = source_coefficient
            for target, source in zip(targets, digits):
                coefficient *= matrix[target][source]
            if coefficient:
                index = tensor_index(targets)
                result[index] = result.get(index, Fraction(0)) + coefficient
    return {index: coefficient for index, coefficient in result.items() if coefficient}


def main():
    matrix = (
        (Fraction(1), Fraction(1), Fraction(0), Fraction(0)),
        (Fraction(0), Fraction(1), Fraction(1), Fraction(0)),
        (Fraction(0), Fraction(0), Fraction(1), Fraction(1)),
        (Fraction(0), Fraction(0), Fraction(0), Fraction(1)),
    )
    alt = alt_basis()

    # det(matrix)=1, so the fourth exterior determinant is unchanged.
    image = apply_tensor(matrix, alt)
    assert image == alt

    # Source projector: one simultaneous average of determinant blocks.
    source_block_coefficient = Fraction(1, ORBIT_LENGTH)
    descended_block_coefficient = Fraction(1, ORBIT_LENGTH**DIMENSION)
    descent_scale = source_block_coefficient / descended_block_coefficient
    assert descent_scale == ORBIT_LENGTH ** (DIMENSION - 1)

    # Each H^1 invariant is the average of its five orbit components. A map
    # identical on each component sends that average to one copy, not five.
    tensor_transfer = sum(
        (source_block_coefficient for _ in range(ORBIT_LENGTH)), Fraction(0)
    )
    exterior_transfer = descent_scale * tensor_transfer
    assert tensor_transfer == 1
    assert exterior_transfer == ORBIT_LENGTH ** (DIMENSION - 1)
    assert tensor_transfer == Fraction(1, descent_scale) * exterior_transfer

    print("determinant_block_dim = 1")
    print("source_projector_average = 1/5")
    print("canonical_descent_relative_scale = 5^3")
    print("small_diagonal_to_descended_exterior = 1/125")
    print("orbit_descent_normalization = PASS")


if __name__ == "__main__":
    main()
