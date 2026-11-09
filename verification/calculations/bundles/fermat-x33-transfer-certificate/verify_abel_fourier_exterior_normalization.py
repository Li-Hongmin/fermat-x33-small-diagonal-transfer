"""Check the canonical exterior/Kunneth normalization.

With Alt_4 = (1/4!) sum_sigma sign(sigma) sigma and the ordered wedge
product, the map wedge^4 V -> V^tensor4 -> wedge^4 V is the identity.
Together with c_J=-1 on each H^1 factor, the four-factor Abel--Fourier
scalar is therefore (-1)^4 = 1.
"""

from fractions import Fraction
from itertools import permutations
from math import factorial


def permutation_sign(permutation):
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(len(permutation))
        for j in range(i + 1, len(permutation))
    )
    return -1 if inversions % 2 else 1


def wedge_basis(tensor):
    """Return (canonical mask, sign) for an ordered basis wedge."""
    mask = 0
    sign = 1
    for basis_index in tensor:
        bit = 1 << basis_index
        if mask & bit:
            return None, 0
        sign *= (-1) ** ((mask >> (basis_index + 1)).bit_count())
        mask |= bit
    return mask, sign


def alternating_wedge_coefficient(dimension=4):
    tensor = tuple(range(dimension))
    coefficient = Fraction(0)
    for permutation in permutations(tensor):
        coefficient += Fraction(
            permutation_sign(permutation)
            * wedge_basis(permutation)[1],
            factorial(dimension),
        )
    return coefficient


def main():
    coefficient = alternating_wedge_coefficient()
    assert coefficient == 1
    c_j = -1
    assert c_j**4 == 1
    print("Alt4_to_wedge4 = 1")
    print("c_J^4 = 1")
    print("lambda_Fourier_standard = 1")
    print("all_checks = PASS")


if __name__ == "__main__":
    main()
