"""Exact coordinate certificate for the X_33 common quotient.

The calculation uses Q[zeta_5] = Q[zeta]/(zeta^4+zeta^3+zeta^2+zeta+1).
For P(s,t)=sum u_k s^(4-k)t^k, it verifies

    Z_i = zeta^i P(1,zeta^i),   Z_5=-5 P(0,1),
    sum_i Z_i = 0,

and checks that the first five Z_i give an invertible coordinate system.
The remaining checks are the exact Kummer exponent action, the permutation
which produces alpha=(7,10,13,19,22,28), and the scalar root used to identify
the two degree-33 function fields.
"""

from itertools import permutations
from math import gcd


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def neg(x):
    return tuple(-value for value in x)


def mul(x, y):
    """Multiply in Q[zeta_5], basis 1,zeta,zeta^2,zeta^3."""
    raw = [0] * 7
    for i, left in enumerate(x):
        for j, right in enumerate(y):
            raw[i + j] += left * right
    # zeta^4 = -(1+zeta+zeta^2+zeta^3), then reduce powers 6,5,4.
    for degree in range(6, 3, -1):
        coefficient = raw[degree]
        if coefficient:
            for lower in range(4):
                raw[degree - 4 + lower] -= coefficient
    return tuple(raw[:4])


def power(base, exponent):
    result = (1, 0, 0, 0)
    for _ in range(exponent):
        result = mul(result, base)
    return result


def determinant(matrix):
    size = len(matrix)
    total = (0, 0, 0, 0)
    for permutation in permutations(range(size)):
        inversions = sum(
            permutation[i] > permutation[j]
            for i in range(size)
            for j in range(i + 1, size)
        )
        term = (1, 0, 0, 0)
        for row, column in enumerate(permutation):
            term = mul(term, matrix[row][column])
        total = add(total, neg(term) if inversions % 2 else term)
    return total


def main():
    zeta = (0, 1, 0, 0)
    powers = [power(zeta, i) for i in range(5)]

    # Coefficient matrix for Z_i = zeta^i P(1,zeta^i), i=0,...,4.
    matrix = [
        [powers[(i * (k + 1)) % 5] for k in range(5)]
        for i in range(5)
    ]
    coordinate_determinant = determinant(matrix)
    assert coordinate_determinant != (0, 0, 0, 0)

    # The sixth coordinate is Z_5=-5u_4; its row is minus the sum of the first
    # five rows, exactly because sum_i zeta^{i(k+1)}=0 for k<4 and =5 for k=4.
    row_sum = [(0, 0, 0, 0) for _ in range(5)]
    for row in matrix:
        for column, value in enumerate(row):
            row_sum[column] = add(row_sum[column], value)
    assert row_sum == [(0, 0, 0, 0)] * 4 + [(5, 0, 0, 0)]

    kummer_exponents = (7, 28, 13, 19, 10, 22)
    alpha = (7, 10, 13, 19, 22, 28)
    reorder = (0, 4, 2, 3, 5, 1)
    assert tuple(kummer_exponents[i] for i in reorder) == alpha
    assert sum(alpha) % 33 == 0
    assert gcd(33, *alpha) == 1

    # Product normalization: Z_5=-5L_infinity contributes (-5)^22=5^22,
    # while the other five coordinates contribute zeta_5^151=zeta_5.
    coordinate_indices = (0, 4, 2, 3, None, 1)
    zeta5_product_exponent = sum(
        index * exponent
        for index, exponent in zip(coordinate_indices, alpha)
        if index is not None
    )
    assert zeta5_product_exponent == 151
    assert alpha[4] % 2 == 0
    kappa = (-zeta5_product_exponent % 5, -alpha[4])
    assert kappa == (4, -22)

    # lambda=zeta_5^3*r^-2 with r^3=5 satisfies lambda^33=kappa.
    lambda_zeta5_exponent, lambda_r_exponent = 3, -2
    assert 33 * lambda_zeta5_exponent % 5 == kappa[0]
    assert 33 * lambda_r_exponent % 3 == 0
    assert 33 * lambda_r_exponent // 3 == kappa[1]

    # The six branch hyperplanes have inertia m/gcd(m, alpha_i).  The
    # z=infinity branch is alpha[4]=22 and has inertia 3.  The auxiliary
    # denominator hyperplane ell=0 has valuation -sum(alpha), divisible by 33.
    local_inertia_orders = tuple(33 // gcd(33, exponent) for exponent in alpha)
    geometric_branch_fiber_cardinalities = tuple(gcd(33, exponent) for exponent in alpha)
    auxiliary_hyperplane_inertia_order = 33 // gcd(33, sum(alpha))
    special_inertia_generator = 11  # <delta^11> has order 3 at alpha_4=22.
    assert (4 * special_inertia_generator) % 33 == special_inertia_generator
    assert local_inertia_orders == (33, 33, 33, 33, 3, 33)
    assert geometric_branch_fiber_cardinalities == (1, 1, 1, 1, 11, 1)
    assert auxiliary_hyperplane_inertia_order == 1
    print("local_inertia_orders =", local_inertia_orders)
    print("geometric_branch_fiber_cardinalities =", geometric_branch_fiber_cardinalities)
    print("auxiliary_hyperplane_inertia_order =", auxiliary_hyperplane_inertia_order)
    print("rho_stable_special_inertia =", True)

    # At the prime divisor Z_0=0 of the common P^4 base, the Kummer
    # function has valuation alpha[0]=7.  Hence it is neither a cube nor
    # an eleventh power in C(P^4)^*, so T^33-f is irreducible over the
    # field containing mu_33 (the standard coprime-prime Kummer criterion).
    valuation_at_Z0 = alpha[0]
    assert valuation_at_Z0 % 3 != 0
    assert valuation_at_Z0 % 11 != 0
    print("valuation_at_Z0 =", valuation_at_Z0)
    print("kummer_degree = 33")

    rho_image = tuple((4 * value) % 33 for value in kummer_exponents)
    assert rho_image == (28, 13, 19, 10, 7, 22)
    assert rho_image[5] == kummer_exponents[5]
    assert tuple((4**5 * value) % 33 for value in kummer_exponents) == kummer_exponents

    print("coordinate_determinant_Q[zeta5] =", coordinate_determinant)
    print("sum_Z_coefficients =", row_sum)
    print("reorder_to_alpha =", reorder)
    print("kappa_exponents_(zeta5,5) =", kappa)
    print("lambda^33_matches_kappa =", True)
    print("rho_exponent_action =", rho_image)
    print("gcd(alpha,33) =", gcd(33, *alpha))
    print("all_checks = PASS")


if __name__ == "__main__":
    main()
