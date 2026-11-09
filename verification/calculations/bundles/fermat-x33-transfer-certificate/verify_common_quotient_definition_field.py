"""Exact definition-field obstruction for the X_33 common quotient.

The marked Kummer model is defined over K0 = Q(zeta_165), but the displayed
common-quotient isomorphism must absorb

    c = zeta_5^4 * 5^(-22)

as a 33rd power.  This script verifies with integer arithmetic that c is not a
33rd power in K0, that adjoining a cube root of 5 supplies an explicit 33rd
root, and that the ordinary cubic trace vanishes on every primitive character
modulo 33.

It certifies only this finite arithmetic.  The existence and field of
definition of the algebraic source cycles remain geometric inputs.
"""

from math import gcd


MODULUS = 33
FINITE_BRANCH_EXPONENTS = (7, 28, 13, 19, 10)
INFINITY_EXPONENT = 22


def zeta3_vector(exponent):
    """Return zeta_3**exponent in the basis (1, zeta_3)."""
    residue = exponent % 3
    if residue == 0:
        return (1, 0)
    if residue == 1:
        return (0, 1)
    return (-1, -1)


def add_vectors(*vectors):
    return tuple(sum(entries) for entries in zip(*vectors))


def main():
    weighted_sum = sum(
        index * exponent
        for index, exponent in enumerate(FINITE_BRANCH_EXPONENTS)
    )
    zeta5_exponent = (-weighted_sum) % 5
    assert weighted_sum == 151
    assert zeta5_exponent == 4

    # In K0, every prime above 5 has v_p(5)=4: the zeta_5 factor is totally
    # ramified of degree four and Q(zeta_33) is unramified at 5.
    valuation_of_5 = 4
    valuation_of_c = -INFINITY_EXPONENT * valuation_of_5
    assert valuation_of_c == -88
    assert valuation_of_c % MODULUS != 0

    # Put r^3=5 and d=zeta_5^3*r^(-2).  Then
    # d^33=zeta_5^(99)*r^(-66)=zeta_5^4*5^(-22)=c.
    root_zeta5_exponent = (3 * MODULUS) % 5
    root_r_exponent = -2 * MODULUS
    assert root_zeta5_exponent == zeta5_exponent
    assert root_r_exponent % 3 == 0
    assert root_r_exponent // 3 == -INFINITY_EXPONENT

    primitive_characters = tuple(
        value for value in range(1, MODULUS) if gcd(value, MODULUS) == 1
    )
    assert len(primitive_characters) == 20
    for character in primitive_characters:
        assert character % 3 != 0
        trace = add_vectors(
            zeta3_vector(0),
            zeta3_vector(character),
            zeta3_vector(2 * character),
        )
        assert trace == (0, 0)

    print("weighted_sum =", weighted_sum)
    print("zeta5_exponent =", zeta5_exponent)
    print("valuation_of_c =", valuation_of_c)
    print("valuation_mod_33 =", valuation_of_c % MODULUS)
    print("explicit_root = zeta5^3 * cuberoot(5)^(-2)")
    print("primitive_characters_checked =", len(primitive_characters))
    print("ordinary_cubic_trace = 0")
    print("all_checks = PASS")


if __name__ == "__main__":
    main()
