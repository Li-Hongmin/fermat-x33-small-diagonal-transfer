"""Exact finite arithmetic for the X_33^4 transfer chain.

No external packages are needed.  The checks are intentionally limited to
discrete data; geometric assertions remain in the accompanying mathematical
notes.
"""

from math import gcd


MODULUS = 33
ALPHA = (7, 10, 13, 19, 22, 28)

# The five finite branch exponents of the Kummer curve C, followed by infinity.
KUMMER_EXPONENTS = (7, 28, 13, 19, 10, 22)


def mul_vector(vector, scalar, modulus=MODULUS):
    return tuple((scalar * value) % modulus for value in vector)


def main():
    assert sum(ALPHA) % MODULUS == 0
    assert gcd(ALPHA[0], MODULUS) == 1

    # z -> zeta_5 z sends the finite branch vector to 4 times itself;
    # infinity is fixed because 4*22 == 22 mod 33.
    finite = KUMMER_EXPONENTS[:5]
    assert mul_vector(finite, 4) == (28, 13, 19, 10, 7)
    assert (4 * KUMMER_EXPONENTS[5]) % MODULUS == KUMMER_EXPONENTS[5]
    assert mul_vector(KUMMER_EXPONENTS, 4) == (
        28,
        13,
        19,
        10,
        7,
        22,
    )

    # Reorder (Z_0,Z_1,Z_2,Z_3,Z_4,Z_5) as
    # (Z_0,Z_4,Z_2,Z_3,Z_5,Z_1).
    reorder = (0, 4, 2, 3, 5, 1)
    assert tuple(KUMMER_EXPONENTS[index] for index in reorder) == ALPHA

    # N = ker(sum:(Z/33)^4 -> Z/33), and A=(Z/33)^6/diagonal.
    n_size = MODULUS**3
    s4_size = 24
    g_size = n_size * s4_size
    a_size = MODULUS**5
    k_size = MODULUS**4
    source_scalar = g_size * k_size

    assert n_size == 35_937
    assert g_size == 862_488
    assert a_size == 33**5
    assert k_size == 1_185_921
    assert source_scalar == 1_022_842_631_448

    # The exact-conductor character orbit under multiplication by 4 has length 5.
    orbit = []
    value = 1
    for _ in range(5):
        orbit.append(value)
        value = (4 * value) % MODULUS
    assert orbit == [1, 4, 16, 31, 25]
    assert value == 1
    assert len(set(orbit)) == 5

    print("alpha_sum_mod_33 =", sum(ALPHA) % MODULUS)
    print("kummer_to_alpha_reorder =", reorder)
    print("rho_exponent_action =", mul_vector(KUMMER_EXPONENTS, 4))
    print("|N| =", n_size)
    print("|N ⋊ S4| =", g_size)
    print("|K_alpha| =", k_size)
    print("source_scalar =", source_scalar)
    print("<4>-orbit(1) =", orbit)
    print("all_checks = PASS")


if __name__ == "__main__":
    main()
