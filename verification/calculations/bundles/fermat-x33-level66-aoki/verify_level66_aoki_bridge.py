"""Exact finite arithmetic for the level-66 Aoki bridge and xi_33 lift.

The geometric implications are the cited Aoki/da Silva results and the
Lefschetz (1,1) theorem, not claims proved by this script.
"""

from collections import Counter
from math import gcd


M = 66
UNITS = tuple(t for t in range(1, M) if gcd(t, M) == 1)

w = (1, 4, 16, 22, 25, 31)
alpha = (7, 10, 13, 19, 22, 28)
a = tuple((2 * x) % M for x in w)
xi33 = (1, 4, 13, 16, 19, 22, 27, 30)
xi33_lift = tuple(2 * x for x in xi33)
Q = (1, 25, 44, 62)
S = (2, 8, 32, 41, 50, 65)
delta0 = (1, 65, 25, 41)
beta = (2, 32, 33, 65)
gamma = (8, 33, 41, 50)
beta_da_silva = (2, 32, 65, 33)
gamma_da_silva = (8, 41, 50, 33)


def hodge_sum(character):
    """Return the numerator of sum_i <t*a_i/M> for every unit t."""
    return tuple(sum((t * x) % M for x in character) for t in UNITS)


def lineal(i):
    return (i, (-i) % M)


def standard(p, i):
    k = M // p
    if p == 2:
        return (i, (i + k) % M, (M - 2 * i) % M, k)
    return tuple((i + j * k) % M for j in range(p)) + ((M - p * i) % M,)


assert len(UNITS) == 20
assert hodge_sum(a) == (3 * M,) * len(UNITS)       # fourfold, (2,2)
for surface_character in (Q, beta, gamma):
    assert hodge_sum(surface_character) == (2 * M,) * len(UNITS)
assert hodge_sum(S) == (3 * M,) * len(UNITS)

assert Counter(a) + Counter(delta0) == Counter(Q) + Counter(S)
assert Counter(beta_da_silva) == Counter(beta)
assert Counter(gamma_da_silva) == Counter(gamma)
assert (beta_da_silva[-1] + gamma_da_silva[-1]) % M == 0
assert Counter(beta_da_silva[:-1] + gamma_da_silva[:-1]) == Counter(S)

# Miranda--Movasati--Rufino--Villaflor's D+S generators give an exact lift
# identity; Q is algebraic by Lefschetz (1,1) after the age check above.
assert xi33_lift == (2, 8, 26, 32, 38, 44, 54, 60)
assert standard(2, 58) == (58, 25, 16, 33)
assert standard(2, 65) == (65, 32, 2, 33)
assert standard(3, 38) == (38, 60, 16, 18)
assert standard(3, 62) == (62, 18, 40, 12)
assert (
    Counter(xi33_lift)
    + Counter(lineal(1))
    + Counter(standard(2, 58))
    + Counter(standard(3, 62))
    == Counter(Q)
    + Counter(lineal(8))
    + Counter(lineal(12))
    + Counter(lineal(26))
    + Counter(standard(2, 65))
    + Counter(standard(3, 38))
)


def is_self_pair(x, y):
    return (x + y) % M == 0


assert is_self_pair(1, 65)
assert is_self_pair(25, 41)

# The target orbit is reached by the Fermat coordinate permutation after 7w.
assert sorted((7 * x) % 33 for x in w) == sorted(alpha)

# The single eigenspace needed by the B_8 source-line criterion is the 59th
# conjugate; its level-66 character descends to -alpha up to permutation.
a59 = tuple((59 * x) % M for x in a)
c = tuple((-7 * x) % 33 for x in w)
minus_alpha = tuple((-x) % 33 for x in alpha)
assert a59 == tuple((-14 * x) % M for x in w)
assert a59 == tuple(2 * x for x in c)
assert sorted(c) == sorted(minus_alpha)
assert all(x % 2 == 0 for x in a59)  # the quotient kernel acts trivially

# For pi : X_66^4 -> X_33^4, the quotient is mu_2^6/mu_2,diag.
quotient_signs = {
    tuple(e ^ bits[0] for e in bits)
    for bits in __import__("itertools").product((0, 1), repeat=6)
}
assert len(quotient_signs) == 32
# A sign change acts on a character by (-1)^(e_i*a_i); a is even coordinatewise.
for signs in quotient_signs:
    assert sum(sign * coordinate for sign, coordinate in zip(signs, a)) % 2 == 0

print("units mod 66:", len(UNITS))
print("Hodge checks: a,S=(2,2) on X_66^4; Q,beta,gamma=(1,1) on surfaces")
print("da Silva contraction: beta',gamma' drop distinguished 33,33 to S")
print("Aoki elimination multiset: a+delta0=Q+S")
print("xi_33 lift identity: exact D+S multiset equality with Q")
print("vanishing pairs in delta0: valid")
print("quotient group size:", len(quotient_signs))
print("target orbit: 7*w = alpha up to permutation")
print("single source line: 59*a = 2*(-7*w), and -7*w = -alpha up to permutation")
print("ALL ASSERTIONS PASSED")
