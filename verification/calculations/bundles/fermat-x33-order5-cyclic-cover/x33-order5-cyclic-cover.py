#!/usr/bin/env python3
"""Exact order-five symmetry attached to the X_33^4 difficult tuple.

The tuple alpha=(7,10,13,19,22,28) is stable as a multiset under
multiplication by 4 modulo 33.  We realize the six entries as the branch
exponents of a cyclic 33-cover of P^1 with finite branch points mu_5 and
infinity, and verify the induced order-five lift on the function field.

The certificate also computes the conductor-33 Hodge structure, its order-five
descent, and the linear hyperplane relation behind the common cyclic quotient
with X_33^4.  It does not prove algebraicity of the resulting generalized Weil
classes and does not invoke Schoen's simple-tuple theorem.
"""

from fractions import Fraction
from math import factorial, gcd


M = 33
ORDER5 = 5
ALPHA = (7, 10, 13, 19, 22, 28)
FINITE_EXPONENTS = (7, 28, 13, 19, 10)
INFINITY_EXPONENT = 22
MULTIPLIER = 4
BRANCH_EXPONENTS = FINITE_EXPONENTS + (INFINITY_EXPONENT,)
STABILIZER = (1, 4, 16, 25, 31)
# The branch order is (1,zeta_5,zeta_5^2,zeta_5^3,zeta_5^4,infinity).
# Reordering by this permutation gives the standard alpha coordinate order.
FERMAT_COORDINATE_PERMUTATION = (0, 4, 2, 3, 5, 1)

# g(z)=prod_i (z-zeta^i)^k_i in the identity
# f(zeta z)=zeta^77*g(z)^33*f(z)^4.
K_EXPONENTS = (0, -3, -1, -2, -1)
WEIGHTS = tuple(MULTIPLIER ** (ORDER5 - 1 - j) for j in range(ORDER5))


def least(value):
    return (value - 1) % M + 1


def multiset_translate(values, unit):
    return tuple(sorted(least(unit * value) for value in values))


def hodge_condition(values):
    expected = M * len(values) // 2
    units = tuple(u for u in range(1, M) if gcd(u, M) == 1)
    return all(sum(least(u * value) for value in values) == expected for u in units)


def permutation_cycles(values, unit):
    remaining = set(values)
    cycles = []
    while remaining:
        start = min(remaining)
        cycle = [start]
        remaining.remove(start)
        value = least(unit * start)
        while value != start:
            cycle.append(value)
            remaining.remove(value)
            value = least(unit * value)
        cycles.append(tuple(cycle))
    return tuple(cycles)


def primitive_units(modulus):
    return tuple(value for value in range(1, modulus) if gcd(value, modulus) == 1)


def character_hodge_split(character):
    """Integral Chevalley--Weil split for a nontrivial deck character.

    For a cyclic cover of P^1 with branch exponents a_i, Schoen's formula is
        -1 + sum_i <character*a_i/M>,
    and the conjugate multiplicity uses -character.  Here fractional parts use
    residues in {0,...,M-1}; this matters for imprimitive characters.
    """
    positive_sum = sum((character * exponent) % M for exponent in BRANCH_EXPONENTS)
    negative_sum = sum((-character * exponent) % M for exponent in BRANCH_EXPONENTS)
    assert positive_sum % M == 0
    assert negative_sum % M == 0
    split = (-1 + positive_sum // M, -1 + negative_sum // M)
    assert min(split) >= 0
    return split


def exact_character_order(character):
    return M // gcd(character, M)


def characters_of_exact_order(order):
    return tuple(
        character
        for character in range(1, M)
        if exact_character_order(character) == order
    )


def character_cycles(characters):
    """Cycles induced by rho on a rho-stable set of deck characters."""
    remaining = set(characters)
    rho_multiplier = pow(MULTIPLIER, -1, M)
    cycles = []
    while remaining:
        start = min(remaining)
        cycle = [start]
        remaining.remove(start)
        value = least(rho_multiplier * start)
        while value != start:
            cycle.append(value)
            remaining.remove(value)
            value = least(rho_multiplier * value)
        cycles.append(tuple(cycle))
    return tuple(cycles)


def character_cycles_under_rho():
    """Cycles on primitive deck characters induced by rho.

    From rho delta rho^{-1}=delta^4, rho sends the delta-character n to
    4^{-1} n = 25 n.  The inverse has the same order-five orbits as 4.
    """
    return character_cycles(primitive_units(M))


def u_mu5_decomposition():
    """Exact abstract decomposition of Schoen's U under the order-five lift.

    For this six-branch cover h=-e(P^1\\B)=4.  Each primitive deck character
    occurs with multiplicity h=4 in H^1(C), with Hodge split (2,2), and the
    corresponding Schoen weight space in U subset H^4(C^4) is det(V_n), hence
    one-dimensional.  rho permutes these one-dimensional spaces in the
    character cycles above.  Since rho^5=1, every five-cycle is a regular
    mu_5 representation (a change of basis removes the nonzero edge scalars).
    """
    cycles = character_cycles_under_rho()
    primitive_splits = tuple(
        character_hodge_split(character) for character in primitive_units(M)
    )
    assert set(primitive_splits) == {(2, 2)}
    assert len(cycles) == 4 and all(len(cycle) == ORDER5 for cycle in cycles)
    # Each regular representation contributes one copy of every fifth root.
    eigenvalue_multiplicities = tuple(len(cycles) for _ in range(ORDER5))
    return cycles, eigenvalue_multiplicities


def h1_conductor_decomposition(quotient_h1_dimension):
    """H^1 dimensions by exact deck-character order and rho-invariant ranks."""
    h1_dimensions = {}
    invariant_dimensions = {}
    cycles_by_order = {}
    for order in (3, 11, 33):
        characters = characters_of_exact_order(order)
        cycles = character_cycles(characters)
        cycles_by_order[order] = cycles
        dimensions = {
            character: sum(character_hodge_split(character))
            for character in characters
        }
        h1_dimensions[order] = sum(dimensions.values())
        if all(len(cycle) == ORDER5 for cycle in cycles):
            assert all(
                len({dimensions[character] for character in cycle}) == 1
                for cycle in cycles
            )
            # rho^5=1 turns each five-cycle of equal-dimensional spaces into
            # that many copies of the regular representation.
            invariant_dimensions[order] = sum(
                dimensions[cycle[0]] for cycle in cycles
            )

    # The order-three characters are fixed by multiplication by 25, so their
    # internal rho eigenvalues are not determined by character permutation.
    # Riemann--Hurwitz gives the total invariant H^1 of C/<rho>; subtract the
    # two long-orbit conductor pieces to determine the remaining rank.
    invariant_dimensions[3] = quotient_h1_dimension - sum(
        invariant_dimensions[order] for order in (11, 33)
    )
    return h1_dimensions, invariant_dimensions, cycles_by_order


def root_of_unity_sum(exponent, order=ORDER5):
    """Return sum_{i=0}^{order-1} zeta_order^(i*exponent)."""
    return order if exponent % order == 0 else 0


def conductor33_projector_value(character):
    """Value of the group-algebra projector e_33 on a deck character."""
    order3_average = int((11 * character) % M == 0)
    order11_average = int((3 * character) % M == 0)
    return (1 - order3_average) * (1 - order11_average)


def rho5_exponent_check():
    """Return finite-branch exponents and the root-of-unity constant in rho^5."""
    # Under z -> zeta^j z, the exponent of (z-zeta^ell) in g is k_{ell+j}.
    branch_exponents = []
    for ell in range(ORDER5):
        exponent = sum(
            WEIGHTS[j] * K_EXPONENTS[(ell + j) % ORDER5]
            for j in range(ORDER5)
        ) + 31 * FINITE_EXPONENTS[ell]
        branch_exponents.append(exponent)

    # zeta_165 is a primitive 165-th root, zeta_5=zeta_165^33, and we choose
    # c=zeta_165^12.  The exponent is the constant in rho^5 after using
    # y^1024=y^(1+31*33)=y*f(z)^31.
    c_exponent = 12
    zeta_exponent = sum(
        WEIGHTS[j] * j * sum(K_EXPONENTS)
        for j in range(ORDER5)
    )
    constant_exponent_mod_165 = (
        c_exponent * sum(WEIGHTS) + 33 * zeta_exponent
    ) % 165
    return tuple(branch_exponents), constant_exponent_mod_165


def affine_fixed_count(modulus, multiplier, translation, power):
    """Count fixed points of a power of x -> multiplier*x+translation."""
    translated = translation * sum(multiplier ** j for j in range(power))
    powered_multiplier = multiplier ** power
    return sum(
        1
        for x in range(modulus)
        if ((powered_multiplier - 1) * x + translated) % modulus == 0
    )


assert tuple(sorted(ALPHA)) == tuple(sorted(FINITE_EXPONENTS + (INFINITY_EXPONENT,)))
assert hodge_condition(ALPHA)
assert multiset_translate(ALPHA, MULTIPLIER) == tuple(sorted(ALPHA))
assert STABILIZER == tuple(sorted({pow(MULTIPLIER, i, M) for i in range(ORDER5)}))
assert all(multiset_translate(ALPHA, unit) == tuple(sorted(ALPHA)) for unit in STABILIZER)
assert all(multiset_translate(ALPHA, unit) != tuple(sorted(ALPHA))
           for unit in range(1, M) if gcd(unit, M) == 1 and unit not in STABILIZER)
assert permutation_cycles(ALPHA, MULTIPLIER) == ((7, 28, 13, 19, 10), (22,))
assert all((FINITE_EXPONENTS[(i + 1) % ORDER5] - MULTIPLIER * FINITE_EXPONENTS[i]) % M == 0
           for i in range(ORDER5))
assert tuple(sum(K_EXPONENTS[(ell + j) % ORDER5] * WEIGHTS[j]
                for j in range(ORDER5)) + 31 * FINITE_EXPONENTS[ell]
             for ell in range(ORDER5)) == (0, 0, 0, 0, 0)
branch_exponents, constant_exponent_mod_165 = rho5_exponent_check()
assert branch_exponents == (0, 0, 0, 0, 0)
assert constant_exponent_mod_165 == 0

# The five finite branch exponents are units, but gcd(33,22)=11 at infinity.
# For y^n=f(z), a branch exponent a contributes n-gcd(n,a) to
# Riemann--Hurwitz and has gcd(n,a) points above the branch point.
branch_gcds = tuple(gcd(M, exponent) for exponent in BRANCH_EXPONENTS)
ramification_contributions = tuple(M - divisor for divisor in branch_gcds)
two_genus_c_minus_two = -2 * M + sum(ramification_contributions)
assert two_genus_c_minus_two % 2 == 0
genus_c = 1 + two_genus_c_minus_two // 2
points_over_infinity = branch_gcds[-1]
infinity_ramification_index = M // points_over_infinity
infinity_affine_multiplier = MULTIPLIER % points_over_infinity
assert branch_gcds == (1, 1, 1, 1, 1, 11)
assert ramification_contributions == (32, 32, 32, 32, 32, 22)
assert genus_c == 59
assert points_over_infinity == 11
assert infinity_ramification_index == 3
assert infinity_affine_multiplier == 4

# Over z=0, rho acts affinely on the unramified Z/33 torsor.  The relation
# rho^5=1 forces 11*b=0 mod 33, hence b is divisible by 3.  Every nontrivial
# power then has exactly three fixed points.  At infinity the 11 points form
# a Z/11 torsor; every affine map with multiplier 4 has one fixed point for
# powers 1,...,4 because 4^k-1 is invertible modulo 11.
affine_order_sum = sum(MULTIPLIER ** j for j in range(ORDER5))
zero_translations = tuple(
    translation
    for translation in range(M)
    if affine_order_sum * translation % M == 0
)
zero_fixed_counts = {
    affine_fixed_count(M, MULTIPLIER, translation, power)
    for translation in zero_translations
    for power in range(1, ORDER5)
}
infinity_fixed_counts = {
    affine_fixed_count(points_over_infinity, MULTIPLIER, translation, power)
    for translation in range(points_over_infinity)
    for power in range(1, ORDER5)
}
assert zero_translations == tuple(range(0, M, 3))
assert zero_fixed_counts == {3}
assert infinity_fixed_counts == {1}
fixed_over_zero_per_power = zero_fixed_counts.pop()
fixed_over_infinity_per_power = infinity_fixed_counts.pop()
fixed_per_power = fixed_over_zero_per_power + fixed_over_infinity_per_power

two_genus_quotient_minus_two_numerator = (
    two_genus_c_minus_two - (ORDER5 - 1) * fixed_per_power
)
assert two_genus_quotient_minus_two_numerator % ORDER5 == 0
two_genus_quotient_minus_two = two_genus_quotient_minus_two_numerator // ORDER5
assert two_genus_quotient_minus_two % 2 == 0
genus_quotient = 1 + two_genus_quotient_minus_two // 2
assert fixed_per_power == 4
assert genus_quotient == 11

# The full semidirect cover C -> C/(mu_33 rtimes mu_5) has three branch
# orders.  The first two are 33 and 5; Riemann--Hurwitz determines the third.
semidirect_group_order = M * ORDER5
third_branch_reciprocal = (
    Fraction(1, 1)
    - Fraction(two_genus_c_minus_two, semidirect_group_order)
    - Fraction(1, M)
    - Fraction(1, ORDER5)
)
assert third_branch_reciprocal.numerator == 1
triangle_signature = (M, ORDER5, third_branch_reciprocal.denominator)
assert triangle_signature == (33, 5, 15)

# Exact conductor decomposition and order-five descent.  The primitive part
# V=e_33 H^1(C,Q) has Q-dimension 80 and K=Q(zeta_33)-rank 4.  Its rho-fixed
# part W has Q-dimension 16 and rank 4 over K^{<4>}, hence corresponds to an
# eight-dimensional abelian factor of J(C/<rho>).
quotient_h1_dimension = 2 * genus_quotient
h1_dimensions_by_order, invariant_dimensions_by_order, cycles_by_order = (
    h1_conductor_decomposition(quotient_h1_dimension)
)
assert h1_dimensions_by_order == {3: 8, 11: 30, 33: 80}
assert sum(h1_dimensions_by_order.values()) == 2 * genus_c
assert invariant_dimensions_by_order == {3: 0, 11: 6, 33: 16}
assert sum(invariant_dimensions_by_order.values()) == quotient_h1_dimension
conductor33_projector_support = tuple(
    character
    for character in range(M)
    if conductor33_projector_value(character) == 1
)
assert conductor33_projector_support == primitive_units(M)

# Schoen's U=exterior_K^4(V) and its order-five descent.  This is still a
# Hodge-structure computation, not an algebraicity claim.
rho_character_cycles, mu5_eigenvalue_multiplicities = u_mu5_decomposition()
u_dimension_over_c = len(primitive_units(M))
u_mu5_invariant_rank = mu5_eigenvalue_multiplicities[0]
fixed_cm_field_degree = len(rho_character_cycles)
primitive_factor_h1_dimension = invariant_dimensions_by_order[33]
primitive_abelian_factor_dimension = primitive_factor_h1_dimension // 2
primitive_rank_over_fixed_field = (
    primitive_factor_h1_dimension // fixed_cm_field_degree
)
fixed_field_hodge_signatures = tuple(
    character_hodge_split(cycle[0]) for cycle in rho_character_cycles
)
generalized_weil_class_dimension = fixed_cm_field_degree
# Each primitive induced representation occurs twice in H^{1,0}; the four
# character orbits form two dual pairs.  Thus the primitive group action leaves
# an eight-dimensional tangent space in Siegel space and does not force CM.
primitive_cycle_sets = tuple(frozenset(cycle) for cycle in rho_character_cycles)
dual_cycle_pairs = {
    frozenset(
        (
            cycle_set,
            frozenset((-character) % M for character in cycle_set),
        )
    )
    for cycle_set in primitive_cycle_sets
}
primitive_siegel_fixed_tangent_dimension = sum(2 * 2 for _ in dual_cycle_pairs)
assert u_dimension_over_c == 20
assert sum(mu5_eigenvalue_multiplicities) == u_dimension_over_c
assert u_mu5_invariant_rank == 4
assert fixed_cm_field_degree == 4
assert primitive_abelian_factor_dimension == 8
assert primitive_rank_over_fixed_field == 4
assert fixed_field_hodge_signatures == ((2, 2),) * fixed_cm_field_degree
assert generalized_weil_class_dimension == u_mu5_invariant_rank
assert len(dual_cycle_pairs) == 2
assert primitive_siegel_fixed_tangent_dimension == 8

# Common cyclic quotient with the Fermat fourfold.  If
# P(s,t)=sum_{k=0}^4 u_k s^(4-k)t^k, the six hyperplanes in Sym^4(P^1)
# associated with b_i=[1:zeta_5^i] and infinity have linear forms
# L_i=sum_k u_k zeta_5^(ik), L_infinity=u_4.  Hence
# sum_i zeta_5^i L_i=5 L_infinity.  After rescaling these are the six
# coordinate hyperplanes in sum Z_i=0.  On the chart L_infinity=1, the
# product of the four curve equations has infinity divisor exponent -77.
# Multiplying its Kummer generator by L_infinity^3 replaces -77 by 22.
# A coordinate permutation then puts the exponents in the standard alpha
# order.  The character alpha is surjective, so both quotient constructions
# have cyclic degree 33.  The two quotient maps also give a nonzero transfer:
# transpose(Gamma) Gamma acts on the G-invariants by |ker chi_alpha| |G|.
hyperplane_fourier_sums = tuple(
    root_of_unity_sum(coefficient_index + 1)
    for coefficient_index in range(ORDER5)
)
raw_infinity_divisor_exponent = -sum(FINITE_EXPONENTS)
infinity_principal_shift = (
    INFINITY_EXPONENT - raw_infinity_divisor_exponent
) // M
normalized_branch_exponents = FINITE_EXPONENTS + (
    raw_infinity_divisor_exponent + infinity_principal_shift * M,
)
fermat_ordered_branch_exponents = tuple(
    normalized_branch_exponents[index]
    for index in FERMAT_COORDINATE_PERMUTATION
)
# A homogeneous monomial is converted into a projective rational function by
# subtracting its total degree from one exponent.  This changes the Kummer
# divisor only by a 33-fold principal divisor.
fermat_projective_exponent_lifts = ALPHA[:-1] + (
    ALPHA[-1] - sum(ALPHA),
)
common_cyclic_cover_degree = M // gcd(M, *BRANCH_EXPONENTS)
schoen_quotient_group_order = factorial(4) * M ** 3
fermat_character_kernel_order = M ** 4
common_quotient_transfer_scalar = (
    schoen_quotient_group_order * fermat_character_kernel_order
)
assert hyperplane_fourier_sums == (0, 0, 0, 0, 5)
assert raw_infinity_divisor_exponent == -77
assert infinity_principal_shift == 3
assert normalized_branch_exponents == BRANCH_EXPONENTS
assert fermat_ordered_branch_exponents == ALPHA
assert sum(fermat_projective_exponent_lifts) == 0
assert tuple(value % M for value in fermat_projective_exponent_lifts) == ALPHA
assert sum(BRANCH_EXPONENTS) % M == 0
assert common_cyclic_cover_degree == 33
assert common_quotient_transfer_scalar > 0

print("alpha:", ALPHA)
print("unit stabilizer generator:", MULTIPLIER)
print("unit stabilizer:", STABILIZER)
print("stabilizer orbit:", permutation_cycles(ALPHA, MULTIPLIER))
print("fixed cyclotomic field:", "Q(zeta_3, sqrt(-11))")
print("finite branch exponents:", FINITE_EXPONENTS)
print("infinity branch exponent:", INFINITY_EXPONENT)
print(f"gcd(33,{INFINITY_EXPONENT}):", points_over_infinity)
print("Kummer correction exponents k_i:", K_EXPONENTS)
print("rho^5 finite-branch exponents:", branch_exponents)
print("rho^5 constant exponent in mu_165:", constant_exponent_mod_165)
print("branch exponent gcds with 33:", branch_gcds)
print("Riemann-Hurwitz contributions:", ramification_contributions)
print("points over infinity:", points_over_infinity)
print("ramification index at infinity:", infinity_ramification_index)
print("infinity affine action on Z/11 multiplier:", infinity_affine_multiplier)
print("fixed points over z=0 per nontrivial lift power:", fixed_over_zero_per_power)
print("fixed points over infinity per nontrivial lift power:", fixed_over_infinity_per_power)
print("fixed points per nontrivial lift power (total):", fixed_per_power)
print("curve genus:", genus_c)
print("quotient genus by order-5 lift:", genus_quotient)
print("triangle signature of full semidirect cover:", triangle_signature)
print("H^1 dimensions by exact character order:", tuple(sorted(h1_dimensions_by_order.items())))
print("rho-invariant H^1 dimensions by exact character order:", tuple(sorted(invariant_dimensions_by_order.items())))
print("conductor-33 projector support size:", len(conductor33_projector_support))
print("Schoen U primitive-character cycles under rho:", rho_character_cycles)
print("Schoen U dimension over C:", u_dimension_over_c)
print("H^1 primitive-character Hodge split (h10,h01):", (2, 2))
print("mu_5 eigenvalue multiplicities on U:", mu5_eigenvalue_multiplicities)
print("mu_5-invariant projector rank on U:", u_mu5_invariant_rank)
print("primitive quotient abelian-factor dimension:", primitive_abelian_factor_dimension)
print("fixed CM field degree and H^1 rank:", (fixed_cm_field_degree, primitive_rank_over_fixed_field))
print("fixed-field Hodge signatures:", fixed_field_hodge_signatures)
print("generalized Weil class dimension:", generalized_weil_class_dimension)
print("primitive Siegel fixed-locus tangent dimension:", primitive_siegel_fixed_tangent_dimension)
print("symmetric hyperplane Fourier sums:", hyperplane_fourier_sums)
print("normalized branch exponents in alpha order:", fermat_ordered_branch_exponents)
print("projective Fermat exponent lifts:", fermat_projective_exponent_lifts)
print("common Fermat cyclic-quotient degree:", common_cyclic_cover_degree)
print("common-quotient transfer scalar:", common_quotient_transfer_scalar)
print("Schoen simple-pair test:", any(a + b == M for i, a in enumerate(ALPHA) for b in ALPHA[i + 1:]))
print("status: Weil descent/common quotient exact")
