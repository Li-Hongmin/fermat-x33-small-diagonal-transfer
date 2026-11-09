"""Exact separator for the quartic-CM determinant projector.

The field is E = Q(sqrt(-3), sqrt(-11)); elements are stored in the basis
1, s3, s11, s3*s11.  The script checks that a0=-1-s3-s11 separates all
35 degree-four monomials in the four conjugates, and prints the pure block
polynomial used in the Chow projector.
"""

from fractions import Fraction
from itertools import product

BASIS = ((0, 0), (1, 0), (0, 1), (1, 1))
ZERO = (0, 0, 0, 0)
ONE = (1, 0, 0, 0)


def add(x, y):
    return tuple(a + b for a, b in zip(x, y))


def neg(x):
    return tuple(-a for a in x)


def mul(x, y):
    out = [0, 0, 0, 0]
    for i, a in enumerate(x):
        for j, b in enumerate(y):
            e, f = BASIS[i]
            g, h = BASIS[j]
            coefficient = a * b
            e += g
            f += h
            if e >= 2:
                e -= 2
                coefficient *= -3
            if f >= 2:
                f -= 2
                coefficient *= -11
            out[BASIS.index((e, f))] += coefficient
    return tuple(out)


def power(x, n):
    out = ONE
    for _ in range(n):
        out = mul(out, x)
    return out


def polynomial_mul(p, q):
    out = [ZERO] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] = add(out[i + j], mul(a, b))
    return out


def polynomial_from_roots(roots):
    out = [ONE]
    for root in roots:
        out = polynomial_mul(out, [neg(root), ONE])
    return out


def scalar_polynomial(coefficients):
    assert all(c[1:] == (0, 0, 0) for c in coefficients)
    return [Fraction(c[0]) for c in coefficients]


def trim(p):
    p = list(p)
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def polynomial_divmod(a, b):
    a = trim(a)
    b = trim(b)
    quotient = [Fraction(0)] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b) and any(a):
        coefficient = a[-1] / b[-1]
        shift = len(a) - len(b)
        quotient[shift] = coefficient
        for i, value in enumerate(b):
            a[i + shift] -= coefficient * value
        a = trim(a)
    return trim(quotient), trim(a)


def polynomial_gcd(a, b):
    while any(b):
        _, remainder = polynomial_divmod(a, b)
        a, b = b, remainder
    leading = a[-1]
    return [value / leading for value in a]


def polynomial_add(a, b):
    out = [Fraction(0)] * max(len(a), len(b))
    for i, value in enumerate(a):
        out[i] += value
    for i, value in enumerate(b):
        out[i] += value
    return trim(out)


def polynomial_sub(a, b):
    return polynomial_add(a, [-value for value in b])


def polynomial_scale(a, scalar):
    return trim([scalar * value for value in a])


def fraction_polynomial_mul(a, b):
    out = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, left in enumerate(a):
        for j, right in enumerate(b):
            out[i + j] += left * right
    return trim(out)


def polynomial_extended_gcd(a, b):
    old_r, r = trim(a), trim(b)
    old_s, s = [Fraction(1)], [Fraction(0)]
    old_t, t = [Fraction(0)], [Fraction(1)]
    while any(r):
        quotient, remainder = polynomial_divmod(old_r, r)
        old_r, r = r, remainder
        old_s, s = s, polynomial_sub(
            old_s, fraction_polynomial_mul(quotient, s)
        )
        old_t, t = t, polynomial_sub(
            old_t, fraction_polynomial_mul(quotient, t)
        )
    scale = Fraction(1, old_r[-1])
    return (
        polynomial_scale(old_r, scale),
        polynomial_scale(old_s, scale),
        polynomial_scale(old_t, scale),
    )


def polynomial_mod(a, modulus):
    return polynomial_divmod(a, modulus)[1]


def main():
    # a0 = -1 - sqrt(-3) - sqrt(-11)
    a0 = (-1, -1, -1, 0)
    conjugates = [
        (a0[0], e * a0[1], f * a0[2], e * f * a0[3])
        for e, f in product((1, -1), repeat=2)
    ]

    values = []
    for k0 in range(5):
        for k1 in range(5 - k0):
            for k2 in range(5 - k0 - k1):
                exponents = (k0, k1, k2, 4 - k0 - k1 - k2)
                value = ONE
                for conjugate, exponent in zip(conjugates, exponents):
                    value = mul(value, power(conjugate, exponent))
                values.append((exponents, value))

    assert len(values) == 35
    assert len({value for _, value in values}) == 35

    pure = polynomial_from_roots(
        [value for exponents, value in values if 4 in exponents]
    )
    mixed = polynomial_from_roots(
        [value for exponents, value in values if 4 not in exponents]
    )
    pure_q = scalar_polynomial(pure)
    mixed_q = scalar_polynomial(mixed)
    assert len(pure_q) == 5
    assert len(mixed_q) == 32
    assert polynomial_gcd(pure_q, mixed_q) == [Fraction(1)]

    # The inverse of P_mix modulo P_pure is the degree-<4 factor in the
    # rational Chow projector Q_E(T)=P_mix(T)*P_mix(T)^(-1) mod P_pure(T).
    gcd_value, inverse, _ = polynomial_extended_gcd(mixed_q, pure_q)
    assert gcd_value == [Fraction(1)]
    inverse_mod_pure = polynomial_mod(inverse, pure_q)
    separator = polynomial_mod(
        fraction_polynomial_mul(mixed_q, inverse_mod_pure), pure_q
    )
    assert separator == [Fraction(1)]

    print("distinct_lambda_count =", len({value for _, value in values}))
    print("pure_polynomial_high_to_low =", list(reversed([int(x) for x in pure_q])))
    print("mixed_degree =", len(mixed_q) - 1)
    print("gcd_degree =", len(polynomial_gcd(pure_q, mixed_q)) - 1)
    print("inverse_mod_pure_low_to_high =", inverse_mod_pure)
    print("projector_remainder =", separator)


if __name__ == "__main__":
    main()
