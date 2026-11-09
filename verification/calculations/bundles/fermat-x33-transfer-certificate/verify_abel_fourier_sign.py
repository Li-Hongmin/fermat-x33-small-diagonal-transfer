"""Exterior-algebra certificate for the Abel--Fourier sign.

With

    theta = sum_i a_i ^ b_i
    c1(P) = sum_i (a_i ^ b'_i - b_i ^ a'_i),

the H^1-output of

    exp(c1(P)) theta^(g-1)/(g-1)! ^ a_0

is exactly -a'_0 for every tested genus 1 <= g <= 11.  Only the linear
term in exp(c1(P)) can contribute to top degree on the first factor.
"""

from fractions import Fraction
from math import factorial


def wedge_mask(left, right):
    if left & right:
        return None, 0
    inversions = 0
    remaining = left
    while remaining:
        bit = remaining & -remaining
        index = bit.bit_length() - 1
        remaining -= bit
        inversions += (right & ((1 << index) - 1)).bit_count()
    return left | right, -1 if inversions % 2 else 1


def add(left, right):
    result = dict(left)
    for mask, coefficient in right.items():
        result[mask] = result.get(mask, Fraction(0)) + coefficient
    return {mask: value for mask, value in result.items() if value}


def multiply(left, right):
    result = {}
    for left_mask, left_coefficient in left.items():
        for right_mask, right_coefficient in right.items():
            mask, sign = wedge_mask(left_mask, right_mask)
            if sign:
                result[mask] = result.get(mask, Fraction(0)) + sign * (
                    left_coefficient * right_coefficient
                )
    return {mask: value for mask, value in result.items() if value}


def power(element, exponent):
    result = {0: Fraction(1)}
    for _ in range(exponent):
        result = multiply(result, element)
    return result


def fourier_output(g):
    theta = {}
    poincare = {}
    for index in range(g):
        left_a = 1 << (2 * index)
        left_b = 1 << (2 * index + 1)
        theta[left_a | left_b] = Fraction(1)

        right_a = 1 << (2 * g + 2 * index)
        right_b = 1 << (2 * g + 2 * index + 1)
        poincare[left_a | right_b] = Fraction(1)
        poincare[left_b | right_a] = Fraction(-1)

    theta_power = power(theta, g - 1)
    theta_power = {
        mask: coefficient / factorial(g - 1)
        for mask, coefficient in theta_power.items()
    }
    alpha = {1: Fraction(1)}  # a_0 on the first factor.
    product = multiply(multiply(theta_power, poincare), alpha)
    left_top = (1 << (2 * g)) - 1
    return {
        mask >> (2 * g): coefficient
        for mask, coefficient in product.items()
        if mask & left_top == left_top
    }


def main():
    outputs = {genus: fourier_output(genus) for genus in range(1, 12)}
    expected = {1: Fraction(-1)}
    assert all(output == expected for output in outputs.values())
    print("genus_range = 1..11")
    print("H1_output =", expected)
    print("c_J = -1")
    print("c_J^4 =", 1)
    print("all_checks = PASS")


if __name__ == "__main__":
    main()
