"""Finite representation certificate for the source projector P_U.

The calculation isolates the only sign-sensitive step: all four factors are
odd, so the geometric permutation action already carries the Koszul sign.
The un-signed graph average is therefore the alternating projector and has
rank one on Lambda^4(Q^4); adding an ordinary sign would incorrectly produce
the symmetric fourth power.
"""

from itertools import permutations


def inversion_parity(permutation):
    inversions = sum(
        permutation[i] > permutation[j]
        for i in range(4)
        for j in range(i + 1, 4)
    )
    return -1 if inversions % 2 else 1


def cycle_count(permutation):
    seen = set()
    count = 0
    for start in range(4):
        if start in seen:
            continue
        count += 1
        current = start
        while current not in seen:
            seen.add(current)
            current = permutation[current]
    return count


def main():
    group = list(permutations(range(4)))
    # Trace of the geometric action on (Q^4)^(tensor 4): every fixed tensor
    # assigns one basis vector per permutation cycle.  The Koszul sign is
    # sign(sigma), since all four inputs have odd degree.
    trace_geometric = sum(
        inversion_parity(sigma) * 4 ** cycle_count(sigma)
        for sigma in group
    )
    trace_wrong_sign = sum(4 ** cycle_count(sigma) for sigma in group)
    rank_geometric_average = trace_geometric // len(group)
    rank_wrong_average = trace_wrong_sign // len(group)

    assert trace_geometric == 24
    assert trace_wrong_sign == 840
    assert rank_geometric_average == 1
    assert rank_wrong_average == 35

    # The annihilator of N={t:sum t_i=0} in (Z/33)^4 is the diagonal:
    # testing t=e_i-e_j forces n_i=n_j for every i,j.
    for labels in ((1, 1, 1, 1), (1, 4, 1, 1), (1, 4, 16, 1)):
        invariant = all(
            sum(t * n for t, n in zip(test, labels)) % 33 == 0
            for test in (
                (1, 32, 0, 0),
                (1, 0, 32, 0),
                (1, 0, 0, 32),
            )
        )
        assert invariant is (len(set(labels)) == 1)

    print("geometric_graph_average_rank =", rank_geometric_average)
    print("ordinary_sign_average_rank =", rank_wrong_average)
    print("N_annihilator = diagonal")
    print("all_checks = PASS")


if __name__ == "__main__":
    main()
