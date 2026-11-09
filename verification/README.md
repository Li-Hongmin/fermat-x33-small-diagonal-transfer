# Supplement

Exact computations used in the paper. Run from this directory.

- Python >= 3.10 with sympy: all `.py` files except the one below.
- SageMath (tested with 10.9): `sage -python calculations/bundles/fermat-x33-frobenius-cm-screen/x33-b-frobenius-cm-screen.py` (Frobenius polynomials at q = 331, 661 with zeta_5 -> 64 mod 331 and zeta_5 -> 247 mod 661).
- PARI/GP (tested with 2.17.3): `gp -q calculations/bundles/fermat-x33-frobenius-cm-screen/verify-b8-end-field.gp` (irreducibility, root ratios, subfields and discriminants).

Each program prints the exact quantities it computes and stops with an error if one of the stated identities fails.

The Sage program fixes the stated fifth roots of unity and checks that it reproduces the polynomials displayed in the paper. The outputs of the two programs are kept in `output-sage.txt` and `output-gp.txt` in the same folder.
