# A small-diagonal transfer from a degree-33 Fermat character to the quartic-CM determinant of a Prym eightfold, and the Hodge conjecture for its powers under explicit quotient hypotheses

[![verify](https://github.com/Li-Hongmin/fermat-x33-small-diagonal-transfer/actions/workflows/verify.yml/badge.svg)](https://github.com/Li-Hongmin/fermat-x33-small-diagonal-transfer/actions/workflows/verify.yml) [![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.24054639-blue)](https://doi.org/10.5281/zenodo.24054639)

**Hongmin Li** (Institute of Science Tokyo; The University of Tokyo) · ORCID [0000-0003-0228-0600](https://orcid.org/0000-0003-0228-0600)

Main result of the manuscript: for the degree-33 Fermat character α = (7,10,13,19,22,28), the Prym eightfold B₈ has End⁰ = Q(ζ₃₃)^⟨4⟩ and special-unitary Hodge group (by exact Frobenius computations at two primes), and a small-diagonal transfer gives the Hodge conjecture for its powers under explicit quotient hypotheses.

- Manuscript on rStation: <https://rstation.ai/problems/hodge_conjecture/papers/20/v5> (version 5)
- Archived on Zenodo: concept DOI [10.5281/zenodo.24054639](https://doi.org/10.5281/zenodo.24054639) (all versions); this repository mirrors version 5, DOI [10.5281/zenodo.24054685](https://doi.org/10.5281/zenodo.24054685)

## Contents

- `paper/paper.md`, `paper/fermat-x33-small-diagonal-transfer.pdf` — the manuscript (version 5).
- `verification/` — the supplement of that version, unpacked unchanged; its own README describes each file.
- `.github/workflows/verify.yml` — reruns every check below on each push.

## What the programs check

| Check | Command (run inside `verification/`) |
|---|---|
| Transfer arithmetic | `python calculations/bundles/fermat-x33-transfer-certificate/verify_transfer_arithmetic.py` |
| Abel-Fourier sign | `python calculations/bundles/fermat-x33-transfer-certificate/verify_abel_fourier_sign.py` |
| Abel-Fourier exterior normalization | `python calculations/bundles/fermat-x33-transfer-certificate/verify_abel_fourier_exterior_normalization.py` |
| Common quotient definition field | `python calculations/bundles/fermat-x33-transfer-certificate/verify_common_quotient_definition_field.py` |
| Common quotient coordinates | `python calculations/bundles/fermat-x33-transfer-certificate/verify_common_quotient_coordinates.py` |
| Klein cubic character map | `python calculations/bundles/fermat-x33-transfer-certificate/verify_klein_cubic_character_map.py` |
| Source projector representation | `python calculations/bundles/fermat-x33-transfer-certificate/verify_source_projector_representation.py` |
| Q* character block | `python calculations/bundles/fermat-x33-transfer-certificate/verify_qstar_character_block.py` |
| Order-5 cyclic cover | `python calculations/bundles/fermat-x33-order5-cyclic-cover/x33-order5-cyclic-cover.py` |
| Explicit separator for the quartic-CM projector | `python calculations/bundles/x33-quartic-cm-projector/verify_explicit_separator.py` |
| Level-66 Aoki bridge | `python calculations/bundles/fermat-x33-level66-aoki/verify_level66_aoki_bridge.py` |
| PARI/GP: irreducibility, root ratios, subfields, discriminants | `gp -q calculations/bundles/fermat-x33-frobenius-cm-screen/verify-b8-end-field.gp < /dev/null` |
| SageMath: Frobenius polynomials at q = 331, 661 | `sage -python calculations/bundles/fermat-x33-frobenius-cm-screen/x33-b-frobenius-cm-screen.py` |

The programs check the finite and exact computations used in the manuscript; the mathematical arguments are those written in the paper.

## Reproducing

```sh
cd verification
pip install sympy==1.14.0
python calculations/bundles/fermat-x33-transfer-certificate/verify_transfer_arithmetic.py
python calculations/bundles/fermat-x33-transfer-certificate/verify_abel_fourier_sign.py
python calculations/bundles/fermat-x33-transfer-certificate/verify_abel_fourier_exterior_normalization.py
python calculations/bundles/fermat-x33-transfer-certificate/verify_common_quotient_definition_field.py
python calculations/bundles/fermat-x33-transfer-certificate/verify_common_quotient_coordinates.py
python calculations/bundles/fermat-x33-transfer-certificate/verify_klein_cubic_character_map.py
python calculations/bundles/fermat-x33-transfer-certificate/verify_source_projector_representation.py
python calculations/bundles/fermat-x33-transfer-certificate/verify_qstar_character_block.py
python calculations/bundles/fermat-x33-order5-cyclic-cover/x33-order5-cyclic-cover.py
python calculations/bundles/x33-quartic-cm-projector/verify_explicit_separator.py
python calculations/bundles/fermat-x33-level66-aoki/verify_level66_aoki_bridge.py
gp -q calculations/bundles/fermat-x33-frobenius-cm-screen/verify-b8-end-field.gp < /dev/null
sage -python calculations/bundles/fermat-x33-frobenius-cm-screen/x33-b-frobenius-cm-screen.py
```

Requirements: Python 3.9+, PARI/GP, SageMath. Run the commands from a copy if you want to keep the recorded outputs untouched.

## Cite

```bibtex
@misc{li2026fermatx33smalldiagonaltr,
  author    = {Li, Hongmin},
  title     = {A small-diagonal transfer from a degree-33 Fermat character to the quartic-CM determinant of a Prym eightfold, and the Hodge conjecture for its powers under explicit quotient hypotheses},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.24054639},
  url       = {https://doi.org/10.5281/zenodo.24054639}
}
```

## License

Manuscript and data: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Code: MIT (see `LICENSE`).
