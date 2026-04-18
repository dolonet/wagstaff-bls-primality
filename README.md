# wagstaff-bls-primality

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19643792.svg)](https://doi.org/10.5281/zenodo.19643792)

Unconditional primality proofs of the Wagstaff numbers
$W_{2617}$, $W_{10501}$, and $W_{12391}$ via the Brillhart–Lehmer–Selfridge
$N-1$ criterion, independent of ECPP.

- **Method.** BLS Theorem 5 \[BLS75\] applied to the cyclotomic
  decomposition $2^{p-1} - 1 = \prod_{d \mid p-1} \Phi_d(2)$, with
  factorisations drawn from the Cunningham project tables and direct
  ECM / Pollard-$\rho$ / $p{-}1$ computation. Every prime of the
  factored part $F$ of $N - 1$ is certified by APR-CL \[APR83, CL84\].
- **Independent sanity.** Each certificate also records the Chua
  congruence $\omega_3^{(W_p + 1)/2} \equiv -1 \pmod{W_p}$, verified in
  $\mathbb{Z}[\sqrt 2] / (W_p)$.
- **No new conjectures.** The proofs do not assume the Chebyshev
  sufficiency conjecture; that conjecture is the subject of a separate
  paper.

## Layout

```
paper/      manuscript sources (Markdown draft + TeX build)
scripts/    primality drivers and the independent re-verifier
data/       JSON certificates (one per W_p)
logs/       (optional) full stdout of each driver run
LICENSE     MIT license for scripts and data
LICENSE-PAPER  CC-BY-4.0 license for paper/ and data/*.json
CITATION.cff   citation metadata
reproducibility.md   how to rebuild every certificate from scratch
```

## Reproducing a certificate

Requirements: Python 3.10+, sympy 1.12+, PARI/GP 2.15+ on `PATH`.

```bash
cd scripts
python3 bls_n_minus_1_w2617.py     # ~10 s, writes ../data/bls_certificate_w2617.json
python3 bls_n_minus_1_w10501.py    # minutes, factored part cached as literals
python3 bls_n_minus_1_w12391.py    # minutes, factored part cached as literals
```

Each driver writes its certificate to `../data/bls_certificate_wXXXX.json`
(falls back to the script directory if `data/` is not present).

## Re-verifying a certificate

`verify_bls.py` treats the certificate as a claim under test. It
re-reads the JSON, rebuilds $F$ from $v_q(N - 1)$ for each listed prime
$q$ (so a forged certificate with bogus factors is caught), re-runs
APR-CL on every $q$, re-checks the BLS witness conditions and the
finite-divisor discriminant argument, and re-runs Condition (II) in
$\mathbb{Z}[\sqrt 2] / (N)$.

```bash
cd scripts
python3 verify_bls.py ../data/bls_certificate_w2617.json
python3 verify_bls.py ../data/bls_certificate_w10501.json
python3 verify_bls.py ../data/bls_certificate_w12391.json
```

A successful run prints `OK  certificate … re-verified from scratch`
and exits 0. Any disagreement exits non-zero.

## Citation

See `CITATION.cff`, or cite as:

> Dolotov, A. *Three Brillhart–Lehmer–Selfridge primality proofs for
> Wagstaff numbers.* 2026. Zenodo. <https://doi.org/10.5281/zenodo.19643792>

## References

- \[APR83\] L. M. Adleman, C. Pomerance, R. S. Rumely. *On distinguishing prime numbers from composite numbers.* Ann. of Math. (2) 117 (1983), 173–206.
- \[BLS75\] J. Brillhart, D. H. Lehmer, J. L. Selfridge. *New primality criteria and factorizations of $2^m \pm 1$.* Math. Comp. 29 (1975), 620–647.
- \[CL84\] H. Cohen, H. W. Lenstra Jr. *Primality testing and Jacobi sums.* Math. Comp. 42 (1984), 297–330.
- \[Chu\] K. S. Chua. *A Chebyshev-like primality criterion for Wagstaff and Mersenne numbers.* Preprint.

## License

- Source code (`scripts/`) and machine-readable certificates
  (`data/*.json`): MIT — see `LICENSE`.
- Manuscript and its build outputs (`paper/`): CC-BY-4.0 — see
  `LICENSE-PAPER`.
