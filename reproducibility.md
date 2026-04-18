# Reproducibility guide

This repository ships three BLS $N-1$ primality certificates for
$W_{2617}$, $W_{10501}$, and $W_{12391}$, plus the drivers that
produced them and an independent re-verifier. Everything listed below
is either computed from scratch on each run or fixed as a literal
in-source (cyclotomic factorisations). No external data files are
downloaded.

## Environment the release was built against

| Tool       | Version                | Source                         |
|------------|------------------------|--------------------------------|
| Python     | 3.12.3                 | Ubuntu 24.04 system interpreter|
| sympy      | 1.14.0                 | `pip install sympy`            |
| PARI/GP    | 2.15.4                 | `apt install pari-gp`          |
| OS kernel  | Linux 6.8.0 x86_64     | Ubuntu 24.04 LTS               |

Exact parity is not required. Any Python $\geq$ 3.10 with sympy
$\geq$ 1.12 and PARI/GP $\geq$ 2.15 should work; nothing in the drivers
depends on CPython or sympy internals.

## What each run does

Each driver `bls_n_minus_1_wXXXX.py`:

1. Defines $N = W_p = (2^p + 1) / 3$ and asserts $N \equiv 3 \pmod 8$.
2. Verifies Chua's Condition (II), $\omega_3^{(N+1)/2} \equiv -1$ in
   $\mathbb{Z}[\sqrt 2] / (N)$, with $\omega_3 = 3 + 2\sqrt 2$. This is
   an independent sanity check on the modular arithmetic; it is not
   used as input to the BLS proof.
3. Produces the cyclotomic decomposition of $2^{p-1} - 1$ and factors
   the small $\Phi_d(2)$ values; large factorisations (the
   Cunningham-project ones and the bespoke factorisations used for
   $W_{10501}$ and $W_{12391}$) are supplied as literals in the
   driver source.
4. APR-CL certifies every prime $q$ of the factored part $F$ of
   $N - 1$ via PARI/GP `isprime(x, 2)` — APR-CL only, no BPSW
   fallback, no size threshold.
5. Checks the BLS Theorem 5 hypothesis $F^3 > N$ with exact integer
   comparison and reports the margin in bits.
6. For each prime $q \mid F$, finds a BLS witness $a$ satisfying
   $a^{N-1} \equiv 1 \pmod N$ and $\gcd(a^{(N-1)/q} - 1, N) = 1$.
7. Runs the finite-divisor discriminant argument (Theorem 5 refinement)
   and concludes PRIME / INCONCLUSIVE.
8. Writes a JSON certificate to `../data/bls_certificate_wXXXX.json`
   (falls back to the script directory if `data/` is absent).

## Expected wall-times

Approximate single-core timings; APR-CL uses only one core. Parallel
runs of distinct drivers on a multi-core machine are safe.

| Driver                      | Machine                         | Wall time |
|-----------------------------|---------------------------------|-----------|
| `bls_n_minus_1_w2617.py`    | Ryzen 5 3600 (Hetzner)          | ~ 10 s    |
| `bls_n_minus_1_w10501.py`   | 256-core cluster (single core)  | minutes   |
| `bls_n_minus_1_w12391.py`   | 256-core cluster (single core)  | minutes   |

The small-$p$ driver is fast enough to run on any workstation. The
two larger ones are memory-bound (up to a few hundred MB peak) due to
`factorint` on $\sim 200{-}400$-digit $\Phi_d(2)$ values.

## Re-verifying an existing certificate

```bash
cd scripts
python3 verify_bls.py ../data/bls_certificate_w2617.json
```

`verify_bls.py` reads the JSON and treats every claim inside it as
untrusted input. It:

* rebuilds $N$ from the claimed formula and checks `digits` and `bits`;
* APR-CL-certifies every prime in `witnesses` (pass `--skip-aprcl` to
  skip, for development use only);
* rebuilds $F$ as $\prod q^{v_q(N-1)}$ from the prime list directly —
  so a certificate that substitutes bogus factors would fail to
  reconstruct the claimed $F$;
* verifies $\gcd(F, R) = 1$ and $F^3 > N$, and re-checks the exact
  margin in bits;
* re-runs the BLS witness check for each $q$;
* reruns the finite-divisor discriminant computation and compares to
  the certificate's `discriminant_sign`;
* re-runs Condition (II) in $\mathbb{Z}[\sqrt 2] / (N)$.

Any disagreement exits non-zero. A successful run prints
`OK  certificate … re-verified from scratch`.

## Certificate schema

Each certificate is a UTF-8 JSON document with the following top-level
fields:

```
number                    string, e.g. "W_2617"
formula                   string, e.g. "(2^2617 + 1) / 3"
digits                    int, len(str(N))
bits                      int, N.bit_length()
result                    "PRIME" | "INCONCLUSIVE"
method                    "BLS Theorem 5 (N-1)"
reference                 BLS75 citation
factored_part             { bits, digits, num_primes }
unfactored_part           { bits, digits }
bls_hypothesis            {
  statement               "F^3 > N",
  satisfied               bool,
  exact_margin_bits       int = (F^3).bit_length() - N.bit_length(),
  definition              string
}
discriminant_sign         "negative" | "non-square" | "square-..."
cyclotomic_decomposition  { base, factorization_of_exponent,
                            num_divisors, num_factored,
                            cunningham_factors_used }
witnesses                 { "q_str": a_int, ... }
aprcl_certification       { tool, scope, num_certified,
                            total_elapsed_seconds,
                            max_elapsed_seconds }
chebyshev_condition_ii    { base, congruence, verified, exponent_bits,
                            elapsed_seconds }
software                  { python, sympy, pari_gp, note }
```

## What is *not* shipped

* Factor tables in bulk. Only the specific primes needed to reach
  $F^3 > N$ for each exponent are baked into the driver source.
* The full Cunningham-project archive. A driver may reference a
  particular Cunningham result as a literal; the full reference list
  is in the manuscript.
* PDFs of the paper — this repository's `paper/` directory contains
  only the TeX and Markdown sources. Build outputs belong in a
  separate archive.

## Reporting a divergence

If a run produces a different JSON than the one in `data/`, or if
`verify_bls.py` reports `FAIL:` against the shipped certificate,
please open an issue including the full stdout, the machine you
ran on, and the versions of Python, sympy, and PARI/GP.
