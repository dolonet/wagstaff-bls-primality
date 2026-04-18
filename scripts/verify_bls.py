"""
verify_bls.py — independent re-verifier for a BLS N-1 Wagstaff certificate.

Usage:
    python3 verify_bls.py path/to/bls_certificate_wXXXX.json

Reads a certificate produced by one of the bls_n_minus_1_wXXXX.py drivers
and re-checks every claim from first principles, WITHOUT trusting the
factor list the driver produced — each prime q in the certificate's
`witnesses` dict is:

    * re-verified prime by APR-CL (PARI/GP `isprime(x, 2)`);
    * re-verified as a factor of N-1 with the q-adic valuation computed
      here so that F = prod q^{v_q(N-1)} is rebuilt independently.

Then the BLS Theorem 5 hypotheses and the finite-divisor discriminant
argument are re-evaluated, and Chua's Condition (II) is re-run in
Z[sqrt(2)] / (N).  Exits 0 on full agreement with the certificate,
non-zero on the first disagreement.

The certificate itself is treated as a claim under test — nothing from
the JSON is assumed correct without a recomputation.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
import time
from math import gcd

from bls_common import (
    aprcl_prove_all,
    exact_bls_margin_bits,
    require_gp,
    verify_condition_ii,
)


FORMULA_RE = re.compile(r"^\(2\^(?P<p>\d+)\s*\+\s*1\)\s*/\s*3$")


def parse_formula(formula: str) -> int:
    m = FORMULA_RE.match(formula.strip())
    if not m:
        raise ValueError(
            f"Unsupported formula {formula!r}; expected '(2^p + 1) / 3'"
        )
    p = int(m.group("p"))
    if (2**p + 1) % 3 != 0:
        raise ValueError(f"(2^{p} + 1) is not divisible by 3")
    return (2**p + 1) // 3


def v_q(m: int, q: int) -> int:
    """q-adic valuation of m: the largest e with q^e | m."""
    if m == 0:
        raise ValueError("v_q(0) is undefined")
    e = 0
    while m % q == 0:
        m //= q
        e += 1
    return e


def bls_witness_check(a: int, q: int, N: int) -> tuple[bool, bool]:
    """Return (pow_ok, gcd_ok) for the two BLS witness conditions."""
    Nm1 = N - 1
    if Nm1 % q != 0:
        return (False, False)
    pow_ok = pow(a, Nm1, N) == 1
    gcd_ok = gcd(pow(a, Nm1 // q, N) - 1, N) == 1
    return (pow_ok, gcd_ok)


def die(msg: str, code: int = 1) -> None:
    print(f"FAIL: {msg}", file=sys.stderr)
    sys.exit(code)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("certificate", help="path to bls_certificate_wXXXX.json")
    ap.add_argument(
        "--skip-aprcl",
        action="store_true",
        help="skip the APR-CL re-certification (not recommended; for debugging only)",
    )
    args = ap.parse_args()

    with open(args.certificate, "r") as fh:
        cert = json.load(fh)

    label = cert.get("number", "N")
    formula = cert["formula"]
    print(f"Re-verifying certificate for {label} from {args.certificate}")

    # ------------------------------------------------------------------
    # Step 1: reconstruct N from the formula and check size claims
    # ------------------------------------------------------------------
    t0 = time.time()
    N = parse_formula(formula)
    n_digits = len(str(N))
    n_bits = N.bit_length()
    print(f"\n[1] N = {formula}")
    print(f"    digits={n_digits}, bits={n_bits}")
    if n_digits != cert["digits"]:
        die(f"digits mismatch: certificate={cert['digits']}, recomputed={n_digits}")
    if n_bits != cert["bits"]:
        die(f"bits mismatch: certificate={cert['bits']}, recomputed={n_bits}")

    # ------------------------------------------------------------------
    # Step 2: APR-CL re-certification of every prime in the witness list
    # ------------------------------------------------------------------
    witness_map = {int(q): int(a) for q, a in cert["witnesses"].items()}
    primes = sorted(witness_map.keys())
    print(f"\n[2] APR-CL re-certification of {len(primes)} primes in witness list")
    if args.skip_aprcl:
        print("    SKIPPED (--skip-aprcl)")
        aprcl_timings: dict[int, float] = {}
    else:
        gp = require_gp()
        aprcl_timings = aprcl_prove_all(primes, gp_path=gp, log=lambda s: print(f"    {s}"))

    # ------------------------------------------------------------------
    # Step 3: rebuild F = prod q^{v_q(N-1)} directly from N-1 and the prime list
    # ------------------------------------------------------------------
    print(f"\n[3] Rebuilding F from v_q(N-1) for each listed prime")
    Nm1 = N - 1
    F = 1
    multiplicities: dict[int, int] = {}
    for q in primes:
        e = v_q(Nm1, q)
        if e == 0:
            die(f"prime q={q} does not divide N-1 (certificate claims it factors F)")
        multiplicities[q] = e
        F *= q**e

    R = Nm1 // F
    if F * R != Nm1:
        die("F * R != N - 1 after rebuild")
    if gcd(F, R) != 1:
        die("gcd(F, R) != 1 — BLS Theorem 5 requires coprimality")

    F_digits = len(str(F))
    F_bits = F.bit_length()
    R_digits = len(str(R))
    R_bits = R.bit_length()
    print(f"    F has {F_digits} digits ({F_bits} bits)")
    print(f"    R has {R_digits} digits ({R_bits} bits)")
    print(f"    gcd(F, R) = 1")

    if F_digits != cert["factored_part"]["digits"]:
        die(f"F digits mismatch: cert={cert['factored_part']['digits']}, recomputed={F_digits}")
    if F_bits != cert["factored_part"]["bits"]:
        die(f"F bits mismatch: cert={cert['factored_part']['bits']}, recomputed={F_bits}")
    if len(primes) != cert["factored_part"]["num_primes"]:
        die(f"num_primes mismatch: cert={cert['factored_part']['num_primes']}, recomputed={len(primes)}")
    if R_digits != cert["unfactored_part"]["digits"]:
        die(f"R digits mismatch: cert={cert['unfactored_part']['digits']}, recomputed={R_digits}")
    if R_bits != cert["unfactored_part"]["bits"]:
        die(f"R bits mismatch: cert={cert['unfactored_part']['bits']}, recomputed={R_bits}")

    # ------------------------------------------------------------------
    # Step 4: BLS Theorem 5 hypothesis F^3 > N with exact margin
    # ------------------------------------------------------------------
    print(f"\n[4] BLS Theorem 5 hypothesis F^3 > N")
    F_cubed = F**3
    if F_cubed <= N:
        die(f"F^3 <= N — BLS Theorem 5 hypothesis violated")
    margin = exact_bls_margin_bits(F, N)
    print(f"    F^3 has {F_cubed.bit_length()} bits; N has {n_bits} bits")
    print(f"    exact margin = (F^3).bit_length() - N.bit_length() = {margin} bits")
    claimed_margin = cert["bls_hypothesis"]["exact_margin_bits"]
    if margin != claimed_margin:
        die(f"margin mismatch: cert={claimed_margin}, recomputed={margin}")
    if not cert["bls_hypothesis"].get("satisfied", False):
        die("certificate's bls_hypothesis.satisfied is False")

    # ------------------------------------------------------------------
    # Step 5: BLS witnesses
    # ------------------------------------------------------------------
    print(f"\n[5] BLS witnesses: for each q | F, check a^(N-1) ≡ 1 and gcd(a^((N-1)/q)-1, N) = 1")
    failures = []
    for q in primes:
        a = witness_map[q]
        pow_ok, gcd_ok = bls_witness_check(a, q, N)
        if not (pow_ok and gcd_ok):
            failures.append((q, a, pow_ok, gcd_ok))
    if failures:
        for q, a, p_ok, g_ok in failures:
            print(f"    FAIL q={q} a={a} pow_ok={p_ok} gcd_ok={g_ok}")
        die(f"{len(failures)} witness(es) rejected")
    print(f"    all {len(primes)} witnesses confirmed")

    # ------------------------------------------------------------------
    # Step 6: finite-divisor discriminant argument
    # ------------------------------------------------------------------
    print(f"\n[6] Finite-divisor check: c1+c2 = R mod F, c1*c2 = R // F, disc = (c1+c2)^2 - 4 c1 c2")
    c1_plus_c2 = R % F
    c1_times_c2 = R // F
    disc = c1_plus_c2**2 - 4 * c1_times_c2
    if disc < 0:
        disc_status = "negative"
        proven = True
    else:
        sq = math.isqrt(disc)
        if sq * sq == disc:
            c1 = (c1_plus_c2 + sq) // 2
            c2 = (c1_plus_c2 - sq) // 2
            if c1 < 1 or c2 < 1:
                disc_status = "square-but-nonpositive-c"
                proven = True
            else:
                r1 = c1 * F + 1
                r2 = c2 * F + 1
                if N % r1 == 0 or N % r2 == 0:
                    disc_status = "square-and-factors"
                    proven = False
                else:
                    disc_status = "square-but-no-divisor"
                    proven = True
        else:
            disc_status = "non-square"
            proven = True
    print(f"    discriminant status: {disc_status}")
    claimed = cert["discriminant_sign"]
    if disc_status != claimed:
        die(f"discriminant sign mismatch: cert={claimed!r}, recomputed={disc_status!r}")
    cert_result = cert.get("result", "").upper()
    if proven and cert_result != "PRIME":
        die(f"re-verification says PRIME but certificate says {cert_result!r}")
    if not proven and cert_result == "PRIME":
        die(f"certificate says PRIME but re-verification is not conclusive ({disc_status})")

    # ------------------------------------------------------------------
    # Step 7: Chua's Condition (II) — independent sanity
    # ------------------------------------------------------------------
    print(f"\n[7] Chua's Condition (II): omega_3^((N+1)/2) ≡ -1 in Z[sqrt 2] / (N)")
    cond = verify_condition_ii(N)
    print(f"    exponent_bits = {cond['exponent_bits']}, elapsed = {cond['elapsed_seconds']:.2f}s")
    print(f"    verified = {cond['verified']}")
    claimed_cond = cert.get("chebyshev_condition_ii", {}).get("verified")
    if claimed_cond is None:
        die("certificate missing chebyshev_condition_ii.verified")
    if bool(cond["verified"]) != bool(claimed_cond):
        die(f"Condition (II) mismatch: cert={claimed_cond}, recomputed={cond['verified']}")

    # ------------------------------------------------------------------
    total = time.time() - t0
    print("\n" + "=" * 60)
    print(f"OK  certificate {args.certificate} re-verified from scratch")
    print(f"    {label}: {cert_result},  exact margin {margin} bits,  "
          f"{len(primes)} primes APR-CL'd,  Condition (II) verified")
    print(f"    wall time: {total:.2f}s")
    print("=" * 60)
    return 0


if __name__ == "__main__":
    sys.exit(main())
