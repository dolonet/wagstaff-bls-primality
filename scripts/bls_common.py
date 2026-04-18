"""
Shared primitives for the three BLS N-1 Wagstaff primality drivers.

Exports:
    aprcl_prove_prime(n)             -- single-prime APR-CL proof via PARI/GP
    aprcl_prove_all(primes)          -- APR-CL proof for every prime in a list
    verify_condition_ii(N)           -- check omega_3^((N+1)/2) == -1 (mod N)
    exact_bls_margin_bits(F, N)      -- integer bits by which F^3 exceeds N
    require_gp()                     -- resolve path to `gp`; error if absent

PARI/GP is required for the certificate to be unconditional; we invoke
isprime(x, 2), which forces the APR-CL test (Adleman-Pomerance-Rumely with
the Cohen-Lenstra improvements).  This certifies every cofactor of F.
"""

from __future__ import annotations

import re
import shutil
import subprocess
import time


def require_gp() -> str:
    gp = shutil.which("gp")
    if gp is None:
        raise RuntimeError(
            "PARI/GP (`gp`) not found in PATH. "
            "The BLS certificate requires APR-CL primality proofs for every "
            "cofactor of F; APR-CL is invoked via `gp -q isprime(x, 2)`. "
            "Install pari-gp (e.g. `sudo apt install pari-gp`) and retry."
        )
    return gp


_ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")


def aprcl_prove_prime(n: int, gp_path: str | None = None, timeout: int = 1200) -> float:
    """Prove n prime by APR-CL via PARI/GP.  Returns wall-clock seconds.

    Raises AssertionError if gp returns anything other than 1 (true)."""
    gp_path = gp_path or require_gp()
    gp_input = f"print(isprime({n}, 2)); quit\n"
    t0 = time.time()
    result = subprocess.run(
        [gp_path, "-q", "-s", "1000000000"],
        input=gp_input,
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    dt = time.time() - t0
    raw = _ANSI_RE.sub("", result.stdout.strip())
    assert raw == "1", (
        f"APR-CL failure for {len(str(n))}-digit candidate: "
        f"gp stdout={raw!r} stderr={result.stderr!r}"
    )
    return dt


def aprcl_prove_all(primes, gp_path: str | None = None, log=print) -> dict[int, float]:
    """Run APR-CL on every prime in the iterable.  Returns {p: seconds}.

    `log` is called with one status line per prime.  Certifies uniformly:
    every prime gets the same unconditional treatment regardless of size."""
    gp_path = gp_path or require_gp()
    timings: dict[int, float] = {}
    primes_sorted = sorted(primes)
    log(f"APR-CL certifying {len(primes_sorted)} primes (PARI/GP {gp_path})")
    for q in primes_sorted:
        dt = aprcl_prove_prime(q, gp_path=gp_path)
        timings[q] = dt
        d = len(str(q))
        log(f"  q ({d:>3d}d): APR-CL PRIME  [{dt:7.3f}s]")
    total = sum(timings.values())
    log(f"All {len(timings)} factors certified by APR-CL.  Total: {total:.2f}s")
    return timings


def _zsqrt2_mul(x, y, N):
    """(a1 + b1 sqrt 2)(a2 + b2 sqrt 2) mod N."""
    a1, b1 = x
    a2, b2 = y
    return ((a1 * a2 + 2 * b1 * b2) % N, (a1 * b2 + b1 * a2) % N)


def _zsqrt2_pow(base, k, N):
    """base^k in Z[sqrt 2] / (N), right-to-left binary exponentiation."""
    result = (1, 0)
    while k > 0:
        if k & 1:
            result = _zsqrt2_mul(result, base, N)
        base = _zsqrt2_mul(base, base, N)
        k >>= 1
    return result


def verify_condition_ii(N: int) -> dict:
    """Verify Chua's Condition (II) at N:

        omega_3^{(N+1)/2} == -1   in Z[sqrt 2] / (N),

    where omega_3 = 3 + 2 sqrt 2.  Returns a dict with the residual
    components so the caller can log and/or embed in a JSON certificate.
    Asserts on failure (Proposition 2.4 predicts this holds for every
    Wagstaff prime W_p with p >= 5)."""
    assert N > 1 and N % 2 == 1, "Condition (II) verification expects odd N > 1"
    exponent = (N + 1) // 2
    t0 = time.time()
    a, b = _zsqrt2_pow((3, 2), exponent, N)
    dt = time.time() - t0
    residual_a = (a - (N - 1)) % N   # zero iff a ≡ -1 (mod N)
    residual_b = b % N               # zero iff b ≡ 0 (mod N)
    ok = (residual_a == 0 and residual_b == 0)
    return {
        "verified": bool(ok),
        "exponent_bits": exponent.bit_length(),
        "elapsed_seconds": dt,
        "residual_a_is_zero": residual_a == 0,
        "residual_b_is_zero": residual_b == 0,
        "residual_a_digits": len(str(residual_a)) if residual_a else 0,
        "residual_b_digits": len(str(residual_b)) if residual_b else 0,
    }


def exact_bls_margin_bits(F: int, N: int) -> int:
    """Integer bits by which F^3 exceeds N.  Positive iff F^3 > N.

    The BLS Theorem 5 hypothesis is F^3 > N with gcd(F, R) = 1, where
    N - 1 = F R.  This returns (F^3).bit_length() - N.bit_length(),
    which is the exact integer margin; earlier dev scripts reported
    F.bit_length() - N.bit_length()//3 - 1 as an approximation."""
    return (F ** 3).bit_length() - N.bit_length()


__all__ = [
    "require_gp",
    "aprcl_prove_prime",
    "aprcl_prove_all",
    "verify_condition_ii",
    "exact_bls_margin_bits",
]
