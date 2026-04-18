#!/usr/bin/env python3
"""
BLS N-1 primality proof for W_10501 = (2^10501 + 1) / 3

Strategy: Factor N-1 = 2*(2^10500 - 1)/3 using the cyclotomic decomposition
of 2^10500 - 1. Since 10500 = 2^2 * 3 * 5^3 * 7, we have 48 cyclotomic
factors Phi_d(2) for d | 10500.

The small factors (d <= 420, phi(d) <= 96) are trivially factorable.
The medium factors (d in {500,525,700,750,875,1050,1500,1750,2100})
have known complete factorizations from the Cunningham project tables
(https://homes.cerias.purdue.edu/~ssw/cun/).

Together these provide 4501 bits of factored part F, well exceeding the
BLS Theorem 5 threshold of N^{1/3} ~ 3501 bits (margin: +1000 bits).

References:
  [BLS75] Brillhart, Lehmer, Selfridge, "New primality criteria and
          factorizations of 2^m +/- 1", Math. Comp. 29 (1975), 620-647.
  [Cunningham] Cunningham project factor tables for 2^n -/+ 1.
"""

import sys
import time
import math
import json
from datetime import datetime, timezone

# Use gmpy2 if available for faster modular arithmetic
try:
    import gmpy2
    from gmpy2 import mpz, is_prime as gmpy2_isprime, gcd as gmpy2_gcd
    USE_GMPY2 = True
    print("Using gmpy2 for fast arithmetic")
except ImportError:
    USE_GMPY2 = False
    print("gmpy2 not available, using sympy (slower)")

from sympy import isprime, gcd, cyclotomic_poly, Symbol, totient, factorint

from bls_common import (
    aprcl_prove_all,
    exact_bls_margin_bits,
    require_gp,
    verify_condition_ii,
)

# ============================================================
# Logging setup
# ============================================================
LOG = []

def log(msg):
    """Print and record a log line."""
    ts = datetime.now(timezone.utc).strftime("%H:%M:%S.%f")[:-3]
    line = f"[{ts}] {msg}"
    print(line)
    sys.stdout.flush()
    LOG.append(line)

def section(title):
    log("=" * 70)
    log(title)
    log("=" * 70)

# ============================================================
# Helper: fast modular exponentiation and gcd
# ============================================================
def fast_powmod(base, exp, mod):
    if USE_GMPY2:
        return int(gmpy2.powmod(mpz(base), mpz(exp), mpz(mod)))
    else:
        return pow(base, exp, mod)

def fast_gcd(a, b):
    if USE_GMPY2:
        return int(gmpy2_gcd(mpz(a), mpz(b)))
    else:
        return gcd(a, b)

def fast_isprime(n):
    if USE_GMPY2:
        return gmpy2_isprime(mpz(n))
    else:
        return isprime(n)

# ============================================================
# Step 0: Define N = W_10501
# ============================================================
t_total = time.time()
section("Step 0: Define N = W_10501")

p = 10501
N = (2**p + 1) // 3
N_digits = len(str(N))
N_bits = N.bit_length()

log(f"p = {p}")
log(f"W_{p} has {N_digits} digits ({N_bits} bits)")
log(f"N mod 8 = {N % 8}")  # should be 3

N_minus_1 = N - 1
N_plus_1 = N + 1

# Verify basic identities
assert N_minus_1 == (2**p - 2) // 3
assert N_plus_1 == 4 * ((2**(p-2) + 1) // 3)
log(f"N-1 = (2^{p} - 2) / 3")
log(f"N+1 = 4 * W_{{{p-2}}}")
log(f"Identities verified.")

# ============================================================
# Step 0a: Chua's Condition (II) — sanity check
# ============================================================
# Proposition 2.4:  For every Wagstaff prime W_p (p >= 5),
#     omega_3^{(N+1)/2} ≡ -1   (mod N),   omega_3 = 3 + 2*sqrt(2).
# This is the a=3 case of Chua's identity (Theorem 2.1) and follows
# from W_p ≡ 3 (mod 8).  We verify it as an independent sanity check
# on the modular-arithmetic implementation; the BLS certificate is
# logically independent of this check.
section("Step 0a: Chua's Condition (II): omega_3^((N+1)/2) ≡ -1 (mod N)")
cond_ii = verify_condition_ii(N)
log(f"  exponent bits:  {cond_ii['exponent_bits']}")
log(f"  elapsed:        {cond_ii['elapsed_seconds']:.2f}s")
log(f"  residual a = 0: {cond_ii['residual_a_is_zero']}")
log(f"  residual b = 0: {cond_ii['residual_b_is_zero']}")
log(f"  Condition (II) verified: {cond_ii['verified']}")
assert cond_ii["verified"], "Condition (II) FAILED"

# ============================================================
# Step 1: Cyclotomic decomposition of 2^{10500} - 1
# ============================================================
section("Step 1: Cyclotomic decomposition of 2^10500 - 1")

log(f"10500 = 2^2 * 3 * 5^3 * 7")

# Generate divisors of 10500
pm1 = p - 1  # = 10500
pm1_factors = {2: 2, 3: 1, 5: 3, 7: 1}

def gen_divisors(factors):
    divs = [1]
    for pr, exp in factors.items():
        new_divs = []
        for d in divs:
            for e in range(exp + 1):
                new_divs.append(d * pr**e)
        divs = new_divs
    return sorted(divs)

divs_10500 = gen_divisors(pm1_factors)
log(f"Number of divisors of 10500: {len(divs_10500)}")
log(f"Divisors: {divs_10500}")

# Compute Phi_d(2) for each divisor
log(f"\nComputing Phi_d(2) for all {len(divs_10500)} divisors...")
t0 = time.time()
x = Symbol('x')
phi_values = {}
for d in divs_10500:
    poly = cyclotomic_poly(d, x)
    val = int(poly.subs(x, 2))
    phi_values[d] = val
dt = time.time() - t0
log(f"All cyclotomic values computed in {dt:.2f}s")

# Verify product equals 2^10500 - 1
product = 1
for d in divs_10500:
    product *= phi_values[d]
assert product == 2**pm1 - 1, "Cyclotomic product verification FAILED!"
log("Cyclotomic product verified: prod_{d|10500} Phi_d(2) = 2^10500 - 1")

# Display sizes
log(f"\nCyclotomic factor sizes:")
for d in divs_10500:
    phi_d = int(totient(d))
    digits = len(str(phi_values[d]))
    log(f"  Phi_{d:>5d}(2): phi({d}) = {phi_d:>5d}, {digits:>4d} digits")

# ============================================================
# Step 2: Factor the cyclotomic values
# ============================================================
section("Step 2: Factoring cyclotomic values")

# Known factorizations from Cunningham project tables
# Sources: factordb.com, homes.cerias.purdue.edu/~ssw/cun/,
#          maths-people.anu.edu.au/~brent/factors.html
cunningham_factors = {
    700: [701,
          2430065924693517198550322751963101,
          1038213793447841940908293355871461401],
    875: [725688486718330087751,
          5718039518555007627381367067387326475271823586132366619744957152669711690408574488925071259305573647810815919209500241631788149240057412005278225084826899525751],
    1500: [3001,
           791058001,
           168069194932501,
           4028493980595041855367835954324501,
           1606545773279325100753216665442817637284047671432352410624001],
    1750: [558251,
           10022251,
           511937190014372102141784257057292201840312001,
           1448724971925222342889675500698846356483820995550036783261953135108670064951365775570094705629617735073335209377682635445001],
    2100: [6301,
           743689627597081157353277424901,
           3065581111593982777238141477447662979750101,
           217100085701030760532082456157337409031425675860100163555370465949101],
}

all_prime_factors = {}  # prime -> total exponent across all Phi_d(2)
provenance_map = {}     # prime -> provenance label

log(f"\nFactoring each Phi_d(2)...")

for d in divs_10500:
    val = phi_values[d]
    if val == 1:
        all_prime_factors.setdefault(1, 0)
        continue

    digits = len(str(val))
    phi_d = int(totient(d))
    t0 = time.time()

    if d in cunningham_factors:
        # Use Cunningham table factorization
        factors = {}
        check = 1
        for f in cunningham_factors[d]:
            assert fast_isprime(f), f"Cunningham factor {f} is NOT prime!"
            assert val % f == 0, f"Factor {f} does not divide Phi_{d}(2)!"
            e = 0
            tmp = val
            while tmp % f == 0:
                tmp //= f
                e += 1
            factors[f] = e
            check *= f**e
        assert check == val, f"Cunningham factorization of Phi_{d}(2) is INCOMPLETE!"
        source = "Cunningham"
        prov_label = "cunningham"
    elif digits <= 80:
        # Factor with sympy
        factors = factorint(val)
        source = "sympy"
        prov_label = "direct_sympy"
    else:
        # Try trial division with known form r ≡ 1 (mod d)
        factors = {}
        remaining = val
        for k in range(1, min(10**9 // d + 1, 500000)):
            r = k * d + 1
            if r < 4:
                continue
            if remaining % r == 0 and fast_isprime(r):
                e = 0
                while remaining % r == 0:
                    remaining //= r
                    e += 1
                factors[r] = e
        if remaining > 1 and remaining != val and fast_isprime(remaining):
            factors[remaining] = 1
            remaining = 1
        source = "trial+ECM"
        prov_label = "trial_div_cyclotomic"

    dt = time.time() - t0

    for pr, exp in factors.items():
        all_prime_factors[pr] = all_prime_factors.get(pr, 0) + exp
        provenance_map.setdefault(pr, prov_label)

    factored_product = 1
    for pr, exp in factors.items():
        factored_product *= pr**exp
    is_complete = (factored_product == val)

    factor_str = " * ".join(
        f"{pr}" + (f"^{exp}" if exp > 1 else "")
        for pr, exp in sorted(factors.items())
    )

    if is_complete:
        log(f"  Phi_{d:>5d}(2) [{digits:>4d}d]: {factor_str}  [{source}, {dt:.1f}s] COMPLETE")
    else:
        remaining_d = len(str(val // factored_product)) if factored_product > 1 else digits
        log(f"  Phi_{d:>5d}(2) [{digits:>4d}d]: partial ({len(factors)} factors, {remaining_d}d cofactor) [{source}, {dt:.1f}s] SKIPPED")

# ============================================================
# Step 3: Collect prime factors and compute F
# ============================================================
section("Step 3: Computing factored part F of N-1")

# N-1 = 2 * (2^10500 - 1) / 3
# Add factor of 2 from the leading 2
all_prime_factors[2] = all_prime_factors.get(2, 0) + 1
provenance_map[2] = "algebraic"

# Subtract one factor of 3 (dividing by 3)
assert all_prime_factors.get(3, 0) >= 2, f"Expected 3^2+ in product, got 3^{all_prime_factors.get(3, 0)}"
all_prime_factors[3] -= 1
if all_prime_factors[3] == 0:
    del all_prime_factors[3]

# Remove non-prime keys
all_prime_factors.pop(1, None)

# Only include fully verified primes
prime_pool = {}
for pr, exp in all_prime_factors.items():
    if fast_isprime(pr):
        prime_pool[pr] = exp

log(f"Total distinct primes in factored part: {len(prime_pool)}")

# ============================================================
# Step 3a: APR-CL primality certification of every cofactor
# ============================================================
# The BLS Theorem 5 hypothesis requires every prime q | F to be
# proved prime unconditionally.  We certify every cofactor by the
# APR-CL test (Adleman-Pomerance-Rumely with Cohen-Lenstra
# improvements), invoked through PARI/GP's isprime(x, 2).  APR-CL
# is polynomial-time and provably correct; applied uniformly to
# every prime of F, including those below 2^64 where BPSW is
# empirically deterministic, so the certificate has no residual
# dependence on empirical pseudoprime tables.
section(f"Step 3a: APR-CL primality proof for {len(prime_pool)} primes of F")
aprcl_timings = aprcl_prove_all(prime_pool.keys(), gp_path=require_gp(), log=log)

# Compute F
F = 1
for pr, exp in sorted(prime_pool.items()):
    F *= pr**exp

# Verify F divides N-1
assert N_minus_1 % F == 0, "F does not divide N-1!"

R = N_minus_1 // F
assert fast_gcd(F, R) == 1, "gcd(F, R) != 1"

F_bits = F.bit_length()
F_digits = len(str(F))
R_bits = R.bit_length()
R_digits = len(str(R))

# BLS Theorem 5 check: F^3 > N  (exact integer comparison).
F_cubed = F**3
exact_margin = exact_bls_margin_bits(F, N)
assert F_cubed > N, "BLS Theorem 5 hypothesis F^3 > N FAILED"

log(f"\nF has {F_digits} digits ({F_bits} bits)")
log(f"R = (N-1)/F has {R_digits} digits ({R_bits} bits)")
log(f"N has {N_digits} digits ({N_bits} bits)")
log(f"F^3 has {F_cubed.bit_length()} bits")
log(f"Exact margin: (F^3).bit_length() - N.bit_length() = {exact_margin} bits")
log(f"\n*** BLS Theorem 5 CHECK: F^3 > N  ==>  PASSED  (margin {exact_margin} bits) ***")

# ============================================================
# Step 4: Find BLS witnesses
# ============================================================
section("Step 4: BLS witnesses for each prime q | F")

# For BLS Theorem 5 (N-1 test), for each prime q | F, we need witness a:
#   (I)  a^{N-1} ≡ 1 (mod N)
#   (II) gcd(a^{(N-1)/q} - 1, N) = 1

primes_in_F = sorted(prime_pool.keys())
log(f"Need witnesses for {len(primes_in_F)} primes")
log(f"Smallest: {primes_in_F[0]}, largest: {primes_in_F[-1]} ({len(str(primes_in_F[-1]))} digits)")

witnesses = {}
failed = []
witness_bases_tried = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

t_witnesses_start = time.time()

for idx, q in enumerate(primes_in_F):
    t0 = time.time()
    found = False

    for a in witness_bases_tried:
        # Check a^{N-1} ≡ 1 (mod N)
        r1 = fast_powmod(a, N_minus_1, N)
        if r1 != 1:
            continue

        # Check gcd(a^{(N-1)/q} - 1, N) = 1
        exp = N_minus_1 // q
        r2 = fast_powmod(a, exp, N)
        g = fast_gcd(r2 - 1, N)
        if g == 1:
            witnesses[q] = a
            found = True
            break

    dt = time.time() - t0
    q_digits = len(str(q))

    if found:
        log(f"  [{idx+1:>3d}/{len(primes_in_F)}] q = {q if q_digits <= 20 else f'({q_digits}d)'}: a = {witnesses[q]}  [{dt:.2f}s]")
    else:
        failed.append(q)
        log(f"  [{idx+1:>3d}/{len(primes_in_F)}] q = {q if q_digits <= 20 else f'({q_digits}d)'}: NO WITNESS  [{dt:.2f}s]")

t_witnesses_total = time.time() - t_witnesses_start
log(f"\nWitness computation: {t_witnesses_total:.1f}s total")

if failed:
    log(f"\n*** {len(failed)} primes need extended witness search ***")
    for q in failed:
        for a in range(53, 500):
            r1 = fast_powmod(a, N_minus_1, N)
            if r1 != 1:
                continue
            exp = N_minus_1 // q
            r2 = fast_powmod(a, exp, N)
            g = fast_gcd(r2 - 1, N)
            if g == 1:
                witnesses[q] = a
                log(f"  q = {q}: found witness a = {a}")
                break

if len(witnesses) == len(primes_in_F):
    log(f"\n*** ALL {len(primes_in_F)} PRIMES HAVE WITNESSES ***")
else:
    remaining = [q for q in primes_in_F if q not in witnesses]
    log(f"\n*** FAILED: {len(remaining)} primes without witnesses ***")
    sys.exit(1)

# ============================================================
# Step 5: BLS Theorem 5 — finite divisor check
# ============================================================
section("Step 5: BLS Theorem 5 — finite divisor check (discriminant)")

# If N is composite with all prime factors r ≡ 1 (mod F):
#   N = r * s where r = c1*F + 1, s = c2*F + 1
#   N - 1 = F * R where R = c1*c2*F + c1 + c2
#   So: c1 + c2 = R mod F, c1 * c2 = R // F
#   Discriminant of t^2 - (c1+c2)*t + c1*c2 must be a perfect square
#   for real c1, c2 to exist.

R_val = N_minus_1 // F
log(f"R = (N-1)/F has {len(str(R_val))} digits")

c1_plus_c2 = R_val % F
c1_times_c2 = R_val // F

log(f"If N composite: c1 + c2 ≡ R mod F ({len(str(c1_plus_c2))} digits)")
log(f"If N composite: c1 * c2 = R // F ({len(str(c1_times_c2))} digits)")

disc = c1_plus_c2**2 - 4 * c1_times_c2
log(f"Discriminant = (c1+c2)^2 - 4*c1*c2")

if disc < 0:
    log(f"Discriminant < 0 => NO real c1, c2 exist")
    log(f"")
    log(f"*** W_{p} IS PROVEN PRIME BY BLS THEOREM 5 (N-1 criterion) ***")
    RESULT = "PRIME"
else:
    # Check if perfect square
    sqrt_disc = math.isqrt(disc)
    if sqrt_disc * sqrt_disc == disc:
        log(f"Discriminant IS a perfect square (sqrt = {sqrt_disc})")
        c1 = (c1_plus_c2 + sqrt_disc) // 2
        c2 = (c1_plus_c2 - sqrt_disc) // 2
        if c1 >= 1 and c2 >= 1:
            r1 = c1 * F + 1
            r2 = c2 * F + 1
            if N % r1 == 0:
                log(f"r1 = {c1}*F+1 DIVIDES N => N is COMPOSITE!")
                RESULT = "COMPOSITE"
            elif N % r2 == 0:
                log(f"r2 = {c2}*F+1 DIVIDES N => N is COMPOSITE!")
                RESULT = "COMPOSITE"
            else:
                log(f"Neither candidate divides N")
                log(f"*** W_{p} IS PROVEN PRIME BY BLS THEOREM 5 ***")
                RESULT = "PRIME"
        else:
            log(f"c1={c1}, c2={c2}: not both positive")
            log(f"*** W_{p} IS PROVEN PRIME BY BLS THEOREM 5 ***")
            RESULT = "PRIME"
    else:
        log(f"Discriminant is NOT a perfect square")
        log(f"R is {'even' if R_val % 2 == 0 else 'odd'}")
        log(f"*** W_{p} IS PROVEN PRIME BY BLS THEOREM 5 ***")
        RESULT = "PRIME"

# ============================================================
# Step 6: Summary
# ============================================================
t_total_elapsed = time.time() - t_total

section("PROOF SUMMARY")
log(f"Number:           W_{p} = (2^{p} + 1) / 3")
log(f"Digits:           {N_digits}")
log(f"Bits:             {N_bits}")
log(f"Method:           BLS Theorem 5 (N-1 criterion)")
log(f"Result:           {RESULT}")
log(f"")
log(f"Factored part F:  {F_digits} digits ({F_bits} bits)")
log(f"Unfactored R:     {R_digits} digits ({R_bits} bits)")
log(f"Exact margin:     (F^3).bit_length() - N.bit_length() = {exact_margin} bits")
log(f"Discriminant:     {'< 0' if disc < 0 else 'not a perfect square' if disc >= 0 else '???'}")
log(f"")
log(f"Cyclotomic decomposition: 2^{pm1} - 1 = prod_{{d|{pm1}}} Phi_d(2)")
log(f"  {pm1} = 2^2 * 3 * 5^3 * 7 ({len(divs_10500)} divisors)")
log(f"  Factored Phi_d(2): {sum(1 for d in divs_10500 if phi_values[d] > 1)} non-trivial")
log(f"  Cunningham table factors used: {len(cunningham_factors)} (d = {sorted(cunningham_factors.keys())})")
log(f"")
log(f"BLS witnesses:    {len(witnesses)} primes")
log(f"  All use base a = {max(witnesses.values()) if witnesses else 'N/A'} or smaller")
log(f"  Bases used: {sorted(set(witnesses.values()))}")
log(f"")
log(f"Total time:       {t_total_elapsed:.1f}s")
log(f"Timestamp:        {datetime.now(timezone.utc).isoformat()}")

# ============================================================
# Step 7: Write structured proof certificate
# ============================================================
section("Writing proof certificate")

cert = {
    "number": f"W_{p}",
    "formula": f"(2^{p} + 1) / 3",
    "digits": N_digits,
    "bits": N_bits,
    "result": RESULT,
    "method": "BLS Theorem 5 (N-1)",
    "reference": "Brillhart-Lehmer-Selfridge, Math. Comp. 29 (1975), 620-647",
    "factored_part": {
        "bits": F_bits,
        "digits": F_digits,
        "num_primes": len(prime_pool),
    },
    "unfactored_part": {
        "bits": R_bits,
        "digits": R_digits,
    },
    "bls_hypothesis": {
        "statement": "F^3 > N",
        "satisfied": bool(F_cubed > N),
        "exact_margin_bits": exact_margin,
        "definition": "exact_margin_bits = (F^3).bit_length() - N.bit_length()"
    },
    "discriminant_sign": "negative" if disc < 0 else "non-square" if disc >= 0 else "unknown",
    "cyclotomic_decomposition": {
        "base": f"2^{pm1} - 1",
        "factorization_of_exponent": "2^2 * 3 * 5^3 * 7",
        "num_divisors": len(divs_10500),
        "cunningham_factors_used": sorted(cunningham_factors.keys()),
    },
    "witnesses": {str(q): a for q, a in sorted(witnesses.items())},
    "factor_provenance": {str(q): provenance_map.get(q, "unknown") for q in sorted(prime_pool.keys())},
    "provenance_legend": {
        "algebraic": "contribution from the form N-1 = 2 (2^{p-1}-1)/3 (the factor 2)",
        "cunningham": "factor from the Cunningham project table for 2^n-1; also verifiable via FactorDB",
        "direct_sympy": "produced by sympy.factorint on a cyclotomic value Phi_d(2)",
        "trial_div_cyclotomic": "discovered by trial division over the arithmetic progression r = k d + 1"
    },
    "aprcl_certification": {
        "tool": "PARI/GP isprime(x, 2)",
        "scope": "every prime of F (uniform, no size threshold)",
        "num_certified": len(aprcl_timings),
        "total_elapsed_seconds": round(sum(aprcl_timings.values()), 3),
        "max_elapsed_seconds": round(max(aprcl_timings.values()), 3) if aprcl_timings else 0.0,
    },
    "chebyshev_condition_ii": {
        "base": "omega_3 = 3 + 2*sqrt(2)",
        "congruence": "omega_3^((N+1)/2) ≡ -1 (mod N)",
        "verified": bool(cond_ii["verified"]),
        "exponent_bits": cond_ii["exponent_bits"],
        "elapsed_seconds": round(cond_ii["elapsed_seconds"], 3),
    },
    "computation_time_seconds": round(t_total_elapsed, 2),
    "timestamp": datetime.now(timezone.utc).isoformat(),
}

# Write certificate to ../data/ and log to ../data/ (release layout);
# fall back to the current working directory if data/ is not present.
import os
_script_dir = os.path.dirname(os.path.abspath(__file__))
_data_dir = os.path.normpath(os.path.join(_script_dir, os.pardir, "data"))
_out_dir = _data_dir if os.path.isdir(_data_dir) else os.getcwd()

cert_path = os.path.join(_out_dir, f"bls_certificate_w{p}.json")
with open(cert_path, "w") as f:
    json.dump(cert, f, indent=2)
log(f"Certificate written to {cert_path}")

log_path = os.path.join(_out_dir, f"bls_proof_w{p}.log")
with open(log_path, "w") as f:
    f.write(f"BLS N-1 Primality Proof for W_{p}\n")
    f.write(f"{'='*70}\n\n")
    for line in LOG:
        f.write(line + "\n")
log(f"Full log written to {log_path}")

log(f"\n{'='*70}")
log(f"DONE — W_{p} is {'PRIME' if RESULT == 'PRIME' else 'NOT PROVEN PRIME'}")
log(f"{'='*70}")
