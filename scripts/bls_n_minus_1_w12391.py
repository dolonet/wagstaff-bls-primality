#!/usr/bin/env python3
"""
BLS N-1 primality proof for W_12391 = (2^12391 + 1) / 3

Strategy: Factor N-1 = 2*(2^12390 - 1)/3 using the cyclotomic decomposition
of 2^12390 - 1. Since 12390 = 2 * 3 * 5 * 7 * 59, we have 32 cyclotomic
factors Phi_d(2) for d | 12390.

All necessary factorizations are known from the Cunningham project tables,
factordb.com, and direct computation. Five large cyclotomic values use
Cunningham/factordb data (d = 885, 1239, 1770, 2065, 2478), and Phi_590(2)
is itself prime (71 digits).

The factored part F has ~5400 bits, exceeding the BLS threshold of 4131 bits.

References:
  [BLS75] Brillhart, Lehmer, Selfridge, "New primality criteria and
          factorizations of 2^m +/- 1", Math. Comp. 29 (1975), 620-647.
"""

import sys
import time
import math
import json
from datetime import datetime, timezone

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

LOG = []

def log(msg):
    ts = datetime.now(timezone.utc).strftime("%H:%M:%S.%f")[:-3]
    line = f"[{ts}] {msg}"
    print(line)
    sys.stdout.flush()
    LOG.append(line)

def section(title):
    log("=" * 70)
    log(title)
    log("=" * 70)

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
# Step 0: Define N = W_12391
# ============================================================
t_total = time.time()
section("Step 0: Define N = W_12391")

p = 12391
N = (2**p + 1) // 3
N_digits = len(str(N))
N_bits = N.bit_length()

log(f"p = {p}")
log(f"W_{p} has {N_digits} digits ({N_bits} bits)")
log(f"N mod 8 = {N % 8}")

N_minus_1 = N - 1
assert N_minus_1 == (2**p - 2) // 3
log(f"N-1 = (2^{p} - 2) / 3")
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
# Step 1: Cyclotomic decomposition of 2^12390 - 1
# ============================================================
section("Step 1: Cyclotomic decomposition of 2^12390 - 1")

log(f"12390 = 2 * 3 * 5 * 7 * 59")

pm1 = p - 1
pm1_factors = {2: 1, 3: 1, 5: 1, 7: 1, 59: 1}

def gen_divisors(factors):
    divs = [1]
    for pr, exp in factors.items():
        new_divs = []
        for d in divs:
            for e in range(exp + 1):
                new_divs.append(d * pr**e)
        divs = new_divs
    return sorted(divs)

divs = gen_divisors(pm1_factors)
log(f"Number of divisors of 12390: {len(divs)}")
log(f"Divisors: {divs}")

log(f"\nComputing Phi_d(2) for all {len(divs)} divisors...")
t0 = time.time()
x = Symbol('x')
phi_values = {}
for d in divs:
    poly = cyclotomic_poly(d, x)
    val = int(poly.subs(x, 2))
    phi_values[d] = val
dt = time.time() - t0
log(f"All cyclotomic values computed in {dt:.2f}s")

product = 1
for d in divs:
    product *= phi_values[d]
assert product == 2**pm1 - 1
log("Cyclotomic product verified: prod_{d|12390} Phi_d(2) = 2^12390 - 1")

log(f"\nCyclotomic factor sizes:")
for d in divs:
    phi_d = int(totient(d))
    digits = len(str(phi_values[d]))
    log(f"  Phi_{d:>5d}(2): phi({d}) = {phi_d:>5d}, {digits:>4d} digits")

# ============================================================
# Step 2: Factor the cyclotomic values
# ============================================================
section("Step 2: Factoring cyclotomic values")

# All known factorizations.
# Sources: FactorDB and the Cunningham project tables (see provenance
# mapping below).  Every factor is re-certified prime by APR-CL (PARI/GP
# isprime(x, 2)) before entering F, so the certificate does not depend
# on the archives' correctness.
known_factors = {
    # d: [list of prime factors]
    # Phi_590(2) = 71 digits — PRIME (handled separately)
    # Phi_295(2) = 70 digits, 7 prime factors (sympy verified)
    295: [4721, 132751, 5794391, 128818831, 3812358161,
          452824604065751, 4410975230650827973711],
    # Phi_413(2) = 105 digits — SKIPPED (not needed for BLS threshold)
    # Phi_826(2) = 105 digits, 3 prime factors
    826: [827,
          170735974773267443,
          6043930497790503973481076813462520042997083539133970912065745573049492802026928038019],
    # Phi_885(2) = 140 digits, 2 prime factors
    885: [2756788662198217256191,
          29293922760297928248078770598052610852712405023442274675016645233570497812348935888398085451640486701886767775742174681],
    # Phi_1239(2) = 210 digits, 4 prime factors
    1239: [263483263,
           1102272524932426318899113975892645895609,
           167087803778100685137282215256230683687068949094424492832211188268061451903148929,
           11763111754911034189958819922265626030162852411746099534272282742979110559598919617],
    # Phi_1770(2) = 140 digits, 4 prime factors
    1770: [516266521,
           873791632531,
           2354488203481,
           34685790485740246824716440792348382055127879712328545535166847444199845117697592973437487413465343206165441],
    # Phi_2065(2) = 420 digits, 5 prime factors (includes p=12391!)
    2065: [12391, 161071, 107429561,
           34612315434702943134556428791471],
    # The 371-digit cofactor will be computed and verified as prime.
    # Phi_2478(2) = 210 digits, 3 prime factors
    2478: [359931645741056789631351742091222797125100809698191524292297,
           616118417295048578293181955408006300135730334724584315102344777,
           1120555975329453797460758793161622336521020113400670993101603640935259508257738747978899],
}

# Per-d provenance labels for known_factors.  These are the original
# sources from which each literal factorization was first obtained
# during development; the literals were then independently re-certified
# by APR-CL.  See "factor_provenance" in the output certificate.
known_factors_source = {
    295: "factordb",
    826: "cunningham",
    885: "cunningham",
    1239: "cunningham",
    1770: "cunningham",
    2065: "factordb",      # small-prime entries 12391, 161071, 107429561 via FactorDB
    2478: "cunningham",
}

# Cyclotomic values we skip factoring — not needed to reach BLS threshold.
# Phi_413(2) is 105 digits; partial factoring found 2006647231 but the
# 96-digit cofactor is composite and we already have +900 bit margin without it.
# Phi_4130, 6195, 12390 are 419-839 digits — completely infeasible to factor.
skip_factoring = {413, 4130, 6195, 12390}

all_prime_factors = {}
provenance_map = {}   # prime -> provenance label

for d in divs:
    val = phi_values[d]
    if val == 1:
        continue

    digits = len(str(val))
    phi_d = int(totient(d))
    t0 = time.time()

    if d in skip_factoring:
        log(f"  Phi_{d:>5d}(2) [{digits:>4d}d]: SKIPPED (not needed for BLS threshold)")
        continue

    computed_cofactor = None  # set when we compute a residual cofactor ourselves

    if d == 590:
        # Phi_590(2) is itself prime
        assert fast_isprime(val), "Phi_590(2) is not prime!"
        factors = {val: 1}
        source = "prime"
        per_prime_source = {val: "cyclotomic_prime"}
    elif d in known_factors and known_factors[d][0] is not None:
        factor_list = known_factors[d]
        factors = {}
        check = 1
        base_source = known_factors_source.get(d, "cunningham")

        if d == 2065:
            # Compute the 371-digit cofactor
            remaining = val
            for f in factor_list:
                while remaining % f == 0:
                    remaining //= f
            assert fast_isprime(remaining), "Phi_2065 cofactor not prime!"
            factor_list = factor_list + [remaining]
            computed_cofactor = remaining
            log(f"    Phi_2065 cofactor ({len(str(remaining))}d): verified prime")

        per_prime_source = {}
        for f in factor_list:
            assert fast_isprime(f), f"Factor {f} is NOT prime!"
            assert val % f == 0, f"Factor {f} does not divide Phi_{d}(2)!"
            e = 0
            tmp = val
            while tmp % f == 0:
                tmp //= f
                e += 1
            factors[f] = e
            check *= f**e
            if f == computed_cofactor:
                per_prime_source[f] = "residual_prime_aprcl"
            else:
                per_prime_source[f] = base_source
        assert check == val, f"Factorization of Phi_{d}(2) incomplete!"
        source = base_source
    elif digits <= 80:
        factors = factorint(val)
        source = "sympy"
        per_prime_source = {f: "direct_sympy" for f in factors}
    else:
        # Try trial division for medium ones, then sympy
        factors = factorint(val)
        source = "sympy"
        per_prime_source = {f: "direct_sympy" for f in factors}

    dt = time.time() - t0

    for pr, exp in factors.items():
        all_prime_factors[pr] = all_prime_factors.get(pr, 0) + exp
        provenance_map.setdefault(pr, per_prime_source.get(pr, "unknown"))

    factored_product = 1
    for pr, exp in factors.items():
        factored_product *= pr**exp
    is_complete = (factored_product == val)

    nf = len(factors)
    log(f"  Phi_{d:>5d}(2) [{digits:>4d}d]: {nf} factors [{source}, {dt:.1f}s] {'COMPLETE' if is_complete else 'PARTIAL'}")

# ============================================================
# Step 3: Compute F
# ============================================================
section("Step 3: Computing factored part F of N-1")

all_prime_factors[2] = all_prime_factors.get(2, 0) + 1
provenance_map[2] = "algebraic"
assert all_prime_factors.get(3, 0) >= 2
all_prime_factors[3] -= 1
if all_prime_factors[3] == 0:
    del all_prime_factors[3]
all_prime_factors.pop(1, None)

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
# every prime of F regardless of size.
section(f"Step 3a: APR-CL primality proof for {len(prime_pool)} primes of F")
aprcl_timings = aprcl_prove_all(prime_pool.keys(), gp_path=require_gp(), log=log)

F = 1
for pr, exp in sorted(prime_pool.items()):
    F *= pr**exp

assert N_minus_1 % F == 0
R = N_minus_1 // F
assert fast_gcd(F, R) == 1

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
# Step 4: BLS witnesses
# ============================================================
section("Step 4: BLS witnesses for each prime q | F")

primes_in_F = sorted(prime_pool.keys())
log(f"Need witnesses for {len(primes_in_F)} primes")
log(f"Smallest: {primes_in_F[0]}, largest: ({len(str(primes_in_F[-1]))} digits)")

witnesses = {}
failed = []
witness_bases = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]

t_witnesses_start = time.time()

for idx, q in enumerate(primes_in_F):
    t0 = time.time()
    found = False

    for a in witness_bases:
        r1 = fast_powmod(a, N_minus_1, N)
        if r1 != 1:
            continue
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
    log(f"\n*** {len(failed)} primes need extended search ***")
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
# Step 5: Discriminant check
# ============================================================
section("Step 5: BLS Theorem 5 — finite divisor check (discriminant)")

R_val = N_minus_1 // F
log(f"R = (N-1)/F has {len(str(R_val))} digits")

c1_plus_c2 = R_val % F
c1_times_c2 = R_val // F

log(f"If N composite: c1 + c2 ({len(str(c1_plus_c2))} digits)")
log(f"If N composite: c1 * c2 ({len(str(c1_times_c2))} digits)")
log(f"Discriminant = (c1+c2)^2 - 4*c1*c2")

disc = c1_plus_c2**2 - 4 * c1_times_c2

if disc < 0:
    log(f"Discriminant < 0 => NO real c1, c2 exist")
    log(f"\n*** W_{p} IS PROVEN PRIME BY BLS THEOREM 5 (N-1 criterion) ***")
    RESULT = "PRIME"
    DISC_STATUS = "negative"
else:
    sqrt_disc = math.isqrt(disc)
    if sqrt_disc * sqrt_disc == disc:
        log(f"Discriminant IS a perfect square")
        c1 = (c1_plus_c2 + sqrt_disc) // 2
        c2 = (c1_plus_c2 - sqrt_disc) // 2
        if c1 >= 1 and c2 >= 1:
            r1 = c1 * F + 1
            r2 = c2 * F + 1
            if N % r1 == 0:
                log(f"r1 DIVIDES N => COMPOSITE!")
                RESULT = "COMPOSITE"
                DISC_STATUS = "square-and-factors"
            elif N % r2 == 0:
                log(f"r2 DIVIDES N => COMPOSITE!")
                RESULT = "COMPOSITE"
                DISC_STATUS = "square-and-factors"
            else:
                log(f"Neither candidate divides N")
                log(f"\n*** W_{p} IS PROVEN PRIME BY BLS THEOREM 5 ***")
                RESULT = "PRIME"
                DISC_STATUS = "square-but-no-divisor"
        else:
            log(f"\n*** W_{p} IS PROVEN PRIME BY BLS THEOREM 5 ***")
            RESULT = "PRIME"
            DISC_STATUS = "square-but-nonpositive-c"
    else:
        log(f"Discriminant is NOT a perfect square")
        log(f"R is {'even' if R_val % 2 == 0 else 'odd'}")
        log(f"\n*** W_{p} IS PROVEN PRIME BY BLS THEOREM 5 ***")
        RESULT = "PRIME"
        DISC_STATUS = "non-square"

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
log(f"")
log(f"Cyclotomic decomposition: 2^{pm1} - 1 = prod_{{d|{pm1}}} Phi_d(2)")
log(f"  {pm1} = 2 * 3 * 5 * 7 * 59 ({len(divs)} divisors)")
log(f"")
log(f"BLS witnesses:    {len(witnesses)} primes")
log(f"  Bases used: {sorted(set(witnesses.values()))}")
log(f"")
log(f"Total time:       {t_total_elapsed:.1f}s")
log(f"Timestamp:        {datetime.now(timezone.utc).isoformat()}")

# ============================================================
# Step 7: Save certificate and log
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
    "discriminant_sign": DISC_STATUS,
    "cyclotomic_decomposition": {
        "base": f"2^{pm1} - 1",
        "factorization_of_exponent": "2 * 3 * 5 * 7 * 59",
        "num_divisors": len(divs),
    },
    "witnesses": {str(q): a for q, a in sorted(witnesses.items())},
    "factor_provenance": {str(q): provenance_map.get(q, "unknown") for q in sorted(prime_pool.keys())},
    "provenance_legend": {
        "algebraic": "contribution from the form N-1 = 2 (2^{p-1}-1)/3 (the factor 2)",
        "cunningham": "factor from the Cunningham project table for 2^n-1",
        "factordb": "factor obtained by FactorDB lookup (used as a discovery aid; re-certified prime by APR-CL before entering F)",
        "cyclotomic_prime": "Phi_d(2) is itself prime (no further factoring needed)",
        "residual_prime_aprcl": "residual cofactor of Phi_d(2) after removing known factors; verified prime by APR-CL before entering F",
        "direct_sympy": "produced by sympy.factorint on a cyclotomic value Phi_d(2)"
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
