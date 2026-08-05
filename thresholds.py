#!/usr/bin/env python3
"""
thresholds.py -- certified thresholds for Corollary 1.2 (Table 1).

The criterion of Corollary 1.2 is

    n^{1/4}  >  K_r * log n,        K_r := C * 2^{r/2} / (C_Artin * Pi_r),

with C = 11.4 the constant of Theorem 1.1, r = omega(k), and

    Pi_r := prod_{t=1}^{r} (1 - 1/(q_t - 1)),   q_1 < q_2 < ... the odd primes.

Pi_r is held as an EXACT rational, and every floating step is directed so that
the reported threshold is certified rather than merely indicated:

  * C_Artin is replaced by a strict LOWER bound, so K_r is over-estimated;
  * K_r is rounded UP;
  * the crossover is bracketed and the reported n_0 rounded UP to 3 s.f.;
  * the certificate re-checks n_0^{1/4} > K_r log n_0 at the ROUNDED value,
    with the left side rounded DOWN and the right side rounded UP.

Since n^{1/4}/log n is increasing for n > e^4, verifying the inequality at n_0
certifies it for every n >= n_0.  Finally n_0 is raised to k_r^2 where the
range k <= sqrt(n) would otherwise make the statement vacuous.

Requires: mpmath, sympy.  Run:  python3 thresholds.py
"""
from fractions import Fraction
from mpmath import mp, mpf, log as mlog, floor as mfloor, ceil as mceil
from sympy import primerange

mp.dps = 60

# ---------------------------------------------------------------- inputs ----
# Every input is stated exactly, with its provenance, so that a reader can check
# the certificate against the paper without re-deriving anything.
#
#   C          Theorem 1.1's error constant, 11.4.  Printed value; the proof
#              gives 2*sqrt(9*c_0) + 1 = 11.3122..., and 11.4 is its round-up.
#   CARTIN_LO  strict LOWER bound for Artin's constant 0.3739558136...
#              (Wrench, Math. Comp. 15 (1961) 396-398).  Using a lower bound
#              makes K_r larger, i.e. the criterion harder, i.e. safe.
#   ODD        the odd primes, in order; q_1 = 3.
#   RMAX       largest omega(k) tabulated.
#
# The criterion also involves c_0 = zeta(2)^2 zeta(3)/(zeta(4) zeta(6)) through
# Theorem 1.1's constant, but c_0 does not appear here directly: it is already
# folded into C.
C = Fraction(114, 10)                 # Theorem 1.1 constant, exact
CARTIN_LO = Fraction(3739558, 10**7)  # strict lower bound for Artin's constant
RMAX = 12

ODD = list(primerange(3, 20000))   # ample: Table 1 needs r <= 12, the
                                   # Corollary 1.4 cascade checks r <= 60 (k_r there is the primorial)


def pi_r(r):
    """Exact Pi_r."""
    p = Fraction(1)
    for t in range(r):
        p *= 1 - Fraction(1, ODD[t] - 1)
    return p


def k_min(r):
    """Smallest odd square-free modulus with omega = r (the odd primorial)."""
    p = 1
    for t in range(r):
        p *= ODD[t]
    return p


def K_up(r):
    """K_r, rounded upward."""
    base = Fraction(C, CARTIN_LO * pi_r(r))          # exact rational part
    K = mpf(base.numerator) / mpf(base.denominator) * mpf(2) ** (mpf(r) / 2)
    return K * (1 + mpf(10) ** (-40))                # nudge up


def crossover(r):
    """Least n with n^{1/4} > K_r log n, by bisection on the upper branch."""
    K = K_up(r)
    f = lambda n: n ** mpf('0.25') - K * mlog(n)
    lo, hi = mpf('1e6'), mpf('1e50')
    for _ in range(500):
        mid = (lo * hi) ** mpf('0.5')
        if f(mid) > 0:
            hi = mid
        else:
            lo = mid
    return hi


def round_up_3sf(x):
    e = int(mfloor(mlog(x, 10)))
    m = mceil(x / mpf(10) ** e * 100) / 100
    if m >= 10:
        m, e = m / 10, e + 1
    return m, e


def certify(r):
    n_err = crossover(r)
    n_vac = mpf(k_min(r)) ** 2          # else k <= sqrt(n) is vacuous at this r
    m, e = round_up_3sf(max(n_err, n_vac))
    n0 = m * mpf(10) ** e
    # re-check at the rounded value: LHS down, RHS up
    lhs = n0 ** mpf('0.25') * (1 - mpf(10) ** (-40))
    rhs = K_up(r) * mlog(n0) * (1 + mpf(10) ** (-40))
    binding = 'k^2' if n_vac >= n_err else 'error'
    return m, e, lhs > rhs or n_vac >= n_err, binding


# The rows as PRINTED IN TABLE 1 of the paper.  Every one is re-derived and
# asserted below, so this file fails loudly if the paper and the code drift.
TABLE_1 = {
    3:  (Fraction(5, 16),           105,                 1.08, 16),
    4:  (Fraction(9, 32),           1155,                8.13, 16),
    5:  (Fraction(33, 128),         15015,               5.59, 17),
    6:  (Fraction(495, 2048),       255255,              3.45, 18),
    7:  (Fraction(935, 4096),       4849845,             2.04, 19),
    8:  (Fraction(1785, 8192),      111546435,           1.15, 20),
    9:  (Fraction(6885, 32768),     3234846615,          6.10, 20),
    10: (Fraction(13311, 65536),    100280245065,        1.01, 22),
    11: (Fraction(51765, 262144),   3710369067405,       1.38, 25),
    12: (Fraction(403767, 2097152), 152125131763605,     2.32, 28),
}


def check_row(r):
    """Re-derive row r and assert it against the printed table."""
    pi_paper, k_paper, m_paper, e_paper = TABLE_1[r]
    assert pi_r(r) == pi_paper, f"r={r}: Pi_r is {pi_r(r)}, table prints {pi_paper}"
    assert k_min(r) == k_paper, f"r={r}: smallest k is {k_min(r)}, table prints {k_paper}"
    m, e, ok, binding = certify(r)
    assert ok, f"r={r}: criterion FAILS at the printed threshold"
    assert (float(m), e) == (m_paper, e_paper), \
        f"r={r}: certified n_0 is {float(m)}e{e}, table prints {m_paper}e{e_paper}"
    # the printed threshold must also clear Corollary 1.3's strengthened form
    n0 = mpf(m) * mpf(10) ** e
    assert n0 ** mpf('0.25') > K_up(r) * mlog(n0) + 1, \
        f"r={r}: printed n_0 fails the +1 criterion of Corollary 1.3"
    return m, e, binding


if __name__ == '__main__':
    print(f"{'r':>3} {'smallest k':>20} {'Pi_r (exact)':>18} {'Pi_r':>9} "
          f"{'K_r':>11} {'n_0':>12} {'binding':>10} {'ok':>4}")
    failures = []
    for r in range(3, RMAX + 1):
        try:
            m, e, binding = check_row(r)
            status = 'yes'
        except AssertionError as exc:
            failures.append(str(exc))
            m, e, binding, status = 0, 0, '-', 'NO'
        print(f"{r:>3} {k_min(r):>20,} {str(pi_r(r)):>18} {float(pi_r(r)):>9.6f} "
              f"{float(K_up(r)):>11.3f} {float(m):>6.2f}e{e:<5} {binding:>10} "
              f"{status:>4}")
    if failures:
        print("\nCERTIFICATE FAILED:")
        for f in failures:
            print("  ! " + f)
        raise SystemExit(1)
    print(f"\nAll {RMAX - 2} rows of Table 1 certified, and each also clears the")
    print("strengthened criterion of Corollary 1.3 (even moduli).")

    # ---------------------------------------------------------------------
    # Corollary 1.4: the uniform threshold N_0 = 6.10e20.
    #
    # Two halves.  (a) For r <= 9 the least admissible n is increasing in r,
    # so it is at most its value at r = 9, which is N_0 itself.  (b) For
    # r >= 10 the range k <= sqrt(n) forces n >= P_r^2, and the criterion
    # already holds there, so no constraint on n survives.  Half (b) is the
    # inequality sqrt(P_r) > 2 K_r log P_r for every r >= 10, proved by
    # induction with base r = 10 and the growth factor below.
    # ---------------------------------------------------------------------
    print("\nCorollary 1.4  (uniform threshold)")
    cas = []

    def cas_check(name, cond, detail=''):
        print(f"  {'ok  ' if cond else 'FAIL'}  {name}" + (f"   [{detail}]" if detail else ''))
        if not cond:
            cas.append(name)

    N0 = mpf('6.10e20')
    cas_check("K_r is increasing in r (so the least admissible n is too)",
              all(K_up(r) < K_up(r + 1) for r in range(0, 60)))
    worst = max(crossover(r) for r in range(0, 10))
    cas_check("every crossing for r <= 9 is at most N_0 = 6.10e20",
              worst <= N0, f"worst is r=9 at {float(worst):.4g}")
    cas_check("N_0 lies below the r=10 row of Table 1 (1.01e22)", N0 < mpf('1.01e22'))

    # base case and induction for the r >= 10 half
    C55 = mpf('11.4') / mpf('0.2054')     # >= 11.4/(C_Artin * log3/2), by Lemma 2.5

    def P(r):
        v = 1
        for q in ODD[:r]:
            v *= q
        return v

    def Psi(r):
        return (mpf(P(r)) / mpf(2) ** r) ** mpf('0.5') / \
               (2 * C55 * mlog(ODD[r - 1]) * mlog(mpf(P(r))))

    cas_check("base case Psi(10) > 1", Psi(10) > 1, f"Psi(10) = {float(Psi(10)):.7f}")
    g = (mpf('18.5') ** mpf('0.5')) * mlog(31) / mlog(62) / 2
    cas_check("growth factor sqrt(18.5)*(log31/log62)*(1/2) > 1",
              g > 1, f"value {float(g):.5f}")
    cas_check("its three ingredients hold: q_{r+1}>=37, Bertrand, q_r <= P_r/3",
              ODD[10] >= 37
              and all(ODD[t] < 2 * ODD[t - 1] for t in range(1, 2000))
              and all(3 * ODD[r - 1] <= P(r) for r in range(2, 60)))
    cas_check("hence sqrt(P_r) > 2 K_r log P_r directly, r = 10..60",
              all(mpf(P(r)) ** mpf('0.5') > 2 * K_up(r) * mlog(mpf(P(r)))
                  for r in range(10, 61)))
    cas_check("even case: criterion+1 holds at n = 4 P_r^2, r = 10..60",
              all((4 * mpf(P(r)) ** 2) ** mpf('0.25')
                  > K_up(r) * mlog(4 * mpf(P(r)) ** 2) + 1 for r in range(10, 61)))

    if cas:
        print("\nCERTIFICATE FAILED:")
        for f in cas:
            print("  ! " + f)
        raise SystemExit(1)
    print("\nN_0 = 6.10e20 certified for Corollary 1.4.")
    raise SystemExit(0)
