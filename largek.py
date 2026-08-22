#!/usr/bin/env python3
"""
largek.py -- certified constants for Proposition 5.2 (large k).

Theorem 1.1 balances the two error terms of the master bound (4.6) by the
choice (4.7) of truncation level z, which forces k <= sqrt(n).  Proposition 5.2
gives up the asymptotic and keeps only positivity, and can then take z far
smaller -- admitting k almost as large as n.  The paper leaves its thresholds
effective but uncomputed; this script computes them.

PART I  (Proposition 5.2(i), fixed r).  Take the CONSTANT level
z = 45 c_0 / delta_r, for which the tail term 9 c_0 n / z of (4.6) is exactly
delta_r n / 5 and the constraint z^2 k <= n reads k <= kappa_r n with
kappa_r := (delta_r / (45 c_0))^2.  Writing delta_r = C_Artin * Pi_r and
tau = 2^r, positivity needs

    (4/5) delta_r n  >  2 n^{5/8} log n + z tau sqrt(n) (log n)^2 + 2 z tau log n,

together with the hypotheses of (4.6) -- n >= 1e5 and 1 <= z <= n^{3/8} -- and
k_r <= kappa_r n so that an admissible k exists at all (k_r = odd primorial).

PART II  (Proposition 5.2(ii), uniform in r).  Take z = (45 c_0 / 0.2054)
log y(k), for which the tail term is delta(k) n / 5 with
delta(k) = 0.2054 / log y(k).  Bounding y(k) by (5.3), tau(k) by the divisor
bound (5.4), and requiring the three remaining terms of (4.6) to fall below
delta(k) n / 5 yields the absolute threshold N; C_3 is then the least constant
for which k <= n / (C_3 (log log n)^2) forces z^2 k <= n throughout n >= N.

Every floating step is directed so the reported constants are certified:

  * C_Artin is replaced by a strict LOWER bound, so delta_r is under-estimated
    (weakening the main term) and kappa_r is under-estimated (shrinking the
    claimed k-range);
  * c_0 is replaced by a strict UPPER bound, so z -- and with it the whole
    error side -- is over-estimated, and kappa_r shrinks further;
  * Pi_r is held as an EXACT rational throughout;
  * every threshold is rounded UP to 3 s.f. and the inequality re-checked at
    the ROUNDED value, left side down and right side up.

In both parts the left side is linear in n and the right side is
O(n^{3/4} log^2 n), so verifying at the threshold certifies every larger n.

Requires: mpmath, sympy.  Run:  python3 largek.py
"""
from mpmath import mp, mpf, log as mlog, floor as mfloor, ceil as mceil, zeta

mp.dps = 60

from thresholds import ODD, pi_r, k_min, CARTIN_LO, TABLE_1, RMAX

# ---------------------------------------------------------------- inputs ----
#   UP / DN    one-ulp nudges at 60 digits, used to direct every rounding.
#   C0_UP      strict UPPER bound for c_0 = zeta(2)^2 zeta(3)/(zeta(4) zeta(6))
#              = 2.953912... (Lemma 2.3).  c_0 enters only through z, and a
#              larger z makes both the error side and admissibility harder.
#   A          the coefficient (45 c_0 / 0.2054) of log y(k) in Part II's z,
#              where 0.2054 is the lower bound of Proposition 3.1(ii).
#   TAU_C      the divisor-bound constant of (5.4): tau(m) <= 8.5 m^{1/4}.
#   THETA_C    Rosser-Schoenfeld: theta(y) > 0.84 y for y >= 101, used in (5.3).
UP = 1 + mpf(10) ** (-40)
DN = 1 - mpf(10) ** (-40)

C0_UP = (zeta(2) ** 2 * zeta(3) / (zeta(4) * zeta(6))) * UP
A = 45 * C0_UP / mpf('0.2054')
TAU_C = mpf('8.5')
BETA = mpf('0.39851')          # constant of Lemma 2.2; scales the GRH remainder in (4.6)
THETA_C = mpf('0.84')

# Constants AS PRINTED IN THE PAPER.  Re-derived and asserted below, so this
# file fails loudly if the paper and the code drift apart.
PROP_52_I = {10: (7.37, 20), 11: (3.77, 21), 12: (1.90, 22)}
PROP_52_II_N = (1.01, 40)
PROP_52_II_C3 = mpf('4.54e5')
PROP_52_II_C2 = mpf('673.2')      # printed C_2', with C_3 = C_2'^2 <= 4.53e5


def round_up_3sf(x):
    e = int(mfloor(mlog(x, 10)))
    m = mceil(x / mpf(10) ** e * 100) / 100
    if m >= 10:
        m, e = m / 10, e + 1
    return m, e


# =========================================================== Part I =========
def delta_lo(r):
    """Strict lower bound for delta_r = C_Artin * Pi_r."""
    f = CARTIN_LO * pi_r(r)
    return mpf(f.numerator) / mpf(f.denominator) * DN


def z_const_up(r):
    """The constant truncation level z = 45 c_0 / delta_r, rounded UP."""
    return 45 * C0_UP / delta_lo(r) * UP


def kappa_lo(r):
    """kappa_r = (delta_r / (45 c_0))^2, rounded DOWN."""
    return (delta_lo(r) / (45 * C0_UP)) ** 2 * DN


def slack_I(n, r):
    """LHS (down) minus RHS (up) of Part I's positivity inequality."""
    z, tau = z_const_up(r), mpf(2) ** r
    lhs = mpf('0.8') * delta_lo(r) * n * DN
    rhs = (2 * n ** mpf('0.625') * mlog(n)
           + BETA * z * tau * n ** mpf('0.5') * mlog(n) ** 2
           + 2 * z * tau * mlog(n)) * UP
    return lhs - rhs


def crossover_I(r):
    lo, hi = mpf('1e5'), mpf('1e60')
    for _ in range(500):
        mid = (lo * hi) ** mpf('0.5')
        if slack_I(mid, r) > 0:
            hi = mid
        else:
            lo = mid
    return hi


def certify_I(r):
    n_err = crossover_I(r)
    n_adm = mpf(k_min(r)) / kappa_lo(r)
    m, e = round_up_3sf(max(n_err, n_adm, mpf('1e5')))
    n_r = m * mpf(10) ** e
    checks = {
        'positivity holds at n_r':        slack_I(n_r, r) > 0,
        'z <= n_r^{3/8}  (Lemma 4.1)':    z_const_up(r) <= n_r ** mpf('0.375'),
        'z >= 1          (Prop. 4.2)':    z_const_up(r) >= 1,
        'n_r >= 1e5      (Lemma 4.1)':    n_r >= mpf('1e5'),
        'an admissible k exists':         mpf(k_min(r)) <= kappa_lo(r) * n_r,
        'kappa_r n_r > sqrt(n_r)':        kappa_lo(r) * n_r > n_r ** mpf('0.5'),
        # surplus above log 2 lets the single term p = 2 be discarded, as in
        # Corollary 1.3, so an even-k analogue is available at the same row.
        'surplus exceeds log 2':          slack_I(n_r, r) > mlog(2),
    }
    return m, e, n_r, checks


# ========================================================== Part II =========
def y_up(n):
    """Upper bound for y(k) over k <= n, from (5.3)."""
    return max(mpf(101), (mlog(n) + mlog(2)) / THETA_C * UP)


def z_var_up(n):
    return A * mlog(y_up(n)) * UP


def delta_k_lo(n):
    """Lower bound for delta(k) = 0.2054 / log y(k), over k <= n."""
    return mpf('0.2054') / (mlog(y_up(n)) * UP) * DN


def slack_II(n):
    z, tau = z_var_up(n), TAU_C * n ** mpf('0.25') * UP
    lhs = delta_k_lo(n) * n / 5 * DN
    rhs = (2 * n ** mpf('0.625') * mlog(n)
           + BETA * z * tau * n ** mpf('0.5') * mlog(n) ** 2
           + 2 * z * tau * mlog(n)) * UP
    return lhs - rhs


def certify_II():
    lo, hi = mpf('1e10'), mpf('1e80')
    for _ in range(600):
        mid = (lo * hi) ** mpf('0.5')
        if slack_II(mid) > 0:
            hi = mid
        else:
            lo = mid
    m, e = round_up_3sf(hi)
    N = m * mpf(10) ** e
    # C_3 must satisfy C_3 (loglog n)^2 >= z(n)^2 for every n >= N.  The ratio
    # z(n)^2/(loglog n)^2 decreases to A^2, so its supremum is at n = N.
    C3 = round_up_3sf(z_var_up(N) ** 2 / mlog(mlog(N)) ** 2 * UP)
    C3 = C3[0] * mpf(10) ** C3[1]
    checks = {
        'positivity holds at N':       slack_II(N) > 0,
        'z(N) <= N^{3/8}':             z_var_up(N) <= N ** mpf('0.375'),
        'z(N) >= 1':                   z_var_up(N) >= 1,
        'C_3 (loglog n)^2 >= z(n)^2':  all(
            C3 * mlog(mlog(N * mpf(10) ** j)) ** 2 >= z_var_up(N * mpf(10) ** j) ** 2
            for j in range(0, 1000, 5)),
        'range exceeds sqrt(n) at N':  N / (C3 * mlog(mlog(N)) ** 2) > N ** mpf('0.5'),
        'surplus exceeds log 2':       slack_II(N) > mlog(2),
    }
    return m, e, N, C3, checks


# ============================================================ main ==========
if __name__ == '__main__':
    failures = []

    print("Proposition 5.2(i)   thresholds n_r, and comparison with Table 1\n")
    print(f"{'r':>3} {'kappa_r':>12} {'z':>9} {'n_r':>12} {'Cor. 1.2':>11} "
          f"{'gain':>11} {'range/sqrt(n)':>14} {'ok':>4}")
    for r in range(3, RMAX + 1):
        m, e, n_r, checks = certify_I(r)
        bad = [k for k, ok in checks.items() if not ok]
        failures += [f"Prop 5.2(i), r={r}: {k}" for k in bad]
        if r in PROP_52_I and (float(m), e) != PROP_52_I[r]:
            failures.append(f"Prop 5.2(i), r={r}: certified n_r is {float(m)}e{e}, "
                            f"paper prints {PROP_52_I[r][0]}e{PROP_52_I[r][1]}")
        t1 = mpf(TABLE_1[r][2]) * mpf(10) ** TABLE_1[r][3]
        print(f"{r:>3} {float(kappa_lo(r)):>12.3e} {float(z_const_up(r)):>9.1f} "
              f"{float(m):>6.2f}e{e:<5} {float(t1):>11.3g} {float(t1 / n_r):>11.2f} "
              f"{float(kappa_lo(r) * n_r / n_r ** mpf('0.5')):>14.3g} "
              f"{'yes' if not bad else 'NO':>4}")

    m, e, N, C3, checks = certify_II()
    bad = [k for k, ok in checks.items() if not ok]
    failures += [f"Prop 5.2(ii): {k}" for k in bad]
    if (float(m), e) != PROP_52_II_N:
        failures.append(f"Prop 5.2(ii): certified N is {float(m)}e{e}, "
                        f"paper prints {PROP_52_II_N[0]}e{PROP_52_II_N[1]}")
    if C3 != PROP_52_II_C3:
        failures.append(f"Prop 5.2(ii): certified C_3 is {float(C3):g}, "
                        f"paper prints {float(PROP_52_II_C3):g}")
    C2 = z_var_up(N) / mlog(mlog(N))
    if not (PROP_52_II_C2 >= C2 and PROP_52_II_C2 ** 2 <= PROP_52_II_C3):
        failures.append(f"Prop 5.2(ii): printed C_2' = {float(PROP_52_II_C2):g} must satisfy "
                        f"{float(C2):.3f} <= C_2' and C_2'^2 <= {float(PROP_52_II_C3):g}")

    print(f"\nProposition 5.2(ii)  uniform in r\n")
    print(f"  N      = {float(m):.2f}e{e}")
    print(f"  C_3    = {float(C3):.3g}")
    print(f"  z(N)   = {float(z_var_up(N)):.2f}   (y(k) <= {float(y_up(N)):.1f})")
    print(f"  k-range at N is {float(N / (C3 * mlog(mlog(N)) ** 2) / N ** mpf('0.5')):.3g} "
          f"times sqrt(N)")
    for name, ok in checks.items():
        print(f"  {'ok  ' if ok else 'FAIL'}  {name}")

    if failures:
        print("\nCERTIFICATE FAILED:")
        for f in failures:
            print("  ! " + f)
        raise SystemExit(1)
    print("\nAll constants of Proposition 5.2 certified.")
    raise SystemExit(0)
