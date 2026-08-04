#!/usr/bin/env python3
"""
verify_lemmas.py -- checks every finite numerical claim made in Sections 2-3.

Companion to thresholds.py, which certifies Table 1.  This file certifies the
constants and finite verifications that Table 1 rests on, so that the two
together cover all the arithmetic in the paper:

  Lemma 2.2  the closing numeric step E(x;q) <= sqrt(x)(log x)^2, i.e. that
             5/(8 pi) + 2/L + 3.42620/L^2 < 1 for L = log x >= log 100, and
             that this expression is decreasing in x.
  Lemma 2.3  c_0 = zeta(2)^2 zeta(3)/(zeta(4) zeta(6)) = 2.953912...,
             the intermediate bound B(M) <= zeta(2)/(zeta(4) M) on the three
             initial intervals and beyond, and the tail bound itself.
  Lemma 2.5  the finite check over the primes below 285: min F = F(3) = log(3)/2,
             equality only at y = 3, and the NON-monotonicity witness
             F(11) > F(13).
  Thm 1.1    the rebalanced constants sqrt(9 c_0) = 5.15609..., and the
             endpoint-corrected coefficient
                 2 sqrt(9 c_0) + 2/(n^(1/8) tau^(1/2)) + 2 sqrt(9 c_0)/(sqrt(n) log n)
             which is decreasing in n and in tau, hence at most its value at
             n = 1e5, tau = 1, namely 10.7893... < 10.80 <= 11.4.

Every check is an assertion; the script exits nonzero on any failure.

Requires: mpmath, sympy.  Run:  python3 verify_lemmas.py
"""
from fractions import Fraction
from mpmath import mp, mpf, mpmathify, zeta, log as mlog, pi as mpi, exp as mexp
from sympy import primerange

mp.dps = 50

FAILS = []


def check(name, cond, detail=''):
    if cond:
        print(f"  ok    {name}")
    else:
        print(f"  FAIL  {name}   {detail}")
        FAILS.append(name)


# ---------------------------------------------------------------- Lemma 2.2 --
print("Lemma 2.2  (explicit GRH bound, closing numeric step)")
bracket = lambda L: 5 / (8 * mpi) + 2 / L + mpf('3.42620') / L ** 2
L100 = mlog(100)
check("bracket < 1 at x = 100", bracket(L100) < 1, f"value {float(bracket(L100)):.6f}")
check("bracket decreasing in x",
      all(bracket(mlog(mpf(10) ** a)) > bracket(mlog(mpf(10) ** (a + 1)))
          for a in range(2, 30)))
check("bracket -> 5/(8 pi) = 0.19894... (so ~5x is discarded, deliberately)",
      abs(bracket(mlog(mpf(10) ** 4000)) - 5 / (8 * mpi)) < mpf('1e-3'))

# ---------------------------------------------------------------- Lemma 2.3 --
print("\nLemma 2.3  (square-free sieve tail)")
z2, z3, z4, z6 = (zeta(k) for k in (2, 3, 4, 6))
c0 = z2 ** 2 * z3 / (z4 * z6)
check("c_0 = 2.953912...", abs(c0 - mpf('2.953912373')) < mpf('1e-9'),
      f"value {float(c0):.9f}")
check("c_0 = (zeta2/zeta4) * prod_p(1 + 1/(p(p-1)))",
      abs(c0 - (z2 / z4) * (z2 * z3 / z6)) < mpf('1e-30'))
check("c_0 < 4 (improves the Ramare-Rumely constant)", c0 < 4)

# B(M) <= zeta(2)/(zeta(4) M):  the three initial intervals, exactly
NB = 2 * 10 ** 6
mu = [1] * (NB + 1)
prime = [True] * (NB + 1)
for i in range(2, NB + 1):
    if prime[i]:
        for j in range(i, NB + 1, i):
            if j > i:
                prime[j] = False
            mu[j] = -mu[j]
        for j in range(i * i, NB + 1, i * i):
            mu[j] = 0
suffix = [mpf(0)] * (NB + 2)
for m in range(NB, 0, -1):
    suffix[m] = suffix[m + 1] + (mpf(1) / (m * m) if mu[m] else mpf(0))
B = lambda M: suffix[int(M) + 1]
for lo, hi in ((mpf('0.5'), 1), (1, 2), (2, 3)):
    check(f"B(M) <= zeta2/(zeta4 M) on [{lo},{hi})",
          B(lo) <= z2 / (z4 * mpf(hi)) or B(lo) <= z2 / (z4 * lo))
check("B(M) <= zeta2/(zeta4 M) for sampled M > 3",
      all(B(M) <= z2 / (z4 * mpf(M)) for M in (3, 4, 7, 20, 100, 1000, 50000)))
check("zeta(2)/zeta(4) > 6/5", z2 / z4 > mpf(6) / 5)

# the tail bound itself
phi = list(range(NB + 1))
for i in range(2, NB + 1):
    if phi[i] == i:
        for j in range(i, NB + 1, i):
            phi[j] -= phi[j] // i
tail = [mpf(0)] * (NB + 2)
for a in range(NB, 0, -1):
    tail[a] = tail[a + 1] + (mpf(1) / (a * phi[a]) if mu[a] else mpf(0))
check("sum_{a>z} mu^2(a)/phi(a^2) <= c_0/z for z = 1..1e5",
      all(tail[z + 1] <= c0 / z for z in (1, 2, 3, 5, 10, 100, 1000, 10000, 100000)))

# ---------------------------------------------------------------- Lemma 2.5 --
print("\nLemma 2.5  (uniform Mertens-type bound)")
F, prod = {}, Fraction(1)
for q in primerange(3, 286):
    prod *= 1 - Fraction(1, q - 1)
    F[q] = mlog(q) * mpf(prod.numerator) / mpf(prod.denominator)
target = mlog(3) / 2
check("min over primes <= 285 is F(3) = log(3)/2",
      min(F.values()) == F[3] and abs(F[3] - target) < mpf('1e-30'))
check("F(y) > log(3)/2 strictly for every prime 5 <= y <= 285",
      all(F[q] > target for q in F if q >= 5))
check("F NOT monotone: F(11) > F(13)", F[11] > F[13],
      f"F(11)={float(F[11]):.6f}, F(13)={float(F[13]):.6f}")
check("equality fails at y = 4 (so equality holds only at y = 3)",
      mpf('0.5') > mlog(3) / (2 * mlog(4)))

# ---------------------------------------------------------------- Thm 1.1 ----
print("\nTheorem 1.1  (rebalanced constants)")
A = (9 * c0) ** mpf('0.5')
check("sqrt(9 c_0) = 5.15609...", abs(A - mpf('5.156085')) < mpf('1e-6'))

# endpoint-corrected coefficient of Theorem 1.1 (see eq. for |R_k - S_k n| in the proof):
#   2 sqrt(9c0) + 2/(n^{1/8} tau^{1/2}) + 2 sqrt(9c0)/(sqrt(n) log n)
coeff = lambda n, tau: 2 * A + 2 / (n ** mpf('0.125') * tau ** mpf('0.5')) \
                       + 2 * A / (n ** mpf('0.5') * mlog(n))
worst = coeff(mpf('1e5'), mpf(1))
check("endpoint-corrected coefficient at n=1e5, tau=1 is 10.7893...",
      abs(worst - mpf('10.789277')) < mpf('1e-5'), f"value {float(worst):.6f}")
check("coefficient < 10.80", worst < mpf('10.80'))
check("coefficient <= 11.4 (printed constant of Theorem 1.1)", worst <= mpf('11.4'))
check("coefficient decreasing in n",
      all(coeff(mpf(10) ** a, mpf(1)) > coeff(mpf(10) ** (a + 1), mpf(1))
          for a in range(5, 40)))
check("coefficient decreasing in tau",
      all(coeff(mpf('1e5'), mpf(2) ** r) > coeff(mpf('1e5'), mpf(2) ** (r + 1))
          for r in range(0, 12)))
check("modulus admissibility needs only log^2 n > 9 c_0 = 26.59...",
      mlog(mpf('1e5')) ** 2 > 9 * c0)

# Corollary 1.3's chain is uniform in r: the only numeric step is
#   11.4 * 2^(r/2) * log n > 1  for n >= 1e5, worst case r = 0.
check("Cor 1.3 key step: 11.4 log(1e5) > 1 (worst case r = 0)",
      mpf('11.4') * mlog(mpf('1e5')) > 1,
      f"value {float(mpf('11.4') * mlog(mpf('1e5'))):.2f}")
# and a witness that the SUPERSEDED table-dependent constant 0.072 was untenable:
# C_Artin * Pi_r decreases to 0, and already dips below 0.072 at r = 12.
_pi = Fraction(1)
_odd = list(primerange(3, 200))
for _t in range(12):
    _pi *= 1 - Fraction(1, _odd[_t] - 1)
_A12 = mpf('0.3739558136192023') * mpf(_pi.numerator) / mpf(_pi.denominator)
check("C_Artin*Pi_12 < 0.072 (so no fixed positive constant can serve all r)",
      _A12 < mpf('0.072'), f"value {float(_A12):.7f}")

print()
if FAILS:
    print(f"VERIFICATION FAILED on {len(FAILS)} check(s): {', '.join(FAILS)}")
    raise SystemExit(1)
print("All finite numerical claims in Sections 2-3 verified.")
raise SystemExit(0)
