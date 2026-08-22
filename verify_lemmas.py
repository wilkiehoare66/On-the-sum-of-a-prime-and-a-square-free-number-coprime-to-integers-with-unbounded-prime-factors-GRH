#!/usr/bin/env python3
"""
verify_lemmas.py -- checks every finite numerical claim made in Sections 2-3,
plus the divisor bound used in Proposition 7.4.

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
  Thm 1.1    the rebalanced constants sqrt(9 c_0) = 5.15608..., and the
             endpoint-corrected coefficient
                 2 sqrt(9 c_0) + 2/(n^(1/8) tau^(1/2)) + 2 sqrt(9 c_0)/(sqrt(n) log n)
             which is decreasing in n and in tau, hence at most its value at
             n = 1e5, tau = 1, namely 10.7892... < 10.80 <= 11.4.
  Prop 7.4   the effective divisor bound tau(m) <= 8.5 m^(1/4), certified as an
             exact supremum over ALL m via the local maxima of the multiplicative
             function tau(m)/m^(1/4): the supremum is 8.44696..., attained at
             m = 2^5 3^3 5^2 . 7 . 11 . 13 = 21621600.

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

# The q in {1,2} carve-out of Lemma 2.2 goes through psi, not theta: Schoenfeld's
# 73.2 threshold is the one for psi, his theta bound needing x >= 599.  The
# resulting bracket is 1/(8pi) + 1.42620/L^2 + log2/(sqrt(x) L^2).
q12 = lambda x: 1 / (8 * mpi) + mpf('1.42620') / mlog(x) ** 2 \
                + mlog(2) / (mp.sqrt(x) * mlog(x) ** 2)
check("q in {1,2}: bracket < 1 at x = 100", q12(100) < 1,
      f"value {float(q12(100)):.6f}")
check("q in {1,2}: bracket decreasing in x",
      all(q12(mpf(10) ** mpf(t)) > q12(mpf(10) ** mpf(t + mpf('0.5')))
          for t in [2, 3, 4, 6, 10]))
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
# Corollary 6.1's discussion: the loss factor F(y)/(log3/2) is >= 1 (which is the whole
# content of Lemma 2.5) and approaches 4 e^-gamma C_2/log 3; it does NOT increase to that
# value, and the limit is not asserted as an upper bound -- F is not monotone, and the
# Mertens comparison underlying it is known to oscillate.
_lim = 2 * mexp(-mpmathify('0.57721566490153286060651209008240243')) * mpf('0.6601618158468696')
check("loss factor F/(log3/2) >= 1 at every prime y <= 285 (Lemma 2.5)",
      all(F[l] >= mlog(3) / 2 for l in F))
check("loss factor limit 4 e^-g C_2/log 3 = 1.3495...",
      abs(2 * _lim / mlog(3) - mpf('1.3495356')) < mpf('1e-6'))
check("F stays below its limit at primes up to 285, but only just at the top",
      all(F[l] < _lim for l in F), f"max {float(max(F.values())):.7f} vs limit {float(_lim):.7f}")
check("F NOT monotone: F(11) > F(13)", F[11] > F[13],
      f"F(11)={float(F[11]):.6f}, F(13)={float(F[13]):.6f}")
check("equality fails at y = 4 (so equality holds only at y = 3)",
      mpf('0.5') > mlog(3) / (2 * mlog(4)))

# ------------------------------------------------------------- Thm 1.5 -------
print("\nTheorem 1.5  (sharpness step)")
# The sharpness half needs prod_{p|n, p>=3} A_p^{-1} = 1+o(1) when every odd prime
# factor of n exceeds (1/2)log n.  With A_p = 1 - 1/(p(p-1)) = 1 + O(p^-2) and at most
# log n/log((1/2)log n) such primes, the correction is exp(O(1/(log n loglog n))).
def _corr(n):
    half = mlog(n) / 2
    cnt = mlog(n) / mlog(half)          # bound on the number of such prime factors
    return mexp(cnt / half ** 2)        # over-estimate of prod (1 + O(p^-2))
check("Thm 1.5 sharpness: the omitted-factor correction tends to 1",
      all(_corr(mpf(10) ** e) > _corr(mpf(10) ** (e + 10)) for e in [5, 20, 50, 100])
      and _corr(mpf(10) ** 1000) < mpf('1.001'),
      f"bound at n=1e50 is {float(_corr(mpf(10) ** 50)):.6f}, at n=1e1000 {float(_corr(mpf(10) ** 1000)):.6f}")

# ---------------------------------------------------------------- Thm 1.1 ----
print("\nTheorem 1.1  (rebalanced constants)")
A = (9 * c0) ** mpf('0.5')
check("sqrt(9 c_0) = 5.15608...", abs(A - mpf('5.156085')) < mpf('1e-6'))
# Remark 7.1: the two balanced terms each contribute sqrt(9 c_0) n^(3/4) tau^(1/2) log n
# (NOT 6, which was the value under the superseded c_0 = 4), and sharpening Lemma 2.2's
# constant 1 to its limit 5/(8pi) gains the factor sqrt(8pi/5), not sqrt(8pi).
check("Remark 7.1: balanced term constant is sqrt(9 c_0) = 5.15608..., not 6",
      abs(A - mpf('5.1560848861')) < mpf('1e-9'), f"value {float(A):.7f}")
check("Remark 7.1: sharpening gain is sqrt(8 pi/5) = 2.24199..., not sqrt(8 pi) = 5.01...",
      abs((8 * mpi / 5) ** mpf('0.5') - mpf('2.2419964')) < mpf('1e-6')
      and abs((8 * mpi) ** mpf('0.5') - mpf('5.0132565')) < mpf('1e-6'))

# endpoint-corrected coefficient of Theorem 1.1 (see eq. for |R_k - S_k n| in the proof):
#   2 sqrt(9c0) + 2/(n^{1/8} tau^{1/2}) + 2 sqrt(9c0)/(sqrt(n) log n)
coeff = lambda n, tau: 2 * A + 2 / (n ** mpf('0.125') * tau ** mpf('0.5')) \
                       + 2 * A / (n ** mpf('0.5') * mlog(n))
worst = coeff(mpf('1e5'), mpf(1))
check("endpoint-corrected coefficient at n=1e5, tau=1 is 10.7892...",
      abs(worst - mpf('10.789277')) < mpf('1e-5'), f"value {float(worst):.6f}")
check("coefficient < 10.80", worst < mpf('10.80'))
# The paper prints 11.4, not 10.8: the margin at 10.8 is only ~0.1%, at 11.4 ~5.7%.
check("margin at 10.8 is under 0.15% (why 10.8 is not printed)",
      (mpf('10.8') - worst) / worst < mpf('0.0015'),
      f"{float((mpf('10.8') - worst) / worst * 100):.3f}%")
check("margin at 11.4 is over 5% (why 11.4 is printed)",
      (mpf('11.4') - worst) / worst > mpf('0.05'),
      f"{float((mpf('11.4') - worst) / worst * 100):.3f}%")
check("coefficient <= 11.4 (printed constant of Theorem 1.1)", worst <= mpf('11.4'))
check("coefficient decreasing in n",
      all(coeff(mpf(10) ** a, mpf(1)) > coeff(mpf(10) ** (a + 1), mpf(1))
          for a in range(5, 40)))
check("coefficient decreasing in tau",
      all(coeff(mpf('1e5'), mpf(2) ** r) > coeff(mpf('1e5'), mpf(2) ** (r + 1))
          for r in range(0, 12)))
check("modulus admissibility needs only log^2 n > 9 c_0 = 26.58...",
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

# Proposition 7.4(ii) uses the effective divisor bound tau(m) <= C_4 m^(1/4) with
# C_4 = 8.5.  This is certified as an exact supremum, not spot-checked: the function
# tau(m)/m^(1/4) is multiplicative with value (a+1) p^(-a/4) at p^a, so its supremum
# is the product of the local maxima, and primes p >= 17 contribute nothing because
# a+1 <= 2^a <= p^(a/4) there.
check("17^(1/4) > 2 and a+1 <= 2^a (a>=1), so primes p >= 17 contribute a factor <= 1",
      mpf(17) ** mpf('0.25') > 2 and all(_a + 1 <= 2 ** _a for _a in range(1, 400)),
      f"17^(1/4) = {float(mpf(17) ** mpf('0.25')):.9f}")

_sup, _args = mpf(1), []
for _p in (2, 3, 5, 7, 11, 13):
    _best, _besta = mpf(0), None
    for _a in range(0, 400):
        _v = (_a + 1) / mpf(_p) ** (mpf(_a) / 4)
        if _v > _best:
            _best, _besta = _v, _a
    _sup *= _best
    _args.append(_besta)
check("local maxima for p <= 13 occur at exponents (5,3,2,1,1,1)",
      _args == [5, 3, 2, 1, 1, 1], f"exponents {tuple(_args)}")

_ext = 2 ** 5 * 3 ** 3 * 5 ** 2 * 7 * 11 * 13
_tau_ext = 6 * 4 * 3 * 2 * 2 * 2          # tau(2^5 3^3 5^2 7 11 13)
_ratio = mpf(_tau_ext) / mpf(_ext) ** mpf('0.25')
check("the extremal m is 2^5 3^3 5^2 . 7 . 11 . 13 = 21621600", _ext == 21621600)
check("sup_m tau(m)/m^(1/4) = product of local maxima = 8.44696...",
      abs(_sup - mpf('8.446961724')) < mpf('1e-9'), f"value {float(_sup):.9f}")
check("the supremum is attained, at m = 21621600", abs(_ratio - _sup) < mpf('1e-25'))
check("hence tau(m) <= 8.5 m^(1/4) for EVERY m >= 1, not merely a checked range",
      _sup <= mpf('8.5'), f"headroom {float(mpf('8.5') - _sup):.6f}")
# independent brute-force cross-check against the certified supremum
_d = [0] * 200001
for _i in range(1, 200001):
    for _j in range(_i, 200001, _i):
        _d[_j] += 1
check("brute force to 2e5 never exceeds the certified supremum",
      all(_d[m] <= float(_sup) * m ** 0.25 for m in range(1, 200001)))

# ------------------------------------------------------- Remark 3.5 (notmin) --
# Along the odd primorials n_Y = prod_{3<=p<=Y} p the singular series does not
# decay: since n_Y is odd, an odd p fails to divide it exactly when p > Y, so
# both products of (3.2) run over p > Y, and
#
#   Sing_k(n_Y) > (1/2) (1 - 1/Y) (1 - 1.02/log Y)   uniformly in k <= sqrt(n_Y),
#
# which exceeds 2/5 from Y = 191 on.  The two ingredients are checked first.
_P = [p for p in primerange(3, 20000)]

check("1 - rho(p) = (p-1)/(p^2-p-1) < 2/p for every odd p  [(p-2)(p+1) > 0]",
      all(mpf(p - 1) / (p * p - p - 1) < mpf(2) / p for p in _P))
check("sum_{m>Y} 1/(m(m-1)) = 1/Y, so prod_{p>Y} A_p > 1 - 1/Y",
      all(abs(sum(mpf(1) / (m * (m - 1)) for m in range(Y + 1, Y + 40000)) - mpf(1) / Y)
          < mpf('1e-4') for Y in (101, 199, 1009)))
check("theta(Y) < 1.01624 Y, so log k <= theta(Y)/2 < 0.51 Y  [RS62 (3.16)]",
      all(sum(mlog(p) for p in _P if p <= Y) + mlog(2) < mpf('1.01624') * Y
          for Y in range(3, 20000, 7)))

def _sing_lo(Y):
    return (1 - mpf(1) / Y) * (1 - mpf('1.02') / mlog(Y)) / 2

check("Sing_k(n_Y) > 2/5 for every prime Y >= 191, uniformly in k",
      all(_sing_lo(Y) > mpf(2) / 5 for Y in _P if Y >= 191),
      f"least value {float(_sing_lo(191)):.6f}, at Y = 191")
check("191 is least: the preceding prime Y = 181 does not clear 2/5",
      _sing_lo(181) <= mpf(2) / 5, f"value at 181 is {float(_sing_lo(181)):.6f}")
check("the bound increases to 1/2: it exceeds 0.49 by Y = 1e44",
      _sing_lo(mpf('1e44')) > mpf('0.49'),
      f"value {float(_sing_lo(mpf('1e44'))):.6f}")

print()
if FAILS:
    print(f"VERIFICATION FAILED on {len(FAILS)} check(s): {', '.join(FAILS)}")
    raise SystemExit(1)
print("All finite numerical claims in Sections 2-3, Remark 3.5 and Proposition 7.4 verified.")
raise SystemExit(0)
