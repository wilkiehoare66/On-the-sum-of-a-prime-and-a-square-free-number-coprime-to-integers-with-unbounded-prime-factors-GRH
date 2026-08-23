# On the sum of a prime and a square-free number coprime to integers with arbitrarily many prime factors, under GRH

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21828256.svg)](https://doi.org/10.5281/zenodo.21828255)

Certificate scripts for the numerical claims in the paper *"On the sum of a
prime and a square-free number coprime to integers with arbitrarily many prime
factors, under GRH"* by W. Hoare, E. S. Lee, and A. Pearce-Crump.

**What this repository is.** Every theorem in the paper is proved analytically;
no result depends on an exhaustive search over a range of `n`. What the paper
does contain is a set of finite numerical claims — the checks that close several
of the lemmas, the ten rows of Table 1, and the constants of `prop:largek` — and
these three scripts certify all of them. Each is assertion-driven and exits
nonzero if any claim fails, so a reader can check every number in the paper with
three commands, in under a minute.

The paper is hosted separately; this repository holds only the certificates.

---

## Requirements

Python 3.8+, with `mpmath` and `sympy`:

```bash
pip install mpmath sympy
```

## Repository layout

```
.
├── thresholds.py      # Table 1 rows 3-9, and the uniform threshold N_0
├── largek.py          # constants of prop:largek: n_r (Table 1 rows 10-12), C_2', C_3, N_1
├── verify_lemmas.py   # finite checks in lem:grh, lem:ramare, lem:mertens,
│                      #   thm:main's assembled constant, eq:taubound, rmk:notmin
└── README.md
```

`largek.py` imports `ODD`, `pi_r`, `k_min`, `CARTIN_LO`, `TABLE_1` and `RMAX`
from `thresholds.py`, so the two cannot disagree about the shared inputs.
`verify_lemmas.py` is independent of both.

## Running

```bash
python3 thresholds.py      # ~1 s
python3 largek.py          # ~1 s
python3 verify_lemmas.py   # ~50 s
echo $?                    # 0 = all claims certified, 1 = at least one failed
```

## A note on cross-references

The scripts refer to the paper's results by **LaTeX label** (`thm:main`,
`cor:positivity`, `eq:master`, ...) rather than by number. Numbers move when the
paper is reorganised; labels do not. To find a result named here, grep the
paper's source for its label.

---

## `thresholds.py` — Table 1 rows 3-9, and `N_0`

`cor:positivity` states that `R_k(n) > 0` — hence that `n` is a sum of a prime
and a square-free number coprime to `k` — as soon as

```
n^(1/4)  >  (6.5 / (C_Artin * Pi_r)) * 2^(r/2) * log n,      r = omega(k),
```

where `Pi_r = prod_{t<=r} (1 - 1/(q_t - 1))` runs over the smallest `r` odd
primes. The criterion involves `k` only through `r`, so it tabulates one row per
`r`. Rows `3 <= r <= 9` of Table 1 are these thresholds; from `r = 10` the table
takes the smaller values supplied by `prop:largek(i)`, which `largek.py`
certifies.

The script re-derives every entry and asserts it: the exact value of `Pi_r`, the
smallest admissible `k`, and the threshold's mantissa and exponent. It solves
against `cor:even`'s strengthened `+1` criterion rather than the plain one,
because at `C = 6.5` the plain crossing rounded up to three significant figures
no longer leaves room for the `+1` at `r = 4` and `r = 5`; solving the weaker
form and asserting the stronger one afterwards would fail there. It also
certifies `N_0 = 5.22e19` for `cor:uniform`: that the least admissible `n`
increases with `r` and is at most `N_0` for `r <= 9`, and that the cascade holds
for `r >= 10`, both through its base case and growth factor and by direct
evaluation for `10 <= r <= 60`.

**Directed rounding.** A floating-point scan is useful for locating a threshold
but is not by itself a certificate, so every step is directed to err on the safe
side:

- `Pi_r` is held as an exact `Fraction`, never as a float, and so is `C`;
- Artin's constant is replaced by a strict **lower** bound, which over-estimates
  `K_r` and so makes the criterion harder;
- `K_r` and `n_0` are both rounded **up**;
- the criterion is re-verified at the **rounded** `n_0`, with the left-hand side
  rounded down and the right-hand side rounded up.

Since `n^(1/4)/log n` is increasing for `n > e^4`, verifying the inequality at
`n_0` certifies it for every `n >= n_0`.

All inputs are documented in the script's header with provenance. The constant
`6.5` is `thm:main`'s error constant, the round-up of the `6.44380...` its proof
assembles; the lower bound for Artin's constant is from Wrench, *Math. Comp.*
**15** (1961), 396-398.

The script's `TABLE_1` dictionary holds the `cor:positivity` crossing at **every**
`r`, including `r >= 10` where the paper prints `prop:largek(i)`'s smaller value
instead. Those rows are retained because `cor:uniform`'s cascade and
`largek.py`'s comparison column both need the `cor:positivity` value at every
`r`. Do not read rows 10-12 of that dictionary as Table 1 entries.

## `largek.py` — the constants of `prop:largek`

`prop:largek` gives up the asymptotic and keeps only positivity, which allows a
far smaller truncation level `z` and hence a much larger range of `k`. Part (i),
at fixed `r`, takes the constant level `z = 30 c_0 / delta_r` and admits
`k <= kappa_r n`; part (ii) takes `z` proportional to `log y(k)` and admits
`k <= n / (C_3 (log log n)^2)`, uniformly in `r`.

The script computes what the paper leaves to it: the thresholds `n_r` of part
(i), which supply rows 10-12 of Table 1, and the constants `C_2' = 448.9`,
`C_3 = 2.02e5` and `N_1 = 1.65e39` of part (ii). It derives them in the same
order the proof does — `N_1` first, from an inequality in `n` alone, then `C_2'`
and `C_3` from `N_1` — so that nothing is used before it exists.

At each threshold it checks not only the governing inequality but every
hypothesis inherited from `eq:master`: that `n >= 1e5`, that `1 <= z <= n^(3/10)`,
that an admissible `k` exists at all, that the claimed range genuinely exceeds
`sqrt(n)`, and that the surplus exceeds `log 2`, so that the term `p = 2` may be
discarded as in `cor:even` and an even-`k` analogue is available at the same row.

## `verify_lemmas.py` — the finite checks

Several of the paper's proofs close on a finite numerical step. These are small
enough to be done by hand, but they are part of the proofs, so they are checked
here explicitly:

| Claim | What is verified |
| --- | --- |
| `lem:grh` | that `beta = 0.39854` dominates the bracket `5/(8 pi) + 2/L + 3.43/L^2` at `x = 1e5` (`L = log x`), that the bracket decreases in `x` so `beta` serves every larger `x`, and the `q in {1,2}` carve-out through `psi`.  The constant `3.43` is that of Lee, Q. J. Math. **74** (2023), App. A, (A7), which `lem:grh` now quotes in place of specialising Grenie-Molteni by hand |
| `lem:ramare` | `c_0 = zeta(2)^2 zeta(3)/(zeta(4) zeta(6)) = 2.953912...`; the intermediate bound `B(M) <= zeta(2)/(zeta(4) M)` on the three initial intervals and beyond; `zeta(2)/zeta(4) > 6/5`; and the tail bound `sum_{a>z} mu^2(a)/phi(a^2) <= c_0/z` itself |
| `lem:mertens` | the finite check over the primes below 285: that `min F = F(3) = log(3)/2`, that `F(l) > log(3)/2` strictly at every prime `5 <= l <= 285`, that equality therefore holds only at `y = 3`, and the non-monotonicity witness `F(11) = 0.674408... > F(13) = 0.661276...` of `rmk:mertens` |
| `thm:main` | the rebalanced constant `sqrt(6 c_0 beta) = 2.65772...` and the assembled coefficient `6.44380... <= 6.5`; that the coefficient decreases in `n` and in `tau(k)`; and that modulus admissibility needs only `log^2 n > 6 c_0 / beta = 44.48...` |
| `lem:tail` | the split point `lambda = 3/10`: that `A(lambda) = 1 + 2/(1-2 lambda) = 6` exactly, that the closing inequality `eq:lamcheck` holds at `n = 1e5` and decreases thereafter, and the optimality claims of `rmk:splitpoint` (minimum `6.43640...` at `lambda = 0.308...`; every `0.29 <= lambda <= 0.33` delivers at most `6.5`; the superseded cut `lambda = 3/8` delivers `6.98884...`) |
| `eq:taubound` | `sup_m tau(m) m^(-1/4) = 8.44696... < 8.5` as an exact supremum over **all** `m`, via the local maxima of the multiplicative function, not over a finite range |
| `rmk:notmin` | that `Sing_k(n_Y) > 2/5` along the odd primorials for every prime `Y >= 191`, uniformly in `k`, and that `Y = 191` is the least prime at which it does |

The `F(11) > F(13)` check is worth singling out: it is what rules out the
tempting (and false) claim that `F` is monotone, and hence why the proof of
`lem:mertens` reduces the range `3 <= y <= 285` to its primes by monotonicity *on
each prime gap* rather than by any global monotonicity.

---

## Notes for reviewers

- The results of record are what the scripts print; the exit status is the
  verdict. No script has a "check" mode distinct from its normal run.
- `thresholds.py` and `largek.py` each hold the paper's printed values in an
  explicit dictionary and assert against them, so the repository fails if the
  paper and the code ever drift apart.
- Every constant in the paper is explicit and asserted by one of these scripts.
  Nothing is left in the effective-but-uncomputed form that a bound of this kind
  usually carries.
- `lem:grh` is not sharp: its bracket is about `0.4` at `x = 1e5` and decreases
  to `5/(8 pi) = 0.19894...`, so carrying it at `x = n` rather than at `x = 1e5`
  would lower every row of Table 1 by a further factor of about `3.7`. The paper
  declines this, for the reason given in `rmk:elementary`: the sharp bracket
  depends on `x`, so `thm:main` would carry no numerical constant at all.
