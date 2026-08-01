# On the sum of a prime and a square-free number coprime to integers with arbitrarily many prime factors, under GRH

Certificate scripts for the numerical claims in the paper

> **On the sum of a prime and a square-free number coprime to integers with
> arbitrarily many prime factors, under GRH**
> W. Hoare, E. S. Lee, and A. Pearce-Crump.

**What this repository is, and is not.** Every theorem in the paper is proved
analytically; unlike the companion work on `omega(k) <= 3`, there is **no
exhaustive computational range** on which any result depends, and so no C++ here.
What the paper does contain is a handful of finite numerical claims — three small
finite checks inside the lemmas, and the ten rows of Table 1 — and these two
scripts certify all of them. Both are assertion-driven and exit nonzero if any
claim fails, so a reader can check every number in the paper with two commands.

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
├── thresholds.py      # certifies Table 1 (Corollary 1.2 thresholds)
├── verify_lemmas.py   # certifies the finite checks in Lemmas 2.2, 2.3, 2.5
└── README.md
```

## Running

```bash
python3 thresholds.py      # ~2 s
python3 verify_lemmas.py   # ~35 s
echo $?                    # 0 = all claims certified, 1 = at least one failed
```

---

## `thresholds.py` — Table 1

Corollary 1.2 states that `R_k(n) > 0` — hence that `n` is a sum of a prime and
a square-free number coprime to `k` — as soon as

```
n^(1/4)  >  (11.4 / (C_Artin * Pi_r)) * 2^(r/2) * log n,      r = omega(k),
```

where `Pi_r = prod_{t<=r} (1 - 1/(q_t - 1))` runs over the smallest `r` odd
primes. The criterion depends on `k` only through `r`, so it tabulates one row
per `r`; Table 1 of the paper records the resulting thresholds `n_0` for
`3 <= r <= 12`.

The script re-derives every printed entry and asserts it. Each row is checked
for: the exact value of `Pi_r`, the smallest admissible modulus, the certified
threshold's mantissa and exponent, and — additionally — that the printed `n_0`
also clears the strengthened `+1` criterion used for even moduli in
Corollary 1.3. Any mismatch between the code and the paper fails loudly.

**Directed rounding.** A floating-point scan is useful for locating a threshold
but is not by itself a certificate, so every step is directed to err on the safe
side:

- `Pi_r` is held as an exact `Fraction`, never as a float;
- Artin's constant is replaced by a strict **lower** bound, which over-estimates
  `K_r` and so makes the criterion harder;
- `K_r` and `n_0` are both rounded **up**;
- the criterion is re-verified at the **rounded** `n_0`, with the left-hand side
  rounded down and the right-hand side rounded up;
- `n_0` is raised to `(prod_{i<=r} q_i)^2` on the rows where the range
  `k <= sqrt(n)` would otherwise make the statement vacuous (this is what binds
  from `r = 10` onwards, and the script reports which constraint fixes each row).

Since `n^(1/4)/log n` is increasing for `n > e^4`, verifying the inequality at
`n_0` certifies it for every `n >= n_0`.

All inputs are documented in the header of the script, with provenance: the
constant `11.4` is Theorem 1.1's error constant (the proof gives
`2*sqrt(9*c_0) + 1 = 11.3122...`), and the lower bound for Artin's constant is
from Wrench, *Math. Comp.* **15** (1961), 396–398.

## `verify_lemmas.py` — the finite checks

Three of the paper's lemmas close on a finite numerical step. These are small
enough to be done by hand, but they are part of the proofs, so they are checked
here explicitly:

| Claim | What is verified |
| --- | --- |
| Lemma 2.2 | `5/(8 pi) + 2/L + 3.42620/L^2 < 1` at `x = 100` (`L = log x`), that it is decreasing in `x`, and that it tends to `5/(8 pi) = 0.19894...` — so the printed constant `1` deliberately discards about a factor of five |
| Lemma 2.3 | `c_0 = zeta(2)^2 zeta(3)/(zeta(4) zeta(6)) = 2.953912...`; the intermediate bound `B(M) <= zeta(2)/(zeta(4) M)` on the three initial intervals and beyond; `zeta(2)/zeta(4) > 6/5`; and the tail bound `sum_{a>z} mu^2(a)/phi(a^2) <= c_0/z` itself |
| Lemma 2.5 | the finite check over the primes below 285: that `min F = F(3) = log(3)/2`, that `F(l) > log(3)/2` strictly at every prime `5 <= l <= 285`, that equality therefore holds only at `y = 3`, and the non-monotonicity witness `F(11) = 0.674408... > F(13) = 0.661276...` |
| Theorem 1.1 | the rebalanced constants `sqrt(9 c_0) = 5.15609...` and `2 sqrt(9 c_0) + 1 = 11.3122... <= 11.4`, and that modulus admissibility needs only `log^2 n > 9 c_0 = 26.59...` |

The `F(11) > F(13)` check is worth singling out: it is what rules out the
tempting (and false) claim that `F` is monotone, and hence why the proof of
Lemma 2.5 reduces the range `3 <= y <= 285` to its primes by monotonicity *on
each prime gap* rather than by any global monotonicity.

---

## Notes for reviewers

- Both scripts are self-contained and import nothing from each other. Together
  they cover every numerical assertion in the paper.
- The results of record are what the scripts print; the exit status is the
  verdict. Neither script has a "check" mode distinct from its normal run.
- `thresholds.py` holds the paper's printed Table 1 in an explicit `TABLE_1`
  dictionary and asserts against it, so the repository fails if the paper and
  the code ever drift apart.
- Table 1's rows for `r >= 10` are fixed by the range `k <= sqrt(n)` rather than
  by the error term. Proposition 7.4 of the paper shows that positivity in fact
  persists for `k` almost as large as `n` itself, but its constants, though
  effective, have not been computed; no entry of Table 1 relies on it.
