#!/usr/bin/env python3
"""Cross-check every constant printed in main.tex against the certificates."""
import re, sys
from mpmath import mp, mpf, log as mlog
mp.dps = 60
tex = open(sys.argv[1] if len(sys.argv) > 1 else 'main.tex').read()

th = {}; exec(compile(open('thresholds.py').read().split("if __name__ == '__main__':")[0],
                      'thresholds.py', 'exec'), th)
lk = {}; exec(compile(open('largek.py').read().split("if __name__ == '__main__':")[0],
                      'largek.py', 'exec'), lk)

bad = []
def want(desc, pattern, expected, count=None):
    hits = re.findall(pattern, tex)
    if not hits:
        bad.append(f"{desc}: pattern not found in main.tex"); return
    if count is not None and len(hits) != count:
        bad.append(f"{desc}: expected {count} occurrences, found {len(hits)}")
    for h in set(hits):
        if h != expected:
            bad.append(f"{desc}: paper prints {h!r}, certificate gives {expected!r}")
    print(f"  ok   {desc:44} {expected}")

# --- Table 1: pull the n_0 cell out of each row and compare with the certificate
rows = dict(re.findall(r"^\$(\d+)\$\s+&.*?&.*?&.*?& \$([0-9.]+\\cdot10\^\{\d+\})\$ &",
                       tex, re.M))
for r in range(3, 13):
    cert = (th['certify'](r)[:2] if r < 10 else lk['certify_I'](r)[:2])
    exp = f"{float(cert[0]):.2f}\\cdot10^{{{cert[1]}}}"
    got = rows.get(str(r))
    if got is None:
        bad.append(f"Table 1 row r={r}: not found")
    elif got != exp:
        bad.append(f"Table 1 row r={r}: paper prints {got}, certificate gives {exp}")
    else:
        print(f"  ok   {'Table 1 row r=' + str(r):44} {exp}")

# --- Table 1 k-range column for r >= 10
kap = dict(re.findall(r"^\$(1[012])\$ &.*?& \$([0-9.]+)\\cdot10\^\{-7\}n\$", tex, re.M))
for r in (10, 11, 12):
    printed = mpf(kap[str(r)]) * mpf('1e-7')
    if printed <= lk['kappa_lo'](r):
        print(f"  ok   {'Table 1 kappa_' + str(r):44} {kap[str(r)]}e-7 <= {float(lk['kappa_lo'](r)):.4e}")
    else:
        bad.append(f"Table 1 kappa_{r}: printed {kap[str(r)]}e-7 exceeds certified {float(lk['kappa_lo'](r)):.6e}")

# --- N_0
N0 = max(th['crossover'](r) for r in range(0, 10))
m, e = th['round_up_3sf'](N0)
want("N_0 (cor:uniform)", r"N_0=(\d\.\d\d\\cdot10\^\{19\})", f"{float(m):.2f}\\cdot10^{{19}}")

# --- prop:largek(ii): N_1, C_3, C_2'
m, e, N, C3, _ = lk['certify_II']()
want("N_1 (prop:largek(ii))", r"(1\.\d\d\\cdot10\^\{39\})", f"{float(m):.2f}\\cdot10^{{39}}")
want("C_3", r"C_3=(\d\.\d\d\\cdot10\^\{5\})", f"{float(C3)/1e5:.2f}\\cdot10^{{5}}")
C2 = lk['z_var_up'](N) / mlog(mlog(N))
printed_C2 = mpf(re.search(r"C_2'=(\d+\.\d)", tex).group(1))
print(f"  {'ok  ' if printed_C2 >= C2 else 'FAIL'} C_2' printed {float(printed_C2)} >= certified {float(C2):.6f}")
if printed_C2 < C2: bad.append("C_2' is not an upper bound")
print(f"  {'ok  ' if printed_C2**2 <= C3 else 'FAIL'} C_2'^2 = {float(printed_C2**2):.2f} <= C_3 = {float(C3):.6g}")
if printed_C2**2 > C3: bad.append("C_2'^2 exceeds C_3")

# --- thm:main's constant, as printed and as assembled
c0 = lk['C0_UP']; BETA = lk['BETA']
coeff = (2*(6*c0*BETA)**mpf('0.5') + 2/mpf('1e5')**mpf('0.05')
         + 2*(6*c0/BETA)**mpf('0.5')/(mpf('1e5')**mpf('0.5')*mlog(mpf('1e5'))))
printed = mpf(re.search(r"\\leq\\ (\d\.\d)\\,n\^\{3/4\}\\,\\tau", tex).group(1))
print(f"  {'ok  ' if coeff <= printed else 'FAIL'} thm:main constant: assembled {float(coeff):.6f} <= printed {float(printed)}")
if coeff > printed: bad.append("thm:main's printed constant is below what the proof assembles")
print(f"  {'ok  ' if float(printed) == float(th['C']) else 'FAIL'} paper's constant {float(printed)} matches thresholds.py C = {float(th['C'])}")
if float(printed) != float(th['C']): bad.append("paper and thresholds.py disagree on C")

print()
if bad:
    print("CROSS-CHECK FAILED:")
    for b in bad: print("  ! " + b)
    raise SystemExit(1)
print("main.tex agrees with the certificates on every constant checked.")
