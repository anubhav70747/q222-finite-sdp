"""nd_core: exact coefficient tensor of the parent form F_c and the native C1/C2 sources.

gram_lift.py (the parent directory of this folder) and boundary/verify_boundary.py (native chart formulas c1_source, c2_source)
are imported by relative path from this file.  T = gram_lift.tensors() is the coefficient tensor of F_c in the 9 source weights.
Only the pieces used by the admission check are kept: fcoef_exact, source, chart_ok.
"""
import os, sys
from fractions import Fraction as Fr
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.append(os.path.normpath(os.path.join(HERE, '..')))                  # gram_lift.py
sys.path.append(os.path.normpath(os.path.join(HERE, '..', 'boundary')))      # verify_boundary.py
import gram_lift                                  # tensors()
from verify_boundary import c1_source, c2_source  # native formulas

NAMES = ['d', 'lam', 'n', 'm', 'C', 'A', 'D', 'alpha', 'beta']

# ---------------------------------------------------------------- exact F_c coefficients
T = gram_lift.tensors()                           # (9, 3, 3, 4, 4); integer-valued entries (0.2 s)
_NZ = [tuple(int(i) for i in ix) for ix in zip(*np.nonzero(T))]
_TFR = {ix: Fr(float(T[ix])) for ix in _NZ}


def fcoef_exact(src):
    """{(a,b,l,m): Fraction}, l<=m, coefficient of u^a v^b psi_l psi_m in F_c (exact zeros preserved)."""
    out = {}
    for (p, a, b, l, m), tv in _TFR.items():
        if l > m: continue
        co = Fr(src[NAMES[p]])*tv*(1 if l == m else 2)
        if co: out[(a, b, l, m)] = out.get((a, b, l, m), Fr(0)) + co
    return {k: x for k, x in out.items() if x != 0}


# ---------------------------------------------------------------- families: source, chart
def source(fam, p):
    return c2_source(*p) if fam == 'C2' else c1_source(*p)


def chart_ok(fam, p):
    """True if the verify_boundary native chart holds AND all source weights are >= 0 (d = 1). Returns (ok, src, ct)."""
    try:
        src, ct, chart = source(fam, p)
    except ZeroDivisionError:
        return False, None, None
    ok = all(chart) and all(src[k] >= 0 for k in NAMES[1:]) and all(x > 0 for x in ct[:2])
    return ok, src, ct



