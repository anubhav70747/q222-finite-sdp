"""exact_check: exact (Fraction) contact-rescaled right-hand side beta(theta) at a rational source theta.
Only exact_beta is used by the admission check."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import nd_core


def exact_beta(fam, th):
    src, ct, chart = nd_core.source(fam, th); assert all(chart)
    fc = nd_core.fcoef_exact(src); ps = [abs(x) for x in ct[2:]]; uc, vc = ct[0], ct[1]; out = {}
    for (a, b, l, m), co in fc.items(): out[(6*a, 6*b, l, m)] = ps[l]*ps[m]*uc**a*vc**b*co
    return out


