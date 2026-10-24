"""Check 3: adj(T) = sum_{|S|=3} (prod a)(prod c) w_S w_S^T ; det T = sum_{|S|=4} (prod a)(prod c) det[h_S]^2.
Symbolic in u, v and fully generic weights a1..a8 (so it holds for any weights)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from sympy.polys.matrices import DomainMatrix

a = list(symbols('a1:9', positive=True))
Tg = T_of(a)
dm = DomainMatrix.from_Matrix(Tg).convert_to(sp.ZZ[[u, v] + a])
adj_direct = dm.adjugate().to_Matrix()
det_direct = dm.domain.to_sympy(dm.det())

adj_cb = zeros(4, 4)
for S in itertools.combinations(range(8), 3):
    w = cross3(*[h[k] for k in S])
    adj_cb += sp.prod([a[k]*c[k] for k in S]) * w * w.T
print('adj(T) - CauchyBinet3 == 0 :', (adj_direct - adj_cb).applyfunc(expand) == zeros(4, 4))

det_cb = sum(sp.prod([a[k]*c[k] for k in S]) * Matrix.hstack(*[h[k] for k in S]).det()**2
             for S in itertools.combinations(range(8), 4))
print('det(T) - CauchyBinet4 == 0 :', expand(det_direct - det_cb) == 0)
