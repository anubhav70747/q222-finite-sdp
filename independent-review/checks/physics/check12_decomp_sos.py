"""Checks 1-2: F = psi^T (T - ll^T) psi, and det(T) F = y^T adj(T) y + det(P) s^2 with y = T psi - s ell."""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from sympy.polys.matrices import DomainMatrix

t0 = time.time()
# Check 1
print('1. psi^T(T-ll^T)psi - F == 0 :', expand((psi.T*P*psi)[0] - F) == 0)

# Check 2 (symbolic).  y = T psi - s*ell = (T - ll^T) psi  because s*ell = ell (ell.psi) = ll^T psi.
# s^2 = r Z^2.  adj/det computed over the polynomial ring ZZ[u,v,params].
gens = [u, v, lam, n, m, C, A, D, d, al, be]
dm = DomainMatrix.from_Matrix(T).convert_to(sp.ZZ[gens])
adjT = dm.adjugate().to_Matrix()
detT = dm.domain.to_sympy(dm.det())
detP = dm.domain.to_sympy(DomainMatrix.from_Matrix(P).convert_to(sp.ZZ[gens]).det())
y = P*psi
lhs = expand(detT*F)
rhs = expand((y.T*adjT*y)[0] + detP*r*Z**2)
print('2. det(T)F - [y^T adj(T) y + det(P) s^2] == 0 (symbolic):', expand(lhs - rhs) == 0)
print('   also det(P) = det(T) - ell^T adj(T) ell :', expand(detP - (detT - r*adjT[3, 3])) == 0)
print('   terms in det(T):', len(expand(detT).args), ' time %.1fs' % (time.time() - t0))
