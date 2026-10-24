"""Check 5: L(U,V) = det(.)/(16 Du^2 Dv^2 U V), U=u^2, V=v^2.  Which (U,V) exponents does each piece reach?
det(P) = det(T) - ell^T adj(T) ell,   ell^T adj(T) ell = 2 d Du Dv det(T_WXY)."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import *
from sympy.polys.matrices import DomainMatrix

U, V = symbols('U V', positive=True)
params = [lam, n, m, C, A, D, d, al, be]
R = sp.ZZ[[u, v] + params]
det = lambda M: R.to_sympy(DomainMatrix.from_Matrix(M).convert_to(R).det())

detT, det3 = det(T), det(T[:3, :3])                       # adj(T)_44 = det of the W,X,Y block
sub_term = expand(r*det3)                                  # ell^T adj(T) ell
print('det(P) == det(T) - 2d Du Dv det(T_WXY):', expand(det(P) - detT + sub_term) == 0)


def laurent(expr):
    """expr(u,v) -> {(i,j): coeff} with expr = 16 Du^2 Dv^2 U V * sum coeff U^i V^j."""
    q = sp.cancel(expand(expr).subs({u**2: U, v**2: V}) / (16*(1 + U)**2*(1 + V)**2))
    assert not q.has(u, v) and sp.denom(q) == 1, 'not even / not divisible'
    pq = sp.Poly(expand(q), U, V)
    return {(i - 1, j - 1): sp.factor(cf) for (i, j), cf in pq.terms()}


LT, LS = laurent(detT), laurent(sub_term)
print('det(T) exponents:', sorted(LT))
print('subtracted exponents:', sorted(LS))
for key in sorted(LS):
    print('  subtracted', key, ':', LS[key])


def nonneg(cf):
    return all(x >= 0 for x in sp.Poly(expand(cf), *params).coeffs())


print('det(T) L-coeffs with a negative monomial (params):', [k for k, cf in LT.items() if not nonneg(cf)])
print('   centre coeff of det(T):', expand(LT[(0, 0)]))
# generic weights a1..a8 as well
a = list(symbols('a1:9', positive=True))
Ra = sp.ZZ[[u, v] + a]
LTa = laurent(Ra.to_sympy(DomainMatrix.from_Matrix(T_of(a)).convert_to(Ra).det()))
print('det(T) L-coeffs with a negative monomial (generic a_k):',
      [k for k, cf in LTa.items() if any(x < 0 for x in sp.Poly(expand(cf), *a).coeffs())])
print('   centre coeff, generic a_k:', expand(LTa[(0, 0)]))
print('net L(det P) coefficients:')
for key in sorted(LT):
    net = expand(LT[key] - LS.get(key, 0))
    negs = [t for t in sp.Add.make_args(net) if t.as_coeff_Mul()[0] < 0]
    print(f'  {key}: hit={key in LS}, negative monomials after cancellation: {negs if negs else "none"}')
