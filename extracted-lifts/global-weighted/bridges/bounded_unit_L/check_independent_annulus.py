"""Fresh exact check of the bounded-unit-L annulus cofactor formula.

Derives the raw4x4 matrix and its determinant/principal cofactors from
the displayed selected parent, without reading author's computed data.
"""
from pathlib import Path
from itertools import permutations
import json
import sympy as S

u,c,s,L,r=S.symbols('u c s L r',real=True)
CC,AA,aa=S.symbols('CC AA aa',real=True)
P,X,Y,T=S.symbols('P X Y T',real=True)
K=(X*X+c*c*(P-Y)**2+s*s*Y*Y-2*s*P*T
   +2*u*c*X*(s*P-s*Y+T)+u*u*s*s*X*X+u*u*(s*Y-T)**2)
F=CC*(P-Y)**2+AA*(1+u*u)*Y*Y+aa*(1+u*u)*(s*Y-T)**2+K/2
v=S.Matrix([P,X,Y,T])
G=S.hessian(F,v)/2

def det_expand(A):
    n=A.rows; out=0
    for perm in permutations(range(n)):
        sign=(-1)**sum(perm[i]>perm[j] for i in range(n) for j in range(i+1,n))
        term=S.Integer(sign)
        for i,j in enumerate(perm): term*=A[i,j]
        out+=term
    return S.expand(out)

def circle_reduce(f):
    return S.Poly(f,c).rem(S.Poly(c*c+s*s-1,c)).as_expr()

det=circle_reduce(det_expand(G))
cof=circle_reduce(sum(det_expand(G.minor_submatrix(i,i)) for i in range(4)))
sub={CC:1/(2*L*(1-r)),AA:(1+L*(1-r))/(2*L*r),aa:L/(2*(1-L*r))}
det_native=S.factor(det.subs(sub))
cof_native=S.factor(cof.subs(sub))
t=u*u; q=1-s*s
N2=(1+r*r+(1-r)*q*(1+r*(1-q))+(1-r)*(2-q*r)*t
    +(1-r)**2*(1-q)*t*t)
N1=2*(1-r)*q+(2-r*q*(3-q))*t+(1-r)*(1-q)*t*t
N0=1+(2-q)*t
claimed=2*(L*L*N2+L*N1+N0)/((1+t)*(L*q+t*(1-q))*(1+L*(1-r)))
assert S.cancel(cof_native/det_native-claimed)==0

# Source-independent event coordinates and all seven coefficient matrices.
e1,e2,e3,e4=S.symbols('e1 e2 e3 e4',real=True)
ev=S.Matrix([e1,e2,e3,e4])
Fe=S.expand(F.subs({P:e1+e2,X:e4,Y:e2,T:s*e2-e3},simultaneous=True))
Ge=S.hessian(Fe,ev)/2
Em=lambda i,j:S.Matrix(4,4,lambda a,b:1 if (a,b)==(i,j) else 0)
arrays=[S.diag(0,AA,aa+S.Rational(1,2),0),Em(0,0)/2,
 -(Em(0,1)+Em(1,0)+Em(1,1))/2,
 (Em(0,2)+Em(2,0)+Em(1,2)+Em(2,1))/2,
 (Em(0,3)+Em(3,0)+Em(1,3)+Em(3,1))/2,
 -(Em(2,3)+Em(3,2))/2,Em(3,3)/2]
mons=[u*u,c*c,s*s,s,u*c*s,u*c,u*u*s*s]
pred=S.diag(CC,AA,aa,S.Rational(1,2))
for mon,A in zip(mons,arrays): pred+=mon*A
assert all(S.expand(a-b)==0 for a,b in zip(Ge,pred))

out={'scope':'Fresh symbolic annulus identity check; inequalities audited analytically.',
 'raw_Gram_from_full_selected_parent':True,
 'trace_inverse_formula_exact':True,'seven_event_coordinate_matrices_exact':True,
 'native_determinant_factorization':str(det_native)}
(Path(__file__).parent/'independent_annulus_checks.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
