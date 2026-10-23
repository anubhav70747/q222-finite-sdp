#!/usr/bin/env python3
"""Exact checks for the universal polynomial-Bob-multiplier obstruction.
Requires Python 3 and SymPy. Run: python3 check_fixed_q_obstruction.py
No numerical SDP, optimization, or floating-point arithmetic is used.
"""
import json
from pathlib import Path
import sympy as S
base=Path(__file__).resolve().parent
data=json.loads((base/'EXACT_SOURCE_AND_FUNCTIONS.json').read_text())
T,t,x,y,u,v,W,X,Y,Z,w,eta,sigma,p,q = S.symbols('T t x y u v W X Y Z w eta sigma p q')
Du,Dv=1+u*u,1+v*v
# Symbolic canonical coefficients from f=w-S-eta(A1+B0).
d=4-w
lo=(w-2-2*eta)/4
hi=(w-2+2*eta)/4
mid=(w-2)/4
lam=m=A=lo
n=D=al=hi
C=be=mid
F=4*lam*Du*(W+v*X)**2+4*n*Du*Dv*X**2+4*m*(v*W-X+u*v*Y-u*Z)**2+4*C*Dv*(u*W-Y)**2+4*A*Du*Dv*Y**2+(4*D+2*d)*(u*W+u*v*X-Y-v*Z)**2+(4*al+2*d)*Du*(v*Y-Z)**2+(4*be+2*d)*Dv*(X+u*Z)**2-2*d*Du*Dv*Z**2
I=S.eye(2); z=S.diag(1,-1); sx=S.Matrix([[0,1],[1,0]])
A1=((1-u*u)*z+2*u*sx)/Du
B1=((1-v*v)*z+2*v*sx)/Dv
kron=S.kronecker_product
Bell=kron(z,z)-kron(z,B1)+kron(A1,z)+kron(A1,B1)+eta*(kron(A1,I)+kron(I,z))
psi=S.Matrix([W,X,Y,Z])
assert S.cancel(F-(psi.T*(w*S.eye(4)-Bell)*psi)[0]*Du*Dv)==0
print('PASS: exact canonical coefficient map (all nine coefficients).')
source_order=['d','lambda','n','m','C','A','D','alpha','beta']
assert data['source_coefficient_order']==source_order
family=[d,lo,hi,lo,mid,lo,hi,hi,mid]
for name,value in zip(source_order,family):
    expected=value.subs({w:4-2*T+4*T**3,eta:1-T})
    assert S.expand(S.sympify(data['source_family'][name],locals={'T':T})-expected)==0
# Standard symmetric doubly tilted CHSH determinant, reduced by circle relations.
c,e,r,s,L=S.symbols('c e r s L')
Omega=1+c+e-c*e
H=S.Matrix([[Omega+2*(1-T),(1-c)*s,(1-e)*r,-r*s],[(1-c)*s,-Omega,-r*s,(e-1)*r],[(1-e)*r,-r*s,-Omega,(c-1)*s],[-r*s,(e-1)*r,(c-1)*s,Omega-2*(1-T)]])
parsed_H=S.Matrix([[S.sympify(z,locals={'Omega':Omega,'c':c,'e':e,'r':r,'s':s,'T':T}) for z in row] for row in data['bell_operator']['H_matrix']])
assert all(S.expand(z)==0 for z in H-parsed_H)
raw=S.Poly(S.expand((L*S.eye(4)-H).det(method='domain-ge')),r,s)
red=0
for (ir,js),a in raw.terms():
    assert ir%2==0 and js%2==0
    red+=a*(1-c*c)**(ir//2)*(1-e*e)**(js//2)
red=S.expand(red.subs({L:4-2*T+sigma,c:1-p,e:1-q},simultaneous=True))
P=sigma**3+(16-8*T)*sigma**2+(84-88*T+20*T*T)*sigma+144-240*T+112*T*T-16*T**3+8*p*q*(1-T)**2
R=(3+2*T-T*T)*(p+q)-2*p*q-18*T+12*T*T-2*T**3
assert S.expand(red-sigma*P-8*p*q*R)==0
print('PASS: exact determinant identity for every angle and both sine signs.')
# Elementary cubic bound used to prove source validity for 0<T<=1/100.
a=S.symbols('a', real=True)
assert S.expand(4*a**3-36*T*a*a+432*T**3-4*(a-6*T)**2*(a+3*T))==0
print('PASS: exact cubic lower-bound identity.')
# Exact polynomial critical test and full state (including all four coordinates).
psit=S.Matrix([1+S.Rational(2,3)*t**4*x*x*y*y,-S.Rational(1,3)*t**3*x*x*y,-t*x+S.Rational(2,3)*t**3*x*y*y,-t*t*x*y-S.Rational(2,3)*t**4*x**3*y])
for name,value in zip(['W','X','Y','Z'],psit):
    assert S.expand(S.sympify(data['critical_substitution'][name],locals={'t':t,'x':x,'y':y})-value)==0
Ft=S.expand(F.subs({w:4-2*t*t+4*t**6,eta:1-t*t,u:-t*x,v:t*y,W:psit[0],X:psit[1],Y:psit[2],Z:psit[3]},simultaneous=True))
pt=S.Poly(Ft,t)
assert all(pt.coeff_monomial(t**j)==0 for j in range(6))
M=x**4*y*y+x*x*y**4-3*x*x*y*y+1
assert S.expand(pt.coeff_monomial(t**6)-S.Rational(4,3)*(M+2))==0
print('PASS: F=t^6*(4/3)*(Motzkin+2)+O(t^8), coefficientwise.')
assert all(j[0]%2==0 for j,a in pt.terms())
# Independent normal directions in the leading full-state calculation.
aa,bb,hh=S.symbols('aa bb hh')
xi=S.Matrix([1,aa*t**3,bb*t**3,hh*t**2])
Vnum=S.Matrix([[1,t*x],[-t*x,1]])
general_psi=kron(Vnum,I)*xi
general_F=S.Poly(S.expand(F.subs({w:4-2*t*t+4*t**6,eta:1-t*t,u:-t*x,v:t*y,W:general_psi[0],X:general_psi[1],Y:general_psi[2],Z:general_psi[3]},simultaneous=True)),t)
assert S.expand(general_F.coeff_monomial(t**4)-4*(hh+x*y)**2)==0
assert S.expand(general_F.coeff_monomial(t**6).subs(hh,-x*y)-(6*(aa-S.Rational(2,3)*x*x*y)**2+6*(bb-S.Rational(2,3)*x*y*y)**2+S.Rational(4,3)*(M+2)))==0
print('PASS: both independent state-normal squares in the critical expansion.')
# The half-line ramification changes the obstruction.
z0=S.symbols('z0')
half=(x*x*z0*z0-1)**2+(x*z0*z0*(z0*z0-1))**2+2*(x*z0*(z0*z0-1))**2
assert S.expand(M.subs(y,z0*z0)-half)==0
halfdata=data['critical_halfline_gram']
features=S.Matrix([S.sympify(f,locals={'x':x,'y':y}).subs(y,z0*z0) for f in halfdata['features']])
features=features.applyfunc(lambda z:z.xreplace({S.sqrt(z0*z0):z0}))
assert S.expand((features.T*S.Matrix(halfdata['gram_matrix'])*features)[0]-(M.subs(y,z0*z0)+2))==0
print('PASS: M(x,z^2) is explicitly SOS; no half-line no-go is claimed.')
# Check every coefficient of the supplied finite expansion.
expanded='F_t coefficients in powers of t; each row is exact.\n'+''.join(f't^{j[0]}: {S.factor(a)}\n' for j,a in sorted(pt.terms(),key=lambda z:z[0]))
assert (base/'CRITICAL_PULLBACK.txt').read_text()==expanded
print('PASS: exact source/function JSON and every supplied pullback coefficient.')

# Check the half-line endpoint/reflection assertions without a positivity guess.
d0,m0,C0,A0,D0,al0,be0=S.symbols('d0 m0 C0 A0 D0 al0 be0')
p1,p2,q2,r0,r1,r2,s1,s2=S.symbols('p1 p2 q2 r0 r1 r2 s1 s2')
P0=S.Matrix([[4*m0*v*v,-4*m0*v],[-4*m0*v,4*m0+(4*be0+2*d0)*(1+v*v)]])
N0=S.Matrix([[4*C0*(1+v*v)+4*D0+2*d0,(4*D0+2*d0)*v],[(4*D0+2*d0)*v,(4*D0+2*d0)*v*v]])
E0=S.Matrix([[4*A0*(1+v*v)+(4*al0+2*d0)*v*v,-(4*al0+2*d0)*v],[-(4*al0+2*d0)*v,4*al0-2*d0*v*v]])
J=S.Matrix([[p1*v+p2*v*v,q2*v*v],[r0+r1*v+r2*v*v,s1*v+s2*v*v]])
selector_symbols={str(z):z for z in [v,p1,p2,q2,r0,r1,r2,s1,s2]}
parsed_J=S.Matrix([[S.sympify(z,locals=selector_symbols) for z in row] for row in data['halfline_selector']['J_plus']])
assert all(S.expand(z)==0 for z in J-parsed_J)
T0=S.diag(1,-1)
parsed_Jminus=S.Matrix([[S.sympify(z,locals=selector_symbols) for z in row] for row in data['halfline_selector']['J_minus_on_v_le_0']])
assert all(S.expand(z)==0 for z in T0*J.subs(v,-v)*T0-parsed_Jminus)
for mat in (P0,N0,E0):
    assert all(S.expand(z)==0 for z in mat.subs(v,-v)-T0*mat*T0)
assert all(S.expand(z)==0 for z in T0*J.subs(v,-v)*T0-S.Matrix([[-p1*v+p2*v*v,-q2*v*v],[-r0+r1*v-r2*v*v,-s1*v+s2*v*v]]))
H0=(P0+E0).row_join(-E0-J).col_join((-E0-J.T).row_join(N0+E0)).subs(v,0)
i,j,k,l=S.symbols('i j k l')
vec=S.Matrix([i,j,k,l])
endpoint=(4*m0+4*be0+2*d0)*j*j+(4*C0+4*D0+2*d0)*k*k-2*r0*j*k+4*A0*(i-k)**2+4*al0*(j-l)**2
assert S.expand((vec.T*H0*vec)[0]-endpoint)==0
print('PASS: half-line endpoint quadratic and sign-sensitive reflection.')
print('No half-line adequacy or global affine lift is asserted.')
