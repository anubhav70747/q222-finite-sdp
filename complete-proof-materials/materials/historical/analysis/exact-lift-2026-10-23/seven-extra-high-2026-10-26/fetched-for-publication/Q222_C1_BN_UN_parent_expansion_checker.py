"""Exact full original reciprocal-C1 parent expansion in both free normals.

Checks the scalar, both diagonal normals, their mixed hard block and
both first-normal coefficients. Independent formal polynomial variables;
no numerical point evaluations, SDP solvers, or sampling.
Requires SymPy.
"""
import sympy as sy
U,V,W,Z,Xe,Ny=sy.symbols('u w W Z Xe Ny')
s,k,l,e,r,rt,m,d,beta,alpha,c0,a0=sy.symbols('s k l e r rt m d beta alpha c0 a0')
S=s**2; a=U**2;c=V**2
D=d-sy.Rational(1,2)
C=c0/(S*(1-r)); A=a0/(S*r)
X=Xe-e*V*W
Y=r*U*W+s*rt*Ny
P=(m*(W-V*X+U*Y-U*V*Z)**2+C*(1+c)*(U*W-Y)**2+A*(1+a)*(1+c)*Y**2+d*(U*V*W+U*X-V*Y-Z)**2 + (alpha+sy.Rational(1,2))*(1+a)*(Y-V*Z)**2+(beta+sy.Rational(1,2))*(1+c)*(X+U*Z)**2-sy.Rational(1,2)*(1+a)*(1+c)*Z**2)
poly=sy.Poly(sy.expand(P),Xe,Ny)
AN=m*(1+r*a+e*c)**2+(beta+sy.Rational(1,2))*e**2*c*(1+c)
AP=(1+c)*(c0*(1-r)+a0*r*(1+a))+S*(d*(1-r-e)**2*c+(alpha+sy.Rational(1,2))*(1+a)*r**2)
K=m*(1+r*a+e*c)+d*(1-r-e)+(alpha+sy.Rational(1,2))*(1+a)*r+(beta+sy.Rational(1,2))*e*(1+c)
FZ=D+alpha*c+beta*a+(m+alpha+beta+sy.Rational(1,2))*a*c
HX=m*c+d*a+(beta+sy.Rational(1,2))*(1+c)
HY=r*c0*(1+c)+(1-r)*a0*(1+a)*(1+c)+S*r*(1-r)*(m*a+d*c+(alpha+sy.Rational(1,2))*(1+a))
HXY=-U*V*s*rt*(m+d)
AX=-m*(1+r*a+e*c)+d*a*(1-r-e)-e*(beta+sy.Rational(1,2))*(1+c)
YP=(1+c)*(-c0+a0*(1+a))+S*(m*(1+r*a+e*c)-d*c*(1-r-e)+(alpha+sy.Rational(1,2))*(1+a)*r)
YC=l*YP-m*a+d-(alpha+sy.Rational(1,2))*(1+a)
# Here YC corresponds to the native coefficient ell ; its relation to Y normal
# includes ell as factor in G, which can be checked via a direct normal expression.

def ck(label,expr):
 v=sy.cancel(expr)
 print(label, "PASS" if v==0 else v, flush=True)
 assert v==0
ck('scalar',poly.coeff_monomial(1)-(AN+a/S*AP)*W**2+2*K*U*V*W*Z-FZ*Z**2)
ck('HX',poly.coeff_monomial(Xe**2)-HX)
ck('HY',poly.coeff_monomial(Ny**2).subs(rt**2,r*(1-r))-HY)
ck('HXY',poly.coeff_monomial(Xe*Ny)-2*HXY)
ck('LX',poly.coeff_monomial(Xe)/2-(V*AX*W+U*(m*c-d+(beta+sy.Rational(1,2))*(1+c))*Z))
# Direct Y normal before its monic contact-conductor substitution:
# missing ell*YP coefficient? Y_P contains source combination; direct L_Y =rt (YP*uW + s*w *[-ma+d-(alpha+.5)(1+a)] Z).
ly_expected=rt*(YP*U*W/s+s*V*(-m*a+d-(alpha+sy.Rational(1,2))*(1+a))*Z)
ck('LY',poly.coeff_monomial(Ny)/2-ly_expected)
print('FULL parent scalar & normals PASS')
