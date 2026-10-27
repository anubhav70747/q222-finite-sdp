"""Bounded exact algebra checks for the R5.5 dependency repair.

No solver, numerical test, source/test identification, or global theorem claim.
"""
import json
from pathlib import Path
import sympy as s

checks = []

def zero(name, expression):
    value = s.cancel(expression)
    assert value == 0, (name, value)
    checks.append(name)

u, v, W, X, Y, Z = s.symbols('u v W X Y Z', real=True)
a, B, h, t, r, m, K, n, A, C, D, alpha = s.symbols(
    'a B h t r m K n A C D alpha', real=True)
g = 1-a*(1-r-h*B)
U, V, Du, Dv = u*u, v*v, 1+u*u, 1+v*v
M = v*W-X+u*v*Y-u*Z
J = u*W+u*v*X-Y-v*Z
P = (m*M*M+n*Du*Dv*X*X+C*Dv*(u*W-Y)**2+A*Du*Dv*Y*Y
     +(D+s.Rational(1,2))*J*J+(alpha+s.Rational(1,2))*Du*(v*Y-Z)**2
     +s.Rational(1,2)*Dv*(X+u*Z)**2-s.Rational(1,2)*Du*Dv*Z*Z)
xi, eta = X+r*u*Z, Y-h*u*W
R1, R2, T = u*W-v*Z/t, v*W-u*Z/a, u*v*W-B*Z
Rmc = -xi+u*v*eta+(1-g)*R2+h*u*T
zero('R5.2 as a full-state row identity', Rmc-(M-g*v*W))
remainder = P.subs(m,K)+2*(m-K)*g*v*W*M-(m-K)*g*g*v*v*W*W
zero('R5.5 against every packet of the four-state parent',
     P-(m-K)*Rmc**2-remainder)

# Exact derivative of the unnormalized WA8 table, with every other actual
# source coefficient fixed. No assumption that a capped source is stationary.
kappa, ell = 1/t, 1/a
b2 = 2*m*h+C*(1-h)**2+(A+alpha+s.Rational(1,2))*h*h
b4 = h*h*(m+A+alpha+s.Rational(1,2))
c4, d2 = n*r*r, (n+s.Rational(1,2))*r*r+2*D*r
d4 = (n+D+s.Rational(1,2))*r*r
beta_s = (B*d4+kappa*c4)/(kappa*ell**2)
gamma_s = b4+d4/(kappa*ell)
j0 = (B*d4-A*h*h*kappa*kappa*ell+d2*ell)/(2*kappa*ell**2)
j1 = d4/(2*kappa*ell)
fs = (b2-c4/ell**2-alpha/B**2-B*d4/(kappa*ell**2))/2
S = s.Matrix([
 [(C+D+s.Rational(1,2))*(1-h)**2+A*h*h+A*h*h*U, u*v*(j0+j1*U), -v*(j0+j1*U)],
 [u*v*(j0+j1*U), m+beta_s*U, u*(fs-j1*U)],
 [-v*(j0+j1*U), u*(fs-j1*U), alpha/B**2+gamma_s*U]])
p0 = (D+r*(n+s.Rational(1,2)))/kappa
p1 = r*(n+D+s.Rational(1,2))/kappa
q1 = (n*r+B*p1)/ell
fx = -m*h+(D+s.Rational(1,2))*(1-h)-p0-q1
Lx = s.Matrix([[u*v*(p0+p1*U), -m+q1*U, u*(fx-p1*U)]])
py = -(C+D+s.Rational(1,2))*(1-h)+A*h
py1 = A*h
kv = -C*(1-h)+A*h+m+(alpha+s.Rational(1,2))*h
kuv = h*(A+m+alpha+s.Rational(1,2))
kz = -m*(1-r)+r*(D+s.Rational(1,2))-alpha-s.Rational(1,2)
qy = (-kz-kappa*py1-B*kuv)/ell
Ly = s.Matrix([[py+py1*U, u*v*qy, v*(kv-qy+kuv*U)]])
H = s.Matrix([
 [n*Du*Dv+m+(D+s.Rational(1,2))*U*V+s.Rational(1,2)*Dv, -u*v*(m+D+s.Rational(1,2))],
 [-u*v*(m+D+s.Rational(1,2)), A*Du*Dv+C*Dv+D+s.Rational(1,2)+m*U*V+(alpha+s.Rational(1,2))*Du*V]])
G = H.row_join(Lx.col_join(Ly)).col_join(Lx.T.row_join(Ly.T).row_join(S))
Lm = s.Matrix([-1,u*v,0,1-g,h*u])
Hm = s.zeros(5)
for i,j,q in [(0,3,-1),(1,4,v),(3,3,2-g),(3,4,h*u)]:
    Hm[i,j] = q
    Hm[j,i] = q
for i in range(5):
    for j in range(i,5):
        zero(f'R5.3 matrix coefficient ({i+1},{j+1})',
             s.diff(G[i,j],m)-Lm[i]*Lm[j]-g*Hm[i,j])

# Contact check without any radical simplification: impose the exact three
# graph-root equations and both normal equations via their consequences.
zero('centered row vanishes at the stationary contact graph',
     (M-g*v*W).subs({X:-r*u*Z,Y:h*u*W}).expand()
     .subs(u*u, a*B).subs(u*Z,a*v*W).expand())

# R6 source cancellation, jointly in s,b,g, before taking residues.
ss, b, gg = s.symbols('s b g', real=True)
rr = 1-ss*(1+h*ss*b-gg)
aa, BB = 1/ss, ss*ss*b
ca = 1+aa*BB*rr
Sc = 1+BB*(1-h+aa*rr)
Delta = (aa+(1+aa*BB)+(1+aa*BB)*BB+2*aa*BB*(1+aa*BB)*rr
         +aa**3*BB**2*rr**2+2*aa**2*BB*rr
         +t*BB*(1+aa+aa*BB+aa**2*BB*rr*(2-rr)))
Ms = (1+aa*BB)*(1+aa*rr)*(h+t)-aa*(1+t*BB)*(1-rr)*(1-h)
source = {
 'mg': aa**2*(1+t*BB)*(1-rr)*ca*Sc/(2*Delta),
 'n/s^2': ((1+aa*rr)*ca*Sc-rr*Delta)/(2*rr*Delta*ss**2),
 'A': t*ca*(1+aa*rr)*Sc/(2*Delta*h),
 'C': ca*Ms/(2*(1-h)*Delta),
 'D': (1+aa*BB)*(1+aa*rr)*(1-t*aa*BB*rr)*Sc/(2*t*(ca-(1-h)/t)*Delta),
 'alpha/(s^3 b^2)': BB**2*(1+aa*rr)*Ms/(2*(1-h*BB)*Delta*ss**3*b*b),
}
expected = {'mg':s.Rational(1,2), 'n/s^2':(1-b*(h+t))/2,
            'A':t/(2*h), 'C':(h+t)/(2*(1-h)),
            'D':1/(2*(t-(1-h))), 'alpha/(s^3 b^2)':(h+t)/2}
for name, expr in source.items():
    cancelled = s.cancel(expr)
    num,den = s.fraction(cancelled)
    den0 = s.factor(den.subs({ss:0,gg:0}))
    assert den0 != 0, (name, den0)
    zero('R6 joint residue '+name, cancelled.subs({ss:0,gg:0})-expected[name])

# Four-state R6 finite remote chart, with n_0=alpha_0=0.
As, Cs, Ds = t/(2*h), (h+t)/(2*(1-h)), 1/(2*(t-(1-h)))
Bfin = P.subs({m:0,n:0,alpha:0,A:As,C:Cs,D:Ds})
kernel = s.expand(Bfin.subs(X,v*W+u*v*Y-u*Z))
vec = s.Matrix([W,Y,Z])
Q = s.hessian(kernel,vec)/2
det_expected = (Du**2*Dv**2*(h+t)/(8*h*(1-h)*(t-(1-h)))
                *(V**2+t*(h+t)*U**2*V+h*U*V**2))
zero('R6 RT10/RB24 determinant on the complete m-event kernel', Q.det()-det_expected)
z0_squares = (Cs*Dv*(u*W-Y)**2+As*Du*Dv*Y*Y
              +(Ds+s.Rational(1,2))*(u*W+u*v*(v*W+u*v*Y)-Y)**2
              +s.Rational(1,2)*Du*V*Y*Y
              +s.Rational(1,2)*Dv*(v*W+u*v*Y)**2)
zero('R6 Z=0 kernel restriction as a positive sum', kernel.subs(Z,0)-z0_squares)

# Homogeneous corrections are polynomials, including literal p=q=0.
pp, qq, mu = s.symbols('p q mu', real=True)
correction = 2*mu*v*W*M-mu*gg*v*v*W*W
cu = s.expand(pp*pp*correction.subs(u,1/pp))
cv = s.expand(qq*qq*correction.subs(v,1/qq))
zero('u-infinity centered correction has zero residue', cu.subs({pp:0,gg:0}))
zero('v-infinity correction is 2 mg W(W+uY)',
     cv.subs({qq:0,gg:0})-2*mu*W*(W+u*Y))
cd = s.expand(pp*pp*qq*qq*correction.subs({u:1/pp,v:1/qq}))
zero('double-infinity centered correction vanishes', cd.subs({pp:0,qq:0,gg:0}))

receipt = {'scope':'R5.2--R5.5 full-parent extraction, R6 source bounds, finite remote kernel, homogeneous correction only',
           'checks':checks,'count':len(checks),
           'limitations':['No global finite-SDP theorem verified.',
                          'No contact conductor or exhaustive chart join verified.',
                          'No inherited source admission transferred from R5 to R6.']}
Path(__file__).with_name('exact-checks.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
print(json.dumps(receipt,indent=2))
