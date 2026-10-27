"""Exact identities supporting the symmetric and homogeneous lift proof.

Run with Python + SymPy.  There are no numerical samples or SDP solves.
The proof of the universal positivity equivalences is in PROOF.md.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).with_name('_deps')))
import sympy as S

a, h, lam, c, d, r, s, u, x, y, z, t = S.symbols(
    'a h lam c d r s u x y z t', real=True)
I2 = S.eye(2)
X = S.Matrix([[0, 1], [1, 0]])
Z = S.diag(1, -1)
Y = S.Matrix([[0, -S.I], [S.I, 0]])
kp = S.kronecker_product
A0, A1 = kp(Z, I2), kp(c*Z+r*X, I2)
B0, B1 = kp(I2, Z), kp(I2, d*Z+s*X)
CHSH = A0*(B0+B1)+A1*(B0-B1)
H = h*CHSH + a*(A0+B0)
Omega = 1+c+d-c*d
H_display = S.Matrix([
    [h*Omega+2*a, h*(1-c)*s, h*(1-d)*r, -h*r*s],
    [h*(1-c)*s, -h*Omega, -h*r*s, h*(d-1)*r],
    [h*(1-d)*r, -h*r*s, -h*Omega, h*(c-1)*s],
    [-h*r*s, h*(d-1)*r, h*(c-1)*s, h*Omega-2*a]])

checks = []
def eq(name, expression):
    if isinstance(expression, S.MatrixBase):
        expression = list(expression)
    elif not isinstance(expression, (list, tuple)):
        expression = [expression]
    for item in expression:
        assert S.cancel(S.expand(item)) == 0, (name, S.factor(item))
    checks.append(name)

def rad_reduce(expression):
    """Polynomial reduction by r^2=1-c^2, s^2=1-d^2."""
    p = S.Poly(S.expand(expression), r, s)
    return S.expand(sum(co*r**(i % 2)*s**(j % 2)*
        (1-c*c)**(i//2)*(1-d*d)**(j//2)
        for (i,j), co in p.terms()))

eq('physical Kronecker operator equals displayed 4x4 block', H-H_display)
eq('CHSH square identity', [rad_reduce(q) for q in
    CHSH**2 - 4*S.eye(4) - 4*r*s*kp(Y,Y)])

coeff = H.charpoly(u).all_coeffs()
poly = S.expand(sum(rad_reduce(q)*(u**(4-i))
    for i,q in enumerate(coeff)))
slack_poly = S.Poly(S.expand(poly.subs(u, lam+u)), u)
e = slack_poly.all_coeffs()
eq('characteristic e1', e[1]-4*lam)
eq('characteristic e2', e[2]-(6*lam**2-4*a**2-8*h**2))
eq('characteristic e3', e[3]-(4*lam**3-8*lam*a**2-
    16*lam*h**2-8*h*a**2*Omega))
D0 = ((lam-2*h)**2-4*a*a)*(lam+2*h)**2
D = S.expand(D0 + 64*h*h*(4*h*h-a*a)*x*y*(x+y) +
    (32*lam*h*a*a+128*h*h*a*a-256*h**4)*x*y-
    256*h**4*x*x*y*y)
eq('complete homogeneous determinant',
    e[4].subs({c:1-2*x,d:1-2*y})-D)
p, q = S.symbols('p q', real=True)
eq('global determinant diagonal difference',
    D.subs({x:p*p,y:q*q})-D.subs({x:p*q,y:p*q})-
    64*h*h*(4*h*h-a*a)*p*p*q*q*(p-q)**2)

b, Habs, v, w = S.symbols('b Habs v w', nonnegative=True)
eq('e2 positivity expansion',
    (6*lam**2-4*b*b-8*Habs*Habs).subs(lam,2*Habs+2*b+v)-
    (16*Habs**2+48*Habs*b+20*b*b+24*(Habs+b)*v+6*v*v))
e3_lower = 4*lam**3-8*lam*b*b-16*lam*Habs*Habs-16*Habs*b*b
eq('e3 positivity expansion',
    e3_lower.subs(lam,2*Habs+2*b+v)-
    (16*b*(b+2*Habs)**2 +
    (32*Habs**2+96*Habs*b+40*b*b)*v +
    24*(Habs+b)*v*v+4*v**3))

rt2 = S.sqrt(2)
Q = S.Matrix([[1,1,0,0], [0,0,1,1], [0,0,1,-1],
    [1,-1,0,0]])/rt2
eq('explicit Bell-basis orthogonality', Q.T*Q-S.eye(4))
H3 = S.Matrix([[2*h*c,2*a,0],
    [2*a,2*h*(1+c-c*c),2*h*(1-c)*r],
    [0,2*h*(1-c)*r,-2*h*(1+c-c*c)]])
Hd = H.subs({d:c,s:r})
B = Q.T*Hd*Q
eq('symmetric 3x3 block and scalar -2hc',
    [S.expand(f).subs(r*r,1-c*c) for f in
     (B-S.diag(H3,S.Matrix([[-2*h*c]])))])

K = S.Matrix([
    [lam-2*h+2*lam*t*t+(lam+2*h)*t**4,
     -2*a*(1+t*t)**2,0],
    [-2*a*(1+t*t)**2,
     lam-2*h+(2*lam-8*h)*t*t+(lam+2*h)*t**4,-8*h*t**3],
    [0,-8*h*t**3,lam+2*h+(2*lam+8*h)*t*t+(lam-2*h)*t**4]])
ct = (1-t*t)/(1+t*t)
rt = 2*t/(1+t*t)
eq('quartic tangent-half-angle identity',
    K-(1+t*t)**2*(lam*S.eye(3)-H3.subs({c:ct,r:rt})))
R = S.diag(1,1,-1)
eq('parity involution', K.subs(t,-t)-R*K*R)
eq('h=0 singular case', K.subs(h,0)-(1+t*t)**2*
    S.Matrix([[lam,-2*a,0],[-2*a,lam,0],[0,0,lam]]))

def sym_matrix(prefix,n):
    cache = {(i,j):S.Symbol(f'{prefix}{i+1}{j+1}')
             for i in range(n) for j in range(i,n)}
    return S.Matrix(n,n,lambda i,j:cache[min(i,j),max(i,j)])
U,V = sym_matrix('u',5),sym_matrix('v',4)
x1,x2,x3=S.symbols('x1 x2 x3')
xx=S.Matrix([x1,x2,x3])
zp=S.Matrix([x1,x2,t*x3,t*t*x1,t*t*x2])
zm=S.Matrix([x3,t*x1,t*x2,t*t*x3])
gram=S.expand((zp.T*U*zp+zm.T*V*zm-xx.T*K*xx)[0])
matching = {
    U[0,0]:lam-2*h, U[0,1]:-2*a, U[1,1]:lam-2*h,
    V[0,0]:lam+2*h, U[0,2]:-V[0,1], U[1,2]:-V[0,2],
    U[0,3]:(2*lam-V[1,1])/2,
    U[0,4]:-4*a-U[1,3]-V[1,2],
    U[1,4]:(2*lam-8*h-V[2,2])/2,
    U[2,2]:2*lam+8*h-2*V[0,3],
    U[2,3]:-V[1,3], U[2,4]:-8*h-V[2,3],
    U[3,3]:lam+2*h, U[3,4]:-2*a, U[4,4]:lam+2*h,
    V[3,3]:lam-2*h}
assert len(S.Poly(gram,t,x1,x2,x3).terms()) == 16
eq('all sixteen affine Gram equations', gram.subs(matching))

CB,DB=(B0+B1)/2,(B0-B1)/2
J=S.eye(4)
SOS20=(J-CB-A1*DB)**2+sum(
    (((J-A0)*(J+By))**2 for By in (B0,B1)), S.zeros(4))/4
eq('critical one-sided SOS', [rad_reduce(f) for f in
    4*J-(CHSH+2*A0)-SOS20])

# The following resultant is an algebraic identity, not a root-selection rule.
Ddiag = D.subs({h:1,x:z,y:z})
F4 = (4*lam**4-4*lam**3*a*a+16*lam**3+11*lam*lam*a**4-
    48*lam*lam*a*a-16*lam*lam+8*lam*a**6-4*lam*a**4-
    80*lam*a*a-128*lam+20*a**6-84*a**4+64*a*a-128)
resultant = S.resultant(Ddiag, S.cancel(S.diff(Ddiag,z)/z), z)
eq('symmetric stationary resultant factorization',
    resultant-2**36*(lam-2)**2*(2*lam+a*a)**2*F4)

print(json.dumps({'status':'all exact identities verified',
    'checks':checks,'count':len(checks),
    'scope':'Named polynomial identities only; universal proof in PROOF.md.'},indent=2))
