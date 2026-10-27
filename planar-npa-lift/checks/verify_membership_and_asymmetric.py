"""Exact identities for the projection lift and signed weighted-ray extension.

Requires SymPy. No numerical optimizers, samples, or repository imports.
The inequalities and closedness arguments are proved in ../proofs/main.tex.
"""
import sympy as s

t, h, a, lam = s.symbols('t h a lam', real=True)
w = s.symbols('w1:4', real=True)
z = s.symbols('z1:17', real=True)
U = s.Matrix([
    [z[0], z[1]/2, z[4]/2, z[6], z[7]/2],
    [z[1]/2, z[2], z[5]/2, z[7]/2, z[8]],
    [z[4]/2, z[5]/2, z[9], z[10]/2, z[11]/2],
    [z[6], z[7]/2, z[10]/2, z[12], z[13]/2],
    [z[7]/2, z[8], z[11]/2, z[13]/2, z[14]],
])
V = s.Matrix([
    [z[3], z[4]/2, z[5]/2, z[9]],
    [z[4]/2, z[6], z[7]/2, z[10]/2],
    [z[5]/2, z[7]/2, z[8], z[11]/2],
    [z[9], z[10]/2, z[11]/2, z[15]],
])

def symvar(n, prefix):
    entries = {(i,j): s.Symbol(f'{prefix}{i+1}{j+1}')
               for i in range(n) for j in range(i,n)}
    return s.Matrix(n,n,lambda i,j: entries[min(i,j),max(i,j)])

G, H = symvar(5,'g'), symvar(4,'q')
left = s.Matrix([
    G[0,0],G[0,1],G[1,1],H[0,0],G[0,2]+H[0,1],G[1,2]+H[0,2],
    2*G[0,3]+H[1,1],G[0,4]+G[1,3]+H[1,2],2*G[1,4]+H[2,2],
    G[2,2]+2*H[0,3],G[2,3]+H[1,3],G[2,4]+H[2,3],
    G[3,3],G[3,4],G[4,4],H[3,3]])
right = s.Matrix([lam-2*h,-2*a,lam-2*h,lam+2*h,0,0,2*lam,-4*a,
                  2*lam-8*h,2*lam+8*h,0,-8*h,lam+2*h,-2*a,lam+2*h,lam-2*h])
assert s.expand((s.Matrix(z).T*left)[0]-s.trace(U*G)-s.trace(V*H)) == 0
T = right.jacobian((lam,a,h))
dual = T.T*s.Matrix(z)
ell = z[0]+z[2]+z[3]+2*z[6]+2*z[8]+2*z[9]+z[12]+z[14]+z[15]
M = 2*z[1]+4*z[7]+2*z[13]
S = 2*z[0]+2*z[2]-2*z[3]+8*z[8]-8*z[9]+8*z[11]-2*z[12]-2*z[14]+2*z[15]
assert dual == s.Matrix([ell,-M,-S])
G0=s.diag(1,1,2,1,1); H0=s.diag(1,2,2,1)
assert left.xreplace(dict(zip(G,G0))|dict(zip(H,H0))) == T*s.Matrix([1,0,0])
assert s.expand(s.trace(U*G0)+s.trace(V*H0)-ell)==0
assert s.Matrix(list(U)+list(V)).jacobian(z).rank()==16
print('PASS sixteen-row adjoint, behavior coordinates, positive anchor and injectivity')

# Independent physical determinant for arbitrary signed tilts.
alpha,beta,c,d,r,q,x,y = s.symbols('alpha beta c d r q x y', real=True)
I=s.eye(2); Z=s.diag(1,-1); X=s.Matrix([[0,1],[1,0]])
kr=s.kronecker_product
Bell=kr(Z,Z)+kr(Z,d*Z+q*X)+kr(c*Z+r*X,Z)-kr(c*Z+r*X,d*Z+q*X)
Bell+=alpha*kr(Z,I)+beta*kr(I,Z)

def circle_reduce(expr):
    expr=s.rem(s.Poly(expr,r),s.Poly(r*r-(1-c*c),r)).as_expr()
    return s.expand(s.rem(s.Poly(expr,q),s.Poly(q*q-(1-d*d),q)).as_expr())

char=Bell.charpoly(lam).all_coeffs()
det=circle_reduce(char[-1]+char[-2]*lam+char[-3]*lam**2+char[-4]*lam**3+lam**4)
A=4-alpha**2; B=4-beta**2
D0=((lam-2)**2-(alpha+beta)**2)*((lam+2)**2-(alpha-beta)**2)
C=32*lam*alpha*beta+64*alpha**2+64*beta**2-256
D=D0+x*y*(C+64*A*x+64*B*y-256*x*y)
assert s.expand(det.subs({c:1-2*x,d:1-2*y})-D)==0
fm=-(lam+alpha)**2+beta**2+4+4*beta*(2*y-1)
fp=(lam-alpha)**2-beta**2-4+4*beta*(2*y-1)
assert s.expand(D.subs(x,1)+fm*fp)==0
Omega=1+c+d-c*d
assert circle_reduce(s.trace(Bell*Bell))-4*(alpha**2+beta**2+4)==0
assert s.expand(circle_reduce(s.trace(Bell**3))-24*alpha*beta*Omega)==0
print('PASS general signed physical determinant, edge factorization and characteristic inputs')

# Radicals are removed by diagonal congruence on the open ray.
rho,zeta,R,Sv=s.symbols('rho zeta R Sv', real=True)
Om=2-4*rho*zeta**2
ray_slack=s.Matrix([
 [lam-Om-alpha-beta,-2*zeta*q,-2*rho*zeta*r,r*q],
 [-2*zeta*q,lam+Om-alpha+beta,r*q,2*rho*zeta*r],
 [-2*rho*zeta*r,r*q,lam+Om+alpha-beta,2*zeta*q],
 [r*q,2*rho*zeta*r,2*zeta*q,lam-Om+alpha+beta]])
Delta=s.diag(r*q,r,q,1)
P=(Delta*ray_slack*Delta).applyfunc(s.expand).subs({r*r:R,q*q:Sv})
assert all(not entry.has(r,q) for entry in P)
P=P.subs({R:4*zeta*(1-zeta),Sv:4*rho*zeta*(1-rho*zeta)})
assert max(s.degree(entry,zeta) for entry in P if entry!=0)==6
assert all(s.diff(entry,lam,2)==0 for entry in P)
print('PASS signed ray congruence, degree six and affinity in the candidate bound')
