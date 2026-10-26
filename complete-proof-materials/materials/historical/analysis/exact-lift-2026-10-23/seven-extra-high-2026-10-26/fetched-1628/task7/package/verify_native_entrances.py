"""Exact symbolic checks supporting the C1/C2 coefficient-density interface.
No numerical source sampling, solver, GitHub access, or Gram search.
"""
import itertools
import sympy as sp


def equal(name, lhs, rhs):
    diff=sp.factor(sp.together(lhs-rhs))
    assert (diff == sp.zeros(*diff.shape) if isinstance(diff, sp.MatrixBase) else diff == 0), (name, diff)
    print('PASS',name)

# C1 forward chart: EN13 / A3, retaining symbolic independent source variables.
k,f,g,a,t = sp.symbols('k f g a t', positive=True)
kap,eta,zeta = k*k,f*f,g*g
Lc=1-a*eta
uc,vc=k,k*f*g
Wc,Xc,Yc,Zc = g/f,-k*(1-t*g*g),k*g*Lc/f,sp.Integer(1)
r=sp.cancel(Yc/(uc*Wc))
angle=1+Xc/uc
eta_inv=vc/(uc*Wc)
zeta_inv=Wc*vc/uc
equal('C1 inverse r=Lc',r,Lc)
equal('C1 inverse eta',eta_inv,eta)
equal('C1 inverse zeta',zeta_inv,zeta)
equal('C1 inverse a',(1-r)/eta_inv,a)
equal('C1 inverse t',angle/zeta_inv,t)
R=1+kap*Lc
T=1+kap-kap*zeta*t
M=R-t
J=T-a
H=1-angle
E=1-kap*eta*H
Q=1-kap*zeta*r
Rm=vc*Wc-Xc+uc*vc*Yc-uc
RD=uc*Wc+uc*vc*Xc-Yc-vc
equal('C1 contact m-gap',Rm,uc*zeta*M)
equal('C1 contact D-gap',RD,-vc*J)
equal('C1 contact E-gap',Wc+vc*Xc,Wc*E)
equal('C1 contact Q-gap',1-vc*Yc,Q)
equal('C1 native E coquotient',1-kap*eta+t*kap*eta*zeta,E)

# C2 forward chart: EN14 / C2 native exact contact.
a,b,c,t,r = sp.symbols('a b c t r', positive=True)
uc,vc=c,a*c
Wc,Xc,Yc,Zc=b/a,-c*r,c*(b/a-a*t),sp.Integer(1)
h=b-a*a*t
g=(1+c*c)*b-(1-r+a*a*c*c*t)
j=1+c*c*r-t
Rm=vc*Wc-Xc+uc*vc*Yc-uc
RD=uc*Wc+uc*vc*Xc-Yc-vc
equal('C2 inverse c',uc,c)
equal('C2 inverse a',vc/uc,a)
equal('C2 inverse b',(vc/uc)*Wc,b)
equal('C2 inverse r',-Xc/uc,r)
equal('C2 inverse t',(uc*Wc-Yc)/vc,t)
equal('C2 inverse h',(vc/uc)*Yc/uc,h)
equal('C2 m-gap',Rm/uc,g)
equal('C2 D-gap',-RD/vc,j)
equal('C2 1-zh-gap',1-vc*Yc,1-c*c*h)
equal('C2 b-a²zr-gap',(vc/uc)*(Wc+vc*Xc),b-a*a*c*c*r)

# The 16 deterministic strategies produce an orthogonal constant/Bell-feature basis.
rows=[]
for A0,A1,B0,B1 in itertools.product((-1,1),repeat=4):
    rows.append([1,A0,A1,B0,B1,A0*B0,A0*B1,A1*B0,A1*B1])
V=sp.Matrix(rows)
equal('deterministic feature orthogonality',V.T*V,16*sp.eye(9))
print('ALL EXACT DENSITY INTERFACE CHECKS PASSED')
