#!/usr/bin/env python3
"""Exact, deterministic symbolic C2 strict-source transfer checks, no solver or samples.

Imports the original verified full native C2 parent generator bundled in this
folder. Checks the unique sixth-root extreme diagonal identities, positive
sign-adaptive contact-vanishing Gram carrier, necessary fractional residual
coefficient equations, full-state 3x3 normal-block identity, and the five-ray
balanced C2 determinant algebra. DOES NOT assert global Z6 admission.
"""
import sympy as sp
import c2_fixed_map as c2

a,b,c,t,r=c2.a,c2.b,c2.c,c2.t,c2.r
n,m,C,A,D,alpha=c2.n,c2.m,c2.C,c2.A,c2.D,c2.al
h=b-a*a*t
wc,gaps=c2.native_weights()

# Source identity and uniquely forced extreme Gram diagonals.
assert sp.cancel((2*A*h-a*a*r*(2*n+1)).subs(wc,simultaneous=True))==0
k=2*c**4*r*(2*n+1)
dY=4*A*c**4*h**2/a**2
dX=2*a*a*c**4*r*r*(2*n+1)
assert sp.cancel((dY-k*h).subs(wc,simultaneous=True))==0
assert sp.expand(dX-k*a*a*r)==0

P=c2.target_template()
psi=[b/a,-c*r,c*h/a,sp.Integer(1)]
def f(key):
    ix,jx,i,j=key
    return sp.factor(P[key]*c**(ix//6)*(a*c)**(jx//6)*psi[i]*psi[j])
assert sp.cancel(f((12,0,2,2))-dY)==0
assert sp.cancel(f((0,12,1,1))-dX)==0
assert (6,0,1,2) not in P or sp.expand(P[(6,0,1,2)])==0
assert (0,6,1,1) not in P or sp.expand(P[(0,6,1,1)])==0
print('PASS 1: exact native source identity and unique positive sixth-root diagonals')

# The 3x3 contact-zero carrier Gram in basis (U phiY,V phiX,phiX).
U,V,x,y=sp.symbols('U V phi_X phi_Y')
Bpos=(dY-dX)*(U*y-x)**2+dX*(U*y-V*x)**2
Bneg=(dX-dY)*(V*x-x)**2+dY*(U*y-V*x)**2
# Both have PSD sums of contact-vanishing squares in the appropriate sign branch.
assert sp.expand(Bpos.subs({U:1,V:1,x:1,y:1}))==0
assert sp.expand(Bneg.subs({U:1,V:1,x:1,y:1}))==0
assert sp.simplify(sp.expand(Bpos).coeff(U,2).coeff(y,2)-(dY))==0
assert sp.simplify(sp.expand(Bpos).coeff(V,2).coeff(x,2)-(dX))==0
assert sp.simplify(sp.expand(Bneg).coeff(U,2).coeff(y,2)-(dY))==0
assert sp.simplify(sp.expand(Bneg).coeff(V,2).coeff(x,2)-(dX))==0
assert sp.simplify(sp.Poly(Bpos,U,V,x,y).coeff_monomial(U*x*y)-(-2*(dY-dX)))==0
assert sp.simplify(sp.Poly(Bneg,U,V,x,y).coeff_monomial(V*x*x)-(-2*(dX-dY)))==0
print('PASS 2: sign-adaptive two-square PSD contact carrier, exact off-support defects')

# Every allowed Z6 Gram sees these two extreme coefficients ONLY as diagonals.
bl=c2.blocks()
features=[z for B in bl for z in B]
for key,want in [((12,0,2,2),(6,0,2)),((0,12,1,1),(0,6,1))]:
    pairs=[(i,j) for block in bl for i in block for j in block
           if c2.product_key(i,j)==key]
    assert pairs==[(want,want)],(key,pairs)
print('PASS 3: both diagonal capacities are uniquely fixed in all four PSD blocks')

# Exact fractional routes available to repair the sign-dependent off-support term
# after the carrier saturates its two extreme diagonal capacities.
for key,ban in [((6,0,1,2),(6,0,2)),((0,6,1,1),(0,6,1))]:
    pairs={(min(i,j),max(i,j)) for block in bl for i in block for j in block
           if c2.product_key(i,j)==key and i!=ban and j!=ban}
    assert pairs, key
    assert any(z[0]%6 not in (0,) or z[1]%6 not in (0,) for ij in pairs for z in ij),key
    print('    remaining polynomial key',key,':',len(pairs),'potential fractional/coupled pair classes')
print('PASS 4: necessary remaining off-support transfer classes identified')

# The complete four-state parent P=psi^T Q psi, and the top-left B=R^TR.
u,v,W,X,Y,Z=c2.u,c2.v,c2.W,c2.X,c2.Y,c2.Z
DU,DV=1+u*u,1+v*v
F=c2.full_parent()
psi4=(W,X,Y,Z)
Q=sp.Matrix(4,4,lambda i,j:sp.diff(F,psi4[i],psi4[j])/2)
rows=[]
for aa in (1,u):
    for bb in (1,v):
        rows.append([0,2*sp.sqrt(n)*aa*bb,0])
rows.append([2*sp.sqrt(m)*v,-2*sp.sqrt(m),2*sp.sqrt(m)*u*v])
for bb in (1,v):
    rows.append([2*sp.sqrt(C)*bb*u,0,-2*sp.sqrt(C)*bb])
for aa in (1,u):
    for bb in (1,v):
        rows.append([0,0,2*sp.sqrt(A)*aa*bb])
rows.append([sp.sqrt(4*D+2)*u,sp.sqrt(4*D+2)*u*v,-sp.sqrt(4*D+2)])
for aa in (1,u):
    rows.append([0,0,sp.sqrt(4*alpha+2)*aa*v])
for bb in (1,v):
    rows.append([0,sp.sqrt(2)*bb,0])
R=sp.Matrix(rows)
assert R.shape==(16,3)
B=Q[:3,:3]
assert all(sp.expand((R.T*R-B)[i,j])==0 for i in range(3) for j in range(3))
assert sp.expand(F-(sp.Matrix(psi4).T*Q*sp.Matrix(psi4))[0])==0
print('PASS 5: exact complete-parent 16-row polynomial Gram of top 3x3 state block')

# 3x3 adjugate completion is checked for an arbitrary symmetric block,
# arbitrary two-normal coupling, and arbitrary fourth-state coefficient.
b00,b01,b02,b11,b12,b22=sp.symbols('b00 b01 b02 b11 b12 b22')
z0,z1,z2,dd=sp.symbols('z0 z1 z2 dd')
x0,x1,x2,z4=sp.symbols('x0 x1 x2 z4')
Bg=sp.Matrix([[b00,b01,b02],[b01,b11,b12],[b02,b12,b22]])
zg=sp.Matrix([z0,z1,z2]);xg=sp.Matrix([x0,x1,x2]);kg=Bg*xg+zg*z4
Pg=Bg.row_join(zg).col_join((zg.T).row_join(sp.Matrix([[dd]])))
full=(xg.T*Bg*xg)[0]+2*(xg.T*zg)[0]*z4+dd*z4*z4
assert sp.expand(Bg.det()*full-(kg.T*Bg.adjugate()*kg)[0]-Pg.det()*z4*z4)==0
print('PASS 6: exact full four-state adjugate/Schur identity (arbitrary entries)')

# A universal five-ray positive decomposition for the C2 six-term
# stationary Laurent determinant. Pairs and triangles, no approximations.
s,tv=sp.symbols('s T',positive=True)
P0,Q0,Z0,K0,L0=sp.symbols('P0 Q0 Z0 K0 L0',nonnegative=True)
weights={( -1,0):P0+K0,(-1,1):Q0+L0,(0,-1):Z0+L0,
         (0,1):Z0+K0,(1,-1):Q0+K0,(1,0):P0+L0}
Laurent=sum(w*(s**i*tv**j-1) for (i,j),w in weights.items())
circuits=[s+1/s-2,s/tv+tv/s-2,tv+1/tv-2,
          1/s+tv+s/tv-3,tv/s+1/tv+s-3]
assert sp.cancel(Laurent-(P0*circuits[0]+Q0*circuits[1]+Z0*circuits[2]+K0*circuits[3]+L0*circuits[4]))==0
assert sum(w*ij[0] for ij,w in weights.items())==0
assert sum(w*ij[1] for ij,w in weights.items())==0
xx,yy,zz=sp.symbols('xx yy zz')
cubes=xx**3+yy**3+zz**3-3*xx*yy*zz
sq=sp.expand(sum(hh*(ii-jj)**2 for hh in (xx,yy,zz) for ii,jj in ((xx,yy),(xx,zz),(yy,zz)))/2)
assert sp.expand(cubes-sq)==0
print('PASS 7: exact five-ray C2 stationary Laurent identity and three-term cubical SOS')
# Explicit monomial transfer from normalized Laurent variables to the
# fixed finite Puiseux root list for UV*L. No independent-angle inverse.
Up,Vp,Uc,Vc=sp.symbols('Up Vp Uc Vc', positive=True)
ss,tt=Up/Uc,Vp/Vc
assert sp.cancel(Up*Vp*(ss+1/ss-2)-Vp*(Up-Uc)**2/Uc)==0
assert sp.cancel(Up*Vp*(ss/tt+tt/ss-2)-(
       sp.sqrt(Uc/Vc)*Vp-sp.sqrt(Vc/Uc)*Up)**2)==0
assert sp.cancel(Up*Vp*(tt+1/tt-2)-Up*(Vp-Vc)**2/Vc)==0
A1,A2,A3=(Uc*Vp)**sp.Rational(1,3),(Up*Vp**2/Vc)**sp.Rational(1,3),(Up**2*Vc/Uc)**sp.Rational(1,3)
B1,B2,B3=(Uc*Vp**2/Vc)**sp.Rational(1,3),(Up*Vc)**sp.Rational(1,3),(Up**2*Vp/Uc)**sp.Rational(1,3)
assert sp.simplify(A1*A2*A3-Up*Vp)==0
assert sp.simplify(B1*B2*B3-Up*Vp)==0
assert sp.simplify((A1**3+A2**3+A3**3)-Up*Vp*(1/ss+tt+ss/tt))==0
assert sp.simplify((B1**3+B2**3+B3**3)-Up*Vp*(tt/ss+1/tt+ss))==0
print('PASS 8: three opposite squares and both nine-square triangle root lists are exact')

print('ALL EIGHT CHECKS PASSED. UNIVERSAL Z6 PSD ADMISSION REMAINS UNPROVED.')
