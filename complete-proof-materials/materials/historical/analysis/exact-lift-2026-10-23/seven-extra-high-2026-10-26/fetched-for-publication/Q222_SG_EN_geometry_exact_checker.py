#!/usr/bin/env python3
"""Exact SG1--SG10 / EN1--EN16 geometry checks.

SymPy only. No SDP or numerical sampling. All equalities are symbolic
polynomial/rational identities. The support check enumerates the 256 finite
sign patterns only to audit the analytical hitting-set lemma.
"""
import sympy as s
from itertools import combinations
from sympy.polys.matrices import DomainMatrix


def zero(expr, message):
    assert s.cancel(expr) == 0, message


def matzero(M, message):
    assert all(s.cancel(a) == 0 for a in M), message


# The eight outer Newton exponents, with positive product packets.
cycle = ('lambda', 'D', 'C', 'A', 'alpha', 'm', 'beta', 'n')
packets = {
    (-1,-1): [('lambda','alpha')],
    (-1, 0): [('lambda','D'),('alpha','m')],
    (-1, 1): [('m','D')],
    ( 0,-1): [('alpha','A'),('lambda','n')],
    ( 0, 1): [('beta','m'),('C','D')],
    ( 1,-1): [('n','A')],
    ( 1, 0): [('beta','n'),('C','A')],
    ( 1, 1): [('beta','C')],
}
E = list(packets)
exclusive_pairs = {
    frozenset((cycle[i],cycle[j]))
    for i,j in combinations(range(8),2)
    if min((i-j)%8,(j-i)%8)>=3
}
assert len(exclusive_pairs) == 12
normal_reps = E
zero_support_patterns = []
for mask in range(256):
    active = {cycle[i] for i in range(8) if mask >> i & 1}
    omitted = set(cycle)-active
    support={e for e,terms in packets.items()
             if any(set(t)<=active for t in terms)}
    coercive = all(any(e[0]*w[0]+e[1]*w[1]>0 for e in support)
                   for w in normal_reps)
    has_exclusive = any(pair<=omitted for pair in exclusive_pairs)
    assert coercive == (not has_exclusive), (active,coercive)
    if not support: zero_support_patterns.append(active)
max_zero = [a for a in zero_support_patterns
            if not any(a<b for b in zero_support_patterns)]
assert len(max_zero)==8 and all(len(a)==3 for a in max_zero)
compatible_triples=[z for z in combinations(range(8),3)
    if all(min((i-j)%8,(j-i)%8)<=2 for i,j in combinations(z,2))]
assert len(compatible_triples)==8
print('SG9: exactly 12 exclusive pairs; eight maximal zero-outer triples; eight compatible omitted triples.')
# Check the three physical event transformations on every actual event label.
events = {
 'lambda':(0,1,1,1), 'n':(0,0,1,-1),
 'm':(1,1,1,-1), 'C':(1,0,-1,1),
 'A':(0,0,-1,1), 'D':(1,1,-1,1),
 'alpha':(0,1,-1,-1), 'beta':(1,0,1,-1)}
def gA(ev):
 x,y,aa,bb=ev;return(1-x,y,aa,(-1)**y*bb)
def gB(ev):
 x,y,aa,bb=ev;return(x,1-y,(-1)**(1-x)*aa,bb)
def gT(ev):
 x,y,aa,bb=ev;return(y,x,(-1)**y*bb,(-1)**x*aa)
for fn,j0 in ((gA,5),(gB,3),(gT,2)):
 for i,name in enumerate(cycle):
  assert fn(events[name]) == events[cycle[(j0-i)%8]], (fn,name)
print('SG8: gA,gB,gT act exactly as the stated physical cycle reflections.')


# The state matrix of SG2.
lam,n,m,C,A,D,alpha,beta,d,u,v=s.symbols('lam n m C A D alpha beta d u v')
Du=1+u*u; Dv=1+v*v
rows=[s.Matrix(k) for k in [
  (1,v,0,0),(0,1,0,0),(v,-1,u*v,-u),(u,0,-1,0),
  (0,0,1,0),(u,u*v,-1,-v),(0,0,v,-1),(0,1,0,u),
  (0,0,0,1)]]
coefficients=[4*lam*Du,4*n*Du*Dv,4*m,4*C*Dv,4*A*Du*Dv,
              4*D+2*d,(4*alpha+2*d)*Du,(4*beta+2*d)*Dv,
              -2*d*Du*Dv]
P=sum((w*x*x.T for w,x in zip(coefficients,rows)),s.zeros(4))

# SG5--SG6, as a polynomial equality for all nine source variables and u,v.
rawdet=DomainMatrix.from_Matrix(P).det().as_expr()
t1=d+2*(lam+n+D);t2=d+2*(alpha+m+A)
t3=d+2*(lam+C+D);t4=d+2*(alpha+beta+m)
t5=d+2*(lam+beta+n);t6=d+2*(alpha+C+A)
t7=d+2*(beta+n+m);t8=d+2*(C+A+D)
q={(0,0):4*lam*alpha*t7*t8,
   (0,1):4*(lam*D*t7*t6+alpha*m*t5*t8),
   (0,2):4*m*D*t5*t6,
   (1,0):4*(alpha*A*t7*t3+lam*n*t8*t4),
   (1,2):4*(beta*m*t1*t6+C*D*t5*t2),
   (2,0):4*n*A*t3*t4,
   (2,1):4*(beta*n*t3*t2+C*A*t1*t4),
   (2,2):4*beta*C*t1*t2}
hcen=rawdet.subs({u:1,v:1})/256-sum(q.values())
assert s.Poly(s.cancel(rawdet-16*Du**2*Dv**2*
      (hcen*u*u*v*v+sum(w*u**(2*i)*v**(2*j)
                          for (i,j),w in q.items()))),u,v).is_zero
print('SG5--SG6: full generic nine-variable determinant identity verified.')

# Signed-core identity SG15, and its nonsingular determinant SG20.
PH=P.subs({lam:0,n:0,m:0,C:0,A:0,D:0,alpha:0,beta:0,d:1})
zero(s.factor(PH.det())+16*Du**2*Dv**2*u*u*v*v,'SG20')
xy=PH.extract([1,2],[1,2]); wz=PH.extract([0,3],[0,3]); cross=PH.extract([0,3],[1,2])
B0=1+v*v+u*u*v*v
Schur=2*Du*Dv/(B0**2-u*u*v*v)*s.Matrix([
  [u*u*v*v,-u*v*B0],[-u*v*B0,u*u*v*v]])
matzero(wz-cross*xy.inv()*cross.T-Schur,'SG15')
print('SG15, SG20: signed-core Schur and determinant verified.')

# EN6 and EN12, with adjugate B44 and B24 computed as explicit minors.
def b44(X):return X.extract([0,1,2],[0,1,2]).det(method='domain-ge')
def b24(X):return X.extract([0,1,2],[0,2,3]).det(method='domain-ge')
Fcof=(d+2*D)*((d+2*alpha)*(m*Dv+C*u*u*v*v)+2*A*Dv*(m+C*u*u))\
    +2*m*C*(d*(u*u*v*v+Dv)+2*alpha+2*D*v*v*Du)
P1=P.subs({lam:0,n:0})
zero(u*b44(P1)+b24(P1)-16*u*v*v*Du**2*Dv*Fcof,'EN6')
P02=P.subs({lam:0,beta:0,n:0})
P2=P.subs({lam:0,beta:0})
zero(b24(P2)-b24(P02),'EN12 B24 independence')
zero(u*(b44(P2)-b44(P02))-4*u*n*Du*Dv*
     P02.extract([0,2],[0,2]).det(method='domain-ge'),'EN12 correction')
print('EN6, EN12: positive-cofactor identities verified.')

# Both exact CHSH-preserving full-parent reflections.
G=s.Matrix([[v,-1,0,0],[-1,-v,0,0],[0,0,-v,1],[0,0,1,v]])
map1={lam:n,n:lam,m:C,C:m,A:alpha,alpha:A,D:beta,beta:D}
matzero(u*u*G.T*P.subs(map1,simultaneous=True).subs(u,1/u)*G-Dv*P,
       'C1 full reflection')
J=s.Matrix([[0,0,0,1],[0,1,0,0],[0,0,1,0],[1,0,0,0]])
map2={lam:beta,beta:lam,m:D,D:m,C:alpha,alpha:C}
matzero(u*u*v*v*J.T*P.subs(map2,simultaneous=True).subs(
    {u:1/v,v:1/u},simultaneous=True)*J-P, 'C2 full reflection')
print('C1/C2 full-parent symmetry identities verified.')

# C1 EN8 selected stationary minor and EN7 Cramer identity.
SS,Y0,th=s.symbols('S Y0 theta')
psi1=s.Matrix([SS,u*(th-1),u*Y0,1])
P1=P.subs({lam:0,n:0,C:2})
def contactcol(Q,psi):
    return s.Matrix(list(Q*psi)+[
        s.expand((psi.T*Q.diff(t)*psi)[0]) for t in (u,v)])
columns=[d,m,A,D,alpha,beta]
M1=s.Matrix.hstack(*[contactcol(P1.diff(t),psi1) for t in columns])
minor=DomainMatrix.from_Matrix(M1).det().as_expr()
Rm=v*SS-u*(th-1)+u*u*v*Y0-u
RD=u*SS+u*u*v*(th-1)-u*Y0-v
Dfac=Y0*Du*(u*Y0+v)*(1-u*v*Y0)*(SS+u*v*(th-1))*Rm*(-RD)
zero(minor+8192*th*u*u*v*Du**2*Dv**2*Dfac,'EN8')
Qcol=contactcol(P.diff(C),psi1)
M1m=M1.copy();M1m[:,1]=-2*Qcol
num=DomainMatrix.from_Matrix(M1m).det().as_expr()
zero(num*(SS+u*v*(th-1))*Rm-2*u**3*(SS-Y0)*(1-th)*Dv*minor,
     'EN7 exact Cramer identity')
print('EN8 selected stationary determinant and EN7 polynomial Cramer identity verified.')

# Two exact positive stationary seeds.
seeds=[(s.Matrix([1,-s.Rational(1,2),s.Rational(1,4),1]),
         {u:1,v:1,d:s.Rational(27,2),lam:0,n:0,m:4,C:2,A:9,
          D:s.Rational(5,4),alpha:s.Rational(13,4),beta:s.Rational(9,4)}),
       (s.Matrix([1,-s.Rational(1,4),s.Rational(1,4),1]),
         {u:1,v:1,d:1,lam:0,beta:0,n:s.Rational(1,8),m:s.Rational(3,8),
          C:s.Rational(1,6),A:s.Rational(5,8),D:s.Rational(3,8),
          alpha:s.Rational(1,6)})]
for psi,point in seeds:
    matzero(P.subs(point)*psi,'seed kernel')
    F=(psi.T*P*psi)[0]
    for t in (u,v):zero(s.diff(F,t).subs(point),'seed angle jet')
print('EN9 and EN C2 seeds: exact state kernel and angle jets verified.')

# EN16, the unweighted C2 six-column Jacobian determinant.
a,c,b,t,r=s.symbols('a c b t r')
z=c*c; h=b-a*a*t;g=(1+z)*b-(1-r+a*a*z*t);j=1+z*r-t
psi2=s.Matrix([b/a,-c*r,c*(b/a-a*t),1]); cp={u:c,v:a*c}
W,X,Y,Z=s.symbols('W X Y Z')
rm=v*W-X+u*v*Y-u*Z;rc=u*W-Y;rd=u*W+u*v*X-Y-v*Z;ral=v*Y-Z
funs=[Du*Dv*X**2,rm**2,Dv*rc**2,Du*Dv*Y**2,rd**2,Du*ral**2]
vals={W:b/a,X:-c*r,Y:c*(b/a-a*t),Z:1,u:c,v:a*c}
M2=s.Matrix([[s.diff(F,xx).subs(vals) for F in funs]
             for xx in [W,X,Y,Z,u,v]])
det2=DomainMatrix.from_Matrix(M2).det().as_expr()
Delta=1+(1+z)*b+z*(1+z)*b*b+2*z*(1+z)*b*r+z*z*r*r+2*z*r\
      +a*a*z*((1+z)*b+1+z*r*(2-r))
zero(det2+64*a*c**6*r*g*t*h*j*(1-z*h)*(1+z)**2*(1+a*a*z)**2*Delta,
     'EN16 positive-factor determinant')
print('EN16: exact C2 six-column determinant, with manifestly positive Delta, verified.')

# C2 native stationary weight identities are in a separate short checker
# (Q222_C2_native_stationarity_check.py) to avoid large expanded forms.

# C1 short native SD9--SD10 identities, proving its stationary certificate.
kap,eta,zeta,aa,tt=s.symbols('kap eta zeta aa tt')
vc=kap*eta*zeta;Lc=1-aa*eta;Hc=1+(1-aa)*eta;Qc=1-kap*zeta+aa*vc
RR=1+kap*(1-aa*eta); TT=1+kap-kap*zeta*tt; KK=1+kap
P4=1+kap+eta*kap*(aa-kap-1)-zeta*kap*tt\
   +vc*(-zeta*kap*tt**2+2*kap*tt+aa*(1-tt)+tt)
Pal=1+kap-aa*eta*kap+zeta*((1+kap)-(2+kap)*tt)\
    +vc*(-zeta*tt**2+aa*(tt-1)+tt)
Pd=(1+kap)**2-aa*eta*kap**2\
  +zeta*kap*((kap+1)-(kap+2)*tt)\
  +vc*(-zeta*kap*tt**2+kap*aa*(tt-1)+kap*tt+kap+1)
c2=2*KK*(1-tt*zeta)*Lc*(1+vc)*Hc*Qc
cb=2*KK*Lc*(1+(1-tt)*zeta)*Hc*Qc
c3=2*KK*Lc*Hc*Qc*(1-kap*eta+tt*vc)
c5=2*KK*Lc*Hc*Qc*(1+vc)
c4=2*Qc*P4;cal=2*Lc*Hc*Pal;cd=4*Lc*Qc*Pd
balance=[cb-c2-zeta*c3,c3-c5+kap*eta*cb,c2-(1-tt*zeta)*c5,
  2*kap*c2+2*(1+vc)*Lc*c4-2*kap**2*zeta*Qc*cal-Qc*cd,
  kap*Lc*c2+c5+(1+vc)*Lc**2*c4+kap*Qc**2*cal\
    +s.Rational(1,2)*(kap*zeta*Lc**2-2*Lc-eta)*cd,
  kap*RR*c2+tt*kap*vc*cb+aa*vc*c3+TT*c5+KK*vc*Lc**2*c4\
   -kap**2*zeta*KK*Qc*Lc*cal-KK*s.Rational(1,2)*(Qc*Lc+eta)*cd]
for index,expr in enumerate(balance):
    assert s.Poly(expr,kap,eta,zeta,aa,tt).is_zero, ('SD9--SD10',index)
print('C1 SD9--SD10: all six exact native polynomial balances verified.')

# SG13: the two 3x3 block identity for the other closed face.
ss,ww,xx,yy,zz,xi=s.symbols('ss ww xx yy zz xi')
p=4*C+4*D+2*d; q0=4*m+4*beta+2*d; rr=4*D+2*d
B1=s.Matrix([[4*m,-4*m,xi],[-4*m,q0+4*n,q0],[xi,q0,q0]])
B2=s.Matrix([[p+4*A,4*D,-p],[4*D,4*D,-(4*m+rr)-xi],
             [-p,-(4*m+rr)-xi,p]])
x1=s.Matrix([ww,xx,ss*zz]);x2=s.Matrix([yy,zz,ss*ww])
G=4*n*xx**2+4*m*(ww-xx-ss*zz)**2+4*C*(ss*ww-yy)**2+4*A*yy**2\
 +rr*(ss*ww-yy-zz)**2+(4*beta+2*d)*(xx+ss*zz)**2\
 -4*d*yy*zz-2*d*zz**2
zero(G-(x1.T*B1*x1)[0]-(x2.T*B2*x2)[0],'SG13 full polynomial')
print('SG13: both 3x3 blocks reproduce the entire limiting parent.')

# The two precise determinant reductions behind the closed-face converses.
Pn=P.subs({lam:0,m:0})
# q_{10}, q_{12}, q_{20}, q_{21}, q_{22} are the only nonzero outer terms.
q10n=q[(1,0)].subs({lam:0,m:0});q12n=q[(1,2)].subs({lam:0,m:0})
hn=Pn.subs({u:1,v:1}).det()/256-sum(val.subs({lam:0,m:0}) for val in q.values())
K0=s.diag(4*C+4*D+2*d,4*A,4*alpha)
K2=s.diag(4*C,4*A+4*alpha+2*d,4*D)
K1=s.Matrix([[0,0,-(4*D+2*d)],[0,0,-(4*alpha+2*d)],
             [-(4*D+2*d),-(4*alpha+2*d),0]])
K=K0+v*K1+v*v*K2
zero(K.det()-(8*Dv/(d+2*(beta+n)))*(q10n+hn*v*v+q12n*v**4),
     'SG12 determinant')

# SG13 limiting 4x4 state matrix (1/2 Hessian), determinant 16(q02+h s^2+q20 s^4).
P13=P.subs({lam:0,alpha:0})
q02p=q[(0,2)].subs({lam:0,alpha:0});q20p=q[(2,0)].subs({lam:0,alpha:0})
h13=P13.subs({u:1,v:1}).det()/256-sum(val.subs({lam:0,alpha:0}) for val in q.values())
yz=[ww,xx,yy,zz]
Gs=s.Matrix([[s.diff(G,z0,z1)/2 for z1 in yz] for z0 in yz])
zero(DomainMatrix.from_Matrix(Gs).det().as_expr()
     -16*(q02p+h13*ss**2+q20p*ss**4),'SG13 determinant')
print('SG12/SG13: closed-face limiting determinant factorizations verified.')

print('ALL CHECKS PASSED')
