"""Independent targeted algebra for strong C2 / IH source and test modules.

This script builds the parent afresh and imports no archive checker.  It checks
identities, not the order/compactness arguments supplied in the audit text.
"""
import sys
import sympy as s
import json
from pathlib import Path

out = Path(__file__).with_suffix('.json')
checks=[]
def zero(name, expr):
    num,den=s.fraction(s.cancel(expr))
    ok=s.expand(num)==0
    checks.append({'name':name,'passed':ok})
    out.write_text(json.dumps({'sympy':s.__version__,'scope':'Fresh algebra only; analytical audit separate','checks':checks,'complete':False},indent=2)+'\n')
    print(name,ok,flush=True)
    if not ok: raise AssertionError((name,s.factor(num)))
def regular(name,expr,subs,expected):
    rat=s.cancel(expr)
    num,den=s.fraction(rat)
    n0=s.factor(num.subs(subs, simultaneous=True)); d0=s.factor(den.subs(subs, simultaneous=True))
    if d0==0: raise AssertionError((name,'singular denominator',s.factor(den)))
    zero(name,s.cancel(n0/d0)-expected)
    checks[-1]['denominator_residue']=str(d0)

# Fresh source, with B=H+L; compensating factors cancel all opening poles.
z,y,H,B,r=s.symbols('z y H B r', positive=True)
L=B-H; ac=1+z*r; srcS=ac+L; j=ac-L/y
delta=z+(1+z)*B+(1+z)*B**2+2*z*(1+z)*B*r+z**3*r**2+2*z**2*r+y*((1+z)*B+z+z**2*r*(2-r))
srcM=(1+z)*(B+z*r)*(H+y)-z*(1+y)*(1-r)*L
G=B+z*(H+r-1)
nv=((B+z*r)*ac*srcS-r*delta)/(2*r*delta)
mv=z**2*(1+y)*(1-r)*ac*srcS/(2*G*delta)
Cv=ac*srcM/(2*L*delta)
Av=y*ac*(B+z*r)*srcS/(2*H*delta)
Dv=(1+z)*(B+z*r)*(B-z*y*r)*srcS/(2*y*j*delta)
alv=(B+z*r)*srcM/(2*(1-H)*delta)
# Cancel the exact m numerator (z*q-B)(z*q+B), with z*q+B=G,
# before multiscale substitution.  This avoids gratuitous expression swell.
qshape=H+r-1
kbase_num=s.expand(H*ac*(B+z*r)*srcS*(y+2)+(1-H)*(B+z*r)*srcM+(1+y)*(1-r)*ac*srcS*(z*qshape-B)+r*(1+2*z)*((B+z*r)*ac*srcS-r*delta)+qshape**2*delta)
kv=2*kbase_num/delta+8*Dv*H*r
zero('k m-term exact cancellation',4*mv*qshape**2-4*B**2*mv/z**2-2*(1+y)*(1-r)*ac*srcS*(z*qshape-B)/delta)
zero('GE positive numerator identity',srcM-((1+z)*(B+z*r)*H+y*G+y*z*j+z*r*L*(1+y)))
zero('GE n numerator at r=1',((B+z*r)*ac*srcS-r*delta).subs(r,1)+(H+y)*(B+z)*(1+z))

rho,ss,tau,ee=s.symbols('rho ss tau ee', positive=True)
strong={y:B/rho,r:ss*rho}; end0={B:1/tau,y:1/(tau*rho),r:ss*rho}
cs=(1+B-H)/B
for name,expr,target in [
 ('nu/rho',(4*nv+2)*r/rho,2*cs),('nr/rho',nv*r/rho,(cs-ss)/2),
 ('D/rho^2',Dv/rho**2,cs/2),('CL',Cv*L,s.Rational(1,2)),
 ('Bm/z^2',B*mv/z**2,cs/2),('AH',Av*H,(1+B-H)/2),
 ('alpha(1-H)',alv*(1-H),B/2),
 ('k/(B rho)',kv/(B*rho),2*(B**2*ss+B*H*ss+B*H-B*ss-H**2+H)/B**2)]:
    regular('strong IH simultaneous '+name,expr.subs(strong,simultaneous=True),{z:0,rho:0},target)
for name,expr,target in [
 ('nu/rho',(4*nv+2)*r/rho,2),('nr/rho',nv*r/rho,(1-ss)/2),
 ('D/rho^2',Dv/rho**2,s.Rational(1,2)),('BC',B*Cv,s.Rational(1,2)),
 ('Bm/z^2',B*mv/z**2,s.Rational(1,2)),('AH/B',Av*H/B,s.Rational(1,2)),
 ('alpha(1-H)/B',alv*(1-H)/B,s.Rational(1,2)),('k/(B rho)',kv/(B*rho),2*ss)]:
    regular('weak IH rho0 simultaneous '+name,expr.subs(end0,simultaneous=True),{z:0,tau:0,rho:0},target)
rho1=(1+z*r-ee)/(1-H*tau)
end1={B:1/tau,y:1/(tau*rho1)}
for name,expr,target in [
 ('De',Dv*ee,s.Rational(1,4)),('nu',(4*nv+2)*r,1),('nr',nv*r,(1-2*r)/4),
 ('BC',B*Cv,s.Rational(1,4)),('Bm/z^2',B*mv/z**2,(1-r)/4),
 ('AH/B',Av*H/B,s.Rational(1,4)),('alpha(1-H)/B',alv*(1-H)/B,s.Rational(1,4)),
 ('kappa/B',(kv-8*Dv*H*r)/B,r)]:
    regular('weak IH rho1 simultaneous '+name,expr.subs(end1,simultaneous=True),{z:0,tau:0,ee:0},target)

# The toric identities underlying the complete local IH Gram.
u,w,W,Z=s.symbols('u w W Z',real=True)
R0=(u*W-B*w*Z)/s.sqrt(B);R1=(z*W-B*u*w*Z)/s.sqrt(B)
R2=u*w*W-B/y*Z;T=(u**2-z)*Z;E=s.sqrt(B*H)*(w**2-1/y)*Z
zero('IC8 toric conductor',B*(u**2-z)*(w**2-1/y)*Z**2-T*R2+z*R0**2-u*R0*R1)
zero('IC8a positive T pivot',B/y*T**2+(u**2-z)*T*R2-u*(z*R0-u*R1)*(u*R0-R1))
zero('SL3 C packet relation',s.sqrt(H)*R2-s.sqrt(B*H)*w*R0-s.sqrt(B)*E)

# Original homogeneous reciprocal-Bob parent, including all state normals.
X,Y,n,m,C,A,D,al=s.symbols('X Y n m C A D al',real=True)
du=1+u**2;dw=1+w**2
ampM=W-w*X+u*Y-u*w*Z;ampJ=u*w*W+u*X-w*Y-Z
parent=4*n*du*dw*X**2+4*m*ampM**2+4*C*dw*(u*W-Y)**2+4*A*du*dw*Y**2+4*D*ampJ**2+4*al*du*(Y-w*Z)**2+2*(ampJ**2+du*(Y-w*Z)**2+dw*(X+u*Z)**2-du*dw*Z**2)
graph={X:-r*u*Z,Y:H*w*Z}
K=4*A*H**2+8*D*H*r+4*al*(1-H)**2+4*m*(H+r-1)**2+2*(H+r-1)**2+4*n*r**2
kk=K-4*B**2*m/z**2+8*z*n*r**2+8*A*H**2/y
graphgram=(4*B*C-kk*z/B)*R0**2+4*B*m/z**2*R1**2+2*(2*C+2*D+1)*R2**2+r**2*(4*n*w**2+2*(2*D+2*n+1))*T**2+H/B*(4*A*u**2+2*(2*A+2*C+2*D+1))*E**2+kk*u/B*R0*R1+(-4*r*(2*D+1)+kk/B)*R2*T-4*s.sqrt(H/B)*(2*C+2*D+1)*R2*E
src={n:nv,m:mv,C:Cv,A:Av,D:Dv,al:alv}
diff=s.Poly(s.expand(parent.subs(graph)-graphgram),u,w,W,Z)
assert set(diff.monoms()) <= {(1,1,1,1),(2,0,0,2),(0,2,0,2),(0,0,0,2)}
checks.append({'name':'IC7 residual has only four contact-jet monomials','passed':True})
# Reconstruct all four residual coefficients from exact native contact jets.
# This proves the identity modulo the native value/state/angle equations,
# without duplicating a huge rational elimination of those same equations.
gc=parent.subs(graph)
cp={u:s.sqrt(z),w:1/s.sqrt(y),W:B/s.sqrt(z*y),Z:1}
ca=s.cancel(s.diff(gc,W).subs(cp)*s.sqrt(y/z))
cb=s.cancel((s.diff(gc,u).subs(cp)-ca*B/(s.sqrt(z)*y))/(2*s.sqrt(z)))
cc=s.cancel((s.diff(gc,w).subs(cp)-ca*B/s.sqrt(y))*s.sqrt(y)/2)
cd=s.cancel(gc.subs(cp)-ca*B/y-cb*z-cc/y)
zero('IC7 exact reconstruction from native contact jets',diff.as_expr()-(ca*u*w*W*Z+cb*u**2*Z**2+cc*w**2*Z**2+cd*Z**2))
PX=(2*D+1)*B/y-2*m*B/z*w**2-2*n*r*(1+u**2)*(1+w**2)-r*(2*D+1)*u**2-2*D-r+w**2*(1-H-r-2*D*H-2*m*(H+r-1))
LX=(2*D+1)*u*R2-2*m*s.sqrt(B)/z*w*R1+u*Z*PX
zero('IC10 full X normal',s.diff(parent,X).subs(graph)/4-LX)
f=2*m-2*C*(1+w**2)-(2*D+1)*w**2
PY=B*f+2*A*H*(1+u**2)*(1+w**2)+2*C*H*(1+w**2)+2*D*H*w**2+2*D*r*u**2+2*D-2*al*(1-H)*(1+u**2)+(2*m+1)*(H+r-1)*u**2+H*(w**2+1)
zero('IC12 full Y normal',s.diff(parent,Y).subs(graph)/4-(s.sqrt(B)*f*R0+w*Z*PY))
zero('IC normal X contact remainder is native state jet',PX.subs(cp)-s.diff(parent,X).subs(graph).subs(cp)/(4*s.sqrt(z)))
zero('IC normal Y contact remainder is native state jet',PY.subs(cp)-s.diff(parent,Y).subs(graph).subs(cp)*s.sqrt(y)/4)

# Independent UI leading model and its full matched Schur Hessian.
aa,bb,rr,vv,xx,zz,uu=s.symbols('aa bb rr vv xx zz uu',positive=True)
az=aa**2; ac0=1+az*rr; gg=(1+az)*bb-1+rr
de=1+(1+az)*bb+az*(1+az)*bb**2+2*az*(1+az)*bb*rr+az**2*rr**2+2*az*rr
n0=((bb+rr)*ac0**2-rr*de)/(2*rr*de)
m0=(1-rr)*ac0**2/(2*gg*de)
al0=az**2*(1+az)*bb*(bb+rr)**2/(2*(1-az*bb)*de)
h0=az*(1+az)*(bb+rr)*bb/(2*de)
aw=ac0**2*(bb+rr)/(2*az*bb*de)
hh,aw0=s.symbols('hh aw0',positive=True)
lead=4*n*(1+uu**2)*xx**2+4*m*(vv*(1+uu**2)-xx-uu*zz)**2+4*hh*vv**2*(uu*xx-zz)**2+4*aw0*(1+uu**2)*uu**2+(4*al+2)*(1+uu**2)*(vv*uu-zz)**2+2*(xx+uu*zz)**2-2*(1+uu**2)*zz**2
lsrc={n:n0,m:m0,al:al0,hh:h0,aw0:aw}
contact={uu:aa,vv:1,xx:-rr/bb,zz:1/(aa*bb)}
for var in [None,xx,zz,vv,uu]:
    ex=lead if var is None else s.diff(lead,var)
    zero('UI leading contact '+str(var),ex.subs(contact).subs(lsrc,simultaneous=True))
hess=s.hessian(lead,(xx,zz,vv,uu))/2
hess=hess.subs(contact).subs(lsrc,simultaneous=True).applyfunc(s.cancel)
hard=hess[:2,:2]; cross=hess[:2,2:]; angular=hess[2:,2:]
harddet=s.cancel(hard.det())
adj=s.Matrix([[hard[1,1],-hard[0,1]],[-hard[1,0],hard[0,0]]])
schur=(angular-cross.T*adj*cross/harddet).applyfunc(s.cancel)
zero('UI4 full Schur uu',schur[1,1]-8*(1-rr)*ac0**2/(az*bb*de))
zero('UI full Schur Vu',schur[0,1]-4*ac0**2*(az*bb+az*rr+2*rr-2)/(aa*bb*de))
zero('UI full Schur VV',schur[0,0]-8*ac0*(az**2*bb**2+az**2*bb*rr+az*bb**2+2*az*bb*rr-az*bb-rr+1)/(bb*de))
zero('UI5 full Schur determinant',schur.det()-16*(bb+rr)**2*ac0**3*(3*az+4-(az+2)**2*rr)/(bb**2*de**2))

# UA formulas checked independently in unnormalized source-weight variables.
U=s.symbols('U',positive=True)
du0=1+U; ns=2*n+1;ts=2*m+2*al+1
P=s.Matrix([[4*n*du0+4*m+2,(4*m+2)*s.sqrt(U)],[(4*m+2)*s.sqrt(U),(4*m+2)*U+4*al*du0]])
ell=s.Matrix([s.sqrt(U),-1]); bc=s.Matrix([-4*m*du0,-du0*s.sqrt(U)*(4*m+4*al+2)])
c2=4*m*du0**2+(4*al+2)*du0*U
pinv=P.inv();qa=al*(2*m+ns)/(ns*ts);tm=m/ts;tam=al*m/ts;iN=1/ns;J0=qa+U*(1-iN)/2
dv=(ell.T*pinv*ell)[0]; ev=(ell.T*pinv*bc)[0];hv=c2-(bc.T*pinv*bc)[0]
zero('UA12 det P normalization',P.det()/(ns*ts)-8*du0*J0)
zero('UA12 d',dv-(U*iN+qa/al)/(4*J0))
zero('UA12 e',ev-s.sqrt(U)*du0/(2*J0))
zero('UA12 h',hv-du0*(4*tam-U)/J0)
# The preceding equalities justify using their reduced expressions here.
# Multiplying the claimed minimizer by P is more economical than expanding
# another inverse expression; P is positive definite on its stated domain.
dv=(U*iN+qa/al)/(4*J0);ev=s.sqrt(U)*du0/(2*J0);hv=du0*(4*tam-U)/J0
zero('UA12 dh+e^2',dv*hv+ev**2-du0*(U+2*tm)/(2*J0))
pbc=s.Matrix([iN*(U-4*tam)/(2*J0),-s.sqrt(U)*(4*tam*iN+(1-iN)*U+1)/(2*J0)])
zero('UA12 minimizer equation1',(P*pbc-bc)[0])
zero('UA12 minimizer equation2',(P*pbc-bc)[1])
fdom=1+rr*(az-1)-bb*rr*az*(az+1)-rr**2*az*(az+2)
tsym=s.symbols('tsym')
bs=tsym*bb+(1-tsym)/(1+az); rs=tsym*rr
zero('UI valid real interior approximation',fdom.subs({bb:bs,rr:rs},simultaneous=True)-((1-tsym)*(1+tsym-tsym*rr)+tsym**2*fdom))

out.write_text(json.dumps({'sympy':s.__version__,'scope':'Fresh algebra only; analytical audit separate','checks':checks,'complete':True,'all_passed':all(x['passed'] for x in checks)},indent=2)+'\n')
print('TOTAL',len(checks),'ALL PASSED',flush=True)
