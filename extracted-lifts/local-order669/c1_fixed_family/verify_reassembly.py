"""Independent full-state symbolic check of the extracted 12-root Gram."""
import json
from pathlib import Path
import sympy as s
import build_core as B

here=Path(__file__).resolve().parent
kk,ll,tt,gg,q,v,wb,zb=s.symbols('kk ll tt gg q v wb zb',positive=True)
ss=s.symbols('s',nonnegative=True)
X,Y=s.symbols('X Y',real=True)
k=kk**2; sig=tt**2; lam=ll**2; rho=lam*k*k/sig
dG=gg**2; S0=dG+rho; q0=lam*k; L0=lam*sig; a0=k/sig
u=ss**2; w=v**2; W=wb**2; Z=zb**2
source={B.k:k,B.sigma:sig,B.lam:lam,B.rho:rho,B.dG:dG,B.S0:S0,
        B.q0:q0,B.L0:L0,B.a0:a0,B.q:q,B.s:ss,B.w:w}
b,p,m,D,C,a,ell,h,nu,R0,Rx,Bx,Ay,chi,KA,c=[B.v[x] for x in
      'b p m D C a ell h nu R0 Rx Bx Ay chi KA c'.split()]
rb=B.r; beta=B.beta
Nx=X+u*h*Z; Ny=(Y-b*w*Z)/q; R=w*Z/q-u*q*W
Gam=sig*Z-u*w*W; At=q*(u*u-k*k)*W; Tt=k*q*(u-k)*W
Bi=[s.sqrt(L0)*u*Z,s.sqrt(rho)*u*w*W,s.sqrt(lam)*k*k*q*W,s.sqrt(lam)*k*w*Z/q]
G=ll*k*ss*v*wb*zb
di=[x-G for x in Bi]
Us=s.sqrt(q0)*q*ss*(u-k)*W; Ks=s.sqrt(q0*k)*q*(u-k)*W
psi=s.Matrix([Nx,Ny,R,gg*Gam,At,Tt]+di+[Us,Ks])
phi=s.Matrix([X,Y,Z,W,u*Z,w*Z,u*W,u*w*W,u*u*W,ss**3*W,ss*W,ss*v*wb*zb])
Am=B.A.subs(source).applyfunc(s.simplify)
rootchecks=[s.expand(x)==0 for x in Am*phi-psi]
assert all(rootchecks)
UG=B.UG.subs(source).applyfunc(s.simplify)
SC=B.full_second.subs(source).applyfunc(s.simplify)
conductors=[s.cancel((UG.T*psi)[0]-u*Gam)==0,
            s.cancel((SC.T*psi)[0]-(R0-8*p*u*u)*q*w*Gam)==0]
assert all(conductors)
JC=Ny+b*R; JL=Ny-(1-b)*R
Eyy=2*q*q*((1+2*m)*u*u+(1+2*D)*w*w)
k1=2*(1+2*m)*(h-(1-b))**2+8*D*h*b+8*p*b*sig+4*beta*(1-h)**2
ERR0=q*q*(8*b*D+8*b*p*sig+4*b*sig*(2*C+2*D+1)*(rb-1))
ERRrest=q*q*(k1*u*u+2*b*b*(2*D+1)*w*w)
EyR=B.EyR.subs(source)
Pbeta=u*Z-k*k/sig*w*W
Scal=L0*Pbeta**2+rho*(u*u-k*k)*(w*w-sig*q*q)*W*W
Ccir=sum(x*x for x in Bi)-4*lam*k*k*u*w*W*Z
Cplus=Ccir-sum(x*x for x in di)
FA=u*(Z-a0*w*W)
epsC=2/lam+u+k; epsG=2*a0*u+(u+k)/L0
Delta=epsC*Ccir+u*u*R*R+epsG*Gam*Gam-At*At
qc=s.Rational(1,2)-c*dG*epsC/S0
target=(2+(2+4*m)*w*w)*Nx*Nx+(4*a+2*q*q*(2*C+1))*JC*JC+4*ell*JL*JL
target+=Eyy*Ny*Ny+(ERR0+ERRrest)*R*R+EyR*Ny*R+dG*Gam*Gam
target+=-4*nu*Nx*u*Gam+(R0-8*p*u*u)*JC*q*w*Gam
target+=-4*u*q*w*(2*m+2*D+1)*Nx*Ny+u*q*w*Rx*Nx*R
target+=q*w*Bx*Nx*At+u*Ay*Ny*At+KA*At*At+chi*u*At*R+Scal
retained=(2+4*m)*w*w*Nx*Nx+(4*a+2*q*q*(2*C+1)-1)*JC*JC+(4*ell-1)*JL*JL
retained+=Eyy*Ny*Ny+ERRrest*R*R+L0*FA*FA/2
retained+=c*dG*Delta/S0+(KA+c)*(At*At-Tt*Tt)/2+qc*Cplus
Q=B.Q.subs(source).applyfunc(s.simplify)
diff=s.expand((psi.T*Q*psi)[0]+retained-target)
res=s.cancel(diff)
assert res==0
report={'source_root_matrix_rows_checked':12,'independent_complete_normal_conductors_checked':2,
 'full_SP6_state_polynomial_reassembly':True,'state_variables':['X','Y','W','Z'],
 'assumptions':'Positive k, sigma, lambda, q; dG>=0 by polynomial extension; c arbitrary; no source sign, state ratio, or determinant-only substitution.',
 'scope':'Checks the complete twelve-root reassembly from SP6, not the upstream P1-to-P4 physical source identity.'}
(here/'REASSEMBLY_VERIFIED.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
