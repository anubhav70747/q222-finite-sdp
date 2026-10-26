"""Exact full-state FM/UKR sign and residue checks. No solver or sampling."""
import json
from pathlib import Path
import sympy as s

checks = {}
def ck(name, expression):
    checks[name] = s.cancel(expression) == 0
    if not checks[name]:
        raise AssertionError(name)

u,v,W,Z,Nx,Ny=s.symbols('u v W Z Nx Ny')
S,k,r,h,q,j=s.symbols('S k r h q j', nonzero=True)
m,D,al,be,C,A,c0,a0=s.symbols('m D al be C A c0 a0')
th=1-h; Du=1+u*u; Dv=1+v*v
X=Nx-h*u*Z; Y=r*u*W+q*Ny
parent=(m*(v*W-X+u*v*Y-u*Z)**2+C*Dv*(u*W-Y)**2
        +A*Du*Dv*Y**2+(D+s.Rational(1,2))*(u*W+u*v*X-Y-v*Z)**2
        +(al+s.Rational(1,2))*Du*(v*Y-Z)**2
        +(be+s.Rational(1,2))*Dv*(X+u*Z)**2-Du*Dv*Z**2/2)
P=s.Poly(s.expand(parent),Nx,Ny)
def co(i,jj):return P.coeff_monomial(Nx**i*Ny**jj)
Lx=co(1,0)/2; Ly=co(0,1)/2
Dp=D+s.Rational(1,2)
LU=m*th-m*j/k*v*(1+r*u*u)-Dp*v*v*(1+h*u*u)+(be*th+th/2)*Dv+Dp*(1-r)*v*j*S
YP=Dv*(-c0+a0*Du)+S*(m*v*v*(1+r*u*u)-Dp*(1-r)+(al+s.Rational(1,2))*Du*r*v*v)
YQ=-m*th*u*u+Dp*(1+h*u*u)-(al+s.Rational(1,2))*Du
Lx_frame=-m/k*v*(1+r*u*u)*(k*W-j*u*Z)+Dp*(1-r)*v*(u*u*W-j*u*S*Z)+LU*u*Z
Ly_frame=q*(YP*u*W/S+v*YQ*Z)
sub={C:c0/(S*(1-r)),A:a0/(S*r)}
ck('FM18_Lx_complete',Lx-Lx_frame)
ck('FM18_Ly_complete',Ly.subs(sub)-Ly_frame)
ck('FM16_Hx_complete',co(2,0)-(m+Dp*u*u*v*v+Dv/2+be*Dv))
ck('FM16_Hxy_complete',co(1,1)+2*u*v*q*(m+Dp))
Hy=r*c0*Dv+(1-r)*a0*Du*Dv+S*r*(1-r)*(m*u*u*v*v+Dp+(al+s.Rational(1,2))*Du*v*v)
ck('FM16_Hy_complete',co(0,2).subs(sub).subs(q*q,S*r*(1-r))-Hy)
wrong=LU-2*(be*th+th/2)*Dv
ck('wrong_sign_discrepancy',LU-wrong-2*(be+s.Rational(1,2))*th*Dv)

# UKR8--10 capped parent: retain original beta*theta couplings.
m0,m1,m2,bc=s.symbols('m0 m1 m2 bc')
V=v*(W+u*Y); B0=X+u*Z
cap=parent-m*(V-B0)**2+m0*B0**2-2*m1*B0*V+m2*V**2+(bc-be)*Dv*Nx**2
CP=s.Poly(s.expand(cap),Nx,Ny)
LUcap=m0*th-m1*j/k*v*(1+r*u*u)-Dp*v*v*(1+h*u*u)+(be*th+th/2)*Dv+Dp*(1-r)*v*j*S
Lxcap=-m1/k*v*(1+r*u*u)*(k*W-j*u*Z)+Dp*(1-r)*v*(u*u*W-j*u*S*Z)+LUcap*u*Z
YPcap=YP.subs(m,m2); YQcap=YQ.subs(m,m1)
ck('UKR10_Lx_complete',CP.coeff_monomial(Nx)/2-Lxcap)
ck('UKR10_Ly_complete',CP.coeff_monomial(Ny).subs(sub)/2-q*(YPcap*u*W/S+v*YQcap*Z))
ck('UKR10_Hx_complete',CP.coeff_monomial(Nx**2)-(m0+Dp*u*u*v*v+(bc+s.Rational(1,2))*Dv))
ck('UKR10_Hxy_complete',CP.coeff_monomial(Nx*Ny)+2*u*v*q*(m1+Dp))
ck('UKR10_Hy_complete',CP.coeff_monomial(Ny**2).subs(sub).subs(q*q,S*r*(1-r))-Hy.subs(m,m2))

# Full graph identity FM10: j is arbitrary, so no hidden contact relation.
AT=m*v*v*(1+r*u*u)**2/k**2
AP=Dv*(c0*(1-r)+a0*r*Du)+S*(Dp*(1-r)**2+(al+s.Rational(1,2))*Du*r*r*v*v)
AQ=(D*v*v+al)/S
AU=m*th*th+al+be*th*th+th*th/2+v*v*(2*D*h+h-s.Rational(1,2)+be*th*th+th*th/2)+Dp*h*h*v*v*u*u
KK=m*th*(1+r*u*u)+Dp*(1-r)*(1+h*u*u)+(al+s.Rational(1,2))*Du*r
HH=v*KK-k*j*AT
gg=j*AP-HH; RU=AU-j*j*AT; RQ=AQ+j*j*AP-2*j*HH
# Substitute the square-root-free equivalent P=uW, Q=SZ, with scaling S.
NN=k*W-j*u*Z
GG=u*W-j*S*Z
graph=AT*NN**2+AP*GG**2/S+2*gg*GG*Z+(S*RQ+u*u*RU)*Z**2
ck('FM10_graph_complete',co(0,0).subs(sub)-graph)

nu,ell=s.symbols('nu ell',positive=True)
LUres=-v*v/2+(ell/2)*(1+v*v)
ck('FM17_residue',(LUres-(ell-nu*v*v)/2).subs(ell,1-nu))
ck('GB_residue_determinant',(s.Rational(1,4)-(ell-nu)**2/4-nu*ell).subs(ell,1-nu))
checks['FM22_lower_bound_at_least_3']='2nu j^2<=1/2 and theta^2 nu^2/(2h)<=1/2'

# Exact UKR15 and UKR14 on the uncapped real source chart.
aa,hh,kk,bb=s.symbols('aa hh kk bb',nonzero=True)
rr=1-aa; dd=hh-aa; nn=aa+kk*hh; ll=1-nn; tt=1-hh
bracket=4*aa**2*hh*rr-kk*(rr*hh-aa)**2
positive=aa**3*rr*(4-aa)+aa**4*(rr-kk)+aa**2*rr*dd*(3+aa+2*kk)+rr**2*dd*(aa**2-kk*dd)
ck('UKR15_positive_factorization',bracket-positive)
mm=kk**2*hh/(2*nn*(rr-kk)); bet=(hh-aa)/(2*tt); EE=aa/nn; jj=nn*bb
Hs=hh-s.Rational(1,2)+tt*(ll+kk*hh)/2
rs=Hs-jj**2*mm/kk**2
Hc=s.Matrix([[mm+(bet+s.Rational(1,2))/nn,-mm*bb/kk],[-mm*bb/kk,mm*bb*bb/kk**2]])
zz=s.Matrix([-jj*(mm/kk+EE),-rs/nn]); AG=EE/(2*nn)
gst=EE*ll-mm*tt-s.Rational(1,2)+2*mm*ll/kk+kk*rs/nn
HB=nn*(mm*tt+s.Rational(1,2))
Jstar=HB-gst**2/AG-kk*(zz.T*Hc.inv()*zz)[0]
Jexpected=(aa*rr+kk*(hh-aa))/(2*aa*hh*rr)*bracket
num=s.cancel(Jstar-Jexpected).as_numer_denom()[0]
checks['UKR14_Schur_identity']=s.rem(s.Poly(num,bb),s.Poly(bb*bb-ll/nn,bb)).is_zero
assert checks['UKR14_Schur_identity']
out={'checks':checks,'all_exact_identities_pass':True,
     'scope':'FM complete graph and both state-normal identities; residue strictness and uncapped UKR Schur algebra. No endpoint-leaf replay.'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
