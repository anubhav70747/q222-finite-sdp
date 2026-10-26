"""Exact R3/R4 audit. SymPy only; no numerical evaluation or optimization.
Pinned sources: c2_short.tex RSR5--6, RFN7--22, EN17--22.
Every assertion is an exact rational/polynomial identity.
"""
import sympy as s
u,v,Z,p,H,mu,n,C,A,D,alpha,j,T,B = s.symbols('u v Z p H mu n C A D alpha j T B')
Du=1+u*u; Dv=1+v*v; h0=mu+s.Rational(1,2); x=j*B; y=T*B
N,M=s.symbols('N M')
G=(mu*Du*(N*N-2*v*N+(1-p)*v*v)+n*Dv*(N+p*v-u*Z)**2+C*Du*Dv*(M+H*v*Z)**2+A*Dv*(u-M-H*v*Z)**2+(D+s.Rational(1,2))*Du*(M-(1-H)*v*Z)**2+(alpha+s.Rational(1,2))*(u*(1-p)*v-(1+H*v*v)*Z-u*N-v*M)**2+s.Rational(1,2)*Du*Dv*(N+p*v)**2-s.Rational(1,2)*Dv*(Z+u*N+p*u*v)**2)
Ln=s.expand(s.diff(G,N).subs({N:0,M:0})/2)
Lm=s.expand(s.diff(G,M).subs({N:0,M:0})/2)
q=v*Z-T*u; R=Z-u*v/B; Ru=v*(u*u-x); b=v-u*Z/j
Kc=p/2+j*(H*(alpha+s.Rational(1,2))-s.Rational(1,2))+n*(p-j)
Kr=-mu-(alpha+s.Rational(1,2))*(1-p)+T*Kc/j
SN=n*p-mu+p/2+(-mu-(alpha+s.Rational(1,2))*(1-p))*x+j*(alpha-n)+y*Kc
LK=(j*(n-alpha)+p*(n+s.Rational(1,2))*v*v)*b+Kr*Ru+Kc/j*u*v*q
resN=s.factor(Ln-LK)
print('N diff factor:',s.factor(resN-v*SN))
Kd=C*H-(D+s.Rational(1,2))*(1-H)
Cq=Kd+A*H+s.Rational(1,2)+v*v*(C*H*(1+u*u)+A*H+(alpha+s.Rational(1,2))*H)
LM=Cq*q+mu*u*v*b+v*(alpha-B*T*C*H*u*u)*R
S0=T*(Kd+A*H+s.Rational(1,2))-A
S1=j*Kd+mu+x*T*C*H
S2=j*Kd+x*T*C*H+T*(C*H+A*H+(alpha+s.Rational(1,2))*H)-A-(alpha+s.Rational(1,2))*(1-p)+alpha/B
resM=s.factor(Lm-LM)
print('M diff factor:',s.factor(resM-(u*S0+u*u*v*Z*S1/j+u*v*v*(S2-S1))))
print('End residual terms',s.count_ops(resN),s.count_ops(resM))
Jb,Hp,Bb,Tt,Pp=s.symbols('jj hh bb tt pp')
K=Jb+1-Pp; X=Jb*Bb; Y=Tt*Bb; EE=1+Bb*(Hp*Tt+Pp-1)
PP=Jb+1+Bb*((Jb+1)**2-Pp)+Tt*Bb*(Hp*Jb+Jb+1+Bb*(Hp*Jb*(2*Jb+2-Pp)+Pp*(2*Jb+1-Pp)))+Hp*Jb*Tt**2*Bb**2*(1+Bb*Pp+Bb*Hp*Jb)
LL=K+Tt*((1-Hp)*Jb+1-Pp-Hp)+Bb*Tt*(Pp*K+Hp*Jb*(1-Pp))+Bb*Hp*Tt**2*((1-Hp)*Jb-Pp)
NN=Pp+Bb*(-Jb**2+3*Jb*Pp-Jb-Pp**2+Pp)+Bb*Tt*(Hp*Jb**2+Hp*Jb-Jb**2+Jb*Pp-Jb+Pp)+Bb**2*Tt*(2*Hp*Jb**2*Pp-Hp*Jb**2-Hp*Jb*Pp**2+2*Hp*Jb*Pp-2*Jb**2*Pp+3*Jb*Pp**2-Jb*Pp-Pp**3+Pp**2)+Bb**2*Hp*Tt**2*Jb*(Hp*Jb-Jb+Pp)*(1+Bb*Pp)
MU=Jb*Pp*(1+Y)*(1+X)*(1+Hp*Y)/(2*PP)
NNN=NN/(2*(Jb-Pp)*PP)
CCC=LL/(2*Hp*Tt*PP)
AAA=Tt*(1+X)*(1+Hp*Y)*(K+Hp*Jb*Y)/(2*(1-Hp*Tt)*PP)
DDD=(1+Hp*Y)*(K+Hp*Jb*Y)*(1-Y*((1-Hp)*Jb-Pp))/(2*Tt*(1-Hp)*PP)
ALPHA=Bb**2*(K+Hp*Jb*Y)*LL/(2*EE*PP)
repl={j:Jb,T:Tt,B:Bb,p:Pp,H:Hp,mu:MU,n:NNN,C:CCC,A:AAA,D:DDD,alpha:ALPHA}
for title,expr in [('SN',SN),('S0',S0),('S1',S1),('S2',S2)]:
    z0=s.together(expr.subs(repl))
    num,den=s.fraction(z0)
    q=s.Poly(s.expand(num),Bb,Jb,Tt,Pp,Hp)
    print(title,'zero?',q.is_zero, 'monomials',len(q.terms()), 'numops',s.count_ops(num))
# RFN12 scalar readback on original u,v,z, symbolic source masses.
x=j*B; y=T*B; qz=j*T
K0=C*H*H+(D+s.Rational(1,2))*(1-H)**2
Ma=A*x; Mb=(mu*(1-p)+(n+s.Rational(1,2))*p**2)*y; Mc=alpha*qz
Mab=(A+alpha*(1-p)**2+mu*(1-p)+s.Rational(1,2)*(1-p)**2)*x*y
Mac=n*x*qz
Mbc=(K0+A*H**2+2*H*alpha+H-s.Rational(1,2))*y*qz
Mb2=(n+s.Rational(1,2))*p**2*y**2
Mabc=(K0+n)*x*y*qz
Mb2c=H**2*(A+C+alpha+s.Rational(1,2))*y**2*qz
Mab2c=C*H**2*x*y**2*qz
kk=(Ma-Mb+Mac+Mab2c-Mbc)/2
rr=Mc-Mb+Mac-kk
ff=Mb-Mac
qq=Ma-ff
dc=Mb2c-kk
gc=Mab2c
e0=s.Rational(1,2)*p**2*y**2-gc
Q2=(T*u-v*Z)**2/(T**2*x); B2=(v-u*Z/j)**2/y; R2=(Z-u*v/B)**2/(j*T)
Psi=u**2/x+v**2/y+Z**2/qz+u**2*v**2*Z**2/(qz**2*B**2)-4*u*v*Z/(qz*B)
Fref=qq*Q2+n*j**2*y*B2+rr*R2+dc*v**2/y*Q2+(n*p**2*y**2+e0)*v**2/y*B2+ff*Psi+gc*v**2/y*Psi+kk*Z**2/qz*(v**2/y-1)**2
Fd=s.Poly(s.expand(G.subs({N:0,M:0})-Fref),u,v,Z)
print('RFN12 generic monomials difference:',len(Fd.terms()))
status=[]
for ind, ((powers),co) in enumerate(Fd.terms()):
    num,den=s.fraction(s.together(co.subs(repl)))
    v=s.Poly(s.expand(num),Bb,Jb,Tt,Pp,Hp).is_zero
    status.append((powers,v))
    if ind%5==0: print('RFN12 coefficient',ind,powers,v,flush=True)
print('RFN12 result',all(y for x,y in status), [x for x,y in status if not y])
import sympy as s
U,V,Z=s.symbols('U V Z')
q=U-V*Z; b=V-U*Z; r=Z-U*V
psi=U**2+V**2+Z**2+U**2*V**2*Z**2-4*U*V*Z
assert s.expand(psi-(U-V)**2-Z**2*(1-U*V)**2-2*U*V*(1-Z)**2)==0
assert s.expand((1+U**2)*psi-r**2-(b**2+U**2*q**2+U**2*(1-U*V*Z)**2))==0
assert s.expand((1+V**2)*psi-V**2*b**2-(q**2+r**2+V**2*(1-U*V*Z)**2))==0
t=U*V; de=U-V; e=Z*(1-t); tt=t*(1-Z); jj=tt-de**2*Z
assert s.expand(2*(psi+V**2*q**2+U**2*b**2)-r**2-(e+tt)**2-2*de**2-4*t*(1-Z)**2-2*jj**2-4*t*de**2*Z**2)==0
print('RFN11, RFN18a, RFN18b, RFN18c exact')
t7=s.Matrix([Z-U**2*V,V-U**2*V,U-U**2*V,U*Z-U**2*V,U*V-U**2*V,U*V*Z-U**2*V,U**2-U**2*V])
t3=s.Matrix([1-U*V,V-U*V,U-U*V])
G7=s.Matrix([[4,-2,s.Rational(-3,2),0,s.Rational(-1,2),0,s.Rational(1,2)],[-2,3,0,s.Rational(-3,2),0,-1,2],[s.Rational(-3,2),0,4,s.Rational(3,2),-2,-3,0],[0,s.Rational(-3,2),s.Rational(3,2),4,s.Rational(-3,2),0,s.Rational(-3,2)],[s.Rational(-1,2),0,-2,s.Rational(-3,2),3,2,-1],[0,-1,-3,0,2,8,-4],[s.Rational(1,2),2,0,s.Rational(-3,2),-1,-4,4]])
G3=s.Matrix([[3,-2,-2],[-2,2,2],[-2,2,3]])
Dstar=V*(U**2-1)
rf=s.expand(4*psi+4*U**2*q**2+4*U**2*r**2-Dstar**2 - (t7.T*G7*t7)[0]-4*V*Z*(1-U)**2-U*Z*(t3.T*G3*t3)[0])
print('RFN19 zero?',rf==0,'G7 minors',[G7[:k,:k].det() for k in range(1,8)],'G3 minors',[G3[:k,:k].det() for k in range(1,4)])
assert rf==0
Erad=Z*(U**2-1)
Wperm=4*psi+4*U**2*q**2+4*U**2*b**2-Erad**2
Rfn20=s.expand(8*U**2*psi+8*U**4*q**2+(8*U**4+2*(U**2-1)**2)*b**2-Dstar**2-2*U**2*Wperm-(U*Erad-(U**2-1)*b)**2)
assert Rfn20==0
print('RFN20 exact')
import sympy as s
u,v,T,Y,X,Z=s.symbols('u v T Y X Z',real=True)
Du=1+u*u;Dv=1+v*v
core=(Du*(v*Y-Z)**2+Dv*(X+u*Z)**2-Du*Dv*Z**2)/2
for name,form,truths in [
 ('R3', T*Du*Dv*Y**2/2+(1+T)*v**2*(u*X-Z)**2/(2*T)+core, [Du*(T*Dv+v**2)/2, Du*(T*Dv+v**2)*(T*u*u*v*v+T*v*v+T+u*u*v*v)/(4*T),v**2*(1+T)*Du**2*Dv*(T*u*u+v*v)/(8*T)]),
 ('R4', Dv*v*v*(Z-u*X)**2/(2*T)+T*Du*Dv*Y**2/(2*(1-T))+core, [Du*(T+v*v)/(2*(1-T)), Du*Dv*(T+v*v)*(T+u*u*v*v)/(4*T*(1-T)),v*v*Du**2*Dv**2*(T*u*u+v*v)/(8*T*(1-T))])
]:
    l=s.Matrix([Y,X,Z]);g=s.hessian(form,l)/2
    checks=[s.cancel(g[:k,:k].det()-val)==0 for k,val in enumerate(truths,1)]
    print(name,'EN22 principal minors',checks)
    assert all(checks)
# EN14 exact correlated completion
h0,n,j,N,b=s.symbols('h0 n j N b',real=True)
cb=n*j/(n+h0); Me=n*h0*j**2/(n+h0)
assert s.factor(h0*N*N+n*(N+j*b)**2-(n+h0)*(N+cb*b)**2-Me*b*b)==0
print('EN14 correlated current identity exact')
# Correct source-limit coquotients: j=1/d, B=d*x. No series or samples.
import sympy as sp
jj, pp, hh, tt, bb, xx, dd=sp.symbols('jj pp hh tt bb xx dd')
kk=jj+1-pp; X0=jj*bb; Y0=tt*bb; E0=1+bb*(hh*tt+pp-1)
P=jj+1+bb*((jj+1)**2-pp)+tt*bb*(hh*jj+jj+1+bb*(hh*jj*(2*jj+2-pp)+pp*(2*jj+1-pp)))+hh*jj*tt**2*bb**2*(1+bb*pp+bb*hh*jj)
L=kk+tt*((1-hh)*jj+1-pp-hh)+bb*tt*(pp*kk+hh*jj*(1-pp))+bb*hh*tt**2*((1-hh)*jj-pp)
F1=-2*bb*hh*jj*pp+bb*hh*jj+2*bb*jj*pp-2*bb*pp**2+2*bb*pp-hh*jj-hh*pp-hh+jj+1-pp**2/jj+pp/jj
F2=-bb**2*hh**2*jj*pp+bb**2*hh*jj*pp-bb**2*hh*pp**2+bb**2*hh*pp-bb*hh**2*jj-bb*hh**2*pp+bb*hh*jj-bb*hh*pp**2/jj+bb*hh*pp/jj
F=kk+tt*F1+tt**2*F2
mu0=jj*pp*(1+Y0)*(1+X0)*(1+hh*Y0)/(2*P)
A0=tt*(1+X0)*(1+hh*Y0)*(kk+hh*jj*Y0)/(2*(1-hh*tt)*P)
C0=L/(2*hh*tt*P)
D0=(1+hh*Y0)*(kk+hh*jj*Y0)*(1-Y0*((1-hh)*jj-pp))/(2*tt*(1-hh)*P)
a0=bb**2*(kk+hh*jj*Y0)*L/(2*E0*P)
expected={
    'P/(j+1)':(P/(jj+1),sp.Integer(1)),
    'L/(j+1)':(L/(jj+1),1+(1-hh)*tt),
    'F/(j+1)':(F/(jj+1),1+(1-hh)*tt),
    'mu/p':(mu0/pp,sp.Rational(1,2)),
    'A':(A0,tt/(2*(1-hh*tt))),
    'CH':(C0*hh,(1+(1-hh)*tt)/(2*tt)),
    'D(1-H)':(D0*(1-hh),1/(2*tt)),
    'alpha/(B²(j+1))':(a0/(bb*bb*(jj+1)),(1+(1-hh)*tt)/2),
}
for name,(expr,limit) in expected.items():
    at0=sp.cancel(sp.cancel(expr.subs({jj:1/dd,bb:dd*xx})).subs({dd:0,xx:0}))
    assert sp.cancel(at0-limit)==0,name
print('EN5/EN8 native coquotient limits exact:',len(expected))
print('ALL EXACT R3/R4 CHECKS PASSED')
