"""Finite reconstruction of all scalar/first-normal inputs to bounded-b endpoint.
Uses only the seven printed factored packets and monic polynomial divisions.
"""
import json,time
from pathlib import Path
import sympy as s
started=time.time(); out={}
e,k,b,a,t,u,v,X,Y,Z,U=s.symbols('e k b a t u v X Y Z U')
K=k*k; eta=e*e*b*b; zeta=k*k*b*b; Vc=K*eta*zeta
Rc=1+K*(1-a*eta); T0=1+K-t*K*zeta
Lc=1-a*eta; Hc=1+(1-a)*eta; Qc=1-K*zeta+a*Vc
P4=1+K+eta*K*(a-K-1)-zeta*K*t+Vc*(-zeta*K*t*t+2*K*t+a*(1-t)+t)
Pal=1+K-a*eta*K+zeta*((1+K)-(2+K)*t)+Vc*(-zeta*t*t+a*(t-1)+t)
Pd=(1+K)**2-a*eta*K*K+zeta*K*((K+1)-(K+2)*t)+Vc*(-zeta*K*t*t+K*a*(t-1)+K*t+K+1)
c2=2*(1+K)*(1-t*zeta)*Lc*(1+Vc)*Hc*Qc
cb=2*(1+K)*Lc*(1+(1-t)*zeta)*Hc*Qc
c3=2*(1+K)*Lc*Hc*Qc*(1-K*eta+t*Vc)
c5=2*(1+K)*Lc*Hc*Qc*(1+Vc)
c4=2*Qc*P4; ca=2*Lc*Hc*Pal; cd=4*Lc*Qc*Pd
vg=e*k*b*b*u; Xg=-u*Z+t*vg; Yg=Lc*u
r2=v-X+u*v*Y-u*Z; rb=X+u*Z; r3=u-Y
r5=u+u*v*X-Y-v*Z; ra=v*Y-Z
q0=vg; z0=eta*u
G=[4*T0*e*e*c2*((r2+t*q0)**2-t*Rc*q0*q0),
4*T0*e*e*cb*(1+v*v)*((rb+(Rc-t)*q0)**2-Rc*(Rc-t)*q0*q0),
4*Rc*c3*(1+v*v)*((r3+(T0-a)*z0)**2-T0*(T0-a)*z0*z0),
4*Rc*c5*((r5-a*z0)**2-a*T0*z0*z0),
4*Rc*T0*eta*eta*c4*(1+u*u)*(1+v*v)*Y*Y,
4*Rc*T0*K*K*eta*zeta*ca*(1+u*u)*ra*ra,
2*Rc*T0*eta*cd*(1+u*u)*(ra*ra-(1+v*v)*Z*Z)]
N2=Yg-Y; S=r5+(T0-a)*z0
G.extend([-4*Rc*c3*(1+v*v)*N2*N2,-4*Rc*c5*S*S])
at={X:Xg,Y:Yg,v:vg}
eg=s.expand(sum(packet.subs(at,simultaneous=True) for packet in G))
def coeffdivide(expr,divisor,name):
    q,r=s.div(expr,divisor,*[e,k,b,a,t,u,Z])
    if r!=0: raise AssertionError(name+' not divisible')
    return q
scaled=coeffdivide(eg,e*e*b**4,'graph')
zp=s.Poly(scaled,Z)
AA=coeffdivide(zp.coeff_monomial(1),e*e,'scalar constant')
BB=coeffdivide(zp.coeff_monomial(Z),e*k,'scalar linear')
CC=coeffdivide(zp.coeff_monomial(Z*Z),k*k,'scalar quadratic')
def even_U(expr):
    p=s.Poly(expr,u); val=0
    for (i,),cf in p.terms():
        if i%2: raise AssertionError('odd scalar parity')
        val+=cf*U**(i//2)
    return s.expand(val)
AA,BB,CC=map(even_U,(AA,BB,CC))
def md(expr,divisor,name):
    q,r=s.div(expr,divisor,U)
    if s.expand(r)!=0: raise AssertionError(name+' remainder')
    return s.expand(q)
NN=md(AA+BB+CC,(U-K)**2,'N')
DD=md(BB+2*CC,U-K,'D')
R0=md(CC+K*DD+K*K*NN,U,'R0')
JJ=md(-DD-2*K*NN,U,'J')
A00=R0.subs(U,u*u); A01=u*JJ.subs(U,u*u); A11=NN.subs(U,u*u)
soft=e*b*b*u*(e-k*Z); qq=e*b*b*(e*u*u-k**3*Z)
out['scalar_multiply_back']=s.expand(eg-A00*soft**2-A01*soft*qq-A11*qq**2)==0
def bucket(expr,kr,es,name):
    p=s.Poly(s.expand(expr),k,e,b,a,t,u,Z); first=0; second=0
    for mon,cf in p.terms():
        powers=list(mon)
        if powers[0]>=kr: powers[0]-=kr; dest=0
        else:
            if powers[1]<es or powers[2]<es: raise AssertionError(name+' invalid bucket '+str(mon))
            powers[1]-=es; powers[2]-=es; dest=1
        trm=cf*s.prod(vv**ii for vv,ii in zip(p.gens,powers))
        if dest==0: first+=trm
        else: second+=trm
    out[name+'_residues']=[str(s.expand(cf).subs({e:0,k:0,u:0,Z:0})) for cf in (first,second)]
bucket(A00,4,2,'A00'); bucket(A01,2,1,'A01')
out['A11_residue']=str(s.expand(A11).subs({e:0,k:0,u:0,Z:0}))
for var,div,name in [(X,e*e*b*b,'X'),(Y,e*b*b,'Y'),(v,e*e*b*b,'v')]:
    raw=s.expand(sum(s.diff(packet,var).subs(at,simultaneous=True) for packet in G))
    rr=coeffdivide(raw,div,'normal '+name)
    Cxi=0 if var!=v else 4*Rc*T0*cd*e*b*b*(1+u*u)
    rt=s.expand(rr-Cxi*Z*u*(e-k*Z))
    rp=s.Poly(rt,Z)
    if rp.degree()>1: raise AssertionError('quadratic normal remains')
    L0=even_U(coeffdivide(rp.coeff_monomial(1),e*u,name+' constant parity'))
    L1=even_U(coeffdivide(rp.coeff_monomial(Z),k*u,name+' linear parity'))
    Jxi=md(L0+L1,U-K,'J_'+name)
    Pxi=s.expand((L0-U*Jxi).subs(U,u*u)+Cxi*Z)
    Qxi=s.expand(u*Jxi.subs(U,u*u))
    out[name+'_multiply_back']=s.expand(rr-Pxi*u*(e-k*Z)-Qxi*(e*u*u-k**3*Z))==0
    bucket(Pxi-(16*k*(T0+u*u) if var==v else 0),2,1,'P_'+name)
    out['Q_'+name+'_residue']=str(Qxi.subs({e:0,k:0,u:0,Z:0}))
out['elapsed_seconds']=round(time.time()-started,3)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
