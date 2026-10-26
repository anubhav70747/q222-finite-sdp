#!/usr/bin/env python3
"""Exact reconstruction of the Q222 CA positive-unit-delta local base.

No SDP, sampled test, or positivity-by-sampling. Uses SymPy rational arithmetic.
Read together with Q222_CA_POSITIVE_UNIT_DELTA_FULL_PROOF_2026-10-26.md.
The finite monic-divison proof is mathematical and is explicitly given there;
this checks its identities and the complete bounded graded parent through CA6.
"""
import sympy as S


def zero(expression, title):
    z = S.cancel(expression)
    assert z == 0, f"{title}: {S.factor(z)}"


def original_native_balances():
    x, eta, zeta, a, t = S.symbols("x eta zeta a t")
    K = 1+x
    r = 1-a*eta
    theta = t*zeta
    g = 1-theta
    Hc = 1+(1-a)*eta
    Qc = 1-x*zeta+a*x*eta*zeta
    V = x*eta*zeta
    R = 1+x*r
    T = K-x*zeta*t
    M = R-t
    EN = 1-x*eta*g
    G = g+zeta
    P4original = 1+x+eta*x*(a-x-1)-zeta*x*t + V*(-zeta*x*t**2+2*x*t+a*(1-t)+t)
    PAoriginal = 1+x-a*eta*x+zeta*(K-(2+x)*t)+V*(-zeta*t**2+a*(t-1)+t)
    Pdoriginal = K**2-a*eta*x**2+zeta*x*(K-(x+2)*t) + V*(-zeta*x*t**2+x*a*(t-1)+x*t+x+1)
    P4compact = T*EN+x*(1-r)*G
    PAcompact = R*g*(1+V)+zeta*M*EN
    Pdcompact = 1+2*x*g+r*x**2*G+x*zeta+eta*x*(K*zeta+x*g*theta)
    for label,a1,b1 in [('P4',P4original,P4compact),('PA',PAoriginal,PAcompact),('Pd',Pdoriginal,Pdcompact)]:
        zero(a1-b1,'native polynomial '+label)
    c2=2*K*g*r*(1+V)*Hc*Qc
    cb=2*K*r*(1+(1-t)*zeta)*Hc*Qc
    c3=2*K*r*Hc*Qc*(1-x*eta+t*V)
    c5=2*K*r*Hc*Qc*(1+V)
    c4=2*Qc*P4original
    ca=2*r*Hc*PAoriginal
    cd=4*r*Qc*Pdoriginal
    identities=[
        cb-c2-zeta*c3,
        c3-c5+x*eta*cb,
        c2-g*c5,
        2*x*c2+2*(1+V)*r*c4-2*x*x*zeta*Qc*ca-Qc*cd,
        x*r*c2+c5+(1+V)*r*r*c4+x*Qc*Qc*ca
          +S.Rational(1,2)*(x*zeta*r*r-2*r-eta)*cd,
        x*R*c2+t*x*V*cb+a*V*c3+T*c5+K*V*r*r*c4
          -x*x*zeta*K*Qc*r*ca-K*S.Rational(1,2)*(Qc*r+eta)*cd,
    ]
    for i,idn in enumerate(identities): zero(idn,'SD9-SD10 '+str(i+1))
    print('PASS original SD1 polynomial coquotients and six SD9–SD10 balances')


def generic_original_packet_subtraction():
    qA,qB,N,M = S.symbols('qA qB N M')
    R,t,T,a = S.symbols('R t T a')
    p=t/R; h=a/T; Mg=R-t; J=T-a
    identities=[
        R/Mg*(N-(1-p)*qA)**2-t/Mg*N**2-(N**2-2*qA*N+(1-p)*qA**2),
        R/t*(N+p*qA)**2-Mg/t*N**2-(N**2+2*qA*N+p*qA**2),
        T/a*(M+h*qB)**2-J/a*M**2-(M**2+2*qB*M+h*qB**2),
        T/J*(M-(1-h)*qB)**2-a/J*M**2-(M**2-2*qB*M+(1-h)*qB**2),
    ]
    for i,idn in enumerate(identities):zero(idn,'SD2 minus SD3 packet '+str(i+1))
    A,r,u,W,Y,Du,Dv=S.symbols('A r u W Y Du Dv')
    zero(4*A*Du*Dv*Y**2-4*(A-1/a)*Du*Dv*(Y-r*u*W)**2
         -(4/a*Du*Dv*(Y-r*u*W)**2+4*A*r*Du*Dv*(2*u*W*Y-r*u*u*W*W)),
         'CA1 centered subtraction')
    print('PASS four original packet subtractions and CA1 centered A identity')


def full_raw_scalar_and_normals():
    u,v,w,N,M=S.symbols('u v w N M')
    b2,bb,b3,b5,L,a,d,alpha,r,p,h,delta,E=S.symbols('b2 bb b3 b5 L a d alpha r p h delta E')
    Du=1+u*u;Dv=1+v*v
    qA=v*Du*w;qB=v*(Du-u*N-p*u*Du*v*w)
    W=w+u*(M+h*qB)/Du
    Y=u*w-(M+h*qB)/Du
    LA=Y-r*u*W
    F=(b2*(N*N-2*qA*N+(1-p)*qA*qA)
      +bb*Dv*(N*N+2*qA*N+p*qA*qA)
      +b3*Dv*(M*M+2*qB*M+h*qB*qB)
      +b5*(M*M-2*qB*M+(1-h)*qB*qB)
      +4/a*Du*Dv*LA**2+4*L*Du*Dv*(2*u*W*Y-r*u*u*W*W)
      +4*alpha*Du*(v*Y-1)**2+2*d*Du*((v*Y-1)**2-Dv))
    V0=1-p*u*v*w
    qB0=Du*v*V0
    W0=w+u*h*v*V0
    Y0=u*w-h*v*V0
    U0=u*W0
    LA0=Y0-r*U0
    A0=w*v;B0=delta*v;C0=delta*E*u*w
    Yraw=C0/(delta*E)-h*B0*V0/delta
    Uraw=C0/(delta*E)+u*u*h*B0*V0/delta
    Sraw=(Du**2*(b2*(1-p)+bb*Dv*p)*A0**2
      +Du**2*V0**2/delta**2*(b3*Dv*h+b5*(1-h))*B0**2
      +4/a*Du*Dv*(Yraw-r*Uraw)**2
      +4*L*Du*Dv*(2*Uraw*Yraw-r*Uraw**2)
      +(4*alpha+2*d)*Du*((B0*Yraw/delta)**2-2*B0*Yraw/delta)
      -2*d*Du*B0**2/delta**2+4*alpha*Du)
    LNraw=(2*qA*(bb*Dv-b2)-2*u*v*qB0*(b3*Dv*h+b5*(1-h))
      +8/a*h*u*v*(1+r*u*u)*Dv*LA0+8*L*h*u*v*Dv*(U0-u*u*LA0)
      +2*(4*alpha+2*d)*h*u*v*v*(v*Y0-1))
    LMraw=(2*qB0*(b3*Dv-b5)-8/a*Dv*(1+r*u*u)*LA0
      +8*L*Dv*(u*u*LA0-U0)-2*(4*alpha+2*d)*v*(v*Y0-1))
    for name,actual,claimed in [
      ('complete scalar',F.subs({N:0,M:0}),Sraw),
      ('complete first normal N',S.diff(F,N).subs({N:0,M:0}),LNraw),
      ('complete first normal M',S.diff(F,M).subs({N:0,M:0}),LMraw),
    ]: zero(actual-claimed,'raw CA5 '+name)
    print('PASS all-order raw scalar and BOTH full normal formulas')


class WSeries:
    """Exact truncated series in an auxiliary WEIGHT parameter (no asymptotic rate assumptions)."""
    n=7
    def __init__(self,x):
        if isinstance(x,WSeries):self.a=x.a;return
        if isinstance(x,list):self.a=(x+[0]*self.n)[:self.n];return
        self.a=[x]+[0]*(self.n-1)
    def __add__(self,y):
        y=WSeries(y);return WSeries([self.a[i]+y.a[i] for i in range(self.n)])
    __radd__=__add__
    def __neg__(self):return WSeries([-x for x in self.a])
    def __sub__(self,y):return self+(-WSeries(y))
    def __rsub__(self,y):return WSeries(y)-self
    def __mul__(self,y):
        y=WSeries(y)
        return WSeries([sum(self.a[j]*y.a[i-j] for j in range(i+1)) for i in range(self.n)])
    __rmul__=__mul__
    def inv(self):
        a=self.a;assert a[0]!=0
        b=[S.S.One/a[0]]
        for i in range(1,self.n):b.append(-sum(a[j]*b[i-j] for j in range(1,i+1))/a[0])
        return WSeries(b)
    def __truediv__(self,y):return self*WSeries(y).inv()
    def __rtruediv__(self,y):return WSeries(y)*self.inv()
    def __pow__(self,k):
        if k<0:return self.inv()**(-k)
        z=WSeries(1)
        for _ in range(k):z=z*self
        return z


def ca6_from_full_native_source():
    k,ell,E,delta,t,u,w,v,N,M=S.symbols('k ell E delta t u w v N M')
    z=lambda val,degree:WSeries([0]*degree+[val])
    ks=z(k,1);ls=z(ell,1);us=z(u,1);ws=z(w,1);vs=z(v,2)
    K=1+ks**2;Q=K**2-ls*ks**2
    a=(K-ls)/E;eta=K*E/Q;r=ls/Q;Delta=delta*Q/K
    zeta=delta**2*ks**2*Q*E/K;R=K**2/Q;T=K-ks**2*t*zeta
    V=ks**2*eta*zeta;Qc=1-ks**2*zeta*r
    Msource=R-t;g=1-t*zeta;Enc=1-ks**2*eta*g;G=g+zeta
    P4=T*Enc+ks**2*(1-r)*G
    Pa=R*g*(1+V)+zeta*Msource*Enc
    Pd=1+2*ks**2*g+r*ks**4*G+ks**2*zeta+eta*ks**2*(K*zeta+ks**2*g*t*zeta)
    om=1-ks**2*eta+t*V
    b2=8*a*g*(1+V)/(R*Delta**2*om)
    bb=8*a*(1+(1-t)*zeta)/(R*Delta**2*om)
    b3=8*a/T
    b5=8*a*(1+V)/(T*om)
    # Exact cancellation of Q in (c3/ell), critical to correct CA6.
    L=2*a*eta**2*P4/(K*(r+eta)*om)
    d=4*a*eta*Pd/(K*(r+eta)*om)
    alpha=2*a*Delta**2*ks**6*eta**2*Pa/(K*Qc*om)
    Du=1+us**2;Dv=1+vs**2;p=t/R;h=a/T
    qA=vs*Du*ws;qB=vs*(Du-us*N-p*us*Du*vs*ws)
    W=ws+us*(M+h*qB)/Du;Y=us*ws-(M+h*qB)/Du
    LA=Y-r*us*W
    F=(b2*(N**2-2*qA*N+(1-p)*qA**2)
       +bb*Dv*(N**2+2*qA*N+p*qA**2)
       +b3*Dv*(M**2+2*qB*M+h*qB**2)
       +b5*(M**2-2*qB*M+(1-h)*qB**2)
       +4/a*Du*Dv*LA**2+4*L*Du*Dv*(2*us*W*Y-r*us**2*W**2)
       +4*alpha*Du*(vs*Y-1)**2+2*d*Du*((vs*Y-1)**2-Dv))
    SC=WSeries([S.expand(q.subs({N:0,M:0})) for q in F.a])
    LN=WSeries([S.expand(S.diff(q,N).subs({N:0,M:0})) for q in F.a])
    LM=WSeries([S.expand(S.diff(q,M).subs({N:0,M:0})) for q in F.a])
    A0=ws*vs;B0=delta*vs;C0=delta*E*us*ws;D0=delta**2*ks**3*E
    x=B0-C0;y=A0-ks*C0;zv=us*B0-ks*C0;h0=ks*C0-D0
    prod=k**2*delta**2*E  # placeholder not used as preassigned leading source
    q1sq=A0*(delta*us*vs)+(ks*C0)*D0-2*(ks**2)*B0*C0
    T4=y**2+zv**2+h0**2-2*y*zv+2*q1sq
    P=(E**2*(E+4)-E*(E**2+6*E+4)*ls+2*(E+1)*(E+2)*ls**2
      +E**2*(8-3*E)*ks**2+4*E**2*ks*us+E**2*(E+4)*us**2)
    S6=4*P/(E**4*delta**2)*x**2-16*(ks+us)/(E**2*delta**2)*x*zv+8/(E*delta**2)*T4
    LN6=WSeries(0)
    LM6=8*(E**2+2*E-2*ls*(E+1))/(E**2*delta)*x
    for label,full,pred,degree in [('scalar',SC,S6,7),('N',LN,LN6,4),('M',LM,LM6,4)]:
        for i in range(degree):zero(full.a[i]-pred.a[i],f'CA6 {label} weight {i}')
    zero(S.diff(F.a[0],N,2).subs({N:0,M:0})/2-16/(E*delta**2),'hard NN residue')
    zero(S.diff(F.a[0],M,2).subs({N:0,M:0})/2-4*(E*E+4)/E,'hard MM residue')
    zero(S.diff(S.diff(F.a[0],N),M).subs({N:0,M:0}),'hard NM residue')
    # Verify α arrives from native formula at exactly grade six, not inserted.
    for i in range(6):zero(alpha.a[i],f'alpha grade {i}')
    zero(alpha.a[6]-2*E*delta**2*k**6,'alpha native grade six')
    print('PASS CA6: full scalar 7 grades, both normals 4 grades, native alpha, full hard residue')


def exact_conductors():
    k,delta,E,u,w,v=S.symbols('k delta E u w v')
    A=w*v;B=delta*v;C=delta*E*u*w;D=delta**2*k**3*E
    x=B-C;y=A-k*C;z=u*B-k*C;h=k*C-D
    q1sq=A*(u*B)+(k*C)*D-2*k**2*B*C
    q2sq=A*(k*C)+(u*B)*D-2*k**2*B*C
    T4=y**2+z**2+h**2-2*y*z+2*q1sq
    checks=[
      C*(u-k)-(z-u*x),
      C*(w-delta*k)-(delta*y-w*x),
      C*(v-delta*k*k*E)-E*(w*z+delta*k*y-k*w*x),
      T4-(A*A+(u*B)**2+(k*C)**2+D**2-4*k*k*B*C),
      T4-(y*y+z*z+h*h-y*z+z*h+q1sq+q2sq),
    ]
    star=y+z-h-2*k*x
    checks.extend([
      star-((v-delta*k*k*E)*((w-delta*k)+delta*(u-k))
          -delta*k*E*(u-k)*(w-delta*k)),
      k*C*star-(q1sq-y*z),
    ])
    Tu,Tw,Tv=S.symbols('Tu Tw Tv')
    Jpoly=(u-k)*Tu+(w-delta*k)*Tw+(v-delta*k*k*E)*Tv
    conductor=(-u*Tu-w*Tw-E*k*w*Tv)*x+(delta*Tw+delta*k*E*Tv)*y+(Tu+E*w*Tv)*z
    checks.append(C*Jpoly-conductor)
    for i,p in enumerate(checks):zero(p,f'CC conductor/repair {i+1}')
    print('PASS CC1–CC6 exact conductor, star, two-pairing and monic-module identities')


if __name__=='__main__':
    original_native_balances()
    generic_original_packet_subtraction()
    full_raw_scalar_and_normals()
    ca6_from_full_native_source()
    exact_conductors()
    print('ALL EXACT IDENTITIES VERIFIED. Finite graded module lift proved in accompanying note.')
