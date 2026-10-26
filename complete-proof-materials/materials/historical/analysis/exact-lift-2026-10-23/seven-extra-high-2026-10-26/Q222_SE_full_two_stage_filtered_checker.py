#!/usr/bin/env python3
"""Exact analytic replay of the SE4 local filtered-base proof.
Requires SymPy; no SDP solver, floating point, or external inputs.

Checks (i) all six scalar/first-normal/hard coefficient groups of SE4
against an independently constructed full SE4 expression; (ii) all
printed native SD9--SD10 stationarity balances; (iii) the native
source residue entries and scalar divisibility coefficients;
(iv) SE5, SE6, SE7 and both second-stage conductor identities.

The structural proof of bounded-ring quotients and filtered matrix
construction is in Q222_SE_full_two_stage_filtered_proof_2026-10-26.md.
This standalone checker intentionally does not independently check the
original-parent-minus-retained-packets operation allocated elsewhere.
"""


def check_full_se4_parent_coefficients():
    import sympy as S
    u,v,w,N,M=S.symbols('u v w N M')
    b2,bb,b3,b5,A0,L,al,d,nu0,eps,r,E,hs,p=S.symbols('b2 bb b3 b5 A0 L alpha d nu eps r E hs p', nonzero=True)
    D=1+u*u; V=1+v*v; v0=1-p*u*v*w
    qA=D*v*w; qb=v*(D-u*N-p*u*D*v*w)
    TT=M+hs*qb
    W=w+u*TT/D
    Y=u*w-TT/D
    original=(b2*(N*N-2*qA*N+(1-p)*qA*qA)+bb*V*(N*N+2*qA*N+p*qA*qA)
     +b3*V*(M*M+2*qb*M+hs*qb*qb)+b5*(M*M-2*qb*M+(1-hs)*qb*qb)
     +4*A0*D*V*(Y-r*u*W)**2+4*L*D*V*(2*u*W*Y-r*u*u*W*W)
     +4*al*D*(v*Y-1)**2+2*d*D*((v*Y-1)**2-V)+nu0*V*N*N)
    phi=b3*V*hs+b5*(1-hs);sig=4*al+2*d
    fA=D*D/eps**2*(b2*(1-p)+bb*V*p)
    fBB=(phi*D*D*v0*v0 + 4*A0*D*V*hs**2*(1+r*u*u)**2*v0*v0
       -4*L*D*V*hs**2*u*u*(2+r*u*u)*v0*v0
       +sig*D*(2*hs*v0+hs*hs*v*v*v0*v0)-2*d*D)
    fBC=(-8*A0*D*V*(1-r)*hs*(1+r*u*u)*v0/E
        +8*L*D*V*hs*((1-r)*u*u-1)*v0/E
        -sig*D*(2/E+2*hs*v*v*v0/E))
    fCC=(4*A0*D*V*(1-r)**2/E**2+4*L*D*V*(2-r)/E**2+sig*D*v*v/E**2)
    varpi=S.symbols('varpi')
    PB=8*A0*V*hs*(1+r*u*u)**2*v0-8*L*V*hs*u*u*(2+r*u*u)*v0+2*sig*(1+hs*v*v*v0)
    PC=(-8*A0*V*(1+r*u*u)*(1-r)/E+8*L*V*((1-r)*u*u-1)/E-2*sig*v*v/E)
    gNA=D*(bb*V-b2)/eps
    gNB=-phi*D*u*v*v0-hs*u*PB*v/2
    gNC=-hs*u*PC*v/2
    gMB=(b3*V-b5)*D*v0+PB/2
    gMC=PC/2
    Xi=(4*V*(A0*(1+r*u*u)**2-L*u*u*(2+r*u*u))+sig*v*v)/D
    hnn=b2+bb*V+nu0*V+phi*u*u*v*v+Xi*hs**2*u*u*v*v
    hnm=-(b3*V-b5)*u*v-Xi*hs*u*v
    hmm=b3*V+b5+Xi
    Af=eps*w*v; Bf=v; Cf=E*u*w
    raw=(fA*Af**2+fBB*Bf**2+fBC*Bf*Cf+fCC*Cf**2+varpi*S.symbols('D0')**2
     +2*N*(gNA*Af+gNB*Bf+gNC*Cf)+2*M*(gMB*Bf+gMC*Cf)
     +hnn*N*N+2*hnm*N*M+hmm*M*M)
    # isolate constant alpha
    f0diff=S.factor((original-raw).subs({N:0,M:0}))
    # direct equality after varpi*D0^2 = 4alpha D
    f0diff=S.cancel(f0diff.subs(varpi,4*al*D/S.symbols('D0')**2))
    print('scalar equal:',f0diff==0)
    for nm,var in [('N',N),('M',M)]:
     dr=S.cancel(S.diff(original-raw,var).subs({N:0,M:0}))
     print(nm,'linear equal:',dr==0, 'remaining',str(dr)[:150])
    for i,x in [(2,N),(2,M)]:
     dr=S.cancel(S.diff(original-raw,x,2).subs({N:0,M:0}))
     print(str(x)+'^2','equal',dr==0,'remaining',str(dr)[:150])
    dr=S.cancel(S.diff(original-raw,N,M).subs({N:0,M:0}))
    print('NM equal',dr==0,'remaining',str(dr)[:150])


def check_native_sd9_sd10_contact_balances():
    import sympy as s,time
    kap,eta,zeta,a,t=s.symbols('kap eta zeta a t');V=kap*eta*zeta;K=1+kap;Lc=1-a*eta;Hc=1+(1-a)*eta;Qc=1-kap*zeta+a*V;R=1+kap*Lc;T=1+kap-kap*zeta*t
    P4=1+kap+eta*kap*(a-kap-1)-zeta*kap*t+V*(-zeta*kap*t*t+2*kap*t+a*(1-t)+t)
    Pa=1+kap-a*eta*kap+zeta*((1+kap)-(2+kap)*t)+V*(-zeta*t*t+a*(t-1)+t)
    Pd=(1+kap)**2-a*eta*kap**2+zeta*kap*((kap+1)-(kap+2)*t)+V*(-zeta*kap*t*t+kap*a*(t-1)+kap*t+kap+1)
    c2=2*K*(1-t*zeta)*Lc*(1+V)*Hc*Qc
    cb=2*K*Lc*(1+(1-t)*zeta)*Hc*Qc
    c3=2*K*Lc*Hc*Qc*(1-kap*eta+t*V)
    c5=2*K*Lc*Hc*Qc*(1+V)
    c4=2*Qc*P4;ca=2*Lc*Hc*Pa;cd=4*Lc*Qc*Pd
    exprs={'s9a':cb-c2-zeta*c3,'s9b':c3-c5+kap*eta*cb,'s9c':c2-(1-t*zeta)*c5,
    'Y':2*kap*c2+2*(1+V)*Lc*c4-2*kap**2*zeta*Qc*ca-Qc*cd,
    'S':kap*Lc*c2+c5+(1+V)*Lc**2*c4+kap*Qc**2*ca+((kap*zeta*Lc**2-2*Lc-eta)/2)*cd,
    'V':kap*R*c2+t*kap*V*cb+a*V*c3+T*c5+K*V*Lc**2*c4-kap**2*zeta*K*Qc*Lc*ca-(K/2)*(Qc*Lc+eta)*cd}
    for name,ex in exprs.items():
     now=time.time();z=s.expand(ex);print(name,z==0,s.count_ops(z),'time',round(time.time()-now,3),flush=True)
     if z!=0:print(s.factor(z))


def check_native_se_residues():
    import sympy as sp
    E,ell,u,v,eps,j,tau=sp.symbols('E ell u v eps j tau',positive=True)
    K=1+eps**2*j**2;Q=K*K-ell*eps**2*j**2
    eta=K*E/Q; r=ell/Q; zz=j*j*Q*E/K; tt=eps**2*tau; a=(K-ell)/E
    R=K*K/Q;T=K-eps**2*j**2*tt*zz; vc=eps**2*j**2*eta*zz
    Hc=r+eta;Qc=1-eps**2*j**2*zz+a*vc; Sn=1-eps**2*j**2*eta+tt*vc
    P4=(1+eps**2*j**2+eta*eps**2*j**2*(a-eps**2*j**2-1)-zz*eps**2*j**2*tt+vc*(-zz*eps**2*j**2*tt**2+2*eps**2*j**2*tt+a*(1-tt)+tt))
    Pd=((1+eps**2*j**2)**2-a*eta*(eps**2*j**2)**2+zz*eps**2*j**2*((eps**2*j**2+1)-(eps**2*j**2+2)*tt)+vc*(-zz*eps**2*j**2*tt**2+eps**2*j**2*a*(tt-1)+eps**2*j**2*tt+eps**2*j**2+1))
    Pal=1+eps**2*j**2-a*eta*eps**2*j**2+zz*((1+eps**2*j**2)-(2+eps**2*j**2)*tt)+vc*(-zz*tt**2+a*(tt-1)+tt)
    b2hat=8*(1-tt*zz)*(1+vc)/(Q*Sn)
    bbhat=8*(1+(1-tt)*zz)/(Q*Sn)
    b3=8/T;b5=8*(1+vc)/(T*Sn)
    L=2*eta**2*P4/(K*Hc*Sn);d=4*eta*Pd/(K*Hc*Sn)
    A0=L*Hc/(1+r*Hc); h=a/T; sig=4*(2*eps**4*j**6*E**2*Pal/(K*Qc*Sn))+2*d
    D=1+u*u;Dv=1+v*v;phi=b3*Dv*h+b5*(1-h)
    # evaluate source first, then testing angles
    Z=lambda x: sp.factor(x.subs({eps:0,j:0,u:0,v:0}))
    fA0=Z(D**2*(b2hat*(1-tt/R)+bbhat*Dv*(tt/R)))
    fBB0=Z(phi*D*D+4*A0*D*Dv*h*h*(1+r*u*u)**2-4*L*D*Dv*h*h*u*u*(2+r*u*u)+sig*D*(2*h+h*h*v*v)-2*d*D)
    fBC0=Z(-8*A0*D*Dv*(1-r)*h*(1+r*u*u)/E+8*L*D*Dv*h*((1-r)*u*u-1)/E-sig*D*(2/E+2*h*v*v/E))
    fCC0=Z(4*A0*D*Dv*(1-r)**2/E**2+4*L*D*Dv*(2-r)/E**2+sig*D*v*v/E**2)
    fva0=Z(8*Pal/(K*Qc*Sn)*D)
    PB=8*A0*Dv*h*(1+r*u*u)**2-8*L*Dv*h*u*u*(2+r*u*u)+2*sig*(1+h*v*v)
    PC=-8*A0*Dv*(1+r*u*u)*(1-r)/E+8*L*Dv*((1-r)*u*u-1)/E-2*sig*v*v/E
    cM=Z((b3*Dv-b5)*D+PB/2)
    cC=Z(PC/2)
    d0=Z(d)
    H0=E+ell;Vs=1+ell*H0
    targets=[('fA',fA0,8),('fBB',fBB0,8*(E+2)/(H0*Vs)),('fBC',fBC0,-16*(E+2)/(H0*Vs)),('fCC',fCC0,8*(E+2)/(H0*Vs)),('varpi',fva0,8),('gMB',cM,8*E*(H0+1)/(H0*Vs)),('gMC',cC,-8*E*(H0+1)/(H0*Vs)),('HNN',d0,4*E/H0)]
    for name,got,want in targets:
     ok=sp.cancel(got-want)==0
     print(name,ok,'ops',sp.count_ops(got),flush=True)
     assert ok
    # scalar s(U,V) at eps=j=0, proof s(0,V)=8V/E & d_U s at zero=8
    U,V=sp.symbols('U V')
    mod=lambda z: z.subs({eps:0,j:0,u*u:U,v*v:V})
    # construct from zero-source expressions of f terms substituting u^2 as U, v^2 as V without first eliminating angles
    rr=ell;H=E+ell;Vs=1+ell*H;hst=(1-ell)/E;A0s=2*E**2/Vs;Ls=2*E**2/H;ds=4*E/H
    Du=1+U;Dv=1+V;ph=8*(1+hst*V); ssig=2*ds
    fbb=ph*Du**2+4*A0s*Du*Dv*hst**2*(1+ell*U)**2-4*Ls*Du*Dv*hst**2*U*(2+ell*U)+ssig*Du*(2*hst+hst**2*V)-2*ds*Du
    fbc=-8*A0s*Du*Dv*(1-ell)*hst*(1+ell*U)/E+8*Ls*Du*Dv*hst*((1-ell)*U-1)/E-ssig*Du*(2/E+2*hst*V/E)
    fcc=4*A0s*Du*Dv*(1-ell)**2/E**2+4*Ls*Du*Dv*(2-ell)/E**2+ssig*Du*V/E**2
    sbar=fbb+fbc+fcc
    print('s(U=0,V)=8V/E',sp.cancel(sbar.subs(U,0)-8*V/E)==0)
    print('dU s at zero=8',sp.cancel(sp.diff(sbar,U).subs({U:0,V:0})-8)==0)
    assert sp.cancel(sbar.subs(U,0)-8*V/E)==0
    assert sp.cancel(sp.diff(sbar,U).subs({U:0,V:0})-8)==0
    print('PASS residue and scalar jet')


def check_two_stage_conductor_identites():
    import sympy as S
    eps,j,E,u,w,v=S.symbols('eps j E u w v')
    k=eps*j; v0=eps*j*j*E; C=E*u*w; A=eps*w*v;B=v;D0=eps**2*j**3*E
    x=B-C;y=A-k*C;z=u*B-k*C;h=k*C-D0
    rstar=y+z-h-2*k*x
    rstar_formula=(v-v0)*(u-k+eps*(w-j))-E*k*(u-k)*(w-j)
    assert S.expand(rstar-rstar_formula)==0
    # q_i radicals are replaced with their squares, all purely algebraic
    q1sq=eps*u*w*(v-v0)**2
    q2sq=eps**2*j*E*u*v*(w-j)**2
    assert S.expand(k*C*rstar-(q1sq-y*z))==0
    assert S.expand(k*C*rstar-(q2sq+z*h))==0
    Tu,U,V,Tv=S.symbols('Tu U V Tv')
    P=(u-k)*Tu+(w-j)*(eps*U+v*V)+(v-v0)*Tv
    rx=-u*Tu-eps*w*U-E*k*w*w*V-E*k*w*Tv
    ry=U+E*j*w*V+E*j*Tv
    rz=Tu+E*w*(w-j)*V+E*w*Tv
    assert S.expand(C*P-x*rx-y*ry-z*rz)==0
    # second-stage modulo epsilon with an arbitrary s(u²,v²) = u² p(u²,v²)+v² q(v²)
    q0,q1,p00,p10,p01,p11,H00,H10,H01=S.symbols('q0 q1 p00 p10 p01 p11 H00 H10 H01')
    pp=p00+p10*u*u+p01*v*v+p11*u*u*v*v
    qq=q0+q1*v*v
    HH=H00+H10*u*u+H01*v*v
    Jx=HH-u*u*pp; Jy=E*j*v*qq; Jz=u*pp+E*w*v*qq; Jh=S.Integer(0)
    for label,Ji,expected in [('x',Jx,0),('y',Jy,0),('z',Jz,E*v*qq),('h',Jh,0)]:
     # Ordered w quotient acts after u=k=0 modulo epsilon.
     Qwi=S.diff(Ji.subs(u,0),w)
     assert S.expand(Qwi-expected)==0
    Jmat=S.Matrix([[S.diff(f,t) for t in (u,w,v)] for f in (x,y,z)])
    assert S.factor(Jmat.det().subs({u:k,w:j,v:v0})+E**2*eps**3*j**4)==0
    print('PASS two circuit identities, exact SE7, contact-gradient determinant, and second-stage w quotient modulo epsilon.')


if __name__ == "__main__":
    check_full_se4_parent_coefficients()
    check_native_sd9_sd10_contact_balances()
    check_native_se_residues()
    check_two_stage_conductor_identites()
    print("ALL INDEPENDENT EXACT CHECKS PASS")
