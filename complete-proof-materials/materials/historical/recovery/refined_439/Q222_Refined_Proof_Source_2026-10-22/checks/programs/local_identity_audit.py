"""Independent exact algebra audit of the recovered R1 proof (2026-10-22).

The formulas below are transcribed from the TeX proof, not imported from its
checkers. SymPy is used only as a diagnostic expansion engine.
"""
import sys
import sympy as s

j,B,T,p,H=s.symbols('j B T p H', positive=True)
u,v,z,N,M=s.symbols('u v z N M')
U,V,Z=s.symbols('U V Z')
half=s.Rational(1,2)
x,y,qz=j*B,T*B,j*T
K=j+1-p
Pc=j+1+B*((j+1)**2-p)+T*B*(H*j+j+1+B*(H*j*(2*j+2-p)+p*(2*j+1-p)))+H*j*T**2*B**2*(1+B*p+B*H*j)
Wc=(1-H)*j+1-p
E=1+B*(H*T+p-1)
Lc=K+T*(Wc-H)+B*T*(p*K+H*j*(1-p))+B*H*T**2*(Wc-1)
Nc=p+B*(-j**2+3*j*p-j-p**2+p)+B*T*(H*j**2+H*j-j**2+j*p-j+p)+B**2*T*(2*H*j**2*p-H*j**2-H*j*p**2+2*H*j*p-2*j**2*p+3*j*p**2-j*p-p**3+p**2)+B**2*H*T**2*j*(H*j-j+p)*(1+B*p)
mu=j*p*(1+y)*(1+x)*(1+H*y)/(2*Pc)
m=mu/(1-p)
n=Nc/(2*(j-p)*Pc)
C=Lc/(2*H*T*Pc)
A=T*(1+x)*(1+H*y)*(K+H*j*y)/(2*(1-H*T)*Pc)
D=(1+H*y)*(K+H*j*y)*(1-y*((1-H)*j-p))/(2*T*(1-H)*Pc)
alpha=B**2*(K+H*j*y)*Lc/(2*E*Pc)

audit_passes=[]

def check(name,expr):
    num=s.together(expr).as_numer_denom()[0]
    out=s.expand(num)
    if out != 0:
        print('FAIL',name, s.factor(out),flush=True)
        raise AssertionError(name)
    audit_passes.append(name)
    print('PASS',name,flush=True)

def coefficient_checks(name,expr,vars):
    for mon,coeff in s.Poly(s.expand(expr),*vars).terms():
        check(name+' coefficient '+str(mon),coeff)

if '--native' in sys.argv or '--all' in sys.argv:
    a=(j+p*x)/(1+H*T*x)
    b=B*(1+H*T*x)/(1+p*B)
    t=T*(1+p*B)/(1+H*T*x)
    r=(j-p)/(j+p*x)
    h=(1-H*T)/(1+H*T*x)
    ca=s.cancel(1+a*b*r)
    sc=s.cancel(1+b*(1-h+a*r))
    gm=s.cancel(1+a*(h*b+r-1))
    gd=s.cancel(ca-(1-h)/t)
    de=s.cancel(a+(1+a*b)+(1+a*b)*b+2*a*b*(1+a*b)*r+a**3*b**2*r**2+2*a**2*b*r+t*b*(1+a+a*b+a**2*b*r*(2-r)))
    mass=s.cancel((1+a*b)*(1+a*r)*(h+t)-a*(1+t*b)*(1-r)*(1-h))
    native={
        'n':((1+a*r)*ca*sc-r*de)/(2*r*de),
        'm':a*a*(1+t*b)*(1-r)*ca*sc/(2*gm*de),
        'C':ca*mass/(2*(1-h)*de),
        'A':t*ca*(1+a*r)*sc/(2*h*de),
        'D':(1+a*b)*(1+a*r)*(1-a*t*b*r)*sc/(2*t*gd*de),
        'alpha':b*b*(1+a*r)*mass/(2*(1-h*b)*de)}
    for name,val in native.items():
        check('native source '+name, val-globals()[name])

if '--normal' in sys.argv or '--all' in sys.argv:
    ma,na,ca,aa,da,al=s.symbols('mu n C A D alpha')
    Du,Dv=1+u*u,1+v*v
    ga=ma*Du*(N*N-2*v*N+(1-p)*v*v)+na*Dv*(N+p*v-u*z)**2+ca*Du*Dv*(M+H*v*z)**2+aa*Dv*(u-M-H*v*z)**2+(da+half)*Du*(M-(1-H)*v*z)**2+(al+half)*(u*(1-p)*v-(1+H*v*v)*z-u*N-v*M)**2+half*Du*Dv*(N+p*v)**2-half*Dv*(z+u*N+p*u*v)**2
    W=1+u*(M+H*v*z); XX=N+p*v-u*z
    YY=u-M-H*v*z; ZZ=z+u*N+p*u*v
    parent=ma/(1-p)*(v*W-XX+u*v*YY-u*ZZ)**2+na*Du*Dv*XX**2+ca*Dv*(u*W-YY)**2+aa*Du*Dv*YY**2+(da+half)*(u*W+u*v*XX-YY-v*ZZ)**2+(al+half)*Du*(v*YY-ZZ)**2+half*Dv*(XX+u*ZZ)**2-half*Du*Dv*ZZ**2
    check('entire native parent to full normal G',(parent-ma*p*Du**2*N**2/(1-p))/Du-ga)
    kk=ca*H**2+(da+half)*(1-H)**2
    ftable=aa*u**2+(ma*(1-p)+(na+half)*p**2)*v**2+al*z**2+(aa+al*(1-p)**2+ma*(1-p)+half*(1-p)**2)*u**2*v**2+na*u**2*z**2+(kk+aa*H**2+2*H*al+H-half)*v**2*z**2+(na+half)*p**2*v**4+(kk+na)*u**2*v**2*z**2+H**2*(aa+ca+al+half)*v**4*z**2+ca*H**2*u**2*v**4*z**2-2*(aa*H+al*(1-p)+na*p+half)*u*v*z-2*(aa*H+H*al*(1-p)+na*p+half*(H*(1-p)+p))*u*v**3*z
    check('complete scalar coefficient table',ga.subs({N:0,M:0})-ftable)
    hn=(ma+half)*Du+half*v*v+al*u*u+na*Dv
    hm=Du*(ca*Dv+da+half)+aa*Dv+(al+half)*v*v
    check('full HN',s.diff(ga,N,2)/2-hn)
    check('full HM',s.diff(ga,M,2)/2-hm)
    check('full NM cross',s.diff(ga,N,M)/2-(al+half)*u*v)
    b=v-u*z/j
    q=v*z-T*u
    R=z-u*v/B
    Ru=v*(u*u-x)
    kc=p/2+j*(H*(al+half)-half)+na*(p-j)
    kr=-ma-(al+half)*(1-p)+T*kc/j
    sn=na*p-ma+p/2+(-ma-(al+half)*(1-p))*x+j*(al-na)+y*kc
    ln=(j*(na-al)+p*(na+half)*v*v)*b+kr*Ru+kc*u*v*q/j
    # The source text specifies derivative minus displayed formula.
    check('full first normal symbolic cancellation',s.diff(ga,N).subs({N:0,M:0})/2-ln-v*sn)
    kd=ca*H-(da+half)*(1-H)
    cq=kd+aa*H+half+v*v*(ca*H*Du+aa*H+(al+half)*H)
    lm=cq*q+ma*u*v*b+v*(al-B*T*ca*H*u*u)*R
    s0=T*(kd+aa*H+half)-aa
    s1=j*kd+ma+x*T*ca*H
    s2=j*kd+x*T*ca*H+T*(ca*H+aa*H+(al+half)*H)-aa-(al+half)*(1-p)+al/B
    check('full second normal symbolic cancellation',s.diff(ga,M).subs({N:0,M:0})/2-lm-u*s0-u*u*v*z*s1/j-u*v*v*(s2-s1))
    sub={ma:mu,na:n,ca:C,aa:A,da:D,al:alpha}
    for name,expr in [('SN',sn),('S0',s0),('S1',s1),('S2',s2)]:
        check('actual source '+name,expr.subs(sub))

if '--scalar' in sys.argv or '--all' in sys.argv:
    K0=C*H**2+(D+half)*(1-H)**2
    Ma=A*x
    Mb=(mu*(1-p)+(n+half)*p*p)*y
    Mc=alpha*qz
    Mab=(A+alpha*(1-p)**2+mu*(1-p)+half*(1-p)**2)*x*y
    Mac=n*x*qz
    Mbc=(K0+A*H*H+2*H*alpha+H-half)*y*qz
    Mb2=(n+half)*p*p*y*y
    Mabc=(K0+n)*x*y*qz
    Mb2c=H*H*(A+C+alpha+half)*y*y*qz
    Mab2c=C*H*H*x*y*y*qz
    O1=(A*H+alpha*(1-p)+n*p+half)*j*B*T
    O3=(A*H+H*alpha*(1-p)+n*p+half*(H*(1-p)+p))*j*B*B*T*T
    kc=(Ma-Mb+Mac+Mab2c-Mbc)/2
    rc=Mc-Mb+Mac-kc
    fc=Mb-Mac
    qc=Ma-fc
    ec=Mb2-Mab2c
    dc=Mb2c-kc
    gc=Mab2c
    Q=U-V*Z; BB=V-U*Z; R=Z-U*V
    psi=U*U+V*V+Z*Z+U*U*V*V*Z*Z-4*U*V*Z
    # Compare only the non-tautological linear coefficient identities.
    check('scalar UV contact',qc+Mac+rc+2*fc-O1)
    check('scalar UV^3Z contact',dc+ec+2*gc-O3)
    check('scalar U^2V^2 balance',rc+dc+gc-Mab)
    check('scalar U^2V^2Z^2 balance',ec+fc-Mabc)
    f1=-2*B*H*j*p+B*H*j+2*B*j*p-2*B*p*p+2*B*p-H*j-H*p-H+j+1-p*p/j+p/j
    f2=-B*B*H*H*j*p+B*B*H*j*p-B*B*H*p*p+B*B*H*p-B*H*H*j-B*H*H*p+B*H*j-B*H*p*p/j+B*H*p/j
    fpoly=K+T*f1+T*T*f2
    kpoly=p*Wc+B*(1-p)*(j*(H*(p+1)-p)-(1-p)*(p-H))+T*B*H*((p+1)*(1-H)*j+(1-p)*(B*H*j+p+1-H))+T*T*B*B*H*H*j*(1-H)
    check('closed fc',fc-T*x*x*fpoly/(2*Pc))
    check('closed rc',rc-B*B*j*T*(1-p)*(1+x)*(1+H*y)*Lc/(2*E*Pc))
    check('closed gc',gc-H*T*T*j*j*B**3*Lc/(2*Pc))
    check('closed kc',kc+T*T*B*B*j*(1+x)*kpoly/(2*E*Pc))

if '--outer' in sys.argv or '--all' in sys.argv:
    W,X,Y,Zs=s.symbols('W X Y Zs')
    ca,da,aa,g,hc=s.symbols('C D A g h')
    Du,Dv=1+u*u,1+v*v
    R=u*W-Y; Ct=v*(u*X-Zs); Bt=X+u*Zs
    core=((R+Ct)**2+Du*(v*Y-Zs)**2+Dv*Bt**2-Du*Dv*Zs**2)/2
    # Eliminate radicals by expanding the proposed row-Gram evaluations.
    qeval=4*g*Ct**2-2*Dv*Ct**2/(Du*v*v)+4*u*Dv*Ct*Bt/(Du*v)+2*Dv*Bt**2/Du
    rhs=4*(ca+da)*(R+hc*Ct)**2+4*ca*v*v*R**2+2*(R+Ct)**2+2*Du*(v*Y-Zs)**2+qeval
    check('outer R1N13',(4*ca*Dv*R**2+4*da*(R+Ct)**2+4*core-rhs).subs({g:ca*da/(ca+da),hc:da/(ca+da)}))
    # R1N14 expanded using K1, Y1 without square-root symbols.
    B1=Bt-u*v*Y; V0=v*Y
    ge=4*g*Ct**2+4*Ct*Y+Y*Y/g+2*Dv*(B1*B1+V0*V0)+(4*aa*g*Dv-1/Du)*Du*Y*Y/g+4*v*Ct*V0-2*Du*v*v*(v*Y-Zs)**2
    rhs=4*(ca+da)*(R+hc*Ct)**2+4*ca*v*v*R**2+2*(R+Ct)**2+ge
    check('outer R1N14',(4*ca*Dv*R**2+4*da*(R+Ct)**2+4*aa*Du*Dv*Y*Y+4*core-rhs).subs({g:ca*da/(ca+da),hc:da/(ca+da)}))
    # Remaining Q3 source-residue determinants.
    mm,nn,bb=s.symbols('mcap ncap bshift')
    vrow=s.Symbol('Vrow')
    q3=mm*(Du*vrow/u-X-u*Zs)**2+nn*Du*X*X+Du*vrow*vrow/2-Du*vrow*Zs+(X+u*Zs)**2/2-2*bb*vrow*(u*X-Zs)
    mat=s.hessian(q3,(X,Zs,vrow))/2
    m0=j*p/(2*(1-p)*(j+1)); n0=p/(2*(j-p)*(j+1)); k0=K/(j+1)
    mi=mat.subs({mm:m0,nn:n0,bb:k0/2})
    check('R1N22 second minor',mi[:2,:2].det()-p*u*u*Du*K/(4*(j-p)*(j+1)**2*(1-p)))
    check('R1N22 determinant',mi.det()-p*u*u*Du*K*K/(8*(j-p)*(j+1)**3))
    kn=mat[1:,1:].subs({mm:j*j/(2*(1-j)*(j+1)),bb:1/(2*(j+1))})
    check('n-opening kernel determinant',kn.det()-u*u/(4*(j+1)**2))
    km=s.hessian(q3.subs(Zs,Du*vrow/(u*u)-X/u),(X,vrow))/2
    km=km.subs(nn,1/(2*(j*j-1)))
    check('m-opening first diagonal',km[0,0]-Du/(2*(j*j-1)))
    check('R1N23 kernel determinant',km.det()+Du**2*(2*bb*(j-1)-j)*(2*bb*(j+1)-j)/(4*u*u*(j*j-1)))
    check('R1N24 common kernel',q3.subs({X:0,Zs:Du*vrow/(u*u)})-Du*(4*bb-1)*vrow*vrow/(2*u*u))

if '--circuits' in sys.argv or '--all' in sys.argv:
    Q=U-V*Z; BB=V-U*Z; R=Z-U*V
    psi=U*U+V*V+Z*Z+U*U*V*V*Z*Z-4*U*V*Z
    check('explicit Psi SOS',psi-((U-V)**2+Z**2*(1-U*V)**2+2*U*V*(1-Z)**2))
    check('RFN18a',(1+U*U)*psi-R*R-(BB*BB+U*U*Q*Q+U*U*(1-U*V*Z)**2))
    check('RFN18b',(1+V*V)*psi-V*V*BB*BB-(Q*Q+R*R+V*V*(1-U*V*Z)**2))
    tt=U*V; dd=U-V; ee=Z*(1-tt); tstar=tt*(1-Z); js=tstar-dd*dd*Z
    check('RFN18c',2*(psi+V*V*Q*Q+U*U*BB*BB)-R*R-((ee+tstar)**2+2*dd*dd+4*tt*(1-Z)**2+2*js*js+4*tt*dd*dd*Z*Z))
    t7=s.Matrix([Z-U*U*V,V-U*U*V,U-U*U*V,U*Z-U*U*V,U*V-U*U*V,U*V*Z-U*U*V,U*U-U*U*V])
    t3=s.Matrix([1-U*V,V-U*V,U-U*V])
    G7=s.Matrix([[4,-2,-3*half,0,-half,0,half],[-2,3,0,-3*half,0,-1,2],[-3*half,0,4,3*half,-2,-3,0],[0,-3*half,3*half,4,-3*half,0,-3*half],[-half,0,-2,-3*half,3,2,-1],[0,-1,-3,0,2,8,-4],[half,2,0,-3*half,-1,-4,4]])
    G3=s.Matrix([[3,-2,-2],[-2,2,2],[-2,2,3]])
    check('RFN19 radial identity',4*psi+4*U*U*Q*Q+4*U*U*R*R-V*V*(U*U-1)**2-((t7.T*G7*t7)[0]+4*V*Z*(1-U)**2+U*Z*(t3.T*G3*t3)[0]))
    assert [G7[:i,:i].det() for i in range(1,8)]==[4,8,s.Rational(101,4),s.Rational(1049,16),s.Rational(907,16),s.Rational(3279,16),s.Rational(743,16)]
    assert [G3[:i,:i].det() for i in range(1,4)]==[3,2,2]
    print('PASS RFN19 positive principal minors',flush=True)

if '--r34' in sys.argv or '--all' in sys.argv:
    W,X,Y,Zs,VR,Nr,ss,ww=s.symbols('W X Y Zs VR Nr ss ww')
    Du,Dv=1+u*u,1+v*v
    kzero=(Du*(v*Y-Zs)**2+Dv*(X+u*Zs)**2-Du*Dv*Zs**2)/2
    r3=T*Du*Dv*Y*Y/2+(1+T)*v*v*(u*X-Zs)**2/(2*T)+kzero
    r4=Dv*v*v*(Zs-u*X)**2/(2*T)+T*Du*Dv*Y*Y/(2*(1-T))+kzero
    mats=[s.hessian(r3,(Y,X,Zs))/2,s.hessian(r4,(Y,X,Zs))/2]
    mins3=[Du*(T*Dv+v*v)/2,Du*(T*Dv+v*v)*(T*u*u*v*v+T*v*v+T+u*u*v*v)/(4*T),v*v*(1+T)*Du**2*Dv*(T*u*u+v*v)/(8*T)]
    mins4=[Du*(T+v*v)/(2*(1-T)),Du*Dv*(T+v*v)*(T+u*u*v*v)/(4*T*(1-T)),v*v*Du**2*Dv**2*(T*u*u+v*v)/(8*T*(1-T))]
    for label,mat,mins in [('R3',mats[0],mins3),('R4',mats[1],mins4)]:
        for i in range(3):check('EN22 '+label+' minor '+str(i+1),mat[:i+1,:i+1].det()-mins[i])
    origcore=((u*W+u*v*X-Y-v*Zs)**2+(1+u*u)*(v*Y-Zs)**2+(1+v*v)*(X+u*Zs)**2-(1+u*u)*(1+v*v)*Zs**2)/2
    shown=((W+(1+ss*ss)*VR-ss*Y-ss*v*Nr)**2+(1+ss*ss)*v*v*Y*Y-2*(1+ss*ss)*Y*(v*Nr-ss*VR)-(1+ss*ss)*(v*Nr-ss*VR)**2+Dv*Nr*Nr)/2
    check('EN infinite Alice complete signed core',(origcore/u**2).subs(u,1/ss)-shown.subs({Nr:ss*X+Zs,VR:v*X}))

if '--r2' in sys.argv or '--all' in sys.argv:
    apar=s.Symbol('aparam',positive=True)
    kk=apar/(apar-1)
    Q=U-V*Z; BB=V-U*Z; RR=Z-U*V
    psi=U**2+V**2+Z**2+U**2*V**2*Z**2-4*U*V*Z
    dd=U-V; tt=U*V; ee=Z*(1-tt); tstar=tt*(Z-1)
    check('SH1 parameterized opposite identity',apar*psi+kk*V**2*(Q**2+BB**2)-RR**2-(apar*dd**2+(apar-1)*(ee-tstar/(apar-1))**2+(2*apar*tt+kk*V**4)*(1-Z)**2+2*kk*Z*(tt-V**2)**2))
    check('SH4 source cancellation',alpha/y-C*H*x-B*C*H*(1-p)*(1+x)/E)
    ma,na,ca,aa,da,al=s.symbols('mu n C A D alpha')
    Du,Dv=1+u**2,1+v**2
    hn=(ma+half)*Du+half*v*v+al*u*u+na*Dv
    hm=Du*(ca*Dv+da+half)+aa*Dv+(al+half)*v*v
    kd=ca*H-(da+half)*(1-H)
    cq=kd+aa*H+half+v*v*(ca*H*Du+aa*H+(al+half)*H)
    s0=T*(kd+aa*H+half)-aa
    check('RL23 entire shifted second normal',cq-H*hm+da+H*(ca+da+half)*u*u)
    fc,kc=s.symbols('fc kc')
    dc=H**2*(aa+ca+al+half)*y**2*qz-kc
    actual=(aa*x-fc+dc*v*v/y)/(T*x)-2*T*H*cq+T*H**2*hm
    stated=(aa*(1-H*T)**2+(H*T)**2*(ca+da+half)*Du)/T-fc/(T*x)-T*H**2*ca*u*u*v*v-kc*v*v/(y*T*x)
    check('RL14 complete q diagonal modulo S0',actual-stated+2*H*s0)
    bcur=v-u*z/j; qcur=v*z-T*u; Ru=v*(u*u-x)
    kcn=p/2+j*(H*(al+half)-half)+na*(p-j)
    kr=-ma-(al+half)*(1-p)+T*kcn/j
    ln=(j*(na-al)+p*(na+half)*v*v)*bcur+kr*Ru+kcn*u*v*qcur/j
    ks=-ma+(al+half)*(H*T+p-1)
    kn=-(na+half)*(1-p/j)
    check('RL22 full shifted first normal',ln-H*(al+half)*u*v*qcur-((j*(na-al)+p*(na+half)*v*v)*bcur+ks*Ru+kn*v*(u*v*z-j*y)))

if '--r2comparison' in sys.argv or '--all' in sys.argv:
    ss,zz,qq,ww,hh=s.symbols('S z q w h',positive=True)
    ee=zz*ww; rr=1-ee
    aa=1/zz; bb=ss*zz*qq; tt=1/ss
    ca0=1+aa*bb*rr; sc0=1+bb*(1-hh+aa*rr)
    mgap=1+aa*(hh*bb+rr-1); dgap=ca0-(1-hh)/tt
    delta=s.cancel(aa+(1+aa*bb)+(1+aa*bb)*bb+2*aa*bb*(1+aa*bb)*rr+aa**3*bb**2*rr**2+2*aa**2*bb*rr+tt*bb*(1+aa+aa*bb+aa**2*bb*rr*(2-rr)))
    mass=s.cancel((1+aa*bb)*(1+aa*rr)*(hh+tt)-aa*(1+tt*bb)*(1-rr)*(1-hh))
    nact=s.cancel(((1+aa*rr)*ca0*sc0-rr*delta)/(2*rr*delta))
    mact=aa**2*(1+tt*bb)*(1-rr)*ca0*sc0/(2*mgap*delta)
    cact=ca0*mass/(2*(1-hh)*delta)
    aact=tt*ca0*(1+aa*rr)*sc0/(2*hh*delta)
    dact=(1+aa*bb)*(1+aa*rr)*(1-aa*tt*bb*rr)*sc0/(2*tt*dgap*delta)
    alact=bb**2*(1+aa*rr)*mass/(2*(1-hh*bb)*delta)
    xx=ss*qq; eta=1/ss; kk=1+xx; R=1+xx*hh; Tn=1+xx*rr
    GG=rr+zz; MM=R-ww; JJ=Tn-ss*(1-hh); EE=1-qq*rr; QQ=1-xx*zz*hh
    P4=Tn*EE+xx*(1-hh)*GG
    Pal=R*rr*(1+qq*zz)+zz*MM*EE
    pd=1+2*xx*rr+hh*xx**2*GG+xx*zz+eta*xx*(kk*zz+xx*rr*ee)
    cf=kk*(hh+eta)*EE/(2*(1-hh)*pd)
    af=eta*P4/(2*hh*pd)
    mf=kk*xx*(hh+eta)*rr*(1+qq*zz)/(2*MM*zz*pd)
    alf=xx**2*zz*(hh+eta)*Pal/(2*QQ*pd)
    df=(kk*(1+qq*zz)*(hh+eta)-(eta*Tn-1+hh)*pd)/(2*JJ*eta*pd)
    betaf=(kk*xx*(hh+eta)*GG-ee*pd)/(2*ee*pd)
    PP=1+qq*zz+ss*qq*(2+zz-2*ee)+ss*qq**2*(zz+ee-ee**2)+ss**2*hh*qq**2*(1+zz-ee)
    GG2=1+zz-ee
    PC=1+qq*zz+ss*(1+qq*(2+2*zz-ee)+qq**2*zz)+ss**2*(qq*(hh*(1+zz)+1-ee)+qq**2*(hh*ee*zz-ee*(1+zz)+2*zz+1))+ss**3*hh*qq**2*GG2
    PA=1+qq*zz+ss*(qq*(2+zz-2*ee)+qq**2*zz*(1-ee))+ss**2*qq**2*GG2*(GG2-hh*zz)
    PD=1+qq*(zz+ee)+ss*(qq*(hh*(1+zz)+1-ee)+qq**2*(zz+ee-ee**2))+ss**2*hh*qq**2*GG2
    PM=1+qq*zz+ss*(qq*(2+zz-2*ee)+qq**2*zz)+ss**2*qq**2*((1-ee)**2+zz*(1-ee)+hh*ee*zz)
    PB=(1+zz)*(1+qq*zz)+ss*qq*(1+zz)*(2+zz-2*ee)+ss*qq**2*zz*(1+zz-ee**2)+ss**2*qq**2*GG2**2
    check('TU6 common formal polynomial',pd-PP)
    lhs=[(1-hh)*ss*(cact-cf),hh*ss*(aact-af),(dact-df)/ss,(mact-mf)*zz*MM,ee/2-(betaf+half)*ee,(alact-alf)/(ss*zz)]
    rhs=[-nact*rr*PC/(zz*PP),-nact*rr*PA/(zz*PP),-nact*rr*(1+qq*zz)*PD/(zz*JJ*PP),nact*rr*(1+qq*zz)*PM/(zz*PP),nact*rr*PB/(zz*PP),nact*rr*qq**2*(1+qq*zz)*((1+ss)*ee-ss*(1-hh)*(1+zz))/(zz*QQ*PP)]
    for label,left,right in zip(['C','A','D','m','beta','alpha'],lhs,rhs):
        check('TU9 native comparison '+label,left-right)
    check('TU10 exact z-divisibility numerator', (PB-(1+qq*zz)*PM).subs(zz,0))
    # UL7 and UL13 with source openings retained before specialization.
    UU,VV=s.symbols('Ut Vt'); Du,Dv=1+UU,1+VV
    dc=-rr*PC/((1-hh)*PP); da=-rr*PA/(hh*PP)
    dd=-rr*(1+qq*zz)*PD/(JJ*PP); dm=rr*(1+qq*zz)*PM/PP; db=rr*PB/PP
    dal=rr*qq**2*(1+qq*zz)*((1+ss)*ee-ss*(1-hh)*(1+zz))/(QQ*PP)
    rho=kk*hh/R; at=1-rho+(ww/R)*VV
    FA=4*Dv*R**2*(dc*(1-rho)**2+da*rho**2/Du)+4*dal*zz*ss**2*VV*R**2*rho**2/Du+4*dd*ss**2*R**2*at**2
    FB=-8*dal*zz*ss*R*rho-8*dd*ss*R*at*Du
    FC=4*dd*Du**2
    FL=4*R*(dm*(1-ww/R)+db*(ww/R)*Dv)
    zero={qq:0,UU:0,VV:0}
    check('UL7 first specialized cancellation',(FA-FC).subs(zero))
    check('UL7 second specialized cancellation',(FB+2*FC+8*Du*Dv*rr).subs(zero))
    check('UL7 third specialized cancellation',(FL-4*Du*Dv*(1-ee**2)).subs(zero))
    FH=Dv*(rho*da-(1-rho)*dc)-Dv*da*UU*rho/Du+dal*zz*ss**2*VV*rho/Du-dd*ss**2*at
    FV=dd*ss*Du-dal*zz*ss
    check('UL13 H specialized cancellation',FH.subs(zero)-rr*ss/(1-ss*(1-hh)))
    check('UL13 v specialized cancellation',FV.subs(zero)+rr*ss/(1-ss*(1-hh)))
    FN=(PB*Dv-(1+qq*zz)*PM)/PP-zz*Du*Dv
    check('UL10 complete N polynomial at q0',FN.subs(qq,0)-VV+zz*UU*Dv)

if '--r67' in sys.argv or '--all' in sys.argv:
    om,pp,cc,ee,es,ds=s.symbols('omega p c epsilon E Drow')
    mu1,mu2,n1,n2,ca,aa,dd,al,acap,ah,bd,hh=s.symbols('mu1 mu2 n1 n2 C A D alpha Kcap Ah betaD h')
    Du,Dv=1+u*u,1+v*v
    qa=Du*v*om; XX=N+pp*qa-u; qb=v*(1-u*XX)
    WW=om+u*(M+cc*qb)/Du; YY=u*om-(M+cc*qb)/Du
    JJ=M-(1-cc)*qb
    base=-2*mu1*qa*N+mu2*qa**2-2*n1*Du*Dv*XX*v*om-n2*Du*Dv*(v*om)**2+(al+half)*Du*(v*YY-1)**2+half*Dv*(XX+u)**2-half*Du*Dv
    gri=acap*N*N+base+ca*Dv*(M+cc*qb)**2+aa*Du*Dv*YY**2+(dd+half)*JJ**2
    AR=ee*om*v; BR=v; CR=es*u*om; V0=1-pp*u*v*om
    Hh=al+half; LL=aa*Dv+Hh*v*v
    F0=pp*pp/2-n2-2*n1*pp+mu2; AN=pp-2*mu1-2*n1
    fA=Du**2*Dv*F0/ee**2
    fBB=(ca*cc**2*Du**2*Dv+(dd+half)*(1-cc)**2*Du**2+Du*LL*cc**2)*V0**2+2*Hh*Du*cc*V0-Du/2-Du**2*mu2*v*v*om*om
    fBC=(-2*Du*LL*cc*V0-2*Hh*Du+2*n1*Du*Dv)/es
    fCC=(Du*LL+Du*Dv*n2*v*v)/es**2
    ftab=fA*AR**2+fBB*BR**2+fBC*BR*CR+fCC*CR**2+al*Du
    check('RI8 RI9 complete raw scalar',gri.subs({N:0,M:0})-ftab)
    RA=Du*AN/ee
    RB=Du*(pp-2*n1)*v*v*om-2*u*v*(ca*cc**2*Du*Dv+(dd+half)*(1-cc)**2*Du+aa*cc**2*Dv)*V0-2*Hh*cc*u*v*(1+cc*v*v*V0)
    RC=2*cc*u*v*LL/es
    SB=2*(Du*(ca*cc*Dv-(dd+half)*(1-cc))+cc*LL)*V0+2*Hh
    SC=-2*LL/es
    check('RI11 RI12 entire first N normal',s.diff(gri,N).subs({N:0,M:0})-(RA*AR+RB*BR+RC*CR))
    check('RI11 RI12 entire first M normal',s.diff(gri,M).subs({N:0,M:0})-(SB*BR+SC*CR))
    check('RI13 full hard N',s.diff(gri,N,2).subs({u:0,v:0})/2-acap-half)
    check('RI13 full hard M',s.diff(gri,M,2).subs({u:0,v:0})/2-aa-ca-dd-half)
    gru=base+ca*Dv*(M+cc*qb)**2+2*ah*Du*Dv*u*WW*YY-ah*hh*Du*Dv*(u*WW)**2-2*bd*M*qb+bd*(1-cc)*qb**2+JJ**2/2
    fBB2=(ca*cc**2*Du**2*Dv+(bd*(1-cc)+(1-cc)**2/2)*Du**2+Hh*Du*v*v*cc**2-ah*Du*Dv*(2*u*u+hh*u**4)*cc**2)*V0**2+2*Hh*Du*cc*V0-Du/2-Du**2*mu2*v*v*om*om
    fBC2=(2*ah*Du*Dv*((1-hh)*u*u-1)*cc*V0-2*Hh*Du*v*v*cc*V0-2*Hh*Du+2*n1*Du*Dv)/es
    fCC2=(ah*(2-hh)*Du*Dv+Hh*Du*v*v+Du*Dv*n2*v*v)/es**2
    check('RU13 RU14 complete raw scalar',gru.subs({N:0,M:0})-(fA*AR**2+fBB2*BR**2+fBC2*BR*CR+fCC2*CR**2+al*Du))
    RB2=Du*(pp-2*n1)*v*v*om-2*u*v*(ca*cc**2*Du*Dv+(bd*(1-cc)+(1-cc)**2/2)*Du)*V0+2*ah*cc**2*u*v*Dv*u*u*(2+hh*u*u)*V0-2*Hh*cc*u*v*(1+cc*v*v*V0)
    RC2=(2*ah*cc*u*v*Dv*(1-(1-hh)*u*u)+2*Hh*cc*u*v**3)/es
    SB2=2*(ca*cc*Du*Dv-(bd+(1-cc)/2)*Du-ah*cc*Dv*u*u*(2+hh*u*u))*V0+2*Hh*(1+cc*v*v*V0)
    SC2=(2*ah*Dv*((1-hh)*u*u-1)-2*Hh*v*v)/es
    check('RU15 RU16 entire first N normal',s.diff(gru,N).subs({N:0,M:0})-(RA*AR+RB2*BR+RC2*CR))
    check('RU15 RU16 entire first M normal',s.diff(gru,M).subs({N:0,M:0})-(SB2*BR+SC2*CR))
    check('RU full hard N',s.diff(gru,N,2).subs({u:0,v:0})/2-half)
    check('RU full hard M',s.diff(gru,M,2).subs({u:0,v:0})/2-ca-half)

if '--r67source' in sys.argv or '--all' in sys.argv:
    ss,bs,gs,hs,ts=s.symbols('s b g h t')
    av=1/ss; bv=ss**2*bs; xs=ss*bs; ws=1+hs*xs-gs; rv=1-ss*ws
    KS=1+xs; RS=1+hs*xs; TS=1+xs*rv; ls=1-hs
    SC=1+bv*(ls+av*rv)
    DE=s.cancel(av+(1+xs)+(1+xs)*bv+2*xs*(1+xs)*rv+av**3*bv**2*rv**2+2*av**2*bv*rv+ts*bv*(1+av+xs+av**2*bv*rv*(2-rv)))
    MASS=s.cancel((1+xs)*(1+av*rv)*(hs+ts)-av*(1+ts*bv)*(1-rv)*ls)
    NUM=(1+av*rv)*TS*SC-rv*DE
    MG=s.cancel(av**2*(1+ts*bv)*(1-rv)*TS*SC/(2*DE))
    NN=s.cancel(NUM/(2*rv*DE))
    NG=s.cancel(NUM*KS/(2*ss*RS*DE))
    NG2=s.cancel(NUM*rv*KS**2/(2*ss**2*RS**2*DE))
    CS=s.cancel(TS*MASS/(2*ls*DE))
    AS=s.cancel(ts*TS*(1+av*rv)*SC/(2*hs*DE))
    DS=s.cancel(KS*(1+av*rv)*(1-av*ts*bv*rv)*SC/(2*ts*(TS-ls/ts)*DE))
    AL=s.cancel(bv**2*(1+av*rv)*MASS/(2*(1-hs*bv)*DE))
    ccS=ls/(ts*TS); ppS=ws/RS
    BD=s.cancel(DS*(1-ccS))
    AN=s.cancel(ppS-2*MG/RS-2*NG)
    FZERO=s.cancel(ppS**2/2-NG2-2*NG*ppS+MG*gs/RS**2)
    # Exact denominator-safe divisions: numerator divisibility is polynomial.
    anquot=s.cancel(AN/(bs*ss)); fquot=s.cancel(FZERO/bs); alqu=s.cancel(AL/(ss**3*bs**2))
    for label,expr,expect in [('AN/(bs)',anquot,hs+ts),('F0/b',fquot,(hs+ts)/2),('alpha normalized',alqu,(hs+ts)/2)]:
        num,den=expr.as_numer_denom()
        assert den.subs({ss:0,bs:0,gs:0})!=0
        check('R6/R7 source '+label+' leading',expr.subs({ss:0,bs:0,gs:0})-expect)
    check('RB n/s2 physical residual',s.cancel(NN/ss**2).subs({ss:0,gs:0})-(1-bs*(hs+ts))/2)
    def taylor2(expr):
        num,den=s.fraction(s.cancel(expr.subs(gs,0)))
        pn,pd=s.Poly(num,ss,bs),s.Poly(den,ss,bs)
        out={}; d00=pd.coeff_monomial(1)
        for total in range(3):
            for i in range(total+1):
                j0=total-i
                val=pn.coeff_monomial(ss**i*bs**j0)
                for (a0,b0),dval in pd.terms():
                    if a0+b0 and a0<=i and b0<=j0:
                        val-=dval*out.get((i-a0,j0-b0),0)
                out[(i,j0)]=s.cancel(val/d00)
        return sum(val*ss**i*bs**j0 for (i,j0),val in out.items())
    jets={
       'A':(AS,ts*(1-ss+ss**2)/(2*hs)),
       'C':(CS,(hs+ts-(ts+1)*ss+(ts+1)*ss**2)/(2*ls)),
       'D':(DS,(1-ss+ss**2)/(2*(ts-ls))-ss*bs*ts*(hs+ts)/(2*(ts-ls)**2)),
       'ngamma':(NG,(ss-ss**2-ss*bs*(hs+ts))/2),
       'Ah at h0':(s.cancel(hs*AS).subs(hs,0),ts*(1-ss+ss**2)/2),
       'betaD at h0':(BD.subs(hs,0),(1-ss+ss**2-(ts+1)*ss*bs)/(2*ts)),
    }
    for label,(expr,expect) in jets.items():check('RI18 RU22 exact Taylor '+label,taylor2(expr)-expect)
    def low_degree_zero(name,expr,variables,degree):
        num=s.together(expr).as_numer_denom()[0]
        terms=s.Poly(s.expand(num),*variables).terms()
        low=sum(coef*s.prod(var**power for var,power in zip(variables,mon)) for mon,coef in terms if sum(mon)<=degree)
        check(name,low)
    # RI chosen raw allocation, through its full required angular weight four.
    # Each of s,b,u has weight two; omitted v and u*v*omega have weight >=6.
    Aj=jets['A'][1]; Cj=jets['C'][1]; Dj=jets['D'][1]; Nj=jets['ngamma'][1]
    ccj=ls*(1-ss*bs)/ts; Ej=ts*(1+ls*ss*bs); duj=1+u*u
    BBj=Cj*ccj**2*duj**2+(Dj+half)*(1-ccj)**2*duj**2+duj*Aj*ccj**2+duj*ccj-duj/2
    BCj=(-2*duj*Aj*ccj-duj+2*Nj*duj)/Ej
    CCj=duj*Aj/Ej**2
    kap=(hs+ts)/(2*ts)
    Lam=(1-ss+ss**2+u*u-(2*ls+hs*(hs+ts))*ss*bs)/(2*hs*ts)
    low_degree_zero('RI chosen fBB filtration',BBj-Lam-kap*u*u,(ss,bs,u),2)
    low_degree_zero('RI chosen fBC filtration',BCj+2*Lam+4*kap*ss*bs,(ss,bs,u),2)
    low_degree_zero('RI chosen fCC filtration',CCj-Lam-kap*ss*bs,(ss,bs,u),2)
    # RU chosen raw allocation at h=g=0, with the two divergent packets removed.
    Cj=jets['C'][1].subs(hs,0); Ahj=jets['Ah at h0'][1]; Bj=jets['betaD at h0'][1]
    Nj=jets['ngamma'][1].subs(hs,0); ccj=(1-ss*bs)/ts; Ej=ts*(1+ss*bs); qj=1-ccj
    BBj=Cj*ccj**2*duj**2+(Bj*qj+qj**2/2)*duj**2-2*Ahj*duj*u*u*ccj**2+duj*ccj-duj/2
    BCj=(2*Ahj*duj*(u*u-1)*ccj-duj+2*Nj*duj)/Ej
    CCj=2*Ahj*duj/Ej**2
    Lam=1/(ts*(1+ss)); ksq=ss*bs
    low_degree_zero('RU small-b chosen fBB filtration',BBj-Lam-u*u/2+(1/ts+half)*ksq,(ss,bs,u),2)
    low_degree_zero('RU small-b chosen fBC filtration',BCj+2*Lam+u*u/ts-(3/ts-1)*ksq,(ss,bs,u),2)
    low_degree_zero('RU small-b chosen fCC filtration',CCj-Lam-u*u/ts+2*ksq/ts,(ss,bs,u),2)

    # Complementary positive-unit b chart: weight(s)=2, weight(u)=1.
    ER=ts*KS/RS; invd2=bs*ts*KS**2/RS**2; cR=ccS
    BRB=(CS*cR**2*(1+u*u)**2+(DS+half)*(1-cR)**2*(1+u*u)**2+(1+u*u)*AS*cR**2+(1+u*u)*cR-(1+u*u)/2)*invd2
    BCRB=(-2*(1+u*u)*AS*cR-(1+u*u)+2*NG*(1+u*u))*invd2/ER
    CCRB=(1+u*u)*AS*invd2/ER**2
    kaprb=bs*(hs+ts)/2; lamrb=bs*(1+u*u)/(2*hs)-ss*bs*(1+bs*hs*(hs+ts))/(2*hs)
    def rbweighted(name,expr):
        num=s.fraction(s.cancel(expr.subs(gs,0)))[0]
        low=sum(coef*ss**mon[0]*u**mon[1] for mon,coef in s.Poly(num,ss,u).terms() if 2*mon[0]+mon[1]<=2)
        check(name,low)
    rbweighted('RB unit-b chosen fBB filtration',BRB-lamrb-kaprb*u*u)
    rbweighted('RB unit-b chosen fBC filtration',BCRB+2*lamrb+4*kaprb*ss*bs)
    rbweighted('RB unit-b chosen fCC filtration',CCRB-lamrb-kaprb*ss*bs)

if '--wa' in sys.argv or '--all' in sys.argv:
    aa,bb,tt,rr,hh=s.symbols('as bs ts rs hs')
    mw,nw,cw,aw,dw,alw=s.symbols('mw nw cw aw dw alw')
    WW,XX,YY,ZZ=s.symbols('WW XX YY ZZ')
    UU,VV=u*u,v*v; du,dv=1+u*u,1+v*v
    ll=1-hh; gapj=1+aa*bb*rr-ll/tt
    ca=1+aa*bb*rr; sc=1+bb*(1-hh+aa*rr); gm=1+aa*(hh*bb+rr-1)
    de=aa+(1+aa*bb)+(1+aa*bb)*bb+2*aa*bb*(1+aa*bb)*rr+aa**3*bb**2*rr**2+2*aa**2*bb*rr+tt*bb*(1+aa+aa*bb+aa**2*bb*rr*(2-rr))
    mass=(1+aa*bb)*(1+aa*rr)*(hh+tt)-aa*(1+tt*bb)*(1-rr)*(1-hh)
    nwactual=((1+aa*rr)*ca*sc-rr*de)/(2*rr*de)
    muactual=(1+tt*bb)*(1-rr)*ca*sc/(2*gm*de)
    Cactual=ca*mass/(2*(1-hh)*de)
    Ahactual=tt*ca*(1+aa*rr)*sc/(2*de)
    Dactual=(1+aa*bb)*(1+aa*rr)*(1-aa*tt*bb*rr)*sc/(2*tt*gapj*de)
    asactual=(1+aa*rr)*mass/(2*(1-hh*bb)*de)
    actual={nw:nwactual,mw:aa**2*muactual,cw:Cactual,aw:Ahactual/hh,dw:Dactual,alw:bb*bb*asactual}
    root=s.Matrix([XX+rr*u*ZZ,YY-hh*u*WW,u*WW-v*ZZ/tt,aa*v*WW-u*ZZ,u*v*WW-bb*ZZ])
    ast=(cw+dw+half)*ll**2+aw*hh**2; a4=aw*hh**2
    b2=2*mw*hh+cw*ll**2+(aw+alw+half)*hh**2
    b4=hh**2*(mw+aw+alw+half)
    c4=nw*rr**2; d2=(nw+half)*rr**2+2*dw*rr; d4=(nw+dw+half)*rr**2
    J=(aa*bb*d4-a4/tt**2+d2+d4*UU)*tt/2
    beta=(bb*d4+c4/tt)*tt
    fs=(b2-alw/bb**2-aa*aa*c4-aa*aa*bb*tt*d4)/(2*aa)
    j1=d4*tt/2
    gamma=b4+d4*aa*tt
    SS=s.Matrix([[ast+a4*UU,u*v*J,-aa*v*J],[u*v*J,mw/aa**2+beta*UU,u*(fs-j1*UU)],[-aa*v*J,u*(fs-j1*UU),alw/bb**2+gamma*UU]])
    p0=tt*(dw+rr*(nw+half)); p1=tt*rr*(nw+dw+half)
    q1=aa*(nw*rr+bb*p1); fx=-mw*hh+(dw+half)*ll-p0-q1
    LX=s.Matrix([[u*v*(p0+p1*UU),(-mw+q1*UU)/aa,u*(fx-p1*UU)]])
    py=-(cw+dw+half)*ll+aw*hh; pyp=aw*hh
    kv=-cw*ll+aw*hh+mw+(alw+half)*hh; kuv=hh*(aw+mw+alw+half)
    kz=-mw*(1-rr)+rr*(dw+half)-alw-half
    qy=aa*(-kz-pyp/tt-bb*kuv); fy=kv-qy
    LY=s.Matrix([[py+pyp*UU,u*v*qy/aa,v*(fy+kuv*UU)]])
    hxx=nw*du*dv+mw+(dw+half)*UU*VV+half*dv
    hxy=-u*v*(mw+dw+half)
    hyy=aw*du*dv+cw*dv+dw+half+mw*UU*VV+(alw+half)*du*VV
    GG=s.Matrix.vstack(s.Matrix.hstack(s.Matrix([[hxx,hxy],[hxy,hyy]]),s.Matrix.vstack(LX,LY)),s.Matrix.hstack(s.Matrix.hstack(LX.T,LY.T),SS))
    pp=mw*(v*WW-XX+u*v*YY-u*ZZ)**2+nw*du*dv*XX**2+cw*dv*(u*WW-YY)**2+aw*du*dv*YY**2+(dw+half)*(u*WW+u*v*XX-YY-v*ZZ)**2+(alw+half)*du*(v*YY-ZZ)**2+half*dv*(XX+u*ZZ)**2-half*du*dv*ZZ**2
    diff=s.Poly(s.expand(pp-(root.T*GG*root)[0]),u,v,WW,XX,YY,ZZ)
    uniques={s.factor(cc) for cc in diff.coeffs()}
    print('WA generic residual coefficient count',len(uniques),flush=True)
    for ii,cc in enumerate(uniques): check('WA complete parent exact source residual '+str(ii),cc.subs(actual))
    dr=s.Matrix([u*v,-1,ll,rr*u*v,-aa*rr*v])
    sy=s.Matrix([0,0,(aa*bb-UU)*v,-(bb-VV/tt)*u,UU-aa*VV/tt])
    kr=s.Matrix([-2*rr*tt*u,0,0,-rr*rr*tt*u,aa*rr*rr*tt])
    HD=s.zeros(5)
    HD[0,2]=HD[2,0]=tt*u*v; HD[0,4]=HD[4,0]=-tt*u
    HD[2,3]=HD[3,2]=rr*tt*u*v; HD[2,4]=HD[4,2]=-aa*rr*tt*v
    # R5 exact centered m gauge on the unscaled second scalar root.
    scaleM=s.diag(1,1,1,aa,1); GM=scaleM*GG*scaleM
    lm=s.Matrix([-1,u*v,0,1-gm,hh*u]); HM=s.zeros(5)
    HM[0,3]=HM[3,0]=-1; HM[1,4]=HM[4,1]=v
    HM[3,3]=2-gm; HM[3,4]=HM[4,3]=hh*u
    check('R5 complete centered m matrix gauge',sum(s.expand(e)**2 for e in GM.diff(mw)-lm*lm.T-gm*HM))
    check('RD physical syzygy',(sy.T*root)[0])
    check('RD centered D row',(dr.T*root)[0]-(u*WW+u*v*XX-YY-v*ZZ+gapj*v*ZZ))
    md=GG.diff(dw)-dr*dr.T-(sy*kr.T+kr*sy.T)/2-gapj*HD
    check('RD entire D gauge sum of entry squares',sum(e*e for e in md))
    bsactual=2*aa**2*muactual*hh+Cactual*ll**2+Ahactual*hh+(bb*bb*asactual+half)*hh**2
    QQ=s.cancel((bsactual-asactual)/aa)
    check('WE exact monic a divisibility',s.together(bsactual-asactual).as_numer_denom()[0].subs(aa,0))
    check('SL RD divided a source corner',QQ.subs({bb:0,hh:0,rr:0}))
    nn0=(1-rr)/(2*rr*(1+aa)); mm0=aa**2*(1-rr)/(2*(1-aa*(1-rr))*(1+aa))
    cc0=(tt*(1+aa*rr)-aa*(1-rr))/(2*(1+aa)); dd0=(1+aa*rr)/(2*(1+aa)*(tt-1))
    p3=s.expand(pp.subs({YY:0,nw:nn0,mw:mm0,cw:cc0,dw:dd0,alw:0,aw:0}))
    P3=s.hessian(p3,(WW,XX,ZZ))/2
    adet=dd0*mm0*(2*nn0+1)/2; cdet=cc0*(2*mm0+1)*(2*dd0+2*nn0+1)/4; ddet=cc0*dd0*(2*nn0+1)/2
    check('WR complete remote determinant',(P3[0,0]*P3[1,1]*P3[2,2]+2*P3[0,1]*P3[0,2]*P3[1,2]-P3[0,0]*P3[1,2]**2-P3[1,1]*P3[0,2]**2-P3[2,2]*P3[0,1]**2)-du*dv*(adet*(VV-tt*UU/aa)**2+cdet*UU**2*VV+ddet*UU*VV**2))
    nrd=(1-rr)/(2*rr*(1+aa)); mrd=aa**2*(1-rr)/(2*(1-aa*(1-rr))*(1+aa)); crd=(1+aa*(2*rr-1))/(2*(1+aa))
    # JD20: restriction after Y=J=0, no hidden normal relaxation.
    p2=(nrd*du*dv+mrd*dv**2+crd*dv*UU*VV+dv/2)*XX**2+2*(-mrd*dv*(VV-UU)/u-crd*dv*u*VV+u*dv/2)*XX*ZZ+(mrd*(VV-UU)**2/UU+crd*dv*VV+(UU-VV)/2)*ZZ**2
    p2orig=pp.subs({YY:0,WW:v*ZZ/u-v*XX,nw:nrd,mw:mrd,cw:crd,alw:0,aw:0})
    check('JD20 full restricted parent',p2orig-p2)
    gmr=1-aa*(1-rr)
    detjd=du*dv*(1+aa*rr)*( (1-rr)*(UU-aa*VV)**2+rr*(gmr+aa*rr)*UU**2*VV+gmr*(gmr+aa*rr)*UU*VV**2 )/(4*rr*(1+aa)**2*gmr)
    check('JD21 complete remote determinant',UU*(s.diff(p2,XX,2)*s.diff(p2,ZZ,2)-s.diff(p2,XX,ZZ)**2)/4-detjd)
    # R5 remote three-state determinant at arbitrary bounded-a residue.
    m5=aa**2/(2*(1-aa**2)); a5=tt/(2*hh*(1+aa)); c5=(hh+tt-aa*(1-hh))/(2*(1-hh)*(1+aa)); d5=1/(2*(1+aa)*(tt-1+hh))
    p5=pp.subs({XX:0,nw:0,mw:m5,aw:a5,cw:c5,dw:d5,alw:0})
    P5=s.hessian(p5,(WW,YY,ZZ))/2
    det5=P5[0,0]*P5[1,1]*P5[2,2]+2*P5[0,1]*P5[0,2]*P5[1,2]-P5[0,0]*P5[1,2]**2-P5[1,1]*P5[0,2]**2-P5[2,2]*P5[0,1]**2
    expected5=du*dv*((hh+tt)*(tt*UU-aa*VV)**2+(hh+tt-aa*(1-hh))*(tt*(tt-1+hh)*UU**2*VV+(hh+(1-aa)*tt)*UU*VV**2))/(8*hh*(1-aa)*(1+aa)**3*(1-hh)*(tt-1+hh))
    check('R5.8 complete remote determinant',det5-expected5)
    # Entire Bell reflection used for R2 outer tests, including both zero events.
    la,be=s.symbols('lambda beta')
    general=pp+la*du*(WW+v*XX)**2+be*dv*(XX+u*ZZ)**2
    w2=u*v*WW-u*XX-v*YY+ZZ; x2=-v*WW+XX-u*v*YY+u*ZZ
    y2=-u*WW-u*v*XX+YY+v*ZZ; z2=WW+v*XX+u*YY+u*v*ZZ
    reflect={u:v,v:u,WW:w2,XX:x2,YY:y2,ZZ:z2,la:alw,alw:la,mw:nw,nw:mw,aw:dw,dw:aw}
    check('R2 entire Bell reflection',general.xreplace(reflect)-du*dv*general)


if '--rc' in sys.argv or '--all' in sys.argv:
    aa,bb,tt,rr,hh=s.symbols('as bs ts rs hs')
    mw,nw,cw,aw,dw,alw=s.symbols('mw nw cw aw dw alw')
    WW,XX,YY,ZZ,EE,NN,XXI=s.symbols('WW XX YY ZZ EE NN XXI')
    UU,VV=u*u,v*v; du,dv=1+u*u,1+v*v
    ev=1-rr; xs,ys=aa*bb,tt*bb
    ca=1+aa*bb*rr; sc=1+bb*(1-hh+aa*rr); gm=1+aa*(hh*bb+rr-1); gapj=ca-(1-hh)/tt
    de=aa+(1+aa*bb)+(1+aa*bb)*bb+2*aa*bb*(1+aa*bb)*rr+aa**3*bb**2*rr**2+2*aa**2*bb*rr+tt*bb*(1+aa+aa*bb+aa**2*bb*rr*(2-rr))
    mass=(1+aa*bb)*(1+aa*rr)*(hh+tt)-aa*(1+tt*bb)*(1-rr)*(1-hh)
    nwactual=((1+aa*rr)*ca*sc-rr*de)/(2*rr*de)
    muactual=(1+tt*bb)*(1-rr)*ca*sc/(2*gm*de)
    Cactual=ca*mass/(2*(1-hh)*de)
    Ahactual=tt*ca*(1+aa*rr)*sc/(2*de)
    Dactual=(1+aa*bb)*(1+aa*rr)*(1-aa*tt*bb*rr)*sc/(2*tt*gapj*de)
    asactual=(1+aa*rr)*mass/(2*(1-hh*bb)*de)
    actual={nw:nwactual,mw:aa**2*muactual,cw:Cactual,aw:Ahactual/hh,dw:Dactual,alw:bb*bb*asactual}
    qs=nw*rr**2+(mw+half)*ev**2; qq=qs+alw; muw=mw-aa**2*qq; lam=muw/(tt*aa**2*bb)
    ast=(cw+dw+half)*(1-hh)**2+aw*hh**2
    b2=2*mw*hh+cw*(1-hh)**2+(aw+alw+half)*hh**2
    b4=hh**2*(mw+aw+alw+half)
    d2=(nw+half)*rr**2+2*dw*rr; d4=(nw+dw+half)*rr**2
    e0=mw*ev+dw*(1-hh)+alw*hh+half
    e1=mw*hh*ev+rr*(dw+half)*(1-hh)+(alw+half)*hh
    cF=(b2-aa**2*nw*rr**2)/tt**2
    fA0=aw*hh**2/tt**2+d2-2*e1/tt+2*aa*nw*rr**2/tt
    fC0=(alw-cF*tt**2*bb**2)/(tt*aa*bb**2)
    fG=-(ast/tt**2+dw-2*e0/tt+2*aa*qq/tt+2*cF*tt*bb)/(aa*bb)
    aE=2*aw*hh**2/tt**2-2*e1/tt; bE=2*b2/tt**2; cE=2*b4/tt**2
    aY=2*(aw*hh/tt-mw*ev+rr*(dw+half)-alw-half)
    bY=2*(-cw*(1-hh)+aw*hh+mw+(alw+half)*hh)/tt
    cY=2*(aw*hh+hh*(mw+alw+half))/tt
    x3=2*((dw+half)*((1-hh)/tt-1)+ev/2-nw*rr*(1+aa/tt)-mw*hh/tt)
    x5=-2*rr*(nw+dw+half)
    KX=-nw*rr+(mw+half)*ev-aa*qq
    KE=ast/tt**2-e0/tt+aa*qq/tt
    KY=(-(cw+dw+half)*(1-hh)+aw*hh)/tt+dw-alw
    balances=[fA0-lam+2*d4*xs+b4*ys/tt**2,fC0-lam-b4*ys/tt**2,fG-3*lam+d4*xs-b4*ys/tt**2,2*KE+aE*xs+bE*ys+cE*xs*ys+2*lam*xs,2*KX+x3*ys+x5*xs*ys-2*tt*lam*xs,2*KY+aY*xs+bY*ys+cY*xs*ys]
    for ii,cc in enumerate(balances): check('RC1 RC6 exact contact balance '+str(ii),cc.subs(actual))
    lamactual=s.cancel(lam.subs(actual)); cfactual=s.cancel(cF.subs(actual))
    check('bounded-a lambda joint residue',lamactual.subs({bb:0,hh:0,rr:1})-half)
    check('bounded-a cF joint residue',cfactual.subs({bb:0,hh:0,rr:1})-1/(2*tt))
    check('bounded-a alpha normalized residue',asactual.subs({bb:0,hh:0,rr:1})-tt*(1+aa)/2)
    print('RC normalized lambda denominator at bounded corner',s.factor(s.denom(lamactual).subs({bb:0,hh:0,rr:1})),flush=True)
    pp=mw*(v*WW-XX+u*v*YY-u*ZZ)**2+nw*du*dv*XX**2+cw*dv*(u*WW-YY)**2+aw*du*dv*YY**2+(dw+half)*(u*WW+u*v*XX-YY-v*ZZ)**2+(alw+half)*du*(v*YY-ZZ)**2+half*dv*(XX+u*ZZ)**2-half*du*dv*ZZ**2
    frame={WW:(EE+v*ZZ)/(tt*u),XX:XXI-rr*u*ZZ,YY:hh*(EE+v*ZZ)/tt+NN}
    parent=s.expand(pp.subs(frame))
    QR=aa*v*frame[WW]-u*ZZ-aa*XXI; NR=rr*(aa*VV/tt-UU)*ZZ+u*XXI
    remainder=parent-qq*QR**2-nw*NR**2-muw*(v*EE/(tt*u)-XXI)**2
    H0=s.expand(remainder.subs(ZZ,0))
    assert not s.denom(s.together(H0)).has(u,v,aa)
    print('RC full H0 coefficients have no testing inverse or a denominator',flush=True)
    # Use c^2=a/t to remove square roots in the full generic identity.
    cc=s.symbols('c'); arel=tt*cc**2
    RD=v*(UU-xs)*ZZ; FV=(VV-ys)*ZZ
    JR=(u*v-xs/cc)*ZZ
    R2=lam*tt**2*xs*(cc*v-u)**2*frame[WW]**2
    T2=2*lam*xs*u*v*(tt*cc*frame[WW]-ZZ)**2/cc
    circuitSum=lam*(tt**2*xs*(UU+aa*VV/tt)*frame[WW]**2-4*tt*xs*u*v*frame[WW]*ZZ+(UU*VV+xs*ys)*ZZ**2)
    check('RC positive circuit sum',(T2+R2+lam*JR**2-circuitSum).subs(aa,arel))
    rhs=qq*QR**2+nw*NR**2+circuitSum+d4*RD**2+(cF+b4*xs/tt**2)*FV**2+b4*v*RD*FV/tt**2+H0-lam*xs*EE**2-muw*XXI**2-2*tt*lam*xs*QR*XXI+EE*((aE+cE*VV)*RD+(bE+cE*xs)*v*FV)+NN*((aY+cY*VV)*RD+(bY+cY*xs)*v*FV)+XXI*(x5*u*v*RD+u*(x3+x5*xs)*FV)
    residual=s.cancel(parent-rhs)
    pol=s.Poly(s.expand(s.fraction(residual)[0]),EE,NN,XXI,ZZ,u,v)
    uniques={s.factor(co) for co in pol.coeffs()}
    print('R10.9 residual coefficient count',len(uniques),flush=True)
    actrel=actual
    for ii,coef in enumerate(uniques):check('R10.9 full parent source residual '+str(ii),coef.subs(actrel))
    # Entire divergent-D channel before source specialization.
    R10form=rhs-qq*QR**2-nw*NR**2-T2
    LD=(1-hh)*EE/tt-rr*RD-NN+u*v*XXI
    check('JD4 complete D channel',s.diff(R10form,dw)-LD**2+2*gapj*u*FV*XXI)
    # Source polynomial positivity and new compact bounded-a Schur lower bound.
    Qsrc=1+bb*(1+tt+aa*hh+2*aa*rr+aa*tt)+bb**2*(aa*(hh+tt)+aa**2*(hh*(rr+1)+tt*(1+rr-rr**2)))
    Ksrc=(1+aa)*(hh+tt)*(1+aa*bb)
    check('BS1 exact physical source scale',2*rr*de*nwactual-ev*Qsrc+bb*Ksrc)
    bound=half-tt*aa**2/(2*(1+aa)*(1+tt*(1+aa)))
    check('new bounded-a full-parent Schur bound',bound-(1+aa+tt*(1+2*aa))/(2*(1+aa)*(1+tt*(1+aa))))

if '--tl' in sys.argv or '--all' in sys.argv:
    ss,bb,ww,hh,tt=s.symbols('s b w h t')
    ee=ss*ww; rr=1-ee; aa=1/ss; bv=ss**2*bb; xx=ss*bb
    kk=1+xx; ca=1+xx*rr; sc=1+bv*(1-hh)+xx*rr; mgap=1-ww+hh*xx; gj=ca-(1-hh)/tt
    delt=s.cancel(aa+(1+xx)+(1+xx)*bv+2*xx*(1+xx)*rr+aa**3*bv**2*rr**2+2*aa**2*bv*rr+tt*bv*(1+aa+xx+aa**2*bv*rr*(2-rr)))
    mass=s.cancel((1+xx)*(1+aa*rr)*(hh+tt)-aa*(1+tt*bv)*ee*(1-hh))
    nact=s.cancel(((1+aa*rr)*ca*sc-rr*delt)/(2*rr*delt))
    mact=aa**2*(1+tt*bv)*ee*ca*sc/(2*mgap*delt)
    cact=ca*mass/(2*(1-hh)*delt); ahact=tt*ca*(1+aa*rr)*sc/(2*delt)
    dact=kk*(1+aa*rr)*(1-tt*xx*rr)*sc/(2*tt*gj*delt)
    alact=bv**2*(1+aa*rr)*mass/(2*(1-hh*bv)*delt)
    R=1+xx*hh; GG=rr+ss; MM=R-ww; QQ=1-xx*ss*hh; EJ=1-tt*xx*rr
    P4=ca*EJ+xx*(1-hh)*GG
    Pal=R*rr*(1+tt*xx*ss)+ss*MM*EJ
    pd=1+2*xx*rr+hh*xx**2*GG+xx*ss+tt*xx*(kk*ss+xx*rr*ee)
    cf=kk*(hh+tt)*EJ/(2*(1-hh)*pd); ahf=tt*P4/(2*pd)
    mf=kk*xx*(hh+tt)*rr*(1+tt*xx*ss)/(2*MM*ss*pd)
    alf=xx**2*ss*(hh+tt)*Pal/(2*QQ*pd)
    df=(kk*(1+tt*xx*ss)*(hh+tt)-(tt*ca-1+hh)*pd)/(2*gj*tt*pd)
    betaf=(kk*xx*(hh+tt)*GG-ee*pd)/(2*ee*pd)
    kap=s.cancel((1+2*nact)/(1+2*betaf)); theta=nact/bv
    Pkap=ca**2+ss*xx*(1+tt)+ss*xx**2*(rr+hh*ee+tt)
    PD=tt*xx**2*rr**2+xx*(tt-1+ee*(1-hh-2*tt))-1
    Palpha=1+ss-2*ee+xx*(rr**2+ss*(rr+hh*ee))
    expected=[0,-theta*ca/(1-hh),-theta*tt*ca*sc/(kk*(hh+tt)),theta*sc*PD/(kk*(hh+tt)*tt*gj),-nact*xx*Palpha/(kk*QQ),-2*theta*Pkap/(kk*(hh+tt)),-nact]
    actual=[mact-kap*mf,cact-kap*cf,ahact-kap*ahf,dact-kap*df,alact-kap*alf,1-kap,-kap*betaf+(1-kap)/2]
    for label,lhs,rhs in zip(['m','C','Ah','D','alpha','d','paired-beta'],actual,expected):check('TL3 exact native comparison '+label,lhs-rhs)
    bD=gj/ca; etaD=(actual[3]+actual[5]/2)*bD
    check('JD16 exact paired D normalized correction',etaD-theta*(sc*PD/tt-gj*Pkap)/(kk*(hh+tt)*ca))
    MB,qB,db=s.symbols('MB qB bD')
    check('JD15 whole D repair',(MB-db*qB)**2-(1-db)*MB**2-db*(MB**2-2*MB*qB+db*qB**2))
    om,pp,gamma,NN=s.symbols('omega p gamma N')
    du,dv=1+u*u,1+v*v; qA=v*om; XN=NN+pp*du*qA-u
    check('TL5 entire extra n row',XN**2-(XN+gamma*qA)**2-(2*gamma*u*qA-(2*gamma*pp*du+gamma**2)*qA**2-2*gamma*qA*NN))
    # UW2 is a polynomial divisibility proof before any residue is taken.
    zr=s.symbols('zratio'); xu=ss*ww*zr; bu=ss**2*ww*zr; ru=1-ss*ww
    cu=1+xu*ru; su=1+bu*(1-hh)+xu*ru; gu=1-ww+hh*xu; Ju=tt*cu-(1-hh)
    Deu=s.cancel((ss*delt).subs(bb,ww*zr))
    Mu=(1+xu)*(ss+ru)*(hh+tt)-ss*ww*(1+tt*bu)*(1-hh)
    NC=cu*Mu; DC=2*(1-hh)*Deu; NA=tt*cu*(ss+ru)*su; DA=2*Deu
    ND=(1+xu)*(ss+ru)*(1-tt*xu*ru)*su; DD=2*Ju*Deu
    NM1=(1+tt*bu)*cu*su; DM=2*gu*Deu
    Puw=(hh*DA+2*NA)*(DD+2*ND)*DC-2*NC*(4*NA*ND-hh*DA*DD)
    Quw=2*(4*NA*ND-hh*DA*DD)*DC
    check('UW3 exact s divisibility',Puw.subs(ss,0))
    check('UW3 exact w divisibility',Puw.subs(ww,0))
    check('UW3 denominator at joint corner',(Quw*NM1).subs({ss:0,hh:0})-16*tt)

if '--ba' in sys.argv or '--all' in sys.argv:
    tt,aa,hh,rr,mw,nw,cw,aw,dw,alw=s.symbols('t a h r m n C A D alpha')
    WW,XX,YY,ZZ,EE,ET,XXI=s.symbols('W X Y Z E eta xi')
    du,dv=1+u*u,1+v*v; ee=1-rr
    qq=nw*rr**2+(mw+half)*ee**2+alw; muw=mw-aa**2*qq
    pp=mw*(v*WW-XX+u*v*YY-u*ZZ)**2+nw*du*dv*XX**2+cw*dv*(u*WW-YY)**2+aw*du*dv*YY**2+(dw+half)*(u*WW+u*v*XX-YY-v*ZZ)**2+(alw+half)*du*(v*YY-ZZ)**2+half*dv*(XX+u*ZZ)**2-half*du*dv*ZZ**2
    Ystar=hh*EE/tt+ET
    explicit=cw*dv*((1-hh)*EE/tt-ET)**2+aw*du*dv*Ystar**2+(dw+half)*((1-hh)*EE/tt-ET+u*v*XXI)**2+(alw+half)*du*v*v*Ystar**2+(nw*(1+v*v+u*u*v*v)+dv/2)*XXI**2+2*mw*(v*v*EE/tt-u*v*XXI)*Ystar+mw*u*u*v*v*Ystar**2
    hard=pp.subs({ZZ:0,WW:EE/(tt*u),XX:XXI,YY:Ystar})-qq*aa**2*(v*EE/(tt*u)-XXI)**2-nw*u*u*XXI**2-muw*(v*EE/(tt*u)-XXI)**2
    check('new BA6 explicit complete H0',hard-explicit)
    ppv,qqv=s.symbols('p q')
    dp=1+ppv**2; K=ppv*XX+ZZ; VX=v*XX; JI=WW+dp*VX-ppv*YY-ppv*v*K
    rhs=dp*dv*nw*XX**2+mw*(ppv*v*WW+v*YY-K)**2+cw*dv*(WW-ppv*YY)**2+dp*dv*aw*YY**2+dw*JI**2+dp*alw*(ppv*XX+v*YY-K)**2+half*(JI**2+(1-ppv**2*v*v)*K**2+dp*(v*v*YY**2-2*v*YY*K+2*ppv*YY*VX+2*ppv*v*K*VX-ppv**2*VX**2))
    check('new BA11 complete six-root remote identity',ppv**2*pp.subs(u,1/ppv)-rhs)
    PI=ppv*ZZ; QI=qqv*ZZ
    rhs2=nw*(1+ppv**2)*(1+qqv**2)*XX**2+mw*(ppv*WW-ppv*qqv*XX+YY-QI)**2+cw*(1+qqv**2)*(WW-ppv*YY)**2+aw*(1+ppv**2)*(1+qqv**2)*YY**2+dw*(qqv*WW+XX-ppv*qqv*YY-PI)**2+alw*(1+ppv**2)*(YY-QI)**2+half*((qqv*WW+XX-ppv*qqv*YY-PI)**2+(1+ppv**2)*(YY-QI)**2+(1+qqv**2)*(ppv**2*XX**2+2*XX*PI-PI**2))
    check('new BA15 complete double-infinity identity',ppv**2*qqv**2*pp.subs({u:1/ppv,v:1/qqv})-rhs2)
    p3=pp.subs({YY:0,mw:0,nw:0,alw:0,cw:tt/2,dw:1/(2*(tt-1))})
    P=s.hessian(p3,(WW,XX,ZZ))/2
    det=P[0,0]*P[1,1]*P[2,2]+2*P[0,1]*P[0,2]*P[1,2]-P[0,0]*P[1,2]**2-P[1,1]*P[0,2]**2-P[2,2]*P[0,1]**2
    check('new BA14 finite-unit strict determinant',det-du*dv*u*u*v*v*(tt*(1/(tt-1)+1)*u*u/8+tt*v*v/(8*(tt-1))))
    pu=nw*dv*XX**2+mw*(v*WW-XX)**2+(cw+aw+half)*dv*YY**2+dw*(YY+v*ZZ)**2+alw*(v*YY-ZZ)**2+dv*XX**2/2
    pv=nw*du*XX**2+(mw+half)*(XX+u*ZZ)**2+(cw+dw+half)*(u*WW-YY)**2+aw*du*YY**2+alw*du*ZZ**2
    check('new BA literal u0 complete parent exit',pp.subs(u,0)-pu)
    check('new BA literal v0 complete parent exit',pp.subs(v,0)-pv)

if '--log-json' in sys.argv:
    import json
    outfile=sys.argv[sys.argv.index('--log-json')+1]
    json.dump({'checked_identities':len(audit_passes),'checks':audit_passes,'engine':'SymPy 1.14.0 from trusted PyPI installation /tmp/q222-refinement-trusted-deps','source_formula_policy':'Independently transcribed TeX formulas; no archived checker or Python dependency imported'},open(outfile,'w'),indent=2)
print('ALL REQUESTED DIAGNOSTIC CHECKS PASSED',len(audit_passes),flush=True)
