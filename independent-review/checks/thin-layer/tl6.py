exec(open('base.py').read())
Du=1+u**2; Dv=1+v**2
Rm=v*W-X+u*v*Y-u*Z; RD=u*W+u*v*X-Y-v*Z
F2=m2*Rm**2+n2*Du*Dv*X**2+C2*Dv*(u*W-Y)**2+A2*Du*Dv*Y**2+(D2+Rational(1,2))*RD**2+(al2+Rational(1,2))*Du*(v*Y-Z)**2+Rational(1,2)*Dv*(X+u*Z)**2-Rational(1,2)*Du*Dv*Z**2
F1=m1n*Rm**2+bt1n*Dv*(X+u*Z)**2+C1n*Dv*(u*W-Y)**2+A1n*Du*Dv*Y**2+Dt1n*RD**2+al1n*Du*(v*Y-Z)**2+Rational(1,2)*Du*((v*Y-Z)**2-Dv*Z**2)
# TL2 check
Dstar=s*Dl; rho=n2/e**2
Qs=1+B*(1+t+a*h+2*a*r+a*t)+B**2*(a*(h+t)+a**2*(h*(r+1)+t*(1+r-r**2)))
Ks=(1+s)*(h+t)*(1+x)
for sub in [{a:4,B:Rational(1,9),t:4,h:Rational(1,7),e:Rational(1,11)},{a:9,B:Rational(1,4),t:4,h:Rational(2,13),e:Rational(1,5)},{a:Rational(25,4),B:Rational(4,25)*Rational(9,16),t:Rational(16,9),h:Rational(1,17),e:Rational(2,23)}]:
    print('R10.2',simplify((2*r*Dstar*n2/s**2-(w*Qs-b*Ks)).subs(sub)))
    print('TL2 pi',simplify((pi_-(1-2*r*Dstar*rho*w/Pd)).subs(sub)), 'BPd-Pn?')
    xv=x.subs(sub); kk=sqrt(xv); e1=sqrt((t*b).subs(sub)); Dl1=1/e1
    etav=t.subs(sub); Lcv=h.subs(sub); t1v=t1.subs(sub); zv=zeta.subs(sub)
    cont={u:kk,v:Dl1*kk**2*etav,W:Dl1*kk,X:-kk*(1-t1v*zv),Y:Dl1*kk**2*Lcv,Z:1}
    f2=F2.subs(sub); f1=F1.subs(sub); kp=kappa.subs(sub)
    print(' contact', [simplify(f.subs(cont)) for f in (f2,f1)])
    print(' jet2',[simplify(diff(f2,q).subs(cont)) for q in (u,v,W,X,Y,Z)])
    print(' jet1',[simplify(diff(f1,q).subs(cont)) for q in (u,v,W,X,Y,Z)])
    # TL4/TL5
    frakw=(W+u*Y)/Du; At=v*frakw
    Atc=At.subs(cont); Kc=simplify(Atc/kk); print(' Xc+rk',simplify(cont[X]+r.subs(sub)*kk))
    gam=r.subs(sub)/Kc; nn=n2.subs(sub)
    p=t1v/R1.subs(sub)
    N=X-p*Du*At+u*Z
    Nsq=nn*Du*Dv*(X+gam*At)**2
    TL5=2*nn*gam*Du*Dv*u*At*Z-(2*nn*gam*p*Du+nn*gam**2)*Du*Dv*At**2-2*nn*gam*Du*Dv*At*N
    print(' TL5',simplify(expand(nn*Du*Dv*X**2-Nsq-TL5)))
    hv=h.subs(sub)
    dC=(C2-kappa*C1n).subs(sub); dD=(D2-kappa*(Dt1n-Rational(1,2))).subs(sub); dd=(1-kappa).subs(sub)
    dal=(al2-kappa*al1n).subs(sub); dAh=((A2-kappa*A1n)*h).subs(sub); dA=(A2-kappa*A1n).subs(sub)
    E=(nn*Du*Dv*X**2-Nsq)-nn*Dv*(X+u*Z)**2+dC*Dv*(u*W-Y)**2+(dD+dd/2)*RD**2+dal*Du*(v*Y-Z)**2+(dd/2)*Du*((v*Y-Z)**2-Dv*Z**2)+dAh*Du*Dv*(2*u*W*Y-hv*u**2*W**2)
    diffE=expand(f2-kp*f1-Nsq-dA*Du*Dv*(Y-hv*u*W)**2-E)
    print(' TL6 identity',simplify(diffE))
    print(' E contact',simplify(E.subs(cont)),[simplify(diff(E,q).subs(cont)) for q in (u,v,W,X,Y,Z)])
