"""Reconstruct JF1 from the complete original quarter-parent, before source entry."""
import json
from pathlib import Path
import sympy as s
u,v,EE,QQ,xi=s.symbols('u v EE QQ xi',nonzero=True)
x,eta,zeta,r,h=s.symbols('x eta zeta r h',nonzero=True)
m,C,A,D,al,be=s.symbols('m C A D al be')
a=1/zeta; t=eta; y=t*x/a; U=u*u; V=v*v; B=be+s.Rational(1,2); Dp=D+s.Rational(1,2)
j=t*(1+x*h)-1+r
W=(v+EE)/(t*u); Y=QQ+r*v/t; X=xi-h*u
raw=m*(v*W-X+u*v*Y-u)**2+C*(1+V)*(u*W-Y)**2+A*(1+U)*(1+V)*Y**2
raw+=Dp*(u*W+u*v*X-Y-v)**2+(al+s.Rational(1,2))*(1+U)*(v*Y-1)**2+B*(1+V)*(X+u)**2-s.Rational(1,2)*(1+U)*(1+V)
P=EE/t-QQ; RD=v*(U-x); TD=RD-u*v*xi/h
q=al+(m+B)*(1-h)**2; mu=m-a*a*q; lam=mu/(t*a*x); HH=a*v*W-u-a*xi
F0=V-y; Fb=F0+v*EE; Ft=F0+(t/r)*v*QQ
T0=A+al+s.Rational(1,2); Tm=T0+m; e4=r*r*Tm/t**2
bb=A*r/t-m*(1-h)-al-s.Rational(1,2); bprime=bb+e4*y*t/(2*r)
kap=2*(B*(1-h)-Dp*j/t-m*r/t)
gam=2*bprime/h+2*Dp*j/r-2*t*B*(1-h)/r
LL=U*V+(a*x/t)*V*V/U+t*x*x/a-3*x*V
circuit=lam*LL+mu*V/(t*t*U)*EE**2+2*(mu*v**3/(t*t*U)-lam*x*v)*EE
rad=C/t**2*(Fb-r*Ft)**2+r*r*T0/t**2*Ft**2+(2*m*r/t**2+e4*x)*Fb*Ft
jf=q*HH**2+circuit+C*P**2+A*QQ**2+Dp*(P-h*TD)**2+2*bprime*QQ*TD
jf+=rad+e4*U*Ft**2-e4*y*u*(HH+a*xi)*Ft+A*U*QQ**2+B*(1+V)*xi**2-2*(mu/a)*HH*xi-mu*xi**2+kap*u*Ft*xi+gam*u*v*QQ*xi
pol=s.Poly(s.expand(raw-jf),EE,QQ,xi,u,v)
K=1+x; R=1+x*r; J=1+x*h-(1-r)/eta; MM=R-(1-h)/zeta
Es=1-x*eta*h; Qs=1-x*zeta*r; yy=x*eta*zeta
P4=(1+x*h)*Es+x*(1-r)*(h+zeta)
Pal=R*h*(1+yy)+zeta*MM*Es
Pd=1+2*x*h+r*x*x*(h+zeta)+x*zeta+eta*x*(K*zeta+x*h*(1-h))
src={C:K*(r+eta)*Es/(2*(1-r)*Pd),A:eta*P4/(2*r*Pd),m:K*x*(r+eta)*h*(1+yy)/(2*MM*zeta*Pd),al:x*x*zeta*(r+eta)*Pal/(2*Qs*Pd),D:(K*(1+yy)*(r+eta)-(eta*(1+x*h)-1+r)*Pd)/(2*J*eta*Pd),be:(K*x*(r+eta)*(h+zeta)-(1-h)*Pd)/(2*(1-h)*Pd)}
out={'residual_support':[list(mon) for mon,c in pol.terms()],'coefficient_checks':{}}
for mon,cf in pol.terms():
    vv=s.factor(cf.subs(src,simultaneous=True))
    out['coefficient_checks'][str(mon)]={'zero':vv==0}
    if vv!=0: out['coefficient_checks'][str(mon)]['residual']=str(vv)
out['all_zero']=all(v['zero'] for v in out['coefficient_checks'].values())
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
