"""Exact scaled-source formulas showing simultaneous endpoint residues."""
import sympy as s
z,tau,rho,H,sh=s.symbols('z tau rho H sh')
B=1/tau; y=1/(tau*rho); r=sh*rho; R=z*r; L=B-H
S0=1+R+L;G=B+R-z*(1-H)
delta=(1+z)*B**2+2*(1+z)*B*R+(1+z)*(1+y)*B+(z-y)*R**2+2*z*(1+y)*R+z*(1+y)
NH=(1+z)*(B+R)*(H+y)-(1+y)*(z-R)*L
Dbar=s.cancel(delta*tau**2*rho/(1+rho))
Nbar=s.cancel(NH*tau**2*rho)
Sb=s.cancel(S0*tau);Gb=s.cancel(G*tau)
print('Dbar:',s.factor(Dbar)); print('Nbar:',s.factor(Nbar));print('Sb:',Sb,'Gb:',Gb)
T1=(1+R); T2=1+R*tau; b1=(1+tau*rho)/(1+rho)
Cbar=T1*Nbar/(2*(1-H*tau)*(1+rho)*Dbar) # B*C
Mbar=b1*(1-r)*T1*Sb/(2*Gb*Dbar)    # B*m/z^2
Abar=T1*T2*Sb/(2*(1+rho)*Dbar)   # A*H/B
AlphaBar=T2*Nbar/(2*(1+rho)*Dbar) # alpha*(1-H)/B
nu=2*rho*T2*T1*Sb/((1+rho)*Dbar)
nr=nu/4-r/2
E=1+R-rho+H*tau*rho
D=rho**2*(1+z)*T2*(1-z*sh)*Sb/(2*E*(1+rho)*Dbar)
# define exact k/B using normalized source factors
m=z*z*tau*Mbar
koverB=(4*Abar*H+8*D*H*r*tau+4*AlphaBar*(1-H)+4*m*tau*(H+r-1)**2 +2*tau*(H+r-1)**2+4*nr*r*tau
 -4*Mbar+8*z*nr*r*tau+8*Abar*H*tau*rho)
print('koverB at rho=0:',flush=True)
print(s.factor(koverB.subs(rho,0)))
assert s.cancel(koverB.subs(rho,0))==0
# to compute derivative at origin use substitution first z=tau=0 after quotient derivative
print('derivative rho at z=tau=0:',s.factor(s.diff(koverB,rho).subs({z:0,tau:0,rho:0})))
print('rho quotient residue:',s.factor((s.diff(koverB,rho)).subs({z:0,tau:0,rho:0})))
assert s.simplify(s.diff(koverB,rho).subs({z:0,tau:0,rho:0})-2*sh)==0
print('PASS simultaneous endpoint rho=0 normalized k bound')
# The other singular shape: rho -> 1; e=1+zr-rho+rho H tau is the native j>0 opening.
ra,e=s.symbols('ra e')
GDb=s.factor(Dbar.subs(sh,ra/rho)); GNbar=s.factor(Nbar.subs(sh,ra/rho))
GSb=Sb.subs(sh,ra/rho)
GT1=1+z*ra;GT2=1+tau*z*ra
GAb=GT1*GT2*GSb/(2*(1+rho)*GDb)
GCb=GT1*GNbar/(2*(1-H*tau)*(1+rho)*GDb)
GMb=((1+tau*rho)/(1+rho))*(1-ra)*GT1*GSb/(2*Gb.subs(sh,ra/rho)*GDb)
GAlphab=GT2*GNbar/(2*(1+rho)*GDb)
Gnu=2*rho*GT2*GT1*GSb/((1+rho)*GDb)
Gnr=Gnu/4-ra/2
GDe=rho**2*(1+z)*GT2*(1-z*ra/rho)*GSb/(2*(1+rho)*GDb)
GkapB=(4*GAb*H+4*GAlphab*(1-H)+4*z*z*tau*tau*GMb*tau*(H+ra-1)**2
+2*tau*(H+ra-1)**2+4*Gnr*ra*tau-4*GMb+8*z*Gnr*ra*tau+8*GAb*H*tau*rho)
origin={rho:1,z:0,tau:0}
checks={
 'BC':(GCb,s.Rational(1,4)), 'Bm/z^2':(GMb,(1-ra)/4),
 'AH/B':(GAb,s.Rational(1,4)),'alpha(1-H)/B':(GAlphab,s.Rational(1,4)),
 'nu':(Gnu,s.Integer(1)),'nr':(Gnr,(1-2*ra)/4),
 'De':(GDe,s.Rational(1,4)),'(k-8DHr)/B':(GkapB,ra)}
for label,(expr,target) in checks.items():
    val=s.cancel(expr.subs(origin));assert s.cancel(val-target)==0,(label,val,target)
    print('rho=1:',label,'PASS')
print('PASS both singular source-shape residues')
