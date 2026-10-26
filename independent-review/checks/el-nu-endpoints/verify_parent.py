import sympy as sp
u,w,W,Z,Xe,Ny=sp.symbols('u w W Z Xe Ny')
r,e,ell,zeta,S,t,gam=sp.symbols('r e ell zeta S t gam')
c0,a0,a1,m0,m1,m2,D,beta=sp.symbols('c0 a0 a1 m0 m1 m2 D beta')
Du=1+u*u;Dw=1+w*w; a=u*u;c=w*w
C=c0/(S*(1-r));A=a0/(S*r);alpha=S*a1;d=D+sp.Rational(1,2)
X=Xe-e*w*W;Y=r*u*W+t*gam*Ny
V=W+u*Y;Bc=w*(X+u*Z)
Mpart=m0*Bc**2-2*m1*Bc*V+m2*V**2
Pcap=Mpart+C*Dw*(u*W-Y)**2+A*Du*Dw*Y**2+d*(u*w*W+u*X-w*Y-Z)**2+(alpha+sp.Rational(1,2))*Du*(Y-w*Z)**2+(beta+sp.Rational(1,2))*Dw*(X+u*Z)**2-sp.Rational(1,2)*Du*Dw*Z**2
# apply gam^2 = r*(1-r), t^2=S, for algebraic checks

def norm(expr):
    expr=sp.expand(expr)
    expr=expr.subs(S,t**2).subs(gam**2,r*(1-r))
    expr=sp.cancel(expr)
    return sp.factor(expr)

def chk(tag,expr):
    x=norm(expr)
    if x==0: print(tag,'OK')
    else:print(tag,'FAIL:',str(x)[:1000])

AN=m2*(1+r*a)**2+2*m1*e*c*(1+r*a)+m0*e**2*c**2+(beta+sp.Rational(1,2))*e**2*c*(1+c)-d*e*(2*(1-r)-e)*a*c
AP=(1+c)*(c0*(1-r)+a0*r*(1+a))+S*(d*(1-r)**2*c+(alpha+sp.Rational(1,2))*(1+a)*r**2)
K=m1*(1+r*a)+m0*e*c+d*(1-r-e)+(alpha+sp.Rational(1,2))*(1+a)*r+(beta+sp.Rational(1,2))*e*(1+c)
FZ=D+alpha*c+beta*a+(m0+alpha+beta+sp.Rational(1,2))*a*c
J1=K-ell*AP;L1=J1-zeta*AN
H=FZ-ell**2*S*c*AP+a*c*(zeta**2*AN-2*zeta*J1)
N=W-zeta*u*w*Z;G=u*W/t-ell*t*w*Z
Graph=AN*N**2+AP*G**2-2*u*w*Z*L1*N+H*Z**2
F0=Pcap.subs({Xe:0,Ny:0})
chk('NU12 scalar graph',F0-Graph)

Ax=-m0*e*c-m1*(1+r*a)+d*a*(1-r-e)-e*(beta+sp.Rational(1,2))*(1+c)
Tx=zeta*c*Ax+beta-D+c*(m0+beta+sp.Rational(1,2))
YP=(1+c)*(-c0+a0*(1+a))+S*(m2*(1+r*a)+m1*e*c-d*c*(1-r-e)+(alpha+sp.Rational(1,2))*(1+a)*r)
YQ=-m1*a+d-(alpha+sp.Rational(1,2))*(1+a)
YC=ell*YP+YQ
Lx=w*Ax*N+u*Z*Tx
Ly=gam*(YP*G+YC*t*w*Z)
HN=m0*c+d*a+(beta+sp.Rational(1,2))*(1+c)
Hy=r*c0*(1+c)+(1-r)*a0*(1+a)*(1+c)+S*r*(1-r)*(m2*a+d*c+(alpha+sp.Rational(1,2))*(1+a))
Hxy=-u*w*t*gam*(m1+d)
px=sp.expand(Pcap).coeff(Xe,1).subs(Ny,0)/2
py=sp.expand(Pcap).coeff(Ny,1).subs(Xe,0)/2
chk('NU19 first X normal',px-Lx)
chk('NU19 first Y normal',py-Ly)
chk('NU19 hard XX',sp.expand(Pcap).coeff(Xe,2)-HN)
chk('NU19 hard YY',sp.expand(Pcap).coeff(Ny,2)-Hy)
chk('NU19 hard XY',sp.expand(Pcap).coeff(Xe*Ny)-2*Hxy)
print('done')
