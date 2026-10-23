import sympy as S
import pickle,json
from pathlib import Path

OUT=Path(__file__).resolve().parents[2]/'data'
u,v,W,N,xi,ze,Rh,R,T=S.symbols('u v W N xi ze Rh R T')
z,k,H,B,r,cp,dp,g,q,ss,ts,js,ls,als,al,nu,mu,aa=S.symbols('z k H B r cp dp g q ss ts js ls als al nu mu aa')
U,V=S.symbols('U V')
y=k*k
eta=k*(N+q*Rh)
X=xi-z*r*v/B*W
Y=eta+H*u/B*W
Z=ze+u*v/B*W
K=u*k*W
Qn=u*xi-ze
J=(u*W-Y)+v*(u*X-Z)
M=v*W-X+u*v*Y-u*Z
du=1+u*u; dv=1+v*v
G=4*ls*N*N+4*g*Rh*Rh+8*(dp-cp)/B*K*N-8*g*ss*K*Rh
G+=4*(cp*ts+dp*js)/B**2*K*K+4*cp*ts/B**2*v*v*K*K-8*cp/B*v*v*K*(N+q*Rh)
G+=4*nu*du*dv*(z*z*r*v*v/B**2*W*W-2*z*v/B*W*xi)
G+=4*als*du*(v*eta-ze)**2+4*al*du*((1-H)*u*u*v*v/B**2*W*W-2*u*v/B*W*(v*eta-ze))
G+=4*z*z*mu*M*M+4*y*aa*du*dv*Y*Y
G+=2*(J*J+du*(v*Y-Z)**2+dv*(X+u*Z)**2-du*dv*Z*Z)
G=S.Poly(S.expand(G),N,xi,ze,Rh,W)
states=[N,xi,ze,Rh]
def cf(x,y):
    return G.coeff_monomial(x*y)
Q0=S.zeros(4)
for i,x in enumerate(states):
    for j in range(i,4):
        yy=states[j]
        Q0[i,j]=Q0[j,i]=cf(x,yy)/(1 if i==j else 2)
rho=S.cancel(cf(Rh,W)/(u*k))
pxi=S.cancel(cf(xi,W)/v+u*u*rho)
pze=S.cancel(cf(ze,W)/(u*v)-rho)
pn=S.cancel(cf(N,W)/(u*k))
phi=S.cancel(G.coeff_monomial(W*W)-ss*u*u*(v*v-y)*rho)
def uv(p):
    p=S.Poly(S.expand(p),u,v)
    assert all(i%2==j%2==0 for i,j in p.monoms())
    return S.Add(*(a*U**(i//2)*V**(j//2) for (i,j),a in p.terms()))
pxi,pze,pn,phi=map(uv,[pxi,pze,pn,phi])
def divvar(p,x,a):
    quot,rem=S.div(p,x-a,x)
    return S.factor(quot),S.factor(rem)
Axi,pxiz=divvar(pxi,U,z)
Bxi,rx=divvar(pxiz/z,V,y)
Aze,pzez=divvar(pze,U,z)
Bze,rz=divvar(pzez,V,y)
An,pnz=divvar(pn,U,z)
Bn,rn=divvar(pnz,V,y)
ds=4/B**2*(z*z*mu*(B+z*r)**2+(nu*r+S.Rational(1,2)*r*r)*z*z)
es=4/B**2*(nu*r+S.Rational(1,2)*r*r)*z*z
phis=S.Poly(S.expand(phi-(ds+es*V)/z**2*V*(U-z)**2),U)
assert S.factor(phis.coeff_monomial(1))==0
As,ra=divvar(phis.coeff_monomial(U*U),V,y)
Bs,rb=divvar(phis.coeff_monomial(U),V,y)
Cs,rc=divvar(Bs+z*As,V,y)
def tv(p):return p.subs({U:u*u,V:v*v})
Qn=u*xi-ze
f=S.Matrix([N,xi,ze,Rh,R,T])
out=(S.Matrix(states).T*Q0*S.Matrix(states))[0]
out+=(tv(Axi)-(v*v-y)*tv(Bxi))*R*xi+tv(Bxi)/ss*(u*v*v*Qn-u*v*k*Rh)*xi
out+=u*tv(Aze)*R*ze+tv(Bze)/ss*(v*v*Qn-v*k*Rh)*ze
out+=tv(An)*T*N+tv(Bn)/ss*(v*k*Qn-y*Rh)*N
out+=(ds+es*v*v)/z**2*R*R+tv(As)/ss*(-Rh*T+u*R*Qn)+tv(Cs)/ss**2*(v*v*Qn**2-2*v*k*Rh*Qn+y*Rh*Rh)
out=S.Poly(S.expand(out),*list(f))
Q=S.zeros(6)
for i,x in enumerate(f):
    for j in range(i,6):
        yy=f[j]
        Q[i,j]=Q[j,i]=S.factor(out.coeff_monomial(x*yy)/(1 if i==j else 2))
dat={'Q':Q,'Q0':Q0,'rho':rho,'Axi':Axi,'Bxi':Bxi,'Aze':Aze,'Bze':Bze,'An':An,'Bn':Bn,'As':As,'Bs':Bs,'Cs':Cs,'contact_residuals':[rx,rz,rn,ra,rb,rc],'phi':phi,'normal_polynomials':[pxi,pze,pn], 'symbols':list(locals().get('states',[]))}
with (OUT/'lII_generic.pkl').open('wb') as fp:pickle.dump(dat,fp)
summary={name:str(S.factor(dat[name])) for name in ['rho','Axi','Bxi','Aze','Bze','An','Bn','As','Bs','Cs']}
summary['Gram']=[[str(x) for x in Q.row(i)] for i in range(6)]
(OUT/'lII_generic.json').write_text(json.dumps(summary,indent=2))
print('Saved generic exact Gram and quotients.')
