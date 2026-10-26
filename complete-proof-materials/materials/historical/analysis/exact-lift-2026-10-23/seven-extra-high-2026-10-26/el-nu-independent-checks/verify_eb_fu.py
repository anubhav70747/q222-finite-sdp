import sympy as sp
u,v,W,X,Y,Z=sp.symbols('u v W X Y Z');rr,S,kappa,lam,aa,bb=sp.symbols('rr S kappa ell aa bb'); c,s=sp.symbols('c s');t=sp.symbols('t',nonzero=True)
C,A,m,D,alpha,beta=sp.symbols('C A m D alpha beta');Du=1+u**2;Au=1+u*u*s*s
# simplified EB core identity EB5, strip pure-field square roots via factor algebra
Ru,Ra=sp.symbols('Ru Ra',nonzero=True)
pp=u*W/t;yy=Y/t;qq=t*Z
Nb=Ra*(X+u*c*c/Au*Z); Nr=Ru*u*s*c/Ra*Z
Kh=(c*(u*W-Y)+s*(u*X-Z))**2+Du*(s*Y-c*Z)**2+(X+u*Z)**2-Du*Z**2
Kexp=Nb**2+Nr**2-2*Du*s*c*pp*qq +t**2*c**2*(pp-yy)**2+t**2*Du*s**2*yy**2+2*t*u*s*c/Ra*(pp-yy)*Nb +2*t*u*s**2*Ru/Ra*(pp-yy)*Nr
# remove Ru etc with substitution of square relations; norm as polynomial in c,s, W X Y Z, reduce c2+s2=1

def norm_eb(z):
 z=sp.cancel(sp.together(z))
 num,den=sp.fraction(z)
 num=sp.expand(num).subs(Ru**2,Du).subs(Ra**2,Au)
 num=sp.cancel(num)
 num=sp.Poly(num,c).rem(sp.Poly(c**2+s**2-1,c)).as_expr()
 return sp.factor(num)
print('EB5 identity:',norm_eb(Kh-Kexp)==0,flush=True)
# EB6 normalized original event roots with homogeneous coordinates
MB=sp.sqrt(m)*(s*(W+u*Y)-c*(X+u*Z))
# equality checked without radicals by squared amplitude simplification with substitutions m->1, D->1, alpha->1
# m RHS Nm -sqrt m c/Ra Nb -sqrt(m Du) s/Ra Nr +sqrt(m S)us y
Nm=s*W
m_r=Nm-c/Ra*Nb-Ru*s/Ra*Nr+t*u*s*yy
print('EB6 m amplitude:',norm_eb(s*(W+u*Y)-c*(X+u*Z)-m_r)==0,flush=True)
d_r=t*c*(pp-yy)+u*s/Ra*Nb-Du*s/(Au*t)*qq
print('EB6 D amplitude:',norm_eb(c*(u*W-Y)+s*(u*X-Z)-d_r)==0,flush=True)
a_r=Du*(s*Y-c*Z)
print('EB6 alpha amplitude:',sp.expand(s*Y-c*Z-(t*s*yy-c*qq/t))==0,flush=True)
# NU25 factor exact. Compute FU15 expression at nu=0 and 0<kappa<1
k,y,z=sp.symbols('k y z',positive=True)
E=1/(1+k);m0=k*k/(2*(1-k*k));q0=sp.Rational(1,2);a1=E/2
p21=E*(m0+q0);p11=-m0-q0+2*E*q0*a1;p10=2*(m0+q0)*E*a1;p01=4*m0*q0*a1
P=p21*y*y*z+p11*y*z+p10*y+p01*z
Pgood=E*(m0+sp.Rational(1,2))*z*(y-k)**2+(m0+sp.Rational(1,2))*E*E*y
print('NU25 FU factor:',sp.factor(P-Pgood)==0,flush=True)
print('NU26 m limit:',sp.limit(sp.expand(P/m0).subs(y,y).subs(z,z),k,1).factor(),flush=True)
# FU full normal verification, symbolically no native source specialization
B0,Ny=sp.symbols('B0 Ny');gam=sp.symbols('gam');mc,betac,c0,a0,a1,ds=sp.symbols('mc betac c0 a0 a1 ds')
Yf=rr*u*W+t*gam*Ny;Xf=B0-u*Z
cf=c0/(t*t*(1-rr));af=a0/(t*t*rr);alf=t*t*a1;Df=t*t*ds
Pfull=mc*(v*W-Xf+u*v*Yf-u*Z)**2+cf*(1+v*v)*(u*W-Yf)**2+af*Du*(1+v*v)*Yf**2+(Df+sp.Rational(1,2))*(u*W+u*v*Xf-Yf-v*Z)**2+(alf+sp.Rational(1,2))*Du*(v*Yf-Z)**2+(betac+sp.Rational(1,2))*(1+v*v)*(Xf+u*Z)**2-sp.Rational(1,2)*Du*(1+v*v)*Z**2
M0=(v/kappa)*(1+rr*u*u)*(kappa*W)-B0
AP=(1+v*v)*(c0*(1-rr)+a0*rr*(1+u*u)) +t*t*((Df+sp.Rational(1,2))*(1-rr)**2+(alf+sp.Rational(1,2))*(1+u*u)*rr**2*v*v)
P=u*W/t;Q=t*Z;T=kappa*W;U=u*Z
F0=mc*M0**2+AP*P**2+Du*(ds*Du*v*v+a1)*Q**2+sp.Rational(1,2)*Du*v*v*U**2-2*(Df+sp.Rational(1,2))*Du*v*v*B0*U-Du*v*(1+2*Df*(1-rr)+2*alf*rr)*P*Q+((betac+sp.Rational(1,2))*(1+v*v)+(Df+sp.Rational(1,2))*u*u*v*v)*B0**2+2*(Df+sp.Rational(1,2))*(1-rr)*u*v*t*B0*P
Ly=mc*u*v*t*gam*M0+gam*(1+v*v)*(-c0+a0*Du)*P+t*t*gam*(-(Df+sp.Rational(1,2))*(1-rr)+(alf+sp.Rational(1,2))*Du*rr*v*v)*P-(Df+sp.Rational(1,2))*u*v*t*gam*B0+gam*v*Du*(Df-alf)*Q
Hy=rr*c0*(1+v*v)+(1-rr)*a0*Du*(1+v*v)+t*t*gam**2*(mc*u*u*v*v+Df+sp.Rational(1,2)+(alf+sp.Rational(1,2))*Du*v*v)

def normal(e):
 z=sp.cancel(sp.together(e.subs(gam**2,rr*(1-rr))))
 z=sp.cancel(sp.together(z))
 return sp.factor(z)
print('FU8 all scalar:',normal(Pfull.subs(Ny,0)-F0)==0,flush=True)
print('FU9 full first normal:',normal(sp.expand(Pfull).coeff(Ny,1)/2-Ly)==0,flush=True)
print('FU9 full hard:',normal(sp.expand(Pfull).coeff(Ny,2)-Hy)==0,flush=True)
