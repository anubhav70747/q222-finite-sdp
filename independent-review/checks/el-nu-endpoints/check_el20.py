import sympy as s
a,c,x,q,kappa,nu,ell,S,r,p,mu,beta,alpha,D,c0,a0,la0,la1=s.symbols('a c x q kappa nu ell S r p mu beta alpha D c0 a0 la0 la1')
d=D+s.Rational(1,2)
L0=mu*(1-p)+(beta+s.Rational(1,2))*p*p
L1=(beta+s.Rational(1,2))*p*p
Bc=1-r+p*c*(1+r*a)
AP=(1+c)*(c0*(1-r)+a0*r*(1+a))+S*(d*Bc**2+(alpha+s.Rational(1,2))*(1+a)*r*r*c)
K=(1+a)*(d*Bc+(alpha+s.Rational(1,2))*r)
FZ=(1+a)*(alpha+D*c+d*a*c)
lam=(la0+la1*c)*(1+r*a)**2
Ahat=AP-kappa*lam
Khat=K-2*kappa*nu*lam
Fhat=FZ-lam*(nu**2*a*c+kappa*nu*ell*S)
J=Khat-nu*Ahat
H=s.expand(Fhat+nu**2*S*c*Ahat-2*nu*S*c*Khat)
print('H(a,0):',s.factor(H.subs(c,0)),flush=True)
print('H_aa(a,0)/2:',s.factor(s.diff(H,a,2).subs(c,0)/2),flush=True)
Hc3=s.expand(H).coeff(c,3)
print('H cubic:',s.factor(Hc3),flush=True)
R=(Hc3-nu**2*S**2*d*p**2*(1+r*a)**2)
print('difference cubic:',s.factor(R),flush=True)
assert s.expand(H.subs(c,0) - ((1+a)*alpha-S*ell*kappa*nu*la0*(1+r*a)**2)) == 0
assert s.expand(s.diff(H,a,2).subs(c,0)/2 + S*ell*kappa*la0*nu*r**2) == 0
assert s.expand(R) == 0
