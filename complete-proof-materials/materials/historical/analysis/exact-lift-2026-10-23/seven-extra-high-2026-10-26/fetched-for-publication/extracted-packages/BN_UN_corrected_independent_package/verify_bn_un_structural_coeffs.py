"""Mechanical SN/BN/UN coefficient cancellations before source substitution.
All equations are identities in polynomial coefficients (contact jets only establish divisibility)."""
import sympy as s
k,L,S,r,tau,m,be,D,alpha,c0,a0=s.symbols('k L S r tau m be D alpha c0 a0',nonzero=True)
x=k*S; zeta=L/k;e=k*tau/L; q=s.symbols('q',nonzero=True)
a,c=s.symbols('a c');d=D+s.Rational(1,2);th=s.symbols('th')

def chk(label,z):
  z=s.cancel(z);assert z==0,(label,s.factor(z));print('PASS',label)

def polycoef(f,a,c,i,j):return s.Poly(s.expand(f),a,c).coeff_monomial(a**i*c**j)
# SN dictionary before monic identities: no rational native source expressions substituted
AP=(1+c)*(c0*(1-r)+a0*r*(1+a))+S*(d*(1-r-e)**2*c+(alpha+s.Rational(1,2))*(1+a)*r**2)
AN=m*(1+r*a+e*c)**2+(be+s.Rational(1,2))*e**2*c*(1+c)
K=m*(1+r*a+e*c)+d*(1-r-e)+(alpha+s.Rational(1,2))*(1+a)*r+(be+s.Rational(1,2))*e*(1+c)
FZ=D+alpha*c+be*a+(m+alpha+be+s.Rational(1,2))*a*c
J1=K-L*AP;L1=J1-zeta*AN
H=FZ-L**2*S*c*AP+a*c*(zeta**2*AN-2*zeta*J1)
# polynomial quotients exactly via coefficients in a, avoiding huge expansions
def second_a_quot(f):
 P=s.Poly(s.expand(f),a)
 return s.expand(sum(coeff*(sum((power-1-i)*a**i*x**(power-2-i) for i in range(power-1))) for (power,),coeff in P.terms() if power>=2))
# Equivalent: (f(a)-f(x)-(a-x)f'(x))/(a-x)^2
La=s.cancel((L1-L1.subs(a,x))/(a-x))
Lc=s.cancel((L1.subs(a,x)-L1.subs({a:x,c:q}))/(c-q))
H20=s.cancel(second_a_quot(H)/c)
CR=s.cancel(H20+2*zeta*La)
chk('SN12 CR exact',CR+L*(m/k**2)*r*(a*L*r+2*c*k*tau+2*L))
# H02 is c3 coefficient of H at a=x divided by S
H02_from_c3=s.Poly(s.expand(H.subs(a,x)),c).coeff_monomial(c**3)/S
chk('SN12 H02 exact',H02_from_c3-k*tau**2*(be+m+s.Rational(1,2)))
# For BN/UN scalar after circuit, exact polynomial shifts by subtraction
# BN: b=beta/q, zn=k beta/q^2
b=be/q;zn=k*be/q**2
Ap=AP-zn*c/L**2;An=zeta**2*AN-b;Kp=K-2*be/(zeta*q);Jp=Kp-L*Ap;Lp=Jp-An/zeta
Hprime=FZ-be*a-be*x*c/q-L**2*S*c*Ap+a*c*(An-2*zeta*Jp)
# compare BNshift to original L1,H polynomial
chk('BN L prime - L',Lp-L1 - k*(be/q**2)/L*(c-q))
chk('BN H prime shift',Hprime-H + S*(c-q)**2*k*(be/q**2)+ 2*c*(a-x)*(c-q)*(be/q**2) - (a-x)*(c-q)*(be/q))
# above relationship is equivalent to Hprime = H minus S k ns delta² minus 2 ns c (a-x)delta
# UN: source polynomial with independent theta, h (theta=1-h)
h=s.symbols('h');theta=1-h;n=s.symbols('n',nonzero=True)
APu=(1+c)*(c0*(1-r)+a0*r*(1+a))+S*(d*(1-r)**2*c+(alpha+s.Rational(1,2))*(1+a)*r**2)
Ku=m*theta*(1+r*a)+d*(1-r)*(1+h*a)+(alpha+s.Rational(1,2))*(1+a)*r
FZu=D+alpha*c+n*a+(m*theta**2+alpha+(be+s.Rational(1,2))*theta**2)*a*c+d*h**2*a**2
A=m*zeta**2*(1+r*a)**2-n/q
z_n=k*n/q**2;APup=APu-z_n*c/L**2;Kup=Ku-2*n/(zeta*q);Jup=Kup-L*APup;Lup=Jup-A/zeta
Hup=FZu-n*a-n*x*c/q-L**2*S*c*APup+a*c*(A-2*zeta*Jup)
Lau=s.cancel((Lup-Lup.subs(a,x))/(a-x))
# c2 coeff of partial a H at a=x equals H11, if double jet holds
H11=s.Poly(s.diff(Hup,a).subs(a,x).expand(),c).coeff_monomial(c**2)
Lcu=s.Poly(s.expand(Lup.subs(a,x)),c).coeff_monomial(c) # linear c
P1=s.Poly(s.expand(APu),c).coeff_monomial(c)
chk('UN13 CG exact', H11/zeta- (2*L*P1.subs(a,x)+L*x*a0*r-2*z_n/L))
chk('UN13 CN exact', -H11-2*zeta*Lcu +S*a0*L**2*r)
H20u=s.cancel(second_a_quot(Hup)/c-d*h**2/c) # subtract dh² before dividing by c
CRu=s.cancel(H20u+2*zeta*Lau)
chk('UN13 CR exact',CRu+m*zeta**2*r*(2+r*a))

# BN12 dangerous C_N cancellation, computed only from c^j coefficients.
Ha=s.Poly(s.expand(s.diff(H,a).subs(a,x)),c)
h2=Ha.coeff_monomial(c**2);h3=Ha.coeff_monomial(c**3)
H11bn=h2+q*h3+h3*c
L0=s.Poly(s.expand(L1.subs(a,x)),c)
Lcbn=L0.coeff_monomial(c)+L0.coeff_monomial(c**2)*(c+q)
CNbn=-H11bn-2*L/k*Lcbn
CNbnexpected=-S*a0*L**2*r+(be+m+s.Rational(1,2))*tau**2*(c+1/L)+2*L*k*(m/k**2)*tau-m*tau**2
chk('BN12 source-safe CN cancellation',(CNbn-CNbnexpected).subs(q,(1-L)/L))


P1bn=s.Poly(s.expand(AP),c).coeff_monomial(c)
SN13=2*L**2*P1bn.subs(a,x)+L**2*x*a0*r+2*L**2*k*(m/k**2)*e*(1+2*r*x)+L**2*(be+s.Rational(1,2))*e**2/k-2*L*(m+be+s.Rational(1,2))*e+(c+q)*L**2*(m+be+s.Rational(1,2))*e**2/k
chk('SN13 full contact coquotient identity',(k*H11bn-SN13).subs(q,(1-L)/L))

print('All universal structural identities proved exactly')
