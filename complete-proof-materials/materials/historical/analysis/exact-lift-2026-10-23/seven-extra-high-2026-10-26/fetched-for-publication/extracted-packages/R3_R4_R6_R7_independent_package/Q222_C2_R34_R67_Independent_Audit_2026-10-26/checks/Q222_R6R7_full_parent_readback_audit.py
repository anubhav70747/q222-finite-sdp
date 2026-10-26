import sympy as s
u,v,o,N,MB=s.symbols('u v o N MB',real=True)
p,c,E,eps,j,mu1,mu2,n1,n2,Km,C,A,D,al,Ah,betaD,h,l=s.symbols('p c E eps j mu1 mu2 n1 n2 Km C A D al Ah betaD h l',real=True)
Du=1+u*u;Dv=1+v*v; q=1-c;H=al+s.Rational(1,2)
qA=Du*v*o
X=N+p*qA-u
qB=v*(Du-u*N-p*u*qA)
Tb=MB+c*qB
W=o+u*Tb/Du
Y=u*o-Tb/Du
J=MB-q*qB
F0=p*p/2-n2-2*n1*p+mu2
AN=p-2*mu1-2*n1
V0=1-p*u*v*o

# helpers
AB=eps*o*v;BB=v;CC=E*u*o;DD=eps**2*j**3*E
varrho=al*Du/DD**2
fA=Du**2*Dv*F0/eps**2

def cr(expr,name):
  expr=s.together(expr)
  out=s.cancel(expr)
  assert out==0,(name,out)
  print('PASS',name,flush=True)

# R6 source-excess remainder
GR6=Km*N**2-2*mu1*qA*N+mu2*qA**2-2*n1*Du*Dv*X*v*o-n2*Du*Dv*(v*o)**2+C*Dv*Tb**2+A*Du*Dv*Y**2+(D+s.Rational(1,2))*J**2+H*Du*(v*Y-1)**2+s.Rational(1,2)*Dv*(X+u)**2-s.Rational(1,2)*Du*Dv
L=A*Dv+H*v*v
fBB=(C*c*c*Du**2*Dv+(D+s.Rational(1,2))*q*q*Du**2+Du*L*c*c)*V0**2+2*H*Du*c*V0-Du/2-Du**2*mu2*v*v*o*o
fBC=(-2*Du*L*c*V0-2*H*Du+2*n1*Du*Dv)/E
fCC=(Du*L+Du*Dv*n2*v*v)/(E*E)
R_A=Du*AN/eps
R_C=2*c*u*v*L/E
R_B=Du*(p-2*n1)*v*v*o-2*u*v*(C*c*c*Du*Dv+A*c*c*Dv+(D+s.Rational(1,2))*q*q*Du+H*c*c*v*v)*V0-2*H*c*u*v
S_B=(2*C*c*Du*Dv-2*(D+s.Rational(1,2))*q*Du+2*A*c*Dv+2*H*c*v*v)*V0+2*H
S_C=-2*L/E
HN=Km+Dv/2+C*c*c*Dv*u*u*v*v+(A*c*c*Dv/Du)*u*u*v*v+(D+s.Rational(1,2))*q*q*u*u*v*v+(H*c*c/Du)*u*u*v**4
HNM=-C*c*Dv*u*v-(A*c*Dv/Du)*u*v+(D+s.Rational(1,2))*q*u*v-(H*c/Du)*u*v**3
HMM=C*Dv+A*Dv/Du+D+s.Rational(1,2)+H*v*v/Du
base0=GR6.subs({N:0,MB:0})
targ=fA*AB**2+fBB*BB**2+fBC*BB*CC+fCC*CC**2+varrho*DD**2
cr(base0-targ,'R6 scalar five packet identity')
cr(s.diff(GR6,N).subs({N:0,MB:0})-(R_A*AB+R_B*BB+R_C*CC),'R6 first N normal')
cr(s.diff(GR6,MB).subs({N:0,MB:0})-(S_B*BB+S_C*CC),'R6 first MB normal')
cr(s.diff(GR6,N,2)/2-HN,'R6 hard NN')
cr(s.diff(GR6,MB,N)/2-HNM,'R6 hard NM')
cr(s.diff(GR6,MB,2)/2-HMM,'R6 hard MM')
# R7 use original parent, subtract 4 positive squares and exact canceled A,D formulas, direct generic
# Need source-bound matching: D*q = betaD; Ah=A*h, q=1-c. Treat A,D split as independent bounded Ah,betaD.
# Reconstruct expression RU7 from direct selected polynomial identities with A=Ah/h,D=betaD/q.
# Avoid h inverse in target: identity A(Y^2-(Y-h*u*W)^2)=Ah(2uWY-hu^2W^2).
# D(J^2-MB^2)=betaD(-2MB*qB+q*qB^2).
GR7=-2*mu1*qA*N+mu2*qA**2-2*n1*Du*Dv*X*v*o-n2*Du*Dv*(v*o)**2+C*Dv*Tb**2+2*Ah*Du*Dv*u*W*Y-Ah*h*Du*Dv*(u*W)**2-2*betaD*MB*qB+betaD*q*qB**2+s.Rational(1,2)*J**2+H*Du*(v*Y-1)**2+s.Rational(1,2)*Dv*(X+u)**2-s.Rational(1,2)*Du*Dv
fBB7=(C*c*c*Du**2*Dv+(betaD*q+q*q/2)*Du**2+H*Du*v*v*c*c-Ah*Du*Dv*(2*u*u+h*u**4)*c*c)*V0**2+2*H*Du*c*V0-Du/2-Du**2*mu2*v*v*o*o
fBC7=(2*Ah*Du*Dv*(l*u*u-1)*c*V0-2*H*Du*v*v*c*V0-2*H*Du+2*n1*Du*Dv)/E
fCC7=(Ah*(2-h)*Du*Dv+H*Du*v*v+Du*Dv*n2*v*v)/(E*E)
RB7=Du*(p-2*n1)*v*v*o-2*u*v*(C*c*c*Du*Dv+(betaD*q+q*q/2)*Du+H*c*c*v*v-Ah*c*c*Dv*(2*u*u+h*u**4))*V0-2*H*c*u*v
RC7=2*c*u*v*(Ah*Dv*(1-l*u*u)+H*v*v)/E
SB7=(2*C*c*Du*Dv-(2*betaD+q)*Du-2*Ah*c*Dv*(2*u*u+h*u**4)+2*H*c*v*v)*V0+2*H
SC7=2*(Ah*Dv*(l*u*u-1)-H*v*v)/E
HNN7=Dv/2+C*c*c*Dv*u*u*v*v-(Ah*c*c*Dv/Du)*(2*u*u+h*u**4)*u*u*v*v+(betaD*q+q*q/2)*u*u*v*v+(H*c*c/Du)*u*u*v**4
HNM7=-C*c*Dv*u*v+(Ah*c*Dv/Du)*(2*u*u+h*u**4)*u*v+(betaD+q/2)*u*v-(H*c/Du)*u*v**3
HMM7=C*Dv-(Ah*Dv/Du)*(2*u*u+h*u**4)+s.Rational(1,2)+H*v*v/Du
# Identity uses l=1-h, replace before compare
cr(GR7.subs({N:0,MB:0})-(fA*AB**2+fBB7*BB**2+fBC7*BB*CC+fCC7*CC**2+varrho*DD**2).subs(l,1-h),'R7 scalar five packet identity')
cr(s.diff(GR7,N).subs({N:0,MB:0})-(R_A*AB+RB7*BB+RC7*CC).subs(l,1-h),'R7 first N normal')
cr(s.diff(GR7,MB).subs({N:0,MB:0})-(SB7*BB+SC7*CC).subs(l,1-h),'R7 first MB normal')
cr(s.diff(GR7,N,2)/2-HNN7,'R7 hard NN')
cr(s.diff(GR7,MB,N)/2-HNM7,'R7 hard NM')
cr(s.diff(GR7,MB,2)/2-HMM7,'R7 hard MM')
