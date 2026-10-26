import sympy as s,sys
sys.path.insert(0,'/mnt/data/Q222_C1_fixed_admission_2026_10_26')
from c1_native import NATIVE,CONTACT,k,f,g,a,t
vals={k:4,f:s.Rational(1,1000),g:s.Rational(1,8),a:1,t:1}
w=[s.factor(t.subs(vals)) for t in NATIVE]; c=[s.factor(t.subs(vals)) for t in CONTACT]
d,la,n,m,C,A,D,al,be=w
# Event reflection: lambda<->n, m<->C, A<->alpha, D<->beta; CHSH invariant
wr=[d,0,0,C,m,al,be,A,D]
uc,vc,wc,xc,yc,zc=c
S0=wc;Y0=yc/uc;theta=t*g*g
th=s.factor(theta.subs(vals)); u1=1/uc;v1=vc
w1=s.factor(vc*S0-uc*(th-1)); x1=s.factor(-S0-uc*vc*(th-1)); y1=s.factor(1-uc*vc*Y0); z1=s.factor(uc*Y0+vc)
cp=[s.factor(z/z1) for z in (w1,x1,y1,z1)]
# same feature diagonal scaling iff contact normalized Z=1
import math
B=2*wr[4]+2*wr[6]+wr[0]
Sb=s.factor(4*wr[5]*u1**2*cp[2]**2)
Sc=s.factor((4*wr[8]+2*wr[0])*v1**2*cp[1]**2)
Se=s.factor(4*wr[8]*u1**2*v1**2)
Ub=s.factor(64*wr[3]*wr[5]*(2*wr[8]+wr[0])*(2*wr[4]+2*wr[6]+wr[0])*u1**2*v1**2*cp[0]**4/((2*wr[8]+wr[0]+2*wr[3])*(2*wr[5]+2*wr[4]+2*wr[6]+wr[0])))
print('reflected source weights:',wr)
print('contact u,v psi:',u1,v1,cp)
print('signs positive:',all(x>0 for x in wr if x!=0),'W>0 X<0 Y>0',cp[0]>0,cp[1]<0,cp[2]>0)
print('necessary squared triangle:',[float(y) for y in (Sb,Sc,Se)],'upper gamma^2',float(Ub))
print('sqrt bounds approximate',float(s.sqrt(Sb)-s.sqrt(Sc)-s.sqrt(Se)),float(s.sqrt(s.sqrt(Ub))))
# Test formula for reflected contact using full original parent
from c1_small_map import sym_parent,u,v,W,X,Y,Z
expr=s.factor(sym_parent(wr).subs({u:u1,v:v1,W:cp[0],X:cp[1],Y:cp[2],Z:1}));print('F(contact) exact',expr)
