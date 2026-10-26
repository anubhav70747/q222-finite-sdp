import sympy as S
s,b,h,t,u=S.symbols('s b h t u');l=1-h
C=(h+t)/(2*l)-(t+1)/(2*l)*s+(t+1)/(2*l)*s**2
Ah=t*(1-s+s**2)/2
bd=(1-s+s**2-(t+1)*s*b)/(2*t)
n1=(s-s**2-(h+t)*s*b)/2
c=l*(1-s*b)/t;q=1-c;E=t*(1+l*s*b);H=S.Rational(1,2);A=Ah/h
D0=lambda X:S.expand(X)
def trunc(ex):
  ex=S.cancel(ex)
  num,den=S.fraction(ex)
  # computed directly taylor total s,b<=2 with successive b and s derivatives
  r0=S.cancel(ex.subs({s:0,b:0}))
  ds=S.cancel(S.diff(ex,s).subs({s:0,b:0}));db=S.cancel(S.diff(ex,b).subs({s:0,b:0}))
  dss=S.cancel(S.diff(ex,s,2).subs({s:0,b:0})/2)
  dsb=S.cancel(S.diff(ex,s,b).subs({s:0,b:0}))
  dbb=S.cancel(S.diff(ex,b,2).subs({s:0,b:0})/2)
  return S.factor(r0+ds*s+db*b+dss*s*s+dsb*s*b+dbb*b*b)
B0=C*c*c+bd*q+q*q/2+A*c*c+2*H*c-S.Rational(1,2)
BC0=(-2*A*c-2*H+2*n1)/E
CC0=A/E**2
B_u2=2*C*c*c+2*(bd*q+q*q/2)+A*c*c+2*H*c-S.Rational(1,2)
BC_u2=BC0;CC_u2=CC0
Lam=(1-s+s*s-(2*l+h*(h+t))*s*b)/(2*h*t);kap=(h+t)/(2*t)
for name,ex,model in [('RI_BB0',B0,Lam),('RI_BC0',BC0,-2*Lam-4*kap*s*b),('RI_CC0',CC0,Lam+kap*s*b),('RI_BB_u2',B_u2,1/(2*h*t)+kap),('RI_BC_u2',BC_u2,-1/(h*t)),('RI_CC_u2',CC_u2,1/(2*h*t))]:
  z=S.factor(trunc(ex-model))
  if "u2" in name: z=S.factor(z.subs({s:0,b:0}))
  print(name, z)
  assert z==0
# RU h=0, Ah remains finite; source jets
rC=C.subs(h,0);rAh=Ah;rbd=bd.subs(h,0);rn1=n1.subs(h,0);rc=c.subs(h,0);rq=1-rc;rE=E.subs(h,0)
R_BB0=rC*rc**2+rbd*rq+rq**2/2+2*H*rc-S.Rational(1,2)
R_BC0=(-2*rAh*rc-2*H+2*rn1)/rE
R_CC0=2*rAh/rE**2
R_BB_u2=2*rC*rc**2+2*(rbd*rq+rq**2/2)-2*rAh*rc**2+2*H*rc-S.Rational(1,2)
R_BC_u2=(-2*H+2*rn1)/rE
R_CC_u2=2*rAh/rE**2
RLam=1/(t*(1+s)); RLamt= (1-s+s*s)/t
for name,ex,model in [('RU_BB0',R_BB0,RLamt-(1/t+S.Rational(1,2))*s*b),('RU_BC0',R_BC0,-2*RLamt+(3/t-1)*s*b),('RU_CC0',R_CC0,RLamt-2*s*b/t),('RU_BB_u2',R_BB_u2,S.Rational(1,2)),('RU_BC_u2',R_BC_u2,-1/t),('RU_CC_u2',R_CC_u2,1/t)]:
  z=S.factor(trunc(ex-model))
  if "u2" in name: z=S.factor(z.subs({s:0,b:0}))
  print(name,z)
  assert z==0
print('12 exact jet identities verified')
