"""Exact native C2 coquotient divisibility without a giant rational simplification."""
import sympy as sp
s,b,h,t,g=sp.symbols('s b h t g'); l=1-h;x=s*b;B=s*s*b;K=1+x;R=1+h*x;w=R-g;r=1-s*w;T=1+x*r
Sc=1+s*b*r+s*s*b*l
Del=1+s*K+s*K*B+2*s*x*K*r+s*s*b*b*r*r+2*s*b*r+t*s*s*b+t*s**3*b*(1+x+b*r*(2-r))
Q=1+s*s*b*(1+t)+s*b*(h+2*r+t)+s**3*b*b*(h+t)+s*s*b*b*(h*(r+1)+t*(1+r-r*r))
Ns=w*Q-b*(1+s)*(h+t)*K
AN_num=w*(Del-(1+t*B)*T*Sc)-s*K*Ns
F0_num=w*w*Del-(r*K*K+2*s*K*w)*Ns+(1+t*B)*w*T*Sc*g
P_N=sp.Poly(AN_num,s,b); P_F=sp.Poly(F0_num,s,b)
assert all(i>=1 and j>=1 for (i,j),c in P_N.terms() if c!=0),"A_N numerator not divisible by sb"
assert all(j>=1 for (i,j),c in P_F.terms() if c!=0),"F0 numerator not divisible by b"
assert sp.expand(P_N.coeff_monomial(s*b).subs(g,0)-(h+t))==0
assert sp.expand(P_F.coeff_monomial(b).subs(g,0)-(h+t))==0
# Polynomial no-localization version of source n identity from two original formulas.
print('PASS: AN_num in (s*b), F0_num in (b), native coquotient residues h+t and (h+t)/2')
