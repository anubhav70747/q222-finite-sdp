"""A single formal rational Puiseux specialization shows R8 is nonempty; not a sweep."""
import sympy as s
h=s.symbols('h',positive=True)
a=1/h;B=h*h/4;x=a*B;r=1-h/2;e=1-r;l=1-h;t=l/(1+x*r-h)
y=t*B;g=1-a*e+h*x;j=1+x*r-l/t
Q=1+B+y+x*(h+t+2*r)+x*x*(t*(1+e*r)+h*(1+r))+B*x*(h+t)
K=(a+1)*(h+t)*(1+x)
expr=s.factor(e*Q-B*K)
assert s.simplify(j-h)==0
assert s.limit(expr/h,h,0,dir='+')==s.Rational(1,4)
assert s.limit(g,h,0,dir='+')==s.Rational(1,2)
assert s.limit(t,h,0,dir='+')==1
print('Formal R8 exact opening j=h; g->1/2, t->1, (eQ-BK)/h->1/4; n>0 for small h.')
