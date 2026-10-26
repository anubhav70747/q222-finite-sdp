import sympy as S
x,h,t=S.symbols('x h t',positive=True)
G=(1-h)/t
K=1+x;T=K-x*(1-h);M=K-t
Pd=S.factor(1+x*K*G+x*(x+2)*h)
be=S.factor(x*(h+G)/(G*t)-Pd/K)
# comparison: lower bound coefficient > gamma upper leading coefficient condition
E=(h+G)
Diff=S.factor(x*T*(E*M+h*t)-4*K*h*E)
print('beta/(2a f^2)=',S.factor(be))
print('violating difference=',Diff)
print('numerator factor=',S.factor(S.together(Diff).as_numer_denom()[0]))
print('k=1 =',S.factor(Diff.subs(x,1)))
print('k=1 beta =',S.factor(be.subs(x,1)))
print('k=1 t=1 diff=',S.factor(Diff.subs({x:1,t:1})))
print('k=1 t=1 beta=',S.factor(be.subs({x:1,t:1})))
