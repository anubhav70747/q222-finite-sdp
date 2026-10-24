import sympy as sp
half=sp.Rational(1,2)
u,v,W,X,Y,Z,q,p=sp.symbols('u v W X Y Z q p',real=True)
K=sp.symbols('K_m',nonnegative=True)
h,t=sp.symbols('h t',positive=True)
l=1-h
Cs=(h+t)/(2*l); As=t/(2*h); Ds=1/(2*(t-l))
# R6.6 residue (q=0, n=alpha=0) written directly from the chart:
Q6=(K*(W+u*Y)**2 + Cs*(u*W-Y)**2 + As*(1+u**2)*Y**2 + (Ds+half)*(u*X-Z)**2
    + half*(1+u**2)*Y**2 + half*(X+u*Z)**2 - half*(1+u**2)*Z**2)
st=[W,Y,X,Z]
G=sp.Matrix(4,4,lambda i,j: sp.expand(Q6).coeff(st[i],2) if i==j else half*sp.expand(Q6).coeff(st[i],1).coeff(st[j],1))
# block structure: (W,Y) and (X,Z) decouple?
print('cross (W,Y)x(X,Z) zero:', all(sp.simplify(G[i,j])==0 for i in (0,1) for j in (2,3)))
print('(W,Y) block det:', sp.factor(G[:2,:2].det()), ' W pivot:', sp.factor(G[0,0]))
print('(X,Z) block det:', sp.factor(G[2:,2:].det()), ' X pivot:', sp.factor(G[2,2]))
# R6.5 residue at v=0 degenerates?
Q5=lambda vv: (Cs*(1+vv**2)*W**2 + As*(1+vv**2)*Y**2 + (Ds+half)*(W+vv*X)**2 + (K+half)*(vv*Y-Z)**2)
print('R6.5 residue at v0=0 has X absent:', sp.expand(Q5(0)).coeff(X,2)==0)
