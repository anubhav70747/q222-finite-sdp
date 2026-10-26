import sympy as s,time
now=time.time()
a,b,c,t,r=s.symbols('a b c t r');u,v=s.symbols('u v'); d,n,m,C,A,D,al,be,L=s.symbols('d n m C A D al be L')
z=c*c; Ac=1+z*r;h=b-a*a*t; g=(1+z)*b-(1-r+a*a*z*t);j=1+z*r-t
De=(1+z)*z*b*b+2*z*(1+z)*b*r+(1+z)*(1+a*a*z)*b+z*z*(1-a*a)*r*r+2*z*(1+a*a*z)*r+1+a*a*z
N=(1+z)*(b+r)*(b+a*a*(1-t))-a*a*(1+a*a*z)*(1-r)*t
weights={d:1,L:0,be:0,n:((b+r)*Ac*(Ac+a*a*z*t)-r*De)/(2*r*De),m:(1+a*a*z)*(1-r)*Ac*(Ac+a*a*z*t)/(2*g*De),C:Ac*N/(2*a*a*t*De),A:a*a*Ac*(b+r)*(Ac+a*a*z*t)/(2*h*De),D:(1+z)*(b+r)*(b-a*a*z*r)*(Ac+a*a*z*t)/(2*a*a*j*De),al:z*z*(b+r)*N/(2*(1-z*h)*De)}
Du=1+u*u;Dv=1+v*v
rows=[s.Matrix(zz) for zz in [(1,v,0,0),(0,1,0,0),(v,-1,u*v,-u),(u,0,-1,0),(0,0,1,0),(u,u*v,-1,-v),(0,0,v,-1),(0,1,0,u),(0,0,0,1)]]
coefs=[4*L*Du,4*n*Du*Dv,4*m,4*C*Dv,4*A*Du*Dv,4*D+2*d,(4*al+2*d)*Du,(4*be+2*d)*Dv,-2*d*Du*Dv]
P=sum([w*q*q.T for w,q in zip(coefs,rows)],s.zeros(4))
psi=s.Matrix([b/a,-c*r,c*(b/a-a*t),1]); angles={u:c,v:a*c}
Pt=s.Matrix([(s.factor(i.subs(angles).subs(weights))) for i in (P*psi)])
print('check P psi=',list(Pt),'time',time.time()-now,flush=True)
F=(psi.T*P*psi)[0]
print('jet u=',s.factor(s.diff(F,u).subs(angles).subs(weights)),'time',time.time()-now,flush=True)
print('jet v=',s.factor(s.diff(F,v).subs(angles).subs(weights)),'time',time.time()-now,flush=True)
