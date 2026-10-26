import sympy as S
u,v,h,t=S.symbols('u v h t');W,Y,Z=S.symbols('W Y Z');U=u*u;V=v*v;du=1+U;dv=1+V;l=1-h
C=(h+t)/(2*l);A=t/(2*h);D=1/(2*(t-l))
X=v*W+u*v*Y-u*Z
J=u*W+u*v*X-Y-v*Z
P=C*dv*(u*W-Y)**2+A*du*dv*Y**2+(D+S.Rational(1,2))*J**2+S.Rational(1,2)*du*(v*Y-Z)**2+S.Rational(1,2)*dv*(X+u*Z)**2-S.Rational(1,2)*du*dv*Z**2
M=S.hessian(P,(W,Y,Z))/2
print('det computing')
det=S.factor(M.det())
claimed=(du**2*dv**2*(h+t)/(8*h*l*(t-l)))*(V**2+t*(h+t)*U**2*V+h*U*V**2)
print('K det',det)
print('RB24 claimed diff',S.factor(det-claimed))
