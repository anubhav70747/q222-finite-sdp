import sympy as S
u,w,v,eps,j,E=S.symbols('u w v eps j E');k=eps*j
A=eps*w*v;B=v;C=E*u*w;D=eps**2*j**3*E
x=B-C;y=A-k*C;z=u*B-k*C;h=k*C-D
T_u,U,V,T_v=S.symbols('Tu U V Tv')
relations=[C*(u-k)-(z-u*x),eps*C*(w-j)-(y-eps*w*x),C*(v-eps*j*j*E)-E*(w*z+j*y-k*w*x),C*v*(w-j)-E*w*(w-j)*z-E*j*w*y+E*k*w*w*x]
for i,q in enumerate(relations):print('RI14',i,S.factor(q)==0)
Sleft=C*((u-k)*T_u+(w-j)*(eps*U+v*V)+(v-eps*j*j*E)*T_v)
Sright=(x*(-u*T_u-eps*w*U-E*k*w*w*V-E*k*w*T_v)+y*(U+E*j*w*V+E*j*T_v)+z*(T_u+E*w*(w-j)*V+E*w*T_v))
print('RI15',S.factor(Sleft-Sright)==0)
A1,A2,A3,A4=S.symbols('a1 a2 a3 a4')
Q1=A1*A2-A3*A4;Q2=A1*A3-A2*A4
YY=A1**2-A3**2;ZZ=A2**2-A3**2;HH=A3**2-A4**2
lhs=A1**4+A2**4+A3**4+A4**4-4*A1*A2*A3*A4
rhs=((YY-ZZ)**2+HH**2+YY**2+(ZZ+HH)**2)/2+Q1**2+Q2**2
print('T4 strict circuit:',S.factor(lhs-rhs)==0)
print('separated monomial product:',S.factor(A*(u*B)*(k*C)*D-k**4*B**2*C**2)==0)
print('gradient determinant candidate:',S.factor(S.det(S.Matrix([[S.diff(q,c) for c in (u,w,v)] for q in (x,y,z)]).subs({u:k,w:j,v:eps*j*j*E}))))
