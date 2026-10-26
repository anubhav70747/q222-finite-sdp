import sympy as s
u,v,w,k,ee,j,E=s.symbols('u v w k ee j E',real=True)
A=ee*w*v;B=v;C=E*u*w;D=ee**2*j**3*E
xx=B-C; yy=A-k*C; zz=u*B-k*C; hh=k*C-D
# require k=ee*j
subs={k:ee*j}
for tag,z0 in {
'3.1.a':C*(u-k)-(zz-u*xx),
'3.1.b':ee*C*(w-j)-(yy-ee*w*xx),
'3.1.c':C*(v-ee*j*j*E)-E*(w*zz+j*yy-k*w*xx),
'3.1.d':C*v*(w-j)-E*w*(w-j)*zz-E*j*w*yy+E*k*w*w*xx,
'3.3':k*C*(yy+zz-hh-2*k*xx)-((A*(u*B)-2*ee**2*j*j*B*C+D*(k*C)) - yy*zz),
'3.4':None
}.items():
 if z0 is not None:
  z1=s.expand(z0.subs(subs))
  assert z1==0,(tag,z1)
  print('PASS',tag)
# circuit identity from variances and Q_i squares
q1_sq=A*(u*B)+D*(k*C)-2*ee*ee*j*j*B*C
q2_sq=A*(k*C)+(u*B)*D-2*ee*ee*j*j*B*C
rhs=s.Rational(1,2)*((yy-zz)**2+hh**2+yy**2+(zz+hh)**2)+q1_sq+q2_sq
lhs=A*A+(u*B)**2+(k*C)**2+D*D-4*k*k*B*C
assert s.expand((rhs-lhs).subs(subs))==0
print('PASS', '3.4 circuit with two product squares')

U,V,Z=s.symbols('U V Z',real=True)
Q=U-V*Z;R=Z-U*V; bb=V-U*Z
Psi=U**2+V**2+Z**2+U**2*V**2*Z**2-4*U*V*Z
assert s.expand(Psi-((U-V)**2+Z**2*(1-U*V)**2+2*U*V*(1-Z)**2))==0
print('PASS','RFN12 circuit')
G7=s.Matrix([[4,-2,s.Rational(-3,2),0,s.Rational(-1,2),0,s.Rational(1,2)],[-2,3,0,s.Rational(-3,2),0,-1,2],[s.Rational(-3,2),0,4,s.Rational(3,2),-2,-3,0],[0,s.Rational(-3,2),s.Rational(3,2),4,s.Rational(-3,2),0,s.Rational(-3,2)],[s.Rational(-1,2),0,-2,s.Rational(-3,2),3,2,-1],[0,-1,-3,0,2,8,-4],[s.Rational(1,2),2,0,s.Rational(-3,2),-1,-4,4]])
G3=s.Matrix([[3,-2,-2],[-2,2,2],[-2,2,3]])
T7=s.Matrix([Z-U**2*V,V-U**2*V,U-U**2*V,U*Z-U**2*V,U*V-U**2*V,U*V*Z-U**2*V,U**2-U**2*V]);T3=s.Matrix([1-U*V,V-U*V,U-U*V]);
RHS=(T7.T*G7*T7)[0]+4*V*Z*(1-U)**2+U*Z*(T3.T*G3*T3)[0]
LHS=4*Psi+4*U**2*Q**2+4*U**2*R**2-(V*(U**2-1))**2
assert s.expand(LHS-RHS)==0,'RFN19 identity mismatch'
print('PASS RFN19 7+3 matrix identity')
print('G7 principal minors:',[G7[:i,:i].det() for i in range(1,8)])
print('G3 principal minors:',[G3[:i,:i].det() for i in range(1,4)])
assert all(G7[:i,:i].det()>0 for i in range(1,8))
assert all(G3[:i,:i].det()>0 for i in range(1,4))
expr=8*U**2*Psi+8*U**4*Q**2+(8*U**4+2*(U**2-1)**2)*bb**2-(V*(U**2-1))**2
rhs2=2*U**2*(4*Psi+4*U**2*Q**2+4*U**2*bb**2-(Z*(U**2-1))**2)+(U*Z*(U**2-1)-(U**2-1)*bb)**2
assert s.expand(expr-rhs2)==0
print('PASS RFN20')
