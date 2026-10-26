import sympy as s
x,r,h,z,eta=s.symbols('x r h z eta')
K=1+x;R=1+x*r;T=1+x*h;G=h+z;theta=1-h
E=1-x*eta*h;y=x*eta*z
p0=1+2*x*h+r*x*x*G+x*z
p1=x*(K*z+x*h*theta)
Pd=p0+eta*p1
B5=K*(1+y)*(r+eta)-(eta*T-1+r)*Pd
Bb=K*x*(r+eta)*G-theta*Pd
Acal=theta*T+K*z;Ccal=R*G*(1+r*x*h)
Bcal=h*h*r*x*x+h*r*x*(2-h)+h*h*x+h*r*x*x*z+h*x*z+h*x+h+r*x*z+x*z+2*z+1
assert s.expand(h*B5+Acal*E**2-Bcal*E+Ccal)==0
a2=K*x*z-T*p1
a1=K*(1+x*z*r)-T*p0+(1-r)*p1
a0=K*r+(1-r)*p0
b1=K*x*G-theta*p1
b0=K*x*r*G-theta*p0
lhs=b1*b1*B5-a2*Bb**2-(b1*a1-2*a2*b0)*Bb
rhs=b1*b1*a0-b1*a1*b0+a2*b0*b0
assert s.expand(lhs-rhs)==0
assert s.expand(rhs-x*x*h*h*K*K*R*G*G*(1+2*x*h+r*x*x*h*h+(1-r)*x*z))==0
print('NC6 and NC7 polynomial identities: PASS')
