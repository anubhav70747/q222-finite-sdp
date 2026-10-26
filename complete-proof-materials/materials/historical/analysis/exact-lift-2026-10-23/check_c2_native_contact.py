from pathlib import Path
import sympy as s
z,y,H,L,r=s.symbols('z y H L r')
Bs=H+L; a=1+z*r; T=Bs+z*r; U=Bs-z*y*r; S=a+L; j=a-L/y
delta=z+(1+z)*Bs+(1+z)*Bs**2+2*z*(1+z)*Bs*r+z**3*r**2+2*z**2*r+y*((1+z)*Bs+z+z*z*r*(2-r))
G=Bs+z*(H+r-1); M=(1+z)*T*H+y*G+y*z*j+z*r*L*(1+y)
n=(T*a*S-r*delta)/(2*r*delta)
m=z*z*(1+y)*(1-r)*a*S/(2*G*delta)
C=a*M/(2*L*delta); A=y*a*T*S/(2*H*delta)
D=(1+z)*T*U*S/(2*y*j*delta); al=T*M/(2*(1-H)*delta)
x,w=s.symbols('x w')
P=n*(1+x)*(1+w)*r*r*x*x*w+m*w*(Bs+x*(r+H-1))**2+C*(1+w)*x*L*L+A*(1+x)*(1+w)*H*H*x+(D+s.Rational(1,2))*x*(L-w*(1+r*x))**2+(al+s.Rational(1,2))*(1+x)*x*w*(1-H)**2+s.Rational(1,2)*(1+w)*x*x*w*(1-r)**2-s.Rational(1,2)*(1+x)*(1+w)*x*w
# Cancel the known contact factors before expanding the common numerator.
yj=y*a-L
M=(1+z)*T*H+y*G+z*yj+z*r*L*(1+y)
terms=[
 (T*a*S-r*delta)*(1+z)*(1+y)*r*z*z*y,
 z*z*(1+y)*(1-r)*a*S*y*G,
 a*M*(1+y)*z*L,
 y*a*T*S*(1+z)*(1+y)*H*z,
 (1+z)*T*U*S*z*yj,
 delta*z*yj*yj,
 T*M*(1+z)*z*y*(1-H),
 delta*(1+z)*z*y*(1-H)**2,
 delta*(1+y)*z*z*y*(1-r)**2,
 -delta*(1+z)*(1+y)*z*y]
assert s.Poly(s.expand(sum(terms)),z,y,H,L,r).is_zero
print('Native C2 full-parent contact value is exactly zero.')
