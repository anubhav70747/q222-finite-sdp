from pathlib import Path
import sympy as s
x,w,z,y,H,L,r,n,m,C,A,D,al=s.symbols('x w z y H L r n m C A D alpha', positive=True)
B=H+L
p=(n*(1+x)*(1+w)*r*r*x*x*w+m*w*(B+x*(r+H-1))**2
 +C*(1+w)*x*L*L+A*(1+x)*(1+w)*H*H*x
 +(D+s.Rational(1,2))*x*(L-w*(1+r*x))**2
 +(al+s.Rational(1,2))*(1+x)*x*w*(1-H)**2
 +s.Rational(1,2)*(1+w)*x*x*w*(1-r)**2
 -s.Rational(1,2)*(1+x)*(1+w)*x*w)
q=s.Poly(s.expand(p),w)
lines=['Native C2 contact-normal kernel, unnormalized F/4.',
'Normal kernel vector psi=(H+L,-r*u^2*v,H*u,u*v).',
'For x=u^2,w=v^2, all state normals N1,N2,N3 vanish exactly.',
'Pkernel=p2(x)w^2+p1(x)w+p0(x):']
for power in (2,1,0):
    lines += [f'p{power}(x) = {s.factor(q.coeff_monomial(w**power))}']
# stationary/zero at z,y forces p1(z)=-2y p2(z), p0(z)=y^2p2(z),
# p0'(z)+y p1'(z)+y^2p2'(z)=0.
# exact polynomial identity with divdiff, valid only imposed contact equations.
p0=q.coeff_monomial(1); p1=q.coeff_monomial(w); p2=q.coeff_monomial(w*w)
lines += [
'The nonnegative w^2 coefficient is x times a degree-two polynomial:',
s.sstr(s.Poly(s.expand(p2/x),x)),
'The constant coefficient is x*((C+D+1/2)*L^2 + A*H^2*(1+x)).',
'Both p0,p2 have a common x factor. At x=0 the residual is m*(H+L)^2*w, so the w coefficient alone survives; do not divide by x on that axis.',
'The exact coefficient array is included below; no universal PSD assertion is made.',
s.sstr(s.Poly(s.expand(p),x,w))]
out=Path(__file__).with_name('c2_contact_kernel_reduction.txt')
out.write_text('\n\n'.join(lines)+'\n')
print('\n'.join(lines[:9]))

