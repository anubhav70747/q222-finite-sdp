"""Exact identities for the physical witnesses; the planar chord proof is analytical."""
import sympy as s
u,v,x,y,w=s.symbols('u v x y w',real=True)
z=[u**4*x*x,2*u**4*x*y,u**4*y*y,u**4*w*w,
   2*u**3*v*x*w,2*u**3*v*y*w,u*u*v*v*x*x,2*u*u*v*v*x*y,
   u*u*v*v*y*y,u*u*v*v*w*w,2*u*v**3*x*w,2*u*v**3*y*w,
   v**4*x*x,2*v**4*x*y,v**4*y*y,v**4*w*w]
a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p=z
U=s.Matrix([[a,b/2,e/2,g,h/2],[b/2,c,f/2,h/2,i],
 [e/2,f/2,j,k/2,l/2],[g,h/2,k/2,m,n/2],[h/2,i,l/2,n/2,o]])
V=s.Matrix([[d,e/2,f/2,j],[e/2,g,h/2,k/2],[f/2,h/2,i,l/2],[j,k/2,l/2,p]])
rp=s.Matrix([u*u*x,u*u*y,u*v*w,v*v*x,v*v*y])
rm=s.Matrix([u*u*w,u*v*x,u*v*y,v*v*w])
assert all(s.expand(q)==0 for q in U-rp*rp.T)
assert all(s.expand(q)==0 for q in V-rm*rm.T)
ell=a+c+d+2*g+2*i+2*j+m+o+p
M=2*b+4*h+2*n
S=2*a+2*c-2*d+8*i-8*j+8*l-2*m-2*o+2*p
assert s.expand(ell-(u*u+v*v)**2*(x*x+y*y+w*w))==0
assert s.expand(M-4*x*y*(u*u+v*v)**2)==0
overlap=u*u-v*v
physicalS=2*overlap*x*x+2*(1+overlap-overlap**2)*(y*y-w*w)+8*u*v*(1-overlap)*y*w
assert s.rem(s.Poly(s.expand(S-physicalS),u),s.Poly(u*u+v*v-1,u)).as_expr()==0
assert s.expand(S-(2*(u**4-v**4)*x*x+2*(u**4+4*u*u*v*v-v**4)*(y*y-w*w)+16*u*v**3*y*w))==0
assert S.subs({u:1,v:0,x:0,y:0,w:1})==-2
assert S.subs({u:0,v:1,x:0,y:0,w:1})==2
t,X,Y=s.symbols('t X Y',real=True)
zp=s.Matrix([x,y,t*w,t*t*x,t*t*y]);zm=s.Matrix([w,t*x,t*y,t*t*w])
monomials=set()
for expr in list(zp*zp.T)+list(zm*zm.T):
 monomials.update(s.Poly(expr,t,x,y,w).monoms())
assert len(monomials)==16
motzkin=1+X**4*Y**2+X**2*Y**4-3*X*X*Y*Y
assert s.expand(motzkin.subs(Y,X)-(X*X-1)**2*(2*X*X+1))==0
print('Verified exact shared-entry witnesses, normalization, physical coordinates, both endpoints, 16-dimensional feature span and equal-angle Motzkin factorization.')
print('Scope: algebraic identities only. Convex-hull equality and the two-strategy/rank result use the analytical proofs.')
