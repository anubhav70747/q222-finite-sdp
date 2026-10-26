import sympy as s
u,w,W,Z,z,y,r,B,H,n,m,C,A,D,alpha=s.symbols('u w W Z z y r B H n m C A D alpha',real=True)
b=s.symbols('b',positive=True); aH=s.symbols('aH',positive=True)
Du=1+u*u;Dw=1+w*w
X=-r*u*Z; Y=H*w*Z
Mm=W-w*X+u*Y-u*w*Z
J= u*w*W+u*X-w*Y-Z
Parent=s.expand(4*n*Du*Dw*X**2+4*m*Mm**2+4*C*Dw*(u*W-Y)**2+4*A*Du*Dw*Y**2+4*D*J**2+4*alpha*Du*(Y-w*Z)**2+2*(J**2+Du*(Y-w*Z)**2+Dw*(X+u*Z)**2-Du*Dw*Z**2))
# write sqrt B as b; sqrt H as aH; b^2 = B, aH^2=H; use b as actual sqrtB and aH as actual sqrtH
subsRoot={B:b*b,H:aH*aH}
R0=(u*W-b*b*w*Z)/b
R1=(z*W-b*b*u*w*Z)/b
R2=u*w*W-(b*b/y)*Z
T=(u*u-z)*Z
E=b*aH*(w*w-1/y)*Z
K=4*A*H**2+8*D*H*r+4*alpha*(1-H)**2+4*m*(H+r-1)**2+2*(H+r-1)**2+4*n*r*r
kk=K-4*B*B*m/z**2+8*z*n*r*r+8*A*H*H/y
IC7=(4*B*C-kk*z/B)*R0**2 +4*B*m/z**2*R1**2+2*(2*C+2*D+1)*R2**2+r*r*(4*n*w*w+2*(2*D+2*n+1))*T**2+(H/B)*(4*A*u*u+2*(2*A+2*C+2*D+1))*E**2 +kk*u/B*R0*R1 +(-4*r*(2*D+1)+kk/B)*R2*T -4*aH*(2*C+2*D+1)/b*R2*E
expr=s.cancel(s.expand(Parent.subs(subsRoot)-IC7.subs(subsRoot)))
# reduce aH^2 to H via aH^2 = H ; since use aH as free symbol, replace even powers
# The exact difference has precisely four monomials, independently of source stationarity.
p=s.Poly(expr,W,Z)
allowed={(1,1):{(1,1)},(0,2):{(2,0),(0,2),(0,0)}}
assert set(mon for mon,_ in p.terms()) <= set(allowed)
for mon,coeff in p.terms():
    try:
        exp=s.Poly(s.together(coeff),u,w)
    except Exception as e: raise AssertionError((mon,str(e)))
    assert {t for t,_ in exp.terms()} <=allowed[mon],(mon,exp.terms())
print('IC7 difference support PASS: only uwWZ, u^2 Z^2, w^2 Z^2, Z^2')
