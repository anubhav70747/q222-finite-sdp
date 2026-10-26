"""Independent complete parent coefficient comparisons for BN SN9/SN15/SN18 and UN UN9/UN10.
No tests by sampling; symbolic polynomials in arbitrary source coefficients."""
import sympy as s

def ck(tag,x):
    xx=s.factor(s.cancel(x));assert xx==0,(tag,xx);print('PASS',tag)
# Independent source parameters (without source contact equalities, necessary for original-parent expansion).
S,u,w,W,Z,r,e,h=s.symbols('S u w W Z r e h',nonzero=True)
m,c0,a0,D,alpha,beta=s.symbols('m c0 a0 D alpha beta')
tn=s.symbols('tn',real=True) # pure-source sqrt(S*r*(1-r))
Du=1+u**2;Dw=1+w**2;aa=u**2;cc=w**2;d=D+s.Rational(1,2)
C=c0/(S*(1-r));A=a0/(S*r)

def parent(X,Y):
    return (m*(W-w*X+u*Y-u*w*Z)**2+C*Dw*(u*W-Y)**2+A*Du*Dw*Y**2+
    d*(u*w*W+u*X-w*Y-Z)**2+(alpha+s.Rational(1,2))*Du*(Y-w*Z)**2+
    (beta+s.Rational(1,2))*Dw*(X+u*Z)**2-s.Rational(1,2)*Du*Dw*Z**2)
# BN: X_e and Ny independent
Xe,Ny=s.symbols('Xe Ny')
f=parent(Xe-e*w*W,r*u*W+tn*Ny)
pol=s.Poly(s.expand(f), Xe,Ny)
coef=lambda i,j:pol.coeff_monomial(Xe**i*Ny**j)
red=lambda p:s.cancel(s.together(p).subs(tn**2,S*r*(1-r)))
Hx=m*cc+d*aa+(beta+s.Rational(1,2))*Dw
Hy=r*c0*Dw+(1-r)*a0*Du*Dw+S*r*(1-r)*(m*aa+d*cc+(alpha+s.Rational(1,2))*Du)
Hxy=-u*w*tn*(m+d)
ck('SN18 H_X', red(coef(2,0)-Hx))
ck('SN18 H_Y', red(coef(0,2)-Hy))
ck('SN18 H_XY', red(coef(1,1)-2*Hxy))
AN=m*(1+r*aa+e*cc)**2+(beta+s.Rational(1,2))*e**2*cc*Dw
AP= Dw*(c0*(1-r)+a0*r*Du)+S*(d*(1-r-e)**2*cc+(alpha+s.Rational(1,2))*Du*r**2)
KK=m*(1+r*aa+e*cc)+d*(1-r-e)+(alpha+s.Rational(1,2))*Du*r+(beta+s.Rational(1,2))*e*Dw
FZ=D+alpha*cc+beta*aa+(m+alpha+beta+s.Rational(1,2))*aa*cc
G0=AN*W**2+aa/S*AP*W**2-2*KK*u*w*W*Z+FZ*Z**2
ck('SN9 complete scalar', red(coef(0,0)-G0))
AX=-m*(1+r*aa+e*cc)+d*aa*(1-r-e)-e*(beta+s.Rational(1,2))*Dw
TX=u*(beta-D+cc*(m+beta+s.Rational(1,2)))
ck('SN15 full X first normal', red(coef(1,0)-2*(w*AX*W+TX*Z)))
YP=Dw*(-c0+a0*Du)+S*(m*(1+r*aa+e*cc)-d*cc*(1-r-e)+(alpha+s.Rational(1,2))*Du*r)
YC=s.symbols('ell')*YP-m*aa+d-(alpha+s.Rational(1,2))*Du
LYgeneric=tn*(YP*u*W/S + w*Z*(-m*aa+d-(alpha+s.Rational(1,2))*Du))
ck('SN18 full Y first normal before conductor', red(coef(0,1)-2*LYgeneric))
# UN: X=-huZ+Nx
Nx=s.symbols('Nx'); theta=1-h
ff=parent(Nx-h*u*Z,r*u*W+tn*Ny)
po=s.Poly(s.expand(ff),Nx,Ny)
cs=lambda i,j:po.coeff_monomial(Nx**i*Ny**j)
ck('UN hard Nx',red(cs(2,0)-Hx))
ck('UN hard Ny',red(cs(0,2)-Hy))
ck('UN hard NxNy',red(cs(1,1)-2*Hxy))
APP=Dw*(c0*(1-r)+a0*r*Du)+S*(d*(1-r)**2*cc+(alpha+s.Rational(1,2))*Du*r**2)
KUN=m*theta*(1+r*aa)+d*(1-r)*(1+h*aa)+(alpha+s.Rational(1,2))*Du*r
n=beta*theta**2+2*D*h+h**2/2
FZUN=D+alpha*cc+n*aa+(m*theta**2+alpha+(beta+s.Rational(1,2))*theta**2)*aa*cc+d*h**2*aa**2
FU=m*(1+r*aa)**2*W**2+APP*aa/S*W**2-2*KUN*u*w*W*Z+FZUN*Z**2
ck('UN9 complete scalar',red(cs(0,0)-FU))
LXun=-m*w*((1+r*aa)*W-theta*u*w*Z)+d*u*(u*w*(1-r)*W-(1+h*aa)*Z)+(beta+s.Rational(1,2))*Dw*u*theta*Z
# X derivative includes contribution from m first term and d and beta; X =Nx-huZ
ck('UN X first normal raw',red(cs(1,0)-2*LXun))
Ypun=Dw*(-c0+a0*Du)+S*(m*(1+r*aa)-d*(1-r)*cc+(alpha+s.Rational(1,2))*Du*r)
LYun=tn*(Ypun*u*W/S + w*Z*(-m*theta*aa+d*(1+h*aa)-(alpha+s.Rational(1,2))*Du))
ck('UN17/18 Y first normal raw',red(cs(0,1)-2*LYun))
print('ALL full-parent checks passed')
