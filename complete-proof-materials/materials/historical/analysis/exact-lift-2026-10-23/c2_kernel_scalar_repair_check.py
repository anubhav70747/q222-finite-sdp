"""One exact univariate source family; no numerical solver or parameter sweep."""
from pathlib import Path
import json
import sympy as S

s,x,w=S.symbols("s x w",positive=True)
H=S.Rational(1,2); z=y=r=s; t=S.Rational(1,2); L=y*t; B=H+L
a=1+z*r; T=B+z*r; U=B-z*y*r; j=a-t; sig=a+L
Delta=S.expand(z+(1+z)*B+(1+z)*B**2+2*z*(1+z)*B*r+
               z**3*r**2+2*z**2*r+y*((1+z)*B+z+z**2*r*(2-r)))
G=B+z*(H+r-1)
Ms=(1+z)*T*H+y*G+y*z*j+z*r*L*(1+y)
n=S.factor((T*a*sig-r*Delta)/(2*r*Delta))
m=S.factor(z**2*(1+y)*(1-r)*a*sig/(2*G*Delta))
C=S.factor(a*Ms/(2*L*Delta))
A=S.factor(y*a*T*sig/(2*H*Delta))
D=S.factor((1+z)*T*U*sig/(2*y*j*Delta))
al=S.factor(T*Ms/(2*(1-H)*Delta))
k=S.factor(m*B**2/z**2)
b=A*H**2
nu=n+S.Rational(1,2)
g0=(C+D+S.Rational(1,2))*L**2+b*(1+x)
g2=D+(2*D*r+nu*r**2)*x+(D+nu)*r**2*x**2

# Independently reconstruct the full kernel restriction, before subtraction.
p=(n*(1+x)*(1+w)*r*r*x*x*w+m*w*(B+x*(r+H-1))**2
 +C*(1+w)*x*L*L+A*(1+x)*(1+w)*H*H*x
 +(D+S.Rational(1,2))*x*(L-w*(1+r*x))**2
 +(al+S.Rational(1,2))*(1+x)*x*w*(1-H)**2
 +S.Rational(1,2)*(1+w)*x*x*w*(1-r)**2
 -S.Rational(1,2)*(1+x)*(1+w)*x*w)
p1=S.Poly(p,w).coeff_monomial(w)
g1=S.cancel((p1-k*(x-z)**2)/x)
GG=S.factor(g2.subs(x,z))
gp=S.factor(S.diff(g2,x).subs(x,z))
f=(D+nu)*r**2
c=n*r**2
q0=S.factor(4*y**2*GG*f+4*y*GG*c-(b/y-y*gp)**2)
q1=S.factor(4*b*f+2*c*(b/y+y*gp))
q2=-c**2
Q=q0+q1*(x-z)+q2*(x-z)**2

checks={
 "kernel_identity":S.cancel(p-k*w*(x-z)**2-x*(g2*w*w+g1*w+g0))==0,
 "contact_g0":S.cancel(g0.subs(x,z)-y*y*GG)==0,
 "contact_g1":S.cancel(g1.subs(x,z)+2*y*GG)==0,
 "contact_first_x":S.cancel(S.diff(g0,x).subs(x,z)+
               y*S.diff(g1,x).subs(x,z)+y*y*gp)==0,
 "g1_quadratic":S.cancel(S.Poly(g1,x).coeff_monomial(x*x)-c)==0,
 "det_factor":S.cancel(4*g0*g2-g1*g1-(x-z)**2*Q)==0,
}
det_at_2s=S.factor(s*s*(q0+s*q1-s*s*c*c)/4)
out={"family":{"H":str(H),"z":"s","y":"s","r":"s","t":str(t),
                "L":str(L),"B":str(B),"j":str(j),"Delta":str(Delta),
                "G":str(G),"Ms":str(S.factor(Ms))},
     "coefficients":{name:str(value) for name,value in
        [("n",n),("m",m),("C",C),("A",A),("D",D),("alpha",al),("k",k)]},
     "checks":checks,"q0":str(q0),"q1":str(q1),"q2":str(S.factor(q2)),
     "det_at_x_2s":str(det_at_2s)}
for name,expr in [("q0",q0),("det_at_x_2s",det_at_2s)]:
    num,den=S.fraction(expr)
    out[name+"_num_coeff_ascending"]=[str(v) for v in reversed(S.Poly(num,s).all_coeffs())]
    out[name+"_den"]=str(S.factor(den))
path=Path(__file__).with_suffix(".json")
path.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
assert all(checks.values()),checks
print(json.dumps(out,indent=2))


# One complementary exact small-H family for the proposed linear correction.
h=S.symbols("h",positive=True)
hh=h; zz=yy=h*h; rr=tt=S.Rational(1,2); ll=yy*tt; bb=hh+ll
aa=1+zz*rr; JJ=aa-tt; TT=bb+zz*rr; UU=bb-zz*yy*rr; SS=aa+ll
de=zz+(1+zz)*bb+(1+zz)*bb**2+2*zz*(1+zz)*bb*rr+zz**3*rr**2+2*zz**2*rr+yy*((1+zz)*bb+zz+zz**2*rr*(2-rr))
gg=bb+zz*(hh+rr-1)
mm=(1+zz)*TT*hh+yy*gg+yy*zz*JJ+zz*rr*ll*(1+yy)
nn=S.cancel((TT*aa*SS-rr*de)/(2*rr*de))
mc=S.cancel(zz**2*(1+yy)*(1-rr)*aa*SS/(2*gg*de))
ac=S.cancel(yy*aa*TT*SS/(2*hh*de))
dc=S.cancel((1+zz)*TT*UU*SS/(2*yy*JJ*de))
kc=S.cancel(mc*bb**2/zz**2)
nuu=nn+S.Rational(1,2)
ff=(dc+nuu)*rr**2; cc=nn*rr**2; bbb=ac*hh**2
G2=dc+(2*dc*rr+nuu*rr**2)*x+ff*x*x
ggc=G2.subs(x,zz); gpc=S.diff(G2,x).subs(x,zz)
BT=bbb/yy; VT=yy*gpc; T0=yy*ggc
sc=(BT+VT)**2/T0+cc+yy*ff
q0c=4*yy*yy*ggc*ff-(BT-VT)**2+4*yy*ggc*(cc+sc)
q1c=4*bbb*ff+2*(cc+sc)*(BT+VT)
limits={}
for name,val in [
    ("n",nn),("k0",kc),("A_over_h",ac/h),("hD",h*dc),
    ("B_over_h",BT/h),("T0_over_h",T0/h),("V_over_h",VT/h),
    ("yf_over_h",yy*ff/h),("c",cc),("s_correction",sc),
    ("corrected_q0_over_h",q0c/h),("corrected_q1_over_h",q1c/h)]:
    limits[name]=str(S.limit(S.cancel(val),h,0,dir="+"))
expected={"n":"1/2","k0":"1/4","A_over_h":"1/2","hD":"1",
          "B_over_h":"1/2","T0_over_h":"1","V_over_h":"1",
          "yf_over_h":"1/4","c":"1/8","s_correction":"1/8",
          "corrected_q0_over_h":"1","corrected_q1_over_h":"3/4"}
assert limits==expected,(limits,expected)
out["proposed_correction_small_H_limits"]=limits
path.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
print(json.dumps({"proposed_correction_small_H_limits":limits},indent=2))


