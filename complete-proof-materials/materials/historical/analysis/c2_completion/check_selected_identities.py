"""Independent bounded checks of shared C2 mechanisms; no source checker imports."""
import json
import time
from pathlib import Path
import sympy as S

started = time.monotonic()
results = {}

def equal(name, expression):
    num = S.together(expression).as_numer_denom()[0]
    ok = S.expand(num) == 0
    results[name] = bool(ok)
    if not ok:
        raise AssertionError(name + ': ' + str(S.factor(num)))

U,V,Z,a=S.symbols('U V Z a')
Q=U-V*Z; B=V-U*Z; R=Z-U*V
psi=U**2+V**2+Z**2+U**2*V**2*Z**2-4*U*V*Z
equal('RSR13 circuit',psi-(U-V)**2-Z**2*(1-U*V)**2-2*U*V*(1-Z)**2)
t=U*V; d=U-V; E=Z*(1-t); T=t*(1-Z); J=T-d*d*Z
equal('RFN18c opposite budget',2*(psi+V**2*Q**2+U**2*B**2)-R**2-
      ((E+T)**2+2*d**2+4*t*(1-Z)**2+2*J**2+4*t*d**2*Z**2))
G7=S.Matrix([[4,-2,-S.Rational(3,2),0,-S.Rational(1,2),0,S.Rational(1,2)],
 [-2,3,0,-S.Rational(3,2),0,-1,2],[-S.Rational(3,2),0,4,S.Rational(3,2),-2,-3,0],
 [0,-S.Rational(3,2),S.Rational(3,2),4,-S.Rational(3,2),0,-S.Rational(3,2)],
 [-S.Rational(1,2),0,-2,-S.Rational(3,2),3,2,-1],[0,-1,-3,0,2,8,-4],
 [S.Rational(1,2),2,0,-S.Rational(3,2),-1,-4,4]])
G3=S.Matrix([[3,-2,-2],[-2,2,2],[-2,2,3]])
v7=S.Matrix([Z-U**2*V,V-U**2*V,U-U**2*V,U*Z-U**2*V,U*V-U**2*V,U*V*Z-U**2*V,U**2-U**2*V])
v3=S.Matrix([1-U*V,V-U*V,U-U*V])
equal('RFN19 radial Gram',4*psi+4*U**2*Q**2+4*U**2*R**2-(V*(U**2-1))**2-
      (v7.T*G7*v7)[0]-4*V*Z*(1-U)**2-U*Z*(v3.T*G3*v3)[0])
results['RFN19 G7 leading minors']=[str(G7[:i,:i].det()) for i in range(1,8)]
results['RFN19 G3 leading minors']=[str(G3[:i,:i].det()) for i in range(1,4)]
K=a/(a-1); TT=t*(Z-1)
equal('SH1 parameterized budget',a*psi+K*V**2*(Q**2+B**2)-R**2-
      (a*d**2+(a-1)*(E-TT/(a-1))**2+(2*a*t+K*V**4)*(1-Z)**2+2*K*Z*(t-V**2)**2))

# Reconstruct the WA matrix in unnormalized eta and scaled R2'=a R2.
u,v,aa,bb,tt,hh,rr,nn,mm,cc,AA,DD,al=S.symbols('u v a B t h r n m C A D alpha')
kap=1/tt; ell=1/aa; ll=1-hh; uu=u*u; vv=v*v
aw=(cc+DD+S.Rational(1,2))*ll**2+AA*hh**2
a4=AA*hh**2
b2=2*mm*hh+cc*ll**2+(AA+al+S.Rational(1,2))*hh**2
b4=hh**2*(mm+AA+al+S.Rational(1,2))
c2=nn*rr**2+(mm+S.Rational(1,2))*(1-rr)**2+al
c4=nn*rr**2; d2=(nn+S.Rational(1,2))*rr**2+2*DD*rr
d4=(nn+DD+S.Rational(1,2))*rr**2
bs=(bb*d4+kap*c4)/(kap*ell**2)
gs=b4+d4/(kap*ell)
j0=(bb*d4-a4*kap**2*ell+d2*ell)/(2*kap*ell**2)
j1=d4/(2*kap*ell); js=j0+j1*uu
fs=(b2-c4/ell**2-al/bb**2-bb*d4/(kap*ell**2))/2
SS=S.Matrix([[aw+a4*uu,u*v*js,-v*js],[u*v*js,mm+bs*uu,u*(fs-j1*uu)],[-v*js,u*(fs-j1*uu),al/bb**2+gs*uu]])
p0=(DD+rr*(nn+S.Rational(1,2)))/kap
p1=rr*(nn+DD+S.Rational(1,2))/kap
q1=(nn*rr+bb*p1)/ell
fx=-mm*hh+(DD+S.Rational(1,2))*ll-p0-q1
LX=S.Matrix([[u*v*(p0+p1*uu),-mm+q1*uu,u*(fx-p1*uu)]])
py=-(cc+DD+S.Rational(1,2))*ll+AA*hh; pyp=AA*hh
kv=-cc*ll+AA*hh+mm+(al+S.Rational(1,2))*hh
kuv=hh*(AA+mm+al+S.Rational(1,2))
kz=-mm*(1-rr)+rr*(DD+S.Rational(1,2))-al-S.Rational(1,2)
qy=(-kz-kap*pyp-bb*kuv)/ell
LY=S.Matrix([[py+pyp*uu,u*v*qy,v*(kv-qy+kuv*uu)]])
hxx=nn*(1+uu)*(1+vv)+mm+(DD+S.Rational(1,2))*uu*vv+(1+vv)/2
hxy=-u*v*(mm+DD+S.Rational(1,2))
hyy=AA*(1+uu)*(1+vv)+cc*(1+vv)+DD+S.Rational(1,2)+mm*uu*vv+(al+S.Rational(1,2))*(1+uu)*vv
GG=S.zeros(5); GG[0,0]=hxx; GG[0,1]=GG[1,0]=hxy; GG[1,1]=hyy
for i in range(3):
    GG[0,i+2]=GG[i+2,0]=LX[i]; GG[1,i+2]=GG[i+2,1]=LY[i]
    for j in range(3): GG[i+2,j+2]=SS[i,j]
scale=S.diag(1,1,1,1/aa,1); GP=scale*GG*scale
dcol=S.Matrix([u*v,-1,ll,rr*u*v,-aa*rr*v])
scol=S.Matrix([0,0,(aa*bb-uu)*v,-(bb-vv/tt)*u,uu-aa*vv/tt])
kcol=S.Matrix([-2*rr*tt*u,0,0,-rr**2*tt*u,aa*rr**2*tt])
HD=S.zeros(5)
for i,j,x in [(0,2,tt*u*v),(0,4,-tt*u),(2,3,rr*tt*u*v),(2,4,-aa*rr*tt*v)]: HD[i,j]=HD[j,i]=x
jgap=1+aa*bb*rr-ll/tt
res=GP.diff(DD)-dcol*dcol.T-(scol*kcol.T+kcol*scol.T)/2-jgap*HD
for i in range(5):
    for j in range(i,5): equal(f'RD4 entry {i+1}{j+1}',res[i,j])
g=1-aa*(1-rr-hh*bb)
Lm=S.Matrix([-1,u*v,0,1-g,hh*u]); HM=S.zeros(5)
for i,j,x in [(0,3,-1),(1,4,v),(3,3,2-g),(3,4,hh*u)]: HM[i,j]=HM[j,i]=x
res=GG.diff(mm)-Lm*Lm.T-g*HM
for i in range(5):
    for j in range(i,5): equal(f'R5 centered m entry {i+1}{j+1}',res[i,j])

# Critical-wall polynomial: derive its full four-state determinant and scalar section.
c,A,D,d0=S.symbols('c A D d', nonzero=True)
Bt,Y,Jt,Vt=S.symbols('B_t Y J_t V_t')
Cu=1+u*u; Cv=1+v*v; C4=4*D+2*d0
crit=4*c*Cv*Jt**2+4*A*Cu*Cv*Y**2+C4*(Jt+u*v*Bt-Cu*Vt)**2+2*d0*Cv*Bt**2+2*d0*Cu*(v*v*Y**2-2*Y*Vt-Vt**2)
critc=d0*(d0+2*A)*(d0+2*D)/(2*(4*A*D-d0*d0))
roots=[Bt,Y,Jt,Vt]
M=S.hessian(crit,roots)/2
expected=64*c*d0*Cu**2*Cv**2*(A*(2*D+d0)*u*u+D*(2*A+d0)*v*v)
equal('SEO20 full four-root determinant',(M.det(method='domain-ge')-expected).subs(c,critc))
y=d0/(2*A); bc=(4*A*D-d0*d0)/(A*C4)
sec=S.expand(crit.subs({Bt:0,Y:y*Vt,Jt:(bc+u*u)*Vt}))/Vt**2
gv=d0*(2*A+d0)/A; gam=D*d0*(2*A+d0)**2/(A*A*(2*D+d0))
equal('SEO13 exact zero section',(sec-(gv*u*u+gam*v*v+4*c*u**4+gv*(y+2)*u*u*v*v+4*c*u**4*v*v)).subs(c,critc))
equal('SEO17 exact outer identity',crit.subs(u,0)-crit.subs({u:0,v:0})-v*v*(2*d0*Bt**2+(4*A+2*d0)*Y**2+4*c*Jt**2))
equal('SEO18 exact outer identity',crit.subs(v,0)-crit.subs({u:0,v:0}).subs(Vt,Cu*Vt)-u*u*(4*A*Y**2+2*d0*Cu*Vt**2))

# Actual source in the unbounded-a chart; all poles canceled explicitly.
s,b,h,t,g0=S.symbols('s b h t g')
x=s*b; R0=1+h*x; K0=1+x; ws=R0-g0; r0=1-s*ws; l0=1-h
ca=1+x*r0; Sc=1+s*s*b*l0+x*r0; j0=ca-l0/t
DS=1+s*(1+x)+s**3*b*(1+x)+2*s*x*(1+x)*r0+s*s*b*b*r0*r0+2*s*b*r0+t*s**3*b*(1+x+b*r0*(2-r0))+t*s*s*b
QP=1+s*s*b*(1+t)+s*b*(h+2*r0+t)+s**3*b*b*(h+t)+s*s*b*b*(h*(r0+1)+t*(1+r0-r0*r0))
KP=(1+s)*(h+t)*(1+x); mass=ws*QP-b*KP
mg=(1+t*s*s*b)*ws*ca*Sc/(2*DS)
ns=s*s*mass/(2*r0*DS)
n1=s*K0*mass/(2*R0*DS); n2=r0*K0*K0*mass/(2*R0*R0*DS)
MS=(1+x)*(s+r0)*(h+t)-s*ws*(1+t*s*s*b)*l0
CS=ca*MS/(2*l0*DS)
Ah=t*ca*(s+r0)*Sc/(2*DS)
DSrc=(1+x)*(s+r0)*(1-t*x*r0)*Sc/(2*t*j0*DS)
alp_unit=(s+r0)*MS/(2*(1-h*s*s*b)*DS)
p=ws/R0; mu1=mg/R0; mu2=mg*g0/R0**2
F0=p*p/2-n2-2*n1*p+mu2
AN=p-2*mu1-2*n1
equal('RI5 F0 exact b face',F0.subs(b,0))
equal('RI5 AN exact b face',AN.subs(b,0))
equal('RI6 AN/(bs) g=0 factor s',S.cancel(AN.subs(g0,0)/b).subs(s,0))

def jet(expr, degree=2):
    """Exact Taylor jet in s,b around zero, unit parameters h,t retained."""
    out=0
    for i in range(degree+1):
        for j in range(degree+1-i):
            coeff=S.factor(S.diff(expr,s,i,b,j).subs({s:0,b:0})/S.factorial(i)/S.factorial(j))
            out += coeff*s**i*b**j
    return out

base={g0:0}
source_jets={name:jet(expr.subs(base)) for name,expr in [('A',Ah/h),('C',CS),('D',DSrc),('n_gamma',n1),('beta_D',DSrc*j0/ca)]}
expectA=t/(2*h)*(1-s+s*s)
expectC=(h+t)/(2*l0)-(t+1)*s/(2*l0)+(t+1)*s*s/(2*l0)
expectD=(1-s+s*s)/(2*(t-l0))-s*b*t*(h+t)/(2*(t-l0)**2)
expectn1=(s-s*s-s*b*(h+t))/2
for name,ex in [('A',expectA),('C',expectC),('D',expectD),('n_gamma',expectn1)]: equal('RI18 source jet '+name,source_jets[name]-ex)
for name,ex in [('C',t/2-(t+1)*s/2+(t+1)*s*s/2),('A',t*(1-s+s*s)/2),('beta_D',(1-s+s*s-(t+1)*s*b)/(2*t)),('n_gamma',(s-s*s-t*s*b)/2)]:
    expr=source_jets[name].subs(h,0) if name!='A' else jet(Ah.subs({g0:0,h:0}))
    equal('RU22 source jet '+name,expr-ex)

# Check the low weighted parts of EACH raw angular coefficient, rather than
# relying on cancellation after multiplication by dependent raw monomials.
z2=S.symbols('u_squared')
def low(poly, limit, unit_b=False):
    out=0
    variables=(s,z2) if unit_b else (s,b,z2)
    weights=(2,2) if unit_b else (2,2,4)
    pp=S.Poly(S.expand(poly),*variables)
    for powers,coef in pp.terms():
        if sum(p*w for p,w in zip(powers,weights))<limit:
            mon=coef
            for var,power in zip(variables,powers): mon*=var**power
            out+=mon
    return S.factor(out)
du=1+z2; cj=l0*(1-s*b)/t; ejinv=(1-l0*s*b)/t
aj=source_jets['A']; cjet=source_jets['C']; dj=source_jets['D']; nj=source_jets['n_gamma']
fbb=du**2*(cjet*cj**2+dj*(1-cj)**2+(1-cj)**2/2)+du*aj*cj**2+du*cj-du/2
fbc=(-2*du*aj*cj-du+2*nj*du)*ejinv
fcc=du*aj*ejinv**2
lam=(1-s+s*s+z2-(2*l0+h*(h+t))*s*b)/(2*h*t)
kap0=(h+t)/(2*t)
for name,actual,model in [('BB',fbb,lam+kap0*z2),('BC',fbc,-2*lam-4*kap0*s*b),('CC',fcc,lam+kap0*s*b)]:
    equal('RI raw '+name+' has weight >=6',low(actual-model,6))

ahj=t*(1-s+s*s)/2; ccj=t/2-(t+1)*s/2+(t+1)*s*s/2
bdj=(1-s+s*s-(t+1)*s*b)/(2*t); njr=(s-s*s-t*s*b)/2
cr=(1-s*b)/t; qr=1-cr; eir=(1-s*b)/t
fbbR=du**2*(ccj*cr**2+bdj*qr+qr**2/2)-ahj*du*2*z2*cr**2+du*cr-du/2
fbcR=(2*ahj*du*(z2-1)*cr-du+2*njr*du)*eir
fccR=2*ahj*du*eir**2
lamR=(1-s+s*s)/t
for name,actual,model in [('BB',fbbR,lamR+z2/2-(1/t+S.Rational(1,2))*s*b),('BC',fbcR,-2*lamR-z2/t+(3/t-1)*s*b),('CC',fccR,lamR+z2/t-2*s*b/t)]:
    equal('RU raw '+name+' has weight >=6 small b',low(actual-model,6))
    equal('RU raw '+name+' has weight >=4 unit b',low(actual-model,4,True))

eps2= b*t*(1+2*l0*s*b)
lamB=b*(1+z2)/(2*h)-s*b*(1+b*h*(h+t))/(2*h)
qB=b*(h+t)
for name,actual,model in [('BB',fbb*eps2,lamB+qB*z2/2),('BC',fbc*eps2,-2*lamB-2*qB*s*b),('CC',fcc*eps2,lamB+qB*s*b/2)]:
    equal('RB raw '+name+' has weight >=4 unit b',low(actual-model,4,True))
equal('RI fA leading pivot',S.diff(F0.subs(g0,0),b).subs({s:0,b:0})/t-kap0)
equal('RI varrho leading pivot',alp_unit.subs({s:0,b:0,g0:0})/t-kap0)
equal('RI AN bs leading coefficient',S.diff(AN.subs(g0,0),s,b).subs({s:0,b:0})-(h+t))
SB=2*(cjet*cj-(dj+S.Rational(1,2))*(1-cj)+cj*aj)+1
SC=-2*aj*ejinv
mu=(1-s+s*s)/h
equal('RI full M normal SB weight >=4',low(SB-mu,4))
equal('RI full M normal SC weight >=4',low(SC+mu,4))
SBR=2*(ccj*cr-(bdj+qr/2))+1; SCR=-2*ahj*eir
muR=1-s+s*s
equal('RU full M normal SB weight >=4',low(SBR-muR,4))
equal('RU full M normal SC weight >=4',low(SCR+muR,4))

results['elapsed_seconds']=round(time.monotonic()-started,3)
Path(__file__).with_suffix('.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
print(json.dumps({'identities':sum(v is True for v in results.values()),'all_passed':all(v is not False for v in results.values()),'seconds':results['elapsed_seconds']},indent=2))
