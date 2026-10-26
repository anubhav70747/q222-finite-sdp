import json
from pathlib import Path
import sympy as s

out = {}
def check(name, expr):
    result = s.factor(expr)
    out[name] = result == 0
    if result != 0:
        out[name + '_residual'] = str(result)

u,w,W,X,Y,Z,q,k,h,b,sg = s.symbols('u w W X Y Z q k h b sg', nonzero=True)
m,c,a,ell,D,beta = s.symbols('m c a ell D beta')
nx,ny = s.symbols('nx ny')
z=q*q; x=k*k; du=1+u*u; dw=1+w*w
mw=W-w*X+u*Y-u*w*Z
rw=u*W-Y
jw=u*w*W+u*X-w*Y-Z
bw=X+u*Z
core=u*u*w*w*Z*Z-2*u*u*w*Y*Z-2*u*w*W*Z+2*u*w*w*X*Z
core+=u*u*w*w*W*W+2*u*u*w*W*X-2*u*w*w*W*Y
core+=(du+w*w)*(X*X+Y*Y)-2*u*w*X*Y
parent=4*m*z*mw**2+4*c*z*dw*rw**2+4*a/z*du*dw*Y**2
parent+=4*D*jw**2+4*ell/z*du*(Y-w*Z)**2+4*beta*dw*bw**2+2*core
expanded=s.Poly(s.expand(parent.subs({X:nx-k*h*Z,Y:q*ny+k*b*z*W})),nx,ny)
ju=u-k*b*z; pu=1+u*k*b*z; tu=k*h-u
A=4*m*z*z*tu*tu+4*ell*du+(4*beta+2)*z*tu*tu
B=(4*c*z+4*D+2)*ju*ju+4*a*du*x*b*b*z
C=(4*D*(1+u*k*h)**2+4*beta*(u-k*h)**2+2*du*x*h*h)/sg**2
D1=8*m*z*pu*tu-8*D*ju*(1+u*k*h)-8*ell*du*k*b-4*u*u*k*(h+b*z)-4*u+4*u*x*h*b*z
E=4*m*pu*pu+4*c*ju*ju+4*(a+ell)*du*x*b*b+2*du*x*b*b*z
Ax=-8*m*z*pu+(8*D+4)*u*ju
Bx=-8*D*u*(1+u*k*h)+8*beta*(u-k*h)-4*du*k*h
Cx=-8*m*z*tu+(8*beta+4)*(u-k*h)
Ay=8*m*z*u*pu-8*c*z*ju+8*(a+ell)*du*k*b+4*du*k*b*z
By=-8*c*z*ju+8*a*du*k*b-(8*D+4)*ju
Cy=8*m*z*z*u*tu-8*ell*du+8*D*z*(1+u*k*h)+4*u*z*tu
check('UH7_scalar',expanded.coeff_monomial(1)-(A/z*w*w*Z*Z+B*w*w*W*W+C*sg*sg*Z*Z+D1*w*W*Z+E*z*W*W))
check('UH7_first_Nx',expanded.coeff_monomial(nx)-(Ax*w*W+Bx*Z+Cx*w*w*Z))
check('UH7_first_Ny',expanded.coeff_monomial(ny)-(q*(Ay+By*w*w)*W+Cy*w*Z/q))
check('UH12_Nx2',expanded.coeff_monomial(nx**2)-(4*m*z*w*w+4*D*u*u+4*beta*dw+2*(du+w*w)))
check('UH12_Ny2',expanded.coeff_monomial(ny**2)-(4*m*z*z*u*u+4*c*z*z*dw+4*a*du*dw+4*ell*du+4*D*z*w*w+2*z*(du+w*w)))
check('UH12_NxNy',expanded.coeff_monomial(nx*ny)-(-8*m*z*u*q*w-8*D*u*q*w-4*u*q*w))

U=w*W-sg*Z/k; N=(w*Z-z*W/k)/q
S=q*(u-k)*W; T=sg*(u-k)*Z
check('UH11',T*T-sg*S*S-k*sg*(u-k)*S*N+k*(u-k)*T*U)

# Recompute the claimed matched angular Schur identity in the unit-sigma,
# z=0 finite-interior residue. This uses every printed scalar and first normal.
H=s.symbols('H',nonzero=True)
K=1+x; Es=1-x*H; J=1+sg*(x*H-1); theta=1-sg*H
lim={q:0,h:sg*H,m:K*x*H/(2*(K+sg)),c:K*Es/(2*x*(K+sg)),
     a:1/(2*b*(K+sg)),ell:x/(2*(1-x*b)*(K+sg)),
     D:K/(2*J*(K+sg))-s.Rational(1,2),
     beta:K/(2*theta*(K+sg))-s.Rational(1,2)}
def leading(expr): return s.factor(expr.subs(lim, simultaneous=True))
AA,BB,CC,DD,EE,AX,BX,AY,CY=map(leading,[A,B,C,D1,E,Ax,Bx,Ay,Cy])
L=s.factor(2*AA/k+2*BB*sg/k+DD)
J0=s.factor(CC-BB/k**2)
K0=s.factor(AA/k**2+2*BB*sg/k**2+DD/k+EE)
for name,expr in [('L',L),('J0',J0),('K0',K0)]:
    check('UH9_'+name,expr.subs(u,k))
check('UH9_angular_jet',(s.diff(K0,u)+sg*s.diff(J0,u)).subs(u,k))
check('UH9_first_Nx',(BX+sg*AX/k).subs(u,k))
check('UH9_first_Ny',(AY+CY/k).subs(u,k))
j=s.factor(s.diff(J0,u).subs(u,k))
j2=s.factor(s.diff(J0,u,2).subs(u,k)/2)
k2=s.factor(s.diff(K0,u,2).subs(u,k)/2)
hx=s.factor(s.diff(BX+sg*AX/k,u).subs(u,k)/sg)
hy=s.factor(s.diff(AY+CY/k,u).subs(u,k))
sn=s.factor(s.diff(L,u).subs(u,k)+k*sg*j)
normal_x=leading(4*D*x+4*beta+2*K)
normal_y=leading(4*(a+ell)*K)
def schur(p11,p12,p22,v1,v2):
    return s.factor((p22*v1*v1-2*p12*v1*v2+p11*v2*v2)/(p11*p22-p12*p12))
gamma_s=s.factor(k2-schur(normal_y,CY.subs(u,k)/2,AA.subs(u,k),hy/2,sn/2))
gamma_t=s.factor(j2-schur(normal_x,AX.subs(u,k)/2,BB.subs(u,k),hx/2,-k*j/2))
p=K*x*(K+4*sg)*H**2-2*K*(x+2*sg)*H+x
check('UH19_full_Schur',gamma_s+sg*gamma_t+2*p/(sg*(K+sg)))
def minimizer(p11,p12,p22,v1,v2):
    det=p11*p22-p12*p12
    return s.factor(-(p22*v1-p12*v2)/det),s.factor(-(-p12*v1+p11*v2)/det)
ny_s,n_s=minimizer(normal_y,CY.subs(u,k)/2,AA.subs(u,k),hy/2,sn/2)
nx_t,u_t=minimizer(normal_x,AX.subs(u,k)/2,BB.subs(u,k),hx/2,-k*j/2)
check('UH_unit_sigma_Ny_minimizer',ny_s-b*(1-H*x*K)/K)
check('UH_unit_sigma_N_minimizer',n_s-(1-H*x*K)/(x*K))
check('UH_unit_sigma_Nx_minimizer',nx_t-(sg*H*H*x*K-H*x*K-sg*H*(2*x+1)+x)/(sg*K))
check('UH_unit_sigma_U_minimizer',u_t-(H*x*K*(K+sg)-x*K-sg)/(x*sg*K))
check('UH_unit_sigma_D_row_defect',k*(nx_t+u_t)+Es/k-k*(H*K-1)*J/(sg*K))
check('UH_unit_sigma_beta_row_defect',nx_t+H+x*theta*(H*K-1)/(sg*K))
out['gamma_S']=str(gamma_s)
out['gamma_T']=str(gamma_t)
out['all_checked_identities_zero']=all(v is True for key,v in out.items() if key.startswith('UH') and not key.endswith('_residual'))
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
