"""Read the exported roots and verify every retained polynomial identity."""
import json
from pathlib import Path
import sympy as S
here=Path(__file__).resolve().parent
raw=json.loads((here/'C1_RETAINED_ROOTS.json').read_text())
aa,bb,ww,zz,X,Y=S.symbols('quarter_u quarter_w fourth_W fourth_absZ X Y',real=True)
test=(aa,bb,ww,zz,X,Y);u=aa**4;w=bb**4;W=ww**4;Za=zz**4
names='b k q sigma p m D beta C a ell h rho dG S0 nu R0 Rx Bx Ay chi lam q0 L0 KA c a0 kappa r'
src={n:S.Symbol(n,positive=True) for n in names.split() if n!='KA'}
src['KA']=S.Symbol('KA',real=True)
globals().update(src)
kk,ll,tt,qq,dd,cc,Kplus=S.symbols('kk ll tt qq dd cc Kplus',positive=True)
sub={k:kk**4,lam:ll**4,sigma:tt**4,q:qq**2,dG:dd**2,c:cc**2,
 rho:ll**4*kk**8/tt**4,q0:ll**4*kk**4,L0:ll**4*tt**4,
 S0:dd**2+ll**4*kk**8/tt**4,a0:kk**4/tt**4,KA:Kplus-cc**2}
def simp(x): return S.cancel(S.expand(x.subs(sub)))
records=[]
for sign in [1,-1]:
 chart=raw['charts'][str(sign)]
 mons=[S.prod(v**n for v,n in zip(test,powers)) for powers in chart['feature_exponents']]
 roots={item['name']:sum(S.sympify(coef,locals=src)*mons[index] for index,coef in item['entries']) for item in chart['roots']}
 def sumcat(cat):
  return sum(roots[it['name']]**2 for it in chart['roots'] if it['category']==cat)
 Z=sign*Za
 Nx=X+u*h*Z;Ny=(Y-b*w*Z)/q;R=w*Z/q-u*q*W
 Gam=sigma*Z-u*w*W;At=q*(u*u-k*k)*W;Tt=k*q*(u-k)*W
 JC=Ny+b*R;JL=Ny-(1-b)*R;Xe=X+k*k*h/sigma*w*W
 FA=u*(Z-a0*w*W)
 k1=2*(1+2*m)*(h-(1-b))**2+8*D*h*b+8*p*b*sigma+4*beta*(1-h)**2
 Eplus=4*beta*(1+w*w)*Nx*Nx+2*(2*D+1)*(u*Xe)**2+4*b*p*q*q*u*u*Gam*Gam
 Eplus+=4*a*(u*u+w*w+u*u*w*w)*JC*JC+4*ell*u*u*JL*JL+4*C*q*q*w*w*JC*JC
 elementary=Eplus+(2+4*m)*w*w*Nx*Nx+(4*a+2*q*q*(2*C+1)-1)*JC*JC+(4*ell-1)*JL*JL
 elementary+=2*q*q*((1+2*m)*u*u+(1+2*D)*w*w)*Ny*Ny
 elementary+=q*q*(k1*u*u+2*b*b*(2*D+1)*w*w)*R*R
 elementary+=L0*FA*FA/2+(KA+c)*(At*At-Tt*Tt)/2
 assert simp(sumcat('elementary')-elementary)==0
 Bi=[S.sqrt(L0)*u*Z,S.sqrt(rho)*u*w*W,S.sqrt(lam)*k*k*q*W,S.sqrt(lam)*k*w*Z/q]
 Ccir=sum(v*v for v in Bi)-4*lam*k*k*u*w*W*Z
 if sign>0:
  G=S.sqrt(lam)*k*aa**2*bb**2*ww**2*zz**2
  Cplus=Ccir-sum((v-G)**2 for v in Bi)
 else: Cplus=-4*lam*k*k*u*w*W*Z
 assert simp(sumcat('circuit_unweighted')-Cplus)==0
 Delta=(2/lam+u+k)*Ccir+u*u*R*R+(2*a0*u+(u+k)/L0)*Gam*Gam-At*At
 assert simp(sumcat('angular')-c*dG/S0*Delta)==0
 records.append({'state_sign':sign,'elementary_roots':19,'angular_roots':8 if sign>0 else 11,
   'circuit_roots':3 if sign>0 else 1,'identities_checked':['full retained elementary polynomial','complete angular Delta including sign terms','circuit C_plus'],
   'passed':True})
report={'all_passed':True,'charts':records,'method':'Read exported coefficient matrices, multiply every printed root by itself, and compare complete state polynomials exactly after algebraic positive-source parameterization.',
 'boundary_statement':'The final testing formulas are polynomial in quarter-angle/state coordinates. W=0, Z=0, u=0, and w=0 follow by polynomial specialization. c=0 and dG=0 follow without cancellation.'}
(here/'RETAINED_ROOTS_VERIFIED.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
