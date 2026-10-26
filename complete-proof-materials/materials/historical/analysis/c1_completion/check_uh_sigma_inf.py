import runpy,contextlib,io,json
from pathlib import Path
import sympy as s
with contextlib.redirect_stdout(io.StringIO()):
    d=runpy.run_path(str(Path(__file__).parents[1]/'check_r21_c1_independent.py'))
globals().update({k:v for k,v in d.items() if not k.startswith('__')})
out={}
def ck(name,expr):
    val=s.factor(expr); out[name]=val==0
    if val!=0: out[name+'_residual']=str(val)
zz=q*q; xx=k*k; KK=1+xx; QQ=1-xx*b; RR=1+xx*zz*b
hh=(xx+zz*QQ)/(xx*(KK+zz)); ee=RR/(KK+zz); MM=1-zz*QQ
bs=s.symbols('bs')
lim_inf={h:sg*hh,D:sg**2*(KK*hh*ee/2-xx*bs),beta:sg**2*bs,
 m:xx*hh/(2*MM),c:RR/(2*xx*(1-zz*b)*(KK+zz)),
 a:hh/(2*b),ell:xx*RR/(2*QQ*(KK+zz))}
def lead(expr):
    return s.factor(s.cancel(expr.subs(lim_inf,simultaneous=True)).subs(sg,0))
AA,BB,CC,DD,EE,AX,BX,CX,AY,CY=map(lead,[A,B,C,D1,E,Ax,Bx,Cx,Ay,Cy])
LL=s.factor(2*AA/k+DD)
JJ=s.factor(CC-BB/k**2)
KO=s.factor(AA/k**2+DD/k+EE)
for nm,ex in [('L',LL),('J0',JJ),('K0',KO)]:ck('UH9_sigma_inf_'+nm,ex.subs(u,k))
ck('UH16_Nx',lead(4*D*xx+4*beta+2*KK)-2*KK)
ck('UH16_cross',(AX+zz*CX/k).subs(u,k)-4*KK*xx*hh)
ck('UH16_U',BB.subs(u,k)-2*KK*xx*hh)
HY=lead(4*m*zz**2*xx+4*c*zz**2+4*(a+ell)*KK+2*zz*KK)
hy=s.factor(s.diff(AY+CY/k,u).subs(u,k))
sn=s.factor(s.diff(LL,u).subs(u,k))
ks=s.factor(s.diff(KO,u,2).subs(u,k)/2)
gamma=s.factor(ks-schur(HY,CY.subs(u,k)/2,AA.subs(u,k),hy/2,sn/2))
ck('UH15',gamma-8*RR*hh/KK)
nys,ns=minimizer(HY,CY.subs(u,k)/2,AA.subs(u,k),hy/2,sn/2)
ck('UH_sigma_inf_Ny_minimizer',nys-b*(2*b*xx*zz-xx+1)/KK)
ck('UH_sigma_inf_N_minimizer',ns-(4*b*xx*zz-xx-2*zz+1)/(xx*KK))
et,ft,lt=s.symbols('et ft lt')
gg=ft-lt; du=1+u*u
m0=xx*hh/(2*MM); c0=RR/(2*xx*(1-zz*b)*(KK+zz))
a0=hh/(2*b); l0=xx*RR/(2*QQ*(KK+zz))
raw=4*m0*(et+u*zz*lt)**2+4*c0*(u*et-zz*ft)**2+4*a0*du*ft**2
raw+=4*l0*du*lt**2+2*zz*(u*u*lt**2+ft**2)-4*u*et*gg
mat=s.hessian(raw,(et,ft,lt))/2
ck('UH17_hard_cross',mat[1,2])
sc=s.factor(mat[0,0]-mat[0,1]**2/mat[1,1]-mat[0,2]**2/mat[2,2])
expect=2*(u*u-xx)**2*(1+u*u)*RR*hh/((u*u*(1-zz*b)+RR)*(u*u*(xx+zz*QQ)+xx*MM))
ck('UH17_full',sc-expect)
ck('UH17_minimizer_F',-mat[0,1]/mat[1,1]-u*b*KK/(u*u*(1-zz*b)+RR))
ck('UH17_minimizer_F_minus_G',-mat[0,2]/mat[2,2]+u*KK*QQ/(u*u*(xx+zz*QQ)+xx*MM))
out['all_zero']=all(v is True for key,v in out.items() if not key.endswith('_residual'))
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2))
