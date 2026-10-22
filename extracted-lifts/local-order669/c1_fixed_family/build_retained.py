"""Export all retained elementary, circuit and angular roots on both sign charts."""
import json
from pathlib import Path
import sympy as S
import build_core as B
here=Path(__file__).resolve().parent
aa,bb,ww,zz,X,Y=S.symbols('quarter_u quarter_w fourth_W fourth_absZ X Y',real=True)
test=(aa,bb,ww,zz,X,Y)
u=aa**4;w=bb**4;W=ww**4;Za=zz**4;s=aa**2
src=B.v
globals().update(src)
Qb=1-b;theta=1-h;z=q*q;e=k*k*h/sigma
k1=2*(1+2*m)*(h-Qb)**2+8*D*h*b+8*p*b*sigma+4*beta*theta**2
all_payload={}
for sign in [1,-1]:
 Z=sign*Za
 Nx=X+u*h*Z;Ny=(Y-b*w*Z)/q;R=w*Z/q-u*q*W
 Gam=sigma*Z-u*w*W;At=q*(u*u-k*k)*W;Tt=k*q*(u-k)*W
 JC=Ny+b*R;JL=Ny-Qb*R;Xe=X+e*w*W;FA=u*(Z-a0*w*W)
 roots=[]
 def add(name,value,category): roots.append((name,S.expand(value),category))
 add('Eplus_beta_0',2*S.sqrt(beta)*Nx,'elementary')
 add('Eplus_beta_w',2*S.sqrt(beta)*w*Nx,'elementary')
 add('Eplus_Xe',S.sqrt(2*(2*D+1))*u*Xe,'elementary')
 add('Eplus_Gamma',2*S.sqrt(b*p)*q*u*Gam,'elementary')
 for name,mon in [('u',u),('w',w),('uw',u*w)]:
  add('Eplus_JC_'+name,2*S.sqrt(a)*mon*JC,'elementary')
 add('Eplus_JL_u',2*S.sqrt(ell)*u*JL,'elementary')
 add('Eplus_JC_Cw',2*S.sqrt(C)*q*w*JC,'elementary')
 add('Nx_w_reserve',S.sqrt(2+4*m)*w*Nx,'elementary')
 add('JC_excess',S.sqrt(4*a+2*z*(2*C+1)-1)*JC,'elementary')
 add('JL_excess',S.sqrt(4*ell-1)*JL,'elementary')
 add('Eyy_u',S.sqrt(2*(1+2*m))*q*u*Ny,'elementary')
 add('Eyy_w',S.sqrt(2*(1+2*D))*q*w*Ny,'elementary')
 add('ERR_u',S.sqrt(k1)*q*u*R,'elementary')
 add('ERR_w',S.sqrt(2*(2*D+1))*q*b*w*R,'elementary')
 add('FA',S.sqrt(L0/2)*FA,'elementary')
 add('angular_difference_u',S.sqrt((KA+c)/2)*q*u*(u-k)*W,'elementary')
 add('angular_difference_s',S.sqrt((KA+c)*k)*q*s*(u-k)*W,'elementary')
 # The eight angular roots for |Z|. For the negative chart, append the
 # three exact excess roots instead of dropping the change in state sign.
 Rp=w*Za/q-u*q*W;Gp=sigma*Za-u*w*W;FAp=u*(Za-a0*w*W)
 L=S.sqrt(q0*sigma)*FAp+k/S.sqrt(q0*sigma)*Gp
 Kroot=k*k*q*W-k*w*Za/q
 swWZ=bb**2*ww**2*zz**2
 angular=[At+u*Rp,S.sqrt(2)*k*Rp,L,s/S.sqrt(k)*L,
          S.sqrt(q0)*Kroot,S.sqrt(q0/k)*s*Kroot,
          S.sqrt(2*q0*k)*(u-k)*swWZ,
          S.sqrt(2*q0)*s*(u-k)*swWZ]
 if sign<0:
  angular.extend([2*k*S.sqrt(4+2*q0+1/q0)*s*swWZ,
                  2*S.sqrt(k)*S.sqrt(2*q0+2+1/q0)*u*swWZ,
                  2*u*s*swWZ])
 for i,val in enumerate(angular):
  add('angular_Delta_'+str(i+1),S.sqrt(c*dG/S0)*val,'angular')
 # Circuit roots are stored BEFORE their positive mixed scalar qC.
 # qC=1/2-c*dG/S0*(2/lam+k+u) is split by the fixed-test compiler.
 if sign>0:
  gs=S.sqrt(S.sqrt(lam)*k)
  circuit=[S.sqrt(2)*gs*((L0**S.Rational(1,4))*aa**3*bb*ww*zz**3
                       -(rho**S.Rational(1,4))*aa**3*bb**3*ww**3*zz),
           S.sqrt(2)*gs*((S.sqrt(lam)*k*k*q)**S.Rational(1,2)*aa*bb*ww**3*zz
                       -(S.sqrt(lam)*k/q)**S.Rational(1,2)*aa*bb**3*ww*zz**3),
           2*gs*((lam*k)**S.Rational(1,4)*aa**3*bb**2*ww**2*zz**2
                 -(lam*k**3)**S.Rational(1,4)*aa*bb**2*ww**2*zz**2)]
 else:
  circuit=[2*S.sqrt(lam)*k*s*swWZ]
 for i,val in enumerate(circuit): add('circuit_'+str(i+1),val,'circuit_unweighted')
 features=sorted({powers for _,val,_ in roots for powers,_ in S.Poly(val,*test).terms()})
 index={p:i for i,p in enumerate(features)}
 records=[]
 for name,val,cat in roots:
  records.append({'name':name,'category':cat,'entries':[
     [index[powers],S.sstr(S.factor(coef))] for powers,coef in S.Poly(val,*test).terms()]})
 all_payload[str(sign)]={'physical_state':'(W,X,Y,Z)=(fourth_W^4,X,Y,'+str(sign)+'*fourth_absZ^4)',
      'feature_exponents':[list(x) for x in features],'roots':records,
      'retained_roots_before_circuit_split':len(records),
      'maximum_source_testing_separation_rank':max(len(x['entries']) for x in records)}
data={'schema':'Q222.C1.retained_roots.v1','testing_variables':[str(x) for x in test],
      'testing_substitution':{'u':'quarter_u^4','w':'quarter_w^4','W':'fourth_W^4','absZ':'fourth_absZ^4'},
      'circuit_weight':'qC=1/2-c*dG/S0*(2/lambda+k+u)',
      'source_relations':'source_definitions.json, plus k1 as printed in P4/SP3','charts':all_payload}
(here/'C1_RETAINED_ROOTS.json').write_text(json.dumps(data,indent=2)+'\n')
print(json.dumps({sign:{'features':len(d['feature_exponents']),'roots':len(d['roots']),
 'max_rank':d['maximum_source_testing_separation_rank']} for sign,d in all_payload.items()},indent=2))
